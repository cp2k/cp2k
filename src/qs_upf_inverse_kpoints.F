!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse Bloch density and operator maps on native Quickstep k-point groups.
!>        AO contractions are replicated. Point integrals use the
!>        native PW communicator; group density matrices contribute from one group rank only.
!>        It does not select a representation or replace the SCF operator provider.
! **************************************************************************************************
MODULE qs_upf_inverse_kpoints
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_submatrix
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_prepare
   USE qs_upf_inverse_gaussian,         ONLY: qs_upf_inverse_gaussian_prepare
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_collocate_dense,&
                                              qs_upf_inverse_grid_integrate_dense,&
                                              qs_upf_inverse_grid_metric_dense,&
                                              qs_upf_inverse_grid_prepare
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_gaussian_release,&
                                              qs_upf_inverse_kpoints_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_kpoints_type, qs_upf_inverse_kpoints_prepare, qs_upf_inverse_kpoints_release
   PUBLIC :: qs_upf_inverse_kpoints_density, qs_upf_inverse_kpoints_collocate
   PUBLIC :: qs_upf_inverse_kpoints_integrate, qs_upf_inverse_kpoints_metric, qs_upf_inverse_kpoints_distribute
   PUBLIC :: qs_upf_inverse_kpoints_release_operators

CONTAINS

! **************************************************************************************************
!> \brief Prepare each actual Bloch vector on the same physical-cell PW grid.
!>        All world ranks participate in the same point-integral sequence. A full unreduced
!>        mesh is required until the inverse density's space-group rotations are connected.
!> \param ks_env Live native KS environment with initialized k-point groups and bases.
!> \param grid Native physical-cell grid.
!> \param images Fixed physical-cell images of the Gaussian and original radial fields.
!> \param frame New owned Bloch fields; the native environment is unchanged.
!> \param info Factorization status; nonzero releases the partially prepared fields.
!> \param atom_local Optional nonoverlapping one-center inverse.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kpoints_prepare(ks_env, grid, images, frame, info, atom_local, direct)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      INTEGER, INTENT(IN)                                :: images(:, :)
      TYPE(qs_upf_inverse_kpoints_type), INTENT(OUT)     :: frame
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct

      INTEGER                                            :: ik
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL get_ks_env(ks_env, kpoints=kpoints)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(kpoints%nkp > 0)
      IF (kpoints%symmetry .OR. .NOT. kpoints%full_grid) THEN
         CPABORT("USPP Bloch calculations require a full unreduced mesh")
      END IF
      CPASSERT(ASSOCIATED(kpoints%kp_env))
      CPASSERT(ASSOCIATED(kpoints%para_env_kp))
      CPASSERT(ALL(kpoints%wkp >= 0.0_dp))
      CPASSERT(ABS(SUM(kpoints%wkp) - 1.0_dp) <= 128.0_dp*EPSILON(1.0_dp)*kpoints%nkp)
      frame%kpoints => kpoints
      frame%kp_range = kpoints%kp_range
      ALLOCATE (frame%reduced_k(3, kpoints%nkp), frame%weights(kpoints%nkp))
      frame%reduced_k(:, :) = kpoints%xkp
      frame%weights(:) = kpoints%wkp
      ALLOCATE (frame%frames(kpoints%nkp))
      DO ik = 1, kpoints%nkp
         CALL qs_upf_inverse_grid_prepare(ks_env, grid, images, frame%reduced_k(:, ik), frame%frames(ik), info, &
                                          atom_local, direct)
         IF (info /= 0) THEN
            CALL qs_upf_inverse_kpoints_release(frame)
            RETURN
         END IF
      END DO
      frame%nao = SIZE(frame%frames(1)%ao%soft, 2)
      CALL qs_upf_inverse_augmentation_prepare(ks_env, frame%frames(1)%ao, frame%augmentation)
      CALL qs_upf_inverse_gaussian_prepare(ks_env, frame%gaussian)
      info = 0
   END SUBROUTINE qs_upf_inverse_kpoints_prepare

! **************************************************************************************************
!> \brief Release owned inverse fields; the native group environment remains borrowed.
!> \param frame ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kpoints_release(frame)
      TYPE(qs_upf_inverse_kpoints_type), INTENT(INOUT)   :: frame

      IF (ALLOCATED(frame%frames)) DEALLOCATE (frame%frames)
      IF (ALLOCATED(frame%augmentation%kinds)) DEALLOCATE (frame%augmentation%kinds)
      IF (ALLOCATED(frame%reduced_k)) DEALLOCATE (frame%reduced_k)
      IF (ALLOCATED(frame%weights)) DEALLOCATE (frame%weights)
      CALL qs_upf_inverse_gaussian_release(frame%gaussian)
      NULLIFY (frame%kpoints)
      frame%nao = 0
      frame%kp_range = 0
   END SUBROUTINE qs_upf_inverse_kpoints_release

! **************************************************************************************************
!> \brief Gather native unweighted P(k) with occupations and its actual imaginary sign.
!>        FM gathering is collective within each owning group. Only its source rank
!>        contributes to the following world reduction, independently of BLACS layout.
!> \param frame ...
!> \param spin Native spin index.
!> \param density New replicated AO matrices, indexed by global k; no k weights applied.
!> \param energy_weighted Gather the prepared native W(k), without using the accepted-P cache.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kpoints_density(frame, spin, density, energy_weighted)
      TYPE(qs_upf_inverse_kpoints_type), INTENT(IN)      :: frame
      INTEGER, INTENT(IN)                                :: spin
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: density(:, :, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: energy_weighted

      INTEGER                                            :: ik, local
      LOGICAL                                            :: use_energy_weighted
      REAL(dp)                                           :: defect, tolerance
      REAL(dp), ALLOCATABLE                              :: imaginary(:, :), real_part(:, :)
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL check_live(frame)
      use_energy_weighted = .FALSE.
      IF (PRESENT(energy_weighted)) use_energy_weighted = energy_weighted
      ALLOCATE (density(frame%nao, frame%nao, SIZE(frame%weights)), &
                real_part(frame%nao, frame%nao), imaginary(frame%nao, frame%nao))
      density = CMPLX(0.0_dp, 0.0_dp, dp)
      DO local = 1, SIZE(frame%kpoints%kp_env)
         kp => frame%kpoints%kp_env(local)%kpoint_env
         ik = frame%kp_range(1) + local - 1
         CPASSERT(kp%nkpoint == ik)
         CPASSERT(spin >= 1 .AND. spin <= SIZE(kp%pmat, 2))
         CPASSERT(kp%pmat(1, spin)%matrix_struct%nrow_global == frame%nao)
         CPASSERT(kp%pmat(1, spin)%matrix_struct%ncol_global == frame%nao)
         IF (use_energy_weighted) THEN
            CPASSERT(ASSOCIATED(kp%wmat))
            CPASSERT(ALL(SHAPE(kp%wmat) == SHAPE(kp%pmat)))
            CALL cp_fm_get_submatrix(kp%wmat(1, spin), real_part)
            imaginary = 0.0_dp
            IF (.NOT. frame%kpoints%use_real_wfn) CALL cp_fm_get_submatrix(kp%wmat(2, spin), imaginary)
            IF (frame%kpoints%para_env_kp%is_source()) density(:, :, ik) = CMPLX(real_part, imaginary, dp)
         ELSE IF (ALLOCATED(frame%kpoints%native_density)) THEN
            CALL cp_cfm_get_submatrix(frame%kpoints%native_density(local, spin), density(:, :, ik))
            IF (.NOT. frame%kpoints%para_env_kp%is_source()) density(:, :, ik) = CMPLX(0.0_dp, 0.0_dp, dp)
         ELSE
            CALL cp_fm_get_submatrix(kp%pmat(1, spin), real_part)
            imaginary = 0.0_dp
            IF (.NOT. frame%kpoints%use_real_wfn) CALL cp_fm_get_submatrix(kp%pmat(2, spin), imaginary)
            IF (frame%kpoints%para_env_kp%is_source()) density(:, :, ik) = CMPLX(real_part, imaginary, dp)
         END IF
      END DO
      DO ik = 1, SIZE(frame%weights)
         CALL frame%kpoints%para_env%sum(density(:, :, ik))
         defect = MAXVAL(ABS(density(:, :, ik) - CONJG(TRANSPOSE(density(:, :, ik)))))
         tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(density(:, :, ik))))
         CPASSERT(defect <= tolerance)
      END DO
   END SUBROUTINE qs_upf_inverse_kpoints_density

! **************************************************************************************************
!> \brief Sum native k-point densities with their mesh weights exactly once.
!>        Original Fourier Q is shared and periodic; Bloch dependence remains in W(k).
!>        Tau here is smooth only and does not establish general meta-GGA support.
!> \param frame ...
!> \param spin ...
!> \param density Overwritten native total spin density.
!> \param gradient Overwritten physical Cartesian gradient.
!> \param tau Overwritten smooth kinetic density.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kpoints_collocate(frame, spin, density, gradient, tau)
      TYPE(qs_upf_inverse_kpoints_type), INTENT(INOUT)   :: frame
      INTEGER, INTENT(IN)                                :: spin
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: density, gradient(3), tau

      COMPLEX(dp), ALLOCATABLE                           :: matrices(:, :, :)
      INTEGER                                            :: axis, ik
      TYPE(pw_r3d_rs_type)                               :: dpart(3), part, tpart

      CALL qs_upf_inverse_kpoints_density(frame, spin, matrices)
      CALL part%create(density%pw_grid)
      CALL tpart%create(density%pw_grid)
      density%array = 0.0_dp
      tau%array = 0.0_dp
      DO axis = 1, 3
         CALL dpart(axis)%create(density%pw_grid)
         gradient(axis)%array = 0.0_dp
      END DO
      DO ik = 1, SIZE(frame%weights)
         matrices(:, :, ik) = frame%weights(ik)*matrices(:, :, ik)
         CALL qs_upf_inverse_grid_collocate_dense(frame%frames(ik), matrices(:, :, ik), part, dpart, tpart, &
                                                  augmentation=frame%augmentation, fourier_charge=.TRUE.)
         density%array = density%array + part%array
         tau%array = tau%array + tpart%array
         DO axis = 1, 3
            gradient(axis)%array = gradient(axis)%array + dpart(axis)%array
         END DO
      END DO
      CALL part%release()
      CALL tpart%release()
      DO axis = 1, 3
         CALL dpart(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_kpoints_collocate

! **************************************************************************************************
!> \brief Integrate the matching inverse density adjoint into native complex group matrices.
!> \param frame ...
!> \param vdensity Unweighted potential conjugate to density.
!> \param vgradient Unweighted potentials conjugate to grad(n), not sigma.
!> \param vtau Potential conjugate to smooth tau only.
!> \param operators New owned CFM array, indexed by local k; release after use.
!> \param dense Optional replicated operators indexed by global k, without mesh weights.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kpoints_integrate(frame, vdensity, vgradient, vtau, operators, dense)
      TYPE(qs_upf_inverse_kpoints_type), INTENT(INOUT)   :: frame
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3), vtau
      TYPE(cp_cfm_type), ALLOCATABLE, INTENT(OUT)        :: operators(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: dense(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: all_k(:, :, :), matrix(:, :)
      INTEGER                                            :: ik

      CALL check_live(frame)
      ALLOCATE (all_k(frame%nao, frame%nao, SIZE(frame%weights)))
      DO ik = 1, SIZE(frame%weights)
         CALL qs_upf_inverse_grid_integrate_dense(frame%frames(ik), vdensity, vgradient, vtau, matrix, &
                                                  augmentation=frame%augmentation, fourier_charge=.TRUE.)
         all_k(:, :, ik) = matrix
      END DO
      CALL qs_upf_inverse_kpoints_distribute(frame, all_k, operators)
      IF (PRESENT(dense)) CALL MOVE_ALLOC(all_k, dense)
   END SUBROUTINE qs_upf_inverse_kpoints_integrate

! **************************************************************************************************
!> \brief Build each inverse metric with an explicit choice of original projector metric.
!> \param frame ...
!> \param operators New owned local CFM metrics.
!> \param dense Optional replicated global metrics.
!> \param supplied_q Select supplied PP_Q; default uses the original Fourier-Q radial integral.
!>        A supplied-q defect is retained and never fitted or rescaled.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kpoints_metric(frame, operators, dense, supplied_q)
      TYPE(qs_upf_inverse_kpoints_type), INTENT(IN)      :: frame
      TYPE(cp_cfm_type), ALLOCATABLE, INTENT(OUT)        :: operators(:)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: dense(:, :, :)
      LOGICAL, INTENT(IN), OPTIONAL                      :: supplied_q

      COMPLEX(dp), ALLOCATABLE                           :: all_k(:, :, :), matrix(:, :)
      INTEGER                                            :: ik
      LOGICAL                                            :: use_supplied

      CALL check_live(frame)
      use_supplied = .FALSE.
      IF (PRESENT(supplied_q)) use_supplied = supplied_q
      ALLOCATE (all_k(frame%nao, frame%nao, SIZE(frame%weights)))
      DO ik = 1, SIZE(frame%weights)
         IF (use_supplied) THEN
            CALL qs_upf_inverse_grid_metric_dense(frame%frames(ik), frame%frames(ik)%ao%q, matrix)
         ELSE
            CALL qs_upf_inverse_grid_metric_dense(frame%frames(ik), frame%frames(ik)%ao%radial_q, matrix)
         END IF
         all_k(:, :, ik) = matrix
      END DO
      CALL qs_upf_inverse_kpoints_distribute(frame, all_k, operators)
      IF (PRESENT(dense)) CALL MOVE_ALLOC(all_k, dense)
   END SUBROUTINE qs_upf_inverse_kpoints_metric

! **************************************************************************************************
!> \brief Deliver full complex Bloch operators directly in the owning group's BLACS layout.
!>        No phase sum over truncated real-space images or imaginary-part projection occurs.
!> \param frame ...
!> \param dense Replicated global-k operators, unchanged.
!> \param operators New owned local group matrices.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kpoints_distribute(frame, dense, operators)
      TYPE(qs_upf_inverse_kpoints_type), INTENT(IN)      :: frame
      COMPLEX(dp), INTENT(IN)                            :: dense(:, :, :)
      TYPE(cp_cfm_type), ALLOCATABLE, INTENT(OUT)        :: operators(:)

      INTEGER                                            :: col, ik, local, row
      TYPE(cp_fm_struct_type), POINTER                   :: layout

      CALL check_live(frame)
      CPASSERT(ALL(SHAPE(dense) == [frame%nao, frame%nao, SIZE(frame%weights)]))
      ALLOCATE (operators(SIZE(frame%kpoints%kp_env)))
      DO local = 1, SIZE(operators)
         ik = frame%kp_range(1) + local - 1
         layout => frame%kpoints%kp_env(local)%kpoint_env%pmat(1, 1)%matrix_struct
         CPASSERT(ASSOCIATED(layout%para_env, frame%kpoints%para_env_kp))
         CALL cp_cfm_create(operators(local), layout)
         operators(local)%local_data = CMPLX(0.0_dp, 0.0_dp, dp)
         DO col = 1, SIZE(layout%col_indices)
            DO row = 1, SIZE(layout%row_indices)
               operators(local)%local_data(row, col) = dense(layout%row_indices(row), layout%col_indices(col), ik)
            END DO
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_kpoints_distribute

! **************************************************************************************************
!> \brief Release complex group outputs, including their retained FM structures.
!> \param operators ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_kpoints_release_operators(operators)
      TYPE(cp_cfm_type), ALLOCATABLE, INTENT(INOUT)      :: operators(:)

      INTEGER                                            :: local

      IF (.NOT. ALLOCATED(operators)) RETURN
      DO local = 1, SIZE(operators)
         CALL cp_cfm_release(operators(local))
      END DO
      DEALLOCATE (operators)
   END SUBROUTINE qs_upf_inverse_kpoints_release_operators

! **************************************************************************************************
!> \brief Verify the borrowed native group distribution and fixed Bloch mesh.
!> \param frame ...
! **************************************************************************************************
   SUBROUTINE check_live(frame)
      TYPE(qs_upf_inverse_kpoints_type), INTENT(IN)      :: frame

      INTEGER                                            :: nlocal

      CPASSERT(ASSOCIATED(frame%kpoints))
      CPASSERT(ALLOCATED(frame%frames))
      CPASSERT(frame%kpoints%nkp == SIZE(frame%weights))
      CPASSERT(ALL(frame%kpoints%kp_range == frame%kp_range))
      CPASSERT(ALL(frame%kpoints%xkp == frame%reduced_k))
      CPASSERT(ALL(frame%kpoints%wkp == frame%weights))
      nlocal = MAX(0, frame%kp_range(2) - frame%kp_range(1) + 1)
      CPASSERT(SIZE(frame%kpoints%kp_env) == nlocal)
   END SUBROUTINE check_live
END MODULE qs_upf_inverse_kpoints
