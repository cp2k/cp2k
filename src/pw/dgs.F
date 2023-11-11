!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH (15-Mar-2001) : Update small grid when cell changes
!>                          with dg_grid_change
! **************************************************************************************************
MODULE dgs

   USE fft_tools,                       ONLY: BWFFT,&
                                              FFT_RADIX_CLOSEST,&
                                              FFT_RADIX_NEXT_ODD,&
                                              fft3d,&
                                              fft_radix_operations
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE mathlib,                         ONLY: det_3x3,&
                                              inv_3x3
   USE pw_grid_info,                    ONLY: pw_find_cutoff
   USE pw_grid_types,                   ONLY: HALFSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_change,&
                                              pw_grid_setup
   USE pw_methods,                      ONLY: pw_copy,&
                                              pw_multiply_with,&
                                              pw_zero
   USE pw_types,                        ONLY: COMPLEXDATA3D,&
                                              pw_type
   USE realspace_grid_types,            ONLY: realspace_grid_type
   USE structure_factors,               ONLY: structure_factor_evaluate
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dgs'

   PUBLIC :: dg_get_patch
   PUBLIC :: dg_pme_grid_setup, &
             dg_sum_patch, dg_sum_patch_force_3d, dg_sum_patch_force_1d, &
             dg_get_strucfac, dg_grid_change

   INTERFACE dg_sum_patch
      MODULE PROCEDURE dg_sum_patch_coef, dg_sum_patch_arr
   END INTERFACE

   INTERFACE dg_sum_patch_force_3d
      MODULE PROCEDURE dg_sum_patch_force_coef_3d, dg_sum_patch_force_arr_3d
   END INTERFACE

   INTERFACE dg_sum_patch_force_1d
      MODULE PROCEDURE dg_sum_patch_force_coef_1d, dg_sum_patch_force_arr_1d
   END INTERFACE

   INTERFACE dg_get_patch
      MODULE PROCEDURE dg_get_patch_1, dg_get_patch_2
   END INTERFACE

   INTERFACE dg_add_patch
      MODULE PROCEDURE dg_add_patch_simple, dg_add_patch_folded
   END INTERFACE

   INTERFACE dg_int_patch_3d
      MODULE PROCEDURE dg_int_patch_simple_3d, dg_int_patch_folded_3d
   END INTERFACE

   INTERFACE dg_int_patch_1d
      MODULE PROCEDURE dg_int_patch_simple_1d, dg_int_patch_folded_1d
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param b_cell_hmat ...
!> \param npts_s ...
!> \param cutoff_radius ...
!> \param grid_s ...
!> \param grid_b ...
!> \param grid_ref ...
!> \param rs_dims ...
!> \param iounit ...
!> \param fft_usage ...
! **************************************************************************************************
   SUBROUTINE dg_pme_grid_setup(b_cell_hmat, npts_s, cutoff_radius, grid_s, grid_b, &
                                grid_ref, rs_dims, iounit, fft_usage)

      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: b_cell_hmat
      INTEGER, DIMENSION(:), INTENT(IN)                  :: npts_s
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_radius
      TYPE(pw_grid_type), POINTER                        :: grid_s, grid_b
      TYPE(pw_grid_type), INTENT(IN), OPTIONAL           :: grid_ref
      INTEGER, DIMENSION(2), INTENT(in), OPTIONAL        :: rs_dims
      INTEGER, INTENT(IN), OPTIONAL                      :: iounit
      LOGICAL, INTENT(IN), OPTIONAL                      :: fft_usage

      INTEGER, DIMENSION(2, 3)                           :: bo
      REAL(KIND=dp)                                      :: cutoff, ecut
      REAL(KIND=dp), DIMENSION(3, 3)                     :: s_cell_hmat, unit_cell_hmat

      CALL dg_find_cutoff(b_cell_hmat, npts_s, cutoff_radius, grid_s, grid_b, cutoff)

      ecut = 0.5_dp*cutoff*cutoff
      bo = grid_b%bounds
      IF (PRESENT(grid_ref)) THEN
         CALL pw_grid_setup(b_cell_hmat, grid_b, bounds=bo, cutoff=ecut, spherical=.TRUE., &
                            ref_grid=grid_ref, rs_dims=rs_dims, iounit=iounit, fft_usage=fft_usage)
      ELSE
         CALL pw_grid_setup(b_cell_hmat, grid_b, bounds=bo, cutoff=ecut, spherical=.TRUE., &
                            rs_dims=rs_dims, iounit=iounit, fft_usage=fft_usage)
      END IF

      CALL dg_find_basis(grid_b%npts, b_cell_hmat, unit_cell_hmat)

      CALL dg_set_cell(grid_s%npts, unit_cell_hmat, s_cell_hmat)

      bo = grid_s%bounds
      CALL pw_grid_setup(s_cell_hmat, grid_s, bounds=bo, cutoff=ecut, iounit=iounit, fft_usage=fft_usage)

   END SUBROUTINE dg_pme_grid_setup

! **************************************************************************************************
!> \brief Calculate the lengths of the cell vectors a, b, and c
!> \param cell_hmat ...
!> \return ...
!> \par History 01.2014 copied from cell_types::get_cell()
!> \author Ole Schuett
! **************************************************************************************************
   PURE FUNCTION get_cell_lengths(cell_hmat) RESULT(abc)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: cell_hmat
      REAL(KIND=dp), DIMENSION(3)                        :: abc

      abc(1) = SQRT(cell_hmat(1, 1)*cell_hmat(1, 1) + &
                    cell_hmat(2, 1)*cell_hmat(2, 1) + &
                    cell_hmat(3, 1)*cell_hmat(3, 1))
      abc(2) = SQRT(cell_hmat(1, 2)*cell_hmat(1, 2) + &
                    cell_hmat(2, 2)*cell_hmat(2, 2) + &
                    cell_hmat(3, 2)*cell_hmat(3, 2))
      abc(3) = SQRT(cell_hmat(1, 3)*cell_hmat(1, 3) + &
                    cell_hmat(2, 3)*cell_hmat(2, 3) + &
                    cell_hmat(3, 3)*cell_hmat(3, 3))
   END FUNCTION get_cell_lengths

! **************************************************************************************************
!> \brief ...
!> \param b_cell_hmat ...
!> \param npts_s ...
!> \param cutoff_radius ...
!> \param grid_s ...
!> \param grid_b ...
!> \param cutoff ...
! **************************************************************************************************
   SUBROUTINE dg_find_cutoff(b_cell_hmat, npts_s, cutoff_radius, grid_s, &
                             grid_b, cutoff)

      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: b_cell_hmat
      INTEGER, DIMENSION(:), INTENT(IN)                  :: npts_s
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_radius
      TYPE(pw_grid_type), POINTER                        :: grid_s, grid_b
      REAL(KIND=dp), INTENT(OUT)                         :: cutoff

      INTEGER                                            :: nout(3)
      REAL(KIND=dp)                                      :: b_cell_deth, cell_lengths(3), dr(3)
      REAL(KIND=dp), DIMENSION(3, 3)                     :: b_cell_h_inv

      b_cell_deth = ABS(det_3x3(b_cell_hmat))
      IF (b_cell_deth < 1.0E-10_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "An invalid set of cell vectors was specified. "// &
                       "The determinant det(h) is too small")
      END IF
      b_cell_h_inv = inv_3x3(b_cell_hmat)

      CALL fft_radix_operations(npts_s(1), nout(1), &
                                operation=FFT_RADIX_NEXT_ODD)
      CALL fft_radix_operations(npts_s(1), nout(2), &
                                operation=FFT_RADIX_NEXT_ODD)
      CALL fft_radix_operations(npts_s(1), nout(3), &
                                operation=FFT_RADIX_NEXT_ODD)

      cell_lengths = get_cell_lengths(b_cell_hmat)

      CALL dg_get_spacing(nout, cutoff_radius, dr)
      CALL dg_find_radix(dr, cell_lengths, grid_b%npts)

! In-line code to set grid_b % npts = npts_s if necessary
      IF (nout(1) > grid_b%npts(1)) THEN
         grid_b%npts(1) = nout(1)
         dr(1) = cell_lengths(1)/REAL(nout(1), KIND=dp)
      END IF
      IF (nout(2) > grid_b%npts(2)) THEN
         grid_b%npts(2) = nout(2)
         dr(2) = cell_lengths(2)/REAL(nout(2), KIND=dp)
      END IF
      IF (nout(3) > grid_b%npts(3)) THEN
         grid_b%npts(3) = nout(3)
         dr(3) = cell_lengths(3)/REAL(nout(3), KIND=dp)
      END IF

! big grid bounds
      grid_b%bounds(1, :) = -grid_b%npts/2
      grid_b%bounds(2, :) = +(grid_b%npts - 1)/2
      grid_b%grid_span = HALFSPACE
! small grid bounds
      grid_s%bounds(1, :) = -nout(:)/2
      grid_s%bounds(2, :) = (+nout(:) - 1)/2
      grid_s%grid_span = HALFSPACE
      grid_s%npts = nout

      cutoff = pw_find_cutoff(grid_b%npts, b_cell_h_inv)

   END SUBROUTINE dg_find_cutoff

! **************************************************************************************************
!> \brief ...
!> \param npts ...
!> \param cutoff_radius ...
!> \param dr ...
! **************************************************************************************************
   SUBROUTINE dg_get_spacing(npts, cutoff_radius, dr)

      INTEGER, DIMENSION(:), INTENT(IN)                  :: npts
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_radius
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: dr

      dr(:) = cutoff_radius/(REAL(npts(:), KIND=dp)/2.0_dp)

   END SUBROUTINE dg_get_spacing

! **************************************************************************************************
!> \brief ...
!> \param b_cell_hmat ...
!> \param grid_b ...
!> \param grid_s ...
! **************************************************************************************************
   SUBROUTINE dg_grid_change(b_cell_hmat, grid_b, grid_s)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: b_cell_hmat
      TYPE(pw_grid_type), POINTER                        :: grid_b, grid_s

      REAL(KIND=dp), DIMENSION(3, 3)                     :: s_cell_hmat, unit_cell_hmat

      CALL dg_find_basis(grid_b%npts, b_cell_hmat, unit_cell_hmat)
      CALL dg_set_cell(grid_s%npts, unit_cell_hmat, s_cell_hmat)
      CALL pw_grid_change(s_cell_hmat, grid_s)
   END SUBROUTINE dg_grid_change

! **************************************************************************************************
!> \brief ...
!> \param dr ...
!> \param cell_lengths ...
!> \param npts ...
! **************************************************************************************************
   SUBROUTINE dg_find_radix(dr, cell_lengths, npts)

      REAL(KIND=dp), INTENT(INOUT)                       :: dr(3)
      REAL(KIND=dp), INTENT(IN)                          :: cell_lengths(3)
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: npts

      INTEGER, DIMENSION(3)                              :: nin

      nin(:) = NINT(cell_lengths(:)/dr(:))
      CALL fft_radix_operations(nin(1), npts(1), &
                                operation=FFT_RADIX_CLOSEST)
      CALL fft_radix_operations(nin(2), npts(2), &
                                operation=FFT_RADIX_CLOSEST)
      CALL fft_radix_operations(nin(3), npts(3), &
                                operation=FFT_RADIX_CLOSEST)
      dr(:) = cell_lengths(:)/REAL(npts(:), KIND=dp)

   END SUBROUTINE dg_find_radix

! **************************************************************************************************
!> \brief ...
!> \param npts ...
!> \param cell_hmat ...
!> \param unit_cell_hmat ...
! **************************************************************************************************
   PURE SUBROUTINE dg_find_basis(npts, cell_hmat, unit_cell_hmat)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: npts
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: cell_hmat
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: unit_cell_hmat

      INTEGER                                            :: i

      DO i = 1, 3
         unit_cell_hmat(:, i) = cell_hmat(:, i)/REAL(npts(:), KIND=dp)
      END DO

   END SUBROUTINE dg_find_basis

!! Calculation of the basis on the mesh 'box'

! **************************************************************************************************
!> \brief ...
!> \param npts ...
!> \param unit_cell_hmat ...
!> \param cell_hmat ...
! **************************************************************************************************
   PURE SUBROUTINE dg_set_cell(npts, unit_cell_hmat, cell_hmat)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: npts
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: unit_cell_hmat
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT)        :: cell_hmat

! computing the unit vector along a, b, c and scaling it to length dr:

      cell_hmat(:, 1) = unit_cell_hmat(:, 1)*npts(1)
      cell_hmat(:, 2) = unit_cell_hmat(:, 2)*npts(2)
      cell_hmat(:, 3) = unit_cell_hmat(:, 3)*npts(3)

   END SUBROUTINE dg_set_cell

! **************************************************************************************************
!> \brief ...
!> \param cell_hmat ...
!> \param r ...
!> \param npts_s ...
!> \param npts_b ...
!> \param centre ...
!> \param lb ...
!> \param ex ...
!> \param ey ...
!> \param ez ...
! **************************************************************************************************
   SUBROUTINE dg_get_strucfac(cell_hmat, r, npts_s, npts_b, centre, lb, ex, ey, ez)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: cell_hmat
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r
      INTEGER, DIMENSION(:), INTENT(IN)                  :: npts_s, npts_b
      INTEGER, INTENT(OUT)                               :: centre(3)
      INTEGER, INTENT(IN)                                :: lb(3)
      COMPLEX(KIND=dp), DIMENSION(lb(1):), INTENT(OUT)   :: ex
      COMPLEX(KIND=dp), DIMENSION(lb(2):), INTENT(OUT)   :: ey
      COMPLEX(KIND=dp), DIMENSION(lb(3):), INTENT(OUT)   :: ez

      REAL(KIND=dp)                                      :: delta(3)

      CALL dg_get_delta(cell_hmat, r, npts_s, npts_b, centre, delta)

      CALL structure_factor_evaluate(delta, lb, ex, ey, ez)

   END SUBROUTINE dg_get_strucfac

! **************************************************************************************************
!> \brief ...
!> \param cell_hmat ...
!> \param r ...
!> \param npts_s ...
!> \param npts_b ...
!> \param centre ...
!> \param delta ...
! **************************************************************************************************
   SUBROUTINE dg_get_delta(cell_hmat, r, npts_s, npts_b, centre, delta)
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: cell_hmat
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r
      INTEGER, DIMENSION(:), INTENT(IN)                  :: npts_s, npts_b
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: centre
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: delta

      REAL(KIND=dp)                                      :: cell_deth
      REAL(KIND=dp), DIMENSION(3)                        :: grid_i, s
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cell_h_inv

      cell_deth = ABS(det_3x3(cell_hmat))
      IF (cell_deth < 1.0E-10_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "An invalid set of cell vectors was specified. "// &
                       "The determinant det(h) is too small")
      END IF

      cell_h_inv = inv_3x3(cell_hmat)

! compute the scaled coordinate of atomi
      s = MATMUL(cell_h_inv, r)
      s = s - NINT(s)

! find the continuous ``grid'' point (on big grid)
      grid_i(1:3) = REAL(npts_b(1:3), KIND=dp)*s(1:3)

! find the closest grid point (on big grid)
      centre(:) = NINT(grid_i(:))

! find the distance vector
      delta(:) = (grid_i(:) - centre(:))/REAL(npts_s(:), KIND=dp)

      centre(:) = centre(:) + npts_b(:)/2
      centre(:) = MODULO(centre(:), npts_b(:))
      centre(:) = centre(:) - npts_b(:)/2

   END SUBROUTINE dg_get_delta

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param rhos ...
!> \param center ...
! **************************************************************************************************
   PURE SUBROUTINE dg_sum_patch_coef(rs, rhos, center)

      TYPE(realspace_grid_type), INTENT(INOUT)           :: rs
      TYPE(pw_type), INTENT(IN)                          :: rhos
      INTEGER, DIMENSION(3), INTENT(IN)                  :: center

      INTEGER                                            :: i, ia, ii
      INTEGER, DIMENSION(3)                              :: nc
      LOGICAL                                            :: folded

      folded = .FALSE.

      DO i = rhos%pw_grid%bounds(1, 1), rhos%pw_grid%bounds(2, 1)
         ia = i - rhos%pw_grid%bounds(1, 1) + 1
         ii = center(1) + i - rs%lb_local(1)
         IF (ii < 0) THEN
            rs%px(ia) = ii + rs%npts_local(1) + 1
            folded = .TRUE.
         ELSEIF (ii >= rs%npts_local(1)) THEN
            rs%px(ia) = ii - rs%npts_local(1) + 1
            folded = .TRUE.
         ELSE
            rs%px(ia) = ii + 1
         END IF
      END DO
      DO i = rhos%pw_grid%bounds(1, 2), rhos%pw_grid%bounds(2, 2)
         ia = i - rhos%pw_grid%bounds(1, 2) + 1
         ii = center(2) + i - rs%lb_local(2)
         IF (ii < 0) THEN
            rs%py(ia) = ii + rs%npts_local(2) + 1
            folded = .TRUE.
         ELSEIF (ii >= rs%npts_local(2)) THEN
            rs%py(ia) = ii - rs%npts_local(2) + 1
            folded = .TRUE.
         ELSE
            rs%py(ia) = ii + 1
         END IF
      END DO
      DO i = rhos%pw_grid%bounds(1, 3), rhos%pw_grid%bounds(2, 3)
         ia = i - rhos%pw_grid%bounds(1, 3) + 1
         ii = center(3) + i - rs%lb_local(3)
         IF (ii < 0) THEN
            rs%pz(ia) = ii + rs%npts_local(3) + 1
            folded = .TRUE.
         ELSEIF (ii >= rs%npts_local(3)) THEN
            rs%pz(ia) = ii - rs%npts_local(3) + 1
            folded = .TRUE.
         ELSE
            rs%pz(ia) = ii + 1
         END IF
      END DO

      IF (folded) THEN
         CALL dg_add_patch(rs%r, rhos%cr3d, rhos%pw_grid%npts, &
                           rs%px, rs%py, rs%pz)
      ELSE
         nc(1) = rs%px(1) - 1
         nc(2) = rs%py(1) - 1
         nc(3) = rs%pz(1) - 1
         CALL dg_add_patch(rs%r, rhos%cr3d, rhos%pw_grid%npts, nc)
      END IF

   END SUBROUTINE dg_sum_patch_coef

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param rhos ...
!> \param center ...
! **************************************************************************************************
   PURE SUBROUTINE dg_sum_patch_arr(rs, rhos, center)

      TYPE(realspace_grid_type), INTENT(INOUT)           :: rs
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: rhos
      INTEGER, DIMENSION(3), INTENT(IN)                  :: center

      INTEGER                                            :: i, ia, ii
      INTEGER, DIMENSION(3)                              :: lb, nc, ns, ub
      LOGICAL                                            :: folded

      ns(1) = SIZE(rhos, 1)
      ns(2) = SIZE(rhos, 2)
      ns(3) = SIZE(rhos, 3)
      lb = -(ns - 1)/2
      ub = lb + ns - 1
      folded = .FALSE.

      DO i = lb(1), ub(1)
         ia = i - lb(1) + 1
         ii = center(1) + i - rs%lb_local(1)
         IF (ii < 0) THEN
            rs%px(ia) = ii + rs%npts_local(1) + 1
            folded = .TRUE.
         ELSEIF (ii >= rs%npts_local(1)) THEN
            rs%px(ia) = ii - rs%npts_local(1) + 1
            folded = .TRUE.
         ELSE
            rs%px(ia) = ii + 1
         END IF
      END DO
      DO i = lb(2), ub(2)
         ia = i - lb(2) + 1
         ii = center(2) + i - rs%lb_local(2)
         IF (ii < 0) THEN
            rs%py(ia) = ii + rs%npts_local(2) + 1
            folded = .TRUE.
         ELSEIF (ii >= rs%npts_local(2)) THEN
            rs%py(ia) = ii - rs%npts_local(2) + 1
            folded = .TRUE.
         ELSE
            rs%py(ia) = ii + 1
         END IF
      END DO
      DO i = lb(3), ub(3)
         ia = i - lb(3) + 1
         ii = center(3) + i - rs%lb_local(3)
         IF (ii < 0) THEN
            rs%pz(ia) = ii + rs%npts_local(3) + 1
            folded = .TRUE.
         ELSEIF (ii >= rs%npts_local(3)) THEN
            rs%pz(ia) = ii - rs%npts_local(3) + 1
            folded = .TRUE.
         ELSE
            rs%pz(ia) = ii + 1
         END IF
      END DO

      IF (folded) THEN
         CALL dg_add_patch(rs%r, rhos, ns, rs%px, rs%py, rs%pz)
      ELSE
         nc(1) = rs%px(1) - 1
         nc(2) = rs%py(1) - 1
         nc(3) = rs%pz(1) - 1
         CALL dg_add_patch(rs%r, rhos, ns, nc)
      END IF

   END SUBROUTINE dg_sum_patch_arr

! **************************************************************************************************
!> \brief ...
!> \param drpot ...
!> \param rhos ...
!> \param center ...
!> \param force ...
! **************************************************************************************************
   PURE SUBROUTINE dg_sum_patch_force_arr_3d(drpot, rhos, center, force)

      TYPE(realspace_grid_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: drpot
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: rhos
      INTEGER, DIMENSION(3), INTENT(IN)                  :: center
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: force

      INTEGER                                            :: i, ia, ii
      INTEGER, DIMENSION(3)                              :: lb, nc, ns, ub
      LOGICAL                                            :: folded

      ns(1) = SIZE(rhos, 1)
      ns(2) = SIZE(rhos, 2)
      ns(3) = SIZE(rhos, 3)
      lb = -(ns - 1)/2
      ub = lb + ns - 1
      folded = .FALSE.

      DO i = lb(1), ub(1)
         ia = i - lb(1) + 1
         ii = center(1) + i - drpot(1)%lb_local(1)
         IF (ii < 0) THEN
            drpot(1)%px(ia) = ii + drpot(1)%npts_local(1) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot(1)%npts_local(1)) THEN
            drpot(1)%px(ia) = ii - drpot(1)%npts_local(1) + 1
            folded = .TRUE.
         ELSE
            drpot(1)%px(ia) = ii + 1
         END IF
      END DO
      DO i = lb(2), ub(2)
         ia = i - lb(2) + 1
         ii = center(2) + i - drpot(1)%lb_local(2)
         IF (ii < 0) THEN
            drpot(1)%py(ia) = ii + drpot(1)%npts_local(2) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot(1)%npts_local(2)) THEN
            drpot(1)%py(ia) = ii - drpot(1)%npts_local(2) + 1
            folded = .TRUE.
         ELSE
            drpot(1)%py(ia) = ii + 1
         END IF
      END DO
      DO i = lb(3), ub(3)
         ia = i - lb(3) + 1
         ii = center(3) + i - drpot(1)%lb_local(3)
         IF (ii < 0) THEN
            drpot(1)%pz(ia) = ii + drpot(1)%npts_local(3) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot(1)%npts_local(3)) THEN
            drpot(1)%pz(ia) = ii - drpot(1)%npts_local(3) + 1
            folded = .TRUE.
         ELSE
            drpot(1)%pz(ia) = ii + 1
         END IF
      END DO

      IF (folded) THEN
         CALL dg_int_patch_3d(drpot(1)%r, drpot(2)%r, &
                              drpot(3)%r, rhos, force, ns, &
                              drpot(1)%px, drpot(1)%py, drpot(1)%pz)
      ELSE
         nc(1) = drpot(1)%px(1) - 1
         nc(2) = drpot(1)%py(1) - 1
         nc(3) = drpot(1)%pz(1) - 1
         CALL dg_int_patch_3d(drpot(1)%r, drpot(2)%r, &
                              drpot(3)%r, rhos, force, ns, nc)
      END IF

   END SUBROUTINE dg_sum_patch_force_arr_3d

! **************************************************************************************************
!> \brief ...
!> \param drpot ...
!> \param rhos ...
!> \param center ...
!> \param force ...
! **************************************************************************************************
   PURE SUBROUTINE dg_sum_patch_force_arr_1d(drpot, rhos, center, force)

      TYPE(realspace_grid_type), INTENT(INOUT)           :: drpot
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: rhos
      INTEGER, DIMENSION(3), INTENT(IN)                  :: center
      REAL(KIND=dp), INTENT(OUT)                         :: force

      INTEGER                                            :: i, ia, ii
      INTEGER, DIMENSION(3)                              :: lb, nc, ns, ub
      LOGICAL                                            :: folded

      ns(1) = SIZE(rhos, 1)
      ns(2) = SIZE(rhos, 2)
      ns(3) = SIZE(rhos, 3)
      lb = -(ns - 1)/2
      ub = lb + ns - 1
      folded = .FALSE.

      DO i = lb(1), ub(1)
         ia = i - lb(1) + 1
         ii = center(1) + i - drpot%lb_local(1)
         IF (ii < 0) THEN
            drpot%px(ia) = ii + drpot%npts_local(1) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot%desc%npts(1)) THEN
            drpot%px(ia) = ii - drpot%npts_local(1) + 1
            folded = .TRUE.
         ELSE
            drpot%px(ia) = ii + 1
         END IF
      END DO
      DO i = lb(2), ub(2)
         ia = i - lb(2) + 1
         ii = center(2) + i - drpot%lb_local(2)
         IF (ii < 0) THEN
            drpot%py(ia) = ii + drpot%npts_local(2) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot%desc%npts(2)) THEN
            drpot%py(ia) = ii - drpot%npts_local(2) + 1
            folded = .TRUE.
         ELSE
            drpot%py(ia) = ii + 1
         END IF
      END DO
      DO i = lb(3), ub(3)
         ia = i - lb(3) + 1
         ii = center(3) + i - drpot%lb_local(3)
         IF (ii < 0) THEN
            drpot%pz(ia) = ii + drpot%npts_local(3) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot%desc%npts(3)) THEN
            drpot%pz(ia) = ii - drpot%npts_local(3) + 1
            folded = .TRUE.
         ELSE
            drpot%pz(ia) = ii + 1
         END IF
      END DO

      IF (folded) THEN
         CALL dg_int_patch_1d(drpot%r, rhos, force, ns, &
                              drpot%px, drpot%py, drpot%pz)
      ELSE
         nc(1) = drpot%px(1) - 1
         nc(2) = drpot%py(1) - 1
         nc(3) = drpot%pz(1) - 1
         CALL dg_int_patch_1d(drpot%r, rhos, force, ns, nc)
      END IF

   END SUBROUTINE dg_sum_patch_force_arr_1d

! **************************************************************************************************
!> \brief ...
!> \param drpot ...
!> \param rhos ...
!> \param center ...
!> \param force ...
! **************************************************************************************************
   PURE SUBROUTINE dg_sum_patch_force_coef_3d(drpot, rhos, center, force)

      TYPE(realspace_grid_type), DIMENSION(:), &
         INTENT(INOUT)                                   :: drpot
      TYPE(pw_type), INTENT(IN)                          :: rhos
      INTEGER, DIMENSION(3), INTENT(IN)                  :: center
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: force

      INTEGER                                            :: i, ia, ii
      INTEGER, DIMENSION(3)                              :: nc
      LOGICAL                                            :: folded

      folded = .FALSE.

      DO i = rhos%pw_grid%bounds(1, 1), rhos%pw_grid%bounds(2, 1)
         ia = i - rhos%pw_grid%bounds(1, 1) + 1
         ii = center(1) + i - drpot(1)%lb_local(1)
         IF (ii < 0) THEN
            drpot(1)%px(ia) = ii + drpot(1)%desc%npts(1) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot(1)%desc%npts(1)) THEN
            drpot(1)%px(ia) = ii - drpot(1)%desc%npts(1) + 1
            folded = .TRUE.
         ELSE
            drpot(1)%px(ia) = ii + 1
         END IF
      END DO
      DO i = rhos%pw_grid%bounds(1, 2), rhos%pw_grid%bounds(2, 2)
         ia = i - rhos%pw_grid%bounds(1, 2) + 1
         ii = center(2) + i - drpot(1)%lb_local(2)
         IF (ii < 0) THEN
            drpot(1)%py(ia) = ii + drpot(1)%desc%npts(2) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot(1)%desc%npts(2)) THEN
            drpot(1)%py(ia) = ii - drpot(1)%desc%npts(2) + 1
            folded = .TRUE.
         ELSE
            drpot(1)%py(ia) = ii + 1
         END IF
      END DO
      DO i = rhos%pw_grid%bounds(1, 3), rhos%pw_grid%bounds(2, 3)
         ia = i - rhos%pw_grid%bounds(1, 3) + 1
         ii = center(3) + i - drpot(1)%lb_local(3)
         IF (ii < 0) THEN
            drpot(1)%pz(ia) = ii + drpot(1)%desc%npts(3) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot(1)%desc%npts(3)) THEN
            drpot(1)%pz(ia) = ii - drpot(1)%desc%npts(3) + 1
            folded = .TRUE.
         ELSE
            drpot(1)%pz(ia) = ii + 1
         END IF
      END DO

      IF (folded) THEN
         CALL dg_int_patch_3d(drpot(1)%r, drpot(2)%r, &
                              drpot(3)%r, rhos%cr3d, force, rhos%pw_grid%npts, &
                              drpot(1)%px, drpot(1)%py, drpot(1)%pz)
      ELSE
         nc(1) = drpot(1)%px(1) - 1
         nc(2) = drpot(1)%py(1) - 1
         nc(3) = drpot(1)%pz(1) - 1
         CALL dg_int_patch_3d(drpot(1)%r, drpot(2)%r, &
                              drpot(3)%r, rhos%cr3d, force, rhos%pw_grid%npts, nc)
      END IF

   END SUBROUTINE dg_sum_patch_force_coef_3d

! **************************************************************************************************
!> \brief ...
!> \param drpot ...
!> \param rhos ...
!> \param center ...
!> \param force ...
! **************************************************************************************************
   PURE SUBROUTINE dg_sum_patch_force_coef_1d(drpot, rhos, center, force)

      TYPE(realspace_grid_type), INTENT(INOUT)           :: drpot
      TYPE(pw_type), INTENT(IN)                          :: rhos
      INTEGER, DIMENSION(3), INTENT(IN)                  :: center
      REAL(KIND=dp), INTENT(OUT)                         :: force

      INTEGER                                            :: i, ia, ii
      INTEGER, DIMENSION(3)                              :: nc
      LOGICAL                                            :: folded

      folded = .FALSE.

      DO i = rhos%pw_grid%bounds(1, 1), rhos%pw_grid%bounds(2, 1)
         ia = i - rhos%pw_grid%bounds(1, 1) + 1
         ii = center(1) + i - drpot%lb_local(1)
         IF (ii < 0) THEN
            drpot%px(ia) = ii + drpot%desc%npts(1) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot%desc%npts(1)) THEN
            drpot%px(ia) = ii - drpot%desc%npts(1) + 1
            folded = .TRUE.
         ELSE
            drpot%px(ia) = ii + 1
         END IF
      END DO
      DO i = rhos%pw_grid%bounds(1, 2), rhos%pw_grid%bounds(2, 2)
         ia = i - rhos%pw_grid%bounds(1, 2) + 1
         ii = center(2) + i - drpot%lb_local(2)
         IF (ii < 0) THEN
            drpot%py(ia) = ii + drpot%desc%npts(2) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot%desc%npts(2)) THEN
            drpot%py(ia) = ii - drpot%desc%npts(2) + 1
            folded = .TRUE.
         ELSE
            drpot%py(ia) = ii + 1
         END IF
      END DO
      DO i = rhos%pw_grid%bounds(1, 3), rhos%pw_grid%bounds(2, 3)
         ia = i - rhos%pw_grid%bounds(1, 3) + 1
         ii = center(3) + i - drpot%lb_local(3)
         IF (ii < 0) THEN
            drpot%pz(ia) = ii + drpot%desc%npts(3) + 1
            folded = .TRUE.
         ELSEIF (ii >= drpot%desc%npts(3)) THEN
            drpot%pz(ia) = ii - drpot%desc%npts(3) + 1
            folded = .TRUE.
         ELSE
            drpot%pz(ia) = ii + 1
         END IF
      END DO

      IF (folded) THEN
         CALL dg_int_patch_1d(drpot%r, rhos%cr3d, force, &
                              rhos%pw_grid%npts, drpot%px, drpot%py, drpot%pz)
      ELSE
         nc(1) = drpot%px(1) - 1
         nc(2) = drpot%py(1) - 1
         nc(3) = drpot%pz(1) - 1
         CALL dg_int_patch_1d(drpot%r, rhos%cr3d, force, rhos%pw_grid%npts, nc)
      END IF

   END SUBROUTINE dg_sum_patch_force_coef_1d

! **************************************************************************************************
!> \brief ...
!> \param rho0 ...
!> \param rhos1 ...
!> \param charge1 ...
!> \param ex1 ...
!> \param ey1 ...
!> \param ez1 ...
! **************************************************************************************************
   SUBROUTINE dg_get_patch_1(rho0, rhos1, charge1, ex1, ey1, ez1)

      TYPE(pw_type), INTENT(IN)                          :: rho0
      TYPE(pw_type), INTENT(INOUT)                       :: rhos1
      REAL(KIND=dp), INTENT(IN)                          :: charge1
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: ex1, ey1, ez1

      COMPLEX(KIND=dp)                                   :: za, zb
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: zs
      INTEGER                                            :: nd(3)
      TYPE(pw_type)                                      :: cd

      nd = rhos1%pw_grid%npts

      ALLOCATE (zs(nd(1)*nd(2)))
      zs = 0.0_dp
      CALL cd%create(rho0%pw_grid, COMPLEXDATA3D, rho0%in_space)
      CALL pw_zero(cd)

      za = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      zb = CMPLX(charge1, 0.0_dp, KIND=dp)
      CALL rankup(nd, za, cd%cc3d, zb, ex1, ey1, ez1, zs)
      CALL pw_multiply_with(cd, rho0)
      CALL fft3d(BWFFT, nd, cd%cc3d)
      CALL pw_copy(cd, rhos1)

      DEALLOCATE (zs)
      CALL cd%release()

   END SUBROUTINE dg_get_patch_1

! **************************************************************************************************
!> \brief ...
!> \param rho0 ...
!> \param rhos1 ...
!> \param rhos2 ...
!> \param charge1 ...
!> \param charge2 ...
!> \param ex1 ...
!> \param ey1 ...
!> \param ez1 ...
!> \param ex2 ...
!> \param ey2 ...
!> \param ez2 ...
! **************************************************************************************************
   SUBROUTINE dg_get_patch_2(rho0, rhos1, rhos2, charge1, charge2, &
                             ex1, ey1, ez1, ex2, ey2, ez2)

      TYPE(pw_type), INTENT(IN)                          :: rho0
      TYPE(pw_type), INTENT(INOUT)                       :: rhos1, rhos2
      REAL(KIND=dp), INTENT(IN)                          :: charge1, charge2
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: ex1, ey1, ez1, ex2, ey2, ez2

      COMPLEX(KIND=dp)                                   :: za, zb
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: zs
      INTEGER                                            :: nd(3)
      TYPE(pw_type)                                      :: cd

      nd = rhos1%pw_grid%npts

      ALLOCATE (zs(nd(1)*nd(2)))
      zs = 0.0_dp
      CALL cd%create(rhos1%pw_grid, COMPLEXDATA3D, rho0%in_space)
      CALL pw_zero(cd)

      za = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      zb = CMPLX(charge2, 0.0_dp, KIND=dp)
      CALL rankup(nd, za, cd%cc3d, zb, ex2, ey2, ez2, zs)
      za = CMPLX(0.0_dp, 1.0_dp, KIND=dp)
      zb = CMPLX(charge1, 0.0_dp, KIND=dp)
      CALL rankup(nd, za, cd%cc3d, zb, ex1, ey1, ez1, zs)
      CALL pw_multiply_with(cd, rho0)
      CALL fft3d(BWFFT, nd, cd%cc3d)
      CALL copy_cri(cd%cc3d, rhos1%cr3d, rhos2%cr3d)

      DEALLOCATE (zs)
      CALL cd%release()

   END SUBROUTINE dg_get_patch_2

! **************************************************************************************************
!> \brief ...
!> \param rb ...
!> \param rs ...
!> \param ns ...
!> \param nc ...
! **************************************************************************************************
   PURE SUBROUTINE dg_add_patch_simple(rb, rs, ns, nc)

      REAL(KIND=dp), INTENT(INOUT)                       :: rb(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: rs(:, :, :)
      INTEGER, INTENT(IN)                                :: ns(3), nc(3)

      INTEGER                                            :: i, ii, j, jj, k, kk

      DO k = 1, ns(3)
         kk = nc(3) + k
         DO j = 1, ns(2)
            jj = nc(2) + j
            DO i = 1, ns(1)
               ii = nc(1) + i
               rb(ii, jj, kk) = rb(ii, jj, kk) + rs(i, j, k)
            END DO
         END DO
      END DO

   END SUBROUTINE dg_add_patch_simple

! **************************************************************************************************
!> \brief ...
!> \param rb ...
!> \param rs ...
!> \param ns ...
!> \param px ...
!> \param py ...
!> \param pz ...
! **************************************************************************************************
   PURE SUBROUTINE dg_add_patch_folded(rb, rs, ns, px, py, pz)

      REAL(KIND=dp), INTENT(INOUT)                       :: rb(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: rs(:, :, :)
      INTEGER, INTENT(IN)                                :: ns(:)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: px, py, pz

      INTEGER                                            :: i, ii, j, jj, k, kk

      DO k = 1, ns(3)
         kk = pz(k)
         DO j = 1, ns(2)
            jj = py(j)
            DO i = 1, ns(1)
               ii = px(i)
               rb(ii, jj, kk) = rb(ii, jj, kk) + rs(i, j, k)
            END DO
         END DO
      END DO

   END SUBROUTINE dg_add_patch_folded

! **************************************************************************************************
!> \brief ...
!> \param rbx ...
!> \param rby ...
!> \param rbz ...
!> \param rs ...
!> \param f ...
!> \param ns ...
!> \param nc ...
! **************************************************************************************************
   PURE SUBROUTINE dg_int_patch_simple_3d(rbx, rby, rbz, rs, f, ns, nc)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: rbx, rby, rbz, rs
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: f
      INTEGER, INTENT(IN)                                :: ns(3), nc(3)

      INTEGER                                            :: i, ii, j, jj, k, kk
      REAL(KIND=dp)                                      :: s

      f = 0.0_dp
      DO k = 1, ns(3)
         kk = nc(3) + k
         DO j = 1, ns(2)
            jj = nc(2) + j
            DO i = 1, ns(1)
               ii = nc(1) + i
               s = rs(i, j, k)
               f(1) = f(1) + s*rbx(ii, jj, kk)
               f(2) = f(2) + s*rby(ii, jj, kk)
               f(3) = f(3) + s*rbz(ii, jj, kk)
            END DO
         END DO
      END DO

   END SUBROUTINE dg_int_patch_simple_3d

! **************************************************************************************************
!> \brief ...
!> \param rb ...
!> \param rs ...
!> \param f ...
!> \param ns ...
!> \param nc ...
! **************************************************************************************************
   PURE SUBROUTINE dg_int_patch_simple_1d(rb, rs, f, ns, nc)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: rb, rs
      REAL(KIND=dp), INTENT(OUT)                         :: f
      INTEGER, INTENT(IN)                                :: ns(3), nc(3)

      INTEGER                                            :: i, ii, j, jj, k, kk
      REAL(KIND=dp)                                      :: s

      f = 0.0_dp
      DO k = 1, ns(3)
         kk = nc(3) + k
         DO j = 1, ns(2)
            jj = nc(2) + j
            DO i = 1, ns(1)
               ii = nc(1) + i
               s = rs(i, j, k)
               f = f + s*rb(ii, jj, kk)
            END DO
         END DO
      END DO

   END SUBROUTINE dg_int_patch_simple_1d

! **************************************************************************************************
!> \brief ...
!> \param rbx ...
!> \param rby ...
!> \param rbz ...
!> \param rs ...
!> \param f ...
!> \param ns ...
!> \param px ...
!> \param py ...
!> \param pz ...
! **************************************************************************************************
   PURE SUBROUTINE dg_int_patch_folded_3d(rbx, rby, rbz, rs, f, ns, px, py, pz)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: rbx, rby, rbz, rs
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT)         :: f
      INTEGER, INTENT(IN)                                :: ns(3)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: px, py, pz

      INTEGER                                            :: i, ii, j, jj, k, kk
      REAL(KIND=dp)                                      :: s

      f = 0.0_dp
      DO k = 1, ns(3)
         kk = pz(k)
         DO j = 1, ns(2)
            jj = py(j)
            DO i = 1, ns(1)
               ii = px(i)
               s = rs(i, j, k)
               f(1) = f(1) + s*rbx(ii, jj, kk)
               f(2) = f(2) + s*rby(ii, jj, kk)
               f(3) = f(3) + s*rbz(ii, jj, kk)
            END DO
         END DO
      END DO

   END SUBROUTINE dg_int_patch_folded_3d

! **************************************************************************************************
!> \brief ...
!> \param rb ...
!> \param rs ...
!> \param f ...
!> \param ns ...
!> \param px ...
!> \param py ...
!> \param pz ...
! **************************************************************************************************
   PURE SUBROUTINE dg_int_patch_folded_1d(rb, rs, f, ns, px, py, pz)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: rb, rs
      REAL(KIND=dp), INTENT(INOUT)                       :: f
      INTEGER, INTENT(IN)                                :: ns(3)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: px, py, pz

      INTEGER                                            :: i, ii, j, jj, k, kk
      REAL(KIND=dp)                                      :: s

      f = 0.0_dp
      DO k = 1, ns(3)
         kk = pz(k)
         DO j = 1, ns(2)
            jj = py(j)
            DO i = 1, ns(1)
               ii = px(i)
               s = rs(i, j, k)
               f = f + s*rb(ii, jj, kk)
            END DO
         END DO
      END DO

   END SUBROUTINE dg_int_patch_folded_1d

! **************************************************************************************************
!> \brief ...
!> \param n ...
!> \param za ...
!> \param cmat ...
!> \param zb ...
!> \param ex ...
!> \param ey ...
!> \param ez ...
!> \param scr ...
! **************************************************************************************************
   SUBROUTINE rankup(n, za, cmat, zb, ex, ey, ez, scr)
!
! cmat(i,j,k) <- za * cmat(i,j,k) + ex(i) * ey(j) * ez(k)
!

      INTEGER, DIMENSION(3), INTENT(IN)                  :: n
      COMPLEX(KIND=dp), INTENT(IN)                       :: za
      COMPLEX(KIND=dp), DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: cmat
      COMPLEX(KIND=dp), INTENT(IN)                       :: zb
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(IN)         :: ex, ey, ez
      COMPLEX(KIND=dp), DIMENSION(:), INTENT(INOUT)      :: scr

      INTEGER                                            :: n2, n3

      n2 = n(1)*n(2)
      n3 = n2*n(3)
      scr(1:n2) = z_zero
      CALL zgeru(n(1), n(2), zb, ex, 1, ey, 1, scr, n(1))
      CALL zscal(n3, za, cmat, 1)
      CALL zgeru(n2, n(3), z_one, scr, 1, ez, 1, cmat, n2)

   END SUBROUTINE rankup

! **************************************************************************************************
!> \brief Copy a the real and imag. parts of a complex 3D array into two real arrays
!> \param z the complex array
!> \param r1 the real array for the real part
!> \param r2 the real array for the imaginary part
! **************************************************************************************************
   SUBROUTINE copy_cri(z, r1, r2)
!
! r1 = real ( z )
! r2 = imag ( z )
!

      COMPLEX(KIND=dp), INTENT(IN)                       :: z(:, :, :)
      REAL(KIND=dp), INTENT(INOUT)                       :: r1(:, :, :), r2(:, :, :)

!$OMP PARALLEL WORKSHARE DEFAULT(NONE), SHARED(r1,r2,z)
      r1(:, :, :) = REAL(z(:, :, :), KIND=dp)
      r2(:, :, :) = AIMAG(z(:, :, :))
!$OMP END PARALLEL WORKSHARE

   END SUBROUTINE copy_cri

END MODULE dgs
