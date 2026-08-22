!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
PROGRAM kpoint_lattice_fft_unittest
   USE fft_tools,                       ONLY: finalize_fft,&
                                              init_fft
   USE kinds,                           ONLY: dp
   USE kpoint_lattice_fft,              ONLY: cell_to_k_grid_fft,&
                                              regular_kpoint_grid
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init

   IMPLICIT NONE

   INTEGER, PARAMETER                                  :: ncell = 31, nkp = 64, nselected = 3
   COMPLEX(KIND=dp), DIMENSION(2, 3, nkp)              :: direct, transformed
   COMPLEX(KIND=dp), DIMENSION(2, 3, nselected)         :: selected_transformed
   INTEGER                                             :: d, icell, ik, ix, iy, iz
   INTEGER, DIMENSION(3)                               :: nkp_grid
   INTEGER, DIMENSION(3, ncell)                        :: index_to_cell
   INTEGER, DIMENSION(nselected)                       :: selected_kpoints
   LOGICAL                                             :: used_fft
   REAL(KIND=dp), DIMENSION(3, 3)                      :: hmat
   REAL(KIND=dp), DIMENSION(3, nkp)                    :: xkp, xkp_bad
   REAL(KIND=dp), DIMENSION(2, 3, ncell)               :: values_rs
   TYPE(mp_comm_type)                                  :: mp_comm

   CALL mp_world_init(mp_comm)
   CALL init_fft("FFTSG", alltoall=.FALSE., fftsg_sizes=.TRUE., pool_limit=4, &
                 wisdom_file="", plan_style=1)

   nkp_grid = [4, 4, 4]
   hmat = RESHAPE([4.0_dp, 0.1_dp, 0.2_dp, &
                   0.0_dp, 5.0_dp, 0.3_dp, &
                   0.0_dp, 0.0_dp, 6.0_dp], [3, 3])

   ik = 0
   DO iz = 0, nkp_grid(3) - 1
      DO iy = 0, nkp_grid(2) - 1
         DO ix = 0, nkp_grid(1) - 1
            ik = ik + 1
            xkp(:, ik) = [(REAL(ix, dp) + 0.5_dp)/REAL(nkp_grid(1), dp), &
                          REAL(iy, dp)/REAL(nkp_grid(2), dp), &
                          (REAL(iz, dp) + 0.25_dp)/REAL(nkp_grid(3), dp)]
         END DO
      END DO
   END DO
   IF (.NOT. regular_kpoint_grid(xkp, nkp_grid)) ERROR STOP "Regular k-point grid not recognized"

   DO icell = 1, ncell
      index_to_cell(:, icell) = [MODULO(3*icell + 1, 11) - 5, &
                                 MODULO(5*icell + 2, 9) - 4, &
                                 MODULO(7*icell + 3, 7) - 3]
      DO d = 1, 3
         values_rs(1, d, icell) = SIN(0.13_dp*REAL(icell + 2*d, dp))
         values_rs(2, d, icell) = COS(0.17_dp*REAL(2*icell - d, dp))
      END DO
   END DO

   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, transformed, used_fft)
   IF (.NOT. used_fft) ERROR STOP "Complete regular grid did not use lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp, direct)
   IF (MAXVAL(ABS(transformed - direct)) > 2.0E-12_dp) THEN
      ERROR STOP "Shifted lattice FFT differs from direct transform"
   END IF

   selected_kpoints = [2, 31, 64]
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, selected_transformed, &
                           used_fft, selected_kpoints=selected_kpoints)
   IF (.NOT. used_fft) ERROR STOP "Selected distributed output did not use lattice FFT"
   IF (MAXVAL(ABS(selected_transformed - direct(:, :, selected_kpoints))) > 2.0E-12_dp) THEN
      ERROR STOP "Selected lattice-FFT output differs from the full transform"
   END IF

   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp, nkp_grid, transformed, used_fft, 2, hmat)
   IF (.NOT. used_fft) ERROR STOP "Derivative did not use lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp, direct, 2, hmat)
   IF (MAXVAL(ABS(transformed - direct)) > 5.0E-12_dp) THEN
      ERROR STOP "Lattice-FFT derivative differs from direct transform"
   END IF

   nkp_grid = [4, 4, 1]
   ik = 0
   DO iy = 0, nkp_grid(2) - 1
      DO ix = 0, nkp_grid(1) - 1
         ik = ik + 1
         xkp(:, ik) = [(REAL(ix, dp) + 0.5_dp)/REAL(nkp_grid(1), dp), &
                       REAL(iy, dp)/REAL(nkp_grid(2), dp), 0.0_dp]
      END DO
   END DO
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp(:, 1:ik), nkp_grid, &
                           transformed(:, :, 1:ik), used_fft)
   IF (.NOT. used_fft) ERROR STOP "Singleton grid direction did not use padded lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp(:, 1:ik), direct(:, :, 1:ik))
   IF (MAXVAL(ABS(transformed(:, :, 1:ik) - direct(:, :, 1:ik))) > 2.0E-12_dp) THEN
      ERROR STOP "Padded lattice FFT differs from direct transform"
   END IF

   nkp_grid = [4, 4, 4]

   xkp_bad = xkp
   xkp_bad(1, nkp) = xkp_bad(1, nkp) + 0.03125_dp
   CALL cell_to_k_grid_fft(values_rs, index_to_cell, xkp_bad, nkp_grid, transformed, used_fft)
   IF (used_fft) ERROR STOP "Irregular k-point list incorrectly used lattice FFT"
   CALL direct_reference(values_rs, index_to_cell, xkp_bad, direct)
   IF (MAXVAL(ABS(transformed - direct)) > 1.0E-14_dp) THEN
      ERROR STOP "Direct lattice-transform fallback is inconsistent"
   END IF

   CALL finalize_fft(mp_comm, "")
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param values_rs ...
!> \param index_to_cell ...
!> \param xkp ...
!> \param values_k ...
!> \param deriv_direction ...
!> \param hmat ...
! **************************************************************************************************
   SUBROUTINE direct_reference(values_rs, index_to_cell, xkp, values_k, deriv_direction, hmat)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: values_rs
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: xkp
      COMPLEX(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)  :: values_k
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat

      COMPLEX(KIND=dp)                                   :: factor
      INTEGER                                            :: icell, ik
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector

      values_k = z_zero
      DO ik = 1, SIZE(xkp, 2)
         DO icell = 1, SIZE(values_rs, 3)
            factor = EXP(gaussi*twopi*SUM(xkp(:, ik)*REAL(index_to_cell(:, icell), KIND=dp)))
            IF (PRESENT(deriv_direction)) THEN
               cell_vector = MATMUL(hmat, REAL(index_to_cell(:, icell), KIND=dp))
               factor = factor*gaussi*cell_vector(deriv_direction)
            END IF
            values_k(:, :, ik) = values_k(:, :, ik) + factor*values_rs(:, :, icell)
         END DO
      END DO

   END SUBROUTINE direct_reference

END PROGRAM kpoint_lattice_fft_unittest
