!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Three-centre integrals d_lP(T) of periodic RI-RS with k-points, evaluated analytically or
!>        by quadrature (keyword D_LP_COMPUTATION).
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE ri_rs_kp_compute_d_lP
   USE cell_types,                      ONLY: pbc
   USE input_constants,                 ONLY: ri_rs_d_lP_analytic,&
                                              ri_rs_d_lP_auto
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_print_mem_per_rank
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE ri_rs_kp_compute_d_lP_analytic,  ONLY: analytic_d_lP_is_faster,&
                                              compute_d_lP_R_analytic
   USE ri_rs_kp_compute_d_lP_numeric,   ONLY: compute_d_lP_R_numeric,&
                                              compute_d_lP_cutoffs
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ri_rs_kp_compute_d_lP'

   PUBLIC :: compute_d_lP_R

CONTAINS

! **************************************************************************************************
!> \brief d_lP(T) = (f_l | v | φ_P,T) on the cells T within the cutoff R_cut^(lP) and the charges
!>        Q_l at q = 0. D_LP_COMPUTATION AUTO takes the analytic integrals when their estimated
!>        time is smaller than the one of the quadrature.
!> \param bs_env Band-structure environment.
!> \param grid_points Reference-cell RI-RS points r_l.
!> \param cells_R_d_lP Integer cell vectors T of d_lP(T), indexed (3,T).
!> \param d_lP_R d_lP(T) indexed by local point l, RI function P and cell T.
!> \param qpoints Must hold the single point q = 0.
!> \param charge_lq Q_l at q = 0.
! **************************************************************************************************
   SUBROUTINE compute_d_lP_R(bs_env, grid_points, cells_R_d_lP, d_lP_R, qpoints, charge_lq)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: cells_R_d_lP(:, :)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: d_lP_R(:, :, :)
      REAL(KIND=dp), INTENT(IN)                          :: qpoints(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: charge_lq(:, :)

      LOGICAL                                            :: analytic_d_lP

      CALL build_cells_R_for_d_lP_R(bs_env, grid_points, cells_R_d_lP)
      CPASSERT(SIZE(cells_R_d_lP, 2) > 0)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,I9)') 'RI-RS d_lP^R cells R:', SIZE(cells_R_d_lP, 2)
      END IF

      SELECT CASE (bs_env%ri_rs%d_lP_computation)
      CASE (ri_rs_d_lP_auto)
         analytic_d_lP = analytic_d_lP_is_faster(bs_env, grid_points, cells_R_d_lP)
      CASE DEFAULT
         analytic_d_lP = (bs_env%ri_rs%d_lP_computation == ri_rs_d_lP_analytic)
      END SELECT
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T72,A9)') 'RI-RS d_lP computation:', MERGE('ANALYTIC', ' NUMERIC', analytic_d_lP)
         FLUSH (bs_env%unit_nr)
      END IF
      IF (analytic_d_lP) THEN
         CALL compute_d_lP_R_analytic(bs_env, grid_points, cells_R_d_lP, d_lP_R, qpoints, charge_lq)
      ELSE
         CALL compute_d_lP_R_numeric(bs_env, grid_points, cells_R_d_lP, d_lP_R, qpoints, charge_lq)
      END IF
      CALL mp_print_mem_per_rank(bs_env%para_env, bs_env%unit_nr, 'After RI-RS d_lP computation')

   END SUBROUTINE compute_d_lP_R

! **************************************************************************************************
!> \brief Enumerate translations T that survive the R_cut^(lP) criterion for at least one (l,P).
!> \param bs_env ...
!> \param grid_points ...
!> \param cells_R_d_lP ...
! **************************************************************************************************
   SUBROUTINE build_cells_R_for_d_lP_R(bs_env, grid_points, cells_R_d_lP)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: grid_points
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: cells_R_d_lP(:, :)

      INTEGER                                            :: cell_index(3), i_atom, i_l, ix, iy, iz, &
                                                            n_max(3), n_RI, n_work
      INTEGER, ALLOCATABLE                               :: cells_work(:, :)
      LOGICAL                                            :: keep
      REAL(KIND=dp)                                      :: center(3), cutoff, max_cutoff
      REAL(KIND=dp), ALLOCATABLE                         :: automatic_cutoffs(:)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      particle_set => bs_env%ri_rs%particle_set
      CALL compute_d_lP_cutoffs(bs_env, automatic_cutoffs)
      max_cutoff = MAXVAL(automatic_cutoffs)
      DO ix = 1, 3
         IF (bs_env%periodic(ix) == 1) THEN
            ! If |h*T + delta| <= R, then every fractional component obeys
            ! |T_i| <= R*|(h^-1)_i| + |(h^-1 delta)_i|.  Both the atom and the
            ! grid point are mapped into the reference cell, so the last term
            ! is smaller than one.  This reciprocal-vector bound remains safe
            ! for skewed cells, unlike a bound based on direct-vector lengths.
            n_max(ix) = CEILING(max_cutoff*NORM2(bs_env%ri_rs%cell%h_inv(ix, :))) + 1
         ELSE
            n_max(ix) = 0
         END IF
      END DO

      n_work = PRODUCT(2*n_max + 1)
      ALLOCATE (cells_work(3, n_work))
      n_RI = 0
      DO ix = -n_max(1), n_max(1)
         DO iy = -n_max(2), n_max(2)
            DO iz = -n_max(3), n_max(3)
               cell_index = [ix, iy, iz]
               keep = .FALSE.
               DO i_atom = 1, bs_env%n_atom
                  cutoff = automatic_cutoffs(i_atom)
                  center = pbc(particle_set(i_atom)%r, bs_env%ri_rs%cell) + &
                           MATMUL(bs_env%hmat, REAL(cell_index, dp))
                  DO i_l = 1, bs_env%ri_rs%n_grid_points
                     IF (NORM2(center - grid_points(:, i_l)) <= cutoff) THEN
                        keep = .TRUE.
                        EXIT
                     END IF
                  END DO
                  IF (keep) EXIT
               END DO
               IF (.NOT. keep) CYCLE
               n_RI = n_RI + 1
               cells_work(:, n_RI) = cell_index
            END DO
         END DO
      END DO

      ALLOCATE (cells_R_d_lP(3, n_RI))
      cells_R_d_lP(:, :) = cells_work(:, 1:n_RI)
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T69,F10.1,A)') &
            'Periodic RI-RS automatic R_cut^(lP):', &
            max_cutoff*angstrom, ' Å'
      END IF
      DEALLOCATE (automatic_cutoffs, cells_work)

   END SUBROUTINE build_cells_R_for_d_lP_R

END MODULE ri_rs_kp_compute_d_lP
