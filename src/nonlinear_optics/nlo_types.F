!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \brief Nonlinear optical environment: input parameters and data shared by calculation stages.
! **************************************************************************************************
MODULE nlo_types
   USE cell_types,                      ONLY: cell_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_release,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_mo_types,                     ONLY: mo_set_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlo_types'

   PUBLIC                                            :: release_nlo_environment
   ! A frequency is either a regular energy grid or a signed link to an earlier field.
   TYPE, PUBLIC                                      :: nlo_frequency_type
      REAL(KIND=dp)                                  :: start_ev = 0.0_dp, &
                                                        step_ev = 0.0_dp
      INTEGER                                        :: n_values = 0, &
                                                        reference = 0, &
                                                        sign = 1
   END TYPE nlo_frequency_type

   ! One instance, nlo_env, owns the input and the data shared by the five main stages.
   TYPE, PUBLIC                                      :: nlo_type
      INTEGER                                        :: order = 0, &
                                                        k_mesh(3) = 0
      LOGICAL, ALLOCATABLE                           :: elements(:, :)
      REAL(KIND=dp)                                  :: t1_t2_au = 0.0_dp
      TYPE(nlo_frequency_type)                       :: frequency(5)
      ! Number of independent tuples formed from FREQ_1,...,FREQ_ORDER.
      INTEGER                                        :: n_frequency_combinations = 0
      ! Data shared by preparation, Brillouin-zone integration and output.
      INTEGER                                        :: n_bands = 0, &
                                                        n_ao = 0, &
                                                        n_images = 0, &
                                                        n_kpoints = 0, &
                                                        n_spin = 0

      REAL(KIND=dp)                                  :: min_direct_gap = 0.0_dp, &
                                                        projection_error = 0.0_dp
      ! AO matrices at lattice image R: H, S, and dipole d.
      ! Indices: (spin or Cartesian component, AO row, AO column, R).
      REAL(KIND=dp), ALLOCATABLE                     :: h_ao_R(:, :, :, :), &
                                                        s_ao_R(:, :, :, :), &
                                                        d_ao_R(:, :, :, :)

      ! Integrated current: (output direction, input-direction code, frequency combination).
      COMPLEX(dp), ALLOCATABLE                       :: current_spectrum(:, :, :)

      ! Temporary matrices for one k point, allocated once for the integration stage.
      COMPLEX(dp), ALLOCATABLE                       :: band_connection_derivatives(:, :, :, :), &
                                                        band_hamiltonian_derivatives(:, :, :), &
                                                        current_flat(:, :), &
                                                        position_derivatives(:, :, :, :), &
                                                        hamiltonian_derivatives(:, :, :), &
                                                        rho0_derivatives(:, :, :), &
                                                        overlap_derivatives(:, :, :), &
                                                        velocity(:, :, :)
      INTEGER, ALLOCATABLE                           :: occupied(:)
      REAL(KIND=dp), ALLOCATABLE                     :: energies(:)

      ! Number of spin/k-point tasks assigned to the reporting rank.
      INTEGER                                        :: progress_done = 0, &
                                                        progress_next = 1, &
                                                        progress_total = 0
      REAL(KIND=dp)                                  :: progress_start = 0.0_dp

      TYPE(cell_type), POINTER                       :: cell => NULL()
      TYPE(kpoint_type), POINTER                     :: images => NULL()
      TYPE(mo_set_type), POINTER                     :: mos(:) => NULL()
      TYPE(mp_para_env_type), POINTER                :: para_env => NULL()
   END TYPE nlo_type

CONTAINS

! **************************************************************************************************
!> \brief Release NLO-owned matrices and the auxiliary k-point environment, if created.
!> \param nlo_env NLO settings and owned arrays to release
! **************************************************************************************************
   SUBROUTINE release_nlo_environment(nlo_env)
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env

      CALL kpoint_release(nlo_env%images)
      NULLIFY (nlo_env%cell, nlo_env%mos, nlo_env%para_env)
      DEALLOCATE (nlo_env%h_ao_R, nlo_env%s_ao_R, nlo_env%d_ao_R, nlo_env%current_spectrum)
      IF (ALLOCATED(nlo_env%elements)) DEALLOCATE (nlo_env%elements)

   END SUBROUTINE release_nlo_environment

END MODULE nlo_types
