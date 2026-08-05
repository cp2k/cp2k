!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Environment type holding the work and accumulation arrays of the Floquet-Bloch
!>        band-structure calculation (floquet_main), plus the derived sizes.
!> \par History
!> \author Shridhar Shanbhag (27.01.2026)
! **************************************************************************************************
MODULE floquet_types
   USE kinds,                           ONLY: dp
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'floquet_types'

   PUBLIC :: floquet_env_type, &
             floquet_env_create, &
             floquet_env_release

! **************************************************************************************************
!> \brief Work and accumulation arrays for one Floquet-Bloch band-structure run, plus the derived
!>        sizes. The per-(k-point,spin) buffers are overwritten each iteration; the all_* arrays
!>        accumulate the results over all DOS k-points and spins.
!> \param max_f_index Floquet truncation index M
!> \param n_f_size size of the truncated Floquet Hamiltonian, nao*(1 + 2*M)
!> \param n_E number of points on the DOS energy grid
!> \param nkp_only_bs number of band-structure k-points (the k-points the Floquet loop runs over)
!> \param nao number of atomic orbitals
!> \param n_spin number of spin channels
!> \param eigenvalues Floquet eigenvalues for the current (k-point,spin) (n_f_size)
!> \param w0 central-sector weights (n_f_size)
!> \param wE outermost-sector weights (n_f_size)
!> \param a_k spectral function A(k,E) on the energy grid (n_E)
!> \param quasi_energies quasi-energies (nao)
!> \param m0_energies m=0 band energies (nao)
!> \param m0_weights central-sector weight of each m=0 band (nao)
!> \param e_k_kp_spin band energies for the current (k-point,spin) (nao)
!> \param de_dk_kp_spin band-energy k-derivatives for the current (k-point,spin) (3, nao)
!> \param dipole_kp_spin dipole matrix elements for the current (k-point,spin) (3, nao, nao)
!> \param all_quasi_energies quasi-energies for every (band, spin, bs k-point) (nao, n_spin, nkp_only_bs)
!> \param all_a_k spectral function for every (energy, spin, bs k-point) (n_E, n_spin, nkp_only_bs)
!> \param all_m0_energies m=0 band energies for every (band, spin, bs k-point) (nao, n_spin, nkp_only_bs)
!> \param all_m0_weights central-sector weights for every (band, spin, bs k-point) (nao, n_spin, nkp_only_bs)
! **************************************************************************************************
   TYPE floquet_env_type
      INTEGER                                          :: max_f_index = -1, &
                                                          n_f_size = -1, &
                                                          n_E = -1, &
                                                          nkp_only_bs = -1, &
                                                          nao = -1, &
                                                          n_spin = -1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)         :: eigenvalues, w0, wE, a_k, &
                                                          quasi_energies, m0_energies, &
                                                          m0_weights, e_k_kp_spin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)      :: de_dk_kp_spin
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: dipole_kp_spin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)   :: all_quasi_energies, all_a_k, &
                                                          all_m0_energies, all_m0_weights
   END TYPE floquet_env_type

CONTAINS

! **************************************************************************************************
!> \brief Set the derived sizes from bs_env and allocate all work and accumulation arrays.
!> \param floquet_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE floquet_env_create(floquet_env, bs_env)
      TYPE(floquet_env_type), INTENT(OUT)                :: floquet_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'floquet_env_create'

      INTEGER                                            :: handle, nao, nE, nf, nk, ns

      CALL timeset(routineN, handle)

      floquet_env%nao = bs_env%n_ao
      floquet_env%n_spin = bs_env%n_spin
      floquet_env%max_f_index = bs_env%max_floquet_index
      floquet_env%n_f_size = floquet_env%nao*(1 + 2*floquet_env%max_f_index)
      floquet_env%n_E = NINT(2*bs_env%energy_window_floquet/bs_env%energy_step_floquet)
      floquet_env%nkp_only_bs = bs_env%nkp_only_bs

      nao = floquet_env%nao
      nf = floquet_env%n_f_size
      nE = floquet_env%n_E
      ns = floquet_env%n_spin
      nk = floquet_env%nkp_only_bs

      ! per-(k-point,spin) work buffers
      ALLOCATE (floquet_env%eigenvalues(nf), floquet_env%w0(nf), floquet_env%wE(nf))
      ALLOCATE (floquet_env%a_k(nE))
      ALLOCATE (floquet_env%quasi_energies(nao), floquet_env%m0_energies(nao), &
                floquet_env%m0_weights(nao))
      ALLOCATE (floquet_env%e_k_kp_spin(nao), floquet_env%de_dk_kp_spin(3, nao), &
                floquet_env%dipole_kp_spin(3, nao, nao))

      ! per-(spin,k-point) accumulation arrays, filled once by the owning subgroup source and
      ! summed across the global communicator after the k-point loop
      ALLOCATE (floquet_env%all_quasi_energies(nao, ns, nk), source=0.0_dp)
      ALLOCATE (floquet_env%all_a_k(nE, ns, nk), source=0.0_dp)
      ALLOCATE (floquet_env%all_m0_energies(nao, ns, nk), source=0.0_dp)
      ALLOCATE (floquet_env%all_m0_weights(nao, ns, nk), source=0.0_dp)

      CALL timestop(handle)

   END SUBROUTINE floquet_env_create

! **************************************************************************************************
!> \brief Deallocate all arrays held by floquet_env.
!> \param floquet_env ...
! **************************************************************************************************
   SUBROUTINE floquet_env_release(floquet_env)
      TYPE(floquet_env_type), INTENT(INOUT)              :: floquet_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'floquet_env_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ALLOCATED(floquet_env%eigenvalues)) DEALLOCATE (floquet_env%eigenvalues)
      IF (ALLOCATED(floquet_env%w0)) DEALLOCATE (floquet_env%w0)
      IF (ALLOCATED(floquet_env%wE)) DEALLOCATE (floquet_env%wE)
      IF (ALLOCATED(floquet_env%a_k)) DEALLOCATE (floquet_env%a_k)
      IF (ALLOCATED(floquet_env%quasi_energies)) DEALLOCATE (floquet_env%quasi_energies)
      IF (ALLOCATED(floquet_env%m0_energies)) DEALLOCATE (floquet_env%m0_energies)
      IF (ALLOCATED(floquet_env%m0_weights)) DEALLOCATE (floquet_env%m0_weights)
      IF (ALLOCATED(floquet_env%e_k_kp_spin)) DEALLOCATE (floquet_env%e_k_kp_spin)
      IF (ALLOCATED(floquet_env%de_dk_kp_spin)) DEALLOCATE (floquet_env%de_dk_kp_spin)
      IF (ALLOCATED(floquet_env%dipole_kp_spin)) DEALLOCATE (floquet_env%dipole_kp_spin)
      IF (ALLOCATED(floquet_env%all_quasi_energies)) DEALLOCATE (floquet_env%all_quasi_energies)
      IF (ALLOCATED(floquet_env%all_a_k)) DEALLOCATE (floquet_env%all_a_k)
      IF (ALLOCATED(floquet_env%all_m0_energies)) DEALLOCATE (floquet_env%all_m0_energies)
      IF (ALLOCATED(floquet_env%all_m0_weights)) DEALLOCATE (floquet_env%all_m0_weights)

      CALL timestop(handle)

   END SUBROUTINE floquet_env_release

END MODULE floquet_types
