!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Floquet stuff
!> \par History
!> \author Shridhar Shanbhag (27.01.2026)
! **************************************************************************************************
MODULE floquet_main
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_type
   USE floquet_utils,                   ONLY: build_floquet_matrix,&
                                              calculate_floquet_spectral_function,&
                                              check_floquet_conversion,&
                                              write_floquet_spectral_function,&
                                              write_quasi_energy
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE physcon,                         ONLY: a_bohr,&
                                              evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'floquet_main'

   ! Public subroutines
   PUBLIC :: floquet

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE floquet(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'floquet'

      INTEGER                                            :: handle, ikp, ikp_for_file, ispin, &
                                                            max_f_index, n_E, n_f_size, n_fbands, &
                                                            n_spin, nao, nkp, nkp_start, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: indices
      REAL(KIND=dp)                                      :: energy_step, energy_window
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a_k, eigenvalues
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_type)                                  :: cfm_eigenvectors, floquet_copy, &
                                                            floquet_matrix
      TYPE(cp_fm_struct_type), POINTER                   :: floquet_struct
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      para_env=para_env, &
                      mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)

      unit_nr = bs_env%unit_nr
      nkp = bs_env%nkp_bs_and_DOS
      nkp_start = bs_env%nkp_only_DOS
      n_spin = bs_env%n_spin
      max_f_index = bs_env%max_floquet_index
      n_fbands = 1 + 2*max_f_index
      n_f_size = nao*n_fbands

      energy_step = bs_env%energy_step_floquet
      energy_window = 2*bs_env%energy_window_floquet

      n_E = INT(energy_window/energy_step)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(T2,A)') ' '
         WRITE (unit_nr, '(T2,A)') REPEAT('-', 79)
         WRITE (unit_nr, '(T2,A,A78)') '-', '-'
         WRITE (unit_nr, '(T2,A,A51,A27)') '-', 'FLOQUET BANDSTRUCTURE CALCULATION', '-'
         WRITE (unit_nr, '(T2,A,A78)') '-', '-'
         WRITE (unit_nr, '(T2,A)') REPEAT('-', 79)
         WRITE (unit_nr, '(T2,A)') ' '

         WRITE (unit_nr, '(T2,A,T56,ES12.2E2)') &
            "FLOQUET PARAMETERS      Amplitude [V/m]:", bs_env%floquet_amplitude*evolt/a_bohr
         WRITE (unit_nr, '(T2,A,T56,F12.4)') &
            "                        Frequency [eV]:", bs_env%floquet_omega*evolt
         WRITE (unit_nr, '(T2,A)') "                        "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T56,3F8.4)') &
            "                        Polarisation:", bs_env%floquet_polarisation(1:3)
         WRITE (unit_nr, '(T2,A,T56,3F8.4)') &
            "                        Phase offsets:", pi*bs_env%floquet_phi(1:3)
         WRITE (unit_nr, '(T2,A)') "                        "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T56,I12)') &
            "                        Max Floquet index:", bs_env%max_floquet_index
         WRITE (unit_nr, '(T2,A,T56,I12)') &
            "                        Floquet Hamiltonian Size:", n_f_size
         WRITE (unit_nr, '(T2,A)') "                        "//REPEAT("-", 44)
         WRITE (unit_nr, '(T2,A,T56,F12.4)') &
            "                        Energy window [eV]:", bs_env%energy_window_floquet*evolt
         WRITE (unit_nr, '(T2,A,T56,F12.4)') &
            "                        Energy step [eV]:", bs_env%energy_step_floquet*evolt
         WRITE (unit_nr, '(T2,A,T56,F12.4)') &
            "                        Broadening [eV]:", bs_env%broadening_floquet*evolt
         WRITE (unit_nr, '(T2,A)') "                        "//REPEAT("-", 44)
         WRITE (unit_nr, '(A)') ""

         WRITE (unit_nr, '(T2,A)') &
            "We construct the Floquet-Bloch Hamiltonian and diagonalise to obtain the"
         WRITE (unit_nr, '(T2,A)') &
            "eigenvalues which give us the quasi-energies stored in QUASI_ENERGIES.bs"
         WRITE (unit_nr, '(A)') ""
         WRITE (unit_nr, '(T2,A)') &
            "The k-resolved density of states is obtained by computing the trace of"
         WRITE (unit_nr, '(T2,A)') &
            "the retarded Green's function and stored in FLOQUET_DOS.out"
         WRITE (unit_nr, '(T2,A)') &
            "DOS(ω,k) = -1/π*Im[Tr_KS(G^R(ω,k))]"
         WRITE (unit_nr, '(A)') ""
         WRITE (unit_nr, '(T2,A)') "Calculations completed for:"
         WRITE (unit_nr, '(2X,A,T12,A,T22,A)') "KPoint", "Spin", "Coordinate"
      END IF

      CALL cp_fm_struct_create(floquet_struct, &
                               nrow_global=n_f_size, &
                               ncol_global=n_f_size, &
                               context=blacs_env, &
                               para_env=para_env)

      CALL cp_cfm_create(floquet_matrix, floquet_struct, set_zero=.TRUE.)
      CALL cp_cfm_create(floquet_copy, floquet_struct, set_zero=.TRUE.)
      CALL cp_cfm_create(cfm_eigenvectors, floquet_struct, set_zero=.TRUE.)
      ALLOCATE (eigenvalues(n_f_size), indices(nao))
      ALLOCATE (a_k(n_E))
      DO ikp = nkp_start + 1, nkp

         xkp(1:3) = bs_env%kpoints_DOS%xkp(1:3, ikp)
         ikp_for_file = ikp - nkp_start

         DO ispin = 1, n_spin
            ! Build the Floquet-Bloch Hamiltonian matrix H_F(k), with the equilibrium
            ! band energies ε_{nk} on the diagonal (shifted by mħΩ per sector) and
            ! the light-matter coupling blocks A . p_uv on the off-diagonals
            CALL build_floquet_matrix(qs_env, bs_env, xkp, ispin, floquet_matrix)

            ! Use a copy of the Floquet matrix for diagonalization, since
            ! cp_cfm_heevd overwrites it (needed later for the Green's function)
            CALL cp_cfm_to_cfm(floquet_matrix, floquet_copy)

            ! Diagonalize the Floquet-Bloch matrix H_F(k) F_{α,k} = ε_α F_{α,k}
            ! to obtain the quasi-energies ε_α (eigenvalues) and the
            ! Floquet-Bloch coefficients F^{nu}_{α,k} (eigenvectors)
            CALL cp_cfm_set_all(cfm_eigenvectors, z_zero)
            CALL cp_cfm_heevd(floquet_copy, cfm_eigenvectors, eigenvalues)

            ! Run a conversion check to ensure that MAX_FLOQUET_INDEX,
            ! the size of the truncated matrix H_F(k), is sufficiently large.
            CALL check_floquet_conversion(qs_env, bs_env, cfm_eigenvectors)

            ! Write the quasi-energies ε_α folded into the first Floquet Brillouin
            ! zone -ħΩ/2 < ε_α ≤ ħΩ/2 for this k-point and spin channel
            CALL write_quasi_energy(qs_env, bs_env, ispin, ikp_for_file, xkp, eigenvalues)

            ! Compute the physical spectral function
            ! A_phys(k,ω) = -(1/π) Im[ Tr_KS G^R_{00}(k,ω) ]
            ! from the zero-zero Floquet sector of the retarded Green's function
            CALL calculate_floquet_spectral_function(qs_env, bs_env, floquet_matrix, a_k)

            ! Write A_phys(k,ω) (or f(ω,μ)·A_phys for the occupied spectral weight)
            ! to output for this k-point and spin channel
            CALL write_floquet_spectral_function(qs_env, bs_env, ispin, ikp_for_file, xkp, a_k)

            ! Write progress in the main output file
            IF (unit_nr > 0) WRITE (unit_nr, '(2X,I6,T12,I4,T18,3F12.6)') ikp_for_file, ispin, xkp(1:3)

         END DO
      END DO

      IF (unit_nr > 0) WRITE (unit_nr, '(A)') ""
      DEALLOCATE (eigenvalues)
      CALL cp_cfm_release(floquet_matrix)
      CALL cp_cfm_release(floquet_copy)
      CALL cp_cfm_release(cfm_eigenvectors)

      CALL timestop(handle)

   END SUBROUTINE floquet

END MODULE floquet_main
