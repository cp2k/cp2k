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
   USE mathconstants,                   ONLY: z_zero
   USE message_passing,                 ONLY: mp_para_env_type
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

      INTEGER                                            :: handle, ikp, ispin, max_f_index, n_E, &
                                                            n_f_size, n_fbands, n_spin, nao, nkp, &
                                                            nkp_start
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

      nkp = bs_env%nkp_bs_and_DOS
      nkp_start = bs_env%nkp_only_DOS
      n_spin = bs_env%n_spin
      max_f_index = bs_env%max_floquet_index
      n_fbands = 1 + 2*max_f_index
      n_f_size = nao*n_fbands

      energy_step = bs_env%energy_step_floquet
      energy_window = 2*bs_env%energy_window_floquet

      n_E = INT(energy_window/energy_step)

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
         DO ispin = 1, n_spin

            xkp(1:3) = bs_env%kpoints_DOS%xkp(1:3, ikp)

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

            ! TODO make a conversion check with cfm_eigenvectors.
            CALL check_floquet_conversion(qs_env, bs_env, cfm_eigenvectors)

            ! Write the quasi-energies ε_α folded into the first Floquet Brillouin
            ! zone -ħΩ/2 < ε_α ≤ ħΩ/2 for this k-point and spin channel
            CALL write_quasi_energy(qs_env, bs_env, ispin, ikp, xkp, eigenvalues)

            ! Compute the physical spectral function
            ! A_phys(k,ω) = -(1/π) Im[ Tr_KS G^R_{00}(k,ω) ]
            ! from the zero-zero Floquet sector of the retarded Green's function
            CALL calculate_floquet_spectral_function(qs_env, bs_env, floquet_matrix, a_k)

            ! Write A_phys(k,ω) (or f(ω,μ)·A_phys for the occupied spectral weight)
            ! to output for this k-point and spin channel
            CALL write_floquet_spectral_function(qs_env, bs_env, ispin, ikp, xkp, a_k)
         END DO
      END DO
      DEALLOCATE (eigenvalues)
      CALL cp_cfm_release(floquet_matrix)
      CALL cp_cfm_release(floquet_copy)
      CALL cp_cfm_release(cfm_eigenvectors)

      CALL timestop(handle)

   END SUBROUTINE floquet

END MODULE floquet_main
