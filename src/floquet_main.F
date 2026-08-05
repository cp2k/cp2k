!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Floquet-Bloch band-structure calculations: builds and diagonalises the
!>        truncated Floquet-Bloch Hamiltonian H_F(k) for each DOS k-point and writes the m=0
!>        band structure, the quasi-energies and the Floquet density of states.
!> \par History
!> \author Shridhar Shanbhag (27.01.2026)
! **************************************************************************************************
MODULE floquet_main
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE floquet_types,                   ONLY: floquet_env_create,&
                                              floquet_env_release,&
                                              floquet_env_type
   USE floquet_utils,                   ONLY: &
        build_floquet_matrix, calculate_floquet_observables, check_floquet_convergence, &
        compute_e_k_de_dk_dipole, distribute_floquet_kp_data, floquet_sector_weights, &
        make_floquet_subgroups, write_floquet_header, write_floquet_results
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: z_zero
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'floquet_main'

   PUBLIC :: floquet

CONTAINS

! **************************************************************************************************
!> \brief Computes the Floquet-Bloch band structure, density of states and quasi-energies.
!> \param qs_env ...
!> \param bs_env ...
! **************************************************************************************************
   SUBROUTINE floquet(qs_env, bs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'floquet'

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: dipole
      INTEGER                                            :: handle, ikp, ispin, my_group, &
                                                            n_est_done, n_source_done, ngroups
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: group_distribution
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: e_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: de_dk
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      TYPE(cp_cfm_type)                                  :: cfm_eigenvectors, floquet_matrix
      TYPE(cp_fm_struct_type), POINTER                   :: floquet_struct
      TYPE(floquet_env_type)                             :: floquet_env
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_sub

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env)

      ! Floquet-Bloch theory: under a monochromatic electric field drive the Hamiltonian is periodic,
      ! H(k,t) = H(k,t+T) with T = 2π/Ω, where Ω is the angular frequency of the drive.
      ! By Floquet's theorem the solutions are:
      ! ψ_αk(t) = e^(-i λ_αk t) Σ_m e^(-i m Ω t) |F^m_αk>, with quasi-energies λ_αk and
      ! Fourier index m.
      ! In order to obtain the solutions λ_αk, we start with the Floquet Hamiltonian defined as,
      ! H_F(k,t) = H(k,t) - iħ d/dt.
      ! Taking the Fourier transform, the solutions are obtained by diagonalizing the matrix H_F(k)
      !   H_F(k) F_αk = λ_αk F_αk,
      ! where H_F(k) is block-tridiagonal in the sector index m = -M..M (M = MAX_FLOQUET_INDEX)
      ! and each block is nao×nao. We build and diagonalise H_F(k) per DOS k-point k and spin.
      ! The matrix H_F(k) has size (nao*(2*MAX_FLOQUET_INDEX+1) x nao*(2*MAX_FLOQUET_INDEX+1)).

      CALL write_floquet_header(bs_env)

      CALL floquet_env_create(floquet_env, bs_env)

      ! Equilibrium band data at every band-structure k-point k from H(k)C(k) = S(k)C(k)ε(k):
      !   band energies    ε_nk
      !   k-derivatives    ∇_k ε_nk = C_n(k)^† ∇_k H(k) C_n(k) - ε_nk C_n(k)^† ∇_k S(k) C_n(k)
      !   dipoles          d_nm(k) = <ψ_nk| r |ψ_mk>
      ! Results are distributed across ranks (k-point ikp stored on rank MOD(ikp-1,num_pe)).
      CALL compute_e_k_de_dk_dipole(qs_env, bs_env, e_k, de_dk, dipole)

      ! Split the global communicator into Floquet subgroups (Hamiltonian distributed within each
      ! subgroup, k-points distributed across subgroups) and create the subgroup BLACS context.
      CALL make_floquet_subgroups(qs_env, bs_env, para_env_sub, blacs_env_sub, &
                                  group_distribution, ngroups)

      !get inside make_floquet_subroutine
      my_group = group_distribution(para_env%mepos)

      CALL cp_fm_struct_create(floquet_struct, context=blacs_env_sub, para_env=para_env_sub, &
                               nrow_global=floquet_env%n_f_size, ncol_global=floquet_env%n_f_size)
      CALL cp_cfm_create(floquet_matrix, floquet_struct, set_zero=.TRUE.)
      CALL cp_cfm_create(cfm_eigenvectors, floquet_struct, set_zero=.TRUE.)

      n_source_done = 0

      DO ikp = 1, floquet_env%nkp_only_bs
         ! Only the subgroup that owns this k-point's data works on it
         IF (group_distribution(MOD(ikp - 1, para_env%num_pe)) /= my_group) CYCLE

         DO ispin = 1, floquet_env%n_spin
            ! Broadcast this k-point/spin's ε_nk, ∇_k ε_nk and d_nm(k) to every rank of the
            ! owning subgroup (its owner rank holds them from the distributed precompute above).
            CALL distribute_floquet_kp_data(para_env, para_env_sub, ispin, ikp, e_k, de_dk, &
                                            dipole, floquet_env)

            ! Assemble H_F(k) in Sambe space (blocks m,m' = -M..M, each nao×nao)
            !  diagonal      (H_F)^{m,m}_{nn'}   = (ε_nk + m ħΩ) δ_{nn'}
            !  off diagonal  (H_F)^{m,m±1}_{nn'} = Σ_a p^a_{nn'}(k) (i E_a e^{±iφ_a})/(2Ω),  a=x,y,z
            ! with momentum   p^a_{nn'}(k) = i d^a_{nn'}(k)(ε_nk - ε_n'k) + ∇_{k_a} ε_nk δ_{nn'}
            ! and field       E_a = |E| pol_a.
            CALL build_floquet_matrix(bs_env, floquet_env, floquet_matrix)

            ! Diagonalise H_F(k) F_αk = λ_αk F_αk -> quasi-energies λ_αk and eigenvector components
            ! F^{n,m}_αk = <n,m|αk> (band index n, sector m).
            CALL cp_cfm_set_all(cfm_eigenvectors, z_zero)
            CALL cp_cfm_heevd(floquet_matrix, cfm_eigenvectors, floquet_env%eigenvalues)

            ! Calculate sector weights,
            ! w_αk^(m) = Σ_{n=1}^{nao} |<n,m|αk>|^2 (with Σ_m w_αk^(m) = 1).
            ! Keep only the central weight w0_αk = w_αk^(0) and the outermost-rung weight
            ! wE_αk = w_αk^(-M) + w_αk^(+M)
            CALL floquet_sector_weights(floquet_env, cfm_eigenvectors)

            ! Truncation check: the m=0-dominated states (w0_αk > 1/2) must leak negligibly onto the
            ! outermost rungs, max_{α: w0_αk > 1/2} wE_αk ≤ EPS_FLOQUET, else M is too small.
            CALL check_floquet_convergence(bs_env, floquet_env)

            ! From λ_αk and w0_αk obtain, and store on the subgroup source for (ispin,ikp):
            !   DOS:
            !      A(k,E) = -(1/π) Im Tr_KS G^R_00(k,E) = -(1/π) Σ_α w0_αk Im[1/(E+iη-λ_αk)]
            !      with G^R_00 = [(E+iη) - H_F(k)]^{-1}_{m=0,m=0}, η = broadening/2.
            !      The central (m=0) block of the retarded Floquet Green's function gives the
            !      physical, time-averaged spectral function / density of states; see
            !      T. Oka and S. Kitamura, Annu. Rev. Condens. Matter Phys. 10, 387 (2019), and
            !      N. Tsuji, T. Oka, and H. Aoki, Phys. Rev. B 78, 235124 (2008).
            !   m=0 bands:
            !      ε^{0}_ik = λ_{α(i)}k, with band i placed at the eigenvalue where the cumulative
            !      m=0 weight Σ_{λ_αk<E} w0_αk first crosses i - 1/2.
            !      We know that: Σ_α w0_αk = nao. So we choose the matrix eigenvalue such that:
            !      N(ε) = Σ_{α : λ_αk ≤ ε} w0_αk first crosses i - 1/2 to find the m=0 band.
            !   Quasi-energy:
            !      ε̃_ik = ε^{0}_ik - mħΩ (choose m such that ε̃_ik is folded to (-ħΩ/2, ħΩ/2]).
            CALL calculate_floquet_observables(bs_env, para_env_sub, ispin, ikp, floquet_env)
         END DO

         ! Rough live K-Point progress
         IF (bs_env%unit_nr > 0) THEN
            n_source_done = n_source_done + 1
            n_est_done = MIN(n_source_done*ngroups, floquet_env%nkp_only_bs)
            WRITE (bs_env%unit_nr, '(T2,A,I6,A,F6.1,A)') &
               "FLOQUET PROGRESS | Estimated no. of k-points done: ", n_est_done, "  (", &
               100.0_dp*REAL(n_est_done, dp)/REAL(floquet_env%nkp_only_bs, dp), " %)"
         END IF
      END DO

      ! Write the m=0 bands ε^{0}_ik, the quasi-energies ε̃_ik and the DOS A(k,E) in their
      ! respective files
      CALL write_floquet_results(bs_env, floquet_env)

      DEALLOCATE (e_k, de_dk, dipole)
      DEALLOCATE (group_distribution)
      CALL floquet_env_release(floquet_env)
      CALL cp_cfm_release(floquet_matrix)
      CALL cp_cfm_release(cfm_eigenvectors)
      CALL cp_fm_struct_release(floquet_struct)
      CALL cp_blacs_env_release(blacs_env_sub)
      CALL mp_para_env_release(para_env_sub)

      CALL timestop(handle)

   END SUBROUTINE floquet

END MODULE floquet_main
