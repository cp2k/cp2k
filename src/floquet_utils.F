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
MODULE floquet_utils
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_invert,&
                                              cp_cfm_scale_and_add
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_submatrix,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_set_submatrix,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_k_r_trafo_simple,         ONLY: replicate_rs_matrices,&
                                              rs_to_kp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              pi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: geeig_right,&
                                              gemm_square
   USE message_passing,                 ONLY: mp_para_env_type
   USE physcon,                         ONLY: evolt
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_moments,                      ONLY: qs_moment_kpoints_deep
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'floquet_utils'

   ! Public subroutines
   PUBLIC :: build_floquet_matrix, &
             calculate_epsilon_derivative, &
             build_momentum_matrix, &
             check_floquet_conversion, &
             calculate_floquet_spectral_function, &
             write_quasi_energy, &
             write_floquet_spectral_function

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param xkp ...
!> \param e_k ...
!> \param de_dk ...
! **************************************************************************************************
   SUBROUTINE calculate_epsilon_derivative(qs_env, xkp, e_k, de_dk)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         OPTIONAL                                        :: e_k
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :), OPTIONAL                 :: de_dk

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_epsilon_derivative'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: C_dH_C, C_dS_C, C_k, dH_dk_i, dS_dk_i, &
                                                            H_k, S_k
      INTEGER                                            :: handle, i_dir, ikp, ispin, n, n_img_all, &
                                                            n_spin, nao, nkp
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_all
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_all
      LOGICAL                                            :: present_dedk, present_ek
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvals
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: H_rs, S_rs
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp, matrix_s_kp
      TYPE(kpoint_type), POINTER                         :: kpoints_all, kpoints_scf
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all

      CALL timeset(routineN, handle)

      present_ek = PRESENT(e_k)  ! calculate band energies, ε_k for all kpoints xkp
      present_dedk = PRESENT(de_dk)  ! calculate derivative, ∇_k ε_k of band energies for all kpoints xkp

      IF (.NOT. (present_ek .OR. present_dedk)) CPABORT("Subroutine needs either e_k or de_dk")

      CALL get_qs_env(qs_env, &
                      matrix_ks_kp=matrix_ks_kp, &
                      matrix_s_kp=matrix_s_kp, &
                      sab_all=sab_all, &
                      cell=cell, &
                      kpoints=kpoints_scf, &
                      para_env=para_env, &
                      mos=mos)

      CALL get_mo_set(mo_set=mos(1), nao=nao)
      CALL get_cell(cell=cell, h=hmat)

      n_spin = SIZE(matrix_ks_kp, 1)
      nkp = SIZE(xkp, 2)

      ! create kpoint environment kpoints_all which contains all neighbor cells R
      ! without considering any lattice symmetry
      NULLIFY (kpoints_all)
      CALL kpoint_create(kpoints_all)
      CALL kpoint_init_cell_index(kpoints_all, sab_all, para_env, n_img_all)
      CALL get_kpoint_info(kpoints_all, cell_to_index=cell_to_index_all, index_to_cell=index_to_cell_all)

      ALLOCATE (S_rs(1, nao, nao, n_img_all), H_rs(n_spin, nao, nao, n_img_all), source=0.0_dp)

      ! Convert real-space dbcsr matrices into arrays
      CALL replicate_rs_matrices(matrix_s_kp, kpoints_scf, S_rs, cell_to_index_all)
      CALL replicate_rs_matrices(matrix_ks_kp, kpoints_scf, H_rs, cell_to_index_all)

      IF (present_dedk) ALLOCATE (de_dk(n_spin, nkp, 3, nao), source=0.0_dp)
      IF (present_ek) ALLOCATE (e_k(n_spin, nkp, nao), source=0.0_dp)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(ikp, ispin, S_k, H_k, eigenvals, C_k, dS_dk_i, dH_dk_i, C_dS_C, C_dH_C) &
!$OMP SHARED(nao, n_spin, de_dk, e_k, present_ek, present_dedk, &
!$OMP nkp, xkp, S_rs, H_rs, index_to_cell_all, hmat)
      IF (present_dedk) ALLOCATE (dS_dk_i(nao, nao), C_dS_C(nao, nao), &
                                  dH_dk_i(nao, nao), C_dH_C(nao, nao), source=z_zero)
      ALLOCATE (C_k(nao, nao), S_k(nao, nao), H_k(nao, nao), source=z_zero)
      ALLOCATE (eigenvals(nao), source=0.0_dp)
!$OMP DO COLLAPSE(2)
      DO ispin = 1, n_spin
         DO ikp = 1, nkp

            ! S^R -> S(k), H^R -> H(k)
            S_k = 0
            H_k = 0
            CALL rs_to_kp(S_rs(1, :, :, :), S_k, index_to_cell_all, xkp(:, ikp))
            CALL rs_to_kp(H_rs(ispin, :, :, :), H_k, index_to_cell_all, xkp(:, ikp))

            ! Diagonalize H(k)C(k) = S(k)C(k)ε(k)
            CALL geeig_right(H_k, S_k, eigenvals, C_k)
            IF (present_ek) e_k(ispin, ikp, :) = eigenvals(:)

            IF (present_dedk) THEN
               ! Evaluate the derivatives using
               ! ∇ ε_k = C^H(k) ∇ H_k C(k) - ε_k C^H(k) ∇ S_k C(k)
               DO i_dir = 1, 3
                  CALL rs_to_kp(S_rs(1, :, :, :), dS_dk_i, index_to_cell_all, xkp(:, ikp), i_dir, hmat)
                  CALL rs_to_kp(H_rs(ispin, :, :, :), dH_dk_i, index_to_cell_all, xkp(:, ikp), i_dir, hmat)

                  CALL gemm_square(C_k, 'C', dS_dk_i, 'N', C_k, 'N', C_dS_C)
                  CALL gemm_square(C_k, 'C', dH_dk_i, 'N', C_k, 'N', C_dH_C)

                  DO n = 1, nao
                     de_dk(ispin, ikp, i_dir, n) = DBLE(C_dH_C(n, n)) - DBLE(eigenvals(n)*C_dS_C(n, n))
                  END DO
               END DO
            END IF
         END DO
      END DO
!$OMP END DO
      IF (present_dedk) DEALLOCATE (dS_dk_i, C_dS_C, dH_dk_i, C_dH_C)
      DEALLOCATE (S_k, H_k, C_k, eigenvals)
!$OMP END PARALLEL
      DEALLOCATE (S_rs, H_rs)
      CALL kpoint_release(kpoints_all)

      CALL timestop(handle)

   END SUBROUTINE calculate_epsilon_derivative

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param xkp ...
!> \param momentum ...
! **************************************************************************************************
   SUBROUTINE build_momentum_matrix(qs_env, xkp, momentum)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: momentum

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_momentum_matrix'

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :, :)                        :: dipole
      INTEGER                                            :: handle, i, i_dir, ispin, j, n_spin, nao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: e_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: de_dk
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      ! We calculate momentum matrix elements p_nm = <ψ_n|-iħ∇_r|ψ_m>
      ! p_nm = <u_n|(ħk - iħ∇_r)|u_m>
      ! p_nm = i d_nm (ε_n - ε_m) + ∇_k ε_n δ_nm

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)

      ALLOCATE (xkp1(3, 1), source=0.0_dp)
      xkp1(:, 1) = xkp(:)
      CALL calculate_epsilon_derivative(qs_env, xkp1, e_k, de_dk)
      n_spin = dft_control%nspins
      CALL qs_moment_kpoints_deep(qs_env, xkp1, dipole)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(ispin, i_dir, i, j) &
!$OMP SHARED(nao, n_spin, momentum, dipole, e_k, de_dk)
!$OMP DO COLLAPSE(4)
      DO ispin = 1, n_spin
         DO i_dir = 1, 3
            DO i = 1, nao
               DO j = 1, nao
                  IF (j == i) THEN
                     momentum(ispin, i_dir, i, j) = de_dk(ispin, 1, i_dir, i)
                  ELSE
                     momentum(ispin, i_dir, i, j) = gaussi*(e_k(ispin, 1, i) - e_k(ispin, 1, j))* &
                                                    dipole(ispin, 1, i_dir, i, j)
                  END IF
               END DO
            END DO
         END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL
      DEALLOCATE (xkp1, dipole, e_k, de_dk)

      CALL timestop(handle)
   END SUBROUTINE build_momentum_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param xkp ...
!> \param ispin ...
!> \param off_diag_m ...
! **************************************************************************************************
   SUBROUTINE build_off_diagonal_matrix(qs_env, bs_env, xkp, ispin, off_diag_m)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      INTEGER                                            :: ispin
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: off_diag_m

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_off_diagonal_matrix'

      COMPLEX(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :, :, :)                           :: momentum
      COMPLEX(KIND=dp), DIMENSION(3)                     :: efactor
      INTEGER                                            :: handle, i, i_dir, j, n_spin, nao
      REAL(KIND=dp)                                      :: amplitude, omega
      REAL(KIND=dp), DIMENSION(3)                        :: e_vec, phi, polarisation
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      ! Builds the matrix that occupies the off diagonal blocks in the Floquet Matrix H_F

      polarisation(:) = bs_env%floquet_polarisation(:)
      IF (SQRT(SUM(polarisation**2)) < EPSILON(0.0_dp)) THEN
         CPABORT("Invalid (too small) polarisation vector specified for POLARISATION_PUMP")
      END IF

      amplitude = bs_env%floquet_amplitude
      e_vec(:) = amplitude*polarisation
      phi(:) = pi*bs_env%floquet_phi(:)
      omega = bs_env%floquet_omega
      n_spin = bs_env%n_spin

      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)
      ALLOCATE (momentum(n_spin, 3, nao, nao), source=z_zero)
      CALL build_momentum_matrix(qs_env, xkp, momentum)

      ! E_factor(α) = (i E(α) exp(i·φ(α)))/(2ω) where α = x,y,z
      DO i_dir = 1, 3
         efactor(i_dir) = gaussi*e_vec(i_dir)*CMPLX(COS(phi(i_dir)), SIN(phi(i_dir)), KIND=dp)/(2*omega)
      END DO

      DO i_dir = 1, 3
         DO i = 1, nao
            DO j = 1, nao
               off_diag_m(i, j) = off_diag_m(i, j) + &
                                  momentum(ispin, i_dir, i, j)*efactor(i_dir)
            END DO
         END DO
      END DO
      DEALLOCATE (momentum)

      CALL timestop(handle)

   END SUBROUTINE build_off_diagonal_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param xkp ...
!> \param ispin ...
!> \param f_index ...
!> \param diag_e ...
! **************************************************************************************************
   SUBROUTINE build_diagonal_matrix(qs_env, bs_env, xkp, ispin, f_index, diag_e)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      INTEGER                                            :: ispin, f_index
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: diag_e

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_diagonal_matrix'

      INTEGER                                            :: handle, i, n_spin, nao
      REAL(KIND=dp)                                      :: omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: e_k
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      ! Builds the matrix that occupies the off diagonal blocks in the Floquet Matrix H_F

      omega = bs_env%floquet_omega
      n_spin = bs_env%n_spin

      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)
      ALLOCATE (xkp1(3, 1), source=0.0_dp)
      xkp1(:, 1) = xkp(:)
      CALL calculate_epsilon_derivative(qs_env, xkp1, e_k)

      DO i = 1, nao
         diag_e(i, i) = e_k(ispin, 1, i) + f_index*omega
      END DO
      DEALLOCATE (xkp1, e_k)
      CALL timestop(handle)

   END SUBROUTINE build_diagonal_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param xkp ...
!> \param ispin ...
!> \param floquet_matrix ...
! **************************************************************************************************
   SUBROUTINE build_floquet_matrix(qs_env, bs_env, xkp, ispin, floquet_matrix)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      INTEGER                                            :: ispin
      TYPE(cp_cfm_type)                                  :: floquet_matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_floquet_matrix'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: conj_off_diag_m, diag_e, off_diag_m
      INTEGER                                            :: f_index, handle, i, i_f, max_f_index, &
                                                            n_fbands, n_spin, nao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xkp1
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)
      max_f_index = bs_env%max_floquet_index
      n_fbands = 1 + 2*max_f_index
      n_spin = bs_env%n_spin

      ! Creates the floquet matrix H_F by placing the diagonal and off diagonal blocks

      ALLOCATE (xkp1(3, 1), source=0.0_dp)
      ALLOCATE (diag_e(nao, nao), source=z_zero)
      ALLOCATE (off_diag_m(nao, nao), source=z_zero)
      ALLOCATE (conj_off_diag_m(nao, nao), source=z_zero)
      xkp1(:, 1) = xkp(:)
      CALL build_off_diagonal_matrix(qs_env, bs_env, xkp, ispin, off_diag_m)
      conj_off_diag_m(:, :) = CONJG(TRANSPOSE(off_diag_m(:, :)))

      floquet_matrix%local_data(:, :) = z_zero
      DO i = 1, n_fbands
         i_f = 1 + (i - 1)*nao
         f_index = i - max_f_index - 1
         CALL build_diagonal_matrix(qs_env, bs_env, xkp, ispin, f_index, diag_e)
         CALL cp_cfm_set_submatrix(floquet_matrix, diag_e, i_f, i_f)
         IF (i > 1) THEN
            CALL cp_cfm_set_submatrix(floquet_matrix, off_diag_m, i_f, i_f - nao)
            CALL cp_cfm_set_submatrix(floquet_matrix, conj_off_diag_m, i_f - nao, i_f)
         END IF
      END DO
      DEALLOCATE (diag_e, off_diag_m, conj_off_diag_m, xkp1)

      CALL timestop(handle)
   END SUBROUTINE build_floquet_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param cfm_eigenvectors ...
! **************************************************************************************************
   SUBROUTINE check_floquet_conversion(qs_env, bs_env, cfm_eigenvectors)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: cfm_eigenvectors

      CHARACTER(LEN=*), PARAMETER :: routineN = 'check_floquet_conversion'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: vector_block_0, vector_block_1
      INTEGER                                            :: handle, max_f_index, nao, start_row_0, &
                                                            start_row_1
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)

      max_f_index = bs_env%max_floquet_index

      start_row_0 = 1 + nao*max_f_index
      start_row_1 = 1 + nao*(max_f_index + 1)

      ALLOCATE (vector_block_0(3*nao, nao), vector_block_1(3*nao, nao), source=z_zero)
      CALL cp_cfm_get_submatrix(cfm_eigenvectors, vector_block_0, start_row_0, start_row_0)
      CALL cp_cfm_get_submatrix(cfm_eigenvectors, vector_block_1, start_row_1, start_row_1)
      IF (bs_env%eps_floquet > 0) THEN
         IF (ABS(MAXVAL(ABS(vector_block_0)) - MAXVAL(ABS(vector_block_1))) > bs_env%eps_floquet) THEN
            CPABORT("FLOQUET DIAGONALIZATION DID NOT CONVERGE. MAX_FLOQUET_INDEX TOO SMALL")
         END IF
      END IF
      DEALLOCATE (vector_block_0, vector_block_1)
      CALL timestop(handle)

   END SUBROUTINE check_floquet_conversion

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param floquet_matrix ...
!> \param a_k ...
! **************************************************************************************************
   SUBROUTINE calculate_floquet_spectral_function(qs_env, bs_env, floquet_matrix, a_k)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type)                                  :: floquet_matrix
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_floquet_spectral_function'

      COMPLEX(KIND=dp)                                   :: diag
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: g00
      INTEGER                                            :: handle, i, i_E, max_f_index, n_E, nao, &
                                                            start_row_g00
      REAL(KIND=dp)                                      :: broad, E_min, energy, energy_step, mu
      TYPE(cp_cfm_type)                                  :: g_r
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, mos=mos)

      CALL get_mo_set(mo_set=mos(1), mu=mu, nao=nao)

      broad = bs_env%broadening_floquet
      energy_step = bs_env%energy_step_floquet
      E_min = mu - bs_env%energy_window_floquet
      max_f_index = bs_env%max_floquet_index
      n_E = SIZE(a_k)
      start_row_g00 = 1 + nao*max_f_index

      ALLOCATE (g00(nao, nao))
      CALL cp_cfm_create(g_r, floquet_matrix%matrix_struct, set_zero=.TRUE.)

      DO i_E = 1, n_E
         energy = E_min + i_E*energy_step
         diag = energy + gaussi*broad/2.0_dp
         CALL cp_cfm_set_all(g_r, z_zero, diag)
         CALL cp_cfm_scale_and_add(z_one, g_r, -z_one, floquet_matrix)
         CALL cp_cfm_lu_invert(g_r)
         g00(:, :) = 0.0_dp
         CALL cp_cfm_get_submatrix(g_r, g00, start_row_g00, start_row_g00)
         a_k(i_E) = -AIMAG(SUM([(g00(i, i), i=1, nao)]))/pi
      END DO

      DEALLOCATE (g00)
      CALL cp_cfm_release(g_r)

      CALL timestop(handle)

   END SUBROUTINE calculate_floquet_spectral_function

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param ispin ...
!> \param ikp_for_file ...
!> \param xkp ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE write_quasi_energy(qs_env, bs_env, ispin, ikp_for_file, xkp, eigenvalues)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: ispin, ikp_for_file
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_quasi_energy'

      CHARACTER(LEN=default_string_length)               :: fname
      INTEGER                                            :: handle, i, j, n, nao, nkp_start, qunit
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: indices, work
      REAL(KIND=dp)                                      :: fbz_max, fbz_min, omega, qe
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: quasi_energies
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), nao=nao)

      nkp_start = bs_env%nkp_only_DOS
      omega = bs_env%floquet_omega
      fbz_max = omega/2.0_dp
      fbz_min = -omega/2.0_dp

      n = SIZE(eigenvalues)
      ALLOCATE (indices(nao), work(nao), quasi_energies(nao))
      indices(:) = [(i, i=1, nao)]
      DO i = 1, nao
         DO j = 1, n
            IF (ANY(indices == j)) CYCLE
            IF (ABS(eigenvalues(j)) < ABS(eigenvalues(indices(i)))) indices(i) = j
         END DO
      END DO
      CALL sort(indices, nao, work)

      quasi_energies(:) = eigenvalues(indices(:))

      IF (bs_env%para_env%is_source()) THEN
         WRITE (fname, "(2A)") TRIM(bs_env%floquet_qe_file), ".bs"
         IF (ikp_for_file == 1) THEN
            CALL open_file(TRIM(fname), unit_number=qunit, file_status="REPLACE", &
                           file_action="WRITE")
         ELSE
            CALL open_file(TRIM(fname), unit_number=qunit, file_status="OLD", &
                           file_action="WRITE", file_position="APPEND")
         END IF

         WRITE (qunit, "(A)") "# Quasi-energies obtained by diagonalising the Floquet Hamiltonian"
         WRITE (qunit, "(A)") "# (in units of eV, folded to Floquet Brillouin zone)"
         WRITE (qunit, "(A,I0,T10,A,I0,A,T24,3(1X,F14.8))") &
            "#  Spin ", ispin, "  Point ", ikp_for_file, ": ", xkp(1:3)
         WRITE (qunit, "(A)") "#  Floquet band  Quasi-energy [eV]"
         DO i = 1, nao
            qe = quasi_energies(i)
            WRITE (qunit, "(I8,F21.8)") i, qe*evolt
         END DO
         CALL close_file(qunit)
      END IF
      DEALLOCATE (indices, work, quasi_energies)

      CALL timestop(handle)

   END SUBROUTINE write_quasi_energy
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param bs_env ...
!> \param ispin ...
!> \param ikp_for_file ...
!> \param xkp ...
!> \param a_k ...
! **************************************************************************************************
   SUBROUTINE write_floquet_spectral_function(qs_env, bs_env, ispin, ikp_for_file, xkp, a_k)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER                                            :: ispin, ikp_for_file
      REAL(KIND=dp), DIMENSION(3)                        :: xkp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a_k

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_floquet_spectral_function'

      CHARACTER(LEN=default_string_length)               :: fname
      INTEGER                                            :: handle, i_E, n_E, nkp_start, wunit
      REAL(KIND=dp)                                      :: E_min, energy, energy_step, mu
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, mos=mos)

      CALL get_mo_set(mo_set=mos(1), mu=mu)

      energy_step = bs_env%energy_step_floquet
      E_min = mu - bs_env%energy_window_floquet
      n_E = SIZE(a_k)
      nkp_start = bs_env%nkp_only_DOS

      IF (bs_env%para_env%is_source()) THEN
         WRITE (fname, "(2A)") TRIM(bs_env%floquet_dos_file), ".out"
         IF (ikp_for_file == 1) THEN
            CALL open_file(TRIM(fname), unit_number=wunit, file_status="REPLACE", &
                           file_action="WRITE")
         ELSE
            CALL open_file(TRIM(fname), unit_number=wunit, file_status="OLD", &
                           file_action="WRITE", file_position="APPEND")
         END IF

         WRITE (wunit, "(A)") "# Floquet Density of States: D(ω,k) = -1/π*Im[Tr_KS(G^R(ω,k))]"
         WRITE (wunit, "(A,I0,T10,A,I0,A,T24,3(1X,F14.8))") &
            "#  Spin ", ispin, "  Point ", ikp_for_file, ": ", xkp(1:3)
         WRITE (wunit, "(A)") "#Energy-E_F (eV)  A(ω,k) = DOS (1/eV)"
         DO i_E = 1, n_E
            energy = E_min + i_E*energy_step
            WRITE (wunit, "(2X,2G13.4)") energy*evolt, a_k(i_E)/evolt
         END DO
         CALL close_file(wunit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_floquet_spectral_function

END MODULE floquet_utils
