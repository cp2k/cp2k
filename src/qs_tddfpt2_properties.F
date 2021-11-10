!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_properties
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Martin2003,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_solve
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_fm,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_vectorsnorm
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_filter, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE dbcsr_vector,                    ONLY: create_col_vec_from_matrix,&
                                              create_replicated_col_vec_from_matrix,&
                                              create_replicated_row_vec_from_matrix,&
                                              dbcsr_matrix_colvec_multiply
   USE input_constants,                 ONLY: tddfpt_dipole_berry,&
                                              tddfpt_dipole_length,&
                                              tddfpt_dipole_velocity
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kahan_sum,                       ONLY: accurate_dot_product
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              int_8
   USE mathconstants,                   ONLY: twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_irecv,&
                                              mp_isend,&
                                              mp_sum,&
                                              mp_wait
   USE molden_utils,                    ONLY: write_mos_molden
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              set_mo_set
   USE qs_moments,                      ONLY: build_berry_moment_matrix
   USE qs_operators_ao,                 ONLY: rRc_xyz_ao
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE string_utilities,                ONLY: integer_to_string
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_properties'

   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt_dipole_operator, tddfpt_print_summary, tddfpt_print_excitation_analysis, &
             tddfpt_print_nto_analysis

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Compute the action of the dipole operator on the ground state wave function.
!> \param dipole_op_mos_occ  2-D array [x,y,z ; spin] of matrices where to put the computed quantity
!>                           (allocated and initialised on exit)
!> \param tddfpt_control     TDDFPT control parameters
!> \param gs_mos             molecular orbitals optimised for the ground state
!> \param qs_env             Quickstep environment
!> \par History
!>    * 05.2016 created as 'tddfpt_print_summary' [Sergey Chulkov]
!>    * 06.2018 dipole operator based on the Berry-phase formula [Sergey Chulkov]
!>    * 08.2018 splited of from 'tddfpt_print_summary' and merged with code from 'tddfpt'
!>              [Sergey Chulkov]
!> \note \parblock
!>       Adapted version of the subroutine find_contributions() which was originally created
!>       by Thomas Chassaing on 02.2005.
!>
!>       The relation between dipole integrals in velocity and length forms are the following:
!>       \f[<\psi_i|\nabla|\psi_a> = <\psi_i|\vec{r}|\hat{H}\psi_a> - <\hat{H}\psi_i|\vec{r}|\psi_a>
!>                                 = (\epsilon_a - \epsilon_i) <\psi_i|\vec{r}|\psi_a> .\f],
!>       due to the commutation identity:
!>       \f[\vec{r}\hat{H} - \hat{H}\vec{r} = [\vec{r},\hat{H}] = [\vec{r},-1/2 \nabla^2] = \nabla\f] .
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_dipole_operator(dipole_op_mos_occ, tddfpt_control, gs_mos, qs_env)
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(inout)                                   :: dipole_op_mos_occ
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_dipole_operator'

      INTEGER                                            :: handle, i_cos_sin, icol, ideriv, irow, &
                                                            ispin, jderiv, nao, ncols_local, &
                                                            ndim_periodic, nrows_local, nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      REAL(kind=dp)                                      :: eval_occ
      REAL(kind=dp), DIMENSION(3)                        :: kvec, reference_point
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: local_data_ediff, local_data_wfm
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:)     :: gamma_00, gamma_inv_00
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_mos_virt
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: dBerry_mos_occ, gamma_real_imag, opvec
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: ediff_inv, rRc_mos_occ, wfm_ao_ao, &
                                                            wfm_mo_virt_mo_occ
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: berry_cossin_xyz, matrix_s, rRc_xyz
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, cell, matrix_s, pw_env)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, cell=cell, matrix_s=matrix_s, pw_env=pw_env)

      nspins = SIZE(gs_mos)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)
      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
      END DO

      ! +++ allocate dipole operator matrices (must be deallocated elsewhere)
      ALLOCATE (dipole_op_mos_occ(nderivs, nspins))
      DO ispin = 1, nspins
         CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct)

         DO ideriv = 1, nderivs
            NULLIFY (dipole_op_mos_occ(ideriv, ispin)%matrix)
            CALL cp_fm_create(dipole_op_mos_occ(ideriv, ispin)%matrix, fm_struct)
         END DO
      END DO

      ! +++ allocate work matrices
      ALLOCATE (S_mos_virt(nspins))
      DO ispin = 1, nspins
         NULLIFY (S_mos_virt(ispin)%matrix)
         CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, matrix_struct=fm_struct)
         CALL cp_fm_create(S_mos_virt(ispin)%matrix, fm_struct)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, &
                                      gs_mos(ispin)%mos_virt, &
                                      S_mos_virt(ispin)%matrix, &
                                      ncol=nmo_virt(ispin), alpha=1.0_dp, beta=0.0_dp)
      END DO

      ! check that the chosen dipole operator is consistent with the periodic boundary conditions used
      CALL pw_env_get(pw_env, poisson_env=poisson_env)
      ndim_periodic = COUNT(poisson_env%parameters%periodic == 1)

      SELECT CASE (tddfpt_control%dipole_form)
      CASE (tddfpt_dipole_berry)
         IF (ndim_periodic /= 3) THEN
            CALL cp_warn(__LOCATION__, &
                         "Fully periodic Poisson solver (PERIODIC xyz) is needed "// &
                         "for oscillator strengths based on the Berry phase formula")
         END IF

         NULLIFY (berry_cossin_xyz)
         ! index: 1 = Re[exp(-i * G_t * t)],
         !        2 = Im[exp(-i * G_t * t)];
         ! t = x,y,z
         CALL dbcsr_allocate_matrix_set(berry_cossin_xyz, 2)

         DO i_cos_sin = 1, 2
            CALL dbcsr_init_p(berry_cossin_xyz(i_cos_sin)%matrix)
            CALL dbcsr_copy(berry_cossin_xyz(i_cos_sin)%matrix, matrix_s(1)%matrix)
         END DO

         ! +++ allocate berry-phase-related work matrices
         ALLOCATE (gamma_00(nspins), gamma_inv_00(nspins), gamma_real_imag(2, nspins), opvec(2, nspins))
         ALLOCATE (dBerry_mos_occ(nderivs, nspins))
         DO ispin = 1, nspins
            NULLIFY (fm_struct)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_occ(ispin), &
                                     ncol_global=nmo_occ(ispin), context=blacs_env)

            NULLIFY (gamma_00(ispin)%matrix, gamma_inv_00(ispin)%matrix)
            CALL cp_cfm_create(gamma_00(ispin)%matrix, fm_struct)
            CALL cp_cfm_create(gamma_inv_00(ispin)%matrix, fm_struct)

            DO i_cos_sin = 1, 2
               NULLIFY (gamma_real_imag(i_cos_sin, ispin)%matrix)
               CALL cp_fm_create(gamma_real_imag(i_cos_sin, ispin)%matrix, fm_struct)
            END DO
            CALL cp_fm_struct_release(fm_struct)

            ! G_real C_0, G_imag C_0
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct)
            DO i_cos_sin = 1, 2
               NULLIFY (opvec(i_cos_sin, ispin)%matrix)
               CALL cp_fm_create(opvec(i_cos_sin, ispin)%matrix, fm_struct)
            END DO

            ! dBerry * C_0
            DO ideriv = 1, nderivs
               NULLIFY (dBerry_mos_occ(ideriv, ispin)%matrix)
               CALL cp_fm_create(dBerry_mos_occ(ideriv, ispin)%matrix, fm_struct)
               CALL cp_fm_set_all(dBerry_mos_occ(ideriv, ispin)%matrix, 0.0_dp)
            END DO
         END DO

         DO ideriv = 1, nderivs
            kvec(:) = twopi*cell%h_inv(ideriv, :)
            DO i_cos_sin = 1, 2
               CALL dbcsr_set(berry_cossin_xyz(i_cos_sin)%matrix, 0.0_dp)
            END DO
            CALL build_berry_moment_matrix(qs_env, berry_cossin_xyz(1)%matrix, &
                                           berry_cossin_xyz(2)%matrix, kvec)

            DO ispin = 1, nspins
               ! i_cos_sin = 1: cos (real) component; opvec(1) = gamma_real C_0
               ! i_cos_sin = 2: sin (imaginary) component; opvec(2) = gamma_imag C_0
               DO i_cos_sin = 1, 2
                  CALL cp_dbcsr_sm_fm_multiply(berry_cossin_xyz(i_cos_sin)%matrix, &
                                               gs_mos(ispin)%mos_occ, &
                                               opvec(i_cos_sin, ispin)%matrix, &
                                               ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
               END DO

               CALL cp_gemm('T', 'N', nmo_occ(ispin), nmo_occ(ispin), nao, &
                            1.0_dp, gs_mos(ispin)%mos_occ, opvec(1, ispin)%matrix, &
                            0.0_dp, gamma_real_imag(1, ispin)%matrix)

               CALL cp_gemm('T', 'N', nmo_occ(ispin), nmo_occ(ispin), nao, &
                            -1.0_dp, gs_mos(ispin)%mos_occ, opvec(2, ispin)%matrix, &
                            0.0_dp, gamma_real_imag(2, ispin)%matrix)

               CALL cp_fm_to_cfm(msourcer=gamma_real_imag(1, ispin)%matrix, &
                                 msourcei=gamma_real_imag(2, ispin)%matrix, &
                                 mtarget=gamma_00(ispin)%matrix)

               ! gamma_inv_00 = Q = [C_0^T (gamma_real - i gamma_imag) C_0] ^ {-1}
               CALL cp_cfm_set_all(gamma_inv_00(ispin)%matrix, z_zero, z_one)
               CALL cp_cfm_solve(gamma_00(ispin)%matrix, gamma_inv_00(ispin)%matrix)

               CALL cp_cfm_to_fm(msource=gamma_inv_00(ispin)%matrix, &
                                 mtargetr=gamma_real_imag(1, ispin)%matrix, &
                                 mtargeti=gamma_real_imag(2, ispin)%matrix)

               ! dBerry_mos_occ is identical to dBerry_psi0 from qs_linres_op % polar_operators()
               CALL cp_gemm("N", "N", nao, nmo_occ(ispin), nmo_occ(ispin), &
                            1.0_dp, opvec(1, ispin)%matrix, gamma_real_imag(2, ispin)%matrix, &
                            0.0_dp, dipole_op_mos_occ(1, ispin)%matrix)
               CALL cp_gemm("N", "N", nao, nmo_occ(ispin), nmo_occ(ispin), &
                            -1.0_dp, opvec(2, ispin)%matrix, gamma_real_imag(1, ispin)%matrix, &
                            1.0_dp, dipole_op_mos_occ(1, ispin)%matrix)

               DO jderiv = 1, nderivs
                  CALL cp_fm_scale_and_add(1.0_dp, dBerry_mos_occ(jderiv, ispin)%matrix, &
                                           cell%hmat(jderiv, ideriv), dipole_op_mos_occ(1, ispin)%matrix)
               END DO
            END DO
         END DO

         ! --- release berry-phase-related work matrices
         DO ispin = nspins, 1, -1
            DO i_cos_sin = SIZE(opvec, 1), 1, -1
               CALL cp_fm_release(opvec(i_cos_sin, ispin)%matrix)
            END DO

            DO i_cos_sin = SIZE(gamma_real_imag, 1), 1, -1
               CALL cp_fm_release(gamma_real_imag(i_cos_sin, ispin)%matrix)
            END DO

            CALL cp_cfm_release(gamma_inv_00(ispin)%matrix)
            CALL cp_cfm_release(gamma_00(ispin)%matrix)
         END DO
         DEALLOCATE (gamma_00, gamma_inv_00, gamma_real_imag, opvec)
         CALL dbcsr_deallocate_matrix_set(berry_cossin_xyz)

         NULLIFY (fm_struct, wfm_ao_ao)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
         CALL cp_fm_create(wfm_ao_ao, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         ! trans_dipole = 2|e|/|G_mu| * Tr Imag(evects^T * (gamma_real - i gamma_imag) * C_0 * gamma_inv_00) +
         !                2|e|/|G_mu| * Tr Imag(C_0^T * (gamma_real - i gamma_imag) * evects * gamma_inv_00) ,
         !
         ! Taking into account the symmetry of the matrices 'gamma_real' and 'gamma_imag' and the fact
         ! that the response wave-function is a real-valued function, the above expression can be simplified as
         ! trans_dipole = 4|e|/|G_mu| * Tr Imag(evects^T * (gamma_real - i gamma_imag) * C_0 * gamma_inv_00)
         !
         ! 1/|G_mu| = |lattice_vector_mu| / (2*pi) .
         DO ispin = 1, nspins
            ! wfm_ao_ao = S * mos_virt * mos_virt^T
            CALL cp_gemm('N', 'T', nao, nao, nmo_virt(ispin), &
                         1.0_dp/twopi, S_mos_virt(ispin)%matrix, gs_mos(ispin)%mos_virt, &
                         0.0_dp, wfm_ao_ao)

            DO ideriv = 1, nderivs
               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nao, &
                            1.0_dp, wfm_ao_ao, dBerry_mos_occ(ideriv, ispin)%matrix, &
                            0.0_dp, dipole_op_mos_occ(ideriv, ispin)%matrix)
            END DO
         END DO

         CALL cp_fm_release(wfm_ao_ao)

         DO ispin = nspins, 1, -1
            DO ideriv = SIZE(dBerry_mos_occ, 1), 1, -1
               CALL cp_fm_release(dBerry_mos_occ(ideriv, ispin)%matrix)
            END DO
         END DO
         DEALLOCATE (dBerry_mos_occ)

      CASE (tddfpt_dipole_length)
         IF (ndim_periodic /= 0) THEN
            CALL cp_warn(__LOCATION__, &
                         "Non-periodic Poisson solver (PERIODIC none) is needed "// &
                         "for oscillator strengths based on the length operator")
         END IF

         ! compute components of the dipole operator in the length form
         NULLIFY (rRc_xyz)
         CALL dbcsr_allocate_matrix_set(rRc_xyz, nderivs)

         DO ideriv = 1, nderivs
            CALL dbcsr_init_p(rRc_xyz(ideriv)%matrix)
            CALL dbcsr_copy(rRc_xyz(ideriv)%matrix, matrix_s(1)%matrix)
         END DO

         CALL get_reference_point(reference_point, qs_env=qs_env, &
                                  reference=tddfpt_control%dipole_reference, &
                                  ref_point=tddfpt_control%dipole_ref_point)

         CALL rRc_xyz_ao(op=rRc_xyz, qs_env=qs_env, rc=reference_point, order=1, &
                         minimum_image=.FALSE., soft=.FALSE.)

         NULLIFY (fm_struct, wfm_ao_ao)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)
         CALL cp_fm_create(wfm_ao_ao, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         DO ispin = 1, nspins
            NULLIFY (rRc_mos_occ)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, matrix_struct=fm_struct)
            CALL cp_fm_create(rRc_mos_occ, fm_struct)

            ! wfm_ao_ao = S * mos_virt * mos_virt^T
            CALL cp_gemm('N', 'T', nao, nao, nmo_virt(ispin), &
                         1.0_dp, S_mos_virt(ispin)%matrix, gs_mos(ispin)%mos_virt, &
                         0.0_dp, wfm_ao_ao)

            DO ideriv = 1, nderivs
               CALL cp_dbcsr_sm_fm_multiply(rRc_xyz(ideriv)%matrix, &
                                            gs_mos(ispin)%mos_occ, &
                                            rRc_mos_occ, &
                                            ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)

               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nao, &
                            1.0_dp, wfm_ao_ao, rRc_mos_occ, &
                            0.0_dp, dipole_op_mos_occ(ideriv, ispin)%matrix)
            END DO

            CALL cp_fm_release(rRc_mos_occ)
         END DO

         CALL cp_fm_release(wfm_ao_ao)
         CALL dbcsr_deallocate_matrix_set(rRc_xyz)

      CASE (tddfpt_dipole_velocity)
         IF (SIZE(matrix_s) < nderivs + 1) THEN
            CPABORT("Where are the overlap matrix derivatives?")
         END IF

         DO ispin = 1, nspins
            NULLIFY (fm_struct, ediff_inv, wfm_mo_virt_mo_occ)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), &
                                     ncol_global=nmo_occ(ispin), context=blacs_env)
            CALL cp_fm_create(ediff_inv, fm_struct)
            CALL cp_fm_create(wfm_mo_virt_mo_occ, fm_struct)
            CALL cp_fm_struct_release(fm_struct)

            CALL cp_fm_get_info(ediff_inv, nrow_local=nrows_local, ncol_local=ncols_local, &
                                row_indices=row_indices, col_indices=col_indices, local_data=local_data_ediff)
            CALL cp_fm_get_info(wfm_mo_virt_mo_occ, local_data=local_data_wfm)

!$OMP       PARALLEL DO DEFAULT(NONE), &
!$OMP                PRIVATE(eval_occ, icol, irow), &
!$OMP                SHARED(col_indices, gs_mos, ispin, local_data_ediff, ncols_local, nrows_local, row_indices)
            DO icol = 1, ncols_local
               ! E_occ_i ; imo_occ = col_indices(icol)
               eval_occ = gs_mos(ispin)%evals_occ(col_indices(icol))

               DO irow = 1, nrows_local
                  ! ediff_inv_weights(a, i) = 1.0 / (E_virt_a - E_occ_i)
                  ! imo_virt = row_indices(irow)
                  local_data_ediff(irow, icol) = 1.0_dp/(gs_mos(ispin)%evals_virt(row_indices(irow)) - eval_occ)
               END DO
            END DO
!$OMP       END PARALLEL DO

            DO ideriv = 1, nderivs
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(ideriv + 1)%matrix, &
                                            gs_mos(ispin)%mos_occ, &
                                            dipole_op_mos_occ(ideriv, ispin)%matrix, &
                                            ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)

               CALL cp_gemm('T', 'N', nmo_virt(ispin), nmo_occ(ispin), nao, &
                            1.0_dp, gs_mos(ispin)%mos_virt, dipole_op_mos_occ(ideriv, ispin)%matrix, &
                            0.0_dp, wfm_mo_virt_mo_occ)

               ! in-place element-wise (Schur) product;
               ! avoid allocation of a temporary [nmo_virt x nmo_occ] matrix which is needed
               ! for cp_fm_schur_product() subroutine call

!$OMP          PARALLEL DO DEFAULT(NONE), &
!$OMP                   PRIVATE(icol, irow), &
!$OMP                   SHARED(ispin, local_data_ediff, local_data_wfm, ncols_local, nrows_local)
               DO icol = 1, ncols_local
                  DO irow = 1, nrows_local
                     local_data_wfm(irow, icol) = local_data_wfm(irow, icol)*local_data_ediff(irow, icol)
                  END DO
               END DO
!$OMP          END PARALLEL DO

               CALL cp_gemm('N', 'N', nao, nmo_occ(ispin), nmo_virt(ispin), &
                            1.0_dp, S_mos_virt(ispin)%matrix, wfm_mo_virt_mo_occ, &
                            0.0_dp, dipole_op_mos_occ(ideriv, ispin)%matrix)
            END DO

            CALL cp_fm_release(wfm_mo_virt_mo_occ)
            CALL cp_fm_release(ediff_inv)
         END DO

      CASE DEFAULT
         CPABORT("Unimplemented form of the dipole operator")
      END SELECT

      ! --- release work matrices
      DO ispin = nspins, 1, -1
         CALL cp_fm_release(S_mos_virt(ispin)%matrix)
      END DO
      DEALLOCATE (S_mos_virt)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_dipole_operator

! **************************************************************************************************
!> \brief Print final TDDFPT excitation energies and oscillator strengths.
!> \param log_unit           output unit
!> \param evects             TDDFPT trial vectors (SIZE(evects,1) -- number of spins;
!>                           SIZE(evects,2) -- number of excited states to print)
!> \param evals              TDDFPT eigenvalues
!> \param mult               multiplicity
!> \param dipole_op_mos_occ  action of the dipole operator on the ground state wave function
!>                           [x,y,z ; spin]
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 06.2016 transition dipole moments and oscillator strengths [Sergey Chulkov]
!>    * 07.2016 spin-unpolarised electron density [Sergey Chulkov]
!>    * 08.2018 compute 'dipole_op_mos_occ' in a separate subroutine [Sergey Chulkov]
!> \note \parblock
!>       Adapted version of the subroutine find_contributions() which was originally created
!>       by Thomas Chassaing on 02.2005.
!>
!>       Transition dipole moment along direction 'd' is computed as following:
!>       \f[ t_d(spin) = Tr[evects^T dipole\_op\_mos\_occ(d, spin)] .\f]
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_print_summary(log_unit, evects, evals, mult, dipole_op_mos_occ)
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      INTEGER, INTENT(in)                                :: mult
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: dipole_op_mos_occ

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_summary'

      CHARACTER(len=1)                                   :: lsd_str
      CHARACTER(len=20)                                  :: mult_str
      INTEGER                                            :: handle, ideriv, ispin, istate, nspins, &
                                                            nstates
      REAL(kind=dp)                                      :: osc_strength
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: trans_dipoles

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)

      IF (nspins > 1) THEN
         lsd_str = 'U'
      ELSE
         lsd_str = 'R'
      END IF

      ! *** summary header ***
      IF (log_unit > 0) THEN
         CALL integer_to_string(mult, mult_str)
         WRITE (log_unit, '(/,1X,A1,A,1X,A,/)') lsd_str, "-TDDFPT states of multiplicity", TRIM(mult_str)

         WRITE (log_unit, '(T10,A,T19,A,T37,A,T69,A)') "State", "Excitation", &
            "Transition dipole (a.u.)", "Oscillator"
         WRITE (log_unit, '(T10,A,T19,A,T37,A,T49,A,T61,A,T67,A)') "number", "energy (eV)", &
            "x", "y", "z", "strength (a.u.)"
         WRITE (log_unit, '(T10,72("-"))')
      END IF

      ! transition dipole moment
      ALLOCATE (trans_dipoles(nstates, nderivs, nspins))
      trans_dipoles(:, :, :) = 0.0_dp

      ! nspins == 1 .AND. mult == 3 : spin-flip transitions are forbidden due to symmetry reasons
      IF (nspins > 1 .OR. mult == 1) THEN
         DO ispin = 1, nspins
            DO ideriv = 1, nderivs
               CALL cp_fm_trace(evects(ispin, :), dipole_op_mos_occ(ideriv, ispin)%matrix, &
                                trans_dipoles(:, ideriv, ispin))
            END DO
         END DO

         IF (nspins == 1) THEN
            trans_dipoles(:, :, 1) = SQRT(2.0_dp)*trans_dipoles(:, :, 1)
         ELSE
            trans_dipoles(:, :, 1) = SQRT(trans_dipoles(:, :, 1)**2 + trans_dipoles(:, :, 2)**2)
         END IF
      END IF

      ! *** summary information ***
      DO istate = 1, nstates
         IF (log_unit > 0) THEN
            osc_strength = 2.0_dp/3.0_dp*evals(istate)* &
                           accurate_dot_product(trans_dipoles(istate, :, 1), trans_dipoles(istate, :, 1))

            WRITE (log_unit, '(1X,A,T9,I7,T19,F11.5,T31,3(1X,ES11.4E2),T69,ES12.5E2)') &
               "TDDFPT|", istate, evals(istate)*evolt, trans_dipoles(istate, 1:nderivs, 1), osc_strength
         END IF
      END DO

      ! punch a checksum for the regs
      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A,E14.6)') 'TDDFPT : CheckSum  =', SQRT(SUM(evals**2))
      END IF

      DEALLOCATE (trans_dipoles)

      CALL timestop(handle)
   END SUBROUTINE tddfpt_print_summary

! **************************************************************************************************
!> \brief Print excitation analysis.
!> \param log_unit           output unit
!> \param evects             TDDFPT trial vectors (SIZE(evects,1) -- number of spins;
!>                           SIZE(evects,2) -- number of excited states to print)
!> \param evals              TDDFPT eigenvalues
!> \param gs_mos             molecular orbitals optimised for the ground state
!> \param matrix_s           overlap matrix
!> \param min_amplitude      the smallest excitation amplitude to print
!> \par History
!>    * 05.2016 created as 'tddfpt_print_summary' [Sergey Chulkov]
!>    * 08.2018 splited of from 'tddfpt_print_summary' [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_print_excitation_analysis(log_unit, evects, evals, gs_mos, matrix_s, min_amplitude)
      INTEGER, INTENT(in)                                :: log_unit
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      REAL(kind=dp), INTENT(in)                          :: min_amplitude

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_excitation_analysis'

      CHARACTER(len=5)                                   :: spin_label
      INTEGER :: handle, icol, iproc, irow, ispin, istate, nao, ncols_local, nrows_local, nspins, &
         nstates, send_handler, send_handler2, state_spin
      INTEGER(kind=int_8)                                :: iexc, imo_occ, imo_virt, ind, nexcs, &
                                                            nexcs_local, nexcs_max_local, &
                                                            nmo_virt_occ, nmo_virt_occ_alpha
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:)     :: inds_local, inds_recv, nexcs_recv
      INTEGER(kind=int_8), DIMENSION(1)                  :: nexcs_send
      INTEGER(kind=int_8), DIMENSION(maxspins)           :: nmo_occ8, nmo_virt8
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inds, recv_handlers, recv_handlers2
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      LOGICAL                                            :: do_exc_analysis
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: weights_local, weights_neg_abs_recv, &
                                                            weights_recv
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: local_data
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: S_mos_virt, weights_fm
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)
      do_exc_analysis = min_amplitude < 1.0_dp

      CALL cp_fm_get_info(gs_mos(1)%mos_occ, context=blacs_env, para_env=para_env)
      CALL dbcsr_get_info(matrix_s, nfullrows_total=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
         nmo_occ8(ispin) = SIZE(gs_mos(ispin)%evals_occ, kind=int_8)
         nmo_virt8(ispin) = SIZE(gs_mos(ispin)%evals_virt, kind=int_8)
      END DO

      ! *** excitation analysis ***
      IF (do_exc_analysis) THEN
         CPASSERT(log_unit <= 0 .OR. para_env%mepos == para_env%source)
         nmo_virt_occ_alpha = INT(nmo_virt(1), int_8)*INT(nmo_occ(1), int_8)

         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                            Excitation analysis                              -", &
               "-------------------------------------------------------------------------------"
            WRITE (log_unit, '(8X,A,T27,A,T49,A,T69,A)') "State", "Occupied", "Virtual", "Excitation"
            WRITE (log_unit, '(8X,A,T28,A,T49,A,T69,A)') "number", "orbital", "orbital", "amplitude"
            WRITE (log_unit, '(1X,79("-"))')

            IF (nspins == 1) THEN
               state_spin = 1
               spin_label = '     '
            END IF
         END IF

         ALLOCATE (S_mos_virt(nspins), weights_fm(nspins))
         DO ispin = 1, nspins
            NULLIFY (S_mos_virt(ispin)%matrix)
            CALL cp_fm_get_info(gs_mos(ispin)%mos_virt, matrix_struct=fm_struct)
            CALL cp_fm_create(S_mos_virt(ispin)%matrix, fm_struct)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, &
                                         gs_mos(ispin)%mos_virt, &
                                         S_mos_virt(ispin)%matrix, &
                                         ncol=nmo_virt(ispin), alpha=1.0_dp, beta=0.0_dp)

            NULLIFY (fm_struct, weights_fm(ispin)%matrix)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin), &
                                     context=blacs_env)
            CALL cp_fm_create(weights_fm(ispin)%matrix, fm_struct)
            CALL cp_fm_struct_release(fm_struct)
         END DO

         nexcs_max_local = 0
         DO ispin = 1, nspins
            CALL cp_fm_get_info(weights_fm(ispin)%matrix, nrow_local=nrows_local, ncol_local=ncols_local)
            nexcs_max_local = nexcs_max_local + INT(nrows_local, int_8)*INT(ncols_local, int_8)
         END DO

         ALLOCATE (weights_local(nexcs_max_local), inds_local(nexcs_max_local))

         DO istate = 1, nstates
            nexcs_local = 0
            nmo_virt_occ = 0

            ! analyse matrix elements locally and transfer only significant
            ! excitations to the master node for subsequent ordering
            DO ispin = 1, nspins
               ! compute excitation amplitudes
               CALL cp_gemm('T', 'N', nmo_virt(ispin), nmo_occ(ispin), nao, 1.0_dp, S_mos_virt(ispin)%matrix, &
                            evects(ispin, istate)%matrix, 0.0_dp, weights_fm(ispin)%matrix)

               CALL cp_fm_get_info(weights_fm(ispin)%matrix, nrow_local=nrows_local, ncol_local=ncols_local, &
                                   row_indices=row_indices, col_indices=col_indices, local_data=local_data)

               ! locate single excitations with significant amplitudes (>= min_amplitude)
               DO icol = 1, ncols_local
                  DO irow = 1, nrows_local
                     IF (ABS(local_data(irow, icol)) >= min_amplitude) THEN
                        ! number of non-negligible excitations
                        nexcs_local = nexcs_local + 1
                        ! excitation amplitude
                        weights_local(nexcs_local) = local_data(irow, icol)
                        ! index of single excitation (ivirt, iocc, ispin) in compressed form
                        inds_local(nexcs_local) = nmo_virt_occ + INT(row_indices(irow), int_8) + &
                                                  INT(col_indices(icol) - 1, int_8)*nmo_virt8(ispin)
                     END IF
                  END DO
               END DO

               nmo_virt_occ = nmo_virt_occ + nmo_virt8(ispin)*nmo_occ8(ispin)
            END DO

            IF (para_env%mepos == para_env%source) THEN
               ! master node
               ALLOCATE (nexcs_recv(para_env%num_pe), recv_handlers(para_env%num_pe), recv_handlers2(para_env%num_pe))

               ! collect number of non-negligible excitations from other nodes
               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos) THEN
                     CALL mp_irecv(nexcs_recv(iproc:iproc), iproc - 1, para_env%group, recv_handlers(iproc), 0)
                  ELSE
                     nexcs_recv(iproc) = nexcs_local
                  END IF
               END DO

               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos) &
                     CALL mp_wait(recv_handlers(iproc))
               END DO

               ! compute total number of non-negligible excitations
               nexcs = 0
               DO iproc = 1, para_env%num_pe
                  nexcs = nexcs + nexcs_recv(iproc)
               END DO

               ! receive indices and amplitudes of selected excitations
               ALLOCATE (weights_recv(nexcs), weights_neg_abs_recv(nexcs))
               ALLOCATE (inds_recv(nexcs), inds(nexcs))

               nmo_virt_occ = 0
               DO iproc = 1, para_env%num_pe
                  IF (nexcs_recv(iproc) > 0) THEN
                     IF (iproc - 1 /= para_env%mepos) THEN
                        ! excitation amplitudes
                        CALL mp_irecv(weights_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)), &
                                      iproc - 1, para_env%group, recv_handlers(iproc), 1)
                        ! compressed indices
                        CALL mp_irecv(inds_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)), &
                                      iproc - 1, para_env%group, recv_handlers2(iproc), 2)
                     ELSE
                        ! data on master node
                        weights_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)) = weights_local(1:nexcs_recv(iproc))
                        inds_recv(nmo_virt_occ + 1:nmo_virt_occ + nexcs_recv(iproc)) = inds_local(1:nexcs_recv(iproc))
                     END IF

                     nmo_virt_occ = nmo_virt_occ + nexcs_recv(iproc)
                  END IF
               END DO

               DO iproc = 1, para_env%num_pe
                  IF (iproc - 1 /= para_env%mepos .AND. nexcs_recv(iproc) > 0) THEN
                     CALL mp_wait(recv_handlers(iproc))
                     CALL mp_wait(recv_handlers2(iproc))
                  END IF
               END DO

               DEALLOCATE (nexcs_recv, recv_handlers, recv_handlers2)
            ELSE
               ! working node: send the number of selected excited states to the master node
               nexcs_send(1) = nexcs_local
               CALL mp_isend(nexcs_send, para_env%source, para_env%group, send_handler, 0)
               CALL mp_wait(send_handler)

               IF (nexcs_local > 0) THEN
                  ! send excitation amplitudes
                  CALL mp_isend(weights_local(1:nexcs_local), para_env%source, para_env%group, send_handler, 1)
                  ! send compressed indices
                  CALL mp_isend(inds_local(1:nexcs_local), para_env%source, para_env%group, send_handler2, 2)

                  CALL mp_wait(send_handler)
                  CALL mp_wait(send_handler2)
               END IF
            END IF

            ! sort non-negligible excitations on the master node according to their amplitudes,
            ! uncompress indices and print summary information
            IF (para_env%mepos == para_env%source .AND. log_unit > 0) THEN
               weights_neg_abs_recv(:) = -ABS(weights_recv)
               CALL sort(weights_neg_abs_recv, INT(nexcs), inds)

               WRITE (log_unit, '(T7,I8,F10.5,A)') istate, evals(istate)*evolt, " eV"

               DO iexc = 1, nexcs
                  ind = inds_recv(inds(iexc)) - 1
                  IF (nspins > 1) THEN
                     IF (ind < nmo_virt_occ_alpha) THEN
                        state_spin = 1
                        spin_label = '(alp)'
                     ELSE
                        state_spin = 2
                        ind = ind - nmo_virt_occ_alpha
                        spin_label = '(bet)'
                     END IF
                  END IF

                  imo_occ = ind/nmo_virt8(state_spin) + 1
                  imo_virt = MOD(ind, nmo_virt8(state_spin)) + 1

                  WRITE (log_unit, '(T27,I8,1X,A5,T48,I8,1X,A5,T70,F9.6)') imo_occ, spin_label, &
                     nmo_occ8(state_spin) + imo_virt, spin_label, weights_recv(inds(iexc))
               END DO
            END IF

            ! deallocate temporary arrays
            IF (para_env%mepos == para_env%source) &
               DEALLOCATE (weights_recv, weights_neg_abs_recv, inds_recv, inds)
         END DO

         DEALLOCATE (weights_local, inds_local)
         IF (log_unit > 0) THEN
            WRITE (log_unit, "(1X,A)") &
               "-------------------------------------------------------------------------------"
         END IF
      END IF

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(weights_fm(ispin)%matrix)
         CALL cp_fm_release(S_mos_virt(ispin)%matrix)
      END DO
      DEALLOCATE (S_mos_virt, weights_fm)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_print_excitation_analysis

! **************************************************************************************************
!> \brief Print natural transition orbital analysis.
!> \param qs_env             Information on Kinds and Particles
!> \param evects             TDDFPT trial vectors (SIZE(evects,1) -- number of spins;
!>                           SIZE(evects,2) -- number of excited states to print)
!> \param evals              TDDFPT eigenvalues
!> \param gs_mos             molecular orbitals optimised for the ground state
!> \param matrix_s           overlap matrix
!> \param print_section      ...
!> \par History
!>    * 06.2019 created [JGH]
! **************************************************************************************************
   SUBROUTINE tddfpt_print_nto_analysis(qs_env, evects, evals, gs_mos, matrix_s, print_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_print_nto_analysis'
      INTEGER, PARAMETER                                 :: ntomax = 10

      CHARACTER(LEN=20), DIMENSION(2)                    :: nto_name
      INTEGER                                            :: handle, i, ia, icg, iounit, ispin, &
                                                            istate, nao, nmax, nmo, nnto, nspins, &
                                                            nstates
      INTEGER, DIMENSION(2)                              :: iv
      INTEGER, DIMENSION(2, ntomax)                      :: ia_index
      INTEGER, DIMENSION(:), POINTER                     :: stride
      LOGICAL                                            :: append_cube, cube_file
      REAL(KIND=dp)                                      :: alpha, beta, dkz, rkrk, rr, sume, &
                                                            threshold
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvalues
      REAL(KIND=dp), DIMENSION(ntomax)                   :: ia_eval
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: teig
      TYPE(cp_fm_struct_type), POINTER                   :: fm_mo_struct, fm_struct
      TYPE(cp_fm_type), POINTER                          :: ev, Sev, tmat, wvec
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: dvec, pdens, rep_col_vec, rep_row_vec, &
                                                            rvec, tvec, xvec, zvec
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: nto_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: molden_section, nto_section

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_section, &
                                           "NTO_ANALYSIS"), cp_p_file)) THEN

         CALL cite_reference(Martin2003)

         CALL section_vals_val_get(print_section, "NTO_ANALYSIS%THRESHOLD", r_val=threshold)
         IF (iounit > 0) THEN
            WRITE (iounit, "(1X,A)") "", &
               "-------------------------------------------------------------------------------", &
               "-                            Natural Orbital analysis                         -", &
               "-------------------------------------------------------------------------------"
         END IF

         nspins = SIZE(evects, 1)
         nstates = SIZE(evects, 2)
         CALL dbcsr_get_info(matrix_s, nfullrows_total=nao)

         DO istate = 1, nstates
            IF (iounit > 0) THEN
               WRITE (iounit, "(1X,A,I3,T30,F10.5,A)") "  STATE NR. ", istate, evals(istate)*evolt, " eV"
            END IF
            nmax = 0
            DO ispin = 1, nspins
               ev => evects(ispin, istate)%matrix
               CALL cp_fm_get_info(ev, matrix_struct=fm_struct, ncol_global=nmo)
               nmax = MAX(nmax, nmo)
            END DO
            ALLOCATE (eigenvalues(nmax, nspins))
            eigenvalues = 0.0_dp
            ! SET 1: Hole states
            ! SET 2: Particle states
            nto_name(1) = 'Hole_states'
            nto_name(2) = 'Particle_states'
            ALLOCATE (nto_set(2))
            DO i = 1, 2
               NULLIFY (nto_set(i)%mo_set)
               CALL allocate_mo_set(nto_set(i)%mo_set, nao, ntomax, 0, 0.0_dp, 1.0_dp, 0.0_dp)
               ev => evects(1, istate)%matrix
               CALL cp_fm_get_info(ev, matrix_struct=fm_struct)
               CALL cp_fm_struct_create(fmstruct=fm_mo_struct, template_fmstruct=fm_struct, &
                                        ncol_global=ntomax)
               CALL cp_fm_create(tmat, fm_mo_struct)
               CALL init_mo_set(nto_set(i)%mo_set, fm_ref=tmat, name=nto_name(i))
               CALL cp_fm_release(tmat)
               CALL cp_fm_struct_release(fm_mo_struct)
            END DO
            !
            ALLOCATE (teig(nspins))
            DO ispin = 1, nspins
               NULLIFY (teig(ispin)%matrix)
            END DO
            ! hole states
            DO ispin = 1, nspins
               ev => evects(ispin, istate)%matrix
               CALL cp_fm_get_info(ev, matrix_struct=fm_struct, ncol_global=nmo)
               CALL cp_fm_create(Sev, fm_struct)
               CALL cp_fm_struct_create(fmstruct=fm_mo_struct, template_fmstruct=fm_struct, &
                                        nrow_global=nmo, ncol_global=nmo)
               CALL cp_fm_create(tmat, fm_mo_struct)
               CALL cp_fm_create(teig(ispin)%matrix, fm_mo_struct)
               CALL cp_dbcsr_sm_fm_multiply(matrix_s, ev, Sev, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
               CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, ev, Sev, 0.0_dp, tmat)

               CALL choose_eigv_solver(tmat, teig(ispin)%matrix, eigenvalues(1:nmo, ispin))

               CALL cp_fm_struct_release(fm_mo_struct)
               CALL cp_fm_release(tmat)
               CALL cp_fm_release(Sev)
            END DO
            ! find major determinants i->a
            ia_index = 0
            sume = 0.0_dp
            nnto = 0
            DO i = 1, ntomax
               iv = MAXLOC(eigenvalues)
               ia_eval(i) = eigenvalues(iv(1), iv(2))
               ia_index(1:2, i) = iv(1:2)
               sume = sume + ia_eval(i)
               eigenvalues(iv(1), iv(2)) = 0.0_dp
               nnto = nnto + 1
               IF (sume > threshold) EXIT
            END DO
            ! store hole states
            CALL set_mo_set(nto_set(1)%mo_set, nmo=nnto)
            DO i = 1, nnto
               ia = ia_index(1, i)
               ispin = ia_index(2, i)
               CALL cp_fm_get_info(gs_mos(ispin)%mos_occ, ncol_global=nmo)
               CALL cp_fm_get_info(teig(ispin)%matrix, matrix_struct=fm_struct)
               CALL cp_fm_struct_create(fmstruct=fm_mo_struct, template_fmstruct=fm_struct, &
                                        nrow_global=nmo, ncol_global=1)
               CALL cp_fm_create(tmat, fm_mo_struct)
               CALL cp_fm_struct_release(fm_mo_struct)
               CALL cp_fm_get_info(gs_mos(1)%mos_occ, matrix_struct=fm_struct)
               CALL cp_fm_struct_create(fmstruct=fm_mo_struct, template_fmstruct=fm_struct, &
                                        ncol_global=1)
               CALL cp_fm_create(wvec, fm_mo_struct)
               CALL cp_fm_struct_release(fm_mo_struct)
               CALL cp_fm_to_fm(teig(ispin)%matrix, tmat, 1, ia, 1)
               CALL cp_gemm('N', 'N', nao, 1, nmo, 1.0_dp, gs_mos(ispin)%mos_occ, &
                            tmat, 0.0_dp, wvec)
               CALL cp_fm_to_fm(wvec, nto_set(1)%mo_set%mo_coeff, 1, 1, i)
               CALL cp_fm_release(wvec)
               CALL cp_fm_release(tmat)
            END DO
            ! particle states
            CALL set_mo_set(nto_set(2)%mo_set, nmo=nnto)
            DO i = 1, nnto
               ispin = ia_index(2, i)
               ev => evects(ispin, istate)%matrix
               CALL cp_fm_get_info(ev, matrix_struct=fm_struct, ncol_global=nmo)
               CALL cp_fm_create(Sev, fm_struct)
               CALL cp_dbcsr_sm_fm_multiply(matrix_s, ev, Sev, ncol=nmo, alpha=1.0_dp, beta=0.0_dp)
               CALL dbcsr_create(pdens, template=matrix_s, matrix_type=dbcsr_type_no_symmetry)
               ! P(nu,mu) = SUM_j XT(nu,j)*CT(mu,j)
               CALL cp_dbcsr_plus_fm_fm_t(pdens, ev, Sev, nmo, 1.0_dp, keep_sparsity=.FALSE.)
               CALL dbcsr_filter(pdens, 1.0e-6_dp)
               CALL dbcsr_add_on_diag(pdens, -ia_eval(i))
               ! Conjugate Gradient search for particle states
               CALL create_col_vec_from_matrix(rvec, pdens, 1)
               CALL create_col_vec_from_matrix(dvec, pdens, 1)
               CALL create_col_vec_from_matrix(zvec, pdens, 1)
               CALL create_col_vec_from_matrix(tvec, pdens, 1)
               CALL create_col_vec_from_matrix(xvec, pdens, 1)
               CALL create_replicated_col_vec_from_matrix(rep_col_vec, pdens, 1)
               CALL create_replicated_row_vec_from_matrix(rep_row_vec, pdens, 1)
               !
               CALL dbcsr_set(rvec, 1.0_dp)
               CALL project_vector(rvec, xvec, gs_mos(ispin)%mos_occ, matrix_s)
               CALL dbcsr_matrix_colvec_multiply(pdens, xvec, rvec, -1.0_dp, &
                                                 0.0_dp, rep_row_vec, rep_col_vec)
               CALL dbcsr_copy(dvec, rvec)
               CALL vec_product(rvec, rvec, rkrk)
               DO icg = 1, 20
                  CALL dbcsr_matrix_colvec_multiply(pdens, dvec, tvec, 1.0_dp, &
                                                    0.0_dp, rep_row_vec, rep_col_vec)
                  CALL project_vector(tvec, zvec, gs_mos(ispin)%mos_occ, matrix_s)
                  CALL vec_product(dvec, zvec, dkz)
                  alpha = rkrk/dkz
                  CALL dbcsr_add(xvec, dvec, 1.0_dp, alpha)
                  CALL dbcsr_add(rvec, zvec, 1.0_dp, -alpha)
                  rr = rkrk
                  CALL vec_product(rvec, rvec, rkrk)
                  IF (rkrk < 1.E-8_dp) EXIT
                  beta = rkrk/rr
                  CALL dbcsr_add(dvec, rvec, beta, 1.0_dp)
               END DO
               !
               CALL cp_fm_get_info(gs_mos(1)%mos_occ, matrix_struct=fm_struct)
               CALL cp_fm_struct_create(fmstruct=fm_mo_struct, template_fmstruct=fm_struct, &
                                        nrow_global=nao, ncol_global=1)
               CALL cp_fm_create(wvec, fm_mo_struct)
               CALL cp_fm_struct_release(fm_mo_struct)
               CALL copy_dbcsr_to_fm(xvec, wvec)
               CALL cp_fm_to_fm(wvec, nto_set(2)%mo_set%mo_coeff, 1, 1, i)
               CALL cp_fm_release(wvec)
               !
               CALL dbcsr_release(rep_col_vec)
               CALL dbcsr_release(rep_row_vec)
               CALL dbcsr_release(rvec)
               CALL dbcsr_release(dvec)
               CALL dbcsr_release(zvec)
               CALL dbcsr_release(tvec)
               CALL dbcsr_release(xvec)

               CALL cp_fm_release(Sev)
               CALL dbcsr_release(pdens)
            END DO
            ! print
            IF (iounit > 0) THEN
               sume = 0.0_dp
               DO i = 1, nnto
                  sume = sume + ia_eval(i)
                  WRITE (iounit, "(T6,A,i2,T30,A,i1,T42,A,F8.5,T63,A,F8.5)") &
                     "Particle-Hole state:", i, " Spin:", ia_index(2, i), &
                     "Eigenvalue:", ia_eval(i), " Sum Eigv:", sume
               END DO
            END IF
            ! Cube and Molden files
            nto_section => section_vals_get_subs_vals(print_section, "NTO_ANALYSIS")
            CALL section_vals_val_get(nto_section, "CUBE_FILES", l_val=cube_file)
            CALL section_vals_val_get(nto_section, "STRIDE", i_vals=stride)
            CALL section_vals_val_get(nto_section, "APPEND", l_val=append_cube)
            IF (cube_file) THEN
               CALL print_nto_cubes(qs_env, nto_set, istate, stride, append_cube, nto_section)
            END IF
            CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set)
            molden_section => section_vals_get_subs_vals(print_section, "MOS_MOLDEN")
            CALL write_mos_molden(nto_set, qs_kind_set, particle_set, molden_section)
            !
            DEALLOCATE (eigenvalues)
            DO ispin = 1, nspins
               CALL cp_fm_release(teig(ispin)%matrix)
            END DO
            DEALLOCATE (teig)
            !
            DO i = 1, 2
               CALL deallocate_mo_set(nto_set(i)%mo_set)
            END DO
            DEALLOCATE (nto_set)
         END DO

         IF (iounit > 0) THEN
            WRITE (iounit, "(1X,A)") &
               "-------------------------------------------------------------------------------"
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_print_nto_analysis

! **************************************************************************************************
!> \brief ...
!> \param vin ...
!> \param vout ...
!> \param mos_occ ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE project_vector(vin, vout, mos_occ, matrix_s)
      TYPE(dbcsr_type)                                   :: vin, vout
      TYPE(cp_fm_type), POINTER                          :: mos_occ
      TYPE(dbcsr_type), POINTER                          :: matrix_s

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'project_vector'

      INTEGER                                            :: handle, nao, nmo
      REAL(KIND=dp)                                      :: norm(1)
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_vec_struct
      TYPE(cp_fm_type), POINTER                          :: csvec, svec, vec

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(mos_occ, matrix_struct=fm_struct, nrow_global=nao, ncol_global=nmo)
      CALL cp_fm_struct_create(fmstruct=fm_vec_struct, template_fmstruct=fm_struct, &
                               nrow_global=nao, ncol_global=1)
      CALL cp_fm_create(vec, fm_vec_struct)
      CALL cp_fm_create(svec, fm_vec_struct)
      CALL cp_fm_struct_release(fm_vec_struct)
      CALL cp_fm_struct_create(fmstruct=fm_vec_struct, template_fmstruct=fm_struct, &
                               nrow_global=nmo, ncol_global=1)
      CALL cp_fm_create(csvec, fm_vec_struct)
      CALL cp_fm_struct_release(fm_vec_struct)

      CALL copy_dbcsr_to_fm(vin, vec)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, vec, svec, ncol=1, alpha=1.0_dp, beta=0.0_dp)
      CALL cp_gemm('T', 'N', nmo, 1, nao, 1.0_dp, mos_occ, svec, 0.0_dp, csvec)
      CALL cp_gemm('N', 'N', nao, 1, nmo, -1.0_dp, mos_occ, csvec, 1.0_dp, vec)
      CALL cp_fm_vectorsnorm(vec, norm)
      CPASSERT(norm(1) > 1.e-14_dp)
      norm(1) = SQRT(1._dp/norm(1))
      CALL cp_fm_scale(norm(1), vec)
      CALL copy_fm_to_dbcsr(vec, vout, keep_sparsity=.FALSE.)

      CALL cp_fm_release(csvec)
      CALL cp_fm_release(svec)
      CALL cp_fm_release(vec)

      CALL timestop(handle)

   END SUBROUTINE project_vector

! **************************************************************************************************
!> \brief ...
!> \param va ...
!> \param vb ...
!> \param res ...
! **************************************************************************************************
   SUBROUTINE vec_product(va, vb, res)
      TYPE(dbcsr_type)                                   :: va, vb
      REAL(KIND=dp), INTENT(OUT)                         :: res

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'vec_product'

      INTEGER                                            :: blk, group, handle, icol, irow
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vba, vbb
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      res = 0.0_dp

      CALL dbcsr_get_info(va, group=group)
      CALL dbcsr_iterator_start(iter, va)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, vba, blk)
         CALL dbcsr_get_block_p(vb, row=irow, col=icol, block=vbb, found=found)
         res = res + SUM(vba*vbb)
         CPASSERT(found)
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL mp_sum(res, group)

      CALL timestop(handle)

   END SUBROUTINE vec_product

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mos ...
!> \param istate ...
!> \param stride ...
!> \param append_cube ...
!> \param print_section ...
! **************************************************************************************************
   SUBROUTINE print_nto_cubes(qs_env, mos, istate, stride, append_cube, print_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      INTEGER, INTENT(IN)                                :: istate
      INTEGER, DIMENSION(:), POINTER                     :: stride
      LOGICAL, INTENT(IN)                                :: append_cube
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube, title
      INTEGER                                            :: i, iset, nmo, unit_nr
      LOGICAL                                            :: mpi_io
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys

      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      my_pos_cube = "REWIND"
      IF (append_cube) THEN
         my_pos_cube = "APPEND"
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      particle_set=particle_set)

      DO iset = 1, 2
         CALL get_mo_set(mo_set=mos(iset)%mo_set, mo_coeff=mo_coeff, nmo=nmo)
         DO i = 1, nmo
            CALL calculate_wavefunction(mo_coeff, i, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                        cell, dft_control, particle_set, pw_env)
            IF (iset == 1) THEN
               WRITE (filename, '(a4,I3.3,I2.2,a11)') "NTO_STATE", istate, i, "_Hole_State"
            ELSEIF (iset == 2) THEN
               WRITE (filename, '(a4,I3.3,I2.2,a15)') "NTO_STATE", istate, i, "_Particle_State"
            END IF
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, print_section, '', extension=".cube", &
                                           middle_name=TRIM(filename), file_position=my_pos_cube, &
                                           log_filename=.FALSE., ignore_should_output=.TRUE., mpi_io=mpi_io)
            IF (iset == 1) THEN
               WRITE (title, *) "Natural Transition Orbital Hole State", i
            ELSEIF (iset == 2) THEN
               WRITE (title, *) "Natural Transition Orbital Particle State", i
            END IF
            CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, stride=stride, mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, print_section, '', &
                                              ignore_should_output=.TRUE., mpi_io=mpi_io)
         END DO
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)

   END SUBROUTINE print_nto_cubes

END MODULE qs_tddfpt2_properties
