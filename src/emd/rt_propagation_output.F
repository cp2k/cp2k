!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routine for the real time propagation output.
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

MODULE rt_propagation_output
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_binary_write, dbcsr_convert_offsets_to_sizes, dbcsr_copy, dbcsr_create, &
        dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_distribution_type, dbcsr_filter, &
        dbcsr_get_info, dbcsr_get_occupation, dbcsr_init_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_checksum,&
                                              dbcsr_trace
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_double,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_iter_string,&
                                              cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_printkey_is_on
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE efield_utils,                    ONLY: make_field
   USE greenx_interface,                ONLY: greenx_refine_ft
   USE input_constants,                 ONLY: ehrenfest,&
                                              real_time_propagation
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_comm_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt,&
                                              femtoseconds
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_linres_current,               ONLY: calculate_jrho_resp
   USE qs_linres_types,                 ONLY: current_env_type
   USE qs_mo_io,                        ONLY: write_rt_mos_to_restart
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_operators_ao,                 ONLY: build_lin_mom_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_post_gpw,                 ONLY: qs_scf_post_moments,&
                                              write_mo_dependent_results,&
                                              write_mo_free_results
   USE qs_scf_post_tb,                  ONLY: scf_post_calculation_tb
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE rt_projection_mo_utils,          ONLY: compute_and_write_proj_mo
   USE rt_propagation_ft,               ONLY: multi_fft
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
   USE rt_propagation_utils,            ONLY: write_rtp_mo_cubes,&
                                              write_rtp_mos_to_output_unit
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_output'

   PUBLIC :: rt_prop_output, &
             rt_convergence, &
             rt_convergence_density, &
             report_density_occupation, &
             print_moments, &
             calc_local_moment, &
             print_ft

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param run_type ...
!> \param delta_iter ...
!> \param used_time ...
! **************************************************************************************************
   SUBROUTINE rt_prop_output(qs_env, run_type, delta_iter, used_time)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(in)                                :: run_type
      REAL(dp), INTENT(in), OPTIONAL                     :: delta_iter, used_time

      INTEGER                                            :: i, n_electrons, n_proj, natom, nspin, &
                                                            output_unit, spin, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      LOGICAL                                            :: new_file
      REAL(dp)                                           :: orthonormality, strace, tot_rho_r, trace
      REAL(kind=dp), DIMENSION(3)                        :: field, reference_point
      REAL(KIND=dp), DIMENSION(:), POINTER               :: qs_tot_rho_r
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: j_int
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, P_im, p_xyz, rho_new
      TYPE(dbcsr_type), POINTER                          :: tmp_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(section_vals_type), POINTER                   :: dft_section, input, rtp_section

      NULLIFY (logger, dft_control)

      logger => cp_get_default_logger()
      CALL get_qs_env(qs_env, &
                      rtp=rtp, &
                      matrix_s=matrix_s, &
                      input=input, &
                      rho=rho, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, sab_all=sab_all, sab_orb=sab_orb, &
                      dbcsr_dist=dbcsr_dist, nelectron_spin=nelectron_spin)

      rtp_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION")

      CALL get_qs_kind_set(qs_kind_set, nelectron=n_electrons)
      n_electrons = n_electrons - dft_control%charge

      CALL qs_rho_get(rho_struct=rho, tot_rho_r=qs_tot_rho_r)

      tot_rho_r = accurate_sum(qs_tot_rho_r)

      output_unit = cp_print_key_unit_nr(logger, rtp_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".scfLog")

      IF (output_unit > 0) THEN
         WRITE (output_unit, FMT="(/,(T3,A,T40,I5))") &
            "Information at iteration step:", rtp%iter
         WRITE (UNIT=output_unit, FMT="((T3,A,T41,2F20.10))") &
            "Total electronic density (r-space): ", &
            tot_rho_r, &
            tot_rho_r + &
            REAL(n_electrons, dp)
         WRITE (UNIT=output_unit, FMT="((T3,A,T59,F22.14))") &
            "Total energy:", rtp%energy_new
         IF (run_type == ehrenfest) &
            WRITE (UNIT=output_unit, FMT="((T3,A,T61,F20.14))") &
            "Energy difference to previous iteration step:", rtp%energy_new - rtp%energy_old
         IF (run_type == real_time_propagation) &
            WRITE (UNIT=output_unit, FMT="((T3,A,T61,F20.14))") &
            "Energy difference to initial state:", rtp%energy_new - rtp%energy_old
         IF (PRESENT(delta_iter)) &
            WRITE (UNIT=output_unit, FMT="((T3,A,T61,E20.6))") &
            "Convergence:", delta_iter
         IF (rtp%converged) THEN
            IF (run_type == real_time_propagation) &
               WRITE (UNIT=output_unit, FMT="((T3,A,T61,F12.2))") &
               "Time needed for propagation:", used_time
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,3X,F16.14))") &
               "CONVERGENCE REACHED", rtp%energy_new - rtp%energy_old
         END IF
      END IF

      IF (rtp%converged) THEN
         IF (.NOT. rtp%linear_scaling) THEN
            CALL get_rtp(rtp=rtp, mos_new=mos_new)
            CALL rt_calculate_orthonormality(orthonormality, &
                                             mos_new, matrix_s(1)%matrix)
            IF (output_unit > 0) &
               WRITE (output_unit, FMT="(/,(T3,A,T60,F20.10))") &
               "Max deviation from orthonormalization:", orthonormality
         END IF
      END IF

      IF (output_unit > 0) &
         CALL m_flush(output_unit)
      CALL cp_print_key_finished_output(output_unit, logger, rtp_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      IF (rtp%converged) THEN
         dft_section => section_vals_get_subs_vals(input, "DFT")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              dft_section, "REAL_TIME_PROPAGATION%PRINT%FIELD"), cp_p_file)) &
            CALL print_field_applied(qs_env, dft_section)
         CALL make_moment(qs_env)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              dft_section, "REAL_TIME_PROPAGATION%PRINT%E_CONSTITUENTS"), cp_p_file)) THEN
            CALL print_rtp_energy_components(qs_env, dft_section)
         END IF
         IF (.NOT. dft_control%qs_control%dftb) THEN
            CALL write_available_results(qs_env=qs_env, rtp=rtp)
         END IF

         IF (rtp%linear_scaling) THEN
            CALL get_rtp(rtp=rtp, rho_new=rho_new)

            ! Probably have to rebuild the moment matrix, since atoms can also move, in principle
            IF (dft_control%rtp_control%save_local_moments) THEN
               ! Save the field value
               IF (dft_control%apply_efield_field) THEN
                  CALL make_field(dft_control, field, qs_env%sim_step, qs_env%sim_time)
                  rtp%fields(:, rtp%istep + rtp%i_start + 1) = CMPLX(field(:), 0.0, kind=dp)
               END IF
               IF (.NOT. dft_control%rtp_control%fixed_ions) THEN
                  CALL build_local_moment_matrix(qs_env, rtp%local_moments, 1, reference_point)
               END IF
               ! TODO : Is symmetric rho possible?
               ! Spin + complex parts
               ! Extensions setup
               CALL calc_local_moment(rtp%local_moments, rho_new, &
                                      rtp%local_moments_work, rtp%moments(:, :, rtp%istep + rtp%i_start + 1))
               ! Time 1 is zero (start) time
               rtp%times(rtp%istep + rtp%i_start + 1) = qs_env%sim_time
               output_unit = cp_logger_get_default_io_unit(logger)
               CALL print_moments(section_vals_get_subs_vals(rtp_section, "PRINT%MOMENTS"), output_unit, &
                                  rtp%moments(:, :, rtp%istep + rtp%i_start + 1), qs_env%sim_time, rtp%track_imag_density)
            END IF

            IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                 dft_section, "REAL_TIME_PROPAGATION%PRINT%RESTART"), cp_p_file)) THEN
               CALL write_rt_p_to_restart(rho_new, .FALSE.)
            END IF
            IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                 dft_section, "REAL_TIME_PROPAGATION%PRINT%RESTART_HISTORY"), cp_p_file)) THEN
               CALL write_rt_p_to_restart(rho_new, .TRUE.)
            END IF
            IF (.NOT. dft_control%qs_control%dftb) THEN
               !Not sure if these things could also work with dftb or not
               IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                    dft_section, "REAL_TIME_PROPAGATION%PRINT%CURRENT"), cp_p_file)) THEN
                  DO spin = 1, SIZE(rho_new)/2
                     CALL rt_current(qs_env, rho_new(2*spin)%matrix, dft_section, spin, SIZE(rho_new)/2)
                  END DO
               END IF
            END IF
         ELSE
            CALL get_rtp(rtp=rtp, mos_new=mos_new)
            IF (.NOT. dft_control%qs_control%dftb .AND. .NOT. dft_control%qs_control%xtb) THEN
               IF (rtp%track_imag_density) THEN
                  NULLIFY (P_im, p_xyz)
                  CALL dbcsr_allocate_matrix_set(p_xyz, 3)

! Linear momentum operator
! prepare for allocation
                  natom = SIZE(particle_set, 1)
                  ALLOCATE (first_sgf(natom))
                  ALLOCATE (last_sgf(natom))
                  CALL get_particle_set(particle_set, qs_kind_set, &
                                        first_sgf=first_sgf, &
                                        last_sgf=last_sgf)
                  ALLOCATE (row_blk_sizes(natom))
                  CALL dbcsr_convert_offsets_to_sizes(first_sgf, row_blk_sizes, last_sgf)
                  DEALLOCATE (first_sgf)
                  DEALLOCATE (last_sgf)

                  ALLOCATE (p_xyz(1)%matrix)
                  CALL dbcsr_create(matrix=p_xyz(1)%matrix, &
                                    name="p_xyz", &
                                    dist=dbcsr_dist, matrix_type=dbcsr_type_antisymmetric, &
                                    row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                                    mutable_work=.TRUE.)
                  CALL cp_dbcsr_alloc_block_from_nbl(p_xyz(1)%matrix, sab_orb)
                  CALL dbcsr_set(p_xyz(1)%matrix, 0.0_dp)
                  DO i = 2, 3
                     ALLOCATE (p_xyz(i)%matrix)
                     CALL dbcsr_copy(p_xyz(i)%matrix, p_xyz(1)%matrix, "p_xyz-"//TRIM(ADJUSTL(cp_to_string(i))))
                     CALL dbcsr_set(p_xyz(i)%matrix, 0.0_dp)
                  END DO
                  CALL build_lin_mom_matrix(qs_env, p_xyz)
                  DEALLOCATE (row_blk_sizes)

                  nspin = SIZE(mos_new)/2
                  CALL qs_rho_get(rho, rho_ao_im=P_im)
                  ALLOCATE (j_int(nspin, 3))
                  j_int = 0.0_dp

                  NULLIFY (tmp_ao)
                  CALL dbcsr_init_p(tmp_ao)
                  CALL dbcsr_create(tmp_ao, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry, name="tmp")
                  CALL cp_dbcsr_alloc_block_from_nbl(tmp_ao, sab_all)
                  CALL dbcsr_set(tmp_ao, 0.0_dp)

                  DO i = 1, 3
                     strace = 0.0_dp
                     DO spin = 1, nspin
                        CALL dbcsr_set(tmp_ao, 0.0_dp)
                        CALL dbcsr_multiply("T", "N", 1.0_dp, P_im(spin)%matrix, p_xyz(i)%matrix, &
                                            0.0_dp, tmp_ao)
                        CALL dbcsr_trace(tmp_ao, trace)
                        strace = strace + trace
!     dft_control%rtp_control%vec_pot(1)
                        j_int(spin, i) = -trace + dft_control%rtp_control%vec_pot(i)*nelectron_spin(spin)
!!  j_int(spin, i) = strace
                     END DO
                  END DO
!  PP term missing

                  IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                       dft_section, "REAL_TIME_PROPAGATION%PRINT%CURRENT_INT"), cp_p_file)) THEN

                     output_unit = cp_logger_get_default_io_unit(logger)
                     unit_nr = cp_print_key_unit_nr(logger, dft_section, &
                                                  "REAL_TIME_PROPAGATION%PRINT%CURRENT_INT", extension=".dat", is_new_file=new_file)
                     IF (output_unit > 0) THEN
                        IF (nspin == 2) THEN
                           WRITE (UNIT=unit_nr, FMT="(I10,F16.6,6(F16.8,1X))") qs_env%sim_step, qs_env%sim_time*femtoseconds, &
                              j_int(1, 1:3), j_int(2, 1:3)
                        ELSE
                           WRITE (UNIT=unit_nr, FMT="(I10,F16.6,3(F16.8,1X))") qs_env%sim_step, qs_env%sim_time*femtoseconds, &
                              j_int(1, 1:3)
                        END IF
                     END IF
                     CALL cp_print_key_finished_output(unit_nr, logger, dft_section, &
                                                       "REAL_TIME_PROPAGATION%PRINT%CURRENT_INT")
                  END IF
                  DEALLOCATE (j_int)

                  IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                       dft_section, "REAL_TIME_PROPAGATION%PRINT%CURRENT"), cp_p_file)) THEN
                     DO spin = 1, nspin
                        CALL rt_current(qs_env, P_im(spin)%matrix, dft_section, spin, nspin)
                     END DO
                  END IF
                  CALL dbcsr_deallocate_matrix(tmp_ao)
                  CALL dbcsr_deallocate_matrix_set(p_xyz)
               END IF

! projection of molecular orbitals
               IF (dft_control%rtp_control%is_proj_mo) THEN
                  DO n_proj = 1, SIZE(dft_control%rtp_control%proj_mo_list)
                     CALL compute_and_write_proj_mo(qs_env, mos_new, &
                                                    dft_control%rtp_control%proj_mo_list(n_proj)%proj_mo, n_proj)
                  END DO
               END IF
            END IF
            CALL write_rt_mos_to_restart(qs_env%mos, mos_new, particle_set, &
                                         dft_section, qs_kind_set)
         END IF
      END IF

      rtp%energy_old = rtp%energy_new

      IF (.NOT. rtp%converged .AND. rtp%iter >= dft_control%rtp_control%max_iter) &
         CALL cp_abort(__LOCATION__, "EMD did not converge, either increase MAX_ITER "// &
                       "or use a smaller TIMESTEP")

   END SUBROUTINE rt_prop_output

! **************************************************************************************************
!> \brief computes the effective orthonormality of a set of mos given an s-matrix
!>        orthonormality is the max deviation from unity of the C^T S C
!> \param orthonormality ...
!> \param mos_new ...
!> \param matrix_s ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************
   SUBROUTINE rt_calculate_orthonormality(orthonormality, mos_new, matrix_s)
      REAL(KIND=dp), INTENT(out)                         :: orthonormality
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_s

      CHARACTER(len=*), PARAMETER :: routineN = 'rt_calculate_orthonormality'

      INTEGER                                            :: handle, i, im, ispin, j, k, n, &
                                                            ncol_local, nrow_local, nspin, re
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: alpha, max_alpha, max_beta
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type)                                   :: overlap_re, svec_im, svec_re

      NULLIFY (tmp_fm_struct)

      CALL timeset(routineN, handle)

      nspin = SIZE(mos_new)/2
      max_alpha = 0.0_dp
      max_beta = 0.0_dp
      DO ispin = 1, nspin
         re = ispin*2 - 1
         im = ispin*2
         ! get S*C
         CALL cp_fm_create(svec_re, mos_new(im)%matrix_struct)
         CALL cp_fm_create(svec_im, mos_new(im)%matrix_struct)
         CALL cp_fm_get_info(mos_new(im), &
                             nrow_global=n, ncol_global=k)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s, mos_new(re), &
                                      svec_re, k)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s, mos_new(im), &
                                      svec_im, k)

         ! get C^T (S*C)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=k, ncol_global=k, &
                                  para_env=mos_new(re)%matrix_struct%para_env, &
                                  context=mos_new(re)%matrix_struct%context)
         CALL cp_fm_create(overlap_re, tmp_fm_struct)

         CALL cp_fm_struct_release(tmp_fm_struct)

         CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, mos_new(re), &
                            svec_re, 0.0_dp, overlap_re)
         CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, mos_new(im), &
                            svec_im, 1.0_dp, overlap_re)

         CALL cp_fm_release(svec_re)
         CALL cp_fm_release(svec_im)

         CALL cp_fm_get_info(overlap_re, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices)
         DO i = 1, nrow_local
            DO j = 1, ncol_local
               alpha = overlap_re%local_data(i, j)
               IF (row_indices(i) == col_indices(j)) alpha = alpha - 1.0_dp
               max_alpha = MAX(max_alpha, ABS(alpha))
            END DO
         END DO
         CALL cp_fm_release(overlap_re)
      END DO
      CALL mos_new(1)%matrix_struct%para_env%max(max_alpha)
      CALL mos_new(1)%matrix_struct%para_env%max(max_beta)
      orthonormality = max_alpha

      CALL timestop(handle)

   END SUBROUTINE rt_calculate_orthonormality

! **************************************************************************************************
!> \brief computes the convergence criterion for RTP and EMD
!> \param rtp ...
!> \param matrix_s Overlap matrix without the derivatives
!> \param delta_mos ...
!> \param delta_eps ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE rt_convergence(rtp, matrix_s, delta_mos, delta_eps)
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: delta_mos
      REAL(dp), INTENT(out)                              :: delta_eps

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rt_convergence'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: handle, i, icol, im, ispin, j, lcol, &
                                                            lrow, nao, newdim, nmo, nspin, re
      LOGICAL                                            :: double_col, double_row
      REAL(KIND=dp)                                      :: alpha, max_alpha
      TYPE(cp_fm_struct_type), POINTER                   :: newstruct, newstruct1, tmp_fm_struct
      TYPE(cp_fm_type)                                   :: work, work1, work2
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new

      NULLIFY (tmp_fm_struct)

      CALL timeset(routineN, handle)

      CALL get_rtp(rtp=rtp, mos_new=mos_new)

      nspin = SIZE(delta_mos)/2
      max_alpha = 0.0_dp

      DO i = 1, SIZE(mos_new)
         CALL cp_fm_scale_and_add(-one, delta_mos(i), one, mos_new(i))
      END DO

      DO ispin = 1, nspin
         re = ispin*2 - 1
         im = ispin*2

         double_col = .TRUE.
         double_row = .FALSE.
         CALL cp_fm_struct_double(newstruct, &
                                  delta_mos(re)%matrix_struct, &
                                  delta_mos(re)%matrix_struct%context, &
                                  double_col, &
                                  double_row)

         CALL cp_fm_create(work, matrix_struct=newstruct)
         CALL cp_fm_create(work1, matrix_struct=newstruct)

         CALL cp_fm_get_info(delta_mos(re), ncol_local=lcol, ncol_global=nmo, &
                             nrow_global=nao)
         CALL cp_fm_get_info(work, ncol_global=newdim)

         CALL cp_fm_set_all(work, zero, zero)

         DO icol = 1, lcol
            work%local_data(:, icol) = delta_mos(re)%local_data(:, icol)
            work%local_data(:, icol + lcol) = delta_mos(im)%local_data(:, icol)
         END DO

         CALL cp_dbcsr_sm_fm_multiply(matrix_s, work, work1, ncol=newdim)

         CALL cp_fm_release(work)

         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmo, ncol_global=nmo, &
                                  para_env=delta_mos(re)%matrix_struct%para_env, &
                                  context=delta_mos(re)%matrix_struct%context)
         CALL cp_fm_struct_double(newstruct1, &
                                  tmp_fm_struct, &
                                  delta_mos(re)%matrix_struct%context, &
                                  double_col, &
                                  double_row)

         CALL cp_fm_create(work, matrix_struct=newstruct1)
         CALL cp_fm_create(work2, matrix_struct=newstruct1)

         CALL parallel_gemm("T", "N", nmo, newdim, nao, one, delta_mos(re), &
                            work1, zero, work)

         CALL parallel_gemm("T", "N", nmo, newdim, nao, one, delta_mos(im), &
                            work1, zero, work2)

         CALL cp_fm_get_info(work, nrow_local=lrow)
         DO i = 1, lrow
            DO j = 1, lcol
               alpha = SQRT((work%local_data(i, j) + work2%local_data(i, j + lcol))**2 + &
                            (work%local_data(i, j + lcol) - work2%local_data(i, j))**2)
               max_alpha = MAX(max_alpha, ABS(alpha))
            END DO
         END DO

         CALL cp_fm_release(work)
         CALL cp_fm_release(work1)
         CALL cp_fm_release(work2)
         CALL cp_fm_struct_release(tmp_fm_struct)
         CALL cp_fm_struct_release(newstruct)
         CALL cp_fm_struct_release(newstruct1)

      END DO

      CALL delta_mos(1)%matrix_struct%para_env%max(max_alpha)
      delta_eps = SQRT(max_alpha)

      CALL timestop(handle)

   END SUBROUTINE rt_convergence

! **************************************************************************************************
!> \brief computes the convergence criterion for RTP and EMD based on the density matrix
!> \param rtp ...
!> \param delta_P ...
!> \param delta_eps ...
!> \author Samuel Andermatt (02.14)
! **************************************************************************************************

   SUBROUTINE rt_convergence_density(rtp, delta_P, delta_eps)

      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: delta_P
      REAL(dp), INTENT(out)                              :: delta_eps

      CHARACTER(len=*), PARAMETER :: routineN = 'rt_convergence_density'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: col_atom, handle, i, ispin, row_atom
      REAL(dp)                                           :: alpha, max_alpha
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_values
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_new
      TYPE(dbcsr_type), POINTER                          :: tmp
      TYPE(mp_comm_type)                                 :: group

      CALL timeset(routineN, handle)

      CALL get_rtp(rtp=rtp, rho_new=rho_new)

      DO i = 1, SIZE(rho_new)
         CALL dbcsr_add(delta_P(i)%matrix, rho_new(i)%matrix, one, -one)
      END DO
      !get the maximum value of delta_P
      DO i = 1, SIZE(delta_P)
         !square all entries of both matrices
         CALL dbcsr_iterator_start(iter, delta_P(i)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
            block_values = block_values*block_values
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO
      NULLIFY (tmp)
      ALLOCATE (tmp)
      CALL dbcsr_create(tmp, template=delta_P(1)%matrix, matrix_type="N")
      DO ispin = 1, SIZE(delta_P)/2
         CALL dbcsr_desymmetrize(delta_P(2*ispin - 1)%matrix, tmp)
         CALL dbcsr_add(delta_P(2*ispin)%matrix, tmp, one, one)
      END DO
      !the absolute values are now in the even entries of delta_P
      max_alpha = zero
      DO ispin = 1, SIZE(delta_P)/2
         CALL dbcsr_iterator_start(iter, delta_P(2*ispin)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
            alpha = MAXVAL(block_values)
            IF (alpha > max_alpha) max_alpha = alpha
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO
      CALL dbcsr_get_info(delta_P(1)%matrix, group=group)
      CALL group%max(max_alpha)
      delta_eps = SQRT(max_alpha)
      CALL dbcsr_deallocate_matrix(tmp)
      CALL timestop(handle)

   END SUBROUTINE rt_convergence_density

! **************************************************************************************************
!> \brief interface to qs_moments. Does only work for nonperiodic dipole
!> \param qs_env ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE make_moment(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'make_moment'

      INTEGER                                            :: handle, output_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%qs_control%dftb) THEN
         CALL scf_post_calculation_tb(qs_env, "DFTB", .FALSE.)
      ELSE IF (dft_control%qs_control%xtb) THEN
         CALL scf_post_calculation_tb(qs_env, "xTB", .FALSE.)
      ELSE
         CALL qs_scf_post_moments(qs_env%input, logger, qs_env, output_unit)
      END IF
      CALL timestop(handle)

   END SUBROUTINE make_moment

! **************************************************************************************************
!> \brief Reports the sparsity pattern of the complex density matrix
!> \param filter_eps ...
!> \param rho ...
!> \author Samuel Andermatt (09.14)
! **************************************************************************************************

   SUBROUTINE report_density_occupation(filter_eps, rho)

      REAL(KIND=dp)                                      :: filter_eps
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho

      CHARACTER(len=*), PARAMETER :: routineN = 'report_density_occupation'

      INTEGER                                            :: handle, i, im, ispin, re, unit_nr
      REAL(KIND=dp)                                      :: eps, occ
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: tmp

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_io_unit(logger)
      NULLIFY (tmp)
      CALL dbcsr_allocate_matrix_set(tmp, SIZE(rho))
      DO i = 1, SIZE(rho)
         CALL dbcsr_init_p(tmp(i)%matrix)
         CALL dbcsr_create(tmp(i)%matrix, template=rho(i)%matrix)
         CALL dbcsr_copy(tmp(i)%matrix, rho(i)%matrix)
      END DO
      DO ispin = 1, SIZE(rho)/2
         re = 2*ispin - 1
         im = 2*ispin
         eps = MAX(filter_eps, 1.0E-11_dp)
         DO WHILE (eps < 1.1_dp)
            CALL dbcsr_filter(tmp(re)%matrix, eps)
            occ = dbcsr_get_occupation(tmp(re)%matrix)
            IF (unit_nr > 0) WRITE (unit_nr, FMT="((T3,A,I1,A,F15.12,A,T61,F20.10))") "Occupation of rho spin ", &
               ispin, " eps ", eps, " real: ", occ
            eps = eps*10
         END DO
         eps = MAX(filter_eps, 1.0E-11_dp)
         DO WHILE (eps < 1.1_dp)
            CALL dbcsr_filter(tmp(im)%matrix, eps)
            occ = dbcsr_get_occupation(tmp(im)%matrix)
            IF (unit_nr > 0) WRITE (unit_nr, FMT="((T3,A,I1,A,F15.12,A,T61,F20.10))") "Occupation of rho spin ", &
               ispin, " eps ", eps, " imag: ", occ
            eps = eps*10.0_dp
         END DO
      END DO
      CALL dbcsr_deallocate_matrix_set(tmp)
      CALL timestop(handle)

   END SUBROUTINE report_density_occupation

! **************************************************************************************************
!> \brief Writes the density matrix and the atomic positions to a restart file
!> \param rho_new ...
!> \param history ...
!> \author Samuel Andermatt (09.14)
! **************************************************************************************************

   SUBROUTINE write_rt_p_to_restart(rho_new, history)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_new
      LOGICAL                                            :: history

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_rt_p_to_restart'

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: handle, im, ispin, re, unit_nr
      REAL(KIND=dp)                                      :: cs_pos
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      project_name = logger%iter_info%project_name
      DO ispin = 1, SIZE(rho_new)/2
         re = 2*ispin - 1
         im = 2*ispin
         IF (history) THEN
            WRITE (file_name, '(A,I0,A)') &
               TRIM(project_name)//"_LS_DM_SPIN_RE", ispin, "_"//TRIM(cp_iter_string(logger%iter_info))//"_RESTART.dm"
         ELSE
            WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_RE", ispin, "_RESTART.dm"
         END IF
         cs_pos = dbcsr_checksum(rho_new(re)%matrix, pos=.TRUE.)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,E20.8)') "Writing restart DM "//TRIM(file_name)//" with checksum: ", cs_pos
         END IF
         CALL dbcsr_binary_write(rho_new(re)%matrix, file_name)
         IF (history) THEN
            WRITE (file_name, '(A,I0,A)') &
               TRIM(project_name)//"_LS_DM_SPIN_IM", ispin, "_"//TRIM(cp_iter_string(logger%iter_info))//"_RESTART.dm"
         ELSE
            WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_IM", ispin, "_RESTART.dm"
         END IF
         cs_pos = dbcsr_checksum(rho_new(im)%matrix, pos=.TRUE.)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,E20.8)') "Writing restart DM "//TRIM(file_name)//" with checksum: ", cs_pos
         END IF
         CALL dbcsr_binary_write(rho_new(im)%matrix, file_name)
      END DO

      CALL timestop(handle)

   END SUBROUTINE write_rt_p_to_restart

! **************************************************************************************************
!> \brief Collocation of the current and printing of it in a cube file
!> \param qs_env ...
!> \param P_im ...
!> \param dft_section ...
!> \param spin ...
!> \param nspin ...
!> \author Samuel Andermatt (06.15)
! **************************************************************************************************
   SUBROUTINE rt_current(qs_env, P_im, dft_section, spin, nspin)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: P_im
      TYPE(section_vals_type), POINTER                   :: dft_section
      INTEGER                                            :: spin, nspin

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rt_current'

      CHARACTER(len=1)                                   :: char_spin
      CHARACTER(len=14)                                  :: ext
      CHARACTER(len=2)                                   :: sdir
      INTEGER                                            :: dir, handle, print_unit
      INTEGER, DIMENSION(:), POINTER                     :: stride(:)
      LOGICAL                                            :: mpi_io
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(current_env_type)                             :: current_env
      TYPE(dbcsr_type), POINTER                          :: tmp, zero
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_c1d_gs_type)                               :: gs
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: rs
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      CALL get_qs_env(qs_env=qs_env, subsys=subsys, pw_env=pw_env)
      CALL qs_subsys_get(subsys, particles=particles)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      NULLIFY (zero, tmp)
      ALLOCATE (zero, tmp)
      CALL dbcsr_create(zero, template=P_im)
      CALL dbcsr_copy(zero, P_im)
      CALL dbcsr_set(zero, 0.0_dp)
      CALL dbcsr_create(tmp, template=P_im)
      CALL dbcsr_copy(tmp, P_im)
      IF (nspin == 1) THEN
         CALL dbcsr_scale(tmp, 0.5_dp)
      END IF
      current_env%gauge = -1
      current_env%gauge_init = .FALSE.
      CALL auxbas_pw_pool%create_pw(rs)
      CALL auxbas_pw_pool%create_pw(gs)

      NULLIFY (stride)
      ALLOCATE (stride(3))

      DO dir = 1, 3

         CALL pw_zero(rs)
         CALL pw_zero(gs)

         CALL calculate_jrho_resp(zero, tmp, zero, zero, dir, dir, rs, gs, qs_env, current_env, retain_rsgrid=.TRUE.)

         stride = section_get_ivals(dft_section, "REAL_TIME_PROPAGATION%PRINT%CURRENT%STRIDE")

         IF (dir == 1) THEN
            sdir = "-x"
         ELSEIF (dir == 2) THEN
            sdir = "-y"
         ELSE
            sdir = "-z"
         END IF
         WRITE (char_spin, "(I1)") spin

         ext = "-SPIN-"//char_spin//sdir//".cube"
         mpi_io = .TRUE.
         print_unit = cp_print_key_unit_nr(logger, dft_section, "REAL_TIME_PROPAGATION%PRINT%CURRENT", &
                                           extension=ext, file_status="REPLACE", file_action="WRITE", &
                                           log_filename=.FALSE., mpi_io=mpi_io)

         CALL cp_pw_to_cube(rs, print_unit, "EMD current", particles=particles, stride=stride, &
                            mpi_io=mpi_io)

         CALL cp_print_key_finished_output(print_unit, logger, dft_section, "REAL_TIME_PROPAGATION%PRINT%CURRENT", &
                                           mpi_io=mpi_io)

      END DO

      CALL auxbas_pw_pool%give_back_pw(rs)
      CALL auxbas_pw_pool%give_back_pw(gs)

      CALL dbcsr_deallocate_matrix(zero)
      CALL dbcsr_deallocate_matrix(tmp)

      DEALLOCATE (stride)

      CALL timestop(handle)

   END SUBROUTINE rt_current

! **************************************************************************************************
!> \brief Interface routine to trigger writing of results available from normal
!>        SCF. Can write MO-dependent and MO free results (needed for call from
!>        the linear scaling code)
!>        Update: trigger also some of prints for time-dependent runs
!> \param qs_env ...
!> \param rtp ...
!> \par History
!>      2022-11 Update [Guillaume Le Breton]
! **************************************************************************************************
   SUBROUTINE write_available_results(qs_env, rtp)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rt_prop_type), POINTER                        :: rtp

      CHARACTER(len=*), PARAMETER :: routineN = 'write_available_results'

      INTEGER                                            :: handle
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, scf_env=scf_env)
      IF (rtp%linear_scaling) THEN
         CALL write_mo_free_results(qs_env)
      ELSE
         CALL write_mo_free_results(qs_env)
         CALL write_mo_dependent_results(qs_env, scf_env)
         ! Time-dependent MO print
         CALL write_rtp_mos_to_output_unit(qs_env, rtp)
         CALL write_rtp_mo_cubes(qs_env, rtp)
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_available_results

! **************************************************************************************************
!> \brief Print the field applied to the system. Either the electric
!>        field or the vector potential depending on the gauge used
!> \param qs_env ...
!> \param dft_section ...
!> \par History
!>      2023-01  Created [Guillaume Le Breton]
! **************************************************************************************************
   SUBROUTINE print_field_applied(qs_env, dft_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(LEN=3), DIMENSION(3)                     :: rlab
      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: i, i_step, output_unit, unit_nr
      LOGICAL                                            :: new_file
      REAL(kind=dp)                                      :: field(3)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(rt_prop_type), POINTER                        :: rtp

      NULLIFY (dft_control)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, dft_control=dft_control, rtp=rtp)

      i_step = rtp%istep

      unit_nr = cp_print_key_unit_nr(logger, dft_section, &
                                     "REAL_TIME_PROPAGATION%PRINT%FIELD", extension=".dat", is_new_file=new_file)

      IF (output_unit > 0) THEN
         rlab = [CHARACTER(LEN=3) :: "X", "Y", "Z"]
         IF (unit_nr /= output_unit) THEN
            INQUIRE (UNIT=unit_nr, NAME=filename)
            WRITE (UNIT=output_unit, FMT="(/,T2,A,2(/,T3,A),/)") &
               "FIELD", "The field applied is written to the file:", &
               TRIM(filename)
         ELSE
            WRITE (UNIT=output_unit, FMT="(/,T2,A)") "FIELD APPLIED [a.u.]"
            WRITE (UNIT=output_unit, FMT="(T5,3(A,A,E16.8,1X))") &
               (TRIM(rlab(i)), "=", dft_control%rtp_control%field(i), i=1, 3)
         END IF

         IF (new_file) THEN
            IF (dft_control%apply_efield_field) THEN
               WRITE (UNIT=unit_nr, FMT='("#",5X,A,8X,A,3(6X,A))') "Step Nr.", "Time[fs]", " Field X", "   Field Y", "   Field Z"
            ELSE IF (dft_control%apply_vector_potential) THEN
               WRITE (UNIT=unit_nr, FMT='("#",5X,A,8X,A,6(6X,A))') "Step Nr.", "Time[fs]", " Field X", "   Field Y", "   Field Z", &
                  "  Vec. Pot. X", "  Vec. Pot. Y", "    Vec. Pot. Z"
            END IF
         END IF

         field = 0.0_dp
         IF (dft_control%apply_efield_field) THEN
            CALL make_field(dft_control, field, qs_env%sim_step, qs_env%sim_time)
            WRITE (UNIT=unit_nr, FMT="(I10,F16.6,3(F16.8,1X))") qs_env%sim_step, qs_env%sim_time*femtoseconds, &
               field(1), field(2), field(3)
!            DO i=1,3
!               IF (ABS(field(i))< 10E-10) field(i) = 0.0_dp
!            END IF
         ELSE IF (dft_control%apply_vector_potential) THEN
            WRITE (UNIT=unit_nr, FMT="(I10,F16.6,6(F16.8,1X))") qs_env%sim_step, qs_env%sim_time*femtoseconds, &
               dft_control%rtp_control%field(1), dft_control%rtp_control%field(2), dft_control%rtp_control%field(3), &
               dft_control%rtp_control%vec_pot(1), dft_control%rtp_control%vec_pot(2), dft_control%rtp_control%vec_pot(3)
         END IF

      END IF

      CALL cp_print_key_finished_output(unit_nr, logger, dft_section, &
                                        "REAL_TIME_PROPAGATION%PRINT%FIELD")

   END SUBROUTINE print_field_applied

! **************************************************************************************************
!> \brief Print the components of the total energy used in an RTP calculation
!> \param qs_env ...
!> \param dft_section ...
!> \par History
!>      2024-02  Created [ANB]
! **************************************************************************************************
   SUBROUTINE print_rtp_energy_components(qs_env, dft_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: i_step, output_unit, unit_nr
      LOGICAL                                            :: new_file
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(rt_prop_type), POINTER                        :: rtp

      NULLIFY (dft_control, energy, rtp)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, dft_control=dft_control, rtp=rtp, energy=energy)
      i_step = rtp%istep

      unit_nr = cp_print_key_unit_nr(logger, dft_section, &
                                     "REAL_TIME_PROPAGATION%PRINT%E_CONSTITUENTS", extension=".ener", &
                                     file_action="WRITE", is_new_file=new_file)

      IF (output_unit > 0) THEN
         IF (unit_nr /= output_unit) THEN
            INQUIRE (UNIT=unit_nr, NAME=filename)
            WRITE (UNIT=output_unit, FMT="(/,T2,A,2(/,T3,A),/)") &
               "ENERGY_CONSTITUENTS", "Total Energy constituents written to file:", &
               TRIM(filename)
         ELSE
            WRITE (UNIT=output_unit, FMT="(/,T2,A)") "ENERGY_CONSTITUENTS"
         END IF

         IF (new_file) THEN
            ! NOTE that these are not all terms contributing to the total energy for RTP, only a selection of those
            ! most significant / impactful. Therefore the printed components likely will not add up to the total energy.
            WRITE (UNIT=unit_nr, FMT='("#",5X,A,8X,A,10(6X,A))') "Step Nr.", "Time[fs]", &
               "Total ener.[a.u.]", "core[a.u.]   ", " overlap [a.u.]", "hartree[a.u.]", " exc. [a.u.] ", &
               " hartree 1c[a.u.]", "exc 1c[a.u.] ", "exc admm[a.u.]", "exc 1c admm[a.u.]", "efield LG"

         END IF
         WRITE (UNIT=unit_nr, FMT="(I10,F20.6,10(F20.9))") &
            qs_env%sim_step, qs_env%sim_time*femtoseconds, &
            energy%total, energy%core, energy%core_overlap, energy%hartree, energy%exc, &
            energy%hartree_1c, energy%exc1, energy%exc_aux_fit, energy%exc1_aux_fit, energy%efield_core

      END IF

      CALL cp_print_key_finished_output(unit_nr, logger, dft_section, &
                                        "REAL_TIME_PROPAGATION%PRINT%E_CONSTITUENTS")

   END SUBROUTINE print_rtp_energy_components

! **************************************************************************************************
!> \brief Print the dipole moments into a file
!> \param moments_section Section of the input defining the file/stream to print the moments to
!> \param info_unit Unit where standard output from the program is written - for add. identifiers
!> \param moments Actual moment values (for specific time step)
!> \param time Current simulation time
!> \param imag_opt Whether to calculate the imaginary part
!> \par History
!>      10.2025  Created [Marek]
! **************************************************************************************************
   SUBROUTINE print_moments(moments_section, info_unit, moments, time, imag_opt)
      TYPE(section_vals_type), POINTER                   :: moments_section
      INTEGER                                            :: info_unit
      COMPLEX(kind=dp), DIMENSION(:, :)                  :: moments
      REAL(kind=dp), OPTIONAL                            :: time
      LOGICAL, OPTIONAL                                  :: imag_opt

      CHARACTER(len=14), DIMENSION(4)                    :: file_extensions
      INTEGER                                            :: i, j, ndir, nspin, print_unit
      LOGICAL                                            :: imaginary
      TYPE(cp_logger_type), POINTER                      :: logger

! Index 1 : spin, Index 2 : direction

      nspin = SIZE(moments, 1)
      ndir = SIZE(moments, 2)

      IF (nspin < 1) CPABORT("Zero spin index size in print moments!")
      IF (ndir < 1) CPABORT("Zero direction index size in print moments!")

      imaginary = .TRUE.
      IF (PRESENT(imag_opt)) imaginary = imag_opt

      ! Get the program run info unit and target unit
      ! If these are the same (most likely the case of __STD_OUT__), add
      ! extra identifier to the printed output
      file_extensions(1) = "_SPIN_A_RE.dat"
      file_extensions(2) = "_SPIN_A_IM.dat"
      file_extensions(3) = "_SPIN_B_RE.dat"
      file_extensions(4) = "_SPIN_B_IM.dat"
      logger => cp_get_default_logger()
      DO i = 1, nspin
         ! Real part
         print_unit = cp_print_key_unit_nr(logger, moments_section, extension=file_extensions(2*i - 1))
         IF (print_unit == info_unit .AND. print_unit > 0) THEN
            ! Do the output with additional suffix
            WRITE (print_unit, "(A18)", advance="no") " MOMENTS_TRACE_RE|"
         END IF
         IF (print_unit > 0) THEN
            IF (PRESENT(time)) WRITE (print_unit, "(E20.8E3)", advance="no") time*femtoseconds
            DO j = 1, ndir - 1
               WRITE (print_unit, "(E20.8E3)", advance="no") REAL(moments(i, j))
            END DO
            ! Write the last direction
            WRITE (print_unit, "(E20.8E3)") REAL(moments(i, ndir))
         END IF
         CALL cp_print_key_finished_output(print_unit, logger, moments_section)
         ! Same for imaginary part
         IF (imaginary) THEN
            print_unit = cp_print_key_unit_nr(logger, moments_section, extension=file_extensions(2*i))
            IF (print_unit == info_unit .AND. print_unit > 0) THEN
               ! Do the output with additional suffix
               WRITE (print_unit, "(A18)", advance="no") " MOMENTS_TRACE_IM|"
            END IF
            IF (print_unit > 0) THEN
               IF (PRESENT(time)) WRITE (print_unit, "(E20.8E3)", advance="no") time*femtoseconds
               DO j = 1, ndir - 1
                  WRITE (print_unit, "(E20.8E3)", advance="no") AIMAG(moments(i, j))
               END DO
               ! Write the last direction
               WRITE (print_unit, "(E20.8E3)") AIMAG(moments(i, ndir))
            END IF
            CALL cp_print_key_finished_output(print_unit, logger, moments_section)
         END IF
      END DO

   END SUBROUTINE print_moments

! **************************************************************************************************
!> \brief Calculate the values of real/imaginary parts of moments in all directions
!> \param moment_matrices Local matrix representations of dipole (position) operator
!> \param density_matrices Density matrices (spin and real+complex parts)
!> \param work Extra dbcsr matrix for work
!> \param moment Resulting moments (spin and direction)
!> \param imag_opt Whether to calculate the imaginary part of the moment
!> \par History
!>      10.2025  Created [Marek]
! **************************************************************************************************
   SUBROUTINE calc_local_moment(moment_matrices, density_matrices, work, moment, imag_opt)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moment_matrices, density_matrices
      TYPE(dbcsr_type)                                   :: work
      COMPLEX(kind=dp), DIMENSION(:, :)                  :: moment
      LOGICAL, OPTIONAL                                  :: imag_opt

      INTEGER                                            :: i, k, nspin
      LOGICAL                                            :: imag
      REAL(kind=dp)                                      :: real_moment

      imag = .FALSE.
      IF (PRESENT(imag_opt)) imag = imag_opt
      nspin = SIZE(density_matrices)/2

      DO i = 1, nspin
         DO k = 1, 3
            CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                density_matrices(2*i - 1)%matrix, moment_matrices(k)%matrix, &
                                0.0_dp, work)
            CALL dbcsr_trace(work, real_moment)
            moment(i, k) = CMPLX(real_moment, 0.0, kind=dp)
            IF (imag) THEN
               CALL dbcsr_multiply("N", "N", -1.0_dp, &
                                   density_matrices(2*i)%matrix, moment_matrices(k)%matrix, &
                                   0.0_dp, work)
               CALL dbcsr_trace(work, real_moment)
               moment(i, k) = moment(i, k) + CMPLX(0.0, real_moment, kind=dp)
            END IF
         END DO
      END DO

   END SUBROUTINE calc_local_moment

! **************************************************************************************************
!> \brief Calculate and print the Fourier transforms + polarizabilites from moment trace
!> \param rtp_section The RTP input section (needed to access PRINT configurations)
!> \param moments Moment trace
!> \param times Corresponding times
!> \param fields Corresponding fields
!> \param rtc rt_control_type that includes metadata
!> \param info_opt ...
!> \param cell If present, used to change the delta peak representation to be in units of reciprocal lattice
!> \par History
!>      10.2025  Created [Marek]
! **************************************************************************************************
   SUBROUTINE print_ft(rtp_section, moments, times, fields, rtc, info_opt, cell)
      TYPE(section_vals_type), POINTER                   :: rtp_section
      COMPLEX(kind=dp), DIMENSION(:, :, :), POINTER      :: moments
      REAL(kind=dp), DIMENSION(:), POINTER               :: times
      COMPLEX(kind=dp), DIMENSION(:, :), POINTER         :: fields
      TYPE(rtp_control_type), POINTER                    :: rtc
      INTEGER, OPTIONAL                                  :: info_opt
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell

      CHARACTER(len=11), DIMENSION(2)                    :: file_extensions
      CHARACTER(len=20), ALLOCATABLE, DIMENSION(:)       :: headers
      CHARACTER(len=21)                                  :: prefix
      CHARACTER(len=5)                                   :: prefix_format
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: omegas_complex, omegas_pade
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :)     :: field_results, field_results_pade, &
                                                            pol_results, pol_results_pade, &
                                                            results, results_pade, value_series
      INTEGER                                            :: ft_unit, i, info_unit, k, n, n_elems, &
                                                            n_pade, nspin
      LOGICAL                                            :: do_moments_ft, do_polarizability
      REAL(kind=dp)                                      :: damping, t0
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: omegas, omegas_pade_real
      REAL(kind=dp), DIMENSION(3)                        :: delta_vec
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: moment_ft_section, pol_section

! For results, using spin * direction for first index, e.g. for nspin = 2
! results(1,:) = (spin=1 and direction=1,:),
! results(5,:) = (spin=2 and direction=2,:)

      logger => cp_get_default_logger()

      moment_ft_section => section_vals_get_subs_vals(rtp_section, "PRINT%MOMENTS_FT")
      pol_section => section_vals_get_subs_vals(rtp_section, "PRINT%POLARIZABILITY")

      nspin = SIZE(moments, 1)
      n = SIZE(times)
      n_elems = SIZE(rtc%print_pol_elements, 1)

      info_unit = -1
      IF (PRESENT(info_opt)) info_unit = info_opt

      ! NOTE : Allows for at most 2 spin species
      file_extensions(1) = "_SPIN_A.dat"
      file_extensions(2) = "_SPIN_B.dat"

      ! Determine whether MOMENTS_FT and/or polarizability needs to be calculated
      do_moments_ft = cp_printkey_is_on(logger%iter_info, moment_ft_section)
      do_polarizability = cp_printkey_is_on(logger%iter_info, pol_section)
      do_polarizability = do_polarizability .AND. (n_elems > 0)

      damping = rtc%ft_damping
      t0 = rtc%ft_t0

      ! Determine field ft if polarizability required
      IF (do_polarizability) THEN
         ALLOCATE (field_results(3, n))
         IF (rtc%apply_delta_pulse) THEN
            ! Constant real FT
            IF (PRESENT(cell)) THEN
               delta_vec(:) = (REAL(rtc%delta_pulse_direction(1), kind=dp)*cell%h_inv(1, :) + &
                               REAL(rtc%delta_pulse_direction(2), kind=dp)*cell%h_inv(2, :) + &
                               REAL(rtc%delta_pulse_direction(3), kind=dp)*cell%h_inv(3, :)) &
                              *twopi*rtc%delta_pulse_scale
            ELSE
               delta_vec(:) = REAL(rtc%delta_pulse_direction(:), kind=dp)*rtc%delta_pulse_scale
            END IF
            DO k = 1, 3
               field_results(k, :) = CMPLX(delta_vec(k), 0.0, kind=dp)
            END DO
         ELSE
            ! Do explicit FT
            CALL multi_fft(times, fields, field_results, &
                           damping_opt=damping, t0_opt=t0, subtract_initial_opt=.TRUE.)
         END IF
      END IF

      IF (do_moments_ft .OR. do_polarizability) THEN
         ! We need to transform at least the moments
         ! NOTE : Might be able to save some memory by only doing FT of actually
         ! required moments, but for now, doing FT of all moment directions
         ALLOCATE (results(3*nspin, n))
         ALLOCATE (omegas(n))
         ALLOCATE (value_series(3*nspin, n))
         DO i = 1, nspin
            DO k = 1, 3
               value_series(3*(i - 1) + k, :) = moments(i, k, :)
            END DO
         END DO
         ! TODO : Choose whether the initial subtraction is applied in &FT section?
         CALL multi_fft(times, value_series, results, omegas, &
                        damping_opt=damping, t0_opt=t0, subtract_initial_opt=.TRUE.)
         DEALLOCATE (value_series)
         DO i = 1, nspin
            ! Output to FT file, if needed
            ft_unit = cp_print_key_unit_nr(logger, moment_ft_section, extension=file_extensions(i), &
                                           file_form="FORMATTED", file_position="REWIND")
            ! Print header
            IF (ft_unit > 0) THEN
               ALLOCATE (headers(7))
               headers(2) = "      x,real [at.u.]"
               headers(3) = "      x,imag [at.u.]"
               headers(4) = "      y,real [at.u.]"
               headers(5) = "      y,imag [at.u.]"
               headers(6) = "      z,real [at.u.]"
               headers(7) = "      z,imag [at.u.]"
               IF (info_unit == ft_unit) THEN
                  headers(1) = "#        Energy [eV]"
                  prefix = " MOMENTS_FT|"
                  prefix_format = "(A12)"
                  CALL print_ft_file(ft_unit, headers, omegas, results(3*(i - 1) + 1:3*(i - 1) + 3, :), &
                                     prefix, prefix_format, evolt)
               ELSE
                  headers(1) = "#      omega [at.u.]"
                  CALL print_ft_file(ft_unit, headers, omegas, results(3*(i - 1) + 1:3*(i - 1) + 3, :))
               END IF
               DEALLOCATE (headers)
            END IF
            CALL cp_print_key_finished_output(ft_unit, logger, moment_ft_section)
         END DO
      END IF

      IF (rtc%pade_requested .AND. (do_moments_ft .OR. do_polarizability)) THEN
         ALLOCATE (omegas_complex(SIZE(omegas)))
         omegas_complex(:) = CMPLX(omegas(:), 0.0, kind=dp)
         n_pade = INT((rtc%pade_e_max - rtc%pade_e_min)/rtc%pade_e_step)
         ALLOCATE (omegas_pade(n_pade))
         ALLOCATE (omegas_pade_real(n_pade))
         ! Construct omegas_pade and omegas_complex
         DO i = 1, n_pade
            omegas_pade_real(i) = (i - 1)*rtc%pade_e_step + rtc%pade_e_min
            omegas_pade(i) = CMPLX(omegas_pade_real(i), 0.0, kind=dp)
         END DO
         ALLOCATE (results_pade(nspin*3, n_pade), source=CMPLX(0.0, 0.0, kind=dp))
         DO i = 1, nspin
            DO k = 1, 3
               CALL greenx_refine_ft(rtc%pade_fit_e_min, rtc%pade_fit_e_max, omegas_complex, results(3*(i - 1) + k, :), &
                                     omegas_pade, results_pade(3*(i - 1) + k, :))
            END DO
            ! Print to a file
            ft_unit = cp_print_key_unit_nr(logger, moment_ft_section, extension="_PADE"//file_extensions(i), &
                                           file_form="FORMATTED", file_position="REWIND")
            IF (ft_unit > 0) THEN
               ALLOCATE (headers(7))
               headers(2) = " x,real,pade [at.u.]"
               headers(3) = " x,imag,pade [at.u.]"
               headers(4) = " y,real,pade [at.u.]"
               headers(5) = " y,imag,pade [at.u.]"
               headers(6) = " z,real,pade [at.u.]"
               headers(7) = " z,imag,pade [at.u.]"
               IF (info_unit == ft_unit) THEN
                  headers(1) = "#        Energy [eV]"
                  prefix = " MOMENTS_FT_PADE|"
                  prefix_format = "(A17)"
                  CALL print_ft_file(ft_unit, headers, omegas_pade_real, results_pade(3*(i - 1) + 1:3*(i - 1) + 3, :), &
                                     prefix, prefix_format, evolt)
               ELSE
                  headers(1) = "#      omega [at.u.]"
                  CALL print_ft_file(ft_unit, headers, omegas_pade_real, results_pade(3*(i - 1) + 1:3*(i - 1) + 3, :))
               END IF
               DEALLOCATE (headers)
            END IF
         END DO
      END IF

      IF (do_polarizability) THEN
         ! get the polarizability elements, as required
         ALLOCATE (pol_results(n_elems, n))
         DO i = 1, nspin
            DO k = 1, n_elems
               ! NOTE - field is regularized to small value
               pol_results(k, :) = results(3*(i - 1) + &
                                           rtc%print_pol_elements(k, 1), :)/ &
                                   (field_results(rtc%print_pol_elements(k, 2), :) + &
                                    1.0e-10*field_results(rtc%print_pol_elements(k, 2), 2))
            END DO
            ! Print to the file
            ft_unit = cp_print_key_unit_nr(logger, pol_section, extension=file_extensions(i), &
                                           file_form="FORMATTED", file_position="REWIND")
            IF (ft_unit > 0) THEN
               ALLOCATE (headers(2*n_elems + 1))
               DO k = 1, n_elems
                  WRITE (headers(2*k), "(A16,I2,I2)") "real pol. elem.", &
                     rtc%print_pol_elements(k, 1), &
                     rtc%print_pol_elements(k, 2)
                  WRITE (headers(2*k + 1), "(A16,I2,I2)") "imag pol. elem.", &
                     rtc%print_pol_elements(k, 1), &
                     rtc%print_pol_elements(k, 2)
               END DO
               ! Write header
               IF (info_unit == ft_unit) THEN
                  headers(1) = "#        Energy [eV]"
                  prefix = " POLARIZABILITY|"
                  prefix_format = "(A16)"
                  CALL print_ft_file(ft_unit, headers, omegas, pol_results, &
                                     prefix, prefix_format, evolt)
               ELSE
                  headers(1) = "#      omega [at.u.]"
                  CALL print_ft_file(ft_unit, headers, omegas, pol_results)
               END IF
               DEALLOCATE (headers)
            END IF
            CALL cp_print_key_finished_output(ft_unit, logger, pol_section)
         END DO
      END IF

      ! Padé polarizability
      IF (rtc%pade_requested .AND. do_polarizability) THEN
         ! Start with the field pade
         ALLOCATE (field_results_pade(3, n_pade))
         IF (rtc%apply_delta_pulse) THEN
            DO k = 1, 3
               field_results_pade(k, :) = CMPLX(delta_vec(k), 0.0, kind=dp)
            END DO
         ELSE
            DO k = 1, 3
               CALL greenx_refine_ft(rtc%pade_fit_e_min, rtc%pade_fit_e_max, &
                                     omegas_complex, field_results(k, :), &
                                     omegas_pade, field_results_pade(k, :))
            END DO
         END IF
         ! Allocate polarisation pade
         ALLOCATE (pol_results_pade(n_elems, n_pade))
         ! Refine
         DO i = 1, nspin
            DO k = 1, n_elems
               ! NOTE : Regularization to small value
               pol_results_pade(k, :) = results_pade(3*(i - 1) + rtc%print_pol_elements(k, 1), :)/( &
                                        field_results_pade(rtc%print_pol_elements(k, 2), :) + &
                                        field_results_pade(rtc%print_pol_elements(k, 2), 2)*1.0e-10_dp)
            END DO
            ! Print to the file
            ft_unit = cp_print_key_unit_nr(logger, pol_section, extension="_PADE"//file_extensions(i), &
                                           file_form="FORMATTED", file_position="REWIND")
            IF (ft_unit > 0) THEN
               ALLOCATE (headers(2*n_elems + 1))
               DO k = 1, n_elems
                  WRITE (headers(2*k), "(A16,I2,I2)") "re,pade,pol.", &
                     rtc%print_pol_elements(k, 1), &
                     rtc%print_pol_elements(k, 2)
                  WRITE (headers(2*k + 1), "(A16,I2,I2)") "im,pade,pol.", &
                     rtc%print_pol_elements(k, 1), &
                     rtc%print_pol_elements(k, 2)
               END DO
               ! Write header
               IF (info_unit == ft_unit) THEN
                  headers(1) = "#        Energy [eV]"
                  prefix = " POLARIZABILITY_PADE|"
                  prefix_format = "(A21)"
                  CALL print_ft_file(ft_unit, headers, omegas_pade_real, pol_results_pade, &
                                     prefix, prefix_format, evolt)
               ELSE
                  headers(1) = "#      omega [at.u.]"
                  CALL print_ft_file(ft_unit, headers, omegas_pade_real, pol_results_pade)
               END IF
               DEALLOCATE (headers)
            END IF
            CALL cp_print_key_finished_output(ft_unit, logger, pol_section)
         END DO
         DEALLOCATE (field_results_pade)
         DEALLOCATE (pol_results_pade)
      END IF

      IF (rtc%pade_requested .AND. (do_moments_ft .OR. do_polarizability)) THEN
         DEALLOCATE (omegas_complex)
         DEALLOCATE (omegas_pade)
         DEALLOCATE (omegas_pade_real)
         DEALLOCATE (results_pade)
      END IF

      IF (do_polarizability) THEN
         DEALLOCATE (pol_results)
         DEALLOCATE (field_results)
      END IF

      IF (do_moments_ft .OR. do_polarizability) THEN
         DEALLOCATE (results)
         DEALLOCATE (omegas)
      END IF
   END SUBROUTINE print_ft

! **************************************************************************************************
!> \brief ...
!> \param ft_unit ...
!> \param headers ...
!> \param xvals ...
!> \param yvals ...
!> \param prefix ...
!> \param prefix_format ...
!> \param xscale_opt ...
! **************************************************************************************************
   SUBROUTINE print_ft_file(ft_unit, headers, xvals, yvals, prefix, prefix_format, xscale_opt)
      INTEGER, INTENT(IN)                                :: ft_unit
      CHARACTER(len=20), DIMENSION(:), INTENT(IN)        :: headers
      REAL(kind=dp), DIMENSION(:), INTENT(IN)            :: xvals
      COMPLEX(kind=dp), DIMENSION(:, :), INTENT(IN)      :: yvals
      CHARACTER(len=21), INTENT(IN), OPTIONAL            :: prefix
      CHARACTER(len=5), INTENT(IN), OPTIONAL             :: prefix_format
      REAL(kind=dp), INTENT(IN), OPTIONAL                :: xscale_opt

      INTEGER                                            :: i, j, ncols, nrows
      LOGICAL                                            :: do_prefix
      REAL(kind=dp)                                      :: xscale

      do_prefix = .FALSE.
      IF (PRESENT(prefix)) THEN
         IF (PRESENT(prefix_format)) THEN
            do_prefix = .TRUE.
         ELSE
            CPABORT("Printing of prefix with missing format!")
         END IF
      END IF

      xscale = 1.0_dp
      IF (PRESENT(xscale_opt)) xscale = xscale_opt

      ncols = SIZE(yvals, 1)
      nrows = SIZE(yvals, 2)

      ! Check whether enough headers for yvals and xvals is present
      IF (SIZE(headers) < 2*ncols + 1) THEN
         CPABORT("Not enough headers to print the file!")
      END IF

      IF (SIZE(xvals) < nrows) THEN
         CPABORT("Not enough xvals to print all yvals!")
      END IF

      IF (ft_unit > 0) THEN
         ! If prefix is present, write prefix
         ! Also write energy in eV instead of at.u.
         IF (do_prefix) THEN
            WRITE (ft_unit, prefix_format, advance="no") prefix
         END IF
         WRITE (ft_unit, "(A20)", advance="no") headers(1)
         ! Print the rest of the headers
         DO j = 1, 2*ncols - 1
            WRITE (ft_unit, "(A20)", advance="no") headers(j + 1)
         END DO
         WRITE (ft_unit, "(A20)") headers(2*ncols + 1)
         ! Done with the headers, print actual data
         DO i = 1, nrows
            IF (do_prefix) THEN
               WRITE (ft_unit, prefix_format, advance="no") prefix
            END IF
            WRITE (ft_unit, "(E20.8E3)", advance="no") xvals(i)*xscale
            DO j = 1, ncols - 1
               WRITE (ft_unit, "(E20.8E3,E20.8E3)", advance="no") &
                  REAL(yvals(j, i)), AIMAG(yvals(j, i))
            END DO
            WRITE (ft_unit, "(E20.8E3,E20.8E3)") &
               REAL(yvals(ncols, i)), AIMAG(yvals(ncols, i))
         END DO
      END IF
   END SUBROUTINE print_ft_file

END MODULE rt_propagation_output
