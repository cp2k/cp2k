!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines needed for EMD
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

MODULE rt_propagation_utils
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_plus_fm_fm_t,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_binary_read, dbcsr_checksum, dbcsr_copy, dbcsr_copy_into_existing, &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_distribution_type, &
        dbcsr_filter, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_scale, &
        dbcsr_set, dbcsr_type
   USE input_constants,                 ONLY: use_restart_wfn,&
                                              use_rt_restart
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: zero
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_multiply,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_dftb_matrices,                ONLY: build_dftb_overlap
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type
   USE qs_mo_io,                        ONLY: read_mo_set_from_restart,&
                                              read_rt_mos_from_restart,&
                                              write_mo_set_to_output_unit
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_wfn_mix,                  ONLY: wfn_mix
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: get_restart_wfn, &
             calc_S_derivs, &
             calc_update_rho, &
             calc_update_rho_sparse, &
             calculate_P_imaginary, &
             write_rtp_mos_to_output_unit, &
             write_rtp_mo_cubes

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation_utils'

CONTAINS

! **************************************************************************************************
!> \brief Calculates dS/dR respectily the velocity weighted derivatves
!>        only needed for ehrenfest MD.
!>
!> \param qs_env the qs environment
!> \par History
!>      02.2009 created [Manuel Guidon]
!>      02.2014 switched to dbcsr matrices [Samuel Andermatt]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE calc_S_derivs(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calc_S_derivs'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: col_atom, handle, i, j, m, maxder, n, &
                                                            nder, row_atom
      INTEGER, DIMENSION(6, 2)                           :: c_map_mat
      LOGICAL                                            :: return_s_derivatives
      REAL(dp), DIMENSION(:), POINTER                    :: block_values
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: C_mat, S_der, s_derivs
      TYPE(dbcsr_type), POINTER                          :: B_mat, tmp_mat, tmp_mat2
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(rt_prop_type), POINTER                        :: rtp

      CALL timeset(routineN, handle)

      return_s_derivatives = .TRUE.

      NULLIFY (particle_set)
      NULLIFY (rtp)
      NULLIFY (s_derivs)
      NULLIFY (dft_control)
      NULLIFY (ks_env)

      CALL get_qs_env(qs_env=qs_env, &
                      rtp=rtp, &
                      particle_set=particle_set, &
                      sab_orb=sab_orb, &
                      dft_control=dft_control, &
                      ks_env=ks_env)

      CALL get_rtp(rtp=rtp, B_mat=B_mat, C_mat=C_mat, S_der=S_der)

      nder = 2
      maxder = ncoset(nder)

      NULLIFY (tmp_mat)
      ALLOCATE (tmp_mat)
      CALL dbcsr_create(tmp_mat, template=S_der(1)%matrix, matrix_type="N")

      IF (rtp%iter < 2) THEN
         ! calculate the overlap derivative matrices
         IF (dft_control%qs_control%dftb) THEN
            CALL build_dftb_overlap(qs_env, nder, s_derivs)
         ELSE
            CALL build_overlap_matrix(ks_env, nderivative=nder, matrix_s=s_derivs, &
                                      basis_type_a="ORB", basis_type_b="ORB", sab_nl=sab_orb)
         END IF

         NULLIFY (tmp_mat2)
         ALLOCATE (tmp_mat2)
         CALL dbcsr_create(tmp_mat2, template=S_der(1)%matrix, matrix_type="S")
         DO m = 1, 9
            CALL dbcsr_copy(tmp_mat2, s_derivs(m + 1)%matrix)
            CALL dbcsr_desymmetrize(tmp_mat2, S_der(m)%matrix)
            CALL dbcsr_scale(S_der(m)%matrix, -one)
            CALL dbcsr_filter(S_der(m)%matrix, rtp%filter_eps)
            !The diagonal should be zero
            CALL dbcsr_iterator_start(iter, S_der(m)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
               IF (row_atom == col_atom) block_values = 0.0_dp
            END DO
            CALL dbcsr_iterator_stop(iter)
         END DO
         CALL dbcsr_deallocate_matrix_set(s_derivs)
         CALL dbcsr_deallocate_matrix(tmp_mat2)
      END IF

      !calculate scalar product v(Rb)*<alpha|d/dRb beta> (B_mat), and store the first derivatives

      CALL dbcsr_set(B_mat, zero)
      DO m = 1, 3
         CALL dbcsr_copy(tmp_mat, S_der(m)%matrix)
         CALL dbcsr_iterator_start(iter, tmp_mat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
            IF (row_atom == col_atom) block_values = 0.0_dp
            block_values = block_values*particle_set(col_atom)%v(m)
         END DO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_add(B_mat, tmp_mat, one, one)
      END DO
      CALL dbcsr_filter(B_mat, rtp%filter_eps)
      !calculate C matrix: v(Rb)*<d/dRa alpha| d/dRb beta>

      c_map_mat = 0
      n = 0
      DO j = 1, 3
         DO m = j, 3
            n = n + 1
            c_map_mat(n, 1) = j
            IF (m == j) CYCLE
            c_map_mat(n, 2) = m
         END DO
      END DO

      DO i = 1, 3
         CALL dbcsr_set(C_mat(i)%matrix, zero)
      END DO
      DO m = 1, 6
         CALL dbcsr_copy(tmp_mat, S_der(m + 3)%matrix)
         DO j = 1, 2
            IF (c_map_mat(m, j) == 0) CYCLE
            CALL dbcsr_add(C_mat(c_map_mat(m, j))%matrix, tmp_mat, one, one)
         END DO
      END DO

      DO m = 1, 3
         CALL dbcsr_iterator_start(iter, C_mat(m)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row_atom, col_atom, block_values)
            block_values = block_values*particle_set(row_atom)%v(m)
         END DO
         CALL dbcsr_iterator_stop(iter)
         CALL dbcsr_filter(C_mat(m)%matrix, rtp%filter_eps)
      END DO

      CALL dbcsr_deallocate_matrix(tmp_mat)
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief reads the restart file. At the moment only SCF (means only real)
!> \param qs_env ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE get_restart_wfn(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=default_path_length)                 :: file_name, project_name
      INTEGER                                            :: i, id_nr, im, ispin, ncol, nspin, &
                                                            output_unit, re, unit_nr
      REAL(KIND=dp)                                      :: alpha, cs_pos
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type)                                   :: mos_occ
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new, mos_old
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_rmpv, rho_new, rho_old
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo_array
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      NULLIFY (atomic_kind_set, qs_kind_set, mo_array, particle_set, rho_struct, para_env)

      CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      mos=mo_array, &
                      input=input, &
                      rtp=rtp, &
                      dft_control=dft_control, &
                      rho=rho_struct, &
                      para_env=para_env)
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      id_nr = 0
      nspin = SIZE(mo_array)
      CALL qs_rho_get(rho_struct, rho_ao=p_rmpv)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      SELECT CASE (dft_control%rtp_control%initial_wfn)
      CASE (use_restart_wfn)
         CALL read_mo_set_from_restart(mo_array, atomic_kind_set, qs_kind_set, particle_set, para_env, &
                                       id_nr=id_nr, multiplicity=dft_control%multiplicity, dft_section=dft_section)
         CALL set_uniform_occupation_mo_array(mo_array, nspin)

         IF (dft_control%rtp_control%apply_wfn_mix_init_restart) &
            CALL wfn_mix(mo_array, particle_set, dft_section, qs_kind_set, para_env, output_unit, &
                         for_rtp=.TRUE.)

         DO ispin = 1, nspin
            CALL calculate_density_matrix(mo_array(ispin), p_rmpv(ispin)%matrix)
         END DO
         IF (rtp%linear_scaling) THEN
            CALL get_rtp(rtp=rtp, rho_old=rho_old, rho_new=rho_new)
            DO ispin = 1, nspin
               re = 2*ispin - 1
               im = 2*ispin
               CALL cp_fm_get_info(mo_array(ispin)%mo_coeff, ncol_global=ncol)
               CALL cp_fm_create(mos_occ, &
                                 matrix_struct=mo_array(ispin)%mo_coeff%matrix_struct, &
                                 name="mos_occ")
               CALL cp_fm_to_fm(mo_array(ispin)%mo_coeff, mos_occ)
               IF (mo_array(ispin)%uniform_occupation) THEN
                  alpha = 3.0_dp - REAL(nspin, dp)
                  CALL cp_fm_column_scale(mos_occ, mo_array(ispin)%occupation_numbers/alpha)
                  CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_old(re)%matrix, &
                                             matrix_v=mos_occ, &
                                             ncol=ncol, &
                                             alpha=alpha, keep_sparsity=.FALSE.)
               ELSE
                  alpha = 1.0_dp
                  CALL cp_fm_column_scale(mos_occ, mo_array(ispin)%occupation_numbers/alpha)
                  CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_old(re)%matrix, &
                                             matrix_v=mo_array(ispin)%mo_coeff, &
                                             matrix_g=mos_occ, &
                                             ncol=ncol, &
                                             alpha=alpha, keep_sparsity=.FALSE.)
               END IF
               CALL dbcsr_filter(rho_old(re)%matrix, rtp%filter_eps)
               CALL dbcsr_copy(rho_new(re)%matrix, rho_old(re)%matrix)
               CALL cp_fm_release(mos_occ)
            END DO
            CALL calc_update_rho_sparse(qs_env)
         ELSE
            CALL get_rtp(rtp=rtp, mos_old=mos_old)
            DO i = 1, SIZE(qs_env%mos)
               CALL cp_fm_to_fm(mo_array(i)%mo_coeff, mos_old(2*i - 1))
               CALL cp_fm_set_all(mos_old(2*i), zero, zero)
            END DO
         END IF
      CASE (use_rt_restart)
         IF (rtp%linear_scaling) THEN
            CALL get_rtp(rtp=rtp, rho_old=rho_old, rho_new=rho_new)
            project_name = logger%iter_info%project_name
            DO ispin = 1, nspin
               re = 2*ispin - 1
               im = 2*ispin
               WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_RE", ispin, "_RESTART.dm"
               CALL dbcsr_get_info(rho_old(re)%matrix, distribution=dist)
               CALL dbcsr_binary_read(file_name, distribution=dist, matrix_new=rho_old(re)%matrix)
               cs_pos = dbcsr_checksum(rho_old(re)%matrix, pos=.TRUE.)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ", cs_pos
               END IF
               WRITE (file_name, '(A,I0,A)') TRIM(project_name)//"_LS_DM_SPIN_IM", ispin, "_RESTART.dm"
               CALL dbcsr_get_info(rho_old(im)%matrix, distribution=dist)
               CALL dbcsr_binary_read(file_name, distribution=dist, matrix_new=rho_old(im)%matrix)
               cs_pos = dbcsr_checksum(rho_old(im)%matrix, pos=.TRUE.)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, '(T2,A,E20.8)') "Read restart DM "//TRIM(file_name)//" with checksum: ", cs_pos
               END IF
            END DO
            DO i = 1, SIZE(rho_new)
               CALL dbcsr_copy(rho_new(i)%matrix, rho_old(i)%matrix)
            END DO
            CALL calc_update_rho_sparse(qs_env)
         ELSE
            CALL get_rtp(rtp=rtp, mos_old=mos_old, mos_new=mos_new)
            CALL read_rt_mos_from_restart(mo_array, mos_old, atomic_kind_set, qs_kind_set, particle_set, para_env, &
                                          id_nr, dft_control%multiplicity, dft_section)
            CALL set_uniform_occupation_mo_array(mo_array, nspin)
            DO ispin = 1, nspin
               CALL calculate_density_matrix(mo_array(ispin), &
                                             p_rmpv(ispin)%matrix)
            END DO
         END IF
      END SELECT

   END SUBROUTINE get_restart_wfn

! **************************************************************************************************
!> \brief Set mo_array(ispin)%uniform_occupation after a restart
!> \param mo_array ...
!> \param nspin ...
!> \author Guillaume Le Breton (03.23)
! **************************************************************************************************

   SUBROUTINE set_uniform_occupation_mo_array(mo_array, nspin)

      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo_array
      INTEGER                                            :: nspin

      INTEGER                                            :: ispin, mo
      LOGICAL                                            :: is_uniform

      DO ispin = 1, nspin
         is_uniform = .TRUE.
         DO mo = 1, mo_array(ispin)%nmo
            IF (mo_array(ispin)%occupation_numbers(mo) /= 0.0 .AND. &
                mo_array(ispin)%occupation_numbers(mo) /= 1.0 .AND. &
                mo_array(ispin)%occupation_numbers(mo) /= 2.0) &
               is_uniform = .FALSE.
         END DO
         mo_array(ispin)%uniform_occupation = is_uniform
      END DO

   END SUBROUTINE set_uniform_occupation_mo_array

! **************************************************************************************************
!> \brief calculates the density from the complex MOs and passes the density to
!>        qs_env.
!> \param qs_env ...
!> \author Florian Schiffmann (02.09)
! **************************************************************************************************

   SUBROUTINE calc_update_rho(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_update_rho'
      REAL(KIND=dp), PARAMETER                           :: one = 1.0_dp, zero = 0.0_dp

      INTEGER                                            :: handle, i, im, ncol, re
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_type)                                   :: mos_occ
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_ao_im
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp

      CALL timeset(routineN, handle)

      NULLIFY (rho, ks_env, mos_new, rtp)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      rho=rho, &
                      rtp=rtp, &
                      mos=mos)
      CALL get_rtp(rtp=rtp, mos_new=mos_new)
      CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao)
      DO i = 1, SIZE(mos_new)/2
         re = 2*i - 1; im = 2*i
         CALL dbcsr_set(rho_ao(i)%matrix, zero)
         CALL cp_fm_get_info(mos_new(re), ncol_global=ncol)
         CALL cp_fm_create(mos_occ, &
                           matrix_struct=mos(i)%mo_coeff%matrix_struct, &
                           name="mos_occ")
         CALL cp_fm_to_fm(mos_new(re), mos_occ)
         IF (mos(i)%uniform_occupation) THEN
            alpha = 3*one - REAL(SIZE(mos_new)/2, dp)
            CALL cp_fm_column_scale(mos_occ, mos(i)%occupation_numbers/alpha)
            CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_ao(i)%matrix, &
                                       matrix_v=mos_occ, &
                                       ncol=ncol, &
                                       alpha=alpha)
         ELSE
            alpha = 1.0_dp
            CALL cp_fm_column_scale(mos_occ, mos(i)%occupation_numbers/alpha)
            CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_ao(i)%matrix, &
                                       matrix_v=mos_new(re), &
                                       matrix_g=mos_occ, &
                                       ncol=ncol, &
                                       alpha=alpha)
         END IF

         ! It is actually complex conjugate but i*i=-1 therefore it must be added
         CALL cp_fm_to_fm(mos_new(im), mos_occ)
         IF (mos(i)%uniform_occupation) THEN
            alpha = 3*one - REAL(SIZE(mos_new)/2, dp)
            CALL cp_fm_column_scale(mos_occ, mos(i)%occupation_numbers/alpha)
            CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_ao(i)%matrix, &
                                       matrix_v=mos_occ, &
                                       ncol=ncol, &
                                       alpha=alpha)
         ELSE
            alpha = 1.0_dp
            CALL cp_fm_column_scale(mos_occ, mos(i)%occupation_numbers/alpha)
            CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=rho_ao(i)%matrix, &
                                       matrix_v=mos_new(im), &
                                       matrix_g=mos_occ, &
                                       ncol=ncol, &
                                       alpha=alpha)
         END IF
         CALL cp_fm_release(mos_occ)
      END DO
      CALL qs_rho_update_rho(rho, qs_env)

      IF (rtp%track_imag_density) THEN
         CALL qs_rho_get(rho_struct=rho, rho_ao_im=rho_ao_im)
         CALL calculate_P_imaginary(qs_env, rtp, rho_ao_im, keep_sparsity=.TRUE.)
         CALL qs_rho_set(rho, rho_ao_im=rho_ao_im)
      END IF

      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE.)

      CALL timestop(handle)

   END SUBROUTINE calc_update_rho

! **************************************************************************************************
!> \brief Copies the density matrix back into the qs_env%rho%rho_ao
!> \param qs_env ...
!> \author Samuel Andermatt (3.14)
! **************************************************************************************************

   SUBROUTINE calc_update_rho_sparse(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_update_rho_sparse'
      REAL(KIND=dp), PARAMETER                           :: zero = 0.0_dp

      INTEGER                                            :: handle, ispin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_ao_im, rho_new
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(rtp_control_type), POINTER                    :: rtp_control

      NULLIFY (rho, ks_env, rtp, dft_control)
      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      rho=rho, &
                      rtp=rtp, &
                      dft_control=dft_control)
      rtp_control => dft_control%rtp_control
      CALL get_rtp(rtp=rtp, rho_new=rho_new)
      CALL qs_rho_get(rho_struct=rho, rho_ao=rho_ao)
      IF (rtp%track_imag_density) CALL qs_rho_get(rho_struct=rho, rho_ao_im=rho_ao_im)
      DO ispin = 1, SIZE(rho_ao)
         CALL dbcsr_set(rho_ao(ispin)%matrix, zero)
         CALL dbcsr_copy_into_existing(rho_ao(ispin)%matrix, rho_new(ispin*2 - 1)%matrix)
         IF (rtp%track_imag_density) THEN
            CALL dbcsr_copy_into_existing(rho_ao_im(ispin)%matrix, rho_new(ispin*2)%matrix)
         END IF
      END DO

      CALL qs_rho_update_rho(rho, qs_env)
      CALL qs_ks_did_change(ks_env, rho_changed=.TRUE.)

      CALL timestop(handle)

   END SUBROUTINE calc_update_rho_sparse

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rtp ...
!> \param matrix_p_im ...
!> \param keep_sparsity ...
! **************************************************************************************************
   SUBROUTINE calculate_P_imaginary(qs_env, rtp, matrix_p_im, keep_sparsity)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rt_prop_type), POINTER                        :: rtp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p_im
      LOGICAL, OPTIONAL                                  :: keep_sparsity

      INTEGER                                            :: i, im, ncol, re
      LOGICAL                                            :: my_keep_sparsity
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new, mos_occ
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL get_rtp(rtp=rtp, mos_new=mos_new)

      my_keep_sparsity = .FALSE.
      IF (PRESENT(keep_sparsity)) my_keep_sparsity = keep_sparsity
      CALL get_qs_env(qs_env, mos=mos)
      ALLOCATE (mos_occ(SIZE(mos)*2))

      DO i = 1, SIZE(mos_new)/2
         re = 2*i - 1; im = 2*i
         alpha = 3.0_dp - REAL(SIZE(matrix_p_im), dp)
         CALL cp_fm_create(mos_occ(re), &
                           matrix_struct=mos(i)%mo_coeff%matrix_struct, &
                           name="mos_occ")
         CALL cp_fm_create(mos_occ(im), &
                           matrix_struct=mos(i)%mo_coeff%matrix_struct, &
                           name="mos_occ")
         CALL dbcsr_set(matrix_p_im(i)%matrix, 0.0_dp)
         CALL cp_fm_get_info(mos_new(re), ncol_global=ncol)
         CALL cp_fm_to_fm(mos_new(re), mos_occ(re))
         CALL cp_fm_column_scale(mos_occ(re), mos(i)%occupation_numbers/alpha)
         CALL cp_fm_to_fm(mos_new(im), mos_occ(im))
         CALL cp_fm_column_scale(mos_occ(im), mos(i)%occupation_numbers/alpha)
         CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=matrix_p_im(i)%matrix, &
                                    matrix_v=mos_occ(im), &
                                    matrix_g=mos_occ(re), &
                                    ncol=ncol, &
                                    keep_sparsity=my_keep_sparsity, &
                                    alpha=2.0_dp*alpha, &
                                    symmetry_mode=-1)
      END DO
      CALL cp_fm_release(mos_occ)

   END SUBROUTINE calculate_P_imaginary

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rtp ...
! **************************************************************************************************
   SUBROUTINE write_rtp_mos_to_output_unit(qs_env, rtp)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rt_prop_type), POINTER                        :: rtp

      CHARACTER(len=*), PARAMETER :: routineN = 'write_rtp_mos_to_output_unit'

      CHARACTER(LEN=10)                                  :: spin
      CHARACTER(LEN=2*default_string_length)             :: name
      INTEGER                                            :: handle, i, ispin, nao, nelectron, nmo, &
                                                            nspins
      LOGICAL                                            :: print_eigvecs, print_mo_info
      REAL(KIND=dp)                                      :: flexible_electron_count, maxocc, n_el_f
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_set_type)                                  :: mo_set_rtp
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      NULLIFY (atomic_kind_set, particle_set, qs_kind_set, input, mos, dft_section)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      input=input, &
                      mos=mos)
      ! Quick return, if no printout of MO information is requested
      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL section_vals_val_get(dft_section, "PRINT%MO%EIGENVECTORS", l_val=print_eigvecs)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      print_mo_info = (cp_print_key_should_output(logger%iter_info, &
                                                  dft_section, "PRINT%MO") /= 0) .OR. &
                      (qs_env%sim_step == 1)

      IF ((.NOT. print_mo_info) .OR. (.NOT. print_eigvecs)) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL get_rtp(rtp=rtp, mos_new=mos_new)

      nspins = SIZE(mos_new)/2

      DO ispin = 1, nspins
         ! initiate mo_set
         CALL get_mo_set(mo_set=mos(ispin), nao=nao, nmo=nmo, nelectron=nelectron, &
                         n_el_f=n_el_f, maxocc=maxocc, flexible_electron_count=flexible_electron_count)

         CALL allocate_mo_set(mo_set_rtp, &
                              nao=nao, &
                              nmo=nmo, &
                              nelectron=nelectron, &
                              n_el_f=n_el_f, &
                              maxocc=maxocc, &
                              flexible_electron_count=flexible_electron_count)

         WRITE (name, FMT="(A,I6)") "RTP MO SET, SPIN ", ispin
         CALL init_mo_set(mo_set_rtp, fm_ref=mos_new(2*ispin - 1), name=name)

         IF (nspins > 1) THEN
            IF (ispin == 1) THEN
               spin = "ALPHA SPIN"
            ELSE
               spin = "BETA SPIN"
            END IF
         ELSE
            spin = ""
         END IF

         mo_set_rtp%occupation_numbers = mos(ispin)%occupation_numbers

         !loop for real (odd) and imaginary (even) parts
         DO i = 1, 0, -1
            CALL cp_fm_to_fm(mos_new(2*ispin - i), mo_set_rtp%mo_coeff)
            CALL write_mo_set_to_output_unit(mo_set_rtp, atomic_kind_set, qs_kind_set, particle_set, &
                                          dft_section, 4, 0, rtp=.TRUE., spin=TRIM(spin), cpart=MOD(i, 2), sim_step=qs_env%sim_step)
         END DO

         CALL deallocate_mo_set(mo_set_rtp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE write_rtp_mos_to_output_unit

! **************************************************************************************************
!> \brief Write the time dependent amplitude of the MOs in real grid.
!>        Very close to qs_scf_post_gpw/qs_scf_post_occ_cubes subroutine.
!> \param qs_env ...
!> \param rtp ...
!> \author Guillaume Le Breton (11.22)
! **************************************************************************************************
   SUBROUTINE write_rtp_mo_cubes(qs_env, rtp)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rt_prop_type), POINTER                        :: rtp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_rtp_mo_cubes'

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube, title
      INTEGER                                            :: handle, homo, i, ir, ispin, ivector, &
                                                            n_rep, nhomo, nlist, nspins, &
                                                            rt_time_step, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: list, list_index
      LOGICAL                                            :: append_cube, do_kpoints, mpi_io
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: density_r, wf_g, wf_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      CALL timeset(routineN, handle)

      NULLIFY (logger, auxbas_pw_pool, pw_pools, pw_env)

      ! Get all the info from qs:
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints, &
                      input=input)

      ! Kill the run in the case of K points
      IF (do_kpoints) THEN
         CPABORT("K points not handled yet for printing MO_CUBE")
      END IF

      dft_section => section_vals_get_subs_vals(input, "DFT")
      logger => cp_get_default_logger()

      ! Quick return if no print required
      IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, dft_section, &
                                                 "PRINT%MO_CUBES"), cp_p_file)) THEN
         CALL timestop(handle)
         RETURN
      END IF

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, &
                      mos=mos, &
                      blacs_env=blacs_env, &
                      qs_kind_set=qs_kind_set, &
                      pw_env=pw_env, &
                      subsys=subsys, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      dft_control=dft_control)
      CALL qs_subsys_get(subsys, particles=particles)

      nspins = dft_control%nspins
      rt_time_step = qs_env%sim_step

      ! Setup the grids needed to compute a wavefunction given a vector
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools)
      CALL auxbas_pw_pool%create_pw(wf_r, &
                                    use_data=REALDATA3D, &
                                    in_space=REALSPACE)
      CALL auxbas_pw_pool%create_pw(wf_g, &
                                    use_data=COMPLEXDATA1D, &
                                    in_space=RECIPROCALSPACE)
      CALL auxbas_pw_pool%create_pw(density_r, &
                                    use_data=REALDATA3D, &
                                    in_space=REALSPACE)
      CALL get_rtp(rtp=rtp, mos_new=mos_new)

      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), homo=homo)

         nhomo = section_get_ival(dft_section, "PRINT%MO_CUBES%NHOMO")
         append_cube = section_get_lval(dft_section, "PRINT%MO_CUBES%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%HOMO_LIST", n_rep_val=n_rep)
         IF (n_rep > 0) THEN ! write the cubes of the list
            nlist = 0
            DO ir = 1, n_rep
               NULLIFY (list)
               CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%HOMO_LIST", i_rep_val=ir, &
                                         i_vals=list)
               IF (ASSOCIATED(list)) THEN
                  CALL reallocate(list_index, 1, nlist + SIZE(list))
                  DO i = 1, SIZE(list)
                     list_index(i + nlist) = list(i)
                  END DO
                  nlist = nlist + SIZE(list)
               END IF
            END DO
         ELSE

            IF (nhomo == -1) nhomo = homo
            nlist = homo - MAX(1, homo - nhomo + 1) + 1
            ALLOCATE (list_index(nlist))
            DO i = 1, nlist
               list_index(i) = MAX(1, homo - nhomo + 1) + i - 1
            END DO
         END IF
         DO i = 1, nlist
            ivector = list_index(i)
            CALL get_qs_env(qs_env=qs_env, &
                            atomic_kind_set=atomic_kind_set, &
                            qs_kind_set=qs_kind_set, &
                            cell=cell, &
                            particle_set=particle_set, &
                            pw_env=pw_env)

            ! density_r contains the density of the MOs
            CALL pw_zero(density_r)
            mo_coeff => mos_new(2*ispin - 1)!Real coeff
            CALL calculate_wavefunction(mo_coeff, ivector, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                        cell, dft_control, particle_set, pw_env)
            ! Adding the real part
            CALL pw_multiply(density_r, wf_r, wf_r, 1.0_dp)

            mo_coeff => mos_new(2*ispin) !Im coeff
            CALL calculate_wavefunction(mo_coeff, ivector, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                        cell, dft_control, particle_set, pw_env)
            ! Adding the im part
            CALL pw_multiply(density_r, wf_r, wf_r, 1.0_dp)

            WRITE (filename, '(a4,I5.5,a1,I1.1)') "WFN_", ivector, "_", ispin
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MO_CUBES", extension=".cube", &
                                           middle_name=TRIM(filename), file_position=my_pos_cube, log_filename=.FALSE., &
                                           mpi_io=mpi_io)
            WRITE (title, *) "DENSITY ", ivector, " spin ", ispin, " i.e. HOMO - ", ivector - homo
            CALL cp_pw_to_cube(density_r, unit_nr, title, particles=particles, &
                               stride=section_get_ivals(dft_section, "PRINT%MO_CUBES%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%MO_CUBES", mpi_io=mpi_io)
         END DO
         IF (ASSOCIATED(list_index)) DEALLOCATE (list_index)
      END DO

      ! Deallocate grids needed to compute wavefunctions
      CALL auxbas_pw_pool%give_back_pw(wf_r)
      CALL auxbas_pw_pool%give_back_pw(wf_g)
      CALL auxbas_pw_pool%give_back_pw(density_r)

      CALL timestop(handle)

   END SUBROUTINE write_rtp_mo_cubes

END MODULE rt_propagation_utils
