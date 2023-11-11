!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines needed for cubic-scaling RPA and SOS-Laplace-MP2 forces
!> \author Augustin Bussy
! **************************************************************************************************
MODULE rpa_im_time_force_methods
   USE admm_methods,                    ONLY: admm_projection_derivative
   USE admm_types,                      ONLY: admm_type,&
                                              get_admm_env
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_power
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_dist2d_to_dist,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_eri_mme_interface,            ONLY: cp_eri_mme_update_local_counts
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_clear, dbcsr_complete_redistribute, dbcsr_copy, &
        dbcsr_create, dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_frobenius_norm, dbcsr_get_block_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE dbt_api,                         ONLY: &
        dbt_batched_contract_finalize, dbt_batched_contract_init, dbt_clear, dbt_contract, &
        dbt_copy, dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, dbt_create, dbt_destroy, &
        dbt_filter, dbt_get_info, dbt_mp_environ_pgrid, dbt_pgrid_create, dbt_pgrid_destroy, &
        dbt_pgrid_type, dbt_scale, dbt_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE ec_methods,                      ONLY: create_kernel
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel
   USE hfx_admm_utils,                  ONLY: tddft_hfx_matrix
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_exx,                         ONLY: add_exx_to_rhs
   USE hfx_ri,                          ONLY: get_2c_der_force,&
                                              get_force_from_3c_trace,&
                                              get_idx_to_atom,&
                                              hfx_ri_update_forces
   USE hfx_types,                       ONLY: alloc_containers,&
                                              block_ind_type,&
                                              dealloc_containers,&
                                              hfx_compression_type,&
                                              hfx_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              do_eri_gpw,&
                                              do_eri_mme,&
                                              do_potential_id,&
                                              ri_rpa_method_gpw
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iterate_matrix,                  ONLY: matrix_exponential
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE mathconstants,                   ONLY: fourpi
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE mp2_eri,                         ONLY: integrate_set_2c
   USE mp2_eri_gpw,                     ONLY: calc_potential_gpw,&
                                              cleanup_gpw,&
                                              prepare_gpw,&
                                              virial_gpw_potential
   USE mp2_types,                       ONLY: mp2_type
   USE orbital_pointers,                ONLY: ncoset
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integral_ab,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                              calculate_wavefunction
   USE qs_density_matrices,             ONLY: calculate_whz_matrix
   USE qs_energy_matrix_w,              ONLY: qs_energies_compute_w
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_integrate_potential,          ONLY: integrate_pgf_product,&
                                              integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_reference,                 ONLY: ks_ref_potential
   USE qs_ks_types,                     ONLY: set_ks_env
   USE qs_linres_types,                 ONLY: linres_control_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_tensors,                      ONLY: &
        build_2c_derivatives, build_2c_integrals, build_2c_neighbor_lists, build_3c_derivatives, &
        build_3c_neighbor_lists, calc_2c_virial, calc_3c_virial, compress_tensor, &
        decompress_tensor, get_tensor_occupancy, neighbor_list_3c_destroy
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              create_tensor_batches,&
                                              distribution_3d_create,&
                                              distribution_3d_type,&
                                              neighbor_list_3c_type
   USE realspace_grid_types,            ONLY: map_gaussian_here,&
                                              realspace_grid_type
   USE response_solver,                 ONLY: response_equation_new
   USE rpa_im_time,                     ONLY: compute_mat_dm_global
   USE rpa_im_time_force_types,         ONLY: im_time_force_type
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE task_list_types,                 ONLY: task_list_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_im_time_force_methods'

   PUBLIC :: init_im_time_forces, calc_laplace_loop_forces, calc_post_loop_forces, &
             keep_initial_quad, calc_rpa_loop_forces

CONTAINS

! **************************************************************************************************
!> \brief Initializes and pre-calculates all needed tensors for the forces
!> \param force_data ...
!> \param fm_matrix_PQ ...
!> \param t_3c_M the 3-center M tensor to be used as a template
!> \param unit_nr ...
!> \param mp2_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE init_im_time_forces(force_data, fm_matrix_PQ, t_3c_M, unit_nr, mp2_env, qs_env)

      TYPE(im_time_force_type), INTENT(INOUT)            :: force_data
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_matrix_PQ
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_M
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'init_im_time_forces'

      INTEGER                                            :: handle, i_mem, i_xyz, ibasis, ispin, &
                                                            n_dependent, n_mem, n_rep, natom, &
                                                            nkind, nspins
      INTEGER(int_8)                                     :: nze, nze_tot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1, dist2, dist_AO_1, dist_AO_2, &
                                                            dist_RI, dummy_end, dummy_start, &
                                                            end_blocks, sizes_AO, sizes_RI, &
                                                            start_blocks
      INTEGER, DIMENSION(2)                              :: pdims_t2c
      INTEGER, DIMENSION(3)                              :: nblks_total, pcoord, pdims, pdims_t3c
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      LOGICAL                                            :: do_periodic, use_virial
      REAL(dp)                                           :: compression_factor, eps_pgf_orb, &
                                                            eps_pgf_orb_old, memory, occ
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dbcsr_type)                                   :: dbcsr_work, dbcsr_work2, dbcsr_work3
      TYPE(dbcsr_type), DIMENSION(1)                     :: t_2c_int_tmp
      TYPE(dbcsr_type), DIMENSION(1, 3)                  :: t_2c_der_tmp
      TYPE(dbt_pgrid_type)                               :: pgrid_t2c, pgrid_t3c
      TYPE(dbt_type)                                     :: t_2c_template, t_2c_tmp, t_3c_template
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:, :, :)    :: t_3c_der_AO_prv, t_3c_der_RI_prv
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d, dist_vir
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_ao, basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis
      TYPE(libint_potential_type)                        :: identity_pot
      TYPE(mp_cart_type)                                 :: mp_comm_t3c, mp_comm_vir
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: qs_section
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (dft_control, para_env, particle_set, qs_kind_set, dist_2d, nl_2c, blacs_env, matrix_s, &
               rho, rho_ao, cell, qs_section, orb_basis, ri_basis, virial)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, natom=natom, nkind=nkind, dft_control=dft_control, para_env=para_env, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set, cell=cell, virial=virial)
      IF (dft_control%qs_control%gapw) THEN
         CPABORT("Low-scaling RPA/SOS-MP2 forces only available with GPW")
      END IF

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      do_periodic = .FALSE.
      IF (ANY(cell%perd == 1)) do_periodic = .TRUE.
      force_data%do_periodic = do_periodic

      !Dealing with the 3-center derivatives
      pdims_t3c = 0
      CALL dbt_pgrid_create(para_env, pdims_t3c, pgrid_t3c)

      !Make sure we use the proper QS EPS_PGF_ORB values
      qs_section => section_vals_get_subs_vals(qs_env%input, "DFT%QS")
      CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", n_rep_val=n_rep)
      IF (n_rep /= 0) THEN
         CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", r_val=eps_pgf_orb)
      ELSE
         CALL section_vals_val_get(qs_section, "EPS_DEFAULT", r_val=eps_pgf_orb)
         eps_pgf_orb = SQRT(eps_pgf_orb)
      END IF
      eps_pgf_orb_old = dft_control%qs_control%eps_pgf_orb

      ALLOCATE (sizes_RI(natom), sizes_AO(natom))
      ALLOCATE (basis_set_ri_aux(nkind), basis_set_ao(nkind))
      CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_ri_aux)
      CALL basis_set_list_setup(basis_set_ao, "ORB", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_ao)

      DO ibasis = 1, SIZE(basis_set_ao)
         orb_basis => basis_set_ao(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(orb_basis, eps_pgf_orb)
         ri_basis => basis_set_ri_aux(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(ri_basis, eps_pgf_orb)
      END DO

      CALL create_3c_tensor(t_3c_template, dist_RI, dist_AO_1, dist_AO_2, pgrid_t3c, &
                            sizes_RI, sizes_AO, sizes_AO, map1=[1], map2=[2, 3], name="der (RI AO | AO)")

      ALLOCATE (t_3c_der_RI_prv(1, 1, 3), t_3c_der_AO_prv(1, 1, 3))
      DO i_xyz = 1, 3
         CALL dbt_create(t_3c_template, t_3c_der_RI_prv(1, 1, i_xyz))
         CALL dbt_create(t_3c_template, t_3c_der_AO_prv(1, 1, i_xyz))
      END DO

      IF (use_virial) THEN
         ALLOCATE (force_data%t_3c_virial, force_data%t_3c_virial_split)
         CALL dbt_create(t_3c_template, force_data%t_3c_virial)
         CALL dbt_create(t_3c_M, force_data%t_3c_virial_split)
      END IF
      CALL dbt_destroy(t_3c_template)

      CALL dbt_mp_environ_pgrid(pgrid_t3c, pdims, pcoord)
      CALL mp_comm_t3c%create(pgrid_t3c%mp_comm_2d, 3, pdims)
      CALL distribution_3d_create(dist_3d, dist_RI, dist_AO_1, dist_AO_2, &
                                  nkind, particle_set, mp_comm_t3c, own_comm=.TRUE.)

      !In case of virial, we need to store the 3c_nl
      IF (use_virial) THEN
         ALLOCATE (force_data%nl_3c)
         CALL mp_comm_vir%create(pgrid_t3c%mp_comm_2d, 3, pdims)
         CALL distribution_3d_create(dist_vir, dist_RI, dist_AO_1, dist_AO_2, &
                                     nkind, particle_set, mp_comm_vir, own_comm=.TRUE.)
         CALL build_3c_neighbor_lists(force_data%nl_3c, basis_set_ri_aux, basis_set_ao, basis_set_ao, &
                                      dist_vir, mp2_env%ri_metric, "RPA_3c_nl", qs_env, op_pos=1, &
                                      sym_jk=.FALSE., own_dist=.TRUE.)
      END IF

      CALL build_3c_neighbor_lists(nl_3c, basis_set_ri_aux, basis_set_ao, basis_set_ao, dist_3d, &
                                   mp2_env%ri_metric, "RPA_3c_nl", qs_env, op_pos=1, sym_jk=.TRUE., &
                                   own_dist=.TRUE.)
      DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)

      !Prepare the resulting 3c tensors in the format of t_3c_M for compatible traces: (RI|AO AO), split blocks
      CALL dbt_get_info(t_3c_M, nblks_total=nblks_total)
      ALLOCATE (force_data%bsizes_RI_split(nblks_total(1)), force_data%bsizes_AO_split(nblks_total(2)))
      CALL dbt_get_info(t_3c_M, blk_size_1=force_data%bsizes_RI_split, blk_size_2=force_data%bsizes_AO_split)
      DO i_xyz = 1, 3
         CALL dbt_create(t_3c_M, force_data%t_3c_der_RI(i_xyz))
         CALL dbt_create(t_3c_M, force_data%t_3c_der_AO(i_xyz))
      END DO

      !Keep track of atom index corresponding to split blocks
      ALLOCATE (force_data%idx_to_at_RI(nblks_total(1)))
      CALL get_idx_to_atom(force_data%idx_to_at_RI, force_data%bsizes_RI_split, sizes_RI)

      ALLOCATE (force_data%idx_to_at_AO(nblks_total(2)))
      CALL get_idx_to_atom(force_data%idx_to_at_AO, force_data%bsizes_AO_split, sizes_AO)

      n_mem = mp2_env%ri_rpa_im_time%cut_memory
      CALL create_tensor_batches(sizes_RI, n_mem, dummy_start, dummy_end, start_blocks, end_blocks)
      DEALLOCATE (dummy_start, dummy_end)

      ALLOCATE (force_data%t_3c_der_AO_comp(n_mem, 3), force_data%t_3c_der_RI_comp(n_mem, 3))
      ALLOCATE (force_data%t_3c_der_AO_ind(n_mem, 3), force_data%t_3c_der_RI_ind(n_mem, 3))

      memory = 0.0_dp
      nze_tot = 0
      DO i_mem = 1, n_mem
         CALL build_3c_derivatives(t_3c_der_RI_prv, t_3c_der_AO_prv, mp2_env%ri_rpa_im_time%eps_filter, &
                                   qs_env, nl_3c, basis_set_ri_aux, basis_set_ao, basis_set_ao, &
                                   mp2_env%ri_metric, der_eps=mp2_env%ri_rpa_im_time%eps_filter, op_pos=1, &
                                   bounds_i=[start_blocks(i_mem), end_blocks(i_mem)])

         DO i_xyz = 1, 3
            CALL dbt_copy(t_3c_der_RI_prv(1, 1, i_xyz), force_data%t_3c_der_RI(i_xyz), move_data=.TRUE.)
            CALL dbt_filter(force_data%t_3c_der_RI(i_xyz), mp2_env%ri_rpa_im_time%eps_filter)
            CALL get_tensor_occupancy(force_data%t_3c_der_RI(i_xyz), nze, occ)
            nze_tot = nze_tot + nze

            CALL alloc_containers(force_data%t_3c_der_RI_comp(i_mem, i_xyz), 1)
            CALL compress_tensor(force_data%t_3c_der_RI(i_xyz), force_data%t_3c_der_RI_ind(i_mem, i_xyz)%ind, &
                                 force_data%t_3c_der_RI_comp(i_mem, i_xyz), mp2_env%ri_rpa_im_time%eps_compress, memory)
            CALL dbt_clear(force_data%t_3c_der_RI(i_xyz))

            CALL dbt_copy(t_3c_der_AO_prv(1, 1, i_xyz), force_data%t_3c_der_AO(i_xyz), move_data=.TRUE.)
            CALL dbt_filter(force_data%t_3c_der_AO(i_xyz), mp2_env%ri_rpa_im_time%eps_filter)
            CALL get_tensor_occupancy(force_data%t_3c_der_AO(i_xyz), nze, occ)
            nze_tot = nze_tot + nze

            CALL alloc_containers(force_data%t_3c_der_AO_comp(i_mem, i_xyz), 1)
            CALL compress_tensor(force_data%t_3c_der_AO(i_xyz), force_data%t_3c_der_AO_ind(i_mem, i_xyz)%ind, &
                                 force_data%t_3c_der_AO_comp(i_mem, i_xyz), mp2_env%ri_rpa_im_time%eps_compress, memory)
            CALL dbt_clear(force_data%t_3c_der_AO(i_xyz))
         END DO
      END DO
      CALL neighbor_list_3c_destroy(nl_3c)
      DO i_xyz = 1, 3
         CALL dbt_destroy(t_3c_der_RI_prv(1, 1, i_xyz))
         CALL dbt_destroy(t_3c_der_AO_prv(1, 1, i_xyz))
      END DO

      CALL para_env%sum(memory)
      compression_factor = REAL(nze_tot, dp)*1.0E-06*8.0_dp/memory
      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="((T3,A,T66,F11.2,A4))") &
            "MEMORY_INFO| Memory for 3-center derivatives (compressed):", memory, ' MiB'

         WRITE (UNIT=unit_nr, FMT="((T3,A,T60,F21.2))") &
            "MEMORY_INFO| Compression factor:                  ", compression_factor
      END IF

      !Dealing with the 2-center derivatives
      CALL get_qs_env(qs_env, distribution_2d=dist_2d, blacs_env=blacs_env, matrix_s=matrix_s)
      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(SIZE(sizes_RI)))
      ALLOCATE (col_bsize(SIZE(sizes_RI)))
      row_bsize(:) = sizes_RI(:)
      col_bsize(:) = sizes_RI(:)

      pdims_t2c = 0
      CALL dbt_pgrid_create(para_env, pdims_t2c, pgrid_t2c)
      CALL create_2c_tensor(t_2c_template, dist1, dist2, pgrid_t2c, force_data%bsizes_RI_split, &
                            force_data%bsizes_RI_split, name='(RI| RI)')
      DEALLOCATE (dist1, dist2)

      CALL dbcsr_create(t_2c_int_tmp(1), "(P|Q) RPA", dbcsr_dist, dbcsr_type_symmetric, row_bsize, col_bsize)
      DO i_xyz = 1, 3
         CALL dbcsr_create(t_2c_der_tmp(1, i_xyz), "(P|Q) RPA der", dbcsr_dist, &
                           dbcsr_type_antisymmetric, row_bsize, col_bsize)
      END DO

      IF (use_virial) THEN
         ALLOCATE (force_data%RI_virial_pot, force_data%RI_virial_met)
         CALL dbcsr_create(force_data%RI_virial_pot, "RI_virial", dbcsr_dist, &
                           dbcsr_type_no_symmetry, row_bsize, col_bsize)
         CALL dbcsr_create(force_data%RI_virial_met, "RI_virial", dbcsr_dist, &
                           dbcsr_type_no_symmetry, row_bsize, col_bsize)
      END IF

      ! Main (P|Q) integrals and derivatives
      ! Integrals are passed as a full matrix => convert to DBCSR
      CALL dbcsr_create(dbcsr_work, template=t_2c_int_tmp(1))
      CALL copy_fm_to_dbcsr(fm_matrix_PQ, dbcsr_work)

      ! We need the  +/- square root of (P|Q)
      CALL dbcsr_create(dbcsr_work2, template=t_2c_int_tmp(1))
      CALL dbcsr_create(dbcsr_work3, template=t_2c_int_tmp(1))
      CALL dbcsr_copy(dbcsr_work2, dbcsr_work)
      CALL cp_dbcsr_power(dbcsr_work, -0.5_dp, 1.0E-7_dp, n_dependent, para_env, blacs_env) !1.0E-7 ev qunenching thresh

      ! Transfer to tensor format with split blocks
      CALL dbt_create(dbcsr_work, t_2c_tmp)
      CALL dbt_copy_matrix_to_tensor(dbcsr_work, t_2c_tmp)
      CALL dbt_create(t_2c_template, force_data%t_2c_pot_msqrt)
      CALL dbt_copy(t_2c_tmp, force_data%t_2c_pot_msqrt, move_data=.TRUE.)
      CALL dbt_filter(force_data%t_2c_pot_msqrt, mp2_env%ri_rpa_im_time%eps_filter)

      CALL dbcsr_multiply('N', 'N', 1.0_dp, dbcsr_work2, dbcsr_work, 0.0_dp, dbcsr_work3)
      CALL dbt_copy_matrix_to_tensor(dbcsr_work3, t_2c_tmp)
      CALL dbt_create(t_2c_template, force_data%t_2c_pot_psqrt)
      CALL dbt_copy(t_2c_tmp, force_data%t_2c_pot_psqrt, move_data=.TRUE.)
      CALL dbt_filter(force_data%t_2c_pot_psqrt, mp2_env%ri_rpa_im_time%eps_filter)
      CALL dbt_destroy(t_2c_tmp)
      CALL dbcsr_release(dbcsr_work2)
      CALL dbcsr_release(dbcsr_work3)
      CALL dbcsr_clear(dbcsr_work)

      ! Deal with the 2c potential derivatives. Only precompute if not in PBCs
      IF (.NOT. do_periodic) THEN
         CALL build_2c_neighbor_lists(nl_2c, basis_set_ri_aux, basis_set_ri_aux, mp2_env%potential_parameter, &
                                      "RPA_2c_nl_pot", qs_env, sym_ij=.TRUE., dist_2d=dist_2d)
         CALL build_2c_derivatives(t_2c_der_tmp, mp2_env%ri_rpa_im_time%eps_filter, qs_env, nl_2c, &
                                   basis_set_ri_aux, basis_set_ri_aux, mp2_env%potential_parameter)
         CALL release_neighbor_list_sets(nl_2c)

         DO i_xyz = 1, 3
            CALL dbt_create(t_2c_der_tmp(1, i_xyz), t_2c_tmp)
            CALL dbt_copy_matrix_to_tensor(t_2c_der_tmp(1, i_xyz), t_2c_tmp)
            CALL dbt_create(t_2c_template, force_data%t_2c_der_pot(i_xyz))
            CALL dbt_copy(t_2c_tmp, force_data%t_2c_der_pot(i_xyz), move_data=.TRUE.)
            CALL dbt_filter(force_data%t_2c_der_pot(i_xyz), mp2_env%ri_rpa_im_time%eps_filter)
            CALL dbt_destroy(t_2c_tmp)
            CALL dbcsr_clear(t_2c_der_tmp(1, i_xyz))
         END DO

         IF (use_virial) THEN
            CALL build_2c_neighbor_lists(force_data%nl_2c_pot, basis_set_ri_aux, basis_set_ri_aux, &
                                         mp2_env%potential_parameter, "RPA_2c_nl_pot", qs_env, &
                                         sym_ij=.FALSE., dist_2d=dist_2d)
         END IF
      END IF
      ! Create a G_PQ matrix to collect the terms for the force trace in the periodic case
      CALL dbcsr_create(force_data%G_PQ, "G_PQ", dbcsr_dist, dbcsr_type_no_symmetry, row_bsize, col_bsize)

      ! we need the RI metric derivatives and the inverse of the integrals
      CALL build_2c_neighbor_lists(nl_2c, basis_set_ri_aux, basis_set_ri_aux, mp2_env%ri_metric, &
                                   "RPA_2c_nl_metric", qs_env, sym_ij=.TRUE., dist_2d=dist_2d)
      CALL build_2c_integrals(t_2c_int_tmp, mp2_env%ri_rpa_im_time%eps_filter, qs_env, nl_2c, &
                              basis_set_ri_aux, basis_set_ri_aux, mp2_env%ri_metric)
      CALL build_2c_derivatives(t_2c_der_tmp, mp2_env%ri_rpa_im_time%eps_filter, qs_env, nl_2c, &
                                basis_set_ri_aux, basis_set_ri_aux, mp2_env%ri_metric)
      CALL release_neighbor_list_sets(nl_2c)

      IF (use_virial) THEN
         CALL build_2c_neighbor_lists(force_data%nl_2c_met, basis_set_ri_aux, basis_set_ri_aux, &
                                      mp2_env%ri_metric, "RPA_2c_nl_metric", qs_env, sym_ij=.FALSE., &
                                      dist_2d=dist_2d)
      END IF

      CALL dbcsr_copy(dbcsr_work, t_2c_int_tmp(1))
      CALL cp_dbcsr_cholesky_decompose(dbcsr_work, para_env=para_env, blacs_env=blacs_env)
      CALL cp_dbcsr_cholesky_invert(dbcsr_work, para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)

      CALL dbt_create(dbcsr_work, t_2c_tmp)
      CALL dbt_copy_matrix_to_tensor(dbcsr_work, t_2c_tmp)
      CALL dbt_create(t_2c_template, force_data%t_2c_inv_metric)
      CALL dbt_copy(t_2c_tmp, force_data%t_2c_inv_metric, move_data=.TRUE.)
      CALL dbt_filter(force_data%t_2c_inv_metric, mp2_env%ri_rpa_im_time%eps_filter)
      CALL dbt_destroy(t_2c_tmp)
      CALL dbcsr_clear(dbcsr_work)
      CALL dbcsr_clear(t_2c_int_tmp(1))

      DO i_xyz = 1, 3
         CALL dbt_create(t_2c_der_tmp(1, i_xyz), t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(t_2c_der_tmp(1, i_xyz), t_2c_tmp)
         CALL dbt_create(t_2c_template, force_data%t_2c_der_metric(i_xyz))
         CALL dbt_copy(t_2c_tmp, force_data%t_2c_der_metric(i_xyz), move_data=.TRUE.)
         CALL dbt_filter(force_data%t_2c_der_metric(i_xyz), mp2_env%ri_rpa_im_time%eps_filter)
         CALL dbt_destroy(t_2c_tmp)
         CALL dbcsr_clear(t_2c_der_tmp(1, i_xyz))
      END DO

      !Pre-calculate matrix K = metric^-1 * V^0.5
      CALL dbt_create(t_2c_template, force_data%t_2c_K)
      CALL dbt_contract(1.0_dp, force_data%t_2c_inv_metric, force_data%t_2c_pot_psqrt, &
                        0.0_dp, force_data%t_2c_K, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[1], notcontract_2=[2], &
                        map_1=[1], map_2=[2], filter_eps=mp2_env%ri_rpa_im_time%eps_filter)

      ! Finally, we need the overlap matrix derivative and the inverse of the integrals
      CALL dbt_destroy(t_2c_template)
      CALL dbcsr_release(dbcsr_work)
      CALL dbcsr_release(t_2c_int_tmp(1))
      DO i_xyz = 1, 3
         CALL dbcsr_release(t_2c_der_tmp(1, i_xyz))
      END DO

      DEALLOCATE (row_bsize, col_bsize)
      ALLOCATE (row_bsize(SIZE(sizes_AO)))
      ALLOCATE (col_bsize(SIZE(sizes_AO)))
      row_bsize(:) = sizes_AO(:)
      col_bsize(:) = sizes_AO(:)

      CALL create_2c_tensor(t_2c_template, dist1, dist2, pgrid_t2c, force_data%bsizes_AO_split, &
                            force_data%bsizes_AO_split, name='(AO| AO)')
      DEALLOCATE (dist1, dist2)

      DO i_xyz = 1, 3
         CALL dbcsr_create(t_2c_der_tmp(1, i_xyz), "(P|Q) RPA der", dbcsr_dist, &
                           dbcsr_type_antisymmetric, row_bsize, col_bsize)
      END DO

      identity_pot%potential_type = do_potential_id
      CALL build_2c_neighbor_lists(nl_2c, basis_set_ao, basis_set_ao, identity_pot, &
                                   "RPA_2c_nl_metric", qs_env, sym_ij=.TRUE., dist_2d=dist_2d)
      CALL build_2c_derivatives(t_2c_der_tmp, mp2_env%ri_rpa_im_time%eps_filter, qs_env, nl_2c, &
                                basis_set_ao, basis_set_ao, identity_pot)
      CALL release_neighbor_list_sets(nl_2c)

      IF (use_virial) THEN
         CALL build_2c_neighbor_lists(force_data%nl_2c_ovlp, basis_set_ao, basis_set_ao, identity_pot, &
                                      "RPA_2c_nl_metric", qs_env, sym_ij=.FALSE., dist_2d=dist_2d)
      END IF

      CALL dbcsr_create(force_data%inv_ovlp, template=matrix_s(1)%matrix)
      CALL dbcsr_copy(force_data%inv_ovlp, matrix_s(1)%matrix)
      CALL cp_dbcsr_cholesky_decompose(force_data%inv_ovlp, para_env=para_env, blacs_env=blacs_env)
      CALL cp_dbcsr_cholesky_invert(force_data%inv_ovlp, para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)

      DO i_xyz = 1, 3
         CALL dbt_create(t_2c_der_tmp(1, i_xyz), t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(t_2c_der_tmp(1, i_xyz), t_2c_tmp)
         CALL dbt_create(t_2c_template, force_data%t_2c_der_ovlp(i_xyz))
         CALL dbt_copy(t_2c_tmp, force_data%t_2c_der_ovlp(i_xyz), move_data=.TRUE.)
         CALL dbt_filter(force_data%t_2c_der_ovlp(i_xyz), mp2_env%ri_rpa_im_time%eps_filter)
         CALL dbt_destroy(t_2c_tmp)
         CALL dbcsr_clear(t_2c_der_tmp(1, i_xyz))
      END DO

      !Create the rest of the 2-center AO tensors
      nspins = dft_control%nspins
      ALLOCATE (force_data%P_virt(nspins), force_data%P_occ(nspins))
      ALLOCATE (force_data%sum_YP_tau(nspins), force_data%sum_O_tau(nspins))
      DO ispin = 1, nspins
         ALLOCATE (force_data%P_virt(ispin)%matrix, force_data%P_occ(ispin)%matrix)
         ALLOCATE (force_data%sum_YP_tau(ispin)%matrix, force_data%sum_O_tau(ispin)%matrix)
         CALL dbcsr_create(force_data%P_virt(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(force_data%P_occ(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(force_data%sum_O_tau(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(force_data%sum_YP_tau(ispin)%matrix, template=matrix_s(1)%matrix)

         CALL dbcsr_copy(force_data%sum_O_tau(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_copy(force_data%sum_YP_tau(ispin)%matrix, matrix_s(1)%matrix)

         CALL dbcsr_set(force_data%sum_O_tau(ispin)%matrix, 0.0_dp)
         CALL dbcsr_set(force_data%sum_YP_tau(ispin)%matrix, 0.0_dp)
      END DO

      !Populate the density matrices: 1 = P_virt*S +P_occ*S ==> P_virt = S^-1 - P_occ
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      CALL dbcsr_copy(force_data%P_occ(1)%matrix, rho_ao(1)%matrix)
      IF (nspins == 1) THEN
         CALL dbcsr_scale(force_data%P_occ(1)%matrix, 0.5_dp) !because double occupency
      ELSE
         CALL dbcsr_copy(force_data%P_occ(2)%matrix, rho_ao(2)%matrix)
      END IF
      DO ispin = 1, nspins
         CALL dbcsr_copy(force_data%P_virt(ispin)%matrix, force_data%inv_ovlp)
         CALL dbcsr_add(force_data%P_virt(ispin)%matrix, force_data%P_occ(ispin)%matrix, 1.0_dp, -1.0_dp)
      END DO

      DO ibasis = 1, SIZE(basis_set_ao)
         orb_basis => basis_set_ao(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(orb_basis, eps_pgf_orb_old)
         ri_basis => basis_set_ri_aux(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(ri_basis, eps_pgf_orb_old)
      END DO

      CALL dbt_destroy(t_2c_template)
      CALL dbcsr_release(dbcsr_work)
      DO i_xyz = 1, 3
         CALL dbcsr_release(t_2c_der_tmp(1, i_xyz))
      END DO
      DEALLOCATE (row_bsize, col_bsize)
      CALL dbt_pgrid_destroy(pgrid_t3c)
      CALL dbt_pgrid_destroy(pgrid_t2c)
      CALL dbcsr_distribution_release(dbcsr_dist)
      CALL timestop(handle)

   END SUBROUTINE init_im_time_forces

! **************************************************************************************************
!> \brief Updates the cubic-scaling SOS-Laplace-MP2 contribution to the forces at each quadrature point
!> \param force_data ...
!> \param mat_P_omega ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param starts_array_mc_block ...
!> \param ends_array_mc_block ...
!> \param num_integ_points ...
!> \param nmo ...
!> \param Eigenval ...
!> \param tau_tj ...
!> \param tau_wj ...
!> \param cut_memory ...
!> \param Pspin ...
!> \param Qspin ...
!> \param open_shell ...
!> \param unit_nr ...
!> \param dbcsr_time ...
!> \param dbcsr_nflop ...
!> \param mp2_env ...
!> \param qs_env ...
!> \note In open-shell, we need to take Q from one spin, and everything from the other
! **************************************************************************************************
   SUBROUTINE calc_laplace_loop_forces(force_data, mat_P_omega, t_3c_M, t_3c_O, t_3c_O_compressed, &
                                       t_3c_O_ind, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                                       fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                       fm_mo_coeff_virt_scaled, starts_array_mc, ends_array_mc, &
                                       starts_array_mc_block, ends_array_mc_block, num_integ_points, &
                                       nmo, Eigenval, tau_tj, tau_wj, cut_memory, Pspin, Qspin, &
                                       open_shell, unit_nr, dbcsr_time, dbcsr_nflop, mp2_env, qs_env)

      TYPE(im_time_force_type), INTENT(INOUT)            :: force_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: mat_P_omega
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_M, t_3c_O
      TYPE(hfx_compression_type), DIMENSION(:)           :: t_3c_O_compressed
      TYPE(block_ind_type), DIMENSION(:), INTENT(INOUT)  :: t_3c_O_ind
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mo_coeff_occ, fm_mo_coeff_virt
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mo_coeff_occ_scaled, &
                                                            fm_mo_coeff_virt_scaled
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc, ends_array_mc, &
                                                            starts_array_mc_block, &
                                                            ends_array_mc_block
      INTEGER, INTENT(IN)                                :: num_integ_points, nmo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tau_wj
      INTEGER, INTENT(IN)                                :: cut_memory, Pspin, Qspin
      LOGICAL, INTENT(IN)                                :: open_shell
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(dp), INTENT(INOUT)                            :: dbcsr_time
      INTEGER(int_8), INTENT(INOUT)                      :: dbcsr_nflop
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_laplace_loop_forces'

      INTEGER :: dummy_int, handle, handle2, i_mem, i_xyz, ibasis, ispin, j_xyz, jquad, k_xyz, &
         n_mem_RI, n_rep, natom, nkind, nspins, unit_nr_dbcsr
      INTEGER(int_8)                                     :: flop, nze, nze_ddint, nze_der_AO, &
                                                            nze_der_RI, nze_KQK
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, batch_blk_end, &
                                                            batch_blk_start, batch_end_RI, &
                                                            batch_start_RI, kind_of, mc_ranges, &
                                                            mc_ranges_RI
      INTEGER, DIMENSION(:, :), POINTER                  :: dummy_ptr
      LOGICAL                                            :: memory_info, use_virial
      REAL(dp)                                           :: eps_filter, eps_pgf_orb, &
                                                            eps_pgf_orb_old, fac, occ, occ_ddint, &
                                                            occ_der_AO, occ_der_RI, occ_KQK, &
                                                            omega, pref, t1, t2, tau
      REAL(dp), DIMENSION(3, 3)                          :: work_virial, work_virial_ovlp
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_occ, mat_dm_virt
      TYPE(dbcsr_type)                                   :: dbcsr_work1, dbcsr_work2, dbcsr_work3, &
                                                            exp_occ, exp_virt, R_occ, R_virt, &
                                                            virial_ovlp, Y_1, Y_2
      TYPE(dbt_type) :: t_2c_AO, t_2c_RI, t_2c_RI_2, t_2c_tmp, t_3c_0, t_3c_1, t_3c_3, t_3c_4, &
         t_3c_5, t_3c_6, t_3c_7, t_3c_8, t_3c_help_1, t_3c_help_2, t_3c_ints, t_3c_sparse, &
         t_3c_work, t_dm_occ, t_dm_virt, t_KQKT, t_M_occ, t_M_virt, t_Q, t_R_occ, t_R_virt
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_P
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_ao, basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis
      TYPE(libint_potential_type)                        :: identity_pot
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: qs_section
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (matrix_s, dummy_ptr, atomic_kind_set, force, matrix_s, matrix_ks, mat_dm_occ, mat_dm_virt)
      NULLIFY (dft_control, virial, particle_set, cell, para_env, orb_basis, ri_basis, qs_section)
      NULLIFY (qs_kind_set)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, matrix_s=matrix_s, natom=natom, atomic_kind_set=atomic_kind_set, &
                      force=force, matrix_ks=matrix_ks, dft_control=dft_control, virial=virial, &
                      particle_set=particle_set, cell=cell, para_env=para_env, nkind=nkind, &
                      qs_kind_set=qs_kind_set)
      eps_filter = mp2_env%ri_rpa_im_time%eps_filter
      nspins = dft_control%nspins

      memory_info = mp2_env%ri_rpa_im_time%memory_info
      IF (memory_info) THEN
         unit_nr_dbcsr = unit_nr
      ELSE
         unit_nr_dbcsr = 0
      END IF

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      IF (use_virial) virial%pv_calculate = .TRUE.

      IF (use_virial) THEN
         qs_section => section_vals_get_subs_vals(qs_env%input, "DFT%QS")
         CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", n_rep_val=n_rep)
         IF (n_rep /= 0) THEN
            CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", r_val=eps_pgf_orb)
         ELSE
            CALL section_vals_val_get(qs_section, "EPS_DEFAULT", r_val=eps_pgf_orb)
            eps_pgf_orb = SQRT(eps_pgf_orb)
         END IF
         eps_pgf_orb_old = dft_control%qs_control%eps_pgf_orb

         ALLOCATE (basis_set_ri_aux(nkind), basis_set_ao(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL basis_set_list_setup(basis_set_ao, "ORB", qs_kind_set)

         DO ibasis = 1, SIZE(basis_set_ao)
            orb_basis => basis_set_ao(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(orb_basis, eps_pgf_orb)
            ri_basis => basis_set_ri_aux(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(ri_basis, eps_pgf_orb)
         END DO
      END IF

      !We follow the general logic of the compute_mat_P_omega routine
      ALLOCATE (t_P(nspins))
      CALL dbt_create(force_data%t_2c_K, t_2c_RI)
      CALL dbt_create(force_data%t_2c_K, t_2c_RI_2)
      CALL dbt_create(force_data%t_2c_der_ovlp(1), t_2c_AO)

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

      ! Always do the batching of the MO on mu and sigma, such that it is consistent between
      ! the occupied and the virtual quantities
      ALLOCATE (mc_ranges(cut_memory + 1))
      mc_ranges(:cut_memory) = starts_array_mc_block(:)
      mc_ranges(cut_memory + 1) = ends_array_mc_block(cut_memory) + 1

      ! Also need some batching on the RI, because it loses sparsity at some point
      n_mem_RI = cut_memory
      CALL create_tensor_batches(force_data%bsizes_RI_split, n_mem_RI, batch_start_RI, batch_end_RI, &
                                 batch_blk_start, batch_blk_end)
      ALLOCATE (mc_ranges_RI(n_mem_RI + 1))
      mc_ranges_RI(1:n_mem_RI) = batch_blk_start(1:n_mem_RI)
      mc_ranges_RI(n_mem_RI + 1) = batch_blk_end(n_mem_RI) + 1
      DEALLOCATE (batch_blk_start, batch_blk_end)

      !Pre-allocate all required tensors and matrices
      DO ispin = 1, nspins
         CALL dbt_create(t_2c_RI, t_P(ispin))
      END DO
      CALL dbt_create(t_2c_RI, t_Q)
      CALL dbt_create(t_2c_RI, t_KQKT)
      CALL dbt_create(t_2c_AO, t_dm_occ)
      CALL dbt_create(t_2c_AO, t_dm_virt)

      !note: t_3c_O and t_3c_M have different mappings (map_1d, map_2d)
      CALL dbt_create(t_3c_O, t_M_occ)
      CALL dbt_create(t_3c_O, t_M_virt)
      CALL dbt_create(t_3c_O, t_3c_0)

      CALL dbt_create(t_3c_O, t_3c_1)
      CALL dbt_create(t_3c_O, t_3c_3)
      CALL dbt_create(t_3c_O, t_3c_4)
      CALL dbt_create(t_3c_O, t_3c_5)
      CALL dbt_create(t_3c_M, t_3c_6)
      CALL dbt_create(t_3c_M, t_3c_7)
      CALL dbt_create(t_3c_M, t_3c_8)
      CALL dbt_create(t_3c_M, t_3c_sparse)
      CALL dbt_create(t_3c_O, t_3c_help_1)
      CALL dbt_create(t_3c_O, t_3c_help_2)
      CALL dbt_create(t_2c_AO, t_R_occ)
      CALL dbt_create(t_2c_AO, t_R_virt)
      CALL dbt_create(t_3c_M, t_3c_ints)
      CALL dbt_create(t_3c_M, t_3c_work)

      !Pre-define the sparsity of t_3c_4 as a function of the derivatives
      occ_der_AO = 0; nze_der_AO = 0
      occ_der_RI = 0; nze_der_RI = 0
      DO i_xyz = 1, 3
         DO i_mem = 1, cut_memory
            CALL decompress_tensor(force_data%t_3c_der_RI(i_xyz), force_data%t_3c_der_RI_ind(i_mem, i_xyz)%ind, &
                                   force_data%t_3c_der_RI_comp(i_mem, i_xyz), mp2_env%ri_rpa_im_time%eps_compress)
            CALL get_tensor_occupancy(force_data%t_3c_der_RI(i_xyz), nze, occ)
            occ_der_RI = occ_der_RI + occ
            nze_der_RI = nze_der_RI + nze
            CALL dbt_copy(force_data%t_3c_der_RI(i_xyz), t_3c_sparse, summation=.TRUE., move_data=.TRUE.)

            CALL decompress_tensor(force_data%t_3c_der_AO(i_xyz), force_data%t_3c_der_AO_ind(i_mem, i_xyz)%ind, &
                                   force_data%t_3c_der_AO_comp(i_mem, i_xyz), mp2_env%ri_rpa_im_time%eps_compress)
            CALL get_tensor_occupancy(force_data%t_3c_der_AO(i_xyz), nze, occ)
            occ_der_AO = occ_der_AO + occ
            nze_der_AO = nze_der_AO + nze
            CALL dbt_copy(force_data%t_3c_der_AO(i_xyz), t_3c_sparse, order=[1, 3, 2], summation=.TRUE.)
            CALL dbt_copy(force_data%t_3c_der_AO(i_xyz), t_3c_sparse, summation=.TRUE., move_data=.TRUE.)
         END DO
      END DO
      occ_der_RI = occ_der_RI/3.0_dp
      occ_der_AO = occ_der_AO/3.0_dp
      nze_der_RI = nze_der_RI/3
      nze_der_AO = nze_der_AO/3

      CALL dbcsr_create(R_occ, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(R_virt, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(dbcsr_work1, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(dbcsr_work2, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(dbcsr_work3, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(exp_occ, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(exp_virt, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      IF (use_virial) CALL dbcsr_create(virial_ovlp, template=dbcsr_work1)

      CALL dbt_batched_contract_init(t_3c_0, batch_range_2=mc_ranges, batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_1, batch_range_2=mc_ranges, batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_3, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)
      CALL dbt_batched_contract_init(t_M_occ, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)
      CALL dbt_batched_contract_init(t_M_virt, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)

      CALL dbt_batched_contract_init(t_3c_ints, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_work, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)

      CALL dbt_batched_contract_init(t_3c_4, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_5, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_6, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_7, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_8, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_sparse, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)

      work_virial = 0.0_dp
      work_virial_ovlp = 0.0_dp
      DO jquad = 1, num_integ_points
         tau = tau_tj(jquad)
         omega = tau_wj(jquad)
         fac = -2.0_dp*omega*mp2_env%scale_S
         IF (open_shell) fac = 0.5_dp*fac
         occ_ddint = 0; nze_ddint = 0

         CALL para_env%sync()
         t1 = m_walltime()

         !Deal with the force contributions where there is no explicit 3-center quantities, i.e. the
         !forces due to the metric and potential derivatives
         DO ispin = 1, nspins
            CALL dbt_create(mat_P_omega(jquad, ispin)%matrix, t_2c_tmp)
            CALL dbt_copy_matrix_to_tensor(mat_P_omega(jquad, ispin)%matrix, t_2c_tmp)
            CALL dbt_copy(t_2c_tmp, t_P(ispin), move_data=.TRUE.)
            CALL dbt_filter(t_P(ispin), eps_filter)
            CALL dbt_destroy(t_2c_tmp)
         END DO

         !Q = K^T*P*K, open-shell: Q is from one spin, everything else from the other
         CALL dbt_contract(1.0_dp, t_P(Qspin), force_data%t_2c_K, 0.0_dp, t_2c_RI, &
                           contract_1=[2], notcontract_1=[1], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=eps_filter, &
                           flop=flop, unit_nr=unit_nr_dbcsr)
         dbcsr_nflop = dbcsr_nflop + flop
         CALL dbt_contract(1.0_dp, force_data%t_2c_K, t_2c_RI, 0.0_dp, t_Q, &
                           contract_1=[1], notcontract_1=[2], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=eps_filter, &
                           flop=flop, unit_nr=unit_nr_dbcsr)
         dbcsr_nflop = dbcsr_nflop + flop
         CALL dbt_clear(t_2c_RI)

         CALL perform_2c_ops(force, t_KQKT, force_data, fac, t_Q, t_P(Pspin), t_2c_RI, t_2c_RI_2, &
                             use_virial, atom_of_kind, kind_of, eps_filter, dbcsr_nflop, unit_nr_dbcsr)
         CALL get_tensor_occupancy(t_KQKT, nze_KQK, occ_KQK)

         !Calculate the pseudo-density matrix in tensor form. There are a few useless arguments for SOS-MP2
         CALL compute_mat_dm_global(fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, tau_tj, num_integ_points, &
                                    nmo, fm_mo_coeff_occ(Pspin), fm_mo_coeff_virt(Pspin), &
                                    fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, mat_dm_occ, mat_dm_virt, &
                                    matrix_s, Pspin, Eigenval(:, Pspin), 0.0_dp, eps_filter, &
                                    mp2_env%ri_rpa_im_time%memory_info, unit_nr, &
                                    jquad, .FALSE., .FALSE., qs_env, dummy_int, dummy_ptr, para_env)

         CALL dbt_create(mat_dm_occ(jquad, 1)%matrix, t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(mat_dm_occ(jquad, 1)%matrix, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_dm_occ, move_data=.TRUE.)
         CALL dbt_filter(t_dm_occ, eps_filter)
         CALL dbt_destroy(t_2c_tmp)

         CALL dbt_create(mat_dm_virt(jquad, 1)%matrix, t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(mat_dm_virt(jquad, 1)%matrix, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_dm_virt, move_data=.TRUE.)
         CALL dbt_filter(t_dm_virt, eps_filter)
         CALL dbt_destroy(t_2c_tmp)

         !Deal with the 3-center quantities.
         CALL perform_3c_ops(force, t_R_occ, t_R_virt, force_data, fac, cut_memory, n_mem_RI, &
                             t_KQKT, t_dm_occ, t_dm_virt, t_3c_O, t_3c_M, t_M_occ, t_M_virt, t_3c_0, t_3c_1, &
                             t_3c_3, t_3c_4, t_3c_5, t_3c_6, t_3c_7, t_3c_8, t_3c_sparse, t_3c_help_1, t_3c_help_2, &
                             t_3c_ints, t_3c_work, starts_array_mc, ends_array_mc, batch_start_RI, &
                             batch_end_RI, t_3c_O_compressed, t_3c_O_ind, use_virial, &
                             atom_of_kind, kind_of, eps_filter, occ_ddint, nze_ddint, dbcsr_nflop, &
                             unit_nr_dbcsr, mp2_env)

         CALL timeset(routineN//"_dbcsr", handle2)
         !We go back to DBCSR matrices from now on
         !Note: R matrices are in fact symmetric, but use a normal type for convenience
         CALL dbt_create(matrix_s(1)%matrix, t_2c_tmp)
         CALL dbt_copy(t_R_occ, t_2c_tmp, move_data=.TRUE.)
         CALL dbt_copy_tensor_to_matrix(t_2c_tmp, R_occ)

         CALL dbt_copy(t_R_virt, t_2c_tmp, move_data=.TRUE.)
         CALL dbt_copy_tensor_to_matrix(t_2c_tmp, R_virt)

         !Iteratively calculate the Y1 and Y2 matrices
         CALL dbcsr_multiply('N', 'N', tau, force_data%P_occ(Pspin)%matrix, &
                             matrix_ks(Pspin)%matrix, 0.0_dp, dbcsr_work1)
         CALL build_Y_matrix(Y_1, dbcsr_work1, force_data%P_occ(Pspin)%matrix, R_virt, eps_filter)
         CALL matrix_exponential(exp_occ, dbcsr_work1, 1.0_dp, 1.0_dp, eps_filter)

         CALL dbcsr_multiply('N', 'N', -tau, force_data%P_virt(Pspin)%matrix, &
                             matrix_ks(Pspin)%matrix, 0.0_dp, dbcsr_work1)
         CALL build_Y_matrix(Y_2, dbcsr_work1, force_data%P_virt(Pspin)%matrix, R_occ, eps_filter)
         CALL matrix_exponential(exp_virt, dbcsr_work1, 1.0_dp, 1.0_dp, eps_filter)

         !The force contribution coming from [-S^-1*(e^-tau*P_virt*F)^T*R_occ*S^-1
         !                                    +tau*S^-1*Y_2^T*F*S^-1] * der_S
         CALL dbcsr_multiply('N', 'N', 1.0_dp, R_occ, force_data%inv_ovlp, 0.0_dp, dbcsr_work1)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, exp_virt, dbcsr_work1, 0.0_dp, dbcsr_work3)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, force_data%inv_ovlp, dbcsr_work3, 0.0_dp, dbcsr_work2)

         CALL dbcsr_multiply('N', 'T', tau, force_data%inv_ovlp, Y_2, 0.0_dp, dbcsr_work3)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, dbcsr_work3, matrix_ks(Pspin)%matrix, 0.0_dp, dbcsr_work1)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, dbcsr_work1, force_data%inv_ovlp, 0.0_dp, dbcsr_work3)

         CALL dbcsr_add(dbcsr_work2, dbcsr_work3, 1.0_dp, -1.0_dp)

         CALL dbt_copy_matrix_to_tensor(dbcsr_work2, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_2c_AO, move_data=.TRUE.)

         pref = -1.0_dp*fac
         CALL get_2c_der_force(force, t_2c_AO, force_data%t_2c_der_ovlp, atom_of_kind, &
                               kind_of, force_data%idx_to_at_AO, pref, do_ovlp=.TRUE.)

         IF (use_virial) CALL dbcsr_add(virial_ovlp, dbcsr_work2, 1.0_dp, pref)

         !The final contribution from Tr[(tau*Y_1*P_occ - tau*Y_2*P_virt) * der_F]
         CALL dbcsr_multiply('N', 'N', tau*fac, Y_1, force_data%P_occ(Pspin)%matrix, 1.0_dp, &
                             force_data%sum_YP_tau(Pspin)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply('N', 'N', -tau*fac, Y_2, force_data%P_virt(Pspin)%matrix, 1.0_dp, &
                             force_data%sum_YP_tau(Pspin)%matrix, retain_sparsity=.TRUE.)

         !Build-up the RHS of the response equation.
         pref = -omega*mp2_env%scale_S
         CALL dbcsr_multiply('N', 'N', pref, R_virt, exp_occ, 1.0_dp, &
                             force_data%sum_O_tau(Pspin)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply('N', 'N', -pref, R_occ, exp_virt, 1.0_dp, &
                             force_data%sum_O_tau(Pspin)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply('N', 'N', pref*tau, matrix_ks(Pspin)%matrix, Y_1, 1.0_dp, &
                             force_data%sum_O_tau(Pspin)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply('N', 'N', pref*tau, matrix_ks(Pspin)%matrix, Y_2, 1.0_dp, &
                             force_data%sum_O_tau(Pspin)%matrix, retain_sparsity=.TRUE.)

         CALL timestop(handle2)

         !Print some info
         CALL para_env%sync()
         t2 = m_walltime()
         dbcsr_time = dbcsr_time + t2 - t1

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/T3,A,1X,I3,A)') &
               'RPA_LOW_SCALING_INFO| Info for time point', jquad, '    (gradients)'
            WRITE (unit_nr, '(T6,A,T56,F25.6)') &
               'Execution time (s):', t2 - t1
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of 3c AO derivs:', REAL(nze_der_AO, dp), '/', occ_der_AO*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of 3c RI derivs:', REAL(nze_der_RI, dp), '/', occ_der_RI*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of the Docc * Dvirt * 3c-int tensor', REAL(nze_ddint, dp), '/', occ_ddint*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of KQK^T 2c-tensor:', REAL(nze_KQK, dp), '/', occ_KQK*100, '%'
            CALL m_flush(unit_nr)
         END IF

         !intermediate clean-up
         CALL dbcsr_release(Y_1)
         CALL dbcsr_release(Y_2)
         CALL dbt_destroy(t_2c_tmp)
      END DO !jquad

      CALL dbt_batched_contract_finalize(t_3c_0)
      CALL dbt_batched_contract_finalize(t_3c_1)
      CALL dbt_batched_contract_finalize(t_3c_3)
      CALL dbt_batched_contract_finalize(t_M_occ)
      CALL dbt_batched_contract_finalize(t_M_virt)

      CALL dbt_batched_contract_finalize(t_3c_ints)
      CALL dbt_batched_contract_finalize(t_3c_work)

      CALL dbt_batched_contract_finalize(t_3c_4)
      CALL dbt_batched_contract_finalize(t_3c_5)
      CALL dbt_batched_contract_finalize(t_3c_6)
      CALL dbt_batched_contract_finalize(t_3c_7)
      CALL dbt_batched_contract_finalize(t_3c_8)
      CALL dbt_batched_contract_finalize(t_3c_sparse)

      !Calculate the 2c and 3c contributions to the virial
      IF (use_virial) THEN
         CALL dbt_copy(force_data%t_3c_virial_split, force_data%t_3c_virial, move_data=.TRUE.)
         CALL calc_3c_virial(work_virial, force_data%t_3c_virial, 1.0_dp, qs_env, force_data%nl_3c, &
                             basis_set_ri_aux, basis_set_ao, basis_set_ao, mp2_env%ri_metric, &
                             der_eps=mp2_env%ri_rpa_im_time%eps_filter, op_pos=1)

         CALL calc_2c_virial(work_virial, force_data%RI_virial_met, 1.0_dp, qs_env, force_data%nl_2c_met, &
                             basis_set_ri_aux, basis_set_ri_aux, mp2_env%ri_metric)
         CALL dbcsr_clear(force_data%RI_virial_met)

         IF (.NOT. force_data%do_periodic) THEN
            CALL calc_2c_virial(work_virial, force_data%RI_virial_pot, 1.0_dp, qs_env, force_data%nl_2c_pot, &
                                basis_set_ri_aux, basis_set_ri_aux, mp2_env%potential_parameter)
            CALL dbcsr_clear(force_data%RI_virial_pot)
         END IF

         identity_pot%potential_type = do_potential_id
         CALL calc_2c_virial(work_virial_ovlp, virial_ovlp, 1.0_dp, qs_env, force_data%nl_2c_ovlp, &
                             basis_set_ao, basis_set_ao, identity_pot)
         CALL dbcsr_release(virial_ovlp)

         DO k_xyz = 1, 3
            DO j_xyz = 1, 3
               DO i_xyz = 1, 3
                  virial%pv_mp2(i_xyz, j_xyz) = virial%pv_mp2(i_xyz, j_xyz) &
                                                - work_virial(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
                  virial%pv_overlap(i_xyz, j_xyz) = virial%pv_overlap(i_xyz, j_xyz) &
                                                    - work_virial_ovlp(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
                  virial%pv_virial(i_xyz, j_xyz) = virial%pv_virial(i_xyz, j_xyz) &
                                                   - work_virial(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz) &
                                                   - work_virial_ovlp(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
               END DO
            END DO
         END DO
      END IF

      !Calculate the periodic contributions of (P|Q) to the force and the virial
      work_virial = 0.0_dp
      IF (force_data%do_periodic) THEN
         IF (mp2_env%eri_method == do_eri_gpw) THEN
            CALL get_2c_gpw_forces(force_data%G_PQ, force, work_virial, use_virial, mp2_env, qs_env)
         ELSE IF (mp2_env%eri_method == do_eri_mme) THEN
            CALL get_2c_mme_forces(force_data%G_PQ, force, mp2_env, qs_env)
            IF (use_virial) CPABORT("Stress tensor not available with MME intrgrals")
         ELSE
            CPABORT("Periodic case not possible with OS integrals")
         END IF
         CALL dbcsr_clear(force_data%G_PQ)
      END IF

      IF (use_virial) THEN
         virial%pv_mp2 = virial%pv_mp2 + work_virial
         virial%pv_virial = virial%pv_virial + work_virial
         virial%pv_calculate = .FALSE.

         DO ibasis = 1, SIZE(basis_set_ao)
            orb_basis => basis_set_ao(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(orb_basis, eps_pgf_orb_old)
            ri_basis => basis_set_ri_aux(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(ri_basis, eps_pgf_orb_old)
         END DO
      END IF

      !clean-up
      IF (ASSOCIATED(dummy_ptr)) DEALLOCATE (dummy_ptr)
      DO ispin = 1, nspins
         CALL dbt_destroy(t_P(ispin))
      END DO
      CALL dbt_destroy(t_3c_0)
      CALL dbt_destroy(t_3c_1)
      CALL dbt_destroy(t_3c_3)
      CALL dbt_destroy(t_3c_4)
      CALL dbt_destroy(t_3c_5)
      CALL dbt_destroy(t_3c_6)
      CALL dbt_destroy(t_3c_7)
      CALL dbt_destroy(t_3c_8)
      CALL dbt_destroy(t_3c_sparse)
      CALL dbt_destroy(t_3c_help_1)
      CALL dbt_destroy(t_3c_help_2)
      CALL dbt_destroy(t_3c_ints)
      CALL dbt_destroy(t_3c_work)
      CALL dbt_destroy(t_R_occ)
      CALL dbt_destroy(t_R_virt)
      CALL dbt_destroy(t_dm_occ)
      CALL dbt_destroy(t_dm_virt)
      CALL dbt_destroy(t_Q)
      CALL dbt_destroy(t_KQKT)
      CALL dbt_destroy(t_M_occ)
      CALL dbt_destroy(t_M_virt)
      CALL dbcsr_release(R_occ)
      CALL dbcsr_release(R_virt)
      CALL dbcsr_release(dbcsr_work1)
      CALL dbcsr_release(dbcsr_work2)
      CALL dbcsr_release(dbcsr_work3)
      CALL dbcsr_release(exp_occ)
      CALL dbcsr_release(exp_virt)

      CALL dbt_destroy(t_2c_RI)
      CALL dbt_destroy(t_2c_RI_2)
      CALL dbt_destroy(t_2c_AO)
      CALL dbcsr_deallocate_matrix_set(mat_dm_occ)
      CALL dbcsr_deallocate_matrix_set(mat_dm_virt)

      CALL timestop(handle)

   END SUBROUTINE calc_laplace_loop_forces

! **************************************************************************************************
!> \brief Updates the cubic-scaling RPA contribution to the forces at each quadrature point. This
!>        routine is adapted from the corresponding Laplace SOS-MP2 loop force one.
!> \param force_data ...
!> \param mat_P_omega ...
!> \param t_3c_M ...
!> \param t_3c_O ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param fm_scaled_dm_occ_tau ...
!> \param fm_scaled_dm_virt_tau ...
!> \param fm_mo_coeff_occ ...
!> \param fm_mo_coeff_virt ...
!> \param fm_mo_coeff_occ_scaled ...
!> \param fm_mo_coeff_virt_scaled ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param starts_array_mc_block ...
!> \param ends_array_mc_block ...
!> \param num_integ_points ...
!> \param nmo ...
!> \param Eigenval ...
!> \param e_fermi ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_cos_tf_w_to_t ...
!> \param tj ...
!> \param wj ...
!> \param tau_tj ...
!> \param cut_memory ...
!> \param ispin ...
!> \param open_shell ...
!> \param unit_nr ...
!> \param dbcsr_time ...
!> \param dbcsr_nflop ...
!> \param mp2_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calc_rpa_loop_forces(force_data, mat_P_omega, t_3c_M, t_3c_O, t_3c_O_compressed, &
                                   t_3c_O_ind, fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, &
                                   fm_mo_coeff_occ, fm_mo_coeff_virt, fm_mo_coeff_occ_scaled, &
                                   fm_mo_coeff_virt_scaled, starts_array_mc, ends_array_mc, &
                                   starts_array_mc_block, ends_array_mc_block, num_integ_points, &
                                   nmo, Eigenval, e_fermi, weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, &
                                   tj, wj, tau_tj, cut_memory, ispin, open_shell, unit_nr, dbcsr_time, &
                                   dbcsr_nflop, mp2_env, qs_env)

      TYPE(im_time_force_type), INTENT(INOUT)            :: force_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: mat_P_omega
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_M, t_3c_O
      TYPE(hfx_compression_type), DIMENSION(:)           :: t_3c_O_compressed
      TYPE(block_ind_type), DIMENSION(:), INTENT(INOUT)  :: t_3c_O_ind
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_scaled_dm_occ_tau, &
                                                            fm_scaled_dm_virt_tau
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mo_coeff_occ, fm_mo_coeff_virt
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mo_coeff_occ_scaled, &
                                                            fm_mo_coeff_virt_scaled
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc, ends_array_mc, &
                                                            starts_array_mc_block, &
                                                            ends_array_mc_block
      INTEGER, INTENT(IN)                                :: num_integ_points, nmo
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: Eigenval
      REAL(KIND=dp), INTENT(IN)                          :: e_fermi
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(IN)                                      :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(IN)                                      :: tj, wj
      REAL(KIND=dp), DIMENSION(0:num_integ_points), &
         INTENT(IN)                                      :: tau_tj
      INTEGER, INTENT(IN)                                :: cut_memory, ispin
      LOGICAL, INTENT(IN)                                :: open_shell
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(dp), INTENT(INOUT)                            :: dbcsr_time
      INTEGER(int_8), INTENT(INOUT)                      :: dbcsr_nflop
      TYPE(mp2_type)                                     :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_rpa_loop_forces'

      INTEGER :: dummy_int, handle, handle2, i_mem, i_xyz, ibasis, iquad, j_xyz, jquad, k_xyz, &
         n_mem_RI, n_rep, natom, nkind, nspins, unit_nr_dbcsr
      INTEGER(int_8)                                     :: flop, nze, nze_ddint, nze_der_AO, &
                                                            nze_der_RI, nze_KBK
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, batch_blk_end, &
                                                            batch_blk_start, batch_end_RI, &
                                                            batch_start_RI, kind_of, mc_ranges, &
                                                            mc_ranges_RI
      INTEGER, DIMENSION(:, :), POINTER                  :: dummy_ptr
      LOGICAL                                            :: memory_info, use_virial
      REAL(dp) :: eps_filter, eps_pgf_orb, eps_pgf_orb_old, fac, occ, occ_ddint, occ_der_AO, &
         occ_der_RI, occ_KBK, omega, pref, spin_fac, t1, t2, tau, weight
      REAL(dp), DIMENSION(3, 3)                          :: work_virial, work_virial_ovlp
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_P_tau, matrix_ks, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_occ, mat_dm_virt
      TYPE(dbcsr_type)                                   :: dbcsr_work1, dbcsr_work2, dbcsr_work3, &
                                                            dbcsr_work_symm, exp_occ, exp_virt, &
                                                            R_occ, R_virt, virial_ovlp, Y_1, Y_2
      TYPE(dbt_type) :: t_2c_AO, t_2c_RI, t_2c_RI_2, t_2c_tmp, t_3c_0, t_3c_1, t_3c_3, t_3c_4, &
         t_3c_5, t_3c_6, t_3c_7, t_3c_8, t_3c_help_1, t_3c_help_2, t_3c_ints, t_3c_sparse, &
         t_3c_work, t_dm_occ, t_dm_virt, t_KBKT, t_M_occ, t_M_virt, t_P, t_R_occ, t_R_virt
      TYPE(dbt_type), ALLOCATABLE, DIMENSION(:)          :: t_B
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_ao, basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis
      TYPE(libint_potential_type)                        :: identity_pot
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: qs_section
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (matrix_s, dummy_ptr, atomic_kind_set, force, matrix_s, matrix_ks, mat_dm_occ, mat_dm_virt)
      NULLIFY (dft_control, virial, particle_set, cell, blacs_env, para_env, orb_basis, ri_basis)
      NULLIFY (qs_kind_set)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, matrix_s=matrix_s, natom=natom, atomic_kind_set=atomic_kind_set, &
                      force=force, matrix_ks=matrix_ks, dft_control=dft_control, virial=virial, &
                      particle_set=particle_set, cell=cell, blacs_env=blacs_env, para_env=para_env, &
                      qs_kind_set=qs_kind_set, nkind=nkind)
      eps_filter = mp2_env%ri_rpa_im_time%eps_filter
      nspins = dft_control%nspins

      memory_info = mp2_env%ri_rpa_im_time%memory_info
      IF (memory_info) THEN
         unit_nr_dbcsr = unit_nr
      ELSE
         unit_nr_dbcsr = 0
      END IF

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      IF (use_virial) virial%pv_calculate = .TRUE.

      IF (use_virial) THEN
         qs_section => section_vals_get_subs_vals(qs_env%input, "DFT%QS")
         CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", n_rep_val=n_rep)
         IF (n_rep /= 0) THEN
            CALL section_vals_val_get(qs_section, "EPS_PGF_ORB", r_val=eps_pgf_orb)
         ELSE
            CALL section_vals_val_get(qs_section, "EPS_DEFAULT", r_val=eps_pgf_orb)
            eps_pgf_orb = SQRT(eps_pgf_orb)
         END IF
         eps_pgf_orb_old = dft_control%qs_control%eps_pgf_orb

         ALLOCATE (basis_set_ri_aux(nkind), basis_set_ao(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL basis_set_list_setup(basis_set_ao, "ORB", qs_kind_set)

         DO ibasis = 1, SIZE(basis_set_ao)
            orb_basis => basis_set_ao(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(orb_basis, eps_pgf_orb)
            ri_basis => basis_set_ri_aux(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(ri_basis, eps_pgf_orb)
         END DO
      END IF

      !We follow the general logic of the compute_mat_P_omega routine
      CALL dbt_create(force_data%t_2c_K, t_2c_RI)
      CALL dbt_create(force_data%t_2c_K, t_2c_RI_2)
      CALL dbt_create(force_data%t_2c_der_ovlp(1), t_2c_AO)

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

      ! Always do the batching of the MO on mu and sigma, such that it is consistent between
      ! the occupied and the virtual quantities
      ALLOCATE (mc_ranges(cut_memory + 1))
      mc_ranges(:cut_memory) = starts_array_mc_block(:)
      mc_ranges(cut_memory + 1) = ends_array_mc_block(cut_memory) + 1

      ! Also need some batching on the RI, because it loses sparsity at some point
      n_mem_RI = cut_memory
      CALL create_tensor_batches(force_data%bsizes_RI_split, n_mem_RI, batch_start_RI, batch_end_RI, &
                                 batch_blk_start, batch_blk_end)
      ALLOCATE (mc_ranges_RI(n_mem_RI + 1))
      mc_ranges_RI(1:n_mem_RI) = batch_blk_start(1:n_mem_RI)
      mc_ranges_RI(n_mem_RI + 1) = batch_blk_end(n_mem_RI) + 1
      DEALLOCATE (batch_blk_start, batch_blk_end)

      !Pre-allocate all required tensors and matrices
      CALL dbt_create(t_2c_RI, t_P)
      CALL dbt_create(t_2c_RI, t_KBKT)
      CALL dbt_create(t_2c_AO, t_dm_occ)
      CALL dbt_create(t_2c_AO, t_dm_virt)

      !note: t_3c_O and t_3c_M have different mappings (map_1d, map_2d)
      CALL dbt_create(t_3c_O, t_M_occ)
      CALL dbt_create(t_3c_O, t_M_virt)
      CALL dbt_create(t_3c_O, t_3c_0)

      CALL dbt_create(t_3c_O, t_3c_1)
      CALL dbt_create(t_3c_O, t_3c_3)
      CALL dbt_create(t_3c_O, t_3c_4)
      CALL dbt_create(t_3c_O, t_3c_5)
      CALL dbt_create(t_3c_M, t_3c_6)
      CALL dbt_create(t_3c_M, t_3c_7)
      CALL dbt_create(t_3c_M, t_3c_8)
      CALL dbt_create(t_3c_M, t_3c_sparse)
      CALL dbt_create(t_3c_O, t_3c_help_1)
      CALL dbt_create(t_3c_O, t_3c_help_2)
      CALL dbt_create(t_2c_AO, t_R_occ)
      CALL dbt_create(t_2c_AO, t_R_virt)
      CALL dbt_create(t_3c_M, t_3c_ints)
      CALL dbt_create(t_3c_M, t_3c_work)

      !Before entring the loop, need to compute the 2c tensors B = (1 + Q(w))^-1 - 1, for each
      !frequency grid point, before doing the transformation to the time grid
      ALLOCATE (t_B(num_integ_points))
      DO jquad = 1, num_integ_points
         CALL dbt_create(t_2c_RI, t_B(jquad))
      END DO

      ALLOCATE (mat_P_tau(num_integ_points))
      DO jquad = 1, num_integ_points
         ALLOCATE (mat_P_tau(jquad)%matrix)
         CALL dbcsr_create(mat_P_tau(jquad)%matrix, template=mat_P_omega(jquad, ispin)%matrix)
      END DO

      CALL dbcsr_create(dbcsr_work_symm, template=force_data%G_PQ, matrix_type=dbcsr_type_symmetric)
      CALL dbt_create(dbcsr_work_symm, t_2c_tmp)

      !loop over freqeuncies
      DO iquad = 1, num_integ_points
         omega = tj(iquad)

         !calculate (1 + Q(w))^-1 - 1 for the given freq.
         !Always take spin alpha (get 2*alpha in closed shell, and alpha+beta in open-shell)
         CALL dbcsr_copy(dbcsr_work_symm, mat_P_omega(iquad, 1)%matrix)
         CALL dbt_copy_matrix_to_tensor(dbcsr_work_symm, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_2c_RI, move_data=.TRUE.)

         CALL dbt_contract(1.0_dp, t_2c_RI, force_data%t_2c_K, 0.0_dp, t_2c_RI_2, &
                           contract_1=[2], notcontract_1=[1], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=eps_filter, &
                           flop=flop, unit_nr=unit_nr_dbcsr)
         dbcsr_nflop = dbcsr_nflop + flop
         CALL dbt_contract(1.0_dp, force_data%t_2c_K, t_2c_RI_2, 0.0_dp, t_2c_RI, &
                           contract_1=[1], notcontract_1=[2], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=eps_filter, &
                           flop=flop, unit_nr=unit_nr_dbcsr)
         CALL dbt_copy(t_2c_RI, t_2c_tmp, move_data=.TRUE.)
         CALL dbt_copy_tensor_to_matrix(t_2c_tmp, dbcsr_work_symm)
         CALL dbcsr_add_on_diag(dbcsr_work_symm, 1.0_dp)

         CALL cp_dbcsr_cholesky_decompose(dbcsr_work_symm, para_env=para_env, blacs_env=blacs_env)
         CALL cp_dbcsr_cholesky_invert(dbcsr_work_symm, para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)

         CALL dbcsr_add_on_diag(dbcsr_work_symm, -1.0_dp)

         DO jquad = 1, num_integ_points
            tau = tau_tj(jquad)

            !the P matrix to time.
            weight = weights_cos_tf_w_to_t(jquad, iquad)*COS(tau*omega)
            IF (open_shell) THEN
               IF (ispin == 1) THEN
                  !mat_P_omega contains the sum of alpha and beta spin => we only want alpha
                  CALL dbcsr_add(mat_P_tau(jquad)%matrix, mat_P_omega(iquad, 1)%matrix, 1.0_dp, weight)
                  CALL dbcsr_add(mat_P_tau(jquad)%matrix, mat_P_omega(iquad, 2)%matrix, 1.0_dp, -weight)
               ELSE
                  CALL dbcsr_add(mat_P_tau(jquad)%matrix, mat_P_omega(iquad, 2)%matrix, 1.0_dp, weight)
               END IF
            ELSE
               !factor 0.5 because originam matrix Q is scaled by 2 in RPA (spin)
               weight = 0.5_dp*weight
               CALL dbcsr_add(mat_P_tau(jquad)%matrix, mat_P_omega(iquad, 1)%matrix, 1.0_dp, weight)
            END IF

            !convert B matrix to time
            weight = weights_cos_tf_t_to_w(iquad, jquad)*COS(tau*omega)*wj(iquad)
            CALL dbt_copy_matrix_to_tensor(dbcsr_work_symm, t_2c_tmp)
            CALL dbt_scale(t_2c_tmp, weight)
            CALL dbt_copy(t_2c_tmp, t_B(jquad), summation=.TRUE., move_data=.TRUE.)
         END DO
      END DO
      CALL dbt_destroy(t_2c_tmp)
      CALL dbcsr_release(dbcsr_work_symm)
      CALL dbt_clear(t_2c_RI)
      CALL dbt_clear(t_2c_RI_2)

      !Pre-define the sparsity of t_3c_4 as a function of the derivatives
      occ_der_AO = 0; nze_der_AO = 0
      occ_der_RI = 0; nze_der_RI = 0
      DO i_xyz = 1, 3
         DO i_mem = 1, cut_memory
            CALL decompress_tensor(force_data%t_3c_der_RI(i_xyz), force_data%t_3c_der_RI_ind(i_mem, i_xyz)%ind, &
                                   force_data%t_3c_der_RI_comp(i_mem, i_xyz), mp2_env%ri_rpa_im_time%eps_compress)
            CALL get_tensor_occupancy(force_data%t_3c_der_RI(i_xyz), nze, occ)
            occ_der_RI = occ_der_RI + occ
            nze_der_RI = nze_der_RI + nze
            CALL dbt_copy(force_data%t_3c_der_RI(i_xyz), t_3c_sparse, summation=.TRUE., move_data=.TRUE.)

            CALL decompress_tensor(force_data%t_3c_der_AO(i_xyz), force_data%t_3c_der_AO_ind(i_mem, i_xyz)%ind, &
                                   force_data%t_3c_der_AO_comp(i_mem, i_xyz), mp2_env%ri_rpa_im_time%eps_compress)
            CALL get_tensor_occupancy(force_data%t_3c_der_AO(i_xyz), nze, occ)
            occ_der_AO = occ_der_AO + occ
            nze_der_AO = nze_der_AO + nze
            CALL dbt_copy(force_data%t_3c_der_AO(i_xyz), t_3c_sparse, order=[1, 3, 2], summation=.TRUE.)
            CALL dbt_copy(force_data%t_3c_der_AO(i_xyz), t_3c_sparse, summation=.TRUE., move_data=.TRUE.)
         END DO
      END DO
      occ_der_RI = occ_der_RI/3.0_dp
      occ_der_AO = occ_der_AO/3.0_dp
      nze_der_RI = nze_der_RI/3
      nze_der_AO = nze_der_AO/3

      CALL dbcsr_create(R_occ, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(R_virt, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(dbcsr_work_symm, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(dbcsr_work1, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(dbcsr_work2, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(dbcsr_work3, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(exp_occ, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(exp_virt, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      IF (use_virial) CALL dbcsr_create(virial_ovlp, template=dbcsr_work1)

      CALL dbt_batched_contract_init(t_3c_0, batch_range_2=mc_ranges, batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_1, batch_range_2=mc_ranges, batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_3, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)
      CALL dbt_batched_contract_init(t_M_occ, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)
      CALL dbt_batched_contract_init(t_M_virt, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)

      CALL dbt_batched_contract_init(t_3c_ints, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_work, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges)

      CALL dbt_batched_contract_init(t_3c_4, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_5, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_6, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_7, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_8, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)
      CALL dbt_batched_contract_init(t_3c_sparse, batch_range_1=mc_ranges_RI, batch_range_2=mc_ranges, &
                                     batch_range_3=mc_ranges)

      fac = 1.0_dp/fourpi*mp2_env%ri_rpa%scale_rpa
      IF (open_shell) fac = 0.5_dp*fac

      work_virial = 0.0_dp
      work_virial_ovlp = 0.0_dp
      DO jquad = 1, num_integ_points
         tau = tau_tj(jquad)
         occ_ddint = 0; nze_ddint = 0

         CALL para_env%sync()
         t1 = m_walltime()

         !Deal with the force contributions where there is no explicit 3-center quantities, i.e. the
         !forces due to the metric and potential derivatives
         CALL dbt_create(mat_P_tau(jquad)%matrix, t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(mat_P_tau(jquad)%matrix, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_P, move_data=.TRUE.)
         CALL dbt_filter(t_P, eps_filter)
         CALL dbt_destroy(t_2c_tmp)

         CALL perform_2c_ops(force, t_KBKT, force_data, fac, t_B(jquad), t_P, t_2c_RI, t_2c_RI_2, &
                             use_virial, atom_of_kind, kind_of, eps_filter, dbcsr_nflop, unit_nr_dbcsr)
         CALL get_tensor_occupancy(t_KBKT, nze_KBK, occ_KBK)

         !Calculate the pseudo-density matrix in tensor form. There are a few useless arguments for SOS-MP2
         CALL compute_mat_dm_global(fm_scaled_dm_occ_tau, fm_scaled_dm_virt_tau, tau_tj, num_integ_points, &
                                    nmo, fm_mo_coeff_occ(ispin), fm_mo_coeff_virt(ispin), &
                                    fm_mo_coeff_occ_scaled, fm_mo_coeff_virt_scaled, mat_dm_occ, mat_dm_virt, &
                                    matrix_s, ispin, Eigenval(:, ispin), e_fermi, eps_filter, &
                                    mp2_env%ri_rpa_im_time%memory_info, unit_nr, &
                                    jquad, .FALSE., .FALSE., qs_env, dummy_int, dummy_ptr, para_env)

         CALL dbt_create(mat_dm_occ(jquad, 1)%matrix, t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(mat_dm_occ(jquad, 1)%matrix, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_dm_occ, move_data=.TRUE.)
         CALL dbt_filter(t_dm_occ, eps_filter)
         CALL dbt_destroy(t_2c_tmp)

         CALL dbt_create(mat_dm_virt(jquad, 1)%matrix, t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(mat_dm_virt(jquad, 1)%matrix, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_dm_virt, move_data=.TRUE.)
         CALL dbt_filter(t_dm_virt, eps_filter)
         CALL dbt_destroy(t_2c_tmp)

         !Deal with the 3-center quantities.
         CALL perform_3c_ops(force, t_R_occ, t_R_virt, force_data, fac, cut_memory, n_mem_RI, &
                             t_KBKT, t_dm_occ, t_dm_virt, t_3c_O, t_3c_M, t_M_occ, t_M_virt, t_3c_0, t_3c_1, &
                             t_3c_3, t_3c_4, t_3c_5, t_3c_6, t_3c_7, t_3c_8, t_3c_sparse, t_3c_help_1, t_3c_help_2, &
                             t_3c_ints, t_3c_work, starts_array_mc, ends_array_mc, batch_start_RI, &
                             batch_end_RI, t_3c_O_compressed, t_3c_O_ind, use_virial, &
                             atom_of_kind, kind_of, eps_filter, occ_ddint, nze_ddint, dbcsr_nflop, &
                             unit_nr_dbcsr, mp2_env)

         CALL timeset(routineN//"_dbcsr", handle2)
         !We go back to DBCSR matrices from now on
         !Note: R matrices are in fact symmetric, but use a normal type for convenience
         CALL dbt_create(matrix_s(1)%matrix, t_2c_tmp)
         CALL dbt_copy(t_R_occ, t_2c_tmp, move_data=.TRUE.)
         CALL dbt_copy_tensor_to_matrix(t_2c_tmp, R_occ)

         CALL dbt_copy(t_R_virt, t_2c_tmp, move_data=.TRUE.)
         CALL dbt_copy_tensor_to_matrix(t_2c_tmp, R_virt)

         !Iteratively calculate the Y1 and Y2 matrices
         CALL dbcsr_copy(dbcsr_work_symm, matrix_ks(ispin)%matrix)
         CALL dbcsr_add(dbcsr_work_symm, matrix_s(1)%matrix, 1.0_dp, -e_fermi)
         CALL dbcsr_multiply('N', 'N', tau, force_data%P_occ(ispin)%matrix, &
                             dbcsr_work_symm, 0.0_dp, dbcsr_work1)
         CALL build_Y_matrix(Y_1, dbcsr_work1, force_data%P_occ(ispin)%matrix, R_virt, eps_filter)
         CALL matrix_exponential(exp_occ, dbcsr_work1, 1.0_dp, 1.0_dp, eps_filter)

         CALL dbcsr_multiply('N', 'N', -tau, force_data%P_virt(ispin)%matrix, &
                             dbcsr_work_symm, 0.0_dp, dbcsr_work1)
         CALL build_Y_matrix(Y_2, dbcsr_work1, force_data%P_virt(ispin)%matrix, R_occ, eps_filter)
         CALL matrix_exponential(exp_virt, dbcsr_work1, 1.0_dp, 1.0_dp, eps_filter)

         !The force contribution coming from [-S^-1*(e^-tau*P_virt*F)^T*R_occ*S^-1
         !                                    +tau*S^-1*Y_2^T*F*S^-1] * der_S
         !as well as -tau*e_fermi*Y_1*P^occ + tau*e_fermi*Y_2*P^virt
         CALL dbcsr_multiply('N', 'N', 1.0_dp, R_occ, force_data%inv_ovlp, 0.0_dp, dbcsr_work1)
         CALL dbcsr_multiply('T', 'N', 1.0_dp, exp_virt, dbcsr_work1, 0.0_dp, dbcsr_work3)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, force_data%inv_ovlp, dbcsr_work3, 0.0_dp, dbcsr_work2)

         CALL dbcsr_multiply('N', 'T', tau, force_data%inv_ovlp, Y_2, 0.0_dp, dbcsr_work3)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, dbcsr_work3, dbcsr_work_symm, 0.0_dp, dbcsr_work1)
         CALL dbcsr_multiply('N', 'N', -1.0_dp, dbcsr_work1, force_data%inv_ovlp, 1.0_dp, dbcsr_work2)

         CALL dbcsr_multiply('N', 'T', tau*e_fermi, force_data%P_occ(ispin)%matrix, Y_1, 1.0_dp, dbcsr_work2)
         CALL dbcsr_multiply('N', 'T', -tau*e_fermi, force_data%P_virt(ispin)%matrix, Y_2, 1.0_dp, dbcsr_work2)

         CALL dbt_copy_matrix_to_tensor(dbcsr_work2, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_2c_AO, move_data=.TRUE.)

         pref = -1.0_dp*fac
         CALL get_2c_der_force(force, t_2c_AO, force_data%t_2c_der_ovlp, atom_of_kind, &
                               kind_of, force_data%idx_to_at_AO, pref, do_ovlp=.TRUE.)

         IF (use_virial) CALL dbcsr_add(virial_ovlp, dbcsr_work2, 1.0_dp, pref)

         !The final contribution from Tr[(tau*Y_1*P_occ - tau*Y_2*P_virt) * der_F]
         CALL dbcsr_multiply('N', 'N', fac*tau, Y_1, force_data%P_occ(ispin)%matrix, 1.0_dp, &
                             force_data%sum_YP_tau(ispin)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply('N', 'N', -fac*tau, Y_2, force_data%P_virt(ispin)%matrix, 1.0_dp, &
                             force_data%sum_YP_tau(ispin)%matrix, retain_sparsity=.TRUE.)

         spin_fac = 0.5_dp*fac
         IF (open_shell) spin_fac = 2.0_dp*spin_fac
         !Build-up the RHS of the response equation.
         CALL dbcsr_multiply('N', 'N', 1.0_dp*spin_fac, R_virt, exp_occ, 1.0_dp, &
                             force_data%sum_O_tau(ispin)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply('N', 'N', -1.0_dp*spin_fac, R_occ, exp_virt, 1.0_dp, &
                             force_data%sum_O_tau(ispin)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply('N', 'N', tau*spin_fac, dbcsr_work_symm, Y_1, 1.0_dp, &
                             force_data%sum_O_tau(ispin)%matrix, retain_sparsity=.TRUE.)
         CALL dbcsr_multiply('N', 'N', tau*spin_fac, dbcsr_work_symm, Y_2, 1.0_dp, &
                             force_data%sum_O_tau(ispin)%matrix, retain_sparsity=.TRUE.)

         CALL timestop(handle2)

         !Print some info
         CALL para_env%sync()
         t2 = m_walltime()
         dbcsr_time = dbcsr_time + t2 - t1

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/T3,A,1X,I3,A)') &
               'RPA_LOW_SCALING_INFO| Info for time point', jquad, '    (gradients)'
            WRITE (unit_nr, '(T6,A,T56,F25.6)') &
               'Time:', t2 - t1
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of 3c AO derivs:', REAL(nze_der_AO, dp), '/', occ_der_AO*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of 3c RI derivs:', REAL(nze_der_RI, dp), '/', occ_der_RI*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of the Docc * Dvirt * 3c-int tensor', REAL(nze_ddint, dp), '/', occ_ddint*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of KBK^T 2c-tensor:', REAL(nze_KBK, dp), '/', occ_KBK*100, '%'
            CALL m_flush(unit_nr)
         END IF

         !intermediate clean-up
         CALL dbcsr_release(Y_1)
         CALL dbcsr_release(Y_2)
         CALL dbt_destroy(t_2c_tmp)

      END DO !jquad

      CALL dbt_batched_contract_finalize(t_3c_0)
      CALL dbt_batched_contract_finalize(t_3c_1)
      CALL dbt_batched_contract_finalize(t_3c_3)
      CALL dbt_batched_contract_finalize(t_M_occ)
      CALL dbt_batched_contract_finalize(t_M_virt)

      CALL dbt_batched_contract_finalize(t_3c_ints)
      CALL dbt_batched_contract_finalize(t_3c_work)

      CALL dbt_batched_contract_finalize(t_3c_4)
      CALL dbt_batched_contract_finalize(t_3c_5)
      CALL dbt_batched_contract_finalize(t_3c_6)
      CALL dbt_batched_contract_finalize(t_3c_7)
      CALL dbt_batched_contract_finalize(t_3c_8)
      CALL dbt_batched_contract_finalize(t_3c_sparse)

      !Calculate the 2c and 3c contributions to the virial
      IF (use_virial) THEN
         CALL dbt_copy(force_data%t_3c_virial_split, force_data%t_3c_virial, move_data=.TRUE.)
         CALL calc_3c_virial(work_virial, force_data%t_3c_virial, 1.0_dp, qs_env, force_data%nl_3c, &
                             basis_set_ri_aux, basis_set_ao, basis_set_ao, mp2_env%ri_metric, &
                             der_eps=mp2_env%ri_rpa_im_time%eps_filter, op_pos=1)

         CALL calc_2c_virial(work_virial, force_data%RI_virial_met, 1.0_dp, qs_env, force_data%nl_2c_met, &
                             basis_set_ri_aux, basis_set_ri_aux, mp2_env%ri_metric)
         CALL dbcsr_clear(force_data%RI_virial_met)

         IF (.NOT. force_data%do_periodic) THEN
            CALL calc_2c_virial(work_virial, force_data%RI_virial_pot, 1.0_dp, qs_env, force_data%nl_2c_pot, &
                                basis_set_ri_aux, basis_set_ri_aux, mp2_env%potential_parameter)
            CALL dbcsr_clear(force_data%RI_virial_pot)
         END IF

         identity_pot%potential_type = do_potential_id
         CALL calc_2c_virial(work_virial_ovlp, virial_ovlp, 1.0_dp, qs_env, force_data%nl_2c_ovlp, &
                             basis_set_ao, basis_set_ao, identity_pot)
         CALL dbcsr_release(virial_ovlp)

         DO k_xyz = 1, 3
            DO j_xyz = 1, 3
               DO i_xyz = 1, 3
                  virial%pv_mp2(i_xyz, j_xyz) = virial%pv_mp2(i_xyz, j_xyz) &
                                                - work_virial(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
                  virial%pv_overlap(i_xyz, j_xyz) = virial%pv_overlap(i_xyz, j_xyz) &
                                                    - work_virial_ovlp(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
                  virial%pv_virial(i_xyz, j_xyz) = virial%pv_virial(i_xyz, j_xyz) &
                                                   - work_virial(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz) &
                                                   - work_virial_ovlp(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
               END DO
            END DO
         END DO
      END IF

      !Calculate the periodic contributions of (P|Q) to the force and the virial
      work_virial = 0.0_dp
      IF (force_data%do_periodic) THEN
         IF (mp2_env%eri_method == do_eri_gpw) THEN
            CALL get_2c_gpw_forces(force_data%G_PQ, force, work_virial, use_virial, mp2_env, qs_env)
         ELSE IF (mp2_env%eri_method == do_eri_mme) THEN
            CALL get_2c_mme_forces(force_data%G_PQ, force, mp2_env, qs_env)
            IF (use_virial) CPABORT("Stress tensor not available with MME intrgrals")
         ELSE
            CPABORT("Periodic case not possible with OS integrals")
         END IF
         CALL dbcsr_clear(force_data%G_PQ)
      END IF

      IF (use_virial) THEN
         virial%pv_mp2 = virial%pv_mp2 + work_virial
         virial%pv_virial = virial%pv_virial + work_virial
         virial%pv_calculate = .FALSE.

         DO ibasis = 1, SIZE(basis_set_ao)
            orb_basis => basis_set_ao(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(orb_basis, eps_pgf_orb_old)
            ri_basis => basis_set_ri_aux(ibasis)%gto_basis_set
            CALL init_interaction_radii_orb_basis(ri_basis, eps_pgf_orb_old)
         END DO
      END IF

      !clean-up
      IF (ASSOCIATED(dummy_ptr)) DEALLOCATE (dummy_ptr)
      DO jquad = 1, num_integ_points
         CALL dbt_destroy(t_B(jquad))
      END DO
      CALL dbt_destroy(t_P)
      CALL dbt_destroy(t_3c_0)
      CALL dbt_destroy(t_3c_1)
      CALL dbt_destroy(t_3c_3)
      CALL dbt_destroy(t_3c_4)
      CALL dbt_destroy(t_3c_5)
      CALL dbt_destroy(t_3c_6)
      CALL dbt_destroy(t_3c_7)
      CALL dbt_destroy(t_3c_8)
      CALL dbt_destroy(t_3c_sparse)
      CALL dbt_destroy(t_3c_help_1)
      CALL dbt_destroy(t_3c_help_2)
      CALL dbt_destroy(t_3c_ints)
      CALL dbt_destroy(t_3c_work)
      CALL dbt_destroy(t_R_occ)
      CALL dbt_destroy(t_R_virt)
      CALL dbt_destroy(t_dm_occ)
      CALL dbt_destroy(t_dm_virt)
      CALL dbt_destroy(t_KBKT)
      CALL dbt_destroy(t_M_occ)
      CALL dbt_destroy(t_M_virt)
      CALL dbcsr_release(R_occ)
      CALL dbcsr_release(R_virt)
      CALL dbcsr_release(dbcsr_work_symm)
      CALL dbcsr_release(dbcsr_work1)
      CALL dbcsr_release(dbcsr_work2)
      CALL dbcsr_release(dbcsr_work3)
      CALL dbcsr_release(exp_occ)
      CALL dbcsr_release(exp_virt)

      CALL dbt_destroy(t_2c_RI)
      CALL dbt_destroy(t_2c_RI_2)
      CALL dbt_destroy(t_2c_AO)
      CALL dbcsr_deallocate_matrix_set(mat_dm_occ)
      CALL dbcsr_deallocate_matrix_set(mat_dm_virt)
      CALL dbcsr_deallocate_matrix_set(mat_P_tau)

      CALL timestop(handle)

   END SUBROUTINE calc_rpa_loop_forces

! **************************************************************************************************
!> \brief This subroutines performs the 2c tensor operations that are common accros low-scaling RPA
!>        and SOS-MP2, including forces and virial
!> \param force ...
!> \param t_KBKT returns the 2c tensor product of K*B*K^T
!> \param force_data ...
!> \param fac ...
!> \param t_B depending on RPA or SOS-MP2, t_B contains (1 + Q)^-1 - 1 or simply Q, respectively
!> \param t_P ...
!> \param t_2c_RI ...
!> \param t_2c_RI_2 ...
!> \param use_virial ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param eps_filter ...
!> \param dbcsr_nflop ...
!> \param unit_nr_dbcsr ...
! **************************************************************************************************
   SUBROUTINE perform_2c_ops(force, t_KBKT, force_data, fac, t_B, t_P, t_2c_RI, t_2c_RI_2, use_virial, &
                             atom_of_kind, kind_of, eps_filter, dbcsr_nflop, unit_nr_dbcsr)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbt_type), INTENT(INOUT)                      :: t_KBKT
      TYPE(im_time_force_type), INTENT(INOUT)            :: force_data
      REAL(dp), INTENT(IN)                               :: fac
      TYPE(dbt_type), INTENT(INOUT)                      :: t_B, t_P, t_2c_RI, t_2c_RI_2
      LOGICAL, INTENT(IN)                                :: use_virial
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of
      REAL(dp), INTENT(IN)                               :: eps_filter
      INTEGER(int_8), INTENT(INOUT)                      :: dbcsr_nflop
      INTEGER, INTENT(IN)                                :: unit_nr_dbcsr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'perform_2c_ops'

      INTEGER                                            :: handle
      INTEGER(int_8)                                     :: flop
      REAL(dp)                                           :: pref
      TYPE(dbt_type)                                     :: t_2c_tmp, t_2c_virial

      CALL timeset(routineN, handle)

      IF (use_virial) CALL dbt_create(force_data%RI_virial_pot, t_2c_virial)

      !P^T*K*B + P*K*B^T (note we calculate and save K*B*K^T for later, and P=P^T)
      CALL dbt_contract(1.0_dp, force_data%t_2c_K, t_B, 0.0_dp, t_2c_RI, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[1], notcontract_2=[2], &
                        map_1=[1], map_2=[2], filter_eps=eps_filter, &
                        flop=flop, unit_nr=unit_nr_dbcsr)
      dbcsr_nflop = dbcsr_nflop + flop

      CALL dbt_contract(1.0_dp, t_2c_RI, force_data%t_2c_K, 0.0_dp, t_KBKT, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[2], notcontract_2=[1], &
                        map_1=[1], map_2=[2], filter_eps=eps_filter, &
                        flop=flop, unit_nr=unit_nr_dbcsr)
      dbcsr_nflop = dbcsr_nflop + flop

      CALL dbt_contract(2.0_dp, t_P, t_2c_RI, 0.0_dp, t_2c_RI_2, & !t_2c_RI_2 holds P^T*K*B
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[1], notcontract_2=[2], &
                        map_1=[1], map_2=[2], filter_eps=eps_filter, &
                        flop=flop, unit_nr=unit_nr_dbcsr)
      dbcsr_nflop = dbcsr_nflop + flop
      CALL dbt_clear(t_2c_RI)
      !t_2c_RI_2 currently holds 2*P^T*K*B = P^T*K*B + P*K*B^T (because of symmetry)

      !For the metric contribution, we need S^-1*(P^T*K*B + P*K*B^T)*K^T
      CALL dbt_contract(1.0_dp, force_data%t_2c_inv_metric, t_2c_RI_2, 0.0_dp, t_2c_RI, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[1], notcontract_2=[2], &
                        map_1=[1], map_2=[2], filter_eps=eps_filter, &
                        flop=flop, unit_nr=unit_nr_dbcsr)
      dbcsr_nflop = dbcsr_nflop + flop

      CALL dbt_contract(1.0_dp, t_2c_RI, force_data%t_2c_K, 0.0_dp, t_2c_RI_2, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[2], notcontract_2=[1], &
                        map_1=[1], map_2=[2], filter_eps=eps_filter, &
                        flop=flop, unit_nr=unit_nr_dbcsr)
      dbcsr_nflop = dbcsr_nflop + flop

      !Here we do the trace for the force
      pref = -1.0_dp*fac
      CALL get_2c_der_force(force, t_2c_RI_2, force_data%t_2c_der_metric, atom_of_kind, &
                            kind_of, force_data%idx_to_at_RI, pref, do_mp2=.TRUE.)
      IF (use_virial) THEN
         CALL dbt_copy(t_2c_RI_2, t_2c_virial)
         CALL dbt_scale(t_2c_virial, pref)
         CALL dbt_copy_tensor_to_matrix(t_2c_virial, force_data%RI_virial_met, summation=.TRUE.)
         CALL dbt_clear(t_2c_virial)
      END IF

      !For the potential contribution, we need S^-1*(P^T*K*B + P*K*B^T)*V^-0.5
      !some of it is still in t_2c_RI: ( S^-1*(P^T*K*B + P*K*B^T) )
      CALL dbt_contract(1.0_dp, t_2c_RI, force_data%t_2c_pot_msqrt, 0.0_dp, t_2c_RI_2, &
                        contract_1=[2], notcontract_1=[1], &
                        contract_2=[1], notcontract_2=[2], &
                        map_1=[1], map_2=[2], filter_eps=eps_filter, &
                        flop=flop, unit_nr=unit_nr_dbcsr)
      dbcsr_nflop = dbcsr_nflop + flop

      !Here we do the trace for the force. In the periodic case, we store the matrix in G_PQ for later
      pref = 0.5_dp*fac
      IF (force_data%do_periodic) THEN
         CALL dbt_scale(t_2c_RI_2, pref)
         CALL dbt_create(force_data%G_PQ, t_2c_tmp)
         CALL dbt_copy(t_2c_RI_2, t_2c_tmp, move_data=.TRUE.)
         CALL dbt_copy_tensor_to_matrix(t_2c_tmp, force_data%G_PQ, summation=.TRUE.)
         CALL dbt_destroy(t_2c_tmp)
      ELSE
         CALL get_2c_der_force(force, t_2c_RI_2, force_data%t_2c_der_pot, atom_of_kind, &
                               kind_of, force_data%idx_to_at_RI, pref, do_mp2=.TRUE.)

         IF (use_virial) THEN
            CALL dbt_copy(t_2c_RI_2, t_2c_virial)
            CALL dbt_scale(t_2c_virial, pref)
            CALL dbt_copy_tensor_to_matrix(t_2c_virial, force_data%RI_virial_pot, summation=.TRUE.)
            CALL dbt_clear(t_2c_virial)
         END IF
      END IF

      CALL dbt_clear(t_2c_RI)
      CALL dbt_clear(t_2c_RI_2)

      IF (use_virial) CALL dbt_destroy(t_2c_virial)

      CALL timestop(handle)

   END SUBROUTINE perform_2c_ops

! **************************************************************************************************
!> \brief This subroutines performs the 3c tensor operations that are common accros low-scaling RPA
!>        and SOS-MP2, including forces and virial
!> \param force ...
!> \param t_R_occ ...
!> \param t_R_virt ...
!> \param force_data ...
!> \param fac ...
!> \param cut_memory ...
!> \param n_mem_RI ...
!> \param t_KBKT ...
!> \param t_dm_occ ...
!> \param t_dm_virt ...
!> \param t_3c_O ...
!> \param t_3c_M ...
!> \param t_M_occ ...
!> \param t_M_virt ...
!> \param t_3c_0 ...
!> \param t_3c_1 ...
!> \param t_3c_3 ...
!> \param t_3c_4 ...
!> \param t_3c_5 ...
!> \param t_3c_6 ...
!> \param t_3c_7 ...
!> \param t_3c_8 ...
!> \param t_3c_sparse ...
!> \param t_3c_help_1 ...
!> \param t_3c_help_2 ...
!> \param t_3c_ints ...
!> \param t_3c_work ...
!> \param starts_array_mc ...
!> \param ends_array_mc ...
!> \param batch_start_RI ...
!> \param batch_end_RI ...
!> \param t_3c_O_compressed ...
!> \param t_3c_O_ind ...
!> \param use_virial ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param eps_filter ...
!> \param occ_ddint ...
!> \param nze_ddint ...
!> \param dbcsr_nflop ...
!> \param unit_nr_dbcsr ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE perform_3c_ops(force, t_R_occ, t_R_virt, force_data, fac, cut_memory, n_mem_RI, &
                             t_KBKT, t_dm_occ, t_dm_virt, t_3c_O, t_3c_M, t_M_occ, t_M_virt, t_3c_0, t_3c_1, &
                             t_3c_3, t_3c_4, t_3c_5, t_3c_6, t_3c_7, t_3c_8, t_3c_sparse, t_3c_help_1, t_3c_help_2, &
                             t_3c_ints, t_3c_work, starts_array_mc, ends_array_mc, batch_start_RI, &
                             batch_end_RI, t_3c_O_compressed, t_3c_O_ind, use_virial, &
                             atom_of_kind, kind_of, eps_filter, occ_ddint, nze_ddint, dbcsr_nflop, &
                             unit_nr_dbcsr, mp2_env)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbt_type), INTENT(INOUT)                      :: t_R_occ, t_R_virt
      TYPE(im_time_force_type), INTENT(INOUT)            :: force_data
      REAL(dp), INTENT(IN)                               :: fac
      INTEGER, INTENT(IN)                                :: cut_memory, n_mem_RI
      TYPE(dbt_type), INTENT(INOUT) :: t_KBKT, t_dm_occ, t_dm_virt, t_3c_O, t_3c_M, t_M_occ, &
         t_M_virt, t_3c_0, t_3c_1, t_3c_3, t_3c_4, t_3c_5, t_3c_6, t_3c_7, t_3c_8, t_3c_sparse, &
         t_3c_help_1, t_3c_help_2, t_3c_ints, t_3c_work
      INTEGER, DIMENSION(:), INTENT(IN)                  :: starts_array_mc, ends_array_mc, &
                                                            batch_start_RI, batch_end_RI
      TYPE(hfx_compression_type), DIMENSION(:)           :: t_3c_O_compressed
      TYPE(block_ind_type), DIMENSION(:), INTENT(INOUT)  :: t_3c_O_ind
      LOGICAL, INTENT(IN)                                :: use_virial
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of
      REAL(dp), INTENT(IN)                               :: eps_filter
      REAL(dp), INTENT(INOUT)                            :: occ_ddint
      INTEGER(int_8), INTENT(INOUT)                      :: nze_ddint, dbcsr_nflop
      INTEGER, INTENT(IN)                                :: unit_nr_dbcsr
      TYPE(mp2_type)                                     :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'perform_3c_ops'

      INTEGER                                            :: dummy_int, handle, handle2, i_mem, &
                                                            i_xyz, j_mem, k_mem
      INTEGER(int_8)                                     :: flop, nze
      INTEGER, DIMENSION(2, 1)                           :: ibounds, jbounds, kbounds
      INTEGER, DIMENSION(2, 2)                           :: bounds_2c
      INTEGER, DIMENSION(2, 3)                           :: bounds_cpy
      INTEGER, DIMENSION(3)                              :: bounds_3c
      REAL(dp)                                           :: memory, occ, pref
      TYPE(block_ind_type), ALLOCATABLE, DIMENSION(:, :) :: blk_indices
      TYPE(hfx_compression_type), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: store_3c

      CALL timeset(routineN, handle)

      CALL dbt_get_info(t_3c_M, nfull_total=bounds_3c)

      !Pre-compute and compress KBK^T * (pq|R)
      ALLOCATE (store_3c(n_mem_RI, cut_memory))
      ALLOCATE (blk_indices(n_mem_RI, cut_memory))
      memory = 0.0_dp
      CALL timeset(routineN//"_pre_3c", handle2)
      !temporarily build the full int 3c tensor
      CALL dbt_copy(t_3c_O, t_3c_0)
      DO i_mem = 1, cut_memory
         CALL decompress_tensor(t_3c_O, t_3c_O_ind(i_mem)%ind, t_3c_O_compressed(i_mem), &
                                mp2_env%ri_rpa_im_time%eps_compress)
         CALL dbt_copy(t_3c_O, t_3c_ints)
         CALL dbt_copy(t_3c_O, t_3c_0, move_data=.TRUE., summation=.TRUE.)

         DO k_mem = 1, n_mem_RI
            kbounds(:, 1) = [batch_start_RI(k_mem), batch_end_RI(k_mem)]

            CALL alloc_containers(store_3c(k_mem, i_mem), 1)

            !contract witht KBK^T over the RI index and store
            CALL dbt_batched_contract_init(t_KBKT)
            CALL dbt_contract(1.0_dp, t_KBKT, t_3c_ints, 0.0_dp, t_3c_work, &
                              contract_1=[2], notcontract_1=[1], &
                              contract_2=[1], notcontract_2=[2, 3], &
                              map_1=[1], map_2=[2, 3], filter_eps=eps_filter, &
                              bounds_2=kbounds, flop=flop, unit_nr=unit_nr_dbcsr)
            CALL dbt_batched_contract_finalize(t_KBKT)
            dbcsr_nflop = dbcsr_nflop + flop

            CALL dbt_copy(t_3c_work, t_3c_M, move_data=.TRUE.)
            CALL compress_tensor(t_3c_M, blk_indices(k_mem, i_mem)%ind, store_3c(k_mem, i_mem), &
                                 mp2_env%ri_rpa_im_time%eps_compress, memory)
         END DO
      END DO !i_mem
      CALL dbt_clear(t_3c_M)
      CALL dbt_copy(t_3c_M, t_3c_ints)
      CALL timestop(handle2)

      CALL dbt_batched_contract_init(t_R_occ)
      CALL dbt_batched_contract_init(t_R_virt)
      DO i_mem = 1, cut_memory
         ibounds(:, 1) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]

         !Compute the matrices M (integrals in t_3c_0)
         CALL timeset(routineN//"_3c_M", handle2)
         CALL dbt_batched_contract_init(t_dm_occ)
         CALL dbt_contract(1.0_dp, t_3c_0, t_dm_occ, 0.0_dp, t_3c_1, &
                           contract_1=[3], notcontract_1=[1, 2], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1, 2], map_2=[3], filter_eps=eps_filter, &
                           bounds_3=ibounds, flop=flop, unit_nr=unit_nr_dbcsr)
         dbcsr_nflop = dbcsr_nflop + flop
         CALL dbt_batched_contract_finalize(t_dm_occ)
         CALL dbt_copy(t_3c_1, t_M_occ, order=[1, 3, 2], move_data=.TRUE.)

         CALL dbt_batched_contract_init(t_dm_virt)
         CALL dbt_contract(1.0_dp, t_3c_0, t_dm_virt, 0.0_dp, t_3c_1, &
                           contract_1=[3], notcontract_1=[1, 2], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1, 2], map_2=[3], filter_eps=eps_filter, &
                           bounds_3=ibounds, flop=flop, unit_nr=unit_nr_dbcsr)
         dbcsr_nflop = dbcsr_nflop + flop
         CALL dbt_batched_contract_finalize(t_dm_virt)
         CALL dbt_copy(t_3c_1, t_M_virt, order=[1, 3, 2], move_data=.TRUE.)
         CALL timestop(handle2)

         !Compute the R matrices
         CALL timeset(routineN//"_3c_R", handle2)
         DO k_mem = 1, n_mem_RI
            CALL decompress_tensor(t_3c_M, blk_indices(k_mem, i_mem)%ind, store_3c(k_mem, i_mem), &
                                   mp2_env%ri_rpa_im_time%eps_compress)
            CALL dbt_copy(t_3c_M, t_3c_3, move_data=.TRUE.)

            CALL dbt_contract(1.0_dp, t_M_occ, t_3c_3, 1.0_dp, t_R_occ, &
                              contract_1=[1, 2], notcontract_1=[3], &
                              contract_2=[1, 2], notcontract_2=[3], &
                              map_1=[1], map_2=[2], filter_eps=eps_filter, &
                              flop=flop, unit_nr=unit_nr_dbcsr)
            dbcsr_nflop = dbcsr_nflop + flop

            CALL dbt_contract(1.0_dp, t_M_virt, t_3c_3, 1.0_dp, t_R_virt, &
                              contract_1=[1, 2], notcontract_1=[3], &
                              contract_2=[1, 2], notcontract_2=[3], &
                              map_1=[1], map_2=[2], filter_eps=eps_filter, &
                              flop=flop, unit_nr=unit_nr_dbcsr)
            dbcsr_nflop = dbcsr_nflop + flop
         END DO
         CALL dbt_copy(t_3c_M, t_3c_3)
         CALL dbt_copy(t_3c_M, t_M_virt)
         CALL timestop(handle2)

         CALL dbt_copy(t_M_occ, t_3c_4, move_data=.TRUE.)

         DO j_mem = 1, cut_memory
            jbounds(:, 1) = [starts_array_mc(j_mem), ends_array_mc(j_mem)]

            bounds_cpy(:, 1) = [1, bounds_3c(1)]
            bounds_cpy(:, 2) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]
            bounds_cpy(:, 3) = [starts_array_mc(j_mem), ends_array_mc(j_mem)]
            CALL dbt_copy(t_3c_sparse, t_3c_7, bounds=bounds_cpy)

            CALL dbt_batched_contract_init(t_dm_virt)
            DO k_mem = 1, n_mem_RI
               bounds_2c(:, 1) = [batch_start_RI(k_mem), batch_end_RI(k_mem)]
               bounds_2c(:, 2) = [starts_array_mc(i_mem), ends_array_mc(i_mem)]

               CALL timeset(routineN//"_3c_dm", handle2)

               !Calculate (mu nu| P) * D_occ * D_virt
               !Note: technically need M_occ*D_virt + M_virt*D_occ, but it is equivalent to 2*M_occ*D_virt
               CALL dbt_contract(2.0_dp, t_3c_4, t_dm_virt, 0.0_dp, t_3c_5, &
                                 contract_1=[3], notcontract_1=[1, 2], &
                                 contract_2=[1], notcontract_2=[2], &
                                 map_1=[1, 2], map_2=[3], filter_eps=eps_filter, &
                                 bounds_2=bounds_2c, bounds_3=jbounds, flop=flop, unit_nr=unit_nr_dbcsr)
               dbcsr_nflop = dbcsr_nflop + flop

               CALL get_tensor_occupancy(t_3c_5, nze, occ)
               nze_ddint = nze_ddint + nze
               occ_ddint = occ_ddint + occ

               CALL dbt_copy(t_3c_5, t_3c_6, move_data=.TRUE.)
               CALL timestop(handle2)

               !Calculate the contraction of the above with K*B*K^T
               CALL timeset(routineN//"_3c_KBK", handle2)
               CALL dbt_batched_contract_init(t_KBKT)
               CALL dbt_contract(1.0_dp, t_KBKT, t_3c_6, 0.0_dp, t_3c_7, &
                                 contract_1=[2], notcontract_1=[1], &
                                 contract_2=[1], notcontract_2=[2, 3], &
                                 map_1=[1], map_2=[2, 3], &
                                 retain_sparsity=.TRUE., flop=flop, unit_nr=unit_nr_dbcsr)
               dbcsr_nflop = dbcsr_nflop + flop
               CALL dbt_batched_contract_finalize(t_KBKT)
               CALL timestop(handle2)
               CALL dbt_copy(t_3c_7, t_3c_8, summation=.TRUE.)

            END DO !k_mem
            CALL dbt_batched_contract_finalize(t_dm_virt)
         END DO !j_mem

         CALL dbt_copy(t_3c_8, t_3c_help_1, move_data=.TRUE.)

         pref = 1.0_dp*fac
         DO k_mem = 1, cut_memory
            DO i_xyz = 1, 3
               CALL dbt_clear(force_data%t_3c_der_RI(i_xyz))
               CALL decompress_tensor(force_data%t_3c_der_RI(i_xyz), force_data%t_3c_der_RI_ind(k_mem, i_xyz)%ind, &
                                      force_data%t_3c_der_RI_comp(k_mem, i_xyz), mp2_env%ri_rpa_im_time%eps_compress)
            END DO
            CALL get_force_from_3c_trace(force, t_3c_help_1, force_data%t_3c_der_RI, atom_of_kind, kind_of, &
                                         force_data%idx_to_at_RI, pref, do_mp2=.TRUE., deriv_dim=1)
         END DO

         IF (use_virial) THEN
            CALL dbt_copy(t_3c_help_1, t_3c_help_2)
            CALL dbt_scale(t_3c_help_2, pref)
            CALL dbt_copy(t_3c_help_2, force_data%t_3c_virial_split, summation=.TRUE., move_data=.TRUE.)
         END IF

         CALL dbt_copy(t_3c_help_1, t_3c_help_2)
         CALL dbt_copy(t_3c_help_1, t_3c_help_2, order=[1, 3, 2], move_data=.TRUE., summation=.TRUE.)
         DO k_mem = 1, cut_memory
            DO i_xyz = 1, 3
               CALL dbt_clear(force_data%t_3c_der_AO(i_xyz))
               CALL decompress_tensor(force_data%t_3c_der_AO(i_xyz), force_data%t_3c_der_AO_ind(k_mem, i_xyz)%ind, &
                                      force_data%t_3c_der_AO_comp(k_mem, i_xyz), mp2_env%ri_rpa_im_time%eps_compress)
            END DO
            CALL get_force_from_3c_trace(force, t_3c_help_2, force_data%t_3c_der_AO, atom_of_kind, kind_of, &
                                         force_data%idx_to_at_AO, pref, do_mp2=.TRUE., deriv_dim=3)
         END DO

         CALL dbt_clear(t_3c_help_2)
      END DO !i_mem
      CALL dbt_batched_contract_finalize(t_R_occ)
      CALL dbt_batched_contract_finalize(t_R_virt)

      DO k_mem = 1, n_mem_RI
         DO i_mem = 1, cut_memory
            CALL dealloc_containers(store_3c(k_mem, i_mem), dummy_int)
         END DO
      END DO
      DEALLOCATE (store_3c, blk_indices)

      CALL timestop(handle)

   END SUBROUTINE perform_3c_ops

! **************************************************************************************************
!> \brief All the forces that can be calculated after the loop on the Laplace quaradture, using
!>        terms collected during the said loop. This inludes the z-vector equation and its reponse
!>        forces, as well as the force coming from the trace with the derivative of the KS matrix
!> \param force_data ...
!> \param unit_nr ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calc_post_loop_forces(force_data, unit_nr, qs_env)

      TYPE(im_time_force_type), INTENT(INOUT)            :: force_data
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calc_post_loop_forces'

      INTEGER                                            :: handle, ispin, nao, nao_aux, nocc, nspins
      LOGICAL                                            :: do_exx
      REAL(dp)                                           :: focc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: cpmos, mo_occ
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_p_work, matrix_p_mp2, &
                                                            matrix_p_mp2_admm, matrix_s, &
                                                            matrix_s_aux, work_admm, YP_admm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(section_vals_type), POINTER                   :: hfx_section, lr_section

      NULLIFY (linres_control, p_env, dft_control, matrix_s, mos, mo_coeff, fm_struct, lr_section, &
               dbcsr_p_work, YP_admm, matrix_p_mp2, admm_env, work_admm, matrix_s_aux, matrix_p_mp2_admm)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, matrix_s=matrix_s, mos=mos)
      nspins = dft_control%nspins

      ! Setting up for the z-vector equation

      ! Initialize linres_control
      lr_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%WF_CORRELATION%LOW_SCALING%CPHF")

      ALLOCATE (linres_control)
      CALL section_vals_val_get(lr_section, "MAX_ITER", i_val=linres_control%max_iter)
      CALL section_vals_val_get(lr_section, "EPS_CONV", r_val=linres_control%eps)
      CALL section_vals_val_get(lr_section, "PRECONDITIONER", i_val=linres_control%preconditioner_type)
      CALL section_vals_val_get(lr_section, "ENERGY_GAP", r_val=linres_control%energy_gap)

      linres_control%do_kernel = .TRUE.
      linres_control%lr_triplet = .FALSE.
      linres_control%converged = .FALSE.
      linres_control%eps_filter = qs_env%mp2_env%ri_rpa_im_time%eps_filter

      CALL set_qs_env(qs_env, linres_control=linres_control)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(T3,A)') 'MP2_CPHF| Iterative solution of Z-Vector equations'
         WRITE (unit_nr, '(T3,A,T45,ES8.1)') 'MP2_CPHF| Convergence threshold:', linres_control%eps
         WRITE (unit_nr, '(T3,A,T45,I8)') 'MP2_CPHF| Maximum number of iterations: ', linres_control%max_iter
      END IF

      ALLOCATE (p_env)
      CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE., linres_control=linres_control)
      CALL p_env_psi0_changed(p_env, qs_env)

      ! Matrix allocation
      CALL dbcsr_allocate_matrix_set(p_env%p1, nspins)
      CALL dbcsr_allocate_matrix_set(p_env%w1, nspins)
      CALL dbcsr_allocate_matrix_set(dbcsr_p_work, nspins)
      DO ispin = 1, nspins
         ALLOCATE (p_env%p1(ispin)%matrix, p_env%w1(ispin)%matrix, dbcsr_p_work(ispin)%matrix)
         CALL dbcsr_create(matrix=p_env%p1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(matrix=p_env%w1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(matrix=dbcsr_p_work(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(p_env%p1(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_copy(p_env%w1(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_copy(dbcsr_p_work(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(p_env%p1(ispin)%matrix, 0.0_dp)
         CALL dbcsr_set(p_env%w1(ispin)%matrix, 0.0_dp)
         CALL dbcsr_set(dbcsr_p_work(ispin)%matrix, 0.0_dp)
      END DO

      IF (dft_control%do_admm) THEN
         CALL get_admm_env(qs_env%admm_env, matrix_s_aux_fit=matrix_s_aux)
         CALL dbcsr_allocate_matrix_set(p_env%p1_admm, nspins)
         CALL dbcsr_allocate_matrix_set(work_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (p_env%p1_admm(ispin)%matrix, work_admm(ispin)%matrix)
            CALL dbcsr_create(p_env%p1_admm(ispin)%matrix, template=matrix_s_aux(1)%matrix)
            CALL dbcsr_copy(p_env%p1_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
            CALL dbcsr_set(p_env%p1_admm(ispin)%matrix, 0.0_dp)
            CALL dbcsr_create(work_admm(ispin)%matrix, template=matrix_s_aux(1)%matrix)
            CALL dbcsr_copy(work_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
            CALL dbcsr_set(work_admm(ispin)%matrix, 0.0_dp)
         END DO
      END IF

      ! Preparing the RHS of the z-vector equation
      CALL prepare_for_response(force_data, qs_env)
      ALLOCATE (cpmos(nspins), mo_occ(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, nao=nao, homo=nocc)
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, ncol_global=nocc, &
                                  template_fmstruct=mo_coeff%matrix_struct)
         CALL cp_fm_create(cpmos(ispin), fm_struct)
         CALL cp_fm_set_all(cpmos(ispin), 0.0_dp)
         CALL cp_fm_create(mo_occ(ispin), fm_struct)
         CALL cp_fm_to_fm(mo_coeff, mo_occ(ispin), nocc)
         CALL cp_fm_struct_release(fm_struct)
      END DO

      ! in case of EXX, need to add the HF Hamiltonian to the RHS of the Z-vector equation
      ! Strategy: we take the ks_matrix, remove the current xc contribution, and then add the RPA HF one
      do_exx = .FALSE.
      IF (qs_env%mp2_env%method == ri_rpa_method_gpw) THEN
         hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
         CALL section_vals_get(hfx_section, explicit=do_exx)
      END IF

      IF (do_exx) THEN
         CALL add_exx_to_rhs(rhs=force_data%sum_O_tau, &
                             qs_env=qs_env, &
                             ext_hfx_section=hfx_section, &
                             x_data=qs_env%mp2_env%ri_rpa%x_data, &
                             recalc_integrals=.FALSE., &
                             do_admm=qs_env%mp2_env%ri_rpa%do_admm, &
                             do_exx=do_exx, &
                             reuse_hfx=qs_env%mp2_env%ri_rpa%reuse_hfx)
      END IF

      focc = 2.0_dp
      IF (nspins == 1) focc = 4.0_dp
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff, homo=nocc)
         CALL cp_dbcsr_sm_fm_multiply(force_data%sum_O_tau(ispin)%matrix, mo_occ(ispin), &
                                      cpmos(ispin), nocc, &
                                      alpha=focc, beta=0.0_dp)
      END DO

      ! The z-vector equation and associated forces
      CALL response_equation_new(qs_env, p_env, cpmos, unit_nr)

      ! Save the mp2 density matrix
      CALL get_qs_env(qs_env, matrix_p_mp2=matrix_p_mp2)
      IF (ASSOCIATED(matrix_p_mp2)) CALL dbcsr_deallocate_matrix_set(matrix_p_mp2)
      DO ispin = 1, nspins
         CALL dbcsr_copy(dbcsr_p_work(ispin)%matrix, p_env%p1(ispin)%matrix)
         CALL dbcsr_add(dbcsr_p_work(ispin)%matrix, force_data%sum_YP_tau(ispin)%matrix, 1.0_dp, 1.0_dp)
      END DO
      CALL set_ks_env(qs_env%ks_env, matrix_p_mp2=dbcsr_p_work)

      IF (dft_control%do_admm) THEN
         CALL dbcsr_allocate_matrix_set(YP_admm, nspins)
         CALL get_qs_env(qs_env, matrix_p_mp2_admm=matrix_p_mp2_admm, admm_env=admm_env)
         nao = admm_env%nao_orb
         nao_aux = admm_env%nao_aux_fit
         IF (ASSOCIATED(matrix_p_mp2_admm)) CALL dbcsr_deallocate_matrix_set(matrix_p_mp2_admm)
         DO ispin = 1, nspins

            !sum_YP_tau in the auxiliary basis
            CALL copy_dbcsr_to_fm(force_data%sum_YP_tau(ispin)%matrix, admm_env%work_orb_orb)
            CALL parallel_gemm('N', 'N', nao_aux, nao, nao, 1.0_dp, admm_env%A, admm_env%work_orb_orb, &
                               0.0_dp, admm_env%work_aux_orb)
            CALL parallel_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, admm_env%work_aux_orb, admm_env%A, &
                               0.0_dp, admm_env%work_aux_aux)
            CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, work_admm(ispin)%matrix, keep_sparsity=.TRUE.)

            !save the admm representation od sum_YP_tau
            ALLOCATE (YP_admm(ispin)%matrix)
            CALL dbcsr_create(YP_admm(ispin)%matrix, template=work_admm(ispin)%matrix)
            CALL dbcsr_copy(YP_admm(ispin)%matrix, work_admm(ispin)%matrix)

            CALL dbcsr_add(work_admm(ispin)%matrix, p_env%p1_admm(ispin)%matrix, 1.0_dp, 1.0_dp)

         END DO
         CALL set_ks_env(qs_env%ks_env, matrix_p_mp2_admm=work_admm)
      END IF

      !Calculate the response force and the force from the trace with F
      CALL update_im_time_forces(p_env, force_data%sum_O_tau, force_data%sum_YP_tau, YP_admm, qs_env)

      !clean-up
      IF (dft_control%do_admm) CALL dbcsr_deallocate_matrix_set(YP_admm)

      CALL cp_fm_release(cpmos)
      CALL cp_fm_release(mo_occ)
      CALL p_env_release(p_env)
      DEALLOCATE (p_env)

      CALL timestop(handle)

   END SUBROUTINE calc_post_loop_forces

! **************************************************************************************************
!> \brief Prepares the RHS of the z-vector equation. Apply the xc and HFX kernel on the previously
!>        stored sum_YP_tau density, and add it to the final force_data%sum_O_tau quantity
!> \param force_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE prepare_for_response(force_data, qs_env)

      TYPE(im_time_force_type), INTENT(INOUT)            :: force_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_for_response'

      INTEGER                                            :: handle, ispin, nao, nao_aux, nspins
      LOGICAL                                            :: do_hfx, do_tau, do_tau_admm
      REAL(dp)                                           :: ehartree
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dbcsr_p_work, ker_tau_admm, matrix_s, &
                                                            matrix_s_aux, work_admm
      TYPE(dbcsr_type)                                   :: dbcsr_work
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: rhoz_tot_gspace, zv_hartree_gspace, &
                                                            zv_hartree_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: rhoz_g, rhoz_r, tauz_r, v_xc, v_xc_tau
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit

      NULLIFY (pw_env, rhoz_g, rhoz_r, tauz_r, v_xc, v_xc_tau, &
               poisson_env, auxbas_pw_pool, dft_control, admm_env, xc_section, rho, rho_aux_fit, &
               task_list_aux_fit, ker_tau_admm, work_admm, dbcsr_p_work, matrix_s, hfx_section)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env, rho=rho, matrix_s=matrix_s)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, poisson_env=poisson_env)
      nspins = dft_control%nspins

      CALL dbcsr_allocate_matrix_set(dbcsr_p_work, nspins)
      DO ispin = 1, nspins
         ALLOCATE (dbcsr_p_work(ispin)%matrix)
         CALL dbcsr_create(matrix=dbcsr_p_work(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(dbcsr_p_work(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(dbcsr_p_work(ispin)%matrix, 0.0_dp)
      END DO

      !Apply the kernel on the density saved in force_data%sum_YP_tau
      ALLOCATE (rhoz_r(nspins), rhoz_g(nspins))
      DO ispin = 1, nspins
         CALL auxbas_pw_pool%create_pw(rhoz_r(ispin), &
                                       use_data=REALDATA3D, in_space=REALSPACE)
         CALL auxbas_pw_pool%create_pw(rhoz_g(ispin), &
                                       use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      END DO
      CALL auxbas_pw_pool%create_pw(rhoz_tot_gspace, &
                                    use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL auxbas_pw_pool%create_pw(zv_hartree_rspace, &
                                    use_data=REALDATA3D, in_space=REALSPACE)
      CALL auxbas_pw_pool%create_pw(zv_hartree_gspace, &
                                    use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      CALL pw_zero(rhoz_tot_gspace)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=qs_env%ks_env, matrix_p=force_data%sum_YP_tau(ispin)%matrix, &
                                 rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin))
         CALL pw_axpy(rhoz_g(ispin), rhoz_tot_gspace)
      END DO

      CALL pw_poisson_solve(poisson_env, rhoz_tot_gspace, ehartree, &
                            zv_hartree_gspace)

      CALL pw_transfer(zv_hartree_gspace, zv_hartree_rspace)
      CALL pw_scale(zv_hartree_rspace, zv_hartree_rspace%pw_grid%dvol)

      CALL qs_rho_get(rho, tau_r_valid=do_tau)
      IF (do_tau) THEN
         BLOCK
            TYPE(pw_type) :: tauz_g
            ALLOCATE (tauz_r(nspins))
            CALL auxbas_pw_pool%create_pw(tauz_g, &
                                          use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%create_pw(tauz_r(ispin), &
                                             use_data=REALDATA3D, in_space=REALSPACE)

               CALL calculate_rho_elec(ks_env=qs_env%ks_env, matrix_p=force_data%sum_YP_tau(ispin)%matrix, &
                                       rho=tauz_r(ispin), rho_gspace=tauz_g, compute_tau=.TRUE.)
            END DO
            CALL auxbas_pw_pool%give_back_pw(tauz_g)
         END BLOCK
      END IF

      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF

      !Primary XC kernel
      CALL create_kernel(qs_env, v_xc, v_xc_tau, rho, rhoz_r, rhoz_g, tauz_r, xc_section)

      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
         CALL pw_axpy(zv_hartree_rspace, v_xc(ispin))
         CALL integrate_v_rspace(qs_env=qs_env, &
                                 v_rspace=v_xc(ispin), &
                                 hmat=dbcsr_p_work(ispin), &
                                 calculate_forces=.FALSE.)

         CALL auxbas_pw_pool%give_back_pw(v_xc(ispin))
      END DO
      CALL auxbas_pw_pool%give_back_pw(rhoz_tot_gspace)
      CALL auxbas_pw_pool%give_back_pw(zv_hartree_rspace)
      CALL auxbas_pw_pool%give_back_pw(zv_hartree_gspace)
      DEALLOCATE (v_xc)

      IF (do_tau) THEN
         DO ispin = 1, nspins
            CALL pw_scale(v_xc_tau(ispin), v_xc_tau(ispin)%pw_grid%dvol)
            CALL integrate_v_rspace(qs_env=qs_env, &
                                    v_rspace=v_xc_tau(ispin), &
                                    hmat=dbcsr_p_work(ispin), &
                                    compute_tau=.TRUE., &
                                    calculate_forces=.FALSE.)
            CALL auxbas_pw_pool%give_back_pw(v_xc_tau(ispin))
         END DO
         DEALLOCATE (v_xc_tau)
      END IF

      !Auxiliary xc kernel (admm)
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CALL get_admm_env(admm_env, matrix_s_aux_fit=matrix_s_aux, &
                           task_list_aux_fit=task_list_aux_fit, rho_aux_fit=rho_aux_fit)

         CALL qs_rho_get(rho_aux_fit, tau_r_valid=do_tau_admm)

         CALL dbcsr_allocate_matrix_set(work_admm, nspins)
         CALL dbcsr_allocate_matrix_set(ker_tau_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (work_admm(ispin)%matrix, ker_tau_admm(ispin)%matrix)
            CALL dbcsr_create(work_admm(ispin)%matrix, template=matrix_s_aux(1)%matrix)
            CALL dbcsr_copy(work_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
            CALL dbcsr_set(work_admm(ispin)%matrix, 0.0_dp)
            CALL dbcsr_create(ker_tau_admm(ispin)%matrix, template=matrix_s_aux(1)%matrix)
            CALL dbcsr_copy(ker_tau_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
            CALL dbcsr_set(ker_tau_admm(ispin)%matrix, 0.0_dp)
         END DO

         !get the density in the auxuliary density
         CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
         CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
         CPASSERT(ASSOCIATED(admm_env%work_aux_aux))
         nao = admm_env%nao_orb
         nao_aux = admm_env%nao_aux_fit
         DO ispin = 1, nspins
            CALL copy_dbcsr_to_fm(force_data%sum_YP_tau(ispin)%matrix, admm_env%work_orb_orb)
            CALL parallel_gemm('N', 'N', nao_aux, nao, nao, 1.0_dp, admm_env%A, admm_env%work_orb_orb, &
                               0.0_dp, admm_env%work_aux_orb)
            CALL parallel_gemm('N', 'T', nao_aux, nao_aux, nao, 1.0_dp, admm_env%work_aux_orb, admm_env%A, &
                               0.0_dp, admm_env%work_aux_aux)
            CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, ker_tau_admm(ispin)%matrix, keep_sparsity=.TRUE.)
         END DO

         IF (.NOT. qs_env%admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
            DO ispin = 1, nspins
               CALL pw_zero(rhoz_r(ispin))
               CALL pw_zero(rhoz_g(ispin))
               CALL calculate_rho_elec(ks_env=qs_env%ks_env, matrix_p=ker_tau_admm(ispin)%matrix, &
                                       rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin), &
                                       basis_type="AUX_FIT", task_list_external=task_list_aux_fit)
            END DO

            IF (do_tau_admm) THEN
               BLOCK
                  TYPE(pw_type) :: tauz_g
                  CALL auxbas_pw_pool%create_pw(tauz_g, &
                                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                  DO ispin = 1, nspins
                     CALL pw_zero(tauz_r(ispin))
                     CALL calculate_rho_elec(ks_env=qs_env%ks_env, matrix_p=ker_tau_admm(ispin)%matrix, &
                                             rho=tauz_r(ispin), rho_gspace=tauz_g, &
                                             basis_type="AUX_FIT", task_list_external=task_list_aux_fit, &
                                             compute_tau=.TRUE.)
                  END DO
                  CALL auxbas_pw_pool%give_back_pw(tauz_g)
               END BLOCK
            END IF

            xc_section => admm_env%xc_section_aux
            CALL create_kernel(qs_env, v_xc, v_xc_tau, rho_aux_fit, rhoz_r, rhoz_g, tauz_r, xc_section)

            DO ispin = 1, nspins
               CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
               CALL integrate_v_rspace(qs_env=qs_env, &
                                       v_rspace=v_xc(ispin), &
                                       hmat=work_admm(ispin), &
                                       calculate_forces=.FALSE., &
                                       basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
               CALL auxbas_pw_pool%give_back_pw(v_xc(ispin))
            END DO
            DEALLOCATE (v_xc)

            IF (do_tau_admm) THEN
               DO ispin = 1, nspins
                  CALL pw_scale(v_xc_tau(ispin), v_xc_tau(ispin)%pw_grid%dvol)
                  CALL integrate_v_rspace(qs_env=qs_env, &
                                          v_rspace=v_xc_tau(ispin), &
                                          hmat=work_admm(ispin), &
                                          calculate_forces=.FALSE., &
                                          basis_type="AUX_FIT", &
                                          task_list_external=task_list_aux_fit, &
                                          compute_tau=.TRUE.)
                  CALL auxbas_pw_pool%give_back_pw(v_xc_tau(ispin))
               END DO
               DEALLOCATE (v_xc_tau)
            END IF
         END IF !admm
      END IF

      DO ispin = 1, nspins
         CALL auxbas_pw_pool%give_back_pw(rhoz_r(ispin))
         CALL auxbas_pw_pool%give_back_pw(rhoz_g(ispin))
      END DO
      DEALLOCATE (rhoz_r, rhoz_g)

      IF (do_tau) THEN
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(tauz_r(ispin))
         END DO
         DEALLOCATE (tauz_r)
      END IF

      !HFX kernel
      hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         IF (dft_control%do_admm) THEN
            CALL tddft_hfx_matrix(work_admm, ker_tau_admm, qs_env, .FALSE., .FALSE.)

            !Going back to primary basis
            CALL dbcsr_create(dbcsr_work, template=dbcsr_p_work(1)%matrix)
            CALL dbcsr_copy(dbcsr_work, dbcsr_p_work(1)%matrix)
            CALL dbcsr_set(dbcsr_work, 0.0_dp)
            DO ispin = 1, nspins
               CALL copy_dbcsr_to_fm(work_admm(ispin)%matrix, admm_env%work_aux_aux)
               CALL parallel_gemm('N', 'N', nao_aux, nao, nao_aux, 1.0_dp, admm_env%work_aux_aux, admm_env%A, &
                                  0.0_dp, admm_env%work_aux_orb)
               CALL parallel_gemm('T', 'N', nao, nao, nao_aux, 1.0_dp, admm_env%A, admm_env%work_aux_orb, &
                                  0.0_dp, admm_env%work_orb_orb)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbcsr_work, keep_sparsity=.TRUE.)
               CALL dbcsr_add(dbcsr_p_work(ispin)%matrix, dbcsr_work, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbcsr_work)
            CALL dbcsr_deallocate_matrix_set(ker_tau_admm)
         ELSE
            CALL tddft_hfx_matrix(dbcsr_p_work, force_data%sum_YP_tau, qs_env, .FALSE., .FALSE.)
         END IF
      END IF

      DO ispin = 1, nspins
         CALL dbcsr_add(force_data%sum_O_tau(ispin)%matrix, dbcsr_p_work(ispin)%matrix, 1.0_dp, 1.0_dp)
      END DO

      CALL dbcsr_deallocate_matrix_set(dbcsr_p_work)
      CALL dbcsr_deallocate_matrix_set(work_admm)

      CALL timestop(handle)

   END SUBROUTINE prepare_for_response

! **************************************************************************************************
!> \brief Calculate the force and virial due to the (P|Q) GPW integral derivatives
!> \param G_PQ ...
!> \param force ...
!> \param h_stress ...
!> \param use_virial ...
!> \param mp2_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_2c_gpw_forces(G_PQ, force, h_stress, use_virial, mp2_env, qs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: G_PQ
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT)           :: h_stress
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_2c_gpw_forces'

      INTEGER :: atom_a, color, handle, i, i_RI, i_xyz, iatom, igrid_level, ikind, ipgf, iset, j, &
         j_RI, jatom, lb_RI, n_RI, natom, ncoa, ncoms, nkind, nproc, nseta, o1, offset, pdims(2), &
         sgfa, ub_RI
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, iproc_map, kind_of, &
                                                            sizes_RI
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, la_max, la_min, npgfa, nsgfa, &
                                                            row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, pgrid
      LOGICAL                                            :: found, one_proc_group
      REAL(dp)                                           :: cutoff_old, radius, relative_cutoff_old
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: e_cutoff_old, wf_vector
      REAL(dp), DIMENSION(3)                             :: force_a, force_b, ra
      REAL(dp), DIMENSION(3, 3)                          :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_tmp, I_ab, pab, pblock, sphi_a, zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_type)                                   :: tmp_G_PQ
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_ext
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env_ext
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: dvg(3), pot_g, psi_L, rho_g, rho_g_copy, &
                                                            rho_r
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), DIMENSION(:), POINTER   :: rs_v
      TYPE(task_list_type), POINTER                      :: task_list_ext

      NULLIFY (sab_orb, task_list_ext, particle_set, qs_kind_set, dft_control, pw_env_ext, auxbas_pw_pool, &
               poisson_env, atomic_kind_set, para_env, cell, rs_v, mos, basis_set_a)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env, sab_orb=sab_orb, &
                      natom=natom, nkind=nkind, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      mos=mos, cell=cell, atomic_kind_set=atomic_kind_set)

      !The idea is to use GPW to compute the integrals and derivatives. Because the potential needs
      !to be calculated for each phi_j (column) of all AO pairs, and because that is expensive, we want
      !to minimize the amount of time we do that. Therefore, we work with a special distribution, where
      !each column of the resulting DBCSR matrix is mapped to a sub-communicator.

      !Try to get the optimal pdims (we want a grid that is flat: many cols, few rows)
      IF (para_env%num_pe <= natom) THEN
         pdims(1) = 1
         pdims(2) = para_env%num_pe
      ELSE
         DO i = natom, 1, -1
            IF (MODULO(para_env%num_pe, i) == 0) THEN
               pdims(1) = para_env%num_pe/i
               pdims(2) = i
               EXIT
            END IF
         END DO
      END IF

      ALLOCATE (row_dist(natom), col_dist(natom))
      DO iatom = 1, natom
         row_dist(iatom) = MODULO(iatom, pdims(1))
      END DO
      DO jatom = 1, natom
         col_dist(jatom) = MODULO(jatom, pdims(2))
      END DO

      ALLOCATE (pgrid(0:pdims(1) - 1, 0:pdims(2) - 1))
      nproc = 0
      DO i = 0, pdims(1) - 1
         DO j = 0, pdims(2) - 1
            pgrid(i, j) = nproc
            nproc = nproc + 1
         END DO
      END DO

      CALL dbcsr_distribution_new(dbcsr_dist, group=para_env%get_handle(), pgrid=pgrid, row_dist=row_dist, col_dist=col_dist)

      !The temporary DBCSR integrals and derivatives matrices in this flat distribution
      CALL dbcsr_create(tmp_G_PQ, template=G_PQ, matrix_type=dbcsr_type_no_symmetry, dist=dbcsr_dist)
      CALL dbcsr_complete_redistribute(G_PQ, tmp_G_PQ)

      ALLOCATE (basis_set_ri_aux(nkind), sizes_RI(natom))
      CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_ri_aux)
      n_RI = SUM(sizes_RI)

      one_proc_group = mp2_env%mp2_num_proc == 1
      ALLOCATE (para_env_ext)
      IF (one_proc_group) THEN
         !one subgroup per proc
         CALL para_env_ext%from_split(para_env, para_env%mepos)
      ELSE
         !Split the communicator accross the columns of the matrix
         ncoms = MIN(pdims(2), para_env%num_pe/mp2_env%mp2_num_proc)
         DO i = 0, pdims(1) - 1
            DO j = 0, pdims(2) - 1
               IF (pgrid(i, j) == para_env%mepos) color = MODULO(j + 1, ncoms)
            END DO
         END DO
         CALL para_env_ext%from_split(para_env, color)
      END IF

      !sab_orb and task_list_ext are essentially dummies
      CALL prepare_gpw(qs_env, dft_control, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_ext, pw_env_ext, &
                       auxbas_pw_pool, poisson_env, task_list_ext, rho_r, rho_g, pot_g, psi_L, sab_orb)

      IF (use_virial) THEN
         CALL auxbas_pw_pool%create_pw(rho_g_copy, &
                                       use_data=COMPLEXDATA1D, &
                                       in_space=RECIPROCALSPACE)
         DO i_xyz = 1, 3
            CALL auxbas_pw_pool%create_pw(dvg(i_xyz), &
                                          use_data=COMPLEXDATA1D, &
                                          in_space=RECIPROCALSPACE)
         END DO
      END IF

      ALLOCATE (wf_vector(n_RI))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

      ALLOCATE (iproc_map(natom))

      !Loop over the atomic blocks
      DO jatom = 1, natom

         !Only calculate if on the correct sub-communicator/proc
         IF (one_proc_group) THEN
            iproc_map = 0
            DO iatom = 1, natom
               IF (pgrid(row_dist(iatom), col_dist(jatom)) == para_env%mepos) iproc_map(iatom) = 1
            END DO
            IF (.NOT. ANY(iproc_map == 1)) CYCLE
         ELSE
            IF (.NOT. MODULO(col_dist(jatom) + 1, ncoms) == color) CYCLE
         END IF

         lb_RI = SUM(sizes_RI(1:jatom - 1))
         ub_RI = lb_RI + sizes_RI(jatom)
         DO j_RI = lb_RI + 1, ub_RI

            wf_vector = 0.0_dp
            wf_vector(j_RI) = 1.0_dp

            CALL calculate_wavefunction(mos(1)%mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env_ext, &
                                        basis_type="RI_AUX", external_vector=wf_vector)

            IF (use_virial) THEN
               CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, mp2_env%potential_parameter, dvg)

               wf_vector = 0.0_dp
               DO iatom = 1, natom
                  !only compute if i,j atom pair on correct proc
                  IF (one_proc_group) THEN
                     IF (.NOT. iproc_map(iatom) == 1) CYCLE
                  END IF

                  CALL dbcsr_get_block_p(tmp_G_PQ, iatom, jatom, pblock, found)
                  IF (.NOT. found) CYCLE

                  i_RI = SUM(sizes_RI(1:iatom - 1))
                  wf_vector(i_RI + 1:i_RI + sizes_RI(iatom)) = pblock(:, j_RI - lb_RI)
               END DO

               CALL pw_copy(rho_g, rho_g_copy)
               CALL calculate_wavefunction(mos(1)%mo_coeff, 1, psi_L, rho_g, atomic_kind_set, &
                                           qs_kind_set, cell, dft_control, particle_set, pw_env_ext, &
                                           basis_type="RI_AUX", external_vector=wf_vector)

               CALL calc_potential_gpw(psi_L, rho_g, poisson_env, pot_g, mp2_env%potential_parameter, &
                                       no_transfer=.TRUE.)
               CALL virial_gpw_potential(rho_g_copy, pot_g, rho_g, dvg, h_stress, &
                                         mp2_env%potential_parameter, para_env_ext)
            ELSE
               CALL calc_potential_gpw(rho_r, rho_g, poisson_env, pot_g, mp2_env%potential_parameter)
            END IF

            NULLIFY (rs_v)
            CALL pw_env_get(pw_env_ext, rs_grids=rs_v)
            CALL potential_pw2rs(rs_v, rho_r, pw_env_ext)

            DO iatom = 1, natom

               !only compute if i,j atom pair on correct proc
               IF (one_proc_group) THEN
                  IF (.NOT. iproc_map(iatom) == 1) CYCLE
               END IF

               force_a(:) = 0.0_dp
               force_b(:) = 0.0_dp
               IF (use_virial) THEN
                  my_virial_a = 0.0_dp
                  my_virial_b = 0.0_dp
               END IF

               ikind = kind_of(iatom)
               atom_a = atom_of_kind(iatom)

               basis_set_a => basis_set_ri_aux(ikind)%gto_basis_set
               first_sgfa => basis_set_a%first_sgf
               la_max => basis_set_a%lmax
               la_min => basis_set_a%lmin
               nseta = basis_set_a%nset
               nsgfa => basis_set_a%nsgf_set
               sphi_a => basis_set_a%sphi
               zeta => basis_set_a%zet
               npgfa => basis_set_a%npgf

               ra(:) = pbc(particle_set(iatom)%r, cell)

               CALL dbcsr_get_block_p(tmp_G_PQ, iatom, jatom, pblock, found)
               IF (.NOT. found) CYCLE

               offset = 0
               DO iset = 1, nseta
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  sgfa = first_sgfa(1, iset)

                  ALLOCATE (h_tmp(ncoa, 1)); h_tmp = 0.0_dp
                  ALLOCATE (I_ab(nsgfa(iset), 1)); I_ab = 0.0_dp
                  ALLOCATE (pab(ncoa, 1)); pab = 0.0_dp

                  I_ab(1:nsgfa(iset), 1) = 2.0_dp*pblock(offset + 1:offset + nsgfa(iset), j_RI - lb_RI)
                  CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), 1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             I_ab(1, 1), nsgfa(iset), 0.0_dp, pab(1, 1), ncoa)

                  igrid_level = gaussian_gridlevel(pw_env_ext%gridlevel_info, MINVAL(zeta(:, iset)))

                  ! The last three parameters are used to check whether a given function is within the own range.
                  ! Here, it is always the case, so let's enforce it because mod(0, 1)==0
                  IF (map_gaussian_here(rs_v(igrid_level), cell%h_inv, ra, 0, 1, 0)) THEN
                     DO ipgf = 1, npgfa(iset)
                        o1 = (ipgf - 1)*ncoset(la_max(iset))
                        igrid_level = gaussian_gridlevel(pw_env_ext%gridlevel_info, zeta(ipgf, iset))

                        radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                                          lb_min=0, lb_max=0, ra=ra, rb=ra, rp=ra, &
                                                          zetp=zeta(ipgf, iset), &
                                                          eps=dft_control%qs_control%eps_gvg_rspace, &
                                                          prefactor=1.0_dp, cutoff=1.0_dp)

                        CALL integrate_pgf_product( &
                           la_max=la_max(iset), zeta=zeta(ipgf, iset), la_min=la_min(iset), &
                           lb_max=0, zetb=0.0_dp, lb_min=0, &
                           ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                           rsgrid=rs_v(igrid_level), &
                           hab=h_tmp, pab=pab, &
                           o1=o1, &
                           o2=0, &
                           radius=radius, &
                           calculate_forces=.TRUE., &
                           force_a=force_a, force_b=force_b, &
                           use_virial=use_virial, my_virial_a=my_virial_a, my_virial_b=my_virial_b)

                     END DO

                  END IF

                  offset = offset + nsgfa(iset)
                  DEALLOCATE (pab, h_tmp, I_ab)
               END DO !iset

               force(ikind)%mp2_non_sep(:, atom_a) = force(ikind)%mp2_non_sep(:, atom_a) + force_a + force_b
               IF (use_virial) h_stress = h_stress + my_virial_a + my_virial_b

            END DO !iatom
         END DO !j_RI
      END DO !jatom

      IF (use_virial) THEN
         CALL auxbas_pw_pool%give_back_pw(rho_g_copy)
         DO i_xyz = 1, 3
            CALL auxbas_pw_pool%give_back_pw(dvg(i_xyz))
         END DO
      END IF

      CALL cleanup_gpw(qs_env, e_cutoff_old, cutoff_old, relative_cutoff_old, para_env_ext, pw_env_ext, &
                       task_list_ext, auxbas_pw_pool, rho_r, rho_g, pot_g, psi_L)

      CALL dbcsr_release(tmp_G_PQ)
      CALL dbcsr_distribution_release(dbcsr_dist)
      DEALLOCATE (col_dist, row_dist, pgrid)

      CALL mp_para_env_release(para_env_ext)

      CALL timestop(handle)

   END SUBROUTINE get_2c_gpw_forces

! **************************************************************************************************
!> \brief Calculate the forces due to the (P|Q) MME integral derivatives
!> \param G_PQ ...
!> \param force ...
!> \param mp2_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_2c_mme_forces(G_PQ, force, mp2_env, qs_env)

      TYPE(dbcsr_type), INTENT(INOUT)                    :: G_PQ
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(mp2_type), INTENT(INOUT)                      :: mp2_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_2c_mme_forces'

      INTEGER :: atom_a, atom_b, blk, G_count, handle, i_xyz, iatom, ikind, iset, jatom, jkind, &
         jset, natom, nkind, nseta, nsetb, offset_hab_a, offset_hab_b, R_count, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(dp)                                           :: new_force, pref
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: hab
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: hdab
      REAL(dp), DIMENSION(:, :), POINTER                 :: pblock
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sphi_a, sphi_b, zeta, zetb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, basis_set_a, basis_set_b, pblock, particle_set, &
               cell, la_max, la_min, lb_min, npgfa, lb_max, npgfb, nsgfa, &
               nsgfb, first_sgfa, first_sgfb, sphi_a, sphi_b, zeta, zetb, &
               atomic_kind_set, para_env)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, nkind=nkind, particle_set=particle_set, &
                      cell=cell, atomic_kind_set=atomic_kind_set, natom=natom, para_env=para_env)

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

      ALLOCATE (basis_set_ri_aux(nkind))
      CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)

      G_count = 0; R_count = 0

      CALL dbcsr_iterator_start(iter, G_PQ)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row=iatom, column=jatom, blk=blk)
         CALL dbcsr_get_block_p(G_PQ, iatom, jatom, pblock, found)
         IF (.NOT. found) CYCLE
         IF (iatom > jatom) CYCLE
         pref = 2.0_dp
         IF (iatom == jatom) pref = 1.0_dp

         ikind = kind_of(iatom)
         jkind = kind_of(jatom)

         atom_a = atom_of_kind(iatom)
         atom_b = atom_of_kind(jatom)

         basis_set_a => basis_set_ri_aux(ikind)%gto_basis_set
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         npgfa => basis_set_a%npgf

         basis_set_b => basis_set_ri_aux(jkind)%gto_basis_set
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet
         npgfb => basis_set_b%npgf

         ra(:) = pbc(particle_set(iatom)%r, cell)
         rb(:) = pbc(particle_set(jatom)%r, cell)

         ALLOCATE (hab(basis_set_a%nsgf, basis_set_b%nsgf))
         ALLOCATE (hdab(3, basis_set_a%nsgf, basis_set_b%nsgf))
         hab(:, :) = 0.0_dp
         hdab(:, :, :) = 0.0_dp

         offset_hab_a = 0
         DO iset = 1, nseta
            sgfa = first_sgfa(1, iset)

            offset_hab_b = 0
            DO jset = 1, nsetb
               sgfb = first_sgfb(1, jset)

               CALL integrate_set_2c(mp2_env%eri_mme_param%par, mp2_env%potential_parameter, la_min(iset), &
                                     la_max(iset), lb_min(jset), lb_max(jset), npgfa(iset), npgfb(jset), &
                                     zeta(:, iset), zetb(:, jset), ra, rb, hab, nsgfa(iset), nsgfb(jset), &
                                     offset_hab_a, offset_hab_b, 0, 0, sphi_a, sphi_b, sgfa, sgfb, &
                                     nsgfa(iset), nsgfb(jset), do_eri_mme, hdab=hdab, &
                                     G_count=G_count, R_count=R_count)

               offset_hab_b = offset_hab_b + nsgfb(jset)
            END DO
            offset_hab_a = offset_hab_a + nsgfa(iset)
         END DO

         DO i_xyz = 1, 3
            new_force = pref*SUM(pblock(:, :)*hdab(i_xyz, :, :))
            force(ikind)%mp2_non_sep(i_xyz, atom_a) = force(ikind)%mp2_non_sep(i_xyz, atom_a) + new_force
            force(jkind)%mp2_non_sep(i_xyz, atom_b) = force(jkind)%mp2_non_sep(i_xyz, atom_b) - new_force
         END DO

         DEALLOCATE (hab, hdab)
      END DO
      CALL dbcsr_iterator_stop(iter)

      CALL cp_eri_mme_update_local_counts(mp2_env%eri_mme_param, para_env, G_count_2c=G_count, R_count_2c=R_count)

      CALL timestop(handle)

   END SUBROUTINE get_2c_mme_forces

! **************************************************************************************************
!> \brief This routines gather all the force updates due to the response density and the trace with F
!>        Also update the forces due to the SCF density for XC and exact exchange
!> \param p_env the p_env coming from the response calculation
!> \param matrix_hz the matrix going into the RHS of the response equation
!> \param matrix_p_F the density matrix with which we evaluate Trace[P*F]
!> \param matrix_p_F_admm ...
!> \param qs_env ...
!> \note very much inspired from the response_force routine in response_solver.F, especially for virial
! **************************************************************************************************
   SUBROUTINE update_im_time_forces(p_env, matrix_hz, matrix_p_F, matrix_p_F_admm, qs_env)

      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz, matrix_p_F, matrix_p_F_admm
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'update_im_time_forces'

      INTEGER                                            :: handle, i, idens, ispin, n_rep_hf, nao, &
                                                            nao_aux, nder, nimages, nocc, nspins
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_exx, do_hfx, do_tau, do_tau_admm, &
                                                            use_virial
      REAL(dp)                                           :: dummy_real1, dummy_real2, ehartree, &
                                                            eps_ppnl, exc, focc
      REAL(dp), DIMENSION(3, 3)                          :: h_stress, pv_loc
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: current_density, current_density_admm, &
         current_mat_h, matrix_p_mp2, matrix_p_mp2_admm, matrix_s, matrix_s_aux_fit, matrix_w, &
         rho_ao, rho_ao_aux, scrm, scrm_admm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: dbcsr_work_h, dbcsr_work_p
      TYPE(dbcsr_type)                                   :: dbcsr_work
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: rho_tot_gspace, rhoz_tot_gspace, &
                                                            vh_rspace, vhxc_rspace, &
                                                            zv_hartree_gspace, zv_hartree_rspace
      TYPE(pw_type), DIMENSION(:), POINTER               :: rhoz_g, rhoz_r, tauz_r, v_xc, v_xc_tau, &
                                                            vadmm_rspace, vtau_rspace, vxc_rspace
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (scrm, rho, dft_control, matrix_p_mp2, matrix_s, matrix_p_mp2_admm, admm_env, sab_orb, &
               cell_to_index, dbcsr_work_p, dbcsr_work_h, sac_ppl, sap_ppnl, force, virial, &
               qs_kind_set, atomic_kind_set, particle_set, pw_env, poisson_env, auxbas_pw_pool, &
               task_list_aux_fit, matrix_s_aux_fit, scrm_admm, rho_aux_fit, rho_ao_aux, x_data, &
               hfx_section, xc_section, para_env, rhoz_g, rhoz_r, tauz_r, v_xc, v_xc_tau, &
               vxc_rspace, vtau_rspace, vadmm_rspace, rho_ao, matrix_w)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, rho=rho, dft_control=dft_control, matrix_s=matrix_s, admm_env=admm_env, &
                      sab_orb=sab_orb, sac_ppl=sac_ppl, sap_ppnl=sap_ppnl, force=force, virial=virial, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                      x_data=x_data, para_env=para_env)
      nspins = dft_control%nspins

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (use_virial) virial%pv_calculate = .TRUE.

      !Whether we replace the force/energy of SCF XC with HF in RPA
      do_exx = .FALSE.
      IF (qs_env%mp2_env%method == ri_rpa_method_gpw) THEN
         hfx_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
         CALL section_vals_get(hfx_section, explicit=do_exx)
      END IF

      !Get the mp2 density matrix which is p_env%p1 + matrix_p_F
      CALL get_qs_env(qs_env, matrix_p_mp2=matrix_p_mp2, matrix_p_mp2_admm=matrix_p_mp2_admm)

      !The kinetic term (only response density)
      NULLIFY (scrm)
      IF (nspins == 2) CALL dbcsr_add(matrix_p_mp2(1)%matrix, matrix_p_mp2(2)%matrix, 1.0_dp, 1.0_dp)
      CALL build_kinetic_matrix(qs_env%ks_env, matrix_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_p_mp2(1)%matrix)
      IF (nspins == 2) CALL dbcsr_add(matrix_p_mp2(1)%matrix, matrix_p_mp2(2)%matrix, 1.0_dp, -1.0_dp)
      CALL dbcsr_deallocate_matrix_set(scrm)

      !The pseudo-potential terms (only reponse density)
      CALL dbcsr_allocate_matrix_set(scrm, nspins)
      DO ispin = 1, nspins
         ALLOCATE (scrm(ispin)%matrix)
         CALL dbcsr_create(scrm(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_copy(scrm(ispin)%matrix, matrix_s(1)%matrix)
         CALL dbcsr_set(scrm(ispin)%matrix, 0.0_dp)
      END DO

      nder = 1
      nimages = 1
      NULLIFY (cell_to_index)
      ALLOCATE (dbcsr_work_p(nspins, 1), dbcsr_work_h(nspins, 1))
      DO ispin = 1, nspins
         dbcsr_work_p(ispin, 1)%matrix => matrix_p_mp2(ispin)%matrix
         dbcsr_work_h(ispin, 1)%matrix => scrm(ispin)%matrix
      END DO

      IF (ASSOCIATED(sac_ppl)) THEN
         CALL build_core_ppl(dbcsr_work_h, dbcsr_work_p, force, &
                             virial, .TRUE., use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, &
                             sab_orb, sac_ppl, nimages, cell_to_index, "ORB")
      END IF

      IF (ASSOCIATED(sap_ppnl)) THEN
         eps_ppnl = dft_control%qs_control%eps_ppnl
         CALL build_core_ppnl(dbcsr_work_h, dbcsr_work_p, force, &
                              virial, .TRUE., use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, &
                              sab_orb, sap_ppnl, eps_ppnl, nimages, cell_to_index, "ORB")
      END IF
      DEALLOCATE (dbcsr_work_p, dbcsr_work_h)

      IF (use_virial) THEN
         h_stress = 0.0_dp
         virial%pv_xc = 0.0_dp
         NULLIFY (vxc_rspace, vtau_rspace, vadmm_rspace)
         CALL ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, &
                               dummy_real1, dummy_real2, h_stress)
         virial%pv_ehartree = virial%pv_ehartree + h_stress/REAL(para_env%num_pe, dp)
         virial%pv_virial = virial%pv_virial + h_stress/REAL(para_env%num_pe, dp)
         IF (.NOT. do_exx) THEN
            !if RPA EXX, then do not consider XC virial (replaced by RPA%HF virial)
            virial%pv_exc = virial%pv_exc - virial%pv_xc
            virial%pv_virial = virial%pv_virial - virial%pv_xc
         END IF
      ELSE
         CALL ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, vadmm_rspace, dummy_real1, dummy_real2)
      END IF
      do_tau = ASSOCIATED(vtau_rspace)

      !Core forces from the SCF
      CALL integrate_v_core_rspace(vh_rspace, qs_env)

      !The Hartree-xc potential term, P*dVHxc (mp2 + SCF density x deriv of the SCF potential)
      !Get the total density
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      DO ispin = 1, nspins
         CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, 1.0_dp)
      END DO

      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL auxbas_pw_pool%create_pw(vhxc_rspace, &
                                    use_data=REALDATA3D, in_space=REALSPACE)

      IF (use_virial) pv_loc = virial%pv_virial

      IF (do_exx) THEN
         !Only want response XC contribution, but SCF+response Hartree contribution
         DO ispin = 1, nspins
            !Hartree
            CALL pw_transfer(vh_rspace, vhxc_rspace)
            CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                    hmat=scrm(ispin), pmat=rho_ao(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE.)
            !XC
            CALL pw_transfer(vxc_rspace(ispin), vhxc_rspace)
            CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                    hmat=scrm(ispin), pmat=matrix_p_mp2(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE.)
            IF (do_tau) THEN
               CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                       hmat=scrm(ispin), pmat=matrix_p_mp2(ispin), &
                                       qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE.)
            END IF
         END DO
      ELSE
         DO ispin = 1, nspins
            CALL pw_transfer(vh_rspace, vhxc_rspace)
            CALL pw_axpy(vxc_rspace(ispin), vhxc_rspace)
            CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                    hmat=scrm(ispin), pmat=rho_ao(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE.)
            IF (do_tau) THEN
               CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                       hmat=scrm(ispin), pmat=rho_ao(ispin), &
                                       qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE.)
            END IF
         END DO
      END IF
      CALL auxbas_pw_pool%give_back_pw(vhxc_rspace)

      IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)

      !The admm projection contribution (mp2 + SCF densities). If EXX, then only mp2 density
      IF (dft_control%do_admm) THEN
         CALL get_admm_env(admm_env, task_list_aux_fit=task_list_aux_fit, rho_aux_fit=rho_aux_fit, &
                           matrix_s_aux_fit=matrix_s_aux_fit)
         CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux)
         CALL dbcsr_allocate_matrix_set(scrm_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (scrm_admm(ispin)%matrix)
            CALL dbcsr_create(scrm_admm(ispin)%matrix, template=matrix_s_aux_fit(1)%matrix)
            CALL dbcsr_copy(scrm_admm(ispin)%matrix, matrix_s_aux_fit(1)%matrix)
            CALL dbcsr_set(scrm_admm(ispin)%matrix, 0.0_dp)
         END DO

         IF (use_virial) pv_loc = virial%pv_virial
         IF (.NOT. qs_env%admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
            DO ispin = 1, nspins
               IF (do_exx) THEN
                  CALL integrate_v_rspace(v_rspace=vadmm_rspace(ispin), &
                                          hmat=scrm_admm(ispin), pmat=matrix_p_mp2_admm(ispin), &
                                          qs_env=qs_env, calculate_forces=.TRUE., &
                                          basis_type="AUX_FIT", task_list_external=task_list_aux_fit)
               ELSE
                  CALL dbcsr_add(rho_ao_aux(ispin)%matrix, matrix_p_mp2_admm(ispin)%matrix, 1.0_dp, 1.0_dp)
                  CALL integrate_v_rspace(v_rspace=vadmm_rspace(ispin), &
                                          hmat=scrm_admm(ispin), pmat=rho_ao_aux(ispin), &
                                          qs_env=qs_env, calculate_forces=.TRUE., &
                                          basis_type="AUX_FIT", task_list_external=task_list_aux_fit)
                  CALL dbcsr_add(rho_ao_aux(ispin)%matrix, matrix_p_mp2_admm(ispin)%matrix, 1.0_dp, -1.0_dp)
               END IF
            END DO
         END IF
         IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)

         CALL tddft_hfx_matrix(scrm_admm, rho_ao_aux, qs_env, .FALSE., .FALSE.)

         IF (do_exx) THEN
            CALL admm_projection_derivative(qs_env, scrm_admm, matrix_p_mp2)
         ELSE
            CALL admm_projection_derivative(qs_env, scrm_admm, rho_ao)
         END IF
      END IF

      !The exact-exchange term (mp2 + SCF densities)
      xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      hfx_section => section_vals_get_subs_vals(xc_section, "HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)

      IF (do_hfx) THEN
         CALL section_vals_get(hfx_section, n_repetition=n_rep_hf)
         CPASSERT(n_rep_hf == 1)
         IF (use_virial) virial%pv_fock_4c = 0.0_dp

         !In case of EXX, only want to response HFX forces, as the SCF will change according to RI_RPA%HF
         IF (do_exx) THEN
            IF (dft_control%do_admm) THEN
               CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit)
               CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux, rho_ao_kp=dbcsr_work_p)
               IF (x_data(1, 1)%do_hfx_ri) THEN

                  CALL hfx_ri_update_forces(qs_env, x_data(1, 1)%ri_data, nspins, &
                                            x_data(1, 1)%general_parameter%fraction, &
                                            rho_ao=dbcsr_work_p, rho_ao_resp=matrix_p_mp2_admm, &
                                            use_virial=use_virial, resp_only=.TRUE.)
               ELSE
                  CALL derivatives_four_center(qs_env, dbcsr_work_p, matrix_p_mp2_admm, hfx_section, para_env, &
                                               1, use_virial, resp_only=.TRUE.)
               END IF
            ELSE
               DO ispin = 1, nspins
                  CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, -1.0_dp)
               END DO
               CALL qs_rho_get(rho, rho_ao_kp=dbcsr_work_p)
               IF (x_data(1, 1)%do_hfx_ri) THEN

                  CALL hfx_ri_update_forces(qs_env, x_data(1, 1)%ri_data, nspins, &
                                            x_data(1, 1)%general_parameter%fraction, &
                                            rho_ao=dbcsr_work_p, rho_ao_resp=matrix_p_mp2, &
                                            use_virial=use_virial, resp_only=.TRUE.)
               ELSE
                  CALL derivatives_four_center(qs_env, dbcsr_work_p, matrix_p_mp2, hfx_section, para_env, &
                                               1, use_virial, resp_only=.TRUE.)
               END IF
               DO ispin = 1, nspins
                  CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, 1.0_dp)
               END DO
            END IF !admm

         ELSE !No Exx
            IF (dft_control%do_admm) THEN
               CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit)
               CALL qs_rho_get(rho_aux_fit, rho_ao=rho_ao_aux, rho_ao_kp=dbcsr_work_p)
               DO ispin = 1, nspins
                  CALL dbcsr_add(rho_ao_aux(ispin)%matrix, matrix_p_mp2_admm(ispin)%matrix, 1.0_dp, 1.0_dp)
               END DO
               IF (x_data(1, 1)%do_hfx_ri) THEN

                  CALL hfx_ri_update_forces(qs_env, x_data(1, 1)%ri_data, nspins, &
                                            x_data(1, 1)%general_parameter%fraction, &
                                            rho_ao=dbcsr_work_p, rho_ao_resp=matrix_p_mp2_admm, &
                                            use_virial=use_virial, resp_only=.FALSE.)
               ELSE
                  CALL derivatives_four_center(qs_env, dbcsr_work_p, matrix_p_mp2_admm, hfx_section, para_env, &
                                               1, use_virial, resp_only=.FALSE.)
               END IF
               DO ispin = 1, nspins
                  CALL dbcsr_add(rho_ao_aux(ispin)%matrix, matrix_p_mp2_admm(ispin)%matrix, 1.0_dp, -1.0_dp)
               END DO
            ELSE
               CALL qs_rho_get(rho, rho_ao_kp=dbcsr_work_p)
               IF (x_data(1, 1)%do_hfx_ri) THEN

                  CALL hfx_ri_update_forces(qs_env, x_data(1, 1)%ri_data, nspins, &
                                            x_data(1, 1)%general_parameter%fraction, &
                                            rho_ao=dbcsr_work_p, rho_ao_resp=matrix_p_mp2, &
                                            use_virial=use_virial, resp_only=.FALSE.)
               ELSE
                  CALL derivatives_four_center(qs_env, dbcsr_work_p, matrix_p_mp2, hfx_section, para_env, &
                                               1, use_virial, resp_only=.FALSE.)
               END IF
            END IF
         END IF !do_exx

         IF (use_virial) THEN
            virial%pv_exx = virial%pv_exx - virial%pv_fock_4c
            virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
         END IF
      END IF

      !retrieve the SCF density
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      DO ispin = 1, nspins
         CALL dbcsr_add(rho_ao(ispin)%matrix, matrix_p_mp2(ispin)%matrix, 1.0_dp, -1.0_dp)
      END DO

      !From here, we need to do everything twice. Once for the response density, and once for the
      !density that is used for the trace Tr[P*F]. The reason is that the former is needed for the
      !eventual overlap contribution from matrix_wz
      !Only with the mp2 density

      ALLOCATE (current_density(nspins), current_mat_h(nspins), current_density_admm(nspins))
      DO idens = 1, 2
         DO ispin = 1, nspins
            IF (idens == 1) THEN
               current_density(ispin)%matrix => matrix_p_F(ispin)%matrix
               current_mat_h(ispin)%matrix => scrm(ispin)%matrix
               IF (dft_control%do_admm) current_density_admm(ispin)%matrix => matrix_p_F_admm(ispin)%matrix
            ELSE
               current_density(ispin)%matrix => p_env%p1(ispin)%matrix
               current_mat_h(ispin)%matrix => matrix_hz(ispin)%matrix
               IF (dft_control%do_admm) current_density_admm(ispin)%matrix => p_env%p1_admm(ispin)%matrix
            END IF
         END DO

         !The core-denstiy derivative
         ALLOCATE (rhoz_r(nspins), rhoz_g(nspins))
         DO ispin = 1, nspins
            CALL auxbas_pw_pool%create_pw(rhoz_r(ispin), &
                                          use_data=REALDATA3D, in_space=REALSPACE)
            CALL auxbas_pw_pool%create_pw(rhoz_g(ispin), &
                                          use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         END DO
         CALL auxbas_pw_pool%create_pw(rhoz_tot_gspace, &
                                       use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL auxbas_pw_pool%create_pw(zv_hartree_rspace, &
                                       use_data=REALDATA3D, in_space=REALSPACE)
         CALL auxbas_pw_pool%create_pw(zv_hartree_gspace, &
                                       use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

         CALL pw_zero(rhoz_tot_gspace)
         DO ispin = 1, nspins
            CALL calculate_rho_elec(ks_env=qs_env%ks_env, matrix_p=current_density(ispin)%matrix, &
                                    rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin))
            CALL pw_axpy(rhoz_g(ispin), rhoz_tot_gspace)
         END DO

         IF (use_virial) THEN

            CALL get_qs_env(qs_env, rho=rho)
            CALL auxbas_pw_pool%create_pw(rho_tot_gspace, &
                                          use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

            CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

            h_stress(:, :) = 0.0_dp
            CALL pw_poisson_solve(poisson_env, &
                                  density=rhoz_tot_gspace, &
                                  ehartree=ehartree, &
                                  vhartree=zv_hartree_gspace, &
                                  h_stress=h_stress, &
                                  aux_density=rho_tot_gspace)

            CALL auxbas_pw_pool%give_back_pw(rho_tot_gspace)

            !Green contribution
            virial%pv_ehartree = virial%pv_ehartree + 2.0_dp*h_stress/REAL(para_env%num_pe, dp)
            virial%pv_virial = virial%pv_virial + 2.0_dp*h_stress/REAL(para_env%num_pe, dp)

         ELSE
            CALL pw_poisson_solve(poisson_env, rhoz_tot_gspace, ehartree, &
                                  zv_hartree_gspace)
         END IF

         CALL pw_transfer(zv_hartree_gspace, zv_hartree_rspace)
         CALL pw_scale(zv_hartree_rspace, zv_hartree_rspace%pw_grid%dvol)
         CALL integrate_v_core_rspace(zv_hartree_rspace, qs_env)

         IF (do_tau) THEN
            BLOCK
               TYPE(pw_type) :: tauz_g
               CALL auxbas_pw_pool%create_pw(tauz_g, &
                                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
               ALLOCATE (tauz_r(nspins))
               DO ispin = 1, nspins
                  CALL auxbas_pw_pool%create_pw(tauz_r(ispin), &
                                                use_data=REALDATA3D, in_space=REALSPACE)

                  CALL calculate_rho_elec(ks_env=qs_env%ks_env, matrix_p=current_density(ispin)%matrix, &
                                          rho=tauz_r(ispin), rho_gspace=tauz_g, compute_tau=.TRUE.)
               END DO
               CALL auxbas_pw_pool%give_back_pw(tauz_g)
            END BLOCK
         END IF

         !Volume contribution to the virial
         IF (use_virial) THEN
            !Volume contribution
            exc = 0.0_dp
            DO ispin = 1, nspins
               exc = exc + pw_integral_ab(rhoz_r(ispin), vxc_rspace(ispin))/ &
                     vxc_rspace(ispin)%pw_grid%dvol
            END DO
            IF (ASSOCIATED(vtau_rspace)) THEN
               DO ispin = 1, nspins
                  exc = exc + pw_integral_ab(tauz_r(ispin), vtau_rspace(ispin))/ &
                        vtau_rspace(ispin)%pw_grid%dvol
               END DO
            END IF
            DO i = 1, 3
               virial%pv_ehartree(i, i) = virial%pv_ehartree(i, i) - 4.0_dp*ehartree/REAL(para_env%num_pe, dp)
               virial%pv_exc(i, i) = virial%pv_exc(i, i) - exc/REAL(para_env%num_pe, dp)
               virial%pv_virial(i, i) = virial%pv_virial(i, i) - 4.0_dp*ehartree/REAL(para_env%num_pe, dp) &
                                        - exc/REAL(para_env%num_pe, dp)
            END DO
         END IF

         !The xc-kernel term.
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            xc_section => admm_env%xc_section_primary
         ELSE
            xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
         END IF

         IF (use_virial) virial%pv_xc = 0.0_dp

         CALL create_kernel(qs_env, &
                            vxc=v_xc, &
                            vxc_tau=v_xc_tau, &
                            rho=rho, &
                            rho1_r=rhoz_r, &
                            rho1_g=rhoz_g, &
                            tau1_r=tauz_r, &
                            xc_section=xc_section, &
                            compute_virial=use_virial, &
                            virial_xc=virial%pv_xc)

         IF (use_virial) THEN
            virial%pv_exc = virial%pv_exc + virial%pv_xc
            virial%pv_virial = virial%pv_virial + virial%pv_xc

            pv_loc = virial%pv_virial
         END IF

         CALL qs_rho_get(rho, rho_ao_kp=dbcsr_work_p)
         DO ispin = 1, nspins
            CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
            CALL pw_axpy(zv_hartree_rspace, v_xc(ispin))
            CALL integrate_v_rspace(qs_env=qs_env, &
                                    v_rspace=v_xc(ispin), &
                                    hmat=current_mat_h(ispin), &
                                    pmat=dbcsr_work_p(ispin, 1), &
                                    calculate_forces=.TRUE.)
            CALL auxbas_pw_pool%give_back_pw(v_xc(ispin))
         END DO
         CALL auxbas_pw_pool%give_back_pw(rhoz_tot_gspace)
         CALL auxbas_pw_pool%give_back_pw(zv_hartree_rspace)
         CALL auxbas_pw_pool%give_back_pw(zv_hartree_gspace)
         DEALLOCATE (v_xc)

         IF (do_tau) THEN
            DO ispin = 1, nspins
               CALL pw_scale(v_xc_tau(ispin), v_xc_tau(ispin)%pw_grid%dvol)
               CALL integrate_v_rspace(qs_env=qs_env, &
                                       v_rspace=v_xc_tau(ispin), &
                                       hmat=current_mat_h(ispin), &
                                       pmat=dbcsr_work_p(ispin, 1), &
                                       compute_tau=.TRUE., &
                                       calculate_forces=.TRUE.)
               CALL auxbas_pw_pool%give_back_pw(v_xc_tau(ispin))
            END DO
            DEALLOCATE (v_xc_tau)
         END IF

         IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)

         IF (do_hfx) THEN
            IF (dft_control%do_admm) THEN
               DO ispin = 1, nspins
                  CALL dbcsr_set(scrm_admm(ispin)%matrix, 0.0_dp)
               END DO
               CALL qs_rho_get(rho_aux_fit, tau_r_valid=do_tau_admm)

               IF (.NOT. admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
                  CALL get_admm_env(admm_env, rho_aux_fit=rho_aux_fit)
                  DO ispin = 1, nspins
                     CALL pw_zero(rhoz_r(ispin))
                     CALL pw_zero(rhoz_g(ispin))
                     CALL calculate_rho_elec(ks_env=qs_env%ks_env, matrix_p=current_density_admm(ispin)%matrix, &
                                             rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin), &
                                             basis_type="AUX_FIT", task_list_external=task_list_aux_fit)
                  END DO

                  IF (do_tau_admm) THEN
                     BLOCK
                        TYPE(pw_type) :: tauz_g
                        CALL auxbas_pw_pool%create_pw(tauz_g, &
                                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
                        DO ispin = 1, nspins
                           CALL pw_zero(tauz_r(ispin))
                           CALL calculate_rho_elec(ks_env=qs_env%ks_env, matrix_p=current_density(ispin)%matrix, &
                                                   rho=tauz_r(ispin), rho_gspace=tauz_g, &
                                                   basis_type="AUX_FIT", task_list_external=task_list_aux_fit, &
                                                   compute_tau=.TRUE.)
                        END DO
                        CALL auxbas_pw_pool%give_back_pw(tauz_g)
                     END BLOCK
                  END IF

                  !Volume contribution to the virial
                  IF (use_virial) THEN
                     exc = 0.0_dp
                     DO ispin = 1, nspins
                        exc = exc + pw_integral_ab(rhoz_r(ispin), vadmm_rspace(ispin))/ &
                              vadmm_rspace(ispin)%pw_grid%dvol
                     END DO
                     DO i = 1, 3
                        virial%pv_exc(i, i) = virial%pv_exc(i, i) - exc/REAL(para_env%num_pe, dp)
                        virial%pv_virial(i, i) = virial%pv_virial(i, i) - exc/REAL(para_env%num_pe, dp)
                     END DO

                     virial%pv_xc = 0.0_dp
                  END IF

                  xc_section => admm_env%xc_section_aux
                  CALL create_kernel(qs_env, v_xc, v_xc_tau, rho_aux_fit, rhoz_r, rhoz_g, tauz_r, xc_section, &
                                     compute_virial=use_virial, virial_xc=virial%pv_xc)

                  IF (use_virial) THEN
                     virial%pv_exc = virial%pv_exc + virial%pv_xc
                     virial%pv_virial = virial%pv_virial + virial%pv_xc

                     pv_loc = virial%pv_virial
                  END IF

                  CALL qs_rho_get(rho_aux_fit, rho_ao_kp=dbcsr_work_p)
                  DO ispin = 1, nspins
                     CALL pw_scale(v_xc(ispin), v_xc(ispin)%pw_grid%dvol)
                     CALL integrate_v_rspace(qs_env=qs_env, &
                                             v_rspace=v_xc(ispin), &
                                             hmat=scrm_admm(ispin), &
                                             pmat=dbcsr_work_p(ispin, 1), &
                                             calculate_forces=.TRUE., &
                                             basis_type="AUX_FIT", &
                                             task_list_external=task_list_aux_fit)
                     CALL auxbas_pw_pool%give_back_pw(v_xc(ispin))
                  END DO
                  DEALLOCATE (v_xc)

                  IF (do_tau_admm) THEN
                     DO ispin = 1, nspins
                        CALL pw_scale(v_xc_tau(ispin), v_xc_tau(ispin)%pw_grid%dvol)
                        CALL integrate_v_rspace(qs_env=qs_env, &
                                                v_rspace=v_xc_tau(ispin), &
                                                hmat=scrm_admm(ispin), &
                                                pmat=dbcsr_work_p(ispin, 1), &
                                                calculate_forces=.TRUE., &
                                                basis_type="AUX_FIT", &
                                                task_list_external=task_list_aux_fit, &
                                                compute_tau=.TRUE.)
                        CALL auxbas_pw_pool%give_back_pw(v_xc_tau(ispin))
                     END DO
                     DEALLOCATE (v_xc_tau)
                  END IF

                  IF (use_virial) virial%pv_ehartree = virial%pv_ehartree + (virial%pv_virial - pv_loc)
               END IF

               CALL tddft_hfx_matrix(scrm_admm, current_density_admm, qs_env, .FALSE., .FALSE.)

               CALL qs_rho_get(rho, rho_ao_kp=dbcsr_work_p)
               CALL admm_projection_derivative(qs_env, scrm_admm, dbcsr_work_p(:, 1))

               !If response density, need to get matrix_hz contribution
               CALL dbcsr_create(dbcsr_work, template=matrix_s(1)%matrix)
               IF (idens == 2) THEN
                  nao = admm_env%nao_orb
                  nao_aux = admm_env%nao_aux_fit
                  DO ispin = 1, nspins
                     CALL dbcsr_copy(dbcsr_work, matrix_hz(ispin)%matrix)
                     CALL dbcsr_set(dbcsr_work, 0.0_dp)

                     CALL cp_dbcsr_sm_fm_multiply(scrm_admm(ispin)%matrix, admm_env%A, &
                                                  admm_env%work_aux_orb, nao)
                     CALL parallel_gemm('T', 'N', nao, nao, nao_aux, &
                                        1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                                        admm_env%work_orb_orb)
                     CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbcsr_work, keep_sparsity=.TRUE.)
                     CALL dbcsr_add(matrix_hz(ispin)%matrix, dbcsr_work, 1.0_dp, 1.0_dp)
                  END DO
               END IF

               CALL dbcsr_release(dbcsr_work)
            ELSE !no admm

               !Need the contribution to matrix_hz as well
               IF (idens == 2) THEN
                  CALL tddft_hfx_matrix(matrix_hz, current_density, qs_env, .FALSE., .FALSE.)
               END IF
            END IF !admm
         END IF !do_hfx

         DO ispin = 1, nspins
            CALL auxbas_pw_pool%give_back_pw(rhoz_r(ispin))
            CALL auxbas_pw_pool%give_back_pw(rhoz_g(ispin))
         END DO
         DEALLOCATE (rhoz_r, rhoz_g)

         IF (do_tau) THEN
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%give_back_pw(tauz_r(ispin))
            END DO
            DEALLOCATE (tauz_r)
         END IF
      END DO !idens
      CALL dbcsr_deallocate_matrix_set(scrm_admm)

      DEALLOCATE (current_density, current_mat_h, current_density_admm)
      CALL dbcsr_deallocate_matrix_set(scrm)

      !The energy weighted and overlap term. ONLY with the response density
      focc = 2.0_dp
      IF (nspins == 2) focc = 1.0_dp
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin), homo=nocc)
         CALL calculate_whz_matrix(mos(ispin)%mo_coeff, matrix_hz(ispin)%matrix, &
                                   p_env%w1(ispin)%matrix, focc, nocc)
      END DO
      IF (nspins == 2) CALL dbcsr_add(p_env%w1(1)%matrix, p_env%w1(2)%matrix, 1.0_dp, 1.0_dp)

      !Add to it the SCF W matrix, except if EXX (because taken care of by HF response)
      IF (.NOT. do_exx) THEN
         CALL qs_energies_compute_w(qs_env, calc_forces=.TRUE.)
         CALL get_qs_env(qs_env, matrix_w=matrix_w)
         CALL dbcsr_add(p_env%w1(1)%matrix, matrix_w(1)%matrix, 1.0_dp, 1.0_dp)
         IF (nspins == 2) CALL dbcsr_add(p_env%w1(1)%matrix, matrix_w(2)%matrix, 1.0_dp, 1.0_dp)
      END IF

      NULLIFY (scrm)
      CALL build_overlap_matrix(qs_env%ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=p_env%w1(1)%matrix)

      IF (.NOT. do_exx) THEN
         CALL dbcsr_add(p_env%w1(1)%matrix, matrix_w(1)%matrix, 1.0_dp, -1.0_dp)
         IF (nspins == 2) CALL dbcsr_add(p_env%w1(1)%matrix, matrix_w(2)%matrix, 1.0_dp, -1.0_dp)
         DO ispin = 1, nspins
            CALL dbcsr_set(matrix_w(ispin)%matrix, 0.0_dp)
         END DO
      END IF

      IF (nspins == 2) CALL dbcsr_add(p_env%w1(1)%matrix, p_env%w1(2)%matrix, 1.0_dp, -1.0_dp)
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (use_virial) virial%pv_calculate = .FALSE.

      !clean-up
      CALL auxbas_pw_pool%give_back_pw(vh_rspace)

      DO ispin = 1, nspins
         CALL auxbas_pw_pool%give_back_pw(vxc_rspace(ispin))
         IF (ASSOCIATED(vtau_rspace)) THEN
            CALL auxbas_pw_pool%give_back_pw(vtau_rspace(ispin))
         END IF
         IF (ASSOCIATED(vadmm_rspace)) THEN
            CALL auxbas_pw_pool%give_back_pw(vadmm_rspace(ispin))
         END IF
      END DO
      DEALLOCATE (vxc_rspace)
      IF (ASSOCIATED(vtau_rspace)) DEALLOCATE (vtau_rspace)
      IF (ASSOCIATED(vadmm_rspace)) DEALLOCATE (vadmm_rspace)

      CALL timestop(handle)

   END SUBROUTINE update_im_time_forces

! **************************************************************************************************
!> \brief Iteratively builds the matrix Y = sum_k Y_k until convergence, where
!>        Y_k = 1/k*2^n (A/2^n) Y_k-1 + 1/k!*2^n * PR(n) * (A/2^n)^(k-1)
!>        n is chosen such that the norm of A is < 1 (and e^A converges fast)
!>        PR(n) =  e^(A/2^n)*PR(n-1) + PR(n-1)*e^(A/2^n), PR(0) = P*R^T
!> \param Y ...
!> \param A ...
!> \param P ...
!> \param R ...
!> \param filter_eps ...
! **************************************************************************************************
   SUBROUTINE build_Y_matrix(Y, A, P, R, filter_eps)

      TYPE(dbcsr_type), INTENT(OUT)                      :: Y
      TYPE(dbcsr_type), INTENT(INOUT)                    :: A, P, R
      REAL(dp), INTENT(IN)                               :: filter_eps

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'build_Y_matrix'

      INTEGER                                            :: handle, k, n
      REAL(dp)                                           :: norm_scalar, threshold
      TYPE(dbcsr_type)                                   :: A2n, exp_A2n, PRn, work, work2, Yk

      CALL timeset(routineN, handle)

      threshold = 1.0E-16_dp

      !Find n such that norm(A) < 1 and we insure convergence of the exponential
      norm_scalar = dbcsr_frobenius_norm(A)

      !checked: result invariant with value of n
      n = 1
      DO
         IF ((norm_scalar/2.0_dp**n) < 1.0_dp) EXIT
         n = n + 1
      END DO

      !Calculate PR(n) recursively
      CALL dbcsr_create(PRn, template=A, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(work, template=A, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_multiply('N', 'N', 1.0_dp, P, R, 0.0_dp, work, filter_eps=filter_eps)
      CALL dbcsr_create(exp_A2n, template=A, matrix_type=dbcsr_type_no_symmetry)

      DO k = 1, n
         CALL matrix_exponential(exp_A2n, A, 1.0_dp, 0.5_dp**k, threshold)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, exp_A2n, work, 0.0_dp, PRn, filter_eps=filter_eps)
         CALL dbcsr_multiply('N', 'N', 1.0_dp, work, exp_A2n, 1.0_dp, PRn, filter_eps=filter_eps)
         CALL dbcsr_copy(work, PRn)
      END DO
      CALL dbcsr_release(exp_A2n)

      !Calculate Y iteratively, until convergence
      CALL dbcsr_create(A2n, template=A, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_copy(A2n, A)
      CALL dbcsr_scale(A2n, 0.5_dp**n)
      CALL dbcsr_create(Y, template=A, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(Yk, template=A, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(work2, template=A, matrix_type=dbcsr_type_no_symmetry)

      !k=1
      CALL dbcsr_scale(PRn, 0.5_dp**n)
      CALL dbcsr_copy(work, PRn)
      CALL dbcsr_copy(work2, PRn)
      CALL dbcsr_add(Y, PRn, 1.0_dp, 1.0_dp)

      k = 1
      DO
         k = k + 1
         CALL dbcsr_multiply('N', 'N', 1.0_dp/REAL(k, dp), A2n, work, 0.0_dp, Yk, filter_eps=filter_eps)
         CALL dbcsr_multiply('N', 'N', 1.0_dp/REAL(k, dp), work2, A2n, 0.0_dp, PRn, filter_eps=filter_eps)

         CALL dbcsr_add(Yk, PRn, 1.0_dp, 1.0_dp)
         CALL dbcsr_add(Y, Yk, 1.0_dp, 1.0_dp)

         IF (dbcsr_frobenius_norm(Yk) < threshold) EXIT
         CALL dbcsr_copy(work, Yk)
         CALL dbcsr_copy(work2, PRn)
      END DO

      CALL dbcsr_release(work)
      CALL dbcsr_release(work2)
      CALL dbcsr_release(PRn)
      CALL dbcsr_release(A2n)
      CALL dbcsr_release(Yk)

      CALL timestop(handle)

   END SUBROUTINE build_Y_matrix

! **************************************************************************************************
!> \brief Overwrites the "optimal" Laplace quadrature with that of the first step
!> \param tj ...
!> \param wj ...
!> \param tau_tj ...
!> \param tau_wj ...
!> \param weights_cos_tf_t_to_w ...
!> \param weights_cos_tf_w_to_t ...
!> \param do_laplace ...
!> \param do_im_time ...
!> \param num_integ_points ...
!> \param unit_nr ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE keep_initial_quad(tj, wj, tau_tj, tau_wj, weights_cos_tf_t_to_w, weights_cos_tf_w_to_t, &
                                do_laplace, do_im_time, num_integ_points, unit_nr, qs_env)

      REAL(dp), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: tj, wj, tau_tj, tau_wj
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: weights_cos_tf_t_to_w, &
                                                            weights_cos_tf_w_to_t
      LOGICAL, INTENT(IN)                                :: do_laplace, do_im_time
      INTEGER, INTENT(IN)                                :: num_integ_points, unit_nr
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: jquad

      IF (do_laplace .OR. do_im_time) THEN
         IF (.NOT. ASSOCIATED(qs_env%mp2_env%ri_rpa_im_time%tau_tj)) THEN
            ALLOCATE (qs_env%mp2_env%ri_rpa_im_time%tau_tj(0:num_integ_points))
            ALLOCATE (qs_env%mp2_env%ri_rpa_im_time%tau_wj(num_integ_points))
            qs_env%mp2_env%ri_rpa_im_time%tau_tj(:) = tau_tj(:)
            qs_env%mp2_env%ri_rpa_im_time%tau_wj(:) = tau_wj(:)
         ELSE
            !If weights already stored, we overwrite the new ones
            tau_tj(:) = qs_env%mp2_env%ri_rpa_im_time%tau_tj(:)
            tau_wj(:) = qs_env%mp2_env%ri_rpa_im_time%tau_wj(:)
         END IF
      END IF
      IF (.NOT. do_laplace) THEN
         IF (.NOT. ASSOCIATED(qs_env%mp2_env%ri_rpa_im_time%tj)) THEN
            ALLOCATE (qs_env%mp2_env%ri_rpa_im_time%tj(num_integ_points))
            ALLOCATE (qs_env%mp2_env%ri_rpa_im_time%wj(num_integ_points))
            qs_env%mp2_env%ri_rpa_im_time%tj(:) = tj(:)
            qs_env%mp2_env%ri_rpa_im_time%wj(:) = wj(:)
            IF (do_im_time) THEN
               ALLOCATE (qs_env%mp2_env%ri_rpa_im_time%weights_cos_tf_t_to_w(num_integ_points, num_integ_points))
               ALLOCATE (qs_env%mp2_env%ri_rpa_im_time%weights_cos_tf_w_to_t(num_integ_points, num_integ_points))
               qs_env%mp2_env%ri_rpa_im_time%weights_cos_tf_t_to_w(:, :) = weights_cos_tf_t_to_w(:, :)
               qs_env%mp2_env%ri_rpa_im_time%weights_cos_tf_w_to_t(:, :) = weights_cos_tf_w_to_t(:, :)
            END IF
         ELSE
            tj(:) = qs_env%mp2_env%ri_rpa_im_time%tj(:)
            wj(:) = qs_env%mp2_env%ri_rpa_im_time%wj(:)
            IF (do_im_time) THEN
               weights_cos_tf_t_to_w(:, :) = qs_env%mp2_env%ri_rpa_im_time%weights_cos_tf_t_to_w(:, :)
               weights_cos_tf_w_to_t(:, :) = qs_env%mp2_env%ri_rpa_im_time%weights_cos_tf_w_to_t(:, :)
            END IF
         END IF
      END IF
      IF (unit_nr > 0) THEN
         !Printing order same as in mp2_grids.F for consistency
         IF (ASSOCIATED(qs_env%mp2_env%ri_rpa_im_time%tj) .AND. (.NOT. do_laplace)) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MINIMAX_INFO| Number of integration points:", num_integ_points
            WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") &
               "MINIMAX_INFO| Minimax params (freq grid, scaled):", "Weights", "Abscissas"
            DO jquad = 1, num_integ_points
               WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") wj(jquad), tj(jquad)
            END DO
            CALL m_flush(unit_nr)
         END IF
         IF (ASSOCIATED(qs_env%mp2_env%ri_rpa_im_time%tau_tj)) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
               "MINIMAX_INFO| Number of integration points:", num_integ_points
            WRITE (UNIT=unit_nr, FMT="(T3,A,T54,A,T72,A)") &
               "MINIMAX_INFO| Minimax params (time grid, scaled):", "Weights", "Abscissas"
            DO jquad = 1, num_integ_points
               WRITE (UNIT=unit_nr, FMT="(T41,F20.10,F20.10)") tau_wj(jquad), tau_tj(jquad)
            END DO
            CALL m_flush(unit_nr)
         END IF
      END IF

   END SUBROUTINE keep_initial_quad

END MODULE rpa_im_time_force_methods
