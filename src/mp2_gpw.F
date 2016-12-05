!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Rountines to calculate MP2 energy using pw
!> \par History
!>      10.2011 created [Joost VandeVondele and Mauro Del Ben]
! **************************************************************************************************
MODULE mp2_gpw
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: BLACS_GRID_COL,&
                                              BLACS_GRID_SQUARE,&
                                              cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type,&
                                              get_blacs_info
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_add_on_diag, cp_dbcsr_allocate_matrix_set, cp_dbcsr_clear_mempools, &
        cp_dbcsr_copy, cp_dbcsr_create, cp_dbcsr_deallocate_matrix_set, cp_dbcsr_distribution_new, &
        cp_dbcsr_distribution_release, cp_dbcsr_filter, cp_dbcsr_get_block_p, cp_dbcsr_get_info, &
        cp_dbcsr_init, cp_dbcsr_init_p, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
        cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, &
        cp_dbcsr_mp_new, cp_dbcsr_mp_release, cp_dbcsr_multiply, cp_dbcsr_p_type, &
        cp_dbcsr_release, cp_dbcsr_reserve_all_blocks, cp_dbcsr_reserve_diag_blocks, cp_dbcsr_set, &
        cp_dbcsr_type, dbcsr_distribution_obj, dbcsr_mp_obj, dbcsr_type_no_symmetry, &
        dbcsr_type_real_default, dbcsr_type_symmetric
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_dist2d_to_dist,&
                                              cp_dbcsr_m_by_n_from_row_template,&
                                              cp_dbcsr_m_by_n_from_template
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_get_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_release,&
                                              cp_logger_set,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                              distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE distribution_methods,            ONLY: distribute_molecules_1d,&
                                              distribute_molecules_2d
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: default_output_unit,&
                                              m_flush,&
                                              m_memory
   USE mao_basis,                       ONLY: mao_generate_basis
   USE message_passing,                 ONLY: mp_comm_split_direct,&
                                              mp_max,&
                                              mp_min,&
                                              mp_sendrecv,&
                                              mp_sum
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types_new,              ONLY: molecule_type
   USE mp2_cphf,                        ONLY: solve_z_vector_eq
   USE mp2_ri_gpw,                      ONLY: mp2_ri_gpw_compute_en,&
                                              mp2_ri_gpw_compute_in
   USE mp2_ri_grad,                     ONLY: calc_ri_mp2_nonsep
   USE mp2_types,                       ONLY: mp2_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_methods,                  ONLY: pw_env_create,&
                                              pw_env_rebuild
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_release,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_scale,&
                                              pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_release
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_interactions,                 ONLY: init_interaction_radii
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE rpa_ri_gpw,                      ONLY: rpa_ri_compute_en
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list,&
                                              task_list_type
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_gpw'

   PUBLIC :: mp2_gpw_main

CONTAINS

! **************************************************************************************************
!> \brief with a big bang to mp2
!> \param qs_env ...
!> \param mp2_env ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_EX ...
!> \param Emp2_S ...
!> \param Emp2_T ...
!> \param mos_mp2 ...
!> \param para_env ...
!> \param unit_nr ...
!> \param calc_forces ...
!> \param calc_ex ...
!> \param do_ri_mp2 ...
!> \param do_ri_rpa ...
!> \param do_ri_sos_laplace_mp2 ...
!> \author Mauro Del Ben and Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE mp2_gpw_main(qs_env, mp2_env, Emp2, Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T, &
                           mos_mp2, para_env, unit_nr, calc_forces, calc_ex, do_ri_mp2, do_ri_rpa, &
                           do_ri_sos_laplace_mp2)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      REAL(KIND=dp)                                      :: Emp2, Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos_mp2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: unit_nr
      LOGICAL, INTENT(IN)                                :: calc_forces
      LOGICAL                                            :: calc_ex
      LOGICAL, OPTIONAL                                  :: do_ri_mp2, do_ri_rpa, &
                                                            do_ri_sos_laplace_mp2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_gpw_main', routineP = moduleN//':'//routineN

      INTEGER :: blacs_grid_layout, color_sub, color_sub_3c, comm_sub, comm_sub_3c, dimen, &
         dimen_RI, group_size_3c, gw_corr_lev_occ, gw_corr_lev_occ_beta, gw_corr_lev_virt, &
         gw_corr_lev_virt_beta, handle, homo, homo_beta, i, i_multigrid, local_unit_nr, &
         my_group_L_end, my_group_L_size, my_group_L_start, n_multigrid, nelectron, &
         nelectron_beta, nmo, nspins, ri_metric, ri_metric_gw
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ends_array, ends_B_all, ends_B_virtual, &
         ends_B_virtual_beta, sizes_array, sizes_B_all, sizes_B_virtual, sizes_B_virtual_beta, &
         starts_array, starts_B_all, starts_B_virtual, starts_B_virtual_beta
      LOGICAL :: blacs_repeatable, do_im_time, do_mao, my_do_gw, my_do_ri_mp2, my_do_ri_rpa, &
         my_do_ri_sos_laplace_mp2, skip_load_balance_distributed
      REAL(KIND=dp) :: cutoff_old, Emp2_AB, Emp2_BB, Emp2_Cou_BB, Emp2_d2_AB, Emp2_d_AB, &
         Emp2_EX_BB, eps_gvg_rspace_old, eps_pgf_orb_old, eps_rho_rspace_old, progression_factor, &
         relative_cutoff_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: e_cutoff_old, Eigenval, Eigenval_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C, BIb_C_beta, BIb_C_gw, &
                                                            BIb_C_gw_beta
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues, mo_eigenvalues_beta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub, blacs_env_sub_im_time_3c, &
                                                            blacs_env_sub_mat_munu
      TYPE(cp_dbcsr_p_type) :: mat_dm_occ_global_mao, mat_dm_occ_local, mat_dm_virt_global_mao, &
         mat_dm_virt_local, mat_M, mat_munu, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, &
         mat_P_global, mat_P_local
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: mao_coeff_occ, mao_coeff_occ_A, &
         mao_coeff_virt, mao_coeff_virt_A, mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
         mat_3c_overl_int_mao_for_virt, matrix_s
      TYPE(cp_dbcsr_type), POINTER :: mo_coeff_all, mo_coeff_all_beta, mo_coeff_gw, &
         mo_coeff_gw_beta, mo_coeff_o, mo_coeff_o_beta, mo_coeff_v, mo_coeff_v_beta
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L_RI_metric, mo_coeff, &
                                                            mo_coeff_beta
      TYPE(cp_logger_type), POINTER                      :: logger, logger_sub
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub, para_env_sub_im_time_3c, &
                                                            para_env_sub_im_time_P, &
                                                            para_env_sub_RPA
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_all, sab_orb_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_p_type)                                    :: pot_g, rho_g, rho_r
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(task_list_type), POINTER                      :: task_list_sub

      CALL timeset(routineN, handle)

      ! check if we want to do ri-mp2
      my_do_ri_mp2 = .FALSE.
      IF (PRESENT(do_ri_mp2)) my_do_ri_mp2 = do_ri_mp2

      ! check if we want to do ri-rpa
      my_do_ri_rpa = .FALSE.
      IF (PRESENT(do_ri_rpa)) my_do_ri_rpa = do_ri_rpa

      ! check if we want to do ri-sos-laplace-mp2
      my_do_ri_sos_laplace_mp2 = .FALSE.
      IF (PRESENT(do_ri_sos_laplace_mp2)) my_do_ri_sos_laplace_mp2 = do_ri_sos_laplace_mp2

      ! check if we want to do imaginary time
      do_im_time = mp2_env%ri_rpa%do_im_time
      do_mao = mp2_env%ri_rpa_im_time%do_mao
      ! Get the number of spins
      nspins = SIZE(mos_mp2)

      ! ... setup needed to be able to qs_integrate in a subgroup.
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, mos=mos)
      CALL get_mo_set(mo_set=mos_mp2(1)%mo_set, nelectron=nelectron, &
                      eigenvalues=mo_eigenvalues, nmo=nmo, homo=homo, &
                      mo_coeff=mo_coeff, nao=dimen)
      IF (nspins == 2) THEN
         CALL get_mo_set(mo_set=mos_mp2(2)%mo_set, nelectron=nelectron_beta, &
                         eigenvalues=mo_eigenvalues_beta, homo=homo_beta, &
                         mo_coeff=mo_coeff_beta)
      ENDIF

      CALL get_mo_set(mo_set=mos(1)%mo_set)

      ! a para_env
      color_sub = para_env%mepos/mp2_env%mp2_num_proc
      CALL mp_comm_split_direct(para_env%group, comm_sub, color_sub)
      NULLIFY (para_env_sub)
      CALL cp_para_env_create(para_env_sub, comm_sub)

      ! each of the sub groups might need to generate output
      logger => cp_get_default_logger()
      IF (para_env%mepos == para_env%source) THEN
         local_unit_nr = cp_logger_get_default_unit_nr(logger, local=.FALSE.)
      ELSE
         local_unit_nr = default_output_unit
      ENDIF

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T71,F10.1)") &
            "GPW_INFO| Density cutoff [a.u.]:", mp2_env%mp2_gpw%cutoff*0.5_dp
         WRITE (UNIT=unit_nr, FMT="(T3,A,T71,F10.1)") &
            "GPW_INFO| Relative density cutoff [a.u.]:", mp2_env%mp2_gpw%relative_cutoff*0.5_dp
         CALL m_flush(unit_nr)
      ENDIF

      ! a logger
      NULLIFY (logger_sub)
      CALL cp_logger_create(logger_sub, para_env=para_env_sub, &
                            default_global_unit_nr=local_unit_nr, close_global_unit_on_dealloc=.FALSE.)
      CALL cp_logger_set(logger_sub, local_filename="MP2_localLog")
      ! set to a custom print level (we could also have a different print level for para_env%source)
      logger_sub%iter_info%print_level = mp2_env%mp2_gpw%print_level
      CALL cp_add_default_logger(logger_sub)

      ! a blacs_env (ignore the globenv stored defaults for now)
      blacs_grid_layout = BLACS_GRID_SQUARE
      blacs_repeatable = .TRUE.
      NULLIFY (blacs_env_sub)
      CALL cp_blacs_env_create(blacs_env_sub, para_env_sub, &
                               blacs_grid_layout, &
                               blacs_repeatable)

      blacs_env_sub_mat_munu => blacs_env_sub

      IF (do_im_time) THEN

         group_size_3c = mp2_env%ri_rpa_im_time%group_size_3c

         IF (group_size_3c > para_env%num_pe) THEN

            group_size_3c = para_env%num_pe
            mp2_env%ri_rpa_im_time%group_size_3c = para_env%num_pe

         END IF

         ! only allow group_size_3c which is a factor of the total number of MPI tasks
         CPASSERT(MODULO(para_env%num_pe, group_size_3c) == 0)

         ! a para_env
         color_sub_3c = para_env%mepos/mp2_env%ri_rpa_im_time%group_size_3c
         CALL mp_comm_split_direct(para_env%group, comm_sub_3c, color_sub_3c)
         NULLIFY (para_env_sub_im_time_3c)
         CALL cp_para_env_create(para_env_sub_im_time_3c, comm_sub_3c)

         blacs_grid_layout = BLACS_GRID_SQUARE
         blacs_repeatable = .TRUE.
         NULLIFY (blacs_env_sub_im_time_3c)
         CALL cp_blacs_env_create(blacs_env_sub_im_time_3c, para_env_sub_im_time_3c, &
                                  blacs_grid_layout, &
                                  blacs_repeatable)

         blacs_env_sub_mat_munu => blacs_env_sub_im_time_3c

      END IF

      ! get stuff
      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      molecule_set=molecule_set, &
                      molecule_kind_set=molecule_kind_set, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s)

      CALL get_eps_old(dft_control, eps_pgf_orb_old, eps_rho_rspace_old, eps_gvg_rspace_old)

      CALL create_mat_munu(mat_munu, qs_env, mp2_env, para_env, dft_control, atomic_kind_set, qs_kind_set, &
                           atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub_mat_munu, &
                           sab_orb_all=sab_orb_all, do_im_time=do_im_time)

      IF (do_im_time) THEN

         CALL create_mao_basis_and_matrices(mat_dm_occ_global_mao, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, &
                                            mat_dm_virt_global_mao, mat_munu, do_mao, qs_env, mp2_env, &
                                            mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, &
                                            matrix_s, mo_coeff, mo_coeff_beta, homo, homo_beta, nmo, nspins, unit_nr, &
                                            mo_eigenvalues, mo_eigenvalues_beta)

      END IF

      ! which RI metric we want to have
      ri_metric = mp2_env%ri_metric

      ! check if we want to do ri-g0w0 on top of ri-rpa
      my_do_gw = mp2_env%ri_rpa%do_ri_g0w0
      gw_corr_lev_occ = mp2_env%ri_g0w0%corr_mos_occ
      gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_virt

      ! metric for GW matrix elements B^nm_P
      ri_metric_gw = mp2_env%ri_g0w0%ri_metric

      ! and the array of mos
      ALLOCATE (Eigenval(dimen))
      Eigenval(:) = mo_eigenvalues(:)
      IF (nspins == 2) THEN
         ALLOCATE (Eigenval_beta(dimen))
         Eigenval_beta(:) = mo_eigenvalues_beta(:)
      ENDIF

      ! for imag. time, we do not need this
      IF (.NOT. do_im_time) THEN

         ! new routine: replicate a full matrix from one para_env to a smaller one
         ! keeping the memory usage as small as possible in this case the
         ! output the two part of the C matrix (virtual, occupied)
         CALL replicate_mat_to_subgroup(mp2_env, para_env, para_env_sub, mo_coeff, dimen, homo, mat_munu, &
                                        mo_coeff_o, mo_coeff_v, mo_coeff_all, mo_coeff_gw, my_do_gw, &
                                        gw_corr_lev_occ, gw_corr_lev_virt)

         ! if open shell case replicate also the coefficient matrix for the beta orbitals
         IF (nspins == 2) THEN

            gw_corr_lev_occ_beta = mp2_env%ri_g0w0%corr_mos_occ_beta
            gw_corr_lev_virt_beta = mp2_env%ri_g0w0%corr_mos_virt_beta

            CALL replicate_mat_to_subgroup(mp2_env, para_env, para_env_sub, mo_coeff_beta, dimen, homo_beta, mat_munu, &
                                           mo_coeff_o_beta, mo_coeff_v_beta, mo_coeff_all_beta, mo_coeff_gw_beta, &
                                           my_do_gw, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta)
         END IF

         ! hack hack hack XXXXXXXXXXXXXXX rebuilds the pw_en with the new cutoffs
         progression_factor = dft_control%qs_control%progression_factor
         n_multigrid = SIZE(dft_control%qs_control%e_cutoff)
         ALLOCATE (e_cutoff_old(n_multigrid))
         e_cutoff_old(:) = dft_control%qs_control%e_cutoff
         cutoff_old = dft_control%qs_control%cutoff

         dft_control%qs_control%cutoff = mp2_env%mp2_gpw%cutoff*0.5_dp
         dft_control%qs_control%e_cutoff(1) = dft_control%qs_control%cutoff
         DO i_multigrid = 2, n_multigrid
            dft_control%qs_control%e_cutoff(i_multigrid) = dft_control%qs_control%e_cutoff(i_multigrid-1) &
                                                           /progression_factor
         END DO

         relative_cutoff_old = dft_control%qs_control%relative_cutoff
         dft_control%qs_control%relative_cutoff = mp2_env%mp2_gpw%relative_cutoff*0.5_dp

         ! a pw_env
         NULLIFY (pw_env_sub)
         CALL pw_env_create(pw_env_sub)
         CALL pw_env_rebuild(pw_env_sub, qs_env, para_env_sub)

         CALL pw_env_get(pw_env_sub, auxbas_pw_pool=auxbas_pw_pool, &
                         poisson_env=poisson_env)
         ! hack hack hack XXXXXXXXXXXXXXX

         ! now we need a task list, hard code skip_load_balance_distributed
         NULLIFY (task_list_sub)
         skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
         CALL allocate_task_list(task_list_sub)
         CALL generate_qs_task_list(ks_env, task_list_sub, &
                                    reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                                    skip_load_balance_distributed=skip_load_balance_distributed, &
                                    pw_env_external=pw_env_sub, sab_orb_external=sab_orb_sub)

         ! get some of the grids ready
         NULLIFY (rho_r%pw, rho_g%pw, pot_g%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, pot_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)

         ! run the FFT once, to set up buffers and to take into account the memory
         rho_r%pw%cr3d = 0.0D0
         CALL pw_transfer(rho_r%pw, rho_g%pw)

      END IF

      ! now we're kind of ready to go....
      Emp2_S = 0.0_dp
      Emp2_T = 0.0_dp
      IF (my_do_ri_mp2 .OR. my_do_ri_rpa .OR. my_do_ri_sos_laplace_mp2) THEN
         ! RI-GPW integrals (same stuff for both RPA and MP2)
         IF (nspins == 2) THEN
            ! open shell case (RI) here the (ia|K) integrals are computed for both the alpha and beta components
            CALL mp2_ri_gpw_compute_in( &
               BIb_C, BIb_C_gw, ends_array, ends_B_virtual, sizes_array, &
               sizes_B_virtual, starts_array, starts_B_virtual, &
               dimen_RI, qs_env, para_env, para_env_sub, color_sub, dft_control, cell, particle_set, &
               atomic_kind_set, qs_kind_set, mo_coeff, fm_matrix_L_RI_metric, nmo, homo, rho_r, rho_g, pot_g, &
               mat_munu, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, sab_orb_sub, sab_orb_all, &
               pw_env_sub, poisson_env, auxbas_pw_pool, &
               task_list_sub, mo_coeff_o, mo_coeff_v, mo_coeff_all, mo_coeff_gw, &
               mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, mp2_env%calc_PQ_cond_num, calc_forces, blacs_env_sub, my_do_gw, &
               starts_B_all, sizes_B_all, ends_B_all, gw_corr_lev_occ, gw_corr_lev_virt, &
               do_im_time, do_mao, mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, &
               mao_coeff_occ, mao_coeff_virt, ri_metric, ri_metric_gw, &
               BIb_C_beta, BIb_C_gw_beta, ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta, &
               homo_beta, mo_coeff_o_beta, mo_coeff_v_beta, mo_coeff_all_beta, mo_coeff_gw_beta)
         ELSE
            ! closed shell case (RI)
            CALL mp2_ri_gpw_compute_in(BIb_C, BIb_C_gw, ends_array, ends_B_virtual, sizes_array, &
                                       sizes_B_virtual, starts_array, starts_B_virtual, &
                                       dimen_RI, qs_env, para_env, para_env_sub, color_sub, dft_control, cell, particle_set, &
                                       atomic_kind_set, qs_kind_set, mo_coeff, fm_matrix_L_RI_metric, &
                                       nmo, homo, rho_r, rho_g, pot_g, &
                                       mat_munu, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, sab_orb_sub, sab_orb_all, &
                                       pw_env_sub, poisson_env, auxbas_pw_pool, &
                                       task_list_sub, mo_coeff_o, mo_coeff_v, mo_coeff_all, mo_coeff_gw, &
                                       mp2_env%mp2_gpw%eps_filter, unit_nr, &
                                       mp2_env%mp2_memory, mp2_env%calc_PQ_cond_num, calc_forces, blacs_env_sub, my_do_gw, &
                                       starts_B_all, sizes_B_all, ends_B_all, gw_corr_lev_occ, gw_corr_lev_virt, &
                                       do_im_time, do_mao, mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, &
                                       mat_3c_overl_int_mao_for_virt, mao_coeff_occ, mao_coeff_virt, ri_metric, ri_metric_gw)

         END IF
      ELSE
         ! Canonical MP2-GPW
         IF (nspins == 2) THEN
            ! alpha-alpha and alpha-beta components
            IF (unit_nr > 0) WRITE (unit_nr, *)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Alpha (ia|'
            CALL mp2_gpw_compute( &
               Emp2, Emp2_Cou, Emp2_EX, qs_env, para_env, para_env_sub, color_sub, dft_control, cell, particle_set, &
               atomic_kind_set, qs_kind_set, mo_coeff, Eigenval, nmo, homo, rho_r, rho_g, pot_g, &
               mat_munu, pw_env_sub, poisson_env, auxbas_pw_pool, task_list_sub, &
               mo_coeff_o, mo_coeff_v, mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, calc_ex, blacs_env_sub, &
               homo_beta, mo_coeff_o_beta, mo_coeff_v_beta, Eigenval_beta, Emp2_AB)

            ! beta-beta component
            IF (unit_nr > 0) WRITE (unit_nr, *)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Beta (ia|'
            CALL mp2_gpw_compute( &
               Emp2_BB, Emp2_Cou_BB, Emp2_EX_BB, qs_env, para_env, para_env_sub, color_sub, dft_control, cell, particle_set, &
               atomic_kind_set, qs_kind_set, mo_coeff_beta, Eigenval_beta, nmo, homo_beta, rho_r, rho_g, pot_g, &
               mat_munu, pw_env_sub, poisson_env, auxbas_pw_pool, task_list_sub, &
               mo_coeff_o_beta, mo_coeff_v_beta, mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, calc_ex, blacs_env_sub)

         ELSE
            ! closed shell case
            CALL mp2_gpw_compute( &
               Emp2, Emp2_Cou, Emp2_EX, qs_env, para_env, para_env_sub, color_sub, dft_control, cell, particle_set, &
               atomic_kind_set, qs_kind_set, mo_coeff, Eigenval, nmo, homo, rho_r, rho_g, pot_g, &
               mat_munu, pw_env_sub, poisson_env, auxbas_pw_pool, task_list_sub, &
               mo_coeff_o, mo_coeff_v, mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, calc_ex, blacs_env_sub)
         END IF
      END IF

      IF (.NOT. do_im_time) THEN

         ! and now free the whole lot
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_g%pw)

      END IF

      ! Free possibly large buffers allocated by dbcsr on the GPU,
      ! large hybrid dgemm/pdgemm's coming later will need the space.
      CALL cp_dbcsr_clear_mempools()

      ! moved down
      ! CALL deallocate_task_list(task_list_sub)
      ! CALL pw_env_release(pw_env_sub)

      IF (calc_forces) THEN
         ! make a copy of mo_coeff_o and mo_coeff_v
         NULLIFY (mp2_env%ri_grad%mo_coeff_o)
         CALL cp_dbcsr_init_p(mp2_env%ri_grad%mo_coeff_o)
         CALL cp_dbcsr_copy(mp2_env%ri_grad%mo_coeff_o, mo_coeff_o, name="mo_coeff_o")
         NULLIFY (mp2_env%ri_grad%mo_coeff_v)
         CALL cp_dbcsr_init_p(mp2_env%ri_grad%mo_coeff_v)
         CALL cp_dbcsr_copy(mp2_env%ri_grad%mo_coeff_v, mo_coeff_v, name="mo_coeff_v")
         IF (nspins == 2) THEN
            NULLIFY (mp2_env%ri_grad%mo_coeff_o_beta)
            CALL cp_dbcsr_init_p(mp2_env%ri_grad%mo_coeff_o_beta)
            CALL cp_dbcsr_copy(mp2_env%ri_grad%mo_coeff_o_beta, mo_coeff_o_beta, name="mo_coeff_o_b")
            NULLIFY (mp2_env%ri_grad%mo_coeff_v_beta)
            CALL cp_dbcsr_init_p(mp2_env%ri_grad%mo_coeff_v_beta)
            CALL cp_dbcsr_copy(mp2_env%ri_grad%mo_coeff_v_beta, mo_coeff_v_beta, name="mo_coeff_v_b")
         ENDIF
         my_group_L_size = sizes_array(color_sub)
         my_group_L_start = starts_array(color_sub)
         my_group_L_end = ends_array(color_sub)
      END IF
      ! Copy mo coeffs for RPA AXK
      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         NULLIFY (mp2_env%ri_rpa%mo_coeff_o)
         CALL cp_dbcsr_init_p(mp2_env%ri_rpa%mo_coeff_o)
         CALL cp_dbcsr_copy(mp2_env%ri_rpa%mo_coeff_o, mo_coeff_o, name="mo_coeff_o")
         NULLIFY (mp2_env%ri_rpa%mo_coeff_v)
         CALL cp_dbcsr_init_p(mp2_env%ri_rpa%mo_coeff_v)
         CALL cp_dbcsr_copy(mp2_env%ri_rpa%mo_coeff_v, mo_coeff_v, name="mo_coeff_v")
      ENDIF

      IF (.NOT. do_im_time) THEN

         CALL cp_dbcsr_release(mo_coeff_o)
         DEALLOCATE (mo_coeff_o)
         CALL cp_dbcsr_release(mo_coeff_v)
         DEALLOCATE (mo_coeff_v)
         IF (my_do_gw) THEN
            CALL cp_dbcsr_release(mo_coeff_all)
            DEALLOCATE (mo_coeff_all)
         END IF

         IF (nspins == 2) THEN
            CALL cp_dbcsr_release(mo_coeff_o_beta)
            DEALLOCATE (mo_coeff_o_beta)
            CALL cp_dbcsr_release(mo_coeff_v_beta)
            DEALLOCATE (mo_coeff_v_beta)
            IF (my_do_gw) THEN
               CALL cp_dbcsr_release(mo_coeff_all_beta)
               DEALLOCATE (mo_coeff_all_beta)
            END IF
         END IF

      END IF

      IF (.NOT. calc_forces) THEN
         IF (.NOT. mp2_env%ri_rpa%do_ri_axk) THEN

            CALL cp_dbcsr_release(mat_munu%matrix)
            DEALLOCATE (mat_munu%matrix)

            DO i = 1, SIZE(sab_orb_sub)
               CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
            END DO
            DEALLOCATE (sab_orb_sub)

         ENDIF

      END IF

      ! decide if to do RI-RPA or RI-MP2
      IF (my_do_ri_rpa .OR. my_do_ri_sos_laplace_mp2) THEN

         IF (do_im_time) THEN

            CALL create_dbcsr_matrices_im_time(mat_munu, mat_P_local, mat_P_global, mat_M, mat_dm_occ_global_mao, &
                                               mat_dm_virt_global_mao, mat_dm_occ_local, mat_dm_virt_local, &
                                               do_mao, qs_env, mp2_env, para_env, dft_control, atomic_kind_set, qs_kind_set, &
                                               atom2d, molecule_kind_set, molecule_set, particle_set, cell, &
                                               para_env_sub_im_time_P, blacs_env_sub_im_time_3c, sab_orb_all)

            para_env_sub_RPA => para_env_sub_im_time_3c

         ELSE

            para_env_sub_RPA => para_env_sub

         END IF

         ! RI-RPA
         IF (nspins == 2) THEN
            CALL rpa_ri_compute_en(qs_env, Emp2, mp2_env, BIb_C, BIb_C_gw, para_env, para_env_sub_RPA, color_sub, &
                                   ends_array, ends_B_virtual, ends_B_all, sizes_array, sizes_B_virtual, sizes_B_all, &
                                   starts_array, starts_B_virtual, starts_B_all, mo_coeff, fm_matrix_L_RI_metric, &
                                   Eigenval, nmo, homo, dimen_RI, gw_corr_lev_occ, gw_corr_lev_virt, &
                                   unit_nr, my_do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_mao, matrix_s, &
                                   mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, &
                                   mat_munu, mat_dm_occ_local, mat_dm_virt_local, &
                                   mat_P_local, mat_P_global, mat_M, &
                                   mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, &
                                   mp2_env%mp2_gpw%eps_filter, BIb_C_beta, homo_beta, Eigenval_beta, &
                                   ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta, &
                                   mo_coeff_beta, BIb_C_gw_beta, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta)
         ELSE
            CALL rpa_ri_compute_en(qs_env, Emp2, mp2_env, BIb_C, BIb_C_gw, para_env, para_env_sub_RPA, color_sub, &
                                   ends_array, ends_B_virtual, ends_B_all, sizes_array, sizes_B_virtual, sizes_B_all, &
                                   starts_array, starts_B_virtual, starts_B_all, mo_coeff, fm_matrix_L_RI_metric, &
                                   Eigenval, nmo, homo, dimen_RI, gw_corr_lev_occ, gw_corr_lev_virt, &
                                   unit_nr, my_do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_mao, matrix_s, &
                                   mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, &
                                   mat_munu, mat_dm_occ_local, mat_dm_virt_local, &
                                   mat_P_local, mat_P_global, mat_M, &
                                   mat_3c_overl_int, mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, &
                                   mp2_env%mp2_gpw%eps_filter)
         END IF

         IF (do_im_time) THEN
            CALL clean_up_im_time(mat_munu, mat_P_local, mat_P_global, mat_M, mat_dm_occ_global_mao, &
                                  mat_dm_virt_global_mao, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, &
                                  mat_dm_occ_local, mat_dm_virt_local, fm_matrix_L_RI_metric, para_env_sub_im_time_3c, &
                                  para_env_sub_im_time_P, mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, mp2_env)
         END IF

         ! Release some memory for AXK
         IF (mp2_env%ri_rpa%do_ri_axk) THEN

            CALL cp_dbcsr_release(mat_munu%matrix)
            DEALLOCATE (mat_munu%matrix)

            DO i = 1, SIZE(sab_orb_sub)
               CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
            END DO
            DEALLOCATE (sab_orb_sub)

         END IF

      ELSE
         IF (my_do_ri_mp2) THEN
            ! RI-MP2-GPW compute energy
            IF (nspins == 2) THEN
               ! alpha-alpha component
               CALL mp2_ri_gpw_compute_en( &
                  Emp2, Emp2_Cou, Emp2_EX, BIb_C, mp2_env, para_env, para_env_sub, color_sub, &
                  ends_array, ends_B_virtual, sizes_array, sizes_B_virtual, starts_array, starts_B_virtual, &
                  Eigenval, nmo, homo, dimen_RI, unit_nr, calc_forces, calc_ex, &
                  open_shell_SS=.TRUE.)

               ! beta-beta component
               CALL mp2_ri_gpw_compute_en( &
                  Emp2_BB, Emp2_Cou_BB, Emp2_EX_BB, BIb_C_beta, mp2_env, para_env, para_env_sub, color_sub, &
                  ends_array, ends_B_virtual_beta, sizes_array, &
                  sizes_B_virtual_beta, starts_array, starts_B_virtual_beta, &
                  Eigenval_beta, nmo, homo_beta, dimen_RI, unit_nr, calc_forces, calc_ex, &
                  open_shell_SS=.TRUE.)

               ! alpha-beta case
               CALL mp2_ri_gpw_compute_en( &
                  Emp2_d_AB, Emp2_AB, Emp2_d2_AB, BIb_C, mp2_env, para_env, para_env_sub, color_sub, &
                  ends_array, ends_B_virtual, sizes_array, sizes_B_virtual, starts_array, starts_B_virtual, &
                  Eigenval, nmo, homo, dimen_RI, unit_nr, calc_forces, .FALSE., &
                  .FALSE., BIb_C_beta, homo_beta, Eigenval_beta, &
                  ends_B_virtual_beta, sizes_B_virtual_beta, starts_B_virtual_beta)

            ELSE
               ! closed shell case
               CALL mp2_ri_gpw_compute_en( &
                  Emp2, Emp2_Cou, Emp2_EX, BIb_C, mp2_env, para_env, para_env_sub, color_sub, &
                  ends_array, ends_B_virtual, sizes_array, sizes_B_virtual, starts_array, starts_B_virtual, &
                  Eigenval, nmo, homo, dimen_RI, unit_nr, calc_forces, calc_ex)
            END IF
            ! if we need forces time to calculate the MP2 non-separable contribution
            ! and start coputing the Largrangian
            IF (calc_forces) THEN
               ! since we have to compute again integrals reinitialize the stuff we need
               ! get some of the grids ready
               NULLIFY (rho_r%pw, rho_g%pw, pot_g%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, rho_r%pw, &
                                      use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_pool_create_pw(auxbas_pw_pool, rho_g%pw, &
                                      use_data=COMPLEXDATA1D, &
                                      in_space=RECIPROCALSPACE)
               CALL pw_pool_create_pw(auxbas_pw_pool, pot_g%pw, &
                                      use_data=COMPLEXDATA1D, &
                                      in_space=RECIPROCALSPACE)

               ! the mu_nu matrix (again)
               ! XXXXXXXXXXXXXXXXXXXXXXXXX
               ! ! build a dbcsr matrix the hard way
               ! CALL get_particle_set(particle_set=particle_set,nsgf=rbs)
               ! CALL array_nullify (row_blk_sizes)
               ! CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
               ! ALLOCATE(mat_munu%matrix, STAT=stat)
               ! CPPostcondition(stat==0,cp_failure_level,routineP,failure)
               ! CALL cp_dbcsr_init(mat_munu%matrix)
               ! CALL cp_dbcsr_create(matrix=mat_munu%matrix,&
               !            name="(ai|munu)",&
               !            dist=dbcsr_dist_sub, &
               !            matrix_type=dbcsr_type_symmetric,&
               !            ! matrix_type=dbcsr_type_no_symmetry,&
               !            row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
               !            nze=0)
               ! CALL cp_dbcsr_alloc_block_from_nbl(mat_munu%matrix,sab_orb_sub)
               ! CALL array_release (row_blk_sizes)

               IF (nspins == 2) THEN ! Open shell
                  CALL calc_ri_mp2_nonsep(qs_env, mp2_env, para_env, para_env_sub, dft_control, cell, &
                                          particle_set, atomic_kind_set, qs_kind_set, mo_coeff, nmo, homo, dimen_RI, Eigenval, &
                                          my_group_L_start, my_group_L_end, my_group_L_size, rho_r, rho_g, pot_g, &
                                          mat_munu, pw_env_sub, poisson_env, auxbas_pw_pool, task_list_sub, &
                                          blacs_env_sub, Eigenval_beta, homo_beta, mo_coeff_beta)
               ELSE ! Closed shell
                  CALL calc_ri_mp2_nonsep(qs_env, mp2_env, para_env, para_env_sub, dft_control, cell, &
                                          particle_set, atomic_kind_set, qs_kind_set, mo_coeff, nmo, homo, dimen_RI, Eigenval, &
                                          my_group_L_start, my_group_L_end, my_group_L_size, rho_r, rho_g, pot_g, &
                                          mat_munu, pw_env_sub, poisson_env, auxbas_pw_pool, task_list_sub, &
                                          blacs_env_sub)
               ENDIF

               ! release
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_r%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_g%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, pot_g%pw)

               CALL cp_dbcsr_release(mp2_env%ri_grad%mo_coeff_o)
               DEALLOCATE (mp2_env%ri_grad%mo_coeff_o)

               CALL cp_dbcsr_release(mp2_env%ri_grad%mo_coeff_v)
               DEALLOCATE (mp2_env%ri_grad%mo_coeff_v)

               IF (nspins == 2) THEN
                  CALL cp_dbcsr_release(mp2_env%ri_grad%mo_coeff_o_beta)
                  DEALLOCATE (mp2_env%ri_grad%mo_coeff_o_beta)

                  CALL cp_dbcsr_release(mp2_env%ri_grad%mo_coeff_v_beta)
                  DEALLOCATE (mp2_env%ri_grad%mo_coeff_v_beta)
               ENDIF

               CALL cp_dbcsr_release(mat_munu%matrix)
               DEALLOCATE (mat_munu%matrix)

               DO i = 1, SIZE(sab_orb_sub)
                  CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
               END DO
               DEALLOCATE (sab_orb_sub)

            END IF
         END IF

         IF (nspins == 2) THEN
            ! make order on the MP2 energy contributions
            Emp2_Cou = Emp2_Cou*0.25_dp
            Emp2_EX = Emp2_EX*0.5_dp

            Emp2_Cou_BB = Emp2_Cou_BB*0.25_dp
            Emp2_EX_BB = Emp2_EX_BB*0.5_dp

            Emp2_S = Emp2_AB
            Emp2_T = Emp2_Cou+Emp2_Cou_BB+Emp2_EX+Emp2_EX_BB

            Emp2_Cou = Emp2_Cou+Emp2_Cou_BB+Emp2_AB
            Emp2_EX = Emp2_EX+Emp2_EX_BB
            Emp2 = Emp2_EX+Emp2_Cou
         END IF

      END IF

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
      ! moved from above
      IF (my_do_gw) THEN
         CALL cp_dbcsr_release(mo_coeff_gw)
         DEALLOCATE (mo_coeff_gw)
         IF (nspins == 2) THEN
            CALL cp_dbcsr_release(mo_coeff_gw_beta)
            DEALLOCATE (mo_coeff_gw_beta)
         END IF
      END IF

      IF (.NOT. do_im_time) THEN
         CALL deallocate_task_list(task_list_sub)
         CALL pw_env_release(pw_env_sub)
      END IF

      ! CALL cp_dbcsr_distribution_release(dbcsr_dist_sub)
      ! DEALLOCATE(dbcsr_dist_sub,STAT=stat)
      ! CPPostcondition(stat==0,cp_failure_level,routineP,failure)

      ! DO i=1,SIZE(sab_orb_sub)
      !    CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
      ! END DO
      ! DEALLOCATE(sab_orb_sub,stat=stat)
      ! CPPostconditionNoFail(stat==0,cp_warning_level,routineP)

      ! CALL distribution_2d_release(distribution_2d_sub)

      ! CALL distribution_1d_release(local_particles_sub)
      ! CALL distribution_1d_release(local_molecules_sub)
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx

      ! re-init the radii to be able to generate pair lists with MP2-appropriate screening
      dft_control%qs_control%eps_pgf_orb = eps_pgf_orb_old
      dft_control%qs_control%eps_rho_rspace = eps_rho_rspace_old
      dft_control%qs_control%eps_gvg_rspace = eps_gvg_rspace_old
      CALL init_interaction_radii(dft_control%qs_control, atomic_kind_set, qs_kind_set)

      IF (.NOT. do_im_time) THEN
         ! restore the initial value of the cutoff
         dft_control%qs_control%e_cutoff = e_cutoff_old
         dft_control%qs_control%cutoff = cutoff_old
         dft_control%qs_control%relative_cutoff = relative_cutoff_old
      END IF

      CALL cp_blacs_env_release(blacs_env_sub)

      CALL cp_rm_default_logger()
      CALL cp_logger_release(logger_sub)

      CALL cp_para_env_release(para_env_sub)

      ! finally solve the z-vector equation if forces are required
      IF (calc_forces) THEN
         IF (nspins == 2) THEN
            CALL solve_z_vector_eq(qs_env, mp2_env, para_env, dft_control, &
                                   atomic_kind_set, mo_coeff, nmo, homo, Eigenval, unit_nr, &
                                   Eigenval_beta, homo_beta, mo_coeff_beta)
         ELSE
            CALL solve_z_vector_eq(qs_env, mp2_env, para_env, dft_control, &
                                   atomic_kind_set, mo_coeff, nmo, homo, Eigenval, unit_nr)
         ENDIF
      END IF

      DEALLOCATE (Eigenval)
      IF (nspins == 2) THEN
         DEALLOCATE (Eigenval_beta)
      END IF

      CALL timestop(handle)

   END SUBROUTINE mp2_gpw_main

! **************************************************************************************************
!> \brief ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_EX ...
!> \param qs_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param color_sub ...
!> \param dft_control ...
!> \param cell ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param mo_coeff ...
!> \param Eigenval ...
!> \param nmo ...
!> \param homo ...
!> \param rho_r ...
!> \param rho_g ...
!> \param pot_g ...
!> \param mat_munu ...
!> \param pw_env_sub ...
!> \param poisson_env ...
!> \param auxbas_pw_pool ...
!> \param task_list_sub ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param eps_filter ...
!> \param unit_nr ...
!> \param mp2_memory ...
!> \param calc_ex ...
!> \param blacs_env_sub ...
!> \param homo_beta ...
!> \param mo_coeff_o_beta ...
!> \param mo_coeff_v_beta ...
!> \param Eigenval_beta ...
!> \param Emp2_AB ...
! **************************************************************************************************
   SUBROUTINE mp2_gpw_compute(Emp2, Emp2_Cou, Emp2_EX, qs_env, para_env, para_env_sub, color_sub, dft_control, &
                              cell, particle_set, atomic_kind_set, qs_kind_set, mo_coeff, Eigenval, nmo, homo, &
                              rho_r, rho_g, pot_g, mat_munu, pw_env_sub, &
                              poisson_env, auxbas_pw_pool, task_list_sub, mo_coeff_o, mo_coeff_v, eps_filter, unit_nr, &
                              mp2_memory, calc_ex, blacs_env_sub, homo_beta, mo_coeff_o_beta, &
                              mo_coeff_v_beta, Eigenval_beta, Emp2_AB)

      REAL(KIND=dp)                                      :: Emp2, Emp2_Cou, Emp2_EX
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      INTEGER                                            :: color_sub
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER                                            :: nmo, homo
      TYPE(pw_p_type)                                    :: rho_r, rho_g, pot_g
      TYPE(cp_dbcsr_p_type)                              :: mat_munu
      TYPE(pw_env_type), POINTER                         :: pw_env_sub
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(task_list_type), POINTER                      :: task_list_sub
      TYPE(cp_dbcsr_type), POINTER                       :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp)                                      :: eps_filter
      INTEGER                                            :: unit_nr
      REAL(KIND=dp)                                      :: mp2_memory
      LOGICAL                                            :: calc_ex
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      INTEGER, OPTIONAL                                  :: homo_beta
      TYPE(cp_dbcsr_type), OPTIONAL, POINTER             :: mo_coeff_o_beta, mo_coeff_v_beta
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: Eigenval_beta
      REAL(KIND=dp), OPTIONAL                            :: Emp2_AB

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_gpw_compute', &
         routineP = moduleN//':'//routineN

      INTEGER :: a, a_group_counter, b, b_global, b_group_counter, blk, col, col_offset, col_size, &
         color_counter, comm_exchange, EX_end, EX_end_send, EX_start, EX_start_send, &
         group_counter, handle, handle2, handle3, i, i_counter, i_group_counter, index_proc_shift, &
         j, max_b_size, max_batch_size_A, max_batch_size_I, max_row_col_local, mepos_in_EX_group, &
         my_A_batch_size, my_A_virtual_end, my_A_virtual_start, my_B_size, my_B_virtual_end, &
         my_B_virtual_start, my_I_batch_size, my_I_occupied_end, my_I_occupied_start, &
         my_q_position, ncol_local, nfullcols_total, nfullrows_total, ngroup, nrow_local, one, p
      INTEGER :: p_best, proc_receive, proc_send, q, q_best, row, row_offset, row_size, size_EX, &
         size_EX_send, size_of_exchange_group, sub_sub_color, virtual, virtual_beta, wfn_calc, &
         wfn_calc_best
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: proc_map, sub_proc_map, vector_B_sizes, &
                                                            vector_batch_A_size_group, &
                                                            vector_batch_I_size_group
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: color_array, exchange_group_sizes, &
                                                            local_col_row_info
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: do_alpha_beta
      REAL(KIND=dp)                                      :: mem_min, mem_real, mem_try, pair_energy, &
                                                            wfn_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: my_Cocc, my_Cvirt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C, BIb_Ex, BIb_send
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_type)                                :: matrix_ia_jb, matrix_ia_jb_beta, &
                                                            matrix_ia_jnu, matrix_ia_jnu_beta
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_jb
      TYPE(cp_para_env_type), POINTER                    :: para_env_exchange
      TYPE(pw_p_type)                                    :: psi_a
      TYPE(pw_p_type), ALLOCATABLE, DIMENSION(:)         :: psi_i

      CALL timeset(routineN, handle)

      do_alpha_beta = .FALSE.
      IF (PRESENT(homo_beta) .AND. &
          PRESENT(mo_coeff_o_beta) .AND. &
          PRESENT(mo_coeff_v_beta) .AND. &
          PRESENT(Eigenval_beta) .AND. &
          PRESENT(Emp2_AB)) do_alpha_beta = .TRUE.

      ! initialize and create the matrix (ia|jnu)
      CALL cp_dbcsr_init(matrix_ia_jnu)
      CALL cp_dbcsr_create(matrix_ia_jnu, template=mo_coeff_o)

      ! Allocate Sparse matrices: (ia|jb)
      CALL cp_dbcsr_init(matrix_ia_jb)
      CALL cp_dbcsr_m_by_n_from_template(matrix_ia_jb, template=mo_coeff_o, m=homo, n=nmo-homo, &
                                         sym=dbcsr_type_no_symmetry)

      ! set all to zero in such a way that the memory is actually allocated
      CALL cp_dbcsr_set(matrix_ia_jnu, 0.0_dp)
      CALL cp_dbcsr_set(matrix_ia_jb, 0.0_dp)
      CALL cp_dbcsr_set(mat_munu%matrix, 0.0_dp)

      IF (calc_ex) THEN
         ! create the analogous of matrix_ia_jb in fm type
         NULLIFY (fm_BIb_jb)
         NULLIFY (fm_struct)
         CALL cp_dbcsr_get_info(matrix_ia_jb, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
         CALL cp_fm_struct_create(fm_struct, context=blacs_env_sub, nrow_global=nfullrows_total, &
                                  ncol_global=nfullcols_total, para_env=para_env_sub)
         CALL cp_fm_create(fm_BIb_jb, fm_struct, name="fm_BIb_jb")

         CALL copy_dbcsr_to_fm(matrix_ia_jb, fm_BIb_jb)
         CALL cp_fm_struct_release(fm_struct)

         CALL cp_fm_get_info(matrix=fm_BIb_jb, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)

         max_row_col_local = MAX(nrow_local, ncol_local)
         CALL mp_max(max_row_col_local, para_env_sub%group)

         ALLOCATE (local_col_row_info(0:max_row_col_local, 2))
         local_col_row_info = 0
         ! 0,1 nrows
         local_col_row_info(0, 1) = nrow_local
         local_col_row_info(1:nrow_local, 1) = row_indices(1:nrow_local)
         ! 0,2 ncols
         local_col_row_info(0, 2) = ncol_local
         local_col_row_info(1:ncol_local, 2) = col_indices(1:ncol_local)
      END IF

      IF (do_alpha_beta) THEN
         ! initialize and create the matrix (ia|jnu)
         CALL cp_dbcsr_init(matrix_ia_jnu_beta)
         CALL cp_dbcsr_create(matrix_ia_jnu_beta, template=mo_coeff_o_beta)

         ! Allocate Sparse matrices: (ia|jb)
         CALL cp_dbcsr_init(matrix_ia_jb_beta)
         CALL cp_dbcsr_m_by_n_from_template(matrix_ia_jb_beta, template=mo_coeff_o_beta, m=homo_beta, n=nmo-homo_beta, &
                                            sym=dbcsr_type_no_symmetry)

         virtual_beta = nmo-homo_beta

         CALL cp_dbcsr_set(matrix_ia_jnu_beta, 0.0_dp)
         CALL cp_dbcsr_set(matrix_ia_jb_beta, 0.0_dp)
      END IF

      CALL m_memory(mem)
      mem_real = (mem+1024*1024-1)/(1024*1024)
      ! mp_min .... a hack.. it should be mp_max, but as it turns out, on some processes the previously freed memory (hfx)
      ! has not been given back to the OS yet.
      CALL mp_min(mem_real, para_env%group)

      virtual = nmo-homo

      wfn_size = REAL(SIZE(rho_r%pw%cr3d), KIND=dp)
      CALL mp_max(wfn_size, para_env%group)

      ngroup = para_env%num_pe/para_env_sub%num_pe

      ! calculate the minimal memory required per MPI task (p=occupied division,q=virtual division)
      p_best = ngroup
      q_best = 1
      mem_min = HUGE(0)
      DO p = 1, ngroup
         q = ngroup/p
         IF (p*q .NE. ngroup) CYCLE

         CALL estimate_memory_usage(wfn_size, p, q, para_env_sub%num_pe, nmo, virtual, homo, calc_ex, mem_try)

         IF (mem_try <= mem_min) THEN
            mem_min = mem_try
            p_best = p
            q_best = q
         END IF
      END DO
      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T69,F9.2,A3)') 'Minimum required memory per MPI process for MP2:', &
         mem_min, ' MB'

      mem_real = mp2_memory-mem_real
      mem_real = MAX(mem_real, mem_min)
      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T69,F9.2,A3)') 'Available memory per MPI process for MP2:', &
         mem_real, ' MB'

      wfn_calc_best = HUGE(wfn_calc_best)
      DO p = 1, ngroup
         q = ngroup/p
         IF (p*q .NE. ngroup) CYCLE

         CALL estimate_memory_usage(wfn_size, p, q, para_env_sub%num_pe, nmo, virtual, homo, calc_ex, mem_try)

         IF (mem_try > mem_real) CYCLE
         wfn_calc = ((homo+p-1)/p)+((virtual+q-1)/q)
         IF (wfn_calc < wfn_calc_best) THEN
            wfn_calc_best = wfn_calc
            p_best = p
            q_best = q
         ENDIF
      ENDDO

      max_batch_size_I = (homo+p_best-1)/p_best
      max_batch_size_A = (virtual+q_best-1)/q_best

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,T77,i4)") &
            "MP2_GPW| max. batch size for the occupied states:", max_batch_size_I
         WRITE (UNIT=unit_nr, FMT="(T3,A,T77,i4)") &
            "MP2_GPW| max. batch size for the virtual states:", max_batch_size_A
      ENDIF

      ALLOCATE (vector_batch_I_size_group(0:p_best-1))
      ALLOCATE (vector_batch_A_size_group(0:q_best-1))

      vector_batch_I_size_group = max_batch_size_I
      IF (SUM(vector_batch_I_size_group) /= homo) THEN
         one = 1
         IF (SUM(vector_batch_I_size_group) > homo) one = -1
         i = -1
         DO
            i = i+1
            vector_batch_I_size_group(i) = vector_batch_I_size_group(i)+one
            IF (SUM(vector_batch_I_size_group) == homo) EXIT
            IF (i == p_best-1) i = -1
         END DO
      END IF

      vector_batch_A_size_group = max_batch_size_A
      IF (SUM(vector_batch_A_size_group) /= virtual) THEN
         one = 1
         IF (SUM(vector_batch_A_size_group) > virtual) one = -1
         i = -1
         DO
            i = i+1
            vector_batch_A_size_group(i) = vector_batch_A_size_group(i)+one
            IF (SUM(vector_batch_A_size_group) == virtual) EXIT
            IF (i == q_best-1) i = -1
         END DO
      END IF

      !XXXXXXXXXXXXX inverse group distribution
      group_counter = 0
      a_group_counter = 0
      my_A_virtual_start = 1
      DO j = 0, q_best-1
         my_I_occupied_start = 1
         i_group_counter = 0
         DO i = 0, p_best-1
            group_counter = group_counter+1
            IF (color_sub == group_counter-1) EXIT
            my_I_occupied_start = my_I_occupied_start+vector_batch_I_size_group(i)
            i_group_counter = i_group_counter+1
         END DO
         my_q_position = j
         IF (color_sub == group_counter-1) EXIT
         my_A_virtual_start = my_A_virtual_start+vector_batch_A_size_group(j)
         a_group_counter = a_group_counter+1
      END DO
      !XXXXXXXXXXXXX inverse group distribution

      my_I_occupied_end = my_I_occupied_start+vector_batch_I_size_group(i_group_counter)-1
      my_I_batch_size = vector_batch_I_size_group(i_group_counter)
      my_A_virtual_end = my_A_virtual_start+vector_batch_A_size_group(a_group_counter)-1
      my_A_batch_size = vector_batch_A_size_group(a_group_counter)

      DEALLOCATE (vector_batch_I_size_group)
      DEALLOCATE (vector_batch_A_size_group)

      ! replicate on a local array on proc 0 the occupied and virtual wavevectior
      ! needed for the calculation of the WF's by calculate_wavefunction
      ! (external vector)
      CALL grep_occ_virt_wavefunc(para_env_sub, nmo, &
                                  my_I_occupied_start, my_I_occupied_end, my_I_batch_size, &
                                  my_A_virtual_start, my_A_virtual_end, my_A_batch_size, &
                                  mo_coeff_o, mo_coeff_v, my_Cocc, my_Cvirt)

      ! divide the b states in the sub_group in such a way to create
      ! b_start and b_end for each proc inside the sub_group
      max_b_size = (virtual+para_env_sub%num_pe-1)/para_env_sub%num_pe
      ALLOCATE (vector_B_sizes(0:para_env_sub%num_pe-1))
      vector_B_sizes = max_b_size
      IF (SUM(vector_B_sizes) /= virtual) THEN
         one = 1
         IF (SUM(vector_B_sizes) > virtual) one = -1
         i = -1
         DO
            i = i+1
            vector_B_sizes(i) = vector_B_sizes(i)+one
            IF (SUM(vector_B_sizes) == virtual) EXIT
            IF (i == para_env_sub%num_pe-1) i = -1
         END DO
      END IF
      ! now give to each proc its b_start and b_end
      b_group_counter = 0
      my_B_virtual_start = 1
      DO j = 0, para_env_sub%num_pe-1
         b_group_counter = b_group_counter+1
         IF (b_group_counter-1 == para_env_sub%mepos) EXIT
         my_B_virtual_start = my_B_virtual_start+vector_B_sizes(j)
      END DO
      my_B_virtual_end = my_B_virtual_start+vector_B_sizes(para_env_sub%mepos)-1
      my_B_size = vector_B_sizes(para_env_sub%mepos)

      DEALLOCATE (vector_B_sizes)

      ! create an array containing a different "color" for each pair of
      ! A_start and B_start, communication will take place only among
      ! those proc that have the same A_start and B_start
      ALLOCATE (color_array(0:para_env_sub%num_pe-1, 0:q_best-1))
      color_array = 0
      color_counter = 0
      DO j = 0, q_best-1
         DO i = 0, para_env_sub%num_pe-1
            color_counter = color_counter+1
            color_array(i, j) = color_counter
         END DO
      END DO
      sub_sub_color = color_array(para_env_sub%mepos, my_q_position)

      DEALLOCATE (color_array)

      ! now create a group that contains all the proc that have the same 2 virtual starting points
      ! in this way it is possible to sum the common integrals needed for the full MP2 energy
      ! in mp_comm_split_direct the color is give by my_a_virtual_start and my_b_virtual_start
      CALL mp_comm_split_direct(para_env%group, comm_exchange, sub_sub_color)
      NULLIFY (para_env_exchange)
      CALL cp_para_env_create(para_env_exchange, comm_exchange)

      ! crate the proc maps
      ALLOCATE (proc_map(-para_env_exchange%num_pe:2*para_env_exchange%num_pe-1))
      DO i = 0, para_env_exchange%num_pe-1
         proc_map(i) = i
         proc_map(-i-1) = para_env_exchange%num_pe-i-1
         proc_map(para_env_exchange%num_pe+i) = i
      END DO

      ALLOCATE (sub_proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
      DO i = 0, para_env_sub%num_pe-1
         sub_proc_map(i) = i
         sub_proc_map(-i-1) = para_env_sub%num_pe-i-1
         sub_proc_map(para_env_sub%num_pe+i) = i
      END DO

      ! create an array containing the information for communication
      ALLOCATE (exchange_group_sizes(0:para_env_exchange%num_pe-1, 3))
      exchange_group_sizes = 0
      exchange_group_sizes(para_env_exchange%mepos, 1) = my_I_occupied_start
      exchange_group_sizes(para_env_exchange%mepos, 2) = my_I_occupied_end
      exchange_group_sizes(para_env_exchange%mepos, 3) = my_I_batch_size
      CALL mp_sum(exchange_group_sizes, para_env_exchange%group)
      mepos_in_EX_group = para_env_exchange%mepos
      size_of_exchange_group = para_env_exchange%num_pe

      NULLIFY (psi_a%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, psi_a%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)

      ALLOCATE (psi_i(my_I_occupied_start:my_I_occupied_end))
      DO i = my_I_occupied_start, my_I_occupied_end
         NULLIFY (psi_i(i)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, psi_i(i)%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL calculate_wavefunction(mo_coeff, i, psi_i(i), rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, &
                                     pw_env_sub, external_vector=my_Cocc(:, i-my_I_occupied_start+1))
      END DO

      Emp2 = 0.0_dp
      Emp2_Cou = 0.0_dp
      Emp2_EX = 0.0_dp
      IF (do_alpha_beta) Emp2_AB = 0.0_dp
      IF (calc_ex) THEN
         ALLOCATE (BIb_C(my_B_size, homo, my_I_batch_size))
      END IF

      CALL timeset(routineN//"_loop", handle2)
      DO a = homo+my_A_virtual_start, homo+my_A_virtual_end

         IF (calc_ex) BIb_C = 0.0_dp

         ! psi_a
         CALL calculate_wavefunction(mo_coeff, a, psi_a, rho_g, atomic_kind_set, &
                                     qs_kind_set, cell, dft_control, particle_set, &
                                     pw_env_sub, external_vector=my_Cvirt(:, a-(homo+my_A_virtual_start)+1))
         i_counter = 0
         DO i = my_I_occupied_start, my_I_occupied_end
            i_counter = i_counter+1

            ! potential
            CALL timeset(routineN//"_pot", handle3)
            rho_r%pw%cr3d = psi_i(i)%pw%cr3d*psi_a%pw%cr3d
            CALL pw_transfer(rho_r%pw, rho_g%pw)
            CALL pw_poisson_solve(poisson_env, rho_g%pw, pair_energy, pot_g%pw)
            CALL pw_transfer(pot_g%pw, rho_r%pw)
            CALL pw_scale(rho_r%pw, rho_r%pw%pw_grid%dvol)
            CALL timestop(handle3)

            ! and finally (ia|munu)
            CALL timeset(routineN//"_int", handle3)
            CALL cp_dbcsr_set(mat_munu%matrix, 0.0_dp)
            CALL integrate_v_rspace(rho_r, hmat=mat_munu, qs_env=qs_env, &
                                    calculate_forces=.FALSE., compute_tau=.FALSE., gapw=.FALSE., &
                                    pw_env_external=pw_env_sub, task_list_external=task_list_sub)
            CALL timestop(handle3)

            ! multiply and goooooooo ...
            CALL timeset(routineN//"_mult_o", handle3)
            CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o, &
                                   0.0_dp, matrix_ia_jnu, filter_eps=eps_filter)
            IF (do_alpha_beta) THEN
               ! transform orbitals using the beta coeff matrix
               CALL cp_dbcsr_multiply("N", "N", 1.0_dp, mat_munu%matrix, mo_coeff_o_beta, &
                                      0.0_dp, matrix_ia_jnu_beta, filter_eps=eps_filter)
            END IF
            CALL timestop(handle3)
            CALL timeset(routineN//"_mult_v", handle3)
            CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_ia_jnu, mo_coeff_v, &
                                   0.0_dp, matrix_ia_jb, filter_eps=eps_filter)
            IF (do_alpha_beta) THEN
               ! transform orbitals using the beta coeff matrix
               CALL cp_dbcsr_multiply("T", "N", 1.0_dp, matrix_ia_jnu_beta, mo_coeff_v_beta, &
                                      0.0_dp, matrix_ia_jb_beta, filter_eps=eps_filter)
            END IF
            CALL timestop(handle3)

            CALL timeset(routineN//"_E_Cou", handle3)
            CALL cp_dbcsr_iterator_start(iter, matrix_ia_jb)
            DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
               CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                                 row_size=row_size, col_size=col_size, &
                                                 row_offset=row_offset, col_offset=col_offset)
               DO b = 1, col_size
               DO j = 1, row_size
                  ! Compute the coulomb MP2 energy
                  Emp2_Cou = Emp2_Cou-2.0_dp*data_block(j, b)**2/ &
                             (Eigenval(a)+Eigenval(homo+col_offset+b-1)-Eigenval(i)-Eigenval(row_offset+j-1))
               ENDDO
               ENDDO
            ENDDO
            CALL cp_dbcsr_iterator_stop(iter)
            IF (do_alpha_beta) THEN
               ! Compute the coulomb only= SO = MP2 alpha-beta  MP2 energy component
               CALL cp_dbcsr_iterator_start(iter, matrix_ia_jb_beta)
               DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
                  CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                                    row_size=row_size, col_size=col_size, &
                                                    row_offset=row_offset, col_offset=col_offset)
                  DO b = 1, col_size
                  DO j = 1, row_size
                     ! Compute the coulomb MP2 energy alpha beta case
                     Emp2_AB = Emp2_AB-data_block(j, b)**2/ &
                               (Eigenval(a)+Eigenval_beta(homo_beta+col_offset+b-1)-Eigenval(i)-Eigenval_beta(row_offset+j-1))
                  ENDDO
                  ENDDO
               ENDDO
               CALL cp_dbcsr_iterator_stop(iter)
            END IF
            CALL timestop(handle3)

            ! now collect my local data from all the other members of the group
            ! b_start, b_end
            IF (calc_ex) THEN
               CALL timeset(routineN//"_E_Ex_1", handle3)
               CALL copy_dbcsr_to_fm(matrix_ia_jb, fm_BIb_jb)
               CALL grep_my_integrals(para_env_sub, fm_BIb_jb, BIb_C(1:my_B_size, 1:homo, i_counter), max_row_col_local, &
                                      sub_proc_map, local_col_row_info, &
                                      my_B_virtual_end, my_B_virtual_start)
               CALL timestop(handle3)
            END IF

         END DO

         IF (calc_ex) THEN
            CALL timeset(routineN//"_E_Ex_2", handle3)
            ! calculate the contribution to MP2 energy for my local data
            DO i = 1, my_I_batch_size
               DO j = my_I_occupied_start, my_I_occupied_end
                  DO b = 1, my_B_size
                     b_global = b-1+my_B_virtual_start
                     Emp2_EX = Emp2_EX+BIb_C(b, j, i)*BIb_C(b, i+my_I_occupied_start-1, j-my_I_occupied_start+1) &
                               /(Eigenval(a)+Eigenval(homo+b_global)-Eigenval(i+my_I_occupied_start-1)-Eigenval(j))
                  END DO
               END DO
            END DO

            ! start communicating and collecting exchange contributions from
            ! other processes in my exchange group
            DO index_proc_shift = 1, size_of_exchange_group-1
               proc_send = proc_map(mepos_in_EX_group+index_proc_shift)
               proc_receive = proc_map(mepos_in_EX_group-index_proc_shift)

               EX_start = exchange_group_sizes(proc_receive, 1)
               EX_end = exchange_group_sizes(proc_receive, 2)
               size_EX = exchange_group_sizes(proc_receive, 3)

               ALLOCATE (BIb_EX(my_B_size, my_I_batch_size, size_EX))
               BIb_EX = 0.0_dp

               EX_start_send = exchange_group_sizes(proc_send, 1)
               EX_end_send = exchange_group_sizes(proc_send, 2)
               size_EX_send = exchange_group_sizes(proc_send, 3)

               ALLOCATE (BIb_send(my_B_size, size_EX_send, my_I_batch_size))
               BIb_send(1:my_B_size, 1:size_EX_send, 1:my_I_batch_size) = &
                  BIb_C(1:my_B_size, EX_start_send:EX_end_send, 1:my_I_batch_size)

               ! send and receive the exchange array
               CALL mp_sendrecv(BIb_send, proc_send, BIb_EX, proc_receive, para_env_exchange%group)

               DO i = 1, my_I_batch_size
                  DO j = 1, size_EX
                     DO b = 1, my_B_size
                        b_global = b-1+my_B_virtual_start
                        Emp2_EX = Emp2_EX+BIb_C(b, j+EX_start-1, i)*BIb_EX(b, i, j) &
                                  /(Eigenval(a)+Eigenval(homo+b_global)-Eigenval(i+my_I_occupied_start-1)-Eigenval(j+EX_start-1))
                     END DO
                  END DO
               END DO

               DEALLOCATE (BIb_EX)
               DEALLOCATE (BIb_send)

            END DO
            CALL timestop(handle3)
         END IF

      ENDDO
      CALL timestop(handle2)

      CALL mp_sum(Emp2_Cou, para_env%group)
      CALL mp_sum(Emp2_EX, para_env%group)
      Emp2 = Emp2_Cou+Emp2_EX
      IF (do_alpha_beta) CALL mp_sum(Emp2_AB, para_env%group)

      DEALLOCATE (my_Cocc)
      DEALLOCATE (my_Cvirt)

      IF (calc_ex) THEN
         CALL cp_fm_release(fm_BIb_jb)
         DEALLOCATE (local_col_row_info)
         DEALLOCATE (BIb_C)
      END IF
      DEALLOCATE (proc_map)
      DEALLOCATE (sub_proc_map)
      DEALLOCATE (exchange_group_sizes)

      CALL cp_para_env_release(para_env_exchange)

      CALL cp_dbcsr_release(matrix_ia_jnu)
      CALL cp_dbcsr_release(matrix_ia_jb)
      IF (do_alpha_beta) THEN
         CALL cp_dbcsr_release(matrix_ia_jnu_beta)
         CALL cp_dbcsr_release(matrix_ia_jb_beta)
      END IF

      DO i = my_I_occupied_start, my_I_occupied_end
         CALL pw_release(psi_i(i)%pw)
      END DO
      DEALLOCATE (psi_i)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, psi_a%pw)

      CALL timestop(handle)

   END SUBROUTINE mp2_gpw_compute

! **************************************************************************************************
!> \brief ...
!> \param wfn_size ...
!> \param p ...
!> \param q ...
!> \param num_w ...
!> \param nmo ...
!> \param virtual ...
!> \param homo ...
!> \param calc_ex ...
!> \param mem_try ...
! **************************************************************************************************
   SUBROUTINE estimate_memory_usage(wfn_size, p, q, num_w, nmo, virtual, homo, calc_ex, mem_try)
      REAL(KIND=dp)                                      :: wfn_size
      INTEGER                                            :: p, q, num_w, nmo, virtual, homo
      LOGICAL                                            :: calc_ex
      REAL(KIND=dp)                                      :: mem_try

      mem_try = 0.0_dp
      ! integrals
      mem_try = mem_try+virtual*REAL(homo, KIND=dp)**2/(p*num_w)
      ! array for the coefficient matrix and wave vectors
      mem_try = mem_try+REAL(homo, KIND=dp)*nmo/p+ &
                REAL(virtual, KIND=dp)*nmo/q+ &
                2.0_dp*MAX(REAL(homo, KIND=dp)*nmo/p, REAL(virtual, KIND=dp)*nmo/q)
      ! temporary array for MO integrals and MO integrals to be exchanged
      IF (calc_ex) THEN
         mem_try = mem_try+2.0_dp*MAX(virtual*REAL(homo, KIND=dp)*MIN(1, num_w-1)/num_w, &
                                      virtual*REAL(homo, KIND=dp)**2/(p*p*num_w))
      ELSE
         mem_try = mem_try+2.0_dp*virtual*REAL(homo, KIND=dp)
      END IF
      ! wfn
      mem_try = mem_try+((homo+p-1)/p)*wfn_size
      ! Mb
      mem_try = mem_try*8.0D+00/1024.0D+00**2

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param para_env_sub ...
!> \param fm_BIb_jb ...
!> \param BIb_jb ...
!> \param max_row_col_local ...
!> \param proc_map ...
!> \param local_col_row_info ...
!> \param my_B_virtual_end ...
!> \param my_B_virtual_start ...
! **************************************************************************************************
   SUBROUTINE grep_my_integrals(para_env_sub, fm_BIb_jb, BIb_jb, max_row_col_local, &
                                proc_map, local_col_row_info, &
                                my_B_virtual_end, my_B_virtual_start)
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      TYPE(cp_fm_type), POINTER                          :: fm_BIb_jb
      REAL(KIND=dp), DIMENSION(:, :)                     :: BIb_jb
      INTEGER                                            :: max_row_col_local
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_col_row_info
      INTEGER                                            :: my_B_virtual_end, my_B_virtual_start

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grep_my_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i_global, iiB, j_global, jjB, ncol_rec, &
                                                            nrow_rec, proc_receive, proc_send, &
                                                            proc_shift
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: rec_col_row_info
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_rec, row_indices_rec
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_BI, rec_BI

      ALLOCATE (rec_col_row_info(0:max_row_col_local, 2))

      rec_col_row_info(:, :) = local_col_row_info

      nrow_rec = rec_col_row_info(0, 1)
      ncol_rec = rec_col_row_info(0, 2)

      ALLOCATE (row_indices_rec(nrow_rec))
      row_indices_rec = rec_col_row_info(1:nrow_rec, 1)

      ALLOCATE (col_indices_rec(ncol_rec))
      col_indices_rec = rec_col_row_info(1:ncol_rec, 2)

      ! accumulate data on BIb_jb buffer starting from myself
      DO jjB = 1, ncol_rec
         j_global = col_indices_rec(jjB)
         IF (j_global >= my_B_virtual_start .AND. j_global <= my_B_virtual_end) THEN
            DO iiB = 1, nrow_rec
               i_global = row_indices_rec(iiB)
               BIb_jb(j_global-my_B_virtual_start+1, i_global) = fm_BIb_jb%local_data(iiB, jjB)
            END DO
         END IF
      END DO

      DEALLOCATE (row_indices_rec)
      DEALLOCATE (col_indices_rec)

      IF (para_env_sub%num_pe > 1) THEN
         ALLOCATE (local_BI(nrow_rec, ncol_rec))
         local_BI(1:nrow_rec, 1:ncol_rec) = fm_BIb_jb%local_data(1:nrow_rec, 1:ncol_rec)

         DO proc_shift = 1, para_env_sub%num_pe-1
            proc_send = proc_map(para_env_sub%mepos+proc_shift)
            proc_receive = proc_map(para_env_sub%mepos-proc_shift)

            ! first exchange information on the local data
            rec_col_row_info = 0
            CALL mp_sendrecv(local_col_row_info, proc_send, rec_col_row_info, proc_receive, para_env_sub%group)
            nrow_rec = rec_col_row_info(0, 1)
            ncol_rec = rec_col_row_info(0, 2)

            ALLOCATE (row_indices_rec(nrow_rec))
            row_indices_rec = rec_col_row_info(1:nrow_rec, 1)

            ALLOCATE (col_indices_rec(ncol_rec))
            col_indices_rec = rec_col_row_info(1:ncol_rec, 2)

            ALLOCATE (rec_BI(nrow_rec, ncol_rec))
            rec_BI = 0.0_dp

            ! then send and receive the real data
            CALL mp_sendrecv(local_BI, proc_send, rec_BI, proc_receive, para_env_sub%group)

            ! accumulate the received data on BIb_jb buffer
            DO jjB = 1, ncol_rec
               j_global = col_indices_rec(jjB)
               IF (j_global >= my_B_virtual_start .AND. j_global <= my_B_virtual_end) THEN
                  DO iiB = 1, nrow_rec
                     i_global = row_indices_rec(iiB)
                     BIb_jb(j_global-my_B_virtual_start+1, i_global) = rec_BI(iiB, jjB)
                  END DO
               END IF
            END DO

            DEALLOCATE (col_indices_rec)
            DEALLOCATE (row_indices_rec)
            DEALLOCATE (rec_BI)
         END DO

         DEALLOCATE (local_BI)
      END IF

      DEALLOCATE (rec_col_row_info)

   END SUBROUTINE grep_my_integrals

! **************************************************************************************************
!> \brief ...
!> \param para_env_sub ...
!> \param dimen ...
!> \param my_I_occupied_start ...
!> \param my_I_occupied_end ...
!> \param my_I_batch_size ...
!> \param my_A_virtual_start ...
!> \param my_A_virtual_end ...
!> \param my_A_batch_size ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param my_Cocc ...
!> \param my_Cvirt ...
! **************************************************************************************************
   SUBROUTINE grep_occ_virt_wavefunc(para_env_sub, dimen, &
                                     my_I_occupied_start, my_I_occupied_end, my_I_batch_size, &
                                     my_A_virtual_start, my_A_virtual_end, my_A_batch_size, &
                                     mo_coeff_o, mo_coeff_v, my_Cocc, my_Cvirt)

      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      INTEGER :: dimen, my_I_occupied_start, my_I_occupied_end, my_I_batch_size, &
         my_A_virtual_start, my_A_virtual_end, my_A_batch_size
      TYPE(cp_dbcsr_type), POINTER                       :: mo_coeff_o, mo_coeff_v
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: my_Cocc, my_Cvirt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'grep_occ_virt_wavefunc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, col, col_offset, col_size, handle, &
                                                            i, i_global, j, j_global, row, &
                                                            row_offset, row_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      ALLOCATE (my_Cocc(dimen, my_I_batch_size))
      my_Cocc = 0.0_dp

      ALLOCATE (my_Cvirt(dimen, my_A_batch_size))
      my_Cvirt = 0.0_dp

      ! accumulate data from mo_coeff_o into Cocc
      CALL cp_dbcsr_iterator_start(iter, mo_coeff_o)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
            j_global = col_offset+j-1
            IF (j_global >= my_I_occupied_start .AND. j_global <= my_I_occupied_end) THEN
               DO i = 1, row_size
                  i_global = row_offset+i-1
                  my_Cocc(i_global, j_global-my_I_occupied_start+1) = data_block(i, j)
               END DO
            END IF
         END DO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      CALL mp_sum(my_Cocc, para_env_sub%group)

      ! accumulate data from mo_coeff_o into Cocc
      CALL cp_dbcsr_iterator_start(iter, mo_coeff_v)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)
         DO j = 1, col_size
            j_global = col_offset+j-1
            IF (j_global >= my_A_virtual_start .AND. j_global <= my_A_virtual_end) THEN
               DO i = 1, row_size
                  i_global = row_offset+i-1
                  my_Cvirt(i_global, j_global-my_A_virtual_start+1) = data_block(i, j)
               END DO
            END IF
         END DO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      CALL mp_sum(my_Cvirt, para_env_sub%group)

      CALL timestop(handle)

   END SUBROUTINE grep_occ_virt_wavefunc

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param mo_coeff ...
!> \param dimen ...
!> \param homo ...
!> \param mat_munu ...
!> \param mo_coeff_o ...
!> \param mo_coeff_v ...
!> \param mo_coeff_all ...
!> \param mo_coeff_gw ...
!> \param my_do_gw ...
!> \param gw_corr_lev_occ ...
!> \param gw_corr_lev_virt ...
!> \param only_mo_coeff_all ...
! **************************************************************************************************
   SUBROUTINE replicate_mat_to_subgroup(mp2_env, para_env, para_env_sub, mo_coeff, dimen, homo, mat_munu, &
                                        mo_coeff_o, mo_coeff_v, mo_coeff_all, mo_coeff_gw, my_do_gw, &
                                        gw_corr_lev_occ, gw_corr_lev_virt, only_mo_coeff_all)
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env, para_env_sub
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      INTEGER                                            :: dimen, homo
      TYPE(cp_dbcsr_p_type)                              :: mat_munu
      TYPE(cp_dbcsr_type), POINTER                       :: mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                                            mo_coeff_gw
      LOGICAL                                            :: my_do_gw
      INTEGER                                            :: gw_corr_lev_occ, gw_corr_lev_virt
      LOGICAL, OPTIONAL                                  :: only_mo_coeff_all

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_mat_to_subgroup', &
         routineP = moduleN//':'//routineN

      INTEGER :: blk, col, col_offset, col_size, handle, i_global, iiB, iproc, itmp(2), j_global, &
         jjB, max_row_col_local, my_mu_end, my_mu_size, my_mu_start, ncol_local, ncol_rec, &
         nrow_local, nrow_rec, proc_receive, proc_receive_static, proc_send, proc_send_static, &
         proc_shift, row, row_offset, row_size, virtual
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_array, proc_map, sizes_array, &
                                                            starts_array
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_col_row_info, rec_col_row_info
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, col_indices_rec, &
                                                            row_indices, row_indices_rec
      LOGICAL                                            :: my_only_mo_coeff_all
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C, Cgw, Cocc, Cvirt, rec_C
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block, local_C, local_C_internal
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_coeff

      my_only_mo_coeff_all = .FALSE.
      IF (PRESENT(only_mo_coeff_all)) my_only_mo_coeff_all = only_mo_coeff_all

      CALL timeset(routineN, handle)

      ALLOCATE (sizes_array(0:para_env_sub%num_pe-1))
      ALLOCATE (starts_array(0:para_env_sub%num_pe-1))
      starts_array = 0
      ALLOCATE (ends_array(0:para_env_sub%num_pe-1))
      ends_array = 0

      DO iproc = 0, para_env_sub%num_pe-1
         itmp = get_limit(dimen, para_env_sub%num_pe, iproc)
         starts_array(iproc) = itmp(1)
         ends_array(iproc) = itmp(2)
         sizes_array(iproc) = itmp(2)-itmp(1)+1
      ENDDO

      my_mu_size = sizes_array(para_env_sub%mepos)
      my_mu_start = starts_array(para_env_sub%mepos)
      my_mu_end = ends_array(para_env_sub%mepos)

      ! local storage for the C matrix
      ALLOCATE (C(my_mu_size, dimen))
      C = 0.0_dp

      ! proc_map, vector that replicate the processor numbers also
      ! for negative and positive number > num_pe
      ! needed to know which is the processor, to respect to another one,
      ! for a given shift
      ALLOCATE (proc_map(-para_env%num_pe:2*para_env%num_pe-1))
      DO iiB = 0, para_env%num_pe-1
         proc_map(iiB) = iiB
         proc_map(-iiB-1) = para_env%num_pe-iiB-1
         proc_map(para_env%num_pe+iiB) = iiB
      END DO

      CALL cp_fm_get_info(matrix=mo_coeff, &
                          matrix_struct=fm_struct_coeff, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices, &
                          local_data=local_C_internal)

      ALLOCATE (local_C(nrow_local, ncol_local))
      local_C = local_C_internal(1:nrow_local, 1:ncol_local)
      NULLIFY (local_C_internal)

      max_row_col_local = MAX(nrow_local, ncol_local)
      CALL mp_max(max_row_col_local, para_env%group)

      ALLOCATE (local_col_row_info(0:max_row_col_local, 2))
      local_col_row_info = 0
      ! 0,1 nrows
      local_col_row_info(0, 1) = nrow_local
      local_col_row_info(1:nrow_local, 1) = row_indices(1:nrow_local)
      ! 0,2 ncols
      local_col_row_info(0, 2) = ncol_local
      local_col_row_info(1:ncol_local, 2) = col_indices(1:ncol_local)

      ALLOCATE (rec_col_row_info(0:max_row_col_local, 2))

      ! accumulate data on C buffer starting from myself
      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)
         IF (i_global >= my_mu_start .AND. i_global <= my_mu_end) THEN
            DO jjB = 1, ncol_local
               j_global = col_indices(jjB)
               C(i_global-my_mu_start+1, j_global) = local_C(iiB, jjB)
            END DO
         END IF
      END DO

      ! start ring communication for collecting the data from the other
      proc_send_static = proc_map(para_env%mepos+1)
      proc_receive_static = proc_map(para_env%mepos-1)
      DO proc_shift = 1, para_env%num_pe-1
         proc_send = proc_map(para_env%mepos+proc_shift)
         proc_receive = proc_map(para_env%mepos-proc_shift)

         ! first exchange information on the local data
         rec_col_row_info = 0
         CALL mp_sendrecv(local_col_row_info, proc_send_static, rec_col_row_info, proc_receive_static, para_env%group)
         nrow_rec = rec_col_row_info(0, 1)
         ncol_rec = rec_col_row_info(0, 2)

         ALLOCATE (row_indices_rec(nrow_rec))
         row_indices_rec = rec_col_row_info(1:nrow_rec, 1)

         ALLOCATE (col_indices_rec(ncol_rec))
         col_indices_rec = rec_col_row_info(1:ncol_rec, 2)

         ALLOCATE (rec_C(nrow_rec, ncol_rec))
         rec_C = 0.0_dp

         ! then send and receive the real data
         CALL mp_sendrecv(local_C, proc_send_static, rec_C, proc_receive_static, para_env%group)

         ! accumulate the received data on C buffer
         DO iiB = 1, nrow_rec
            i_global = row_indices_rec(iiB)
            IF (i_global >= my_mu_start .AND. i_global <= my_mu_end) THEN
               DO jjB = 1, ncol_rec
                  j_global = col_indices_rec(jjB)
                  C(i_global-my_mu_start+1, j_global) = rec_C(iiB, jjB)
               END DO
            END IF
         END DO

         local_col_row_info(:, :) = rec_col_row_info
         DEALLOCATE (local_C)
         ALLOCATE (local_C(nrow_rec, ncol_rec))
         local_C = rec_C

         DEALLOCATE (col_indices_rec)
         DEALLOCATE (row_indices_rec)
         DEALLOCATE (rec_C)
      END DO

      DEALLOCATE (local_col_row_info)
      DEALLOCATE (rec_col_row_info)
      DEALLOCATE (proc_map)

      ! proc_map, for the sub_group
      ALLOCATE (proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe-1))
      DO iiB = 0, para_env_sub%num_pe-1
         proc_map(iiB) = iiB
         proc_map(-iiB-1) = para_env_sub%num_pe-iiB-1
         proc_map(para_env_sub%num_pe+iiB) = iiB
      END DO

      ! split the C matrix into occupied and virtual
      ALLOCATE (Cocc(my_mu_size, homo))
      Cocc(1:my_mu_size, 1:homo) = C(1:my_mu_size, 1:homo)

      virtual = dimen-homo
      ALLOCATE (Cvirt(my_mu_size, virtual))
      Cvirt(1:my_mu_size, 1:virtual) = C(1:my_mu_size, homo+1:dimen)

      IF (.NOT. my_only_mo_coeff_all) THEN
         ! create and fill mo_coeff_o, mo_coeff_v and mo_coeff_all
         CALL build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_o, Cocc, &
                                     homo, blk, row, col, row_size, col_size, row_offset, &
                                     col_offset, i_global, j_global, my_mu_start, my_mu_end, &
                                     mat_munu, iter, &
                                     data_block, ends_array, proc_map, &
                                     sizes_array, starts_array)

         CALL build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_v, Cvirt, &
                                     virtual, blk, row, col, row_size, col_size, row_offset, &
                                     col_offset, i_global, j_global, my_mu_start, my_mu_end, &
                                     mat_munu, iter, &
                                     data_block, ends_array, proc_map, &
                                     sizes_array, starts_array)
      ELSE
         DEALLOCATE (Cocc, Cvirt)
      END IF

      IF (my_do_gw .OR. my_only_mo_coeff_all) THEN

         IF (my_do_gw) THEN
            ! also cut levels homo-gw_corr_lev_occ+1, ..., lumo+gw_corr_lev_virt-1 of C
            ALLOCATE (Cgw(my_mu_size, gw_corr_lev_occ+gw_corr_lev_virt))
            Cgw(1:my_mu_size, 1:(gw_corr_lev_occ+gw_corr_lev_virt)) = &
               C(1:my_mu_size, homo-gw_corr_lev_occ+1:homo+gw_corr_lev_virt)
            CALL build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_gw, Cgw, &
                                        gw_corr_lev_occ+gw_corr_lev_virt, &
                                        blk, row, col, row_size, col_size, row_offset, &
                                        col_offset, i_global, j_global, my_mu_start, my_mu_end, &
                                        mat_munu, iter, &
                                        data_block, ends_array, proc_map, &
                                        sizes_array, starts_array)
         END IF

         ! all levels
         CALL build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_all, C, &
                                     dimen, blk, row, col, row_size, col_size, row_offset, &
                                     col_offset, i_global, j_global, my_mu_start, my_mu_end, &
                                     mat_munu, iter, &
                                     data_block, ends_array, proc_map, &
                                     sizes_array, starts_array)

      ELSE
         DEALLOCATE (C)
      END IF

      DEALLOCATE (proc_map)
      DEALLOCATE (sizes_array)
      DEALLOCATE (starts_array)
      DEALLOCATE (ends_array)
      DEALLOCATE (local_C)

      CALL timestop(handle)

   END SUBROUTINE replicate_mat_to_subgroup

! **************************************************************************************************
!> \brief Encapsulate the building of dbcsr_matrices mo_coeff_(v,o,all)
!> \param mp2_env ...
!> \param para_env_sub ...
!> \param mo_coeff_to_build ...
!> \param Cread ...
!> \param number_of_level ...
!> \param blk ...
!> \param row ...
!> \param col ...
!> \param row_size ...
!> \param col_size ...
!> \param row_offset ...
!> \param col_offset ...
!> \param i_global ...
!> \param j_global ...
!> \param my_mu_start ...
!> \param my_mu_end ...
!> \param mat_munu ...
!> \param iter ...
!> \param data_block ...
!> \param ends_array ...
!> \param proc_map ...
!> \param sizes_array ...
!> \param starts_array ...
!> \author Jan Wilhelm, Code by Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_to_build, Cread, &
                                     number_of_level, blk, row, col, row_size, col_size, row_offset, &
                                     col_offset, i_global, j_global, my_mu_start, my_mu_end, &
                                     mat_munu, iter, &
                                     data_block, ends_array, proc_map, &
                                     sizes_array, starts_array)
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      TYPE(cp_dbcsr_type), POINTER                       :: mo_coeff_to_build
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Cread
      INTEGER                                            :: number_of_level, blk, row, col, &
                                                            row_size, col_size, row_offset, &
                                                            col_offset, i_global, j_global, &
                                                            my_mu_start, my_mu_end
      TYPE(cp_dbcsr_p_type)                              :: mat_munu
      TYPE(cp_dbcsr_iterator)                            :: iter
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_array, proc_map, sizes_array, &
                                                            starts_array

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_mo_coeff_v_o_all', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, j, proc_receive, proc_send, &
                                                            proc_shift, rec_mu_end, rec_mu_size, &
                                                            rec_mu_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rec_C

      CALL timeset(routineN, handle)

      NULLIFY (mo_coeff_to_build)
      CALL cp_dbcsr_init_p(mo_coeff_to_build)
      CALL cp_dbcsr_m_by_n_from_row_template(mo_coeff_to_build, template=mat_munu%matrix, n=number_of_level, &
                                             sym=dbcsr_type_no_symmetry, data_type=dbcsr_type_real_default)
      CALL cp_dbcsr_reserve_all_blocks(mo_coeff_to_build)

      ! accumulate data on mo_coeff_to_build starting from myself
      CALL cp_dbcsr_iterator_start(iter, mo_coeff_to_build)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)
         DO i = 1, row_size
            i_global = row_offset+i-1
            IF (i_global >= my_mu_start .AND. i_global <= my_mu_end) THEN
               DO j = 1, col_size
                  j_global = col_offset+j-1
                  data_block(i, j) = Cread(i_global-my_mu_start+1, col_offset+j-1)
               ENDDO
            END IF
         ENDDO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      ! start ring communication in the subgroup for collecting the data from the other
      ! proc (occupied)
      DO proc_shift = 1, para_env_sub%num_pe-1
         proc_send = proc_map(para_env_sub%mepos+proc_shift)
         proc_receive = proc_map(para_env_sub%mepos-proc_shift)

         rec_mu_start = starts_array(proc_receive)
         rec_mu_end = ends_array(proc_receive)
         rec_mu_size = sizes_array(proc_receive)

         ALLOCATE (rec_C(rec_mu_size, number_of_level))
         rec_C = 0.0_dp

         ! then send and receive the real data
         CALL mp_sendrecv(Cread, proc_send, rec_C, proc_receive, para_env_sub%group)

         ! accumulate data on mo_coeff_to_build the data received from proc_rec
         CALL cp_dbcsr_iterator_start(iter, mo_coeff_to_build)
         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
            CALL cp_dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                              row_size=row_size, col_size=col_size, &
                                              row_offset=row_offset, col_offset=col_offset)
            DO i = 1, row_size
               i_global = row_offset+i-1
               IF (i_global >= rec_mu_start .AND. i_global <= rec_mu_end) THEN
                  DO j = 1, col_size
                     j_global = col_offset+j-1
                     data_block(i, j) = rec_C(i_global-rec_mu_start+1, col_offset+j-1)
                  ENDDO
               END IF
            ENDDO
         ENDDO
         CALL cp_dbcsr_iterator_stop(iter)

         DEALLOCATE (rec_C)

      END DO
      CALL cp_dbcsr_filter(mo_coeff_to_build, mp2_env%mp2_gpw%eps_filter)

      DEALLOCATE (Cread)

      CALL timestop(handle)

   END SUBROUTINE build_mo_coeff_v_o_all

! **************************************************************************************************
!> \brief Encapsulate the building of dbcsr_matrix mat_munu
!> \param mat_munu ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param dft_control ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param atom2d ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param sab_orb_sub ...
!> \param particle_set ...
!> \param cell ...
!> \param blacs_env_sub ...
!> \param do_ri_aux_basis ...
!> \param do_mixed_basis ...
!> \param group_size_prim ...
!> \param sab_orb_all ...
!> \param do_im_time ...
!> \author Jan Wilhelm, code by Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE create_mat_munu(mat_munu, qs_env, mp2_env, para_env, dft_control, atomic_kind_set, qs_kind_set, &
                              atom2d, molecule_kind_set, &
                              molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub, &
                              do_ri_aux_basis, do_mixed_basis, group_size_prim, sab_orb_all, do_im_time)

      TYPE(cp_dbcsr_p_type)                              :: mat_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub
      LOGICAL, OPTIONAL                                  :: do_ri_aux_basis, do_mixed_basis
      INTEGER, OPTIONAL                                  :: group_size_prim
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_all
      LOGICAL, OPTIONAL                                  :: do_im_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mat_munu', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blacs_grid_layout, color_sub_1, &
                                                            comm_sub_1, handle, ikind, natom, nkind
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: blacs_repeatable, my_do_im_time, &
                                                            my_do_mixed_basis, my_do_ri_aux_basis
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: orb_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      REAL(kind=dp)                                      :: subcells
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub_1
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub_im_time_P
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist_sub
      TYPE(distribution_1d_type), POINTER                :: local_molecules_sub, local_particles_sub
      TYPE(distribution_2d_type), POINTER                :: distribution_2d_sub
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set, ri_basis_set

      CALL timeset(routineN, handle)

      NULLIFY (basis_set_ri_aux)

      my_do_ri_aux_basis = .FALSE.
      IF (PRESENT(do_ri_aux_basis)) THEN
         my_do_ri_aux_basis = do_ri_aux_basis
      END IF

      my_do_mixed_basis = .FALSE.
      IF (PRESENT(do_mixed_basis)) THEN
         my_do_mixed_basis = do_mixed_basis
      END IF

      my_do_im_time = .FALSE.
      IF (PRESENT(do_im_time)) THEN
         my_do_im_time = do_im_time
      END IF

      ! hack hack hack XXXXXXXXXXXXXXX ... to be fixed
      dft_control%qs_control%eps_pgf_orb = mp2_env%mp2_gpw%eps_grid
      dft_control%qs_control%eps_rho_rspace = mp2_env%mp2_gpw%eps_grid
      dft_control%qs_control%eps_gvg_rspace = mp2_env%mp2_gpw%eps_grid
      CALL init_interaction_radii(dft_control%qs_control, atomic_kind_set, qs_kind_set)

      ! get a distribution_1d
      NULLIFY (local_particles_sub, local_molecules_sub)
      CALL distribute_molecules_1d(atomic_kind_set=atomic_kind_set, &
                                   particle_set=particle_set, &
                                   local_particles=local_particles_sub, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   local_molecules=local_molecules_sub, &
                                   force_env_section=qs_env%input)

      ! get a distribution_2d
      NULLIFY (distribution_2d_sub)
      CALL distribute_molecules_2d(cell=cell, &
                                   atomic_kind_set=atomic_kind_set, &
                                   qs_kind_set=qs_kind_set, &
                                   particle_set=particle_set, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   distribution_2d=distribution_2d_sub, &
                                   blacs_env=blacs_env_sub, &
                                   force_env_section=qs_env%input)

      ! Build the sub orbital-orbital overlap neighbor lists
      NULLIFY (sab_orb_sub)
      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (atom2d(nkind))

      CALL atom2d_build(atom2d, local_particles_sub, distribution_2d_sub, atomic_kind_set, &
                        molecule_set, molecule_only=.FALSE., particle_set=particle_set)

      ALLOCATE (orb_present(nkind))
      ALLOCATE (orb_radius(nkind))
      ALLOCATE (pair_radius(nkind, nkind))

      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            orb_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, kind_radius=orb_radius(ikind))
         ELSE
            orb_present(ikind) = .FALSE.
            orb_radius(ikind) = 0.0_dp
         ENDIF
      END DO

      CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)

      CALL build_neighbor_lists(sab_orb_sub, particle_set, atom2d, cell, pair_radius, &
                                mic=.FALSE., subcells=subcells, molecular=.FALSE., name="sab_orb_sub")
      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d)
      DEALLOCATE (orb_present, orb_radius, pair_radius)

      ! a dbcsr_dist
      ALLOCATE (dbcsr_dist_sub)
      CALL cp_dbcsr_dist2d_to_dist(distribution_2d_sub, dbcsr_dist_sub)

      ! build a dbcsr matrix the hard way
      natom = SIZE(particle_set)
      ALLOCATE (row_blk_sizes(natom))
      IF (my_do_ri_aux_basis) THEN

         ALLOCATE (basis_set_ri_aux(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, basis=basis_set_ri_aux)
         DEALLOCATE (basis_set_ri_aux)

      ELSE IF (my_do_mixed_basis) THEN

         ALLOCATE (basis_set_ri_aux(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, basis=basis_set_ri_aux)
         DEALLOCATE (basis_set_ri_aux)

         ALLOCATE (col_blk_sizes(natom))

         CALL get_particle_set(particle_set, qs_kind_set, nsgf=col_blk_sizes)
         col_blk_sizes = col_blk_sizes*group_size_prim

      ELSE
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes)
      END IF

      NULLIFY (mat_munu%matrix)
      ALLOCATE (mat_munu%matrix)
      CALL cp_dbcsr_init(mat_munu%matrix)

      IF (my_do_ri_aux_basis) THEN

         CALL cp_dbcsr_create(matrix=mat_munu%matrix, &
                              name="(ai|munu)", &
                              dist=dbcsr_dist_sub, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                              nze=0)

      ELSE IF (my_do_mixed_basis) THEN

         CALL cp_dbcsr_create(matrix=mat_munu%matrix, &
                              name="(ai|munu)", &
                              dist=dbcsr_dist_sub, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
                              nze=0)

      ELSE

         CALL cp_dbcsr_create(matrix=mat_munu%matrix, &
                              name="(ai|munu)", &
                              dist=dbcsr_dist_sub, matrix_type=dbcsr_type_symmetric, &
                              row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                              nze=0)

         IF (.NOT. my_do_im_time) THEN

            CALL cp_dbcsr_alloc_block_from_nbl(mat_munu%matrix, sab_orb_sub)

         END IF

      END IF

      DEALLOCATE (row_blk_sizes)

      IF (my_do_mixed_basis) THEN
         DEALLOCATE (col_blk_sizes)
      END IF

      CALL cp_dbcsr_distribution_release(dbcsr_dist_sub)
      DEALLOCATE (dbcsr_dist_sub)

      CALL distribution_2d_release(distribution_2d_sub)

      IF (PRESENT(sab_orb_all) .AND. my_do_im_time) THEN

         ! a para_env with groups of a single process to get the neighbor list sab_orb_all fully on every process
         color_sub_1 = para_env%mepos
         CALL mp_comm_split_direct(para_env%group, comm_sub_1, color_sub_1)
         NULLIFY (para_env_sub_im_time_P)
         CALL cp_para_env_create(para_env_sub_im_time_P, comm_sub_1)

         ! corresponding blacs_env
         blacs_grid_layout = BLACS_GRID_SQUARE
         blacs_repeatable = .TRUE.
         NULLIFY (blacs_env_sub_1)
         CALL cp_blacs_env_create(blacs_env_sub_1, para_env_sub_im_time_P, &
                                  blacs_grid_layout, &
                                  blacs_repeatable)

         ! get a distribution_2d
         NULLIFY (distribution_2d_sub)
         CALL distribute_molecules_2d(cell=cell, &
                                      atomic_kind_set=atomic_kind_set, &
                                      qs_kind_set=qs_kind_set, &
                                      particle_set=particle_set, &
                                      molecule_kind_set=molecule_kind_set, &
                                      molecule_set=molecule_set, &
                                      distribution_2d=distribution_2d_sub, &
                                      blacs_env=blacs_env_sub_1, &
                                      force_env_section=qs_env%input)

         ! Build the sub orbital-orbital overlap neighbor lists
         NULLIFY (sab_orb_all)
         CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)
         nkind = SIZE(atomic_kind_set)
         ALLOCATE (orb_present(nkind))
         ALLOCATE (orb_radius(nkind))
         ALLOCATE (pair_radius(nkind, nkind))
         ALLOCATE (atom2d(nkind))

!        CALL atom2d_build(atom2d, orb_radius, orb_present, local_particles_sub, distribution_2d_sub, &
!                          atomic_kind_set, qs_kind_set, molecule_set, molecule_only=.FALSE., dftb=.FALSE., &
!                          particle_set=particle_set)

         CALL atom2d_build(atom2d, local_particles_sub, distribution_2d_sub, &
                           atomic_kind_set, molecule_set, molecule_only=.FALSE., &
                           particle_set=particle_set)

         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=ri_basis_set, basis_type="RI_AUX")
            IF (ASSOCIATED(ri_basis_set)) THEN
               orb_present(ikind) = .TRUE.
               CALL get_gto_basis_set(gto_basis_set=ri_basis_set, kind_radius=orb_radius(ikind))
            ELSE
               orb_present(ikind) = .FALSE.
               orb_radius(ikind) = 0.0_dp
            ENDIF
         END DO

         CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)

         CALL build_neighbor_lists(sab_orb_all, particle_set, atom2d, cell, pair_radius, &
                                   mic=.FALSE., subcells=subcells, molecular=.FALSE., name="sab_orb_sub")
         CALL atom2d_cleanup(atom2d)
         DEALLOCATE (atom2d)
         DEALLOCATE (orb_present, orb_radius, pair_radius)

         CALL distribution_2d_release(distribution_2d_sub)

         CALL cp_blacs_env_release(blacs_env_sub_1)

         CALL cp_para_env_release(para_env_sub_im_time_P)

      END IF

      CALL distribution_1d_release(local_particles_sub)
      CALL distribution_1d_release(local_molecules_sub)

      CALL timestop(handle)

   END SUBROUTINE create_mat_munu

! **************************************************************************************************
!> \brief Set up dbcsr matrices for imaginary time
!> \param mat_munu ...
!> \param mat_P_local ...
!> \param mat_P_global ...
!> \param mat_M ...
!> \param mat_dm_occ_global_mao ...
!> \param mat_dm_virt_global_mao ...
!> \param mat_dm_occ_local ...
!> \param mat_dm_virt_local ...
!> \param do_mao ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param para_env ...
!> \param dft_control ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param atom2d ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param particle_set ...
!> \param cell ...
!> \param para_env_sub_im_time_P ...
!> \param blacs_env_sub_im_time_3c ...
!> \param sab_orb_all ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE create_dbcsr_matrices_im_time(mat_munu, mat_P_local, mat_P_global, mat_M, mat_dm_occ_global_mao, &
                                            mat_dm_virt_global_mao, mat_dm_occ_local, mat_dm_virt_local, do_mao, &
                                            qs_env, mp2_env, para_env, dft_control, atomic_kind_set, qs_kind_set, &
                                            atom2d, molecule_kind_set, molecule_set, particle_set, cell, &
                                            para_env_sub_im_time_P, blacs_env_sub_im_time_3c, sab_orb_all)

      TYPE(cp_dbcsr_p_type) :: mat_munu, mat_P_local, mat_P_global, mat_M, mat_dm_occ_global_mao, &
         mat_dm_virt_global_mao, mat_dm_occ_local, mat_dm_virt_local
      LOGICAL                                            :: do_mao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub_im_time_P
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub_im_time_3c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_all

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dbcsr_matrices_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER :: blacs_grid_layout, blk_end, blk_start, col_end_local, col_start_local, &
         color_sub_col, color_sub_P, color_sub_row, comm_sub_P, cut_memory, end_col, &
         end_col_data_block, end_row, end_row_data_block, group_size_P, handle, i, i_mem, icol, &
         igroup, itmp(2), j_mem, n_group_col, n_group_P, n_group_row, n_local_col_prim, &
         n_local_row_prim, nblkrows_RI, nblkrows_total, nblkrows_total_mao_occ, &
         nblkrows_total_mao_virt, nfullcols_to_split, nfullcols_total, nfullrows_to_split, &
         nfullrows_total, ngroup, offset_fullcol, offset_fullrow, ref_col, ref_row, row_end_local, &
         row_start_local, size_col
      INTEGER :: size_mao_occ, size_mao_virt, size_row, start_col, start_col_data_block, &
         start_row, start_row_data_block
      INTEGER, DIMENSION(:), POINTER :: blk_offset_mao_occ, blk_offset_mao_virt, &
         blk_sizes_mao_occ, blk_sizes_mao_virt, blk_sizes_occ, blk_sizes_virt, col_blk_offset, &
         col_blk_sizes_prim, max_col_blk_sizes, row_blk_offset, row_blk_sizes, row_blk_sizes_prim
      INTEGER, DIMENSION(:, :), POINTER                  :: col_blk_sizes_cut_memory
      LOGICAL                                            :: blacs_repeatable
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global, blacs_env_sub_P
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub

      CALL timeset(routineN, handle)

      ! release sab_orb_all
      DO i = 1, SIZE(sab_orb_all)
         CALL deallocate_neighbor_list_set(sab_orb_all(i)%neighbor_list_set)
      END DO
      DEALLOCATE (sab_orb_all)

      CALL create_mat_munu(mat_munu, qs_env, mp2_env, para_env, dft_control, atomic_kind_set, qs_kind_set, &
                           atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub_im_time_3c, &
                           do_im_time=.TRUE.)

      CALL cp_blacs_env_release(blacs_env_sub_im_time_3c)

      group_size_P = mp2_env%ri_rpa_im_time%group_size_P

      IF (group_size_P > para_env%num_pe) THEN

         group_size_P = para_env%num_pe
         mp2_env%ri_rpa_im_time%group_size_P = para_env%num_pe

      END IF

      ! only allow group_size_P which is a factor of the total number of MPI tasks
      CPASSERT(MODULO(para_env%num_pe, group_size_P) == 0)

      n_group_P = para_env%num_pe/group_size_P
      mp2_env%ri_rpa_im_time_util(1)%n_group_P = n_group_P

      ! a para_env with P groups
      color_sub_P = para_env%mepos/group_size_P
      mp2_env%ri_rpa_im_time_util(1)%color_sub_P = color_sub_P

      CALL mp_comm_split_direct(para_env%group, comm_sub_P, color_sub_P)
      NULLIFY (para_env_sub_im_time_P)
      CALL cp_para_env_create(para_env_sub_im_time_P, comm_sub_P)

      ! a blacs_env (ignore the globenv stored defaults for now)
      blacs_grid_layout = BLACS_GRID_COL
      blacs_repeatable = .TRUE.
      NULLIFY (blacs_env_sub_P)
      CALL cp_blacs_env_create(blacs_env_sub_P, para_env_sub_im_time_P, &
                               blacs_grid_layout, &
                               blacs_repeatable)

      DO i = 1, SIZE(sab_orb_sub)
         CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
      END DO
      DEALLOCATE (sab_orb_sub)

      CALL create_mat_munu(mat_P_local, qs_env, mp2_env, para_env, dft_control, &
                           atomic_kind_set, qs_kind_set, atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub_P, &
                           do_ri_aux_basis=.TRUE.)

      ! fragment n_group_P in product of integers of similar size, n_group_row*n_group_col=n_group_P
      ! employing Fermat's factorization method
      CALL generate_integer_product(n_group_P, n_group_row, n_group_col)

      ! get number of row and col blocks
      CALL cp_dbcsr_get_info(mat_munu%matrix, &
                             nblkrows_total=nblkrows_total, &
                             row_blk_size=row_blk_sizes_prim, &
                             col_blk_size=col_blk_sizes_prim, &
                             row_blk_offset=row_blk_offset, &
                             col_blk_offset=col_blk_offset, &
                             nfullrows_total=nfullrows_total, &
                             nfullcols_total=nfullcols_total)

      ! mat_munu has to be square matrix
      CPASSERT(nfullrows_total == nfullcols_total)

      ! cut the memory of mat_M
      cut_memory = mp2_env%ri_rpa_im_time%cut_memory

      ! here: memory_cut
      ALLOCATE (mp2_env%ri_rpa_im_time%sizes_array_cm(cut_memory))
      mp2_env%ri_rpa_im_time%sizes_array_cm(:) = 0
      ALLOCATE (mp2_env%ri_rpa_im_time%starts_array_cm(cut_memory))
      mp2_env%ri_rpa_im_time%starts_array_cm(:) = 0
      ALLOCATE (mp2_env%ri_rpa_im_time%ends_array_cm(cut_memory))
      mp2_env%ri_rpa_im_time%ends_array_cm(:) = 0

      ! when doing maos, we need a memory cut for the occ dm and the virt dm
      ALLOCATE (mp2_env%ri_rpa_im_time%sizes_array_cm_mao_occ(cut_memory))
      mp2_env%ri_rpa_im_time%sizes_array_cm_mao_occ(:) = 0
      ALLOCATE (mp2_env%ri_rpa_im_time%starts_array_cm_mao_occ(cut_memory))
      mp2_env%ri_rpa_im_time%starts_array_cm_mao_occ(:) = 0
      ALLOCATE (mp2_env%ri_rpa_im_time%ends_array_cm_mao_occ(cut_memory))
      mp2_env%ri_rpa_im_time%ends_array_cm_mao_occ(:) = 0

      ALLOCATE (mp2_env%ri_rpa_im_time%sizes_array_cm_mao_virt(cut_memory))
      mp2_env%ri_rpa_im_time%sizes_array_cm_mao_virt(:) = 0
      ALLOCATE (mp2_env%ri_rpa_im_time%starts_array_cm_mao_virt(cut_memory))
      mp2_env%ri_rpa_im_time%starts_array_cm_mao_virt(:) = 0
      ALLOCATE (mp2_env%ri_rpa_im_time%ends_array_cm_mao_virt(cut_memory))
      mp2_env%ri_rpa_im_time%ends_array_cm_mao_virt(:) = 0

      ngroup = cut_memory

      DO igroup = 1, ngroup

         itmp = get_limit(nfullrows_total, ngroup, igroup-1)

         CALL get_start_end_size_indx(mp2_env%ri_rpa_im_time%starts_array_cm(igroup), &
                                      mp2_env%ri_rpa_im_time%ends_array_cm(igroup), &
                                      mp2_env%ri_rpa_im_time%sizes_array_cm(igroup), &
                                      nblkrows_total, itmp(1), itmp(2), row_blk_offset, row_blk_sizes_prim)

      ENDDO

      IF (do_mao) THEN

         CALL cp_dbcsr_get_info(mat_dm_occ_global_mao%matrix, &
                                row_blk_size=blk_sizes_mao_occ, &
                                nfullrows_total=size_mao_occ, &
                                nblkrows_total=nblkrows_total_mao_occ, &
                                row_blk_offset=blk_offset_mao_occ)

         row_blk_sizes_prim => blk_sizes_mao_occ
         row_blk_offset => blk_offset_mao_occ

         CALL cp_dbcsr_get_info(mat_dm_virt_global_mao%matrix, &
                                row_blk_size=blk_sizes_mao_virt, &
                                nfullrows_total=size_mao_virt, &
                                nblkrows_total=nblkrows_total_mao_virt, &
                                row_blk_offset=blk_offset_mao_virt)

         col_blk_sizes_prim => blk_sizes_mao_virt
         col_blk_offset => blk_offset_mao_virt

         ! the same for MAOs for occ and virt
         DO igroup = 1, ngroup

            itmp = get_limit(size_mao_occ, ngroup, igroup-1)

            CALL get_start_end_size_indx(mp2_env%ri_rpa_im_time%starts_array_cm_mao_occ(igroup), &
                                         mp2_env%ri_rpa_im_time%ends_array_cm_mao_occ(igroup), &
                                         mp2_env%ri_rpa_im_time%sizes_array_cm_mao_occ(igroup), &
                                         nblkrows_total_mao_occ, itmp(1), itmp(2), &
                                         blk_offset_mao_occ, blk_sizes_mao_occ)

         ENDDO

         DO igroup = 1, ngroup

            itmp = get_limit(size_mao_virt, ngroup, igroup-1)

            CALL get_start_end_size_indx(mp2_env%ri_rpa_im_time%starts_array_cm_mao_virt(igroup), &
                                         mp2_env%ri_rpa_im_time%ends_array_cm_mao_virt(igroup), &
                                         mp2_env%ri_rpa_im_time%sizes_array_cm_mao_virt(igroup), &
                                         nblkrows_total_mao_virt, itmp(1), itmp(2), &
                                         blk_offset_mao_virt, blk_sizes_mao_virt)

         ENDDO

      ELSE

         mp2_env%ri_rpa_im_time%sizes_array_cm_mao_occ(:) = mp2_env%ri_rpa_im_time%sizes_array_cm(:)
         mp2_env%ri_rpa_im_time%starts_array_cm_mao_occ(:) = mp2_env%ri_rpa_im_time%starts_array_cm(:)
         mp2_env%ri_rpa_im_time%ends_array_cm_mao_occ(:) = mp2_env%ri_rpa_im_time%ends_array_cm(:)

         mp2_env%ri_rpa_im_time%sizes_array_cm_mao_virt(:) = mp2_env%ri_rpa_im_time%sizes_array_cm(:)
         mp2_env%ri_rpa_im_time%starts_array_cm_mao_virt(:) = mp2_env%ri_rpa_im_time%starts_array_cm(:)
         mp2_env%ri_rpa_im_time%ends_array_cm_mao_virt(:) = mp2_env%ri_rpa_im_time%ends_array_cm(:)

      END IF

      ! Cut n_group_row
      ngroup = n_group_row
      mp2_env%ri_rpa_im_time_util(1)%color_sub_row = color_sub_P/n_group_col
      mp2_env%ri_rpa_im_time_util(1)%n_group_row = n_group_row

      DO i_mem = 1, cut_memory

         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row = 0

         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_fullrow(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_fullrow = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow = 0

         DO igroup = 0, ngroup-1

! JW MAO            nfullrows_to_split = mp2_env%ri_rpa_im_time%sizes_array_cm(i_mem)

            nfullrows_to_split = mp2_env%ri_rpa_im_time%sizes_array_cm_mao_occ(i_mem)

            itmp = get_limit(nfullrows_to_split, ngroup, igroup)

            IF (i_mem == 1) THEN

               mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(igroup) = itmp(1)
               mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(igroup) = itmp(2)
               mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_fullrow(igroup) = itmp(2)-itmp(1)+1

               CALL get_blk_from_indx(indx=itmp(1), blk=blk_start, blk_offset=row_blk_offset, blk_sizes=row_blk_sizes_prim)
               CALL get_blk_from_indx(indx=itmp(2), blk=blk_end, blk_offset=row_blk_offset, blk_sizes=row_blk_sizes_prim)

               mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(igroup) = blk_start
               mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row(igroup) = blk_end
               mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(igroup) = blk_end-blk_start+1

            ELSE

               offset_fullrow = mp2_env%ri_rpa_im_time_util(i_mem-1)%ends_array_prim_fullrow(ngroup-1)

               mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(igroup) = itmp(1)+offset_fullrow
               mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(igroup) = itmp(2)+offset_fullrow
               mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_fullrow(igroup) = itmp(2)-itmp(1)+1

               CALL get_blk_from_indx(indx=(itmp(1)+offset_fullrow), blk=blk_start, &
                                      blk_offset=row_blk_offset, blk_sizes=row_blk_sizes_prim)
               CALL get_blk_from_indx(indx=(itmp(2)+offset_fullrow), blk=blk_end, &
                                      blk_offset=row_blk_offset, blk_sizes=row_blk_sizes_prim)

               mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(igroup) = blk_start
               mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row(igroup) = blk_end
               mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(igroup) = blk_end-blk_start+1

            END IF

         ENDDO

      END DO

      ! Cut n_group_col
      ngroup = n_group_col
      mp2_env%ri_rpa_im_time_util(1)%color_sub_col = MODULO(color_sub_P, n_group_col)
      mp2_env%ri_rpa_im_time_util(1)%n_group_col = n_group_col

      DO j_mem = 1, cut_memory

         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col = 0

         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_fullcol(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_fullcol = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(0:ngroup-1))
         mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol = 0

         DO igroup = 0, ngroup-1

! JW MAO            nfullcols_to_split = mp2_env%ri_rpa_im_time%sizes_array_cm(j_mem)

            nfullcols_to_split = mp2_env%ri_rpa_im_time%sizes_array_cm_mao_virt(j_mem)

            itmp = get_limit(nfullcols_to_split, ngroup, igroup)

            IF (j_mem == 1) THEN

               mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(igroup) = itmp(1)
               mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(igroup) = itmp(2)
               mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_fullcol(igroup) = itmp(2)-itmp(1)+1

               CALL get_blk_from_indx(indx=itmp(1), blk=blk_start, blk_offset=col_blk_offset, blk_sizes=col_blk_sizes_prim)
               CALL get_blk_from_indx(indx=itmp(2), blk=blk_end, blk_offset=col_blk_offset, blk_sizes=col_blk_sizes_prim)

               mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(igroup) = blk_start
               mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col(igroup) = blk_end
               mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(igroup) = blk_end-blk_start+1

            ELSE

               offset_fullcol = mp2_env%ri_rpa_im_time_util(j_mem-1)%ends_array_prim_fullcol(ngroup-1)

               mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(igroup) = itmp(1)+offset_fullcol
               mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(igroup) = itmp(2)+offset_fullcol
               mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_fullcol(igroup) = itmp(2)-itmp(1)+1

               CALL get_blk_from_indx(indx=(itmp(1)+offset_fullcol), blk=blk_start, &
                                      blk_offset=col_blk_offset, blk_sizes=col_blk_sizes_prim)
               CALL get_blk_from_indx(indx=(itmp(2)+offset_fullcol), blk=blk_end, &
                                      blk_offset=col_blk_offset, blk_sizes=col_blk_sizes_prim)

               mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(igroup) = blk_start
               mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col(igroup) = blk_end
               mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(igroup) = blk_end-blk_start+1

            END IF

         ENDDO

      END DO

      color_sub_row = mp2_env%ri_rpa_im_time_util(1)%color_sub_row
      color_sub_col = mp2_env%ri_rpa_im_time_util(1)%color_sub_col

      CALL cp_dbcsr_get_info(mat_P_local%matrix, &
                             nblkrows_total=nblkrows_RI, &
                             row_blk_size=row_blk_sizes)

      DO i_mem = 1, cut_memory

         n_local_row_prim = mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(color_sub_row)
         row_start_local = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(color_sub_row)
         row_end_local = row_start_local+n_local_row_prim-1

         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%start_row_data_block(row_start_local:row_start_local+n_local_row_prim-1))
         mp2_env%ri_rpa_im_time_util(i_mem)%start_row_data_block = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%end_row_data_block(row_start_local:row_start_local+n_local_row_prim-1))
         mp2_env%ri_rpa_im_time_util(i_mem)%end_row_data_block = 0

      END DO

      DO j_mem = 1, cut_memory

         n_local_col_prim = mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(color_sub_col)
         col_start_local = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(color_sub_col)
         col_end_local = col_start_local+n_local_col_prim-1

         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%start_col_data_block(col_start_local:col_start_local+n_local_col_prim-1))
         mp2_env%ri_rpa_im_time_util(j_mem)%start_col_data_block = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%end_col_data_block(col_start_local:col_start_local+n_local_col_prim-1))
         mp2_env%ri_rpa_im_time_util(j_mem)%end_col_data_block = 0

      END DO

      ALLOCATE (col_blk_sizes_cut_memory(cut_memory, cut_memory))
      col_blk_sizes_cut_memory(:, :) = 0

      ALLOCATE (mp2_env%ri_rpa_im_time_2d_util(cut_memory, cut_memory))

      DO i_mem = 1, cut_memory

         n_local_row_prim = mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(color_sub_row)
         row_start_local = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(color_sub_row)
         row_end_local = row_start_local+n_local_row_prim-1

         DO j_mem = 1, cut_memory

            n_local_col_prim = mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(color_sub_col)
            col_start_local = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(color_sub_col)
            col_end_local = col_start_local+n_local_col_prim-1

            ALLOCATE (mp2_env%ri_rpa_im_time_2d_util(i_mem, j_mem)%offset_combi_block( &
                      row_start_local:row_start_local+n_local_row_prim-1, &
                      col_start_local:col_start_local+n_local_col_prim-1))
            mp2_env%ri_rpa_im_time_2d_util(i_mem, j_mem)%offset_combi_block = 0

            DO icol = 1, n_local_row_prim*n_local_col_prim

               ref_row = (icol-1)/n_local_col_prim+1+row_start_local-1

               ref_col = MODULO(icol-1, n_local_col_prim)+1+col_start_local-1

               IF (ref_row == row_start_local) THEN
                  start_row = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(color_sub_row)
                  end_row = row_blk_offset(ref_row)+row_blk_sizes_prim(ref_row)-1
                  size_row = end_row-start_row+1
                  end_row_data_block = row_blk_sizes_prim(ref_row)
                  start_row_data_block = end_row_data_block-size_row+1
               ELSE IF (ref_row == row_end_local) THEN
                  start_row = row_blk_offset(ref_row)
                  end_row = mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(color_sub_row)
                  size_row = end_row-start_row+1
                  start_row_data_block = 1
                  end_row_data_block = size_row
               ELSE
                  size_row = row_blk_sizes_prim(ref_row)
                  start_row_data_block = 1
                  end_row_data_block = size_row
               END IF

               ! overwrite the whole stuff when only one block for local indices
               IF (row_start_local == row_end_local) THEN
                  start_row = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(color_sub_row)
                  end_row = mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(color_sub_row)
                  size_row = end_row-start_row+1
                  start_row_data_block = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(color_sub_row)- &
                                         row_blk_offset(ref_row)+1
                  end_row_data_block = start_row_data_block+size_row-1
               END IF

               mp2_env%ri_rpa_im_time_util(i_mem)%start_row_data_block(ref_row) = start_row_data_block
               mp2_env%ri_rpa_im_time_util(i_mem)%end_row_data_block(ref_row) = end_row_data_block

               IF (ref_col == col_start_local) THEN
                  start_col = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(color_sub_col)
                  end_col = col_blk_offset(ref_col)+col_blk_sizes_prim(ref_col)-1
                  size_col = end_col-start_col+1
                  end_col_data_block = col_blk_sizes_prim(ref_col)
                  start_col_data_block = end_col_data_block-size_col+1
               ELSE IF (ref_col == col_end_local) THEN
                  start_col = col_blk_offset(ref_col)
                  end_col = mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(color_sub_col)
                  size_col = end_col-start_col+1
                  start_col_data_block = 1
                  end_col_data_block = size_col
               ELSE
                  size_col = col_blk_sizes_prim(ref_col)
                  start_col_data_block = 1
                  end_col_data_block = size_col
               END IF

               IF (col_start_local == col_end_local) THEN
                  start_col = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(color_sub_col)
                  end_col = mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(color_sub_col)
                  size_col = end_col-start_col+1
                  start_col_data_block = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(color_sub_col)- &
                                         col_blk_offset(ref_col)+1
                  end_col_data_block = start_col_data_block+size_col-1
               END IF

               mp2_env%ri_rpa_im_time_util(j_mem)%start_col_data_block(ref_col) = start_col_data_block
               mp2_env%ri_rpa_im_time_util(j_mem)%end_col_data_block(ref_col) = end_col_data_block

               mp2_env%ri_rpa_im_time_2d_util(i_mem, j_mem)%offset_combi_block(ref_row, ref_col) = &
                  col_blk_sizes_cut_memory(i_mem, j_mem)

               col_blk_sizes_cut_memory(i_mem, j_mem) = col_blk_sizes_cut_memory(i_mem, j_mem)+size_row*size_col

            END DO

         END DO

      END DO

      ALLOCATE (max_col_blk_sizes(1))
      max_col_blk_sizes(1) = MAXVAL(col_blk_sizes_cut_memory)

      CALL create_mat_M(mat_M, blacs_env_sub_P, nblkrows_RI, row_blk_sizes, max_col_blk_sizes, &
                        mp2_env%ri_rpa_im_time_util(1)%mepos_P_from_RI_row)

      DEALLOCATE (col_blk_sizes_cut_memory, max_col_blk_sizes)

      CALL cp_blacs_env_release(blacs_env_sub_P)

      ! a blacs_env (ignore the globenv stored defaults for now)
      blacs_grid_layout = BLACS_GRID_SQUARE
      blacs_repeatable = .TRUE.
      NULLIFY (blacs_env_global)
      CALL cp_blacs_env_create(blacs_env_global, para_env, &
                               blacs_grid_layout, &
                               blacs_repeatable)

      DO i = 1, SIZE(sab_orb_sub)
         CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
      END DO
      DEALLOCATE (sab_orb_sub)

      CALL create_mat_munu(mat_P_global, qs_env, mp2_env, para_env, dft_control, &
                           atomic_kind_set, qs_kind_set, atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_global, &
                           do_ri_aux_basis=.TRUE.)

      CALL cp_dbcsr_reserve_all_blocks(mat_P_global%matrix)

      CALL cp_blacs_env_release(blacs_env_global)

      DO i = 1, SIZE(sab_orb_sub)
         CALL deallocate_neighbor_list_set(sab_orb_sub(i)%neighbor_list_set)
      END DO
      DEALLOCATE (sab_orb_sub)

      IF (do_mao) THEN
         CALL cp_dbcsr_get_info(mat_dm_occ_global_mao%matrix, &
                                row_blk_size=blk_sizes_mao_occ)
         blk_sizes_occ => blk_sizes_mao_occ
         CALL cp_dbcsr_get_info(mat_dm_virt_global_mao%matrix, &
                                row_blk_size=blk_sizes_mao_virt)
         blk_sizes_virt => blk_sizes_mao_virt
      ELSE
         CALL cp_dbcsr_get_info(mat_munu%matrix, &
                                row_blk_size=row_blk_sizes)
         blk_sizes_occ => row_blk_sizes
         blk_sizes_virt => row_blk_sizes
      END IF

      ! create mat_dm_occ/virt_local_mao (in case not using MAOs, the size is the Gauss basis)
      NULLIFY (mat_dm_occ_local%matrix)
      ALLOCATE (mat_dm_occ_local%matrix)
      CALL cp_dbcsr_init(mat_dm_occ_local%matrix)
      CALL cp_dbcsr_create(matrix=mat_dm_occ_local%matrix, template=mat_munu%matrix, &
                           name="mat_dm_occ_local", &
                           row_blk_size=blk_sizes_occ, col_blk_size=blk_sizes_occ)

      NULLIFY (mat_dm_virt_local%matrix)
      ALLOCATE (mat_dm_virt_local%matrix)
      CALL cp_dbcsr_init(mat_dm_virt_local%matrix)
      CALL cp_dbcsr_create(matrix=mat_dm_virt_local%matrix, template=mat_munu%matrix, &
                           name="mat_dm_virt_local", &
                           row_blk_size=blk_sizes_virt, col_blk_size=blk_sizes_virt)

      CALL timestop(handle)

   END SUBROUTINE create_dbcsr_matrices_im_time

! **************************************************************************************************
!> \brief ...
!> \param mat_dm_occ_global_mao ...
!> \param mat_munu_mao_occ_virt ...
!> \param mat_munu_mao_virt_occ ...
!> \param mat_dm_virt_global_mao ...
!> \param mat_munu ...
!> \param do_mao ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param mao_coeff_occ ...
!> \param mao_coeff_virt ...
!> \param mao_coeff_occ_A ...
!> \param mao_coeff_virt_A ...
!> \param matrix_s ...
!> \param mo_coeff ...
!> \param mo_coeff_beta ...
!> \param homo ...
!> \param homo_beta ...
!> \param nmo ...
!> \param nspins ...
!> \param unit_nr ...
!> \param mo_eigenvalues ...
!> \param mo_eigenvalues_beta ...
! **************************************************************************************************
   SUBROUTINE create_mao_basis_and_matrices(mat_dm_occ_global_mao, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, &
                                            mat_dm_virt_global_mao, mat_munu, do_mao, qs_env, mp2_env, &
                                            mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, &
                                            matrix_s, mo_coeff, mo_coeff_beta, homo, homo_beta, nmo, nspins, unit_nr, &
                                            mo_eigenvalues, mo_eigenvalues_beta)

      TYPE(cp_dbcsr_p_type)                              :: mat_dm_occ_global_mao, &
                                                            mat_munu_mao_occ_virt, &
                                                            mat_munu_mao_virt_occ, &
                                                            mat_dm_virt_global_mao, mat_munu
      LOGICAL                                            :: do_mao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_coeff_occ, mao_coeff_virt, &
                                                            mao_coeff_occ_A, mao_coeff_virt_A, &
                                                            matrix_s
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_beta
      INTEGER                                            :: homo, homo_beta, nmo, nspins, unit_nr
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues, mo_eigenvalues_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mao_basis_and_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_mao_occ, blk_sizes_mao_virt
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: dm_for_maos_virt, rho_ao_kp, &
                                                            scaled_dm_for_maos_occ
      TYPE(dbcsr_distribution_obj), POINTER              :: dbcsr_dist
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      IF (do_mao) THEN

         ! get density matrix
         CALL get_qs_env(qs_env, rho=rho)
         NULLIFY (rho_ao_kp)
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

         IF (mp2_env%ri_rpa_im_time%opt_sc_dm_occ) THEN
            CALL build_scaled_dm_occ(scaled_dm_for_maos_occ, rho_ao_kp, mo_coeff, mo_coeff_beta, nspins, homo, homo_beta, nmo, &
                                     mo_eigenvalues, mo_eigenvalues_beta)
         ELSE
            scaled_dm_for_maos_occ => rho_ao_kp
         END IF

         IF (mp2_env%ri_rpa_im_time%opt_sc_dm_virt) THEN

         ELSE
            dm_for_maos_virt => rho_ao_kp
         END IF

         ! get mao transformation matrix
         IF (mp2_env%ri_rpa_im_time%nmao_occ(1) .GE. 0) THEN

            CALL mao_generate_basis(qs_env, mao_coef=mao_coeff_occ, pmat_external=scaled_dm_for_maos_occ, &
                                    max_iter=mp2_env%ri_rpa_im_time%max_iter_occ, &
                                    eps_grad=mp2_env%ri_rpa_im_time%eps_grad_occ, &
                                    nmao_external=mp2_env%ri_rpa_im_time%nmao_occ, unit_nr=unit_nr)
            CALL mao_build_trafo_A(mao_coeff_occ_A, mao_coeff_occ, matrix_s, nspins)

         ELSE

            CALL allocate_and_set_identity_dbscr(mao_coeff_occ, scaled_dm_for_maos_occ, nspins)
            CALL allocate_and_set_identity_dbscr(mao_coeff_occ_A, scaled_dm_for_maos_occ, nspins)

         END IF

         IF (mp2_env%ri_rpa_im_time%nmao_virt(1) .GE. 0) THEN

            ! for the beginning, but mao_coeff_virt is not used, just use full basis for virtuals!
            CALL mao_generate_basis(qs_env, mao_coef=mao_coeff_virt, pmat_external=dm_for_maos_virt, &
                                    max_iter=mp2_env%ri_rpa_im_time%max_iter_virt, &
                                    eps_grad=mp2_env%ri_rpa_im_time%eps_grad_virt, &
                                    nmao_external=mp2_env%ri_rpa_im_time%nmao_virt, unit_nr=unit_nr)
            CALL mao_build_trafo_A(mao_coeff_virt_A, mao_coeff_virt, matrix_s, nspins)

         ELSE

            CALL allocate_and_set_identity_dbscr(mao_coeff_virt, dm_for_maos_virt, nspins)
            CALL allocate_and_set_identity_dbscr(mao_coeff_virt_A, dm_for_maos_virt, nspins)

         END IF

         ! free the scaled density matrices
         IF (mp2_env%ri_rpa_im_time%opt_sc_dm_occ) THEN
            CALL cp_dbcsr_deallocate_matrix_set(scaled_dm_for_maos_occ)
         END IF

         ! the column has the MAO basis
         CALL cp_dbcsr_get_info(mao_coeff_occ(1)%matrix, &
                                col_blk_size=blk_sizes_mao_occ)

         ! the column has the MAO basis
         CALL cp_dbcsr_get_info(mao_coeff_virt(1)%matrix, &
                                col_blk_size=blk_sizes_mao_virt)

         CALL get_qs_env(qs_env=qs_env, dbcsr_dist=dbcsr_dist)

         NULLIFY (mat_dm_occ_global_mao%matrix)
         ALLOCATE (mat_dm_occ_global_mao%matrix)
         CALL cp_dbcsr_init(mat_dm_occ_global_mao%matrix)
         CALL cp_dbcsr_create(matrix=mat_dm_occ_global_mao%matrix, &
                              name="mat_dm_occ_global_mao", &
                              dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=blk_sizes_mao_occ, col_blk_size=blk_sizes_mao_occ)

         NULLIFY (mat_dm_virt_global_mao%matrix)
         ALLOCATE (mat_dm_virt_global_mao%matrix)
         CALL cp_dbcsr_init(mat_dm_virt_global_mao%matrix)
         CALL cp_dbcsr_create(matrix=mat_dm_virt_global_mao%matrix, &
                              name="mat_dm_virt_global_mao", &
                              dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                              row_blk_size=blk_sizes_mao_virt, col_blk_size=blk_sizes_mao_virt)

         NULLIFY (mat_munu_mao_occ_virt%matrix)
         ALLOCATE (mat_munu_mao_occ_virt%matrix)
         CALL cp_dbcsr_init(mat_munu_mao_occ_virt%matrix)
         CALL cp_dbcsr_create(matrix=mat_munu_mao_occ_virt%matrix, template=mat_munu%matrix, &
                              name="mat_munu_mao_occ_virt", &
                              row_blk_size=blk_sizes_mao_occ, col_blk_size=blk_sizes_mao_virt)

         NULLIFY (mat_munu_mao_virt_occ%matrix)
         ALLOCATE (mat_munu_mao_virt_occ%matrix)
         CALL cp_dbcsr_init(mat_munu_mao_virt_occ%matrix)
         CALL cp_dbcsr_create(matrix=mat_munu_mao_virt_occ%matrix, template=mat_munu%matrix, &
                              name="mat_munu_mao_virt_occ", &
                              row_blk_size=blk_sizes_mao_virt, col_blk_size=blk_sizes_mao_occ)

      ELSE

         NULLIFY (mat_dm_occ_global_mao%matrix)
         ALLOCATE (mat_dm_occ_global_mao%matrix)
         CALL cp_dbcsr_init(mat_dm_occ_global_mao%matrix)
         CALL cp_dbcsr_create(mat_dm_occ_global_mao%matrix, &
                              name="mat_dm_occ_global_mao", &
                              template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_dm_virt_global_mao%matrix)
         ALLOCATE (mat_dm_virt_global_mao%matrix)
         CALL cp_dbcsr_init(mat_dm_virt_global_mao%matrix)
         CALL cp_dbcsr_create(mat_dm_virt_global_mao%matrix, &
                              name="mat_dm_virt_global_mao", &
                              template=matrix_s(1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param scaled_dm_for_maos_occ ...
!> \param rho_ao_kp ...
!> \param mo_coeff ...
!> \param mo_coeff_beta ...
!> \param nspins ...
!> \param homo ...
!> \param homo_beta ...
!> \param nmo ...
!> \param mo_eigenvalues ...
!> \param mo_eigenvalues_beta ...
! **************************************************************************************************
   SUBROUTINE build_scaled_dm_occ(scaled_dm_for_maos_occ, rho_ao_kp, mo_coeff, mo_coeff_beta, nspins, homo, homo_beta, nmo, &
                                  mo_eigenvalues, mo_eigenvalues_beta)

      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: scaled_dm_for_maos_occ, rho_ao_kp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_beta
      INTEGER                                            :: nspins, homo, homo_beta, nmo
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues, mo_eigenvalues_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_scaled_dm_occ', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, ispin, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: e_lumo, e_lumo_beta
      TYPE(cp_fm_type), POINTER                          :: fm_dm_occ_scaled, &
                                                            fm_mo_coeff_occ_scaled, &
                                                            fm_mo_coeff_occ_scaled_beta

      CALL timeset(routineN, handle)

      ! get info of fm_mo_coeff
      CALL cp_fm_get_info(matrix=mo_coeff, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)

      CALL cp_fm_create(fm_mo_coeff_occ_scaled, mo_coeff%matrix_struct)
      CALL cp_fm_set_all(fm_mo_coeff_occ_scaled, 0.0_dp)
      CALL cp_fm_to_fm(mo_coeff, fm_mo_coeff_occ_scaled)

      e_lumo = mo_eigenvalues(homo+1)

      DO jjB = 1, nrow_local
         DO iiB = 1, ncol_local

            i_global = col_indices(iiB)

            IF (i_global .LE. homo) THEN

               fm_mo_coeff_occ_scaled%local_data(jjB, iiB) = &
                  mo_coeff%local_data(jjB, iiB)/(e_lumo-mo_eigenvalues(i_global))

            ELSE

               fm_mo_coeff_occ_scaled%local_data(jjB, iiB) = 0.0_dp

            END IF

         END DO
      END DO

      IF (nspins == 2) THEN
         ! get info of fm_mo_coeff
         CALL cp_fm_get_info(matrix=mo_coeff_beta, &
                             nrow_local=nrow_local, &
                             ncol_local=ncol_local, &
                             row_indices=row_indices, &
                             col_indices=col_indices)

         CALL cp_fm_create(fm_mo_coeff_occ_scaled_beta, mo_coeff_beta%matrix_struct)
         CALL cp_fm_set_all(fm_mo_coeff_occ_scaled_beta, 0.0_dp)
         CALL cp_fm_to_fm(mo_coeff_beta, fm_mo_coeff_occ_scaled_beta)

         e_lumo_beta = mo_eigenvalues_beta(homo_beta+1)

         DO jjB = 1, nrow_local
            DO iiB = 1, ncol_local

               i_global = col_indices(iiB)

               IF (i_global .LE. homo) THEN

                  fm_mo_coeff_occ_scaled_beta%local_data(jjB, iiB) = &
                     mo_coeff_beta%local_data(jjB, iiB)/(e_lumo_beta-mo_eigenvalues_beta(i_global))

               ELSE

                  fm_mo_coeff_occ_scaled_beta%local_data(jjB, iiB) = 0.0_dp

               END IF

            END DO
         END DO
      END IF

      CALL cp_fm_create(fm_dm_occ_scaled, mo_coeff%matrix_struct)
      CALL cp_fm_set_all(fm_dm_occ_scaled, 0.0_dp)

      CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                   matrix_a=fm_mo_coeff_occ_scaled, matrix_b=fm_mo_coeff_occ_scaled, beta=0.0_dp, &
                   matrix_c=fm_dm_occ_scaled)

      NULLIFY (scaled_dm_for_maos_occ)
      CALL cp_dbcsr_allocate_matrix_set(scaled_dm_for_maos_occ, nspins, 1)

      DO ispin = 1, nspins
         ALLOCATE (scaled_dm_for_maos_occ(ispin, 1)%matrix)
         CALL cp_dbcsr_init(scaled_dm_for_maos_occ(ispin, 1)%matrix)
         CALL cp_dbcsr_create(matrix=scaled_dm_for_maos_occ(ispin, 1)%matrix, &
                              template=rho_ao_kp(1, 1)%matrix)
      END DO

      CALL copy_fm_to_dbcsr(fm_dm_occ_scaled, &
                            scaled_dm_for_maos_occ(1, 1)%matrix, &
                            keep_sparsity=.FALSE.)

      IF (nspins == 2) THEN
         CALL cp_gemm(transa="N", transb="T", m=nmo, n=nmo, k=nmo, alpha=1.0_dp, &
                      matrix_a=fm_mo_coeff_occ_scaled_beta, matrix_b=fm_mo_coeff_occ_scaled_beta, beta=0.0_dp, &
                      matrix_c=fm_dm_occ_scaled)
         CALL copy_fm_to_dbcsr(fm_dm_occ_scaled, &
                               scaled_dm_for_maos_occ(2, 1)%matrix, &
                               keep_sparsity=.FALSE.)
      END IF

      CALL cp_fm_release(fm_mo_coeff_occ_scaled)

      CALL cp_fm_release(fm_dm_occ_scaled)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param id_mat ...
!> \param rho_ao_kp ...
!> \param nspin ...
! **************************************************************************************************
   SUBROUTINE allocate_and_set_identity_dbscr(id_mat, rho_ao_kp, nspin)
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: id_mat
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: rho_ao_kp
      INTEGER                                            :: nspin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_and_set_identity_dbscr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      NULLIFY (id_mat)
      CALL cp_dbcsr_allocate_matrix_set(id_mat, nspin)

      DO ispin = 1, nspin
         ALLOCATE (id_mat(ispin)%matrix)
         CALL cp_dbcsr_init(id_mat(ispin)%matrix)
         CALL cp_dbcsr_create(matrix=id_mat(ispin)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, &
                              template=rho_ao_kp(1, 1)%matrix)
         CALL cp_dbcsr_reserve_diag_blocks(matrix=id_mat(ispin)%matrix)
         CALL cp_dbcsr_add_on_diag(id_mat(ispin)%matrix, 1.0_dp)
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mao_coeff_A ...
!> \param mao_coeff ...
!> \param matrix_s ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE mao_build_trafo_A(mao_coeff_A, mao_coeff, matrix_s, nspins)
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_coeff_A, mao_coeff, matrix_s
      INTEGER                                            :: nspins

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mao_build_trafo_A', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, handle, ispin, row
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_A, block_B, block_S
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      NULLIFY (mao_coeff_A)
      CALL cp_dbcsr_allocate_matrix_set(mao_coeff_A, nspins)

      DO ispin = 1, nspins
         ALLOCATE (mao_coeff_A(ispin)%matrix)
         CALL cp_dbcsr_init(mao_coeff_A(ispin)%matrix)
         CALL cp_dbcsr_create(matrix=mao_coeff_A(ispin)%matrix, &
                              template=mao_coeff(ispin)%matrix)
         CALL cp_dbcsr_reserve_diag_blocks(matrix=mao_coeff_A(ispin)%matrix)
      END DO

      DO ispin = 1, nspins

!$OMP PARALLEL DEFAULT(NONE) SHARED(mao_coeff, mao_coeff_A, matrix_s, ispin) &
!$OMP PRIVATE(iter,row,col,block_A,block_B,block_S,found)
         CALL cp_dbcsr_iterator_start(iter, mao_coeff_A(ispin)%matrix)
         DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
            CALL cp_dbcsr_iterator_next_block(iter, row, col, block_A)
            CPASSERT(row == col)

            CALL cp_dbcsr_get_block_p(matrix=mao_coeff(ispin)%matrix, row=row, col=col, block=block_B, found=found)
            CPASSERT(ASSOCIATED(block_B))

            CALL cp_dbcsr_get_block_p(matrix=matrix_s(1)%matrix, row=row, col=col, block=block_S, found=found)
            CPASSERT(ASSOCIATED(block_S))

            block_A = MATMUL(block_S, block_B)

         ENDDO
         CALL cp_dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_munu ...
!> \param mat_P_local ...
!> \param mat_P_global ...
!> \param mat_M ...
!> \param mat_dm_occ_global_mao ...
!> \param mat_dm_virt_global_mao ...
!> \param mat_munu_mao_occ_virt ...
!> \param mat_munu_mao_virt_occ ...
!> \param mat_dm_occ_local ...
!> \param mat_dm_virt_local ...
!> \param fm_matrix_L_RI_metric ...
!> \param para_env_sub_im_time_3c ...
!> \param para_env_sub_im_time_P ...
!> \param mao_coeff_occ ...
!> \param mao_coeff_virt ...
!> \param mao_coeff_occ_A ...
!> \param mao_coeff_virt_A ...
!> \param mp2_env ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE clean_up_im_time(mat_munu, mat_P_local, mat_P_global, mat_M, mat_dm_occ_global_mao, &
                               mat_dm_virt_global_mao, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, &
                               mat_dm_occ_local, mat_dm_virt_local, fm_matrix_L_RI_metric, &
                               para_env_sub_im_time_3c, para_env_sub_im_time_P, mao_coeff_occ, &
                               mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, mp2_env)

      TYPE(cp_dbcsr_p_type) :: mat_munu, mat_P_local, mat_P_global, mat_M, mat_dm_occ_global_mao, &
         mat_dm_virt_global_mao, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, mat_dm_occ_local, &
         mat_dm_virt_local
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_L_RI_metric
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub_im_time_3c, &
                                                            para_env_sub_im_time_P
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: mao_coeff_occ, mao_coeff_virt, &
                                                            mao_coeff_occ_A, mao_coeff_virt_A
      TYPE(mp2_type)                                     :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'clean_up_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_mao

      CALL timeset(routineN, handle)

      do_mao = mp2_env%ri_rpa_im_time%do_mao

      CALL cp_dbcsr_release(mat_munu%matrix)
      DEALLOCATE (mat_munu%matrix)

      CALL cp_dbcsr_release(mat_P_local%matrix)
      DEALLOCATE (mat_P_local%matrix)

      CALL cp_dbcsr_release(mat_P_global%matrix)
      DEALLOCATE (mat_P_global%matrix)

      CALL cp_dbcsr_release(mat_M%matrix)
      DEALLOCATE (mat_M%matrix)

      CALL cp_dbcsr_release(mat_dm_occ_local%matrix)
      DEALLOCATE (mat_dm_occ_local%matrix)

      CALL cp_dbcsr_release(mat_dm_virt_local%matrix)
      DEALLOCATE (mat_dm_virt_local%matrix)

      IF (do_mao) THEN
         CALL cp_dbcsr_deallocate_matrix_set(mao_coeff_occ)
         CALL cp_dbcsr_deallocate_matrix_set(mao_coeff_virt)
         CALL cp_dbcsr_deallocate_matrix_set(mao_coeff_occ_A)
         CALL cp_dbcsr_deallocate_matrix_set(mao_coeff_virt_A)
         CALL cp_dbcsr_release(mat_munu_mao_occ_virt%matrix)
         DEALLOCATE (mat_munu_mao_occ_virt%matrix)
         CALL cp_dbcsr_release(mat_munu_mao_virt_occ%matrix)
         DEALLOCATE (mat_munu_mao_virt_occ%matrix)
      END IF

      CALL cp_dbcsr_release(mat_dm_occ_global_mao%matrix)
      DEALLOCATE (mat_dm_occ_global_mao%matrix)
      CALL cp_dbcsr_release(mat_dm_virt_global_mao%matrix)
      DEALLOCATE (mat_dm_virt_global_mao%matrix)

      CALL cp_fm_release(fm_matrix_L_RI_metric)

      CALL cp_para_env_release(para_env_sub_im_time_P)

      CALL cp_para_env_release(para_env_sub_im_time_3c)

      ! should one clear not for all mem_cuts?
      DEALLOCATE (mp2_env%ri_rpa_im_time_util(1)%sizes_array_prim_col, &
                  mp2_env%ri_rpa_im_time_util(1)%starts_array_prim_col, &
                  mp2_env%ri_rpa_im_time_util(1)%ends_array_prim_col, &
                  mp2_env%ri_rpa_im_time_util(1)%sizes_array_prim_fullcol, &
                  mp2_env%ri_rpa_im_time_util(1)%starts_array_prim_fullcol, &
                  mp2_env%ri_rpa_im_time_util(1)%ends_array_prim_fullcol, &
                  mp2_env%ri_rpa_im_time_util(1)%sizes_array_prim_row, &
                  mp2_env%ri_rpa_im_time_util(1)%starts_array_prim_row, &
                  mp2_env%ri_rpa_im_time_util(1)%ends_array_prim_row, &
                  mp2_env%ri_rpa_im_time_util(1)%sizes_array_prim_fullrow, &
                  mp2_env%ri_rpa_im_time_util(1)%starts_array_prim_fullrow, &
                  mp2_env%ri_rpa_im_time_util(1)%ends_array_prim_fullrow, &
                  mp2_env%ri_rpa_im_time_2d_util(1, 1)%offset_combi_block, &
                  mp2_env%ri_rpa_im_time_util(1)%start_row_data_block, &
                  mp2_env%ri_rpa_im_time_util(1)%end_row_data_block, &
                  mp2_env%ri_rpa_im_time_util(1)%start_col_data_block, &
                  mp2_env%ri_rpa_im_time_util(1)%end_col_data_block, &
                  mp2_env%ri_rpa_im_time_util(1)%mepos_P_from_RI_row)

      DEALLOCATE (mp2_env%ri_rpa_im_time_util)

      DEALLOCATE (mp2_env%ri_rpa_im_time%sizes_array_cm, &
                  mp2_env%ri_rpa_im_time%starts_array_cm, &
                  mp2_env%ri_rpa_im_time%ends_array_cm, &
                  mp2_env%ri_rpa_im_time%sizes_array_cm_mao_occ, &
                  mp2_env%ri_rpa_im_time%starts_array_cm_mao_occ, &
                  mp2_env%ri_rpa_im_time%ends_array_cm_mao_occ, &
                  mp2_env%ri_rpa_im_time%sizes_array_cm_mao_virt, &
                  mp2_env%ri_rpa_im_time%starts_array_cm_mao_virt, &
                  mp2_env%ri_rpa_im_time%ends_array_cm_mao_virt)

      CALL timestop(handle)

   END SUBROUTINE clean_up_im_time

! **************************************************************************************************
!> \brief create mat_M, code from cp_dbcsr_dist2d_to_dist
!> \param mat_M ...
!> \param blacs_env_sub_P ...
!> \param nblkrows_total ...
!> \param row_blk_sizes ...
!> \param col_blk_sizes ...
!> \param mepos_P_from_RI_row ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE create_mat_M(mat_M, blacs_env_sub_P, nblkrows_total, row_blk_sizes, col_blk_sizes, mepos_P_from_RI_row)

      TYPE(cp_dbcsr_p_type)                              :: mat_M
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub_P
      INTEGER                                            :: nblkrows_total
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes, col_blk_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: mepos_P_from_RI_row

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mat_M', routineP = moduleN//':'//routineN

      INTEGER                                            :: col_size, handle, icol, irow, mypcol, &
                                                            myproc, myprow, nblkcols_total, &
                                                            numproc, row_size
      INTEGER, DIMENSION(:, :), POINTER                  :: col_dist_data, pgrid, row_dist_data
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_obj)                       :: dist
      TYPE(dbcsr_mp_obj)                                 :: mp_env

      CALL timeset(routineN, handle)

      CALL get_blacs_info(blacs_env_sub_P, para_env=para_env, &
                          my_process_row=myprow, my_process_column=mypcol, &
                          blacs2mpi=pgrid)

      myproc = para_env%mepos
      numproc = para_env%num_pe

      CALL cp_dbcsr_mp_new(mp_env, pgrid, para_env%group, myproc, numproc, &
                           myprow, mypcol)

      ! just round robin for row_dist_data and col_dist_data
      row_size = SIZE(pgrid, 1)
      col_size = SIZE(pgrid, 2)
      ALLOCATE (row_dist_data(nblkrows_total, 2))
      row_dist_data(:, 2) = -1

      nblkcols_total = 1

      ALLOCATE (col_dist_data(nblkcols_total, 2))
      col_dist_data(:, 2) = -1

      ALLOCATE (mepos_P_from_RI_row(nblkrows_total))

      DO irow = 1, nblkrows_total
         row_dist_data(irow, 1) = MODULO(irow, row_size)
         mepos_P_from_RI_row(irow) = MODULO(irow, row_size)
      END DO

      DO icol = 1, nblkcols_total
         col_dist_data(icol, 1) = MODULO(icol, col_size)
      END DO

      CALL cp_dbcsr_distribution_new(dist, mp_env, row_dist_data, col_dist_data)

      NULLIFY (mat_M%matrix)
      ALLOCATE (mat_M%matrix)
      CALL cp_dbcsr_init(mat_M%matrix)

      CALL cp_dbcsr_create(matrix=mat_M%matrix, &
                           name="M_P_alphadelta", &
                           dist=dist, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
                           nze=0)

      CALL cp_dbcsr_distribution_release(dist)
      CALL cp_dbcsr_mp_release(mp_env)

      DEALLOCATE (row_dist_data, col_dist_data)

      CALL timestop(handle)

   END SUBROUTINE create_mat_M

! **************************************************************************************************
!> \brief ...
!> \param num_pe ...
!> \param n_group_row ...
!> \param n_group_col ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE generate_integer_product(num_pe, n_group_row, n_group_col)

      INTEGER                                            :: num_pe, n_group_row, n_group_col

      CHARACTER(LEN=*), PARAMETER :: routineN = 'generate_integer_product', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: a_int, b_int, handle, num_pe_temp_int, &
                                                            sqrt_int, square_int, x_int
      LOGICAL                                            :: stay_while
      REAL(KIND=dp)                                      :: num_pe_real, offset_real, square_real

      CALL timeset(routineN, handle)

      ! check whether num_pe is odd
      num_pe_temp_int = num_pe
      a_int = 1
      b_int = 1

      DO WHILE ((num_pe_temp_int/2)*2 == num_pe_temp_int)

         IF ((num_pe_temp_int/4)*4 == num_pe_temp_int) THEN

            num_pe_temp_int = num_pe_temp_int/4

            a_int = a_int*2
            b_int = b_int*2

         ELSE IF ((num_pe_temp_int/2)*2 == num_pe_temp_int .AND. (num_pe_temp_int/4)*4 .NE. num_pe_temp_int) THEN

            num_pe_temp_int = num_pe_temp_int/2
            a_int = a_int*2

         END IF

      END DO

      num_pe_real = DBLE(num_pe_temp_int)

      offset_real = 0.0_dp

      stay_while = .TRUE.

      DO WHILE (stay_while)

         square_real = (CEILING(SQRT(num_pe_real))+offset_real)**2-num_pe_real

         square_int = NINT(square_real)

         sqrt_int = NINT(SQRT(square_real))

         IF (sqrt_int**2 == square_int) THEN

            stay_while = .FALSE.

         ELSE

            offset_real = offset_real+1.0_dp

         END IF

      END DO

      x_int = NINT(CEILING(SQRT(num_pe_real))+offset_real)

      n_group_row = (x_int+sqrt_int)*b_int
      n_group_col = (x_int-sqrt_int)*a_int

      ! additional balancing
      IF (n_group_row == 2*(n_group_row/2) .AND. n_group_row > 2*n_group_col) THEN

         n_group_row = n_group_row/2
         n_group_col = n_group_col*2

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param start_indx ...
!> \param end_indx ...
!> \param size_indx ...
!> \param nblkrows_total ...
!> \param itmp_1 ...
!> \param itmp_2 ...
!> \param row_blk_offset ...
!> \param row_blk_sizes ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE get_start_end_size_indx(start_indx, end_indx, size_indx, nblkrows_total, itmp_1, itmp_2, &
                                      row_blk_offset, row_blk_sizes)
      INTEGER                                            :: start_indx, end_indx, size_indx, &
                                                            nblkrows_total, itmp_1, itmp_2
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_offset, row_blk_sizes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_start_end_size_indx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk

      DO blk = 1, nblkrows_total

         IF (row_blk_offset(blk) >= itmp_1 .AND. row_blk_offset(blk) <= itmp_2) THEN

            ! check if new first block (on entry, start_block has to be 0)
            IF (start_indx == 0) THEN
               start_indx = row_blk_offset(blk)
            END IF

            end_indx = row_blk_offset(blk)+row_blk_sizes(blk)-1

         END IF

      END DO

      size_indx = end_indx-start_indx+1

      ! have a check that if there is nothing to be done for the specific memory_cut, then we know it
      IF (start_indx == 0 .AND. end_indx == 0) THEN
         size_indx = 0
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param indx ...
!> \param blk ...
!> \param blk_offset ...
!> \param blk_sizes ...
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE get_blk_from_indx(indx, blk, blk_offset, blk_sizes)

      INTEGER                                            :: indx, blk
      INTEGER, DIMENSION(:), POINTER                     :: blk_offset, blk_sizes

      INTEGER                                            :: iblk, nblkrows_total

      nblkrows_total = SIZE(blk_sizes)

      DO iblk = 1, nblkrows_total

         IF (blk_offset(iblk) <= indx .AND. blk_offset(iblk)+blk_sizes(iblk)-1 >= indx) THEN

            blk = iblk

         END IF

      END DO

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param eps_pgf_orb_old ...
!> \param eps_rho_rspace_old ...
!> \param eps_gvg_rspace_old ...
! **************************************************************************************************
   SUBROUTINE get_eps_old(dft_control, eps_pgf_orb_old, eps_rho_rspace_old, eps_gvg_rspace_old)

      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(kind=dp)                                      :: eps_pgf_orb_old, eps_rho_rspace_old, &
                                                            eps_gvg_rspace_old

      ! re-init the radii to be able to generate pair lists with MP2-appropriate screening
      eps_pgf_orb_old = dft_control%qs_control%eps_pgf_orb
      eps_rho_rspace_old = dft_control%qs_control%eps_rho_rspace
      eps_gvg_rspace_old = dft_control%qs_control%eps_gvg_rspace

   END SUBROUTINE

END MODULE mp2_gpw
