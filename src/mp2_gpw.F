!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calls routines to get RI integrals and calculate total energies
!> \par History
!>      10.2011 created [Joost VandeVondele and Mauro Del Ben]
!>      07.2019 split from mp2_gpw.F [Frederick Stein]
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
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_dist2d_to_dist,&
                                              cp_dbcsr_m_by_n_from_row_template,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
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
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_clear_mempools, dbcsr_copy, dbcsr_create, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_get_block_p, &
        dbcsr_get_info, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_reserve_diag_blocks, dbcsr_type, &
        dbcsr_type_no_symmetry, dbcsr_type_real_default, dbcsr_type_symmetric
   USE dbcsr_tensor_api,                ONLY: dbcsr_t_type
   USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                              distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE distribution_methods,            ONLY: distribute_molecules_1d,&
                                              distribute_molecules_2d
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              release_group_dist
   USE input_constants,                 ONLY: do_eri_gpw,&
                                              do_eri_mme,&
                                              do_eri_os,&
                                              do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_long,&
                                              ri_default
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE libint_wrapper,                  ONLY: cp_libint_static_cleanup,&
                                              cp_libint_static_init
   USE machine,                         ONLY: default_output_unit,&
                                              m_flush
   USE mao_basis,                       ONLY: mao_generate_basis
   USE message_passing,                 ONLY: mp_comm_split_direct,&
                                              mp_max,&
                                              mp_sendrecv
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE mp2_cphf,                        ONLY: solve_z_vector_eq
   USE mp2_gpw_method,                  ONLY: mp2_gpw_compute
   USE mp2_integrals,                   ONLY: mp2_ri_gpw_compute_in
   USE mp2_ri_gpw,                      ONLY: mp2_ri_gpw_compute_en
   USE mp2_ri_grad,                     ONLY: calc_ri_mp2_nonsep
   USE mp2_types,                       ONLY: mp2_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE rpa_main,                        ONLY: rpa_ri_compute_en
   USE rpa_rse,                         ONLY: rse_energy
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
      REAL(KIND=dp), INTENT(OUT)                         :: Emp2, Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos_mp2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: calc_forces, calc_ex
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_ri_mp2, do_ri_rpa, &
                                                            do_ri_sos_laplace_mp2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_gpw_main', routineP = moduleN//':'//routineN

      INTEGER :: blacs_grid_layout, color_sub, color_sub_3c, comm_sub, comm_sub_3c, dimen, &
         dimen_RI, dimen_RI_red, eri_method, group_size_3c, group_size_P, gw_corr_lev_occ, &
         gw_corr_lev_occ_beta, gw_corr_lev_virt, gw_corr_lev_virt_beta, handle, homo, homo_beta, &
         local_unit_nr, my_group_L_end, my_group_L_size, my_group_L_start, nelectron, &
         nelectron_beta, nmo, nspins, potential_type, ri_metric_type
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_array_mc_t, starts_array_mc_t
      LOGICAL :: blacs_repeatable, do_bse, do_dbcsr_t, do_gpw, do_im_time, do_kpoints_cubic_RPA, &
         do_mao, my_do_gw, my_do_ri_mp2, my_do_ri_rpa, my_do_ri_sos_laplace_mp2
      REAL(KIND=dp) :: Emp2_AB, Emp2_BB, Emp2_Cou_BB, Emp2_d2_AB, Emp2_d_AB, Emp2_EX_BB, &
         eps_gvg_rspace_old, eps_pgf_orb_old, eps_rho_rspace_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval, Eigenval_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C, BIb_C_beta, BIb_C_bse_ab, &
                                                            BIb_C_bse_ij, BIb_C_gw, BIb_C_gw_beta
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues, mo_eigenvalues_beta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub, blacs_env_sub_im_time_3c, &
                                                            blacs_env_sub_mat_munu, &
                                                            blacs_env_sub_RPA
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_beta
      TYPE(cp_logger_type), POINTER                      :: logger, logger_sub
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub, para_env_sub_im_time_3c, &
                                                            para_env_sub_im_time_P, &
                                                            para_env_sub_RPA
      TYPE(dbcsr_p_type) :: mat_dm_occ_global_mao, mat_dm_occ_local, mat_dm_virt_global_mao, &
         mat_dm_virt_local, mat_M, mat_munu, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, &
         mat_P_global, mat_P_local
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coeff_occ, mao_coeff_occ_A, &
                                                            mao_coeff_virt, mao_coeff_virt_A, &
                                                            matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp
      TYPE(dbcsr_p_type), DIMENSION(:, :, :), POINTER    :: mat_3c_overl_int, &
                                                            mat_3c_overl_int_mao_for_occ, &
                                                            mat_3c_overl_int_mao_for_virt
      TYPE(dbcsr_t_type)                                 :: t_3c_M
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :)   :: t_3c_O, t_3c_overl_int
      TYPE(dbcsr_type), POINTER :: mo_coeff_all, mo_coeff_all_beta, mo_coeff_gw, mo_coeff_gw_beta, &
         mo_coeff_o, mo_coeff_o_beta, mo_coeff_v, mo_coeff_v_beta
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(group_dist_d1_type)                           :: gd_array, gd_B_all, gd_B_occ_bse, &
                                                            gd_B_virt_bse, gd_B_virtual, &
                                                            gd_B_virtual_beta
      TYPE(kpoint_type), POINTER                         :: kpoints, kpoints_from_DFT
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_all, sab_orb_sub
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

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

      ! GW and SOS-MP2 cannot be used together
      IF (my_do_ri_sos_laplace_mp2) THEN
         CPASSERT(.NOT. mp2_env%ri_rpa%do_ri_g0w0)
      END IF

      ! check if we want to do imaginary time
      do_im_time = mp2_env%do_im_time
      do_dbcsr_t = mp2_env%ri_rpa_im_time%do_dbcsr_t
      do_mao = mp2_env%ri_rpa_im_time%do_mao
      do_bse = qs_env%mp2_env%ri_g0w0%do_bse
      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      ! Get the number of spins
      nspins = SIZE(mos_mp2)

      ! ... setup needed to be able to qs_integrate in a subgroup.
      IF (do_kpoints_cubic_RPA) THEN
         CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, kpoints=kpoints_from_DFT)
         mos(1:nspins) => kpoints_from_DFT%kp_env(1)%kpoint_env%mos(1:nspins, 1)
      ELSE
         CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, mos=mos)
      END IF
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
      IF (para_env%ionode) THEN
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

      eri_method = qs_env%mp2_env%eri_method

      IF (do_im_time) THEN

         group_size_3c = mp2_env%ri_rpa_im_time%group_size_3c

         IF (group_size_3c > para_env%num_pe) THEN

            group_size_3c = para_env%num_pe
            mp2_env%ri_rpa_im_time%group_size_3c = para_env%num_pe

         END IF

         IF (mp2_env%ri_rpa_im_time%group_size_P > para_env%num_pe) THEN
            mp2_env%ri_rpa_im_time%group_size_P = para_env%num_pe
         END IF

         group_size_P = mp2_env%ri_rpa_im_time%group_size_P

         ! only allow group_size_P which is a factor of the total number of MPI tasks
         CPASSERT(MODULO(para_env%num_pe, group_size_P) == 0)

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

         IF (eri_method == do_eri_gpw) THEN
            blacs_env_sub_mat_munu => blacs_env_sub
         ELSE
            blacs_env_sub_mat_munu => blacs_env_sub_im_time_3c
         END IF

         IF (mp2_env%ri_metric%potential_type == ri_default) mp2_env%ri_metric%potential_type = do_potential_id

         ! statically initialize libint
         CALL cp_libint_static_init()

      END IF

      IF (mp2_env%ri_metric%potential_type == ri_default) THEN
         mp2_env%ri_metric%potential_type = do_potential_coulomb
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
                      matrix_s_kp=matrix_s_kp)

      matrix_s(1:1) => matrix_s_kp(1:1, 1)

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
      ri_metric_type = mp2_env%ri_metric%potential_type

      ! which interaction potential
      potential_type = mp2_env%potential_parameter%potential_type

      ! whether we need gpw integrals (plus pw stuff)
      do_gpw = (eri_method == do_eri_gpw) .OR. &
               ((potential_type == do_potential_long .OR. ri_metric_type == do_potential_long) &
                .AND. eri_method == do_eri_os) &
               .OR. (ri_metric_type == do_potential_id .AND. eri_method == do_eri_mme)

      ! check if we want to do ri-g0w0 on top of ri-rpa
      my_do_gw = mp2_env%ri_rpa%do_ri_g0w0
      gw_corr_lev_occ = mp2_env%ri_g0w0%corr_mos_occ
      gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_virt
      IF (nspins == 2) THEN
         gw_corr_lev_occ_beta = mp2_env%ri_g0w0%corr_mos_occ_beta
         gw_corr_lev_virt_beta = mp2_env%ri_g0w0%corr_mos_virt_beta
      END IF

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

            CALL replicate_mat_to_subgroup(mp2_env, para_env, para_env_sub, mo_coeff_beta, dimen, homo_beta, mat_munu, &
                                           mo_coeff_o_beta, mo_coeff_v_beta, mo_coeff_all_beta, mo_coeff_gw_beta, &
                                           my_do_gw, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta)
         END IF

      END IF

      ! now we're kind of ready to go....
      Emp2_S = 0.0_dp
      Emp2_T = 0.0_dp
      IF (my_do_ri_mp2 .OR. my_do_ri_rpa .OR. my_do_ri_sos_laplace_mp2) THEN
         ! RI-GPW integrals (same stuff for both RPA and MP2)
         IF (nspins == 2) THEN
            ! open shell case (RI) here the (ia|K) integrals are computed for both the alpha and beta components
            CALL mp2_ri_gpw_compute_in( &
               BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, gd_array, gd_B_virtual, dimen_RI, dimen_RI_red, qs_env, &
               para_env, para_env_sub, color_sub, cell, particle_set, &
               atomic_kind_set, qs_kind_set, mo_coeff, fm_matrix_L_RI_metric, nmo, homo, &
               mat_munu, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, sab_orb_sub, sab_orb_all, &
               mo_coeff_o, mo_coeff_v, mo_coeff_all, mo_coeff_gw, &
               mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, mp2_env%calc_PQ_cond_num, calc_forces, blacs_env_sub, my_do_gw .AND. .NOT. do_im_time, &
               do_bse, gd_B_all, starts_array_mc_t, ends_array_mc_t, gw_corr_lev_occ, gw_corr_lev_virt, &
               do_im_time, do_mao, do_kpoints_cubic_RPA, kpoints, &
               mat_3c_overl_int, do_dbcsr_t, t_3c_overl_int, t_3c_M, t_3c_O, &
               mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, &
               mao_coeff_occ, mao_coeff_virt, mp2_env%ri_metric, &
               gd_B_occ_bse, gd_B_virt_bse, &
               BIb_C_beta, BIb_C_gw_beta, gd_B_virtual_beta, &
               homo_beta, mo_coeff_o_beta, mo_coeff_v_beta, mo_coeff_all_beta, mo_coeff_gw_beta)
         ELSE
            ! closed shell case (RI)
            CALL mp2_ri_gpw_compute_in(BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, gd_array, gd_B_virtual, &
                                       dimen_RI, dimen_RI_red, qs_env, para_env, para_env_sub, &
                                       color_sub, cell, particle_set, &
                                       atomic_kind_set, qs_kind_set, mo_coeff, fm_matrix_L_RI_metric, nmo, homo, &
                                       mat_munu, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, sab_orb_sub, &
                                       sab_orb_all, mo_coeff_o, mo_coeff_v, mo_coeff_all, mo_coeff_gw, &
                                       mp2_env%mp2_gpw%eps_filter, unit_nr, &
                                       mp2_env%mp2_memory, mp2_env%calc_PQ_cond_num, calc_forces, &
                                       blacs_env_sub, my_do_gw .AND. .NOT. do_im_time, do_bse, gd_B_all, &
                                       starts_array_mc_t, ends_array_mc_t, &
                                       gw_corr_lev_occ, gw_corr_lev_virt, &
                                       do_im_time, do_mao, do_kpoints_cubic_RPA, kpoints, &
                                       mat_3c_overl_int, do_dbcsr_t, t_3c_overl_int, t_3c_M, t_3c_O, mat_3c_overl_int_mao_for_occ, &
                                       mat_3c_overl_int_mao_for_virt, mao_coeff_occ, mao_coeff_virt, &
                                       mp2_env%ri_metric, gd_B_occ_bse, gd_B_virt_bse)

         END IF
      ELSE
         ! Canonical MP2-GPW
         IF (nspins == 2) THEN
            ! alpha-alpha and alpha-beta components
            IF (unit_nr > 0) WRITE (unit_nr, *)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Alpha (ia|'
            CALL mp2_gpw_compute( &
               Emp2, Emp2_Cou, Emp2_EX, qs_env, para_env, para_env_sub, color_sub, &
               cell, particle_set, &
               atomic_kind_set, qs_kind_set, mo_coeff, Eigenval, nmo, homo, mat_munu, &
               sab_orb_sub, mo_coeff_o, mo_coeff_v, mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, calc_ex, blacs_env_sub, &
               homo_beta, mo_coeff_o_beta, mo_coeff_v_beta, Eigenval_beta, Emp2_AB)

            ! beta-beta component
            IF (unit_nr > 0) WRITE (unit_nr, *)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Beta (ia|'
            CALL mp2_gpw_compute( &
               Emp2_BB, Emp2_Cou_BB, Emp2_EX_BB, qs_env, para_env, para_env_sub, color_sub, cell, particle_set, &
               atomic_kind_set, qs_kind_set, mo_coeff_beta, Eigenval_beta, nmo, homo_beta, mat_munu, &
               sab_orb_sub, mo_coeff_o_beta, mo_coeff_v_beta, mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, calc_ex, blacs_env_sub)

         ELSE
            ! closed shell case
            CALL mp2_gpw_compute( &
               Emp2, Emp2_Cou, Emp2_EX, qs_env, para_env, para_env_sub, color_sub, cell, particle_set, &
               atomic_kind_set, qs_kind_set, mo_coeff, Eigenval, nmo, homo, mat_munu, &
               sab_orb_sub, mo_coeff_o, mo_coeff_v, mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, calc_ex, blacs_env_sub)
         END IF
      END IF

      ! Free possibly large buffers allocated by dbcsr on the GPU,
      ! large hybrid dgemm/pdgemm's coming later will need the space.
      CALL dbcsr_clear_mempools()

      IF (calc_forces) THEN
         ! make a copy of mo_coeff_o and mo_coeff_v
         NULLIFY (mp2_env%ri_grad%mo_coeff_o)
         CALL dbcsr_init_p(mp2_env%ri_grad%mo_coeff_o)
         CALL dbcsr_copy(mp2_env%ri_grad%mo_coeff_o, mo_coeff_o, name="mo_coeff_o")
         NULLIFY (mp2_env%ri_grad%mo_coeff_v)
         CALL dbcsr_init_p(mp2_env%ri_grad%mo_coeff_v)
         CALL dbcsr_copy(mp2_env%ri_grad%mo_coeff_v, mo_coeff_v, name="mo_coeff_v")
         IF (nspins == 2) THEN
            NULLIFY (mp2_env%ri_grad%mo_coeff_o_beta)
            CALL dbcsr_init_p(mp2_env%ri_grad%mo_coeff_o_beta)
            CALL dbcsr_copy(mp2_env%ri_grad%mo_coeff_o_beta, mo_coeff_o_beta, name="mo_coeff_o_b")
            NULLIFY (mp2_env%ri_grad%mo_coeff_v_beta)
            CALL dbcsr_init_p(mp2_env%ri_grad%mo_coeff_v_beta)
            CALL dbcsr_copy(mp2_env%ri_grad%mo_coeff_v_beta, mo_coeff_v_beta, name="mo_coeff_v_b")
         ENDIF
         CALL get_group_dist(gd_array, color_sub, my_group_L_start, my_group_L_end, my_group_L_size)
      END IF
      ! Copy mo coeffs for RPA AXK
      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         NULLIFY (mp2_env%ri_rpa%mo_coeff_o)
         CALL dbcsr_init_p(mp2_env%ri_rpa%mo_coeff_o)
         CALL dbcsr_copy(mp2_env%ri_rpa%mo_coeff_o, mo_coeff_o, name="mo_coeff_o")
         NULLIFY (mp2_env%ri_rpa%mo_coeff_v)
         CALL dbcsr_init_p(mp2_env%ri_rpa%mo_coeff_v)
         CALL dbcsr_copy(mp2_env%ri_rpa%mo_coeff_v, mo_coeff_v, name="mo_coeff_v")
      ENDIF

      IF (.NOT. do_im_time) THEN

         CALL dbcsr_release(mo_coeff_o)
         DEALLOCATE (mo_coeff_o)
         CALL dbcsr_release(mo_coeff_v)
         DEALLOCATE (mo_coeff_v)
         IF (my_do_gw) THEN
            CALL dbcsr_release(mo_coeff_all)
            DEALLOCATE (mo_coeff_all)
         END IF

         IF (nspins == 2) THEN
            CALL dbcsr_release(mo_coeff_o_beta)
            DEALLOCATE (mo_coeff_o_beta)
            CALL dbcsr_release(mo_coeff_v_beta)
            DEALLOCATE (mo_coeff_v_beta)
            IF (my_do_gw) THEN
               CALL dbcsr_release(mo_coeff_all_beta)
               DEALLOCATE (mo_coeff_all_beta)
            END IF
         END IF

      END IF

      IF (.NOT. calc_forces) THEN
         IF (.NOT. mp2_env%ri_rpa%do_ri_axk) THEN

            CALL dbcsr_release(mat_munu%matrix)
            DEALLOCATE (mat_munu%matrix)

            CALL release_neighbor_list_sets(sab_orb_sub)

         ENDIF

      END IF

      ! decide if to do RI-RPA or RI-MP2
      IF (my_do_ri_rpa .OR. my_do_ri_sos_laplace_mp2) THEN

         IF (do_im_time) THEN

            IF (eri_method == do_eri_gpw) THEN
               para_env_sub_RPA => para_env_sub
               blacs_env_sub_RPA => blacs_env_sub
            ELSE
               para_env_sub_RPA => para_env_sub_im_time_3c
               blacs_env_sub_RPA => blacs_env_sub_im_time_3c
            END IF

            ! release sab_orb_all
            CALL release_neighbor_list_sets(sab_orb_all)

            group_size_P = mp2_env%ri_rpa_im_time%group_size_P

            mp2_env%ri_rpa_im_time_util(1)%n_group_P = para_env%num_pe/group_size_P

            ! a para_env with P groups
            mp2_env%ri_rpa_im_time_util(1)%color_sub_P = para_env%mepos/group_size_P

            IF (do_dbcsr_t) THEN
               CALL create_matrix_P(mat_P_global, qs_env, mp2_env, para_env, dft_control, atomic_kind_set, &
                                    qs_kind_set, atom2d, molecule_kind_set, molecule_set, particle_set, cell)
            ELSE
               CALL create_dbcsr_matrices_im_time(mat_munu, mat_P_local, mat_P_global, mat_M, &
                                                  mat_dm_occ_global_mao, mat_dm_virt_global_mao, &
                                                  mat_dm_occ_local, mat_dm_virt_local, &
                                                  do_mao, qs_env, mp2_env, para_env, dft_control, &
                                                  atomic_kind_set, qs_kind_set, &
                                                  atom2d, molecule_kind_set, molecule_set, particle_set, cell, &
                                                  para_env_sub_im_time_P, blacs_env_sub_RPA)
            ENDIF

            CALL cp_blacs_env_release(blacs_env_sub_im_time_3c)

         ELSE

            para_env_sub_RPA => para_env_sub

         END IF

         ! RI-RPA
         IF (nspins == 2) THEN
            CALL rpa_ri_compute_en(qs_env, Emp2, mp2_env, BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, &
                                   para_env, para_env_sub_RPA, color_sub, &
                                   gd_array, gd_B_virtual, gd_B_all, gd_B_occ_bse, gd_B_virt_bse, &
                                   mo_coeff, fm_matrix_L_RI_metric, kpoints, &
                                   Eigenval, nmo, homo, dimen_RI, dimen_RI_red, gw_corr_lev_occ, gw_corr_lev_virt, &
                                   unit_nr, my_do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_mao, do_bse, matrix_s, &
                                   mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, &
                                   mat_munu, mat_dm_occ_local, mat_dm_virt_local, &
                                   mat_P_local, mat_P_global, mat_M, &
                                   mat_3c_overl_int, do_dbcsr_t, t_3c_overl_int, t_3c_M, t_3c_O, starts_array_mc_t, &
                                   ends_array_mc_t, &
                                   mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, &
                                   mp2_env%mp2_gpw%eps_filter, BIb_C_beta, homo_beta, Eigenval_beta, &
                                   gd_B_virtual_beta, &
                                   mo_coeff_beta, BIb_C_gw_beta, gw_corr_lev_occ_beta, gw_corr_lev_virt_beta)

            IF (mp2_env%ri_rpa%do_rse) &
               CALL rse_energy(qs_env, mp2_env, para_env, dft_control, &
                               mo_coeff, nmo, homo, Eigenval, &
                               Eigenval_beta, homo_beta, mo_coeff_beta)
         ELSE
            CALL rpa_ri_compute_en(qs_env, Emp2, mp2_env, BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, &
                                   para_env, para_env_sub_RPA, color_sub, &
                                   gd_array, gd_B_virtual, gd_B_all, gd_B_occ_bse, gd_B_virt_bse, &
                                   mo_coeff, fm_matrix_L_RI_metric, kpoints, &
                                   Eigenval, nmo, homo, dimen_RI, dimen_RI_red, gw_corr_lev_occ, gw_corr_lev_virt, &
                                   unit_nr, my_do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_mao, do_bse, matrix_s, &
                                   mao_coeff_occ, mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, &
                                   mat_munu, mat_dm_occ_local, mat_dm_virt_local, &
                                   mat_P_local, mat_P_global, mat_M, &
                                   mat_3c_overl_int, do_dbcsr_t, t_3c_overl_int, t_3c_M, t_3c_O, &
                                   starts_array_mc_t, ends_array_mc_t, &
                                   mat_3c_overl_int_mao_for_occ, mat_3c_overl_int_mao_for_virt, &
                                   mp2_env%mp2_gpw%eps_filter)

            IF (mp2_env%ri_rpa%do_rse) &
               CALL rse_energy(qs_env, mp2_env, para_env, dft_control, &
                               mo_coeff, nmo, homo, Eigenval)
         END IF

         IF (do_im_time) THEN
            IF (do_dbcsr_t) THEN
               CALL clean_up_im_time_t(mat_P_global, mat_dm_occ_global_mao, mat_dm_virt_global_mao, &
                                       para_env_sub_im_time_3c)
            ELSE
               CALL clean_up_im_time(mat_munu, mat_P_local, mat_P_global, mat_M, mat_dm_occ_global_mao, &
                                     mat_dm_virt_global_mao, mat_munu_mao_occ_virt, mat_munu_mao_virt_occ, &
                                     mat_dm_occ_local, mat_dm_virt_local, para_env_sub_im_time_3c, &
                                     para_env_sub_im_time_P, mao_coeff_occ, mao_coeff_virt, &
                                     mao_coeff_occ_A, mao_coeff_virt_A, mp2_env)
            ENDIF
            CALL cp_libint_static_cleanup()
         END IF

         ! Release some memory for AXK
         IF (mp2_env%ri_rpa%do_ri_axk) THEN

            CALL dbcsr_release(mat_munu%matrix)
            DEALLOCATE (mat_munu%matrix)

            CALL release_neighbor_list_sets(sab_orb_sub)

         END IF

      ELSE
         IF (my_do_ri_mp2) THEN
            ! RI-MP2-GPW compute energy
            IF (nspins == 2) THEN
               ! alpha-alpha component
               CALL mp2_ri_gpw_compute_en( &
                  Emp2, Emp2_Cou, Emp2_EX, BIb_C, mp2_env, para_env, para_env_sub, color_sub, &
                  gd_array, gd_B_virtual, &
                  Eigenval, nmo, homo, dimen_RI_red, unit_nr, calc_forces, calc_ex, &
                  open_shell_SS=.TRUE.)

               ! beta-beta component
               CALL mp2_ri_gpw_compute_en( &
                  Emp2_BB, Emp2_Cou_BB, Emp2_EX_BB, BIb_C_beta, mp2_env, para_env, para_env_sub, color_sub, &
                  gd_array, gd_B_virtual_beta, &
                  Eigenval_beta, nmo, homo_beta, dimen_RI_red, unit_nr, calc_forces, calc_ex, &
                  open_shell_SS=.TRUE.)

               ! alpha-beta case
               CALL mp2_ri_gpw_compute_en( &
                  Emp2_d_AB, Emp2_AB, Emp2_d2_AB, BIb_C, mp2_env, para_env, para_env_sub, color_sub, &
                  gd_array, gd_B_virtual, &
                  Eigenval, nmo, homo, dimen_RI_red, unit_nr, calc_forces, .FALSE., &
                  .FALSE., BIb_C_beta, homo_beta, Eigenval_beta, &
                  gd_B_virtual_beta)

            ELSE
               ! closed shell case
               CALL mp2_ri_gpw_compute_en( &
                  Emp2, Emp2_Cou, Emp2_EX, BIb_C, mp2_env, para_env, para_env_sub, color_sub, &
                  gd_array, gd_B_virtual, &
                  Eigenval, nmo, homo, dimen_RI_red, unit_nr, calc_forces, calc_ex)
            END IF
            ! if we need forces time to calculate the MP2 non-separable contribution
            ! and start coputing the Lagrangian
            IF (calc_forces) THEN

               IF (nspins == 2) THEN ! Open shell
                  CALL calc_ri_mp2_nonsep(qs_env, mp2_env, para_env, para_env_sub, cell, &
                                          particle_set, atomic_kind_set, qs_kind_set, &
                                          mo_coeff, nmo, homo, dimen_RI, Eigenval, &
                                          my_group_L_start, my_group_L_end, my_group_L_size, &
                                          sab_orb_sub, mat_munu, &
                                          blacs_env_sub, Eigenval_beta, homo_beta, mo_coeff_beta)
               ELSE ! Closed shell
                  CALL calc_ri_mp2_nonsep(qs_env, mp2_env, para_env, para_env_sub, cell, &
                                          particle_set, atomic_kind_set, qs_kind_set, &
                                          mo_coeff, nmo, homo, dimen_RI, Eigenval, &
                                          my_group_L_start, my_group_L_end, my_group_L_size, &
                                          sab_orb_sub, mat_munu, &
                                          blacs_env_sub)
               ENDIF

               CALL dbcsr_release(mp2_env%ri_grad%mo_coeff_o)
               DEALLOCATE (mp2_env%ri_grad%mo_coeff_o)

               CALL dbcsr_release(mp2_env%ri_grad%mo_coeff_v)
               DEALLOCATE (mp2_env%ri_grad%mo_coeff_v)

               IF (nspins == 2) THEN
                  CALL dbcsr_release(mp2_env%ri_grad%mo_coeff_o_beta)
                  DEALLOCATE (mp2_env%ri_grad%mo_coeff_o_beta)

                  CALL dbcsr_release(mp2_env%ri_grad%mo_coeff_v_beta)
                  DEALLOCATE (mp2_env%ri_grad%mo_coeff_v_beta)
               ENDIF

               CALL dbcsr_release(mat_munu%matrix)
               DEALLOCATE (mat_munu%matrix)

               CALL release_neighbor_list_sets(sab_orb_sub)

            END IF
         END IF

         IF (nspins == 2) THEN
            ! make order on the MP2 energy contributions
            Emp2_Cou = Emp2_Cou*0.25_dp
            Emp2_EX = Emp2_EX*0.5_dp

            Emp2_Cou_BB = Emp2_Cou_BB*0.25_dp
            Emp2_EX_BB = Emp2_EX_BB*0.5_dp

            Emp2_S = Emp2_AB
            Emp2_T = Emp2_Cou + Emp2_Cou_BB + Emp2_EX + Emp2_EX_BB

            Emp2_Cou = Emp2_Cou + Emp2_Cou_BB + Emp2_AB
            Emp2_EX = Emp2_EX + Emp2_EX_BB
            Emp2 = Emp2_EX + Emp2_Cou
         END IF

      END IF

      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXx
      ! moved from above
      IF (my_do_gw .AND. .NOT. do_im_time) THEN
         CALL dbcsr_release(mo_coeff_gw)
         DEALLOCATE (mo_coeff_gw)
         IF (nspins == 2) THEN
            CALL dbcsr_release(mo_coeff_gw_beta)
            DEALLOCATE (mo_coeff_gw_beta)
         END IF
      END IF

      ! re-init the radii to be able to generate pair lists with MP2-appropriate screening
      dft_control%qs_control%eps_pgf_orb = eps_pgf_orb_old
      dft_control%qs_control%eps_rho_rspace = eps_rho_rspace_old
      dft_control%qs_control%eps_gvg_rspace = eps_gvg_rspace_old
      CALL init_interaction_radii(dft_control%qs_control, qs_kind_set)

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
      INTEGER, INTENT(IN)                                :: dimen, homo
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_o, mo_coeff_v, mo_coeff_all, &
                                                            mo_coeff_gw
      LOGICAL, INTENT(IN)                                :: my_do_gw
      INTEGER, INTENT(IN)                                :: gw_corr_lev_occ, gw_corr_lev_virt
      LOGICAL, INTENT(IN), OPTIONAL                      :: only_mo_coeff_all

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_mat_to_subgroup', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i_global, iiB, j_global, jjB, max_row_col_local, my_mu_end, my_mu_size, &
         my_mu_start, ncol_local, ncol_rec, nrow_local, nrow_rec, proc_receive, &
         proc_receive_static, proc_send, proc_send_static, proc_shift, virtual
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: proc_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: local_col_row_info, rec_col_row_info
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, col_indices_rec, &
                                                            row_indices, row_indices_rec
      LOGICAL                                            :: my_only_mo_coeff_all
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C, Cgw, Cocc, Cvirt, rec_C
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_C, local_C_internal
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_coeff
      TYPE(group_dist_d1_type)                           :: gd_array

      my_only_mo_coeff_all = .FALSE.
      IF (PRESENT(only_mo_coeff_all)) my_only_mo_coeff_all = only_mo_coeff_all

      CALL timeset(routineN, handle)

      CALL create_group_dist(gd_array, para_env_sub%num_pe, dimen)

      CALL get_group_dist(gd_array, para_env_sub%mepos, my_mu_start, my_mu_end, my_mu_size)

      ! local storage for the C matrix
      ALLOCATE (C(my_mu_size, dimen))
      C = 0.0_dp

      ! proc_map, vector that replicate the processor numbers also
      ! for negative and positive number > num_pe
      ! needed to know which is the processor, to respect to another one,
      ! for a given shift
      ALLOCATE (proc_map(-para_env%num_pe:2*para_env%num_pe - 1))
      DO iiB = 0, para_env%num_pe - 1
         proc_map(iiB) = iiB
         proc_map(-iiB - 1) = para_env%num_pe - iiB - 1
         proc_map(para_env%num_pe + iiB) = iiB
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
               C(i_global - my_mu_start + 1, j_global) = local_C(iiB, jjB)
            END DO
         END IF
      END DO

      ! start ring communication for collecting the data from the other
      proc_send_static = proc_map(para_env%mepos + 1)
      proc_receive_static = proc_map(para_env%mepos - 1)
      DO proc_shift = 1, para_env%num_pe - 1
         proc_send = proc_map(para_env%mepos + proc_shift)
         proc_receive = proc_map(para_env%mepos - proc_shift)

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
                  C(i_global - my_mu_start + 1, j_global) = rec_C(iiB, jjB)
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
      ALLOCATE (proc_map(-para_env_sub%num_pe:2*para_env_sub%num_pe - 1))
      DO iiB = 0, para_env_sub%num_pe - 1
         proc_map(iiB) = iiB
         proc_map(-iiB - 1) = para_env_sub%num_pe - iiB - 1
         proc_map(para_env_sub%num_pe + iiB) = iiB
      END DO

      ! split the C matrix into occupied and virtual
      ALLOCATE (Cocc(my_mu_size, homo))
      Cocc(1:my_mu_size, 1:homo) = C(1:my_mu_size, 1:homo)

      virtual = dimen - homo
      ALLOCATE (Cvirt(my_mu_size, virtual))
      Cvirt(1:my_mu_size, 1:virtual) = C(1:my_mu_size, homo + 1:dimen)

      IF (.NOT. my_only_mo_coeff_all) THEN
         ! create and fill mo_coeff_o, mo_coeff_v and mo_coeff_all
         CALL build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_o, Cocc, &
                                     homo, my_mu_start, my_mu_end, &
                                     mat_munu, gd_array, proc_map)

         CALL build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_v, Cvirt, &
                                     virtual, my_mu_start, my_mu_end, &
                                     mat_munu, gd_array, proc_map)
      ELSE
         DEALLOCATE (Cocc, Cvirt)
      END IF

      IF (my_do_gw .OR. my_only_mo_coeff_all) THEN

         IF (my_do_gw) THEN
            ! also cut levels homo-gw_corr_lev_occ+1, ..., lumo+gw_corr_lev_virt-1 of C
            ALLOCATE (Cgw(my_mu_size, gw_corr_lev_occ + gw_corr_lev_virt))
            Cgw(1:my_mu_size, 1:(gw_corr_lev_occ + gw_corr_lev_virt)) = &
               C(1:my_mu_size, homo - gw_corr_lev_occ + 1:homo + gw_corr_lev_virt)
            CALL build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_gw, Cgw, &
                                        gw_corr_lev_occ + gw_corr_lev_virt, my_mu_start, my_mu_end, &
                                        mat_munu, gd_array, proc_map)
         END IF

         ! all levels
         CALL build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_all, C, &
                                     dimen, my_mu_start, my_mu_end, &
                                     mat_munu, gd_array, proc_map)

      ELSE
         DEALLOCATE (C)
      END IF

      DEALLOCATE (proc_map)
      CALL release_group_dist(gd_array)
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
!> \param my_mu_start ...
!> \param my_mu_end ...
!> \param mat_munu ...
!> \param gd_array ...
!> \param proc_map ...
!> \author Jan Wilhelm, Code by Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE build_mo_coeff_v_o_all(mp2_env, para_env_sub, mo_coeff_to_build, Cread, &
                                     number_of_level, my_mu_start, my_mu_end, &
                                     mat_munu, gd_array, proc_map)
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_to_build
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: Cread
      INTEGER, INTENT(IN)                                :: number_of_level, my_mu_start, my_mu_end
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: proc_map

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_mo_coeff_v_o_all', &
         routineP = moduleN//':'//routineN

      INTEGER :: blk, col, col_offset, col_size, handle, i, i_global, j, j_global, proc_receive, &
         proc_send, proc_shift, rec_mu_end, rec_mu_size, rec_mu_start, row, row_offset, row_size
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rec_C
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      NULLIFY (mo_coeff_to_build)
      CALL dbcsr_init_p(mo_coeff_to_build)
      CALL cp_dbcsr_m_by_n_from_row_template(mo_coeff_to_build, template=mat_munu%matrix, n=number_of_level, &
                                             sym=dbcsr_type_no_symmetry, data_type=dbcsr_type_real_default)
      CALL dbcsr_reserve_all_blocks(mo_coeff_to_build)

      ! accumulate data on mo_coeff_to_build starting from myself
      CALL dbcsr_iterator_start(iter, mo_coeff_to_build)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                        row_size=row_size, col_size=col_size, &
                                        row_offset=row_offset, col_offset=col_offset)
         DO i = 1, row_size
            i_global = row_offset + i - 1
            IF (i_global >= my_mu_start .AND. i_global <= my_mu_end) THEN
               DO j = 1, col_size
                  j_global = col_offset + j - 1
                  data_block(i, j) = Cread(i_global - my_mu_start + 1, col_offset + j - 1)
               ENDDO
            END IF
         ENDDO
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      ! start ring communication in the subgroup for collecting the data from the other
      ! proc (occupied)
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = proc_map(para_env_sub%mepos + proc_shift)
         proc_receive = proc_map(para_env_sub%mepos - proc_shift)

         CALL get_group_dist(gd_array, proc_receive, rec_mu_start, rec_mu_end, rec_mu_size)

         ALLOCATE (rec_C(rec_mu_size, number_of_level))
         rec_C = 0.0_dp

         ! then send and receive the real data
         CALL mp_sendrecv(Cread, proc_send, rec_C, proc_receive, para_env_sub%group)

         ! accumulate data on mo_coeff_to_build the data received from proc_rec
         CALL dbcsr_iterator_start(iter, mo_coeff_to_build)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, data_block, blk, &
                                           row_size=row_size, col_size=col_size, &
                                           row_offset=row_offset, col_offset=col_offset)
            DO i = 1, row_size
               i_global = row_offset + i - 1
               IF (i_global >= rec_mu_start .AND. i_global <= rec_mu_end) THEN
                  DO j = 1, col_size
                     j_global = col_offset + j - 1
                     data_block(i, j) = rec_C(i_global - rec_mu_start + 1, col_offset + j - 1)
                  ENDDO
               END IF
            ENDDO
         ENDDO
         CALL dbcsr_iterator_stop(iter)

         DEALLOCATE (rec_C)

      END DO
      CALL dbcsr_filter(mo_coeff_to_build, mp2_env%mp2_gpw%eps_filter)

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

      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_munu
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
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_ri_aux_basis, do_mixed_basis
      INTEGER, INTENT(IN), OPTIONAL                      :: group_size_prim
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_all
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_im_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mat_munu', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blacs_grid_layout, color_sub_1, &
                                                            comm_sub_1, handle, ikind, natom, nkind
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: blacs_repeatable, do_kpoints, &
                                                            my_do_im_time, my_do_mixed_basis, &
                                                            my_do_ri_aux_basis
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: orb_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      REAL(kind=dp)                                      :: subcells
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub_1
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub_im_time_P
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist_sub
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

!      do_kpoints = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints .OR. &
!                             SUM(qs_env%mp2_env%ri_rpa_im_time%kp_grid) > 0
      do_kpoints = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      IF (do_kpoints) THEN
         ! please choose EPS_PGF_ORB in QS section smaller than EPS_GRID in WFC_GPW section
         CPASSERT(mp2_env%mp2_gpw%eps_grid .GE. dft_control%qs_control%eps_pgf_orb)
      END IF

      ! hack hack hack XXXXXXXXXXXXXXX ... to be fixed
      dft_control%qs_control%eps_pgf_orb = mp2_env%mp2_gpw%eps_grid
      dft_control%qs_control%eps_rho_rspace = mp2_env%mp2_gpw%eps_grid
      dft_control%qs_control%eps_gvg_rspace = mp2_env%mp2_gpw%eps_grid
      CALL init_interaction_radii(dft_control%qs_control, qs_kind_set)

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

      ! for cubic RPA/GW with kpoints, we need all neighbors and not only the symmetric ones
      IF (do_kpoints) THEN
         CALL build_neighbor_lists(sab_orb_sub, particle_set, atom2d, cell, pair_radius, &
                                   mic=.FALSE., subcells=subcells, molecular=.FALSE., nlname="sab_orb_sub", &
                                   symmetric=.FALSE.)
      ELSE
         CALL build_neighbor_lists(sab_orb_sub, particle_set, atom2d, cell, pair_radius, &
                                   mic=.FALSE., subcells=subcells, molecular=.FALSE., nlname="sab_orb_sub")
      END IF
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

      IF (my_do_ri_aux_basis) THEN

         CALL dbcsr_create(matrix=mat_munu%matrix, &
                           name="(ai|munu)", &
                           dist=dbcsr_dist_sub, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0)

      ELSE IF (my_do_mixed_basis) THEN

         CALL dbcsr_create(matrix=mat_munu%matrix, &
                           name="(ai|munu)", &
                           dist=dbcsr_dist_sub, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
                           nze=0)

      ELSE

         CALL dbcsr_create(matrix=mat_munu%matrix, &
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

      CALL dbcsr_distribution_release(dbcsr_dist_sub)
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

         ! for kpoints in cubic RPA/GW we need all neighbor cells
         IF (do_kpoints) THEN
            CALL build_neighbor_lists(sab_orb_all, particle_set, atom2d, cell, pair_radius, &
                                      mic=.FALSE., subcells=subcells, molecular=.FALSE., nlname="sab_orb_sub", &
                                      symmetric=.FALSE.)
         ELSE
            CALL build_neighbor_lists(sab_orb_all, particle_set, atom2d, cell, pair_radius, &
                                      mic=.FALSE., subcells=subcells, molecular=.FALSE., nlname="sab_orb_sub")
         END IF
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
!> \brief ...
!> \param mat_P_global ...
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
! **************************************************************************************************
   SUBROUTINE create_matrix_P(mat_P_global, qs_env, mp2_env, para_env, dft_control, atomic_kind_set, &
                              qs_kind_set, atom2d, molecule_kind_set, molecule_set, particle_set, cell)

      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_P_global
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_matrix_P', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blacs_grid_layout, handle
      LOGICAL                                            :: blacs_repeatable
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb_sub

      CALL timeset(routineN, handle)

      blacs_grid_layout = BLACS_GRID_SQUARE
      blacs_repeatable = .TRUE.
      NULLIFY (blacs_env_global)
      CALL cp_blacs_env_create(blacs_env_global, para_env, &
                               blacs_grid_layout, &
                               blacs_repeatable)

      CALL create_mat_munu(mat_P_global, qs_env, mp2_env, para_env, dft_control, &
                           atomic_kind_set, qs_kind_set, atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_global, &
                           do_ri_aux_basis=.TRUE.)

      CALL dbcsr_reserve_all_blocks(mat_P_global%matrix)
      CALL cp_blacs_env_release(blacs_env_global)

      CALL release_neighbor_list_sets(sab_orb_sub)

      CALL timestop(handle)

   END SUBROUTINE
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
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE create_dbcsr_matrices_im_time(mat_munu, mat_P_local, mat_P_global, mat_M, mat_dm_occ_global_mao, &
                                            mat_dm_virt_global_mao, mat_dm_occ_local, mat_dm_virt_local, do_mao, &
                                            qs_env, mp2_env, para_env, dft_control, atomic_kind_set, qs_kind_set, &
                                            atom2d, molecule_kind_set, molecule_set, particle_set, cell, &
                                            para_env_sub_im_time_P, blacs_env_sub_im_time_3c)

      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_munu, mat_P_local, mat_P_global, &
                                                            mat_M
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_dm_occ_global_mao, &
                                                            mat_dm_virt_global_mao
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_dm_occ_local, mat_dm_virt_local
      LOGICAL, INTENT(IN)                                :: do_mao
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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dbcsr_matrices_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER :: blacs_grid_layout, blk_end, blk_start, col_end_local, col_start_local, &
         color_sub_col, color_sub_P, color_sub_row, comm_sub_P, cut_memory, end_col, &
         end_col_data_block, end_row, end_row_data_block, group_size_P, handle, i_mem, icol, &
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

      CALL create_mat_munu(mat_munu, qs_env, mp2_env, para_env, dft_control, atomic_kind_set, qs_kind_set, &
                           atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub_im_time_3c, &
                           do_im_time=.TRUE.)

      color_sub_P = mp2_env%ri_rpa_im_time_util(1)%color_sub_P
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

      CALL release_neighbor_list_sets(sab_orb_sub)

      CALL create_mat_munu(mat_P_local, qs_env, mp2_env, para_env, dft_control, &
                           atomic_kind_set, qs_kind_set, atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub_P, &
                           do_ri_aux_basis=.TRUE.)

      group_size_P = mp2_env%ri_rpa_im_time%group_size_P
      n_group_P = para_env%num_pe/group_size_P
      ! fragment n_group_P in product of integers of similar size, n_group_row*n_group_col=n_group_P
      ! employing Fermat's factorization method
      CALL generate_integer_product(n_group_P, n_group_row, n_group_col)

      ! get number of row and col blocks
      CALL dbcsr_get_info(mat_munu%matrix, &
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

         itmp = get_limit(nfullrows_total, ngroup, igroup - 1)

         CALL get_start_end_size_indx(mp2_env%ri_rpa_im_time%starts_array_cm(igroup), &
                                      mp2_env%ri_rpa_im_time%ends_array_cm(igroup), &
                                      mp2_env%ri_rpa_im_time%sizes_array_cm(igroup), &
                                      nblkrows_total, itmp(1), itmp(2), row_blk_offset, row_blk_sizes_prim)

      ENDDO

      IF (do_mao) THEN

         CALL dbcsr_get_info(mat_dm_occ_global_mao%matrix, &
                             row_blk_size=blk_sizes_mao_occ, &
                             nfullrows_total=size_mao_occ, &
                             nblkrows_total=nblkrows_total_mao_occ, &
                             row_blk_offset=blk_offset_mao_occ)

         row_blk_sizes_prim => blk_sizes_mao_occ
         row_blk_offset => blk_offset_mao_occ

         CALL dbcsr_get_info(mat_dm_virt_global_mao%matrix, &
                             row_blk_size=blk_sizes_mao_virt, &
                             nfullrows_total=size_mao_virt, &
                             nblkrows_total=nblkrows_total_mao_virt, &
                             row_blk_offset=blk_offset_mao_virt)

         col_blk_sizes_prim => blk_sizes_mao_virt
         col_blk_offset => blk_offset_mao_virt

         ! the same for MAOs for occ and virt
         DO igroup = 1, ngroup

            itmp = get_limit(size_mao_occ, ngroup, igroup - 1)

            CALL get_start_end_size_indx(mp2_env%ri_rpa_im_time%starts_array_cm_mao_occ(igroup), &
                                         mp2_env%ri_rpa_im_time%ends_array_cm_mao_occ(igroup), &
                                         mp2_env%ri_rpa_im_time%sizes_array_cm_mao_occ(igroup), &
                                         nblkrows_total_mao_occ, itmp(1), itmp(2), &
                                         blk_offset_mao_occ, blk_sizes_mao_occ)

         ENDDO

         DO igroup = 1, ngroup

            itmp = get_limit(size_mao_virt, ngroup, igroup - 1)

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

         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row = 0

         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_fullrow(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_fullrow = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow = 0

         DO igroup = 0, ngroup - 1

            nfullrows_to_split = mp2_env%ri_rpa_im_time%sizes_array_cm_mao_occ(i_mem)

            itmp = get_limit(nfullrows_to_split, ngroup, igroup)

            IF (i_mem == 1) THEN

               mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(igroup) = itmp(1)
               mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(igroup) = itmp(2)
               mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_fullrow(igroup) = itmp(2) - itmp(1) + 1

               CALL get_blk_from_indx(indx=itmp(1), blk=blk_start, blk_offset=row_blk_offset, blk_sizes=row_blk_sizes_prim)
               CALL get_blk_from_indx(indx=itmp(2), blk=blk_end, blk_offset=row_blk_offset, blk_sizes=row_blk_sizes_prim)

               mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(igroup) = blk_start
               mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row(igroup) = blk_end
               mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(igroup) = blk_end - blk_start + 1

            ELSE

               offset_fullrow = mp2_env%ri_rpa_im_time_util(i_mem - 1)%ends_array_prim_fullrow(ngroup - 1)

               mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(igroup) = itmp(1) + offset_fullrow
               mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(igroup) = itmp(2) + offset_fullrow
               mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_fullrow(igroup) = itmp(2) - itmp(1) + 1

               CALL get_blk_from_indx(indx=(itmp(1) + offset_fullrow), blk=blk_start, &
                                      blk_offset=row_blk_offset, blk_sizes=row_blk_sizes_prim)
               CALL get_blk_from_indx(indx=(itmp(2) + offset_fullrow), blk=blk_end, &
                                      blk_offset=row_blk_offset, blk_sizes=row_blk_sizes_prim)

               mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(igroup) = blk_start
               mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_row(igroup) = blk_end
               mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(igroup) = blk_end - blk_start + 1

            END IF

         ENDDO

      END DO

      ! Cut n_group_col
      ngroup = n_group_col
      mp2_env%ri_rpa_im_time_util(1)%color_sub_col = MODULO(color_sub_P, n_group_col)
      mp2_env%ri_rpa_im_time_util(1)%n_group_col = n_group_col

      DO j_mem = 1, cut_memory

         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col = 0

         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_fullcol(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_fullcol = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(0:ngroup - 1))
         mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol = 0

         DO igroup = 0, ngroup - 1

            nfullcols_to_split = mp2_env%ri_rpa_im_time%sizes_array_cm_mao_virt(j_mem)

            itmp = get_limit(nfullcols_to_split, ngroup, igroup)

            IF (j_mem == 1) THEN

               mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(igroup) = itmp(1)
               mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(igroup) = itmp(2)
               mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_fullcol(igroup) = itmp(2) - itmp(1) + 1

               CALL get_blk_from_indx(indx=itmp(1), blk=blk_start, blk_offset=col_blk_offset, blk_sizes=col_blk_sizes_prim)
               CALL get_blk_from_indx(indx=itmp(2), blk=blk_end, blk_offset=col_blk_offset, blk_sizes=col_blk_sizes_prim)

               mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(igroup) = blk_start
               mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col(igroup) = blk_end
               mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(igroup) = blk_end - blk_start + 1

            ELSE

               offset_fullcol = mp2_env%ri_rpa_im_time_util(j_mem - 1)%ends_array_prim_fullcol(ngroup - 1)

               mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(igroup) = itmp(1) + offset_fullcol
               mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(igroup) = itmp(2) + offset_fullcol
               mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_fullcol(igroup) = itmp(2) - itmp(1) + 1

               CALL get_blk_from_indx(indx=(itmp(1) + offset_fullcol), blk=blk_start, &
                                      blk_offset=col_blk_offset, blk_sizes=col_blk_sizes_prim)
               CALL get_blk_from_indx(indx=(itmp(2) + offset_fullcol), blk=blk_end, &
                                      blk_offset=col_blk_offset, blk_sizes=col_blk_sizes_prim)

               mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(igroup) = blk_start
               mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_col(igroup) = blk_end
               mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(igroup) = blk_end - blk_start + 1

            END IF

         ENDDO

      END DO

      color_sub_row = mp2_env%ri_rpa_im_time_util(1)%color_sub_row
      color_sub_col = mp2_env%ri_rpa_im_time_util(1)%color_sub_col

      CALL dbcsr_get_info(mat_P_local%matrix, &
                          nblkrows_total=nblkrows_RI, &
                          row_blk_size=row_blk_sizes)

      DO i_mem = 1, cut_memory

         n_local_row_prim = mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(color_sub_row)
         row_start_local = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(color_sub_row)
         row_end_local = row_start_local + n_local_row_prim - 1

         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%start_row_data_block(row_start_local: &
                                                                           row_start_local + n_local_row_prim - 1))
         mp2_env%ri_rpa_im_time_util(i_mem)%start_row_data_block = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(i_mem)%end_row_data_block(row_start_local: &
                                                                         row_start_local + n_local_row_prim - 1))
         mp2_env%ri_rpa_im_time_util(i_mem)%end_row_data_block = 0

      END DO

      DO j_mem = 1, cut_memory

         n_local_col_prim = mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(color_sub_col)
         col_start_local = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(color_sub_col)
         col_end_local = col_start_local + n_local_col_prim - 1

         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%start_col_data_block(col_start_local: &
                                                                           col_start_local + n_local_col_prim - 1))
         mp2_env%ri_rpa_im_time_util(j_mem)%start_col_data_block = 0
         ALLOCATE (mp2_env%ri_rpa_im_time_util(j_mem)%end_col_data_block(col_start_local: &
                                                                         col_start_local + n_local_col_prim - 1))
         mp2_env%ri_rpa_im_time_util(j_mem)%end_col_data_block = 0

      END DO

      ALLOCATE (col_blk_sizes_cut_memory(cut_memory, cut_memory))
      col_blk_sizes_cut_memory(:, :) = 0

      ALLOCATE (mp2_env%ri_rpa_im_time_2d_util(cut_memory, cut_memory))

      DO i_mem = 1, cut_memory

         n_local_row_prim = mp2_env%ri_rpa_im_time_util(i_mem)%sizes_array_prim_row(color_sub_row)
         row_start_local = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_row(color_sub_row)
         row_end_local = row_start_local + n_local_row_prim - 1

         DO j_mem = 1, cut_memory

            n_local_col_prim = mp2_env%ri_rpa_im_time_util(j_mem)%sizes_array_prim_col(color_sub_col)
            col_start_local = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_col(color_sub_col)
            col_end_local = col_start_local + n_local_col_prim - 1

            ALLOCATE (mp2_env%ri_rpa_im_time_2d_util(i_mem, j_mem)%offset_combi_block( &
                      row_start_local:row_start_local + n_local_row_prim - 1, &
                      col_start_local:col_start_local + n_local_col_prim - 1))
            mp2_env%ri_rpa_im_time_2d_util(i_mem, j_mem)%offset_combi_block = 0

            DO icol = 1, n_local_row_prim*n_local_col_prim

               ref_row = (icol - 1)/n_local_col_prim + 1 + row_start_local - 1

               ref_col = MODULO(icol - 1, n_local_col_prim) + 1 + col_start_local - 1

               IF (ref_row == row_start_local) THEN
                  start_row = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(color_sub_row)
                  end_row = row_blk_offset(ref_row) + row_blk_sizes_prim(ref_row) - 1
                  size_row = end_row - start_row + 1
                  end_row_data_block = row_blk_sizes_prim(ref_row)
                  start_row_data_block = end_row_data_block - size_row + 1
               ELSE IF (ref_row == row_end_local) THEN
                  start_row = row_blk_offset(ref_row)
                  end_row = mp2_env%ri_rpa_im_time_util(i_mem)%ends_array_prim_fullrow(color_sub_row)
                  size_row = end_row - start_row + 1
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
                  size_row = end_row - start_row + 1
                  start_row_data_block = mp2_env%ri_rpa_im_time_util(i_mem)%starts_array_prim_fullrow(color_sub_row) - &
                                         row_blk_offset(ref_row) + 1
                  end_row_data_block = start_row_data_block + size_row - 1
               END IF

               mp2_env%ri_rpa_im_time_util(i_mem)%start_row_data_block(ref_row) = start_row_data_block
               mp2_env%ri_rpa_im_time_util(i_mem)%end_row_data_block(ref_row) = end_row_data_block

               IF (ref_col == col_start_local) THEN
                  start_col = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(color_sub_col)
                  end_col = col_blk_offset(ref_col) + col_blk_sizes_prim(ref_col) - 1
                  size_col = end_col - start_col + 1
                  end_col_data_block = col_blk_sizes_prim(ref_col)
                  start_col_data_block = end_col_data_block - size_col + 1
               ELSE IF (ref_col == col_end_local) THEN
                  start_col = col_blk_offset(ref_col)
                  end_col = mp2_env%ri_rpa_im_time_util(j_mem)%ends_array_prim_fullcol(color_sub_col)
                  size_col = end_col - start_col + 1
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
                  size_col = end_col - start_col + 1
                  start_col_data_block = mp2_env%ri_rpa_im_time_util(j_mem)%starts_array_prim_fullcol(color_sub_col) - &
                                         col_blk_offset(ref_col) + 1
                  end_col_data_block = start_col_data_block + size_col - 1
               END IF

               mp2_env%ri_rpa_im_time_util(j_mem)%start_col_data_block(ref_col) = start_col_data_block
               mp2_env%ri_rpa_im_time_util(j_mem)%end_col_data_block(ref_col) = end_col_data_block

               mp2_env%ri_rpa_im_time_2d_util(i_mem, j_mem)%offset_combi_block(ref_row, ref_col) = &
                  col_blk_sizes_cut_memory(i_mem, j_mem)

               col_blk_sizes_cut_memory(i_mem, j_mem) = col_blk_sizes_cut_memory(i_mem, j_mem) + size_row*size_col

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

      CALL release_neighbor_list_sets(sab_orb_sub)

      CALL create_mat_munu(mat_P_global, qs_env, mp2_env, para_env, dft_control, &
                           atomic_kind_set, qs_kind_set, atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_global, &
                           do_ri_aux_basis=.TRUE.)

      CALL dbcsr_reserve_all_blocks(mat_P_global%matrix)

      CALL cp_blacs_env_release(blacs_env_global)

      CALL release_neighbor_list_sets(sab_orb_sub)

      IF (do_mao) THEN
         CALL dbcsr_get_info(mat_dm_occ_global_mao%matrix, &
                             row_blk_size=blk_sizes_mao_occ)
         blk_sizes_occ => blk_sizes_mao_occ
         CALL dbcsr_get_info(mat_dm_virt_global_mao%matrix, &
                             row_blk_size=blk_sizes_mao_virt)
         blk_sizes_virt => blk_sizes_mao_virt
      ELSE
         CALL dbcsr_get_info(mat_munu%matrix, &
                             row_blk_size=row_blk_sizes)
         blk_sizes_occ => row_blk_sizes
         blk_sizes_virt => row_blk_sizes
      END IF

      ! create mat_dm_occ/virt_local_mao (in case not using MAOs, the size is the Gauss basis)
      NULLIFY (mat_dm_occ_local%matrix)
      ALLOCATE (mat_dm_occ_local%matrix)
      CALL dbcsr_create(matrix=mat_dm_occ_local%matrix, template=mat_munu%matrix, &
                        name="mat_dm_occ_local", &
                        row_blk_size=blk_sizes_occ, col_blk_size=blk_sizes_occ)

      NULLIFY (mat_dm_virt_local%matrix)
      ALLOCATE (mat_dm_virt_local%matrix)
      CALL dbcsr_create(matrix=mat_dm_virt_local%matrix, template=mat_munu%matrix, &
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

      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_dm_occ_global_mao, &
                                                            mat_munu_mao_occ_virt, &
                                                            mat_munu_mao_virt_occ, &
                                                            mat_dm_virt_global_mao
      TYPE(dbcsr_p_type), INTENT(IN)                     :: mat_munu
      LOGICAL, INTENT(IN)                                :: do_mao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coeff_occ, mao_coeff_virt, &
                                                            mao_coeff_occ_A, mao_coeff_virt_A, &
                                                            matrix_s
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_beta
      INTEGER, INTENT(IN)                                :: homo, homo_beta, nmo, nspins, unit_nr
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues, mo_eigenvalues_beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mao_basis_and_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_mao_occ, blk_sizes_mao_virt
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: dm_for_maos_virt, rho_ao_kp, &
                                                            scaled_dm_for_maos_occ
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      IF (do_mao) THEN

         ! get density matrix
         CALL get_qs_env(qs_env, rho=rho)
         NULLIFY (rho_ao_kp)
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

         IF (mp2_env%ri_rpa_im_time%opt_sc_dm_occ) THEN
            CALL build_scaled_dm_occ(scaled_dm_for_maos_occ, rho_ao_kp, mo_coeff, &
                                     mo_coeff_beta, nspins, homo, homo_beta, nmo, &
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
            CALL dbcsr_deallocate_matrix_set(scaled_dm_for_maos_occ)
         END IF

         ! the column has the MAO basis
         CALL dbcsr_get_info(mao_coeff_occ(1)%matrix, &
                             col_blk_size=blk_sizes_mao_occ)

         ! the column has the MAO basis
         CALL dbcsr_get_info(mao_coeff_virt(1)%matrix, &
                             col_blk_size=blk_sizes_mao_virt)

         CALL get_qs_env(qs_env=qs_env, dbcsr_dist=dbcsr_dist)

         NULLIFY (mat_dm_occ_global_mao%matrix)
         ALLOCATE (mat_dm_occ_global_mao%matrix)
         CALL dbcsr_create(matrix=mat_dm_occ_global_mao%matrix, &
                           name="mat_dm_occ_global_mao", &
                           dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=blk_sizes_mao_occ, col_blk_size=blk_sizes_mao_occ)

         NULLIFY (mat_dm_virt_global_mao%matrix)
         ALLOCATE (mat_dm_virt_global_mao%matrix)
         CALL dbcsr_create(matrix=mat_dm_virt_global_mao%matrix, &
                           name="mat_dm_virt_global_mao", &
                           dist=dbcsr_dist, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=blk_sizes_mao_virt, col_blk_size=blk_sizes_mao_virt)

         NULLIFY (mat_munu_mao_occ_virt%matrix)
         ALLOCATE (mat_munu_mao_occ_virt%matrix)
         CALL dbcsr_create(matrix=mat_munu_mao_occ_virt%matrix, template=mat_munu%matrix, &
                           name="mat_munu_mao_occ_virt", &
                           row_blk_size=blk_sizes_mao_occ, col_blk_size=blk_sizes_mao_virt)

         NULLIFY (mat_munu_mao_virt_occ%matrix)
         ALLOCATE (mat_munu_mao_virt_occ%matrix)
         CALL dbcsr_create(matrix=mat_munu_mao_virt_occ%matrix, template=mat_munu%matrix, &
                           name="mat_munu_mao_virt_occ", &
                           row_blk_size=blk_sizes_mao_virt, col_blk_size=blk_sizes_mao_occ)

      ELSE

         NULLIFY (mat_dm_occ_global_mao%matrix)
         ALLOCATE (mat_dm_occ_global_mao%matrix)
         CALL dbcsr_create(mat_dm_occ_global_mao%matrix, &
                           name="mat_dm_occ_global_mao", &
                           template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

         NULLIFY (mat_dm_virt_global_mao%matrix)
         ALLOCATE (mat_dm_virt_global_mao%matrix)
         CALL dbcsr_create(mat_dm_virt_global_mao%matrix, &
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
   SUBROUTINE build_scaled_dm_occ(scaled_dm_for_maos_occ, rho_ao_kp, mo_coeff, &
                                  mo_coeff_beta, nspins, homo, homo_beta, nmo, &
                                  mo_eigenvalues, mo_eigenvalues_beta)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: scaled_dm_for_maos_occ, rho_ao_kp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_beta
      INTEGER, INTENT(IN)                                :: nspins, homo, homo_beta, nmo
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

      e_lumo = mo_eigenvalues(homo + 1)

      DO jjB = 1, nrow_local
         DO iiB = 1, ncol_local

            i_global = col_indices(iiB)

            IF (i_global .LE. homo) THEN

               fm_mo_coeff_occ_scaled%local_data(jjB, iiB) = &
                  mo_coeff%local_data(jjB, iiB)/(e_lumo - mo_eigenvalues(i_global))

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

         e_lumo_beta = mo_eigenvalues_beta(homo_beta + 1)

         DO jjB = 1, nrow_local
            DO iiB = 1, ncol_local

               i_global = col_indices(iiB)

               IF (i_global .LE. homo) THEN

                  fm_mo_coeff_occ_scaled_beta%local_data(jjB, iiB) = &
                     mo_coeff_beta%local_data(jjB, iiB)/(e_lumo_beta - mo_eigenvalues_beta(i_global))

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
      CALL dbcsr_allocate_matrix_set(scaled_dm_for_maos_occ, nspins, 1)

      DO ispin = 1, nspins
         ALLOCATE (scaled_dm_for_maos_occ(ispin, 1)%matrix)
         CALL dbcsr_create(matrix=scaled_dm_for_maos_occ(ispin, 1)%matrix, &
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
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: id_mat
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      INTEGER, INTENT(IN)                                :: nspin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_and_set_identity_dbscr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin

      CALL timeset(routineN, handle)

      NULLIFY (id_mat)
      CALL dbcsr_allocate_matrix_set(id_mat, nspin)

      DO ispin = 1, nspin
         ALLOCATE (id_mat(ispin)%matrix)
         CALL dbcsr_create(matrix=id_mat(ispin)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           template=rho_ao_kp(1, 1)%matrix)
         CALL dbcsr_reserve_diag_blocks(matrix=id_mat(ispin)%matrix)
         CALL dbcsr_add_on_diag(id_mat(ispin)%matrix, 1.0_dp)
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
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coeff_A, mao_coeff, matrix_s
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mao_build_trafo_A', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: col, handle, ispin, row
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_A, block_B, block_S
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      NULLIFY (mao_coeff_A)
      CALL dbcsr_allocate_matrix_set(mao_coeff_A, nspins)

      DO ispin = 1, nspins
         ALLOCATE (mao_coeff_A(ispin)%matrix)
         CALL dbcsr_create(matrix=mao_coeff_A(ispin)%matrix, &
                           template=mao_coeff(ispin)%matrix)
         CALL dbcsr_reserve_diag_blocks(matrix=mao_coeff_A(ispin)%matrix)
      END DO

      DO ispin = 1, nspins

!$OMP PARALLEL DEFAULT(NONE) SHARED(mao_coeff, mao_coeff_A, matrix_s, ispin) &
!$OMP PRIVATE(iter,row,col,block_A,block_B,block_S,found)
         CALL dbcsr_iterator_start(iter, mao_coeff_A(ispin)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, row, col, block_A)
            CPASSERT(row == col)

            CALL dbcsr_get_block_p(matrix=mao_coeff(ispin)%matrix, row=row, col=col, block=block_B, found=found)
            CPASSERT(ASSOCIATED(block_B))

            CALL dbcsr_get_block_p(matrix=matrix_s(1)%matrix, row=row, col=col, block=block_S, found=found)
            CPASSERT(ASSOCIATED(block_S))

            block_A = MATMUL(block_S, block_B)

         ENDDO
         CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_P_global ...
!> \param mat_dm_occ_global_mao ...
!> \param mat_dm_virt_global_mao ...
!> \param para_env_sub_im_time_3c ...
! **************************************************************************************************
   SUBROUTINE clean_up_im_time_t(mat_P_global, mat_dm_occ_global_mao, mat_dm_virt_global_mao, &
                                 para_env_sub_im_time_3c)
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_P_global, mat_dm_occ_global_mao, &
                                                            mat_dm_virt_global_mao
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub_im_time_3c

      CALL dbcsr_release(mat_P_global%matrix)
      DEALLOCATE (mat_P_global%matrix)

      ! we are just releasing due to lecacy mao code, this should go away
      CALL dbcsr_release(mat_dm_occ_global_mao%matrix)
      DEALLOCATE (mat_dm_occ_global_mao%matrix)
      CALL dbcsr_release(mat_dm_virt_global_mao%matrix)
      DEALLOCATE (mat_dm_virt_global_mao%matrix)

      CALL cp_para_env_release(para_env_sub_im_time_3c)

   END SUBROUTINE clean_up_im_time_t

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
                               mat_dm_occ_local, mat_dm_virt_local, &
                               para_env_sub_im_time_3c, para_env_sub_im_time_P, mao_coeff_occ, &
                               mao_coeff_virt, mao_coeff_occ_A, mao_coeff_virt_A, mp2_env)

      TYPE(dbcsr_p_type), INTENT(INOUT) :: mat_munu, mat_P_local, mat_P_global, mat_M, &
         mat_dm_occ_global_mao, mat_dm_virt_global_mao, mat_munu_mao_occ_virt, &
         mat_munu_mao_virt_occ, mat_dm_occ_local, mat_dm_virt_local
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub_im_time_3c, &
                                                            para_env_sub_im_time_P
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coeff_occ, mao_coeff_virt, &
                                                            mao_coeff_occ_A, mao_coeff_virt_A
      TYPE(mp2_type), POINTER                            :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'clean_up_im_time', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: do_mao

      CALL timeset(routineN, handle)

      do_mao = mp2_env%ri_rpa_im_time%do_mao

      CALL dbcsr_release(mat_munu%matrix)
      DEALLOCATE (mat_munu%matrix)

      CALL dbcsr_release(mat_P_local%matrix)
      DEALLOCATE (mat_P_local%matrix)

      CALL dbcsr_release(mat_P_global%matrix)
      DEALLOCATE (mat_P_global%matrix)

      CALL dbcsr_release(mat_M%matrix)
      DEALLOCATE (mat_M%matrix)

      CALL dbcsr_release(mat_dm_occ_local%matrix)
      DEALLOCATE (mat_dm_occ_local%matrix)

      CALL dbcsr_release(mat_dm_virt_local%matrix)
      DEALLOCATE (mat_dm_virt_local%matrix)

      IF (do_mao) THEN
         CALL dbcsr_deallocate_matrix_set(mao_coeff_occ)
         CALL dbcsr_deallocate_matrix_set(mao_coeff_virt)
         CALL dbcsr_deallocate_matrix_set(mao_coeff_occ_A)
         CALL dbcsr_deallocate_matrix_set(mao_coeff_virt_A)
         CALL dbcsr_release(mat_munu_mao_occ_virt%matrix)
         DEALLOCATE (mat_munu_mao_occ_virt%matrix)
         CALL dbcsr_release(mat_munu_mao_virt_occ%matrix)
         DEALLOCATE (mat_munu_mao_virt_occ%matrix)
      END IF

      CALL dbcsr_release(mat_dm_occ_global_mao%matrix)
      DEALLOCATE (mat_dm_occ_global_mao%matrix)
      CALL dbcsr_release(mat_dm_virt_global_mao%matrix)
      DEALLOCATE (mat_dm_virt_global_mao%matrix)

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

      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_M
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub_P
      INTEGER, INTENT(IN)                                :: nblkrows_total
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes, col_blk_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: mepos_P_from_RI_row

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_mat_M', routineP = moduleN//':'//routineN

      INTEGER                                            :: col_size, handle, icol, irow, &
                                                            nblkcols_total, row_size
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, row_dist
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL timeset(routineN, handle)

      CALL get_blacs_info(blacs_env_sub_P, para_env=para_env, blacs2mpi=pgrid)

      ! just round robin for row_dist_data and col_dist_data
      row_size = SIZE(pgrid, 1)
      col_size = SIZE(pgrid, 2)
      nblkcols_total = 1
      ALLOCATE (row_dist(nblkrows_total), col_dist(nblkcols_total))
      ALLOCATE (mepos_P_from_RI_row(nblkrows_total))

      DO irow = 1, nblkrows_total
         row_dist(irow) = MODULO(irow, row_size)
         mepos_P_from_RI_row(irow) = MODULO(irow, row_size)
      END DO

      DO icol = 1, nblkcols_total
         col_dist(icol) = MODULO(icol, col_size)
      END DO

      CALL dbcsr_distribution_new(dist, group=para_env%group, pgrid=pgrid, &
                                  row_dist=row_dist, col_dist=col_dist)

      NULLIFY (mat_M%matrix)
      ALLOCATE (mat_M%matrix)

      CALL dbcsr_create(matrix=mat_M%matrix, &
                        name="M_P_alphadelta", &
                        dist=dist, matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
                        nze=0)

      CALL dbcsr_distribution_release(dist)

      DEALLOCATE (row_dist, col_dist)

      CALL timestop(handle)

   END SUBROUTINE create_mat_M

! **************************************************************************************************
!> \brief ...
!> \param num_pe ...
!> \param n_group_row ...
!> \param n_group_col ...
!> \author Jan Wilhelm
! **************************************************************************************************
   PURE SUBROUTINE generate_integer_product(num_pe, n_group_row, n_group_col)

      INTEGER, INTENT(IN)                                :: num_pe
      INTEGER, INTENT(OUT)                               :: n_group_row, n_group_col

      CHARACTER(LEN=*), PARAMETER :: routineN = 'generate_integer_product', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: a_int, b_int, num_pe_temp_int, sqrt_int, &
                                                            square_int, x_int
      LOGICAL                                            :: stay_while
      REAL(KIND=dp)                                      :: num_pe_real, offset_real, square_real

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

      num_pe_real = REAL(num_pe_temp_int, KIND=dp)

      offset_real = 0.0_dp

      stay_while = .TRUE.

      DO WHILE (stay_while)

         square_real = (CEILING(SQRT(num_pe_real)) + offset_real)**2 - num_pe_real

         square_int = NINT(square_real)

         sqrt_int = NINT(SQRT(square_real))

         IF (sqrt_int**2 == square_int) THEN

            stay_while = .FALSE.

         ELSE

            offset_real = offset_real + 1.0_dp

         END IF

      END DO

      x_int = NINT(CEILING(SQRT(num_pe_real)) + offset_real)

      n_group_row = (x_int + sqrt_int)*b_int
      n_group_col = (x_int - sqrt_int)*a_int

      ! additional balancing
      IF (n_group_row == 2*(n_group_row/2) .AND. n_group_row > 2*n_group_col) THEN

         n_group_row = n_group_row/2
         n_group_col = n_group_col*2

      END IF

   END SUBROUTINE generate_integer_product

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
   PURE SUBROUTINE get_start_end_size_indx(start_indx, end_indx, size_indx, nblkrows_total, itmp_1, itmp_2, &
                                           row_blk_offset, row_blk_sizes)
      INTEGER, INTENT(INOUT)                             :: start_indx, end_indx, size_indx
      INTEGER, INTENT(IN)                                :: nblkrows_total, itmp_1, itmp_2
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

            end_indx = row_blk_offset(blk) + row_blk_sizes(blk) - 1

         END IF

      END DO

      size_indx = end_indx - start_indx + 1

      ! have a check that if there is nothing to be done for the specific memory_cut, then we know it
      IF (start_indx == 0 .AND. end_indx == 0) THEN
         size_indx = 0
      END IF

   END SUBROUTINE get_start_end_size_indx

! **************************************************************************************************
!> \brief ...
!> \param indx ...
!> \param blk ...
!> \param blk_offset ...
!> \param blk_sizes ...
!> \author Jan Wilhelm
! **************************************************************************************************
   PURE SUBROUTINE get_blk_from_indx(indx, blk, blk_offset, blk_sizes)

      INTEGER, INTENT(IN)                                :: indx
      INTEGER, INTENT(OUT)                               :: blk
      INTEGER, DIMENSION(:), POINTER                     :: blk_offset, blk_sizes

      INTEGER                                            :: iblk, nblkrows_total

      nblkrows_total = SIZE(blk_sizes)

      DO iblk = 1, nblkrows_total

         IF (blk_offset(iblk) <= indx .AND. blk_offset(iblk) + blk_sizes(iblk) - 1 >= indx) THEN

            blk = iblk

         END IF

      END DO

   END SUBROUTINE get_blk_from_indx

! **************************************************************************************************
!> \brief ...
!> \param dft_control ...
!> \param eps_pgf_orb_old ...
!> \param eps_rho_rspace_old ...
!> \param eps_gvg_rspace_old ...
! **************************************************************************************************
   PURE SUBROUTINE get_eps_old(dft_control, eps_pgf_orb_old, eps_rho_rspace_old, eps_gvg_rspace_old)

      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(kind=dp), INTENT(OUT)                         :: eps_pgf_orb_old, eps_rho_rspace_old, &
                                                            eps_gvg_rspace_old

      ! re-init the radii to be able to generate pair lists with MP2-appropriate screening
      eps_pgf_orb_old = dft_control%qs_control%eps_pgf_orb
      eps_rho_rspace_old = dft_control%qs_control%eps_rho_rspace
      eps_gvg_rspace_old = dft_control%qs_control%eps_gvg_rspace

   END SUBROUTINE get_eps_old

END MODULE mp2_gpw
