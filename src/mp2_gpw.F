!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
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
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_blacs_env,                    ONLY: BLACS_GRID_SQUARE,&
                                              cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist,&
                                              cp_dbcsr_m_by_n_from_row_template
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_type
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
        dbcsr_clear_mempools, dbcsr_copy, dbcsr_create, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_filter, dbcsr_init_p, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_release, dbcsr_reserve_all_blocks, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_type_real_default, dbcsr_type_symmetric
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
   USE input_constants,                 ONLY: &
        do_eri_gpw, do_eri_mme, do_eri_os, do_potential_coulomb, do_potential_id, &
        do_potential_long, do_potential_truncated, eri_default, mp2_method_gpw, ri_default, &
        ri_mp2_method_gpw
   USE input_section_types,             ONLY: section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE libint_wrapper,                  ONLY: cp_libint_static_cleanup,&
                                              cp_libint_static_init
   USE machine,                         ONLY: default_output_unit,&
                                              m_flush
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
   USE rpa_main,                        ONLY: rpa_ri_compute_en
   USE rpa_rse,                         ONLY: rse_energy
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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mp2_gpw_main'

      INTEGER :: blacs_grid_layout, color_sub, comm_sub, dimen, dimen_RI, dimen_RI_red, &
         eri_method, gw_corr_lev_occ, gw_corr_lev_occ_beta, gw_corr_lev_virt, &
         gw_corr_lev_virt_beta, handle, homo, homo_beta, local_unit_nr, my_group_L_end, &
         my_group_L_size, my_group_L_start, nelectron, nelectron_beta, nmo, nspins, &
         potential_type, ri_metric_type
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ends_array_mc, ends_array_mc_block, &
                                                            starts_array_mc, starts_array_mc_block
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL :: blacs_repeatable, do_bse, do_gpw, do_im_time, do_kpoints_cubic_RPA, my_do_gw, &
         my_do_ri_mp2, my_do_ri_rpa, my_do_ri_sos_laplace_mp2
      REAL(KIND=dp) :: Emp2_AB, Emp2_BB, Emp2_Cou_BB, Emp2_d2_AB, Emp2_d_AB, Emp2_EX_BB, &
         eps_gvg_rspace_old, eps_pgf_orb_old, eps_rho_rspace_old
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Eigenval, Eigenval_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C, BIb_C_beta, BIb_C_bse_ab, &
                                                            BIb_C_bse_ij, BIb_C_gw, BIb_C_gw_beta
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues, mo_eigenvalues_beta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_sub, blacs_env_sub_mat_munu
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: fm_matrix_L_RI_metric
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_beta
      TYPE(cp_logger_type), POINTER                      :: logger, logger_sub
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub, para_env_sub_RPA
      TYPE(dbcsr_p_type)                                 :: mat_munu, mat_P_global
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp
      TYPE(dbcsr_t_type)                                 :: t_3c_M
      TYPE(dbcsr_t_type), ALLOCATABLE, DIMENSION(:, :)   :: t_3c_O
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
         POINTER                                         :: sab_orb_sub
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
      do_bse = qs_env%mp2_env%ri_g0w0%do_bse
      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      IF (do_kpoints_cubic_RPA .AND. mp2_env%ri_rpa%do_ri_g0w0) THEN
         CPABORT("Full RPA k-points (DO_KPOINTS in LOW_SCALING section) not implemented with GW")
      ENDIF

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

      CALL get_cell(cell=cell, periodic=periodic)

      IF (do_im_time) THEN

         IF (.NOT. mp2_env%ri_rpa_im_time%group_size_internal) THEN
            ASSOCIATE (group_size_3c=>mp2_env%ri_rpa_im_time%group_size_3c, group_size_P=>mp2_env%ri_rpa_im_time%group_size_P)
               DO WHILE (MOD(para_env%num_pe, group_size_3c) .NE. 0)
                  group_size_3c = group_size_3c - 1
               END DO

               DO WHILE (MOD(para_env%num_pe, group_size_P) .NE. 0)
                  group_size_P = group_size_P - 1
               END DO
            END ASSOCIATE
         ENDIF

         IF (mp2_env%ri_metric%potential_type == ri_default) THEN
            IF (SUM(periodic) > 0) THEN
               mp2_env%ri_metric%potential_type = do_potential_id
            ELSE
               mp2_env%ri_metric%potential_type = do_potential_truncated
            ENDIF
         ENDIF

         ! statically initialize libint
         CALL cp_libint_static_init()

      END IF

      IF (mp2_env%ri_metric%potential_type == ri_default) THEN
         mp2_env%ri_metric%potential_type = do_potential_coulomb
      END IF

      IF (mp2_env%eri_method == eri_default) THEN
         IF (SUM(periodic) > 0) mp2_env%eri_method = do_eri_gpw
         IF (SUM(periodic) == 0) mp2_env%eri_method = do_eri_os
         IF (SUM(mp2_env%ri_rpa_im_time%kp_grid) > 0) mp2_env%eri_method = do_eri_os
         IF (mp2_env%method == mp2_method_gpw) mp2_env%eri_method = do_eri_gpw
         IF (mp2_env%method == ri_mp2_method_gpw) mp2_env%eri_method = do_eri_gpw
         IF (mp2_env%ri_rpa_im_time%do_im_time_kpoints) mp2_env%eri_method = do_eri_os
      END IF
      eri_method = mp2_env%eri_method

      IF (unit_nr > 0 .AND. mp2_env%eri_method == do_eri_gpw) THEN
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

      matrix_s(1:1) => matrix_s_kp(1:1, 1)

      CALL get_eps_old(dft_control, eps_pgf_orb_old, eps_rho_rspace_old, eps_gvg_rspace_old)

      CALL create_mat_munu(mat_munu, qs_env, mp2_env, dft_control, atomic_kind_set, qs_kind_set, &
                           atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub_mat_munu, &
                           do_im_time=do_im_time)

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
               mat_munu, sab_orb_sub, &
               mo_coeff_o, mo_coeff_v, mo_coeff_all, mo_coeff_gw, &
               mp2_env%mp2_gpw%eps_filter, unit_nr, &
               mp2_env%mp2_memory, mp2_env%calc_PQ_cond_num, calc_forces, blacs_env_sub, my_do_gw .AND. .NOT. do_im_time, &
               do_bse, gd_B_all, starts_array_mc, ends_array_mc, starts_array_mc_block, ends_array_mc_block, &
               gw_corr_lev_occ, gw_corr_lev_virt, &
               do_im_time, do_kpoints_cubic_RPA, kpoints, &
               t_3c_M, t_3c_O, &
               mp2_env%ri_metric, &
               gd_B_occ_bse, gd_B_virt_bse, &
               BIb_C_beta, BIb_C_gw_beta, gd_B_virtual_beta, &
               homo_beta, mo_coeff_o_beta, mo_coeff_v_beta, mo_coeff_all_beta, mo_coeff_gw_beta)
         ELSE
            ! closed shell case (RI)
            CALL mp2_ri_gpw_compute_in(BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, gd_array, gd_B_virtual, &
                                       dimen_RI, dimen_RI_red, qs_env, para_env, para_env_sub, &
                                       color_sub, cell, particle_set, &
                                       atomic_kind_set, qs_kind_set, mo_coeff, fm_matrix_L_RI_metric, nmo, homo, &
                                       mat_munu, sab_orb_sub, &
                                       mo_coeff_o, mo_coeff_v, mo_coeff_all, mo_coeff_gw, &
                                       mp2_env%mp2_gpw%eps_filter, unit_nr, &
                                       mp2_env%mp2_memory, mp2_env%calc_PQ_cond_num, calc_forces, &
                                       blacs_env_sub, my_do_gw .AND. .NOT. do_im_time, do_bse, gd_B_all, &
                                       starts_array_mc, ends_array_mc, &
                                       starts_array_mc_block, ends_array_mc_block, &
                                       gw_corr_lev_occ, gw_corr_lev_virt, &
                                       do_im_time, do_kpoints_cubic_RPA, kpoints, &
                                       t_3c_M, t_3c_O, &
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

         para_env_sub_RPA => para_env_sub

         IF (do_im_time) THEN

            CALL create_matrix_P(mat_P_global, qs_env, mp2_env, para_env, dft_control, atomic_kind_set, &
                                 qs_kind_set, atom2d, molecule_kind_set, molecule_set, particle_set, cell)

         END IF

         ! RI-RPA
         IF (nspins == 2) THEN
            CALL rpa_ri_compute_en(qs_env, Emp2, mp2_env, BIb_C, BIb_C_gw, BIb_C_bse_ij, BIb_C_bse_ab, &
                                   para_env, para_env_sub_RPA, color_sub, &
                                   gd_array, gd_B_virtual, gd_B_all, gd_B_occ_bse, gd_B_virt_bse, &
                                   mo_coeff, fm_matrix_L_RI_metric, kpoints, &
                                   Eigenval, nmo, homo, dimen_RI, dimen_RI_red, gw_corr_lev_occ, gw_corr_lev_virt, &
                                   unit_nr, my_do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_bse, matrix_s, &
                                   mat_munu, &
                                   mat_P_global, &
                                   t_3c_M, t_3c_O, &
                                   starts_array_mc, ends_array_mc, &
                                   starts_array_mc_block, ends_array_mc_block, &
                                   BIb_C_beta, homo_beta, Eigenval_beta, &
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
                                   unit_nr, my_do_ri_sos_laplace_mp2, my_do_gw, do_im_time, do_bse, matrix_s, &
                                   mat_munu, &
                                   mat_P_global, &
                                   t_3c_M, t_3c_O, &
                                   starts_array_mc, ends_array_mc, &
                                   starts_array_mc_block, ends_array_mc_block)

            IF (mp2_env%ri_rpa%do_rse) &
               CALL rse_energy(qs_env, mp2_env, para_env, dft_control, &
                               mo_coeff, nmo, homo, Eigenval)
         END IF

         IF (do_im_time) THEN
            CALL cleanup_imtime(mat_P_global)

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'replicate_mat_to_subgroup'

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_mo_coeff_v_o_all'

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
!> \param do_im_time ...
!> \author Jan Wilhelm, code by Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE create_mat_munu(mat_munu, qs_env, mp2_env, dft_control, atomic_kind_set, qs_kind_set, &
                              atom2d, molecule_kind_set, &
                              molecule_set, sab_orb_sub, particle_set, cell, blacs_env_sub, &
                              do_ri_aux_basis, do_mixed_basis, group_size_prim, do_im_time)

      TYPE(dbcsr_p_type), INTENT(OUT)                    :: mat_munu
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
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
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_im_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_mat_munu'

      INTEGER                                            :: handle, ikind, natom, nkind
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes, row_blk_sizes
      LOGICAL                                            :: do_kpoints, my_do_im_time, &
                                                            my_do_mixed_basis, my_do_ri_aux_basis
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: orb_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      REAL(kind=dp)                                      :: subcells
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist_sub
      TYPE(distribution_1d_type), POINTER                :: local_molecules_sub, local_particles_sub
      TYPE(distribution_2d_type), POINTER                :: distribution_2d_sub
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ri_aux
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

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

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_matrix_P'

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

      CALL create_mat_munu(mat_P_global, qs_env, mp2_env, dft_control, &
                           atomic_kind_set, qs_kind_set, atom2d, molecule_kind_set, &
                           molecule_set, sab_orb_sub, particle_set, cell, blacs_env_global, &
                           do_ri_aux_basis=.TRUE.)

      CALL dbcsr_reserve_all_blocks(mat_P_global%matrix)
      CALL cp_blacs_env_release(blacs_env_global)

      CALL release_neighbor_list_sets(sab_orb_sub)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_P_global ...
! **************************************************************************************************
   SUBROUTINE cleanup_imtime(mat_P_global)
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: mat_P_global

      CALL dbcsr_release(mat_P_global%matrix)
      DEALLOCATE (mat_P_global%matrix)

   END SUBROUTINE cleanup_imtime

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
