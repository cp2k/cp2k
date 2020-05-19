!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate MP2 energy
!> \par History
!>      05.2011 created [Mauro Del Ben]
!> \author Mauro Del Ben
! **************************************************************************************************
MODULE mp2
   USE admm_methods,                    ONLY: admm_mo_merge_ks_matrix
   USE admm_types,                      ONLY: admm_type
   USE admm_utils,                      ONLY: admm_correct_for_eigenvalues,&
                                              admm_uncorrect_for_eigenvalues
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE bibliography,                    ONLY: DelBen2012,&
                                              DelBen2015b,&
                                              Rybkin2016,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_get_diag, dbcsr_get_info, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_release_p, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_symmetric
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_types,                       ONLY: &
        alloc_containers, dealloc_containers, hfx_basis_info_type, hfx_basis_type, &
        hfx_container_type, hfx_create_basis_types, hfx_init_container, hfx_release_basis_types, &
        hfx_type
   USE input_constants,                 ONLY: &
        cholesky_inverse, cholesky_off, do_admm_basis_projection, do_admm_purify_none, do_eri_gpw, &
        do_eri_mme, do_potential_TShPSC, gw_print_exx, gw_read_exx, xc_none
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get,&
                                              section_vals_val_set
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_memory,&
                                              m_walltime
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_comm_split_direct,&
                                              mp_max,&
                                              mp_sum,&
                                              mp_sync
   USE mp2_direct_method,               ONLY: mp2_canonical_direct_single_batch
   USE mp2_gpw,                         ONLY: mp2_gpw_main
   USE mp2_optimize_ri_basis,           ONLY: optimize_ri_basis_main
   USE mp2_types,                       ONLY: mp2_biel_type,&
                                              mp2_method_direct,&
                                              mp2_method_gpw,&
                                              mp2_ri_optimize_basis,&
                                              mp2_type,&
                                              ri_mp2_laplace,&
                                              ri_mp2_method_gpw,&
                                              ri_rpa_method_gpw
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_matrix_pools,                 ONLY: mpools_create,&
                                              mpools_rebuild_fm_pools,&
                                              mpools_release,&
                                              qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: eigensolver,&
                                              eigensolver_symm
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE scf_control_types,               ONLY: scf_control_type
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2'

   PUBLIC :: mp2_main

CONTAINS

! **************************************************************************************************
!> \brief the main entry point for MP2 calculations
!> \param qs_env ...
!> \param calc_forces ...
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE mp2_main(qs_env, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'mp2_main', routineP = moduleN//':'//routineN

      INTEGER :: bin, cholesky_method, dimen, handle, handle2, i, i_thread, iatom, ikind, irep, &
         ispin, max_nset, my_bin_size, n_rep_hf, n_threads, nao, natom, ncol_block, ndep, &
         nelec_alpha, nelec_beta, nelectron, nfullcols_total, nfullrows_total, nkind, nrow_block, &
         nspins, unit_nr
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      LOGICAL :: calc_ex, do_admm, do_admm_rpa_exx, do_dynamic_load_balancing, do_exx, do_gw, &
         do_kpoints_cubic_RPA, free_hfx_buffer
      REAL(KIND=dp) :: E_ex_from_GW, Emp2, Emp2_AA, Emp2_AA_Cou, Emp2_AA_ex, Emp2_AB, Emp2_AB_Cou, &
         Emp2_AB_ex, Emp2_BB, Emp2_BB_Cou, Emp2_BB_ex, Emp2_Cou, Emp2_ex, Emp2_S, Emp2_T, maxocc, &
         mem_real, t1, t2, t3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Auto, Auto_alpha, Auto_beta
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: C, C_alpha, C_beta
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_matrix_ks, fm_matrix_s, &
                                                            fm_matrix_work, mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_transl, matrix_s_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_basis_info_type)                          :: basis_info
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      TYPE(hfx_container_type), DIMENSION(:), POINTER    :: integral_containers
      TYPE(hfx_container_type), POINTER                  :: maxval_container
      TYPE(hfx_type), POINTER                            :: actual_x_data
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos, mos_mp2
      TYPE(mp2_biel_type)                                :: mp2_biel
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_matrix_pools_type), POINTER                :: my_mpools
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input
      TYPE(virial_type), POINTER                         :: virial

      ! If SCF has not converged we should abort MP2 calculation
      IF (qs_env%mp2_env%hf_fail) THEN
         CALL cp_abort(__LOCATION__, "SCF not converged: "// &
                       "not possible to run MP2")
      ENDIF

      NULLIFY (virial, dft_control, blacs_env, kpoints)
      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      CALL cite_reference(DelBen2012)

      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      ! for cubic RPA and GW, we have kpoints and therefore, we get other matrices later
      IF (do_kpoints_cubic_RPA) THEN

         CALL get_qs_env(qs_env, &
                         input=input, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         dft_control=dft_control, &
                         particle_set=particle_set, &
                         para_env=para_env, &
                         blacs_env=blacs_env, &
                         energy=energy, &
                         kpoints=kpoints, &
                         scf_env=scf_env, &
                         scf_control=scf_control, &
                         matrix_ks_kp=matrix_ks_transl, &
                         matrix_s_kp=matrix_s_kp, &
                         mp2_env=mp2_env)

         CALL get_gamma(matrix_s, matrix_ks, mos, &
                        matrix_s_kp, matrix_ks_transl, kpoints)

      ELSE

         CALL get_qs_env(qs_env, &
                         input=input, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         dft_control=dft_control, &
                         particle_set=particle_set, &
                         para_env=para_env, &
                         blacs_env=blacs_env, &
                         energy=energy, &
                         mos=mos, &
                         scf_env=scf_env, &
                         scf_control=scf_control, &
                         virial=virial, &
                         matrix_ks=matrix_ks, &
                         matrix_s=matrix_s, &
                         matrix_ks_aux_fit=matrix_ks_aux, &
                         mp2_env=mp2_env, &
                         admm_env=admm_env)

      END IF

      unit_nr = cp_print_key_unit_nr(logger, input, "DFT%XC%WF_CORRELATION%MP2%MP2_INFO", &
                                     extension=".mp2Log")

      IF (unit_nr > 0) THEN
         IF (mp2_env%method .NE. ri_rpa_method_gpw) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T2,A)') 'MP2 section'
            WRITE (unit_nr, '(T2,A)') '-----------'
            WRITE (unit_nr, *)
         ELSE
            WRITE (unit_nr, *)
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T2,A)') 'RI-RPA section'
            WRITE (unit_nr, '(T2,A)') '--------------'
            WRITE (unit_nr, *)
         END IF
      ENDIF

      IF (calc_forces) THEN
         CALL cite_reference(DelBen2015b)
         CALL cite_reference(Rybkin2016)
         IF (mp2_env%method /= ri_mp2_method_gpw) THEN
            CPABORT("No forces/gradients for the selected method.")
         END IF
      ENDIF

      IF (.NOT. do_kpoints_cubic_RPA) THEN
         IF (virial%pv_availability .AND. (.NOT. virial%pv_numer) .AND. mp2_env%eri_method == do_eri_mme) THEN
            CPABORT("analytical stress not implemented with ERI_METHOD = MME")
         ENDIF
      ENDIF

      IF (mp2_env%do_im_time .AND. mp2_env%eri_method .NE. do_eri_gpw) THEN
         mp2_env%mp2_num_proc = 1
      END IF

      IF (mp2_env%mp2_num_proc <= 0 .OR. mp2_env%mp2_num_proc > para_env%num_pe .OR. &
          MOD(para_env%num_pe, mp2_env%mp2_num_proc) .NE. 0) THEN
         IF (unit_nr > 0 .AND. mp2_env%mp2_num_proc .NE. -1) &
            WRITE (unit_nr, '(T3,A,T76,I5)') 'Requested number of processes per group:', mp2_env%mp2_num_proc
         mp2_env%mp2_num_proc = para_env%num_pe
      ENDIF
      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T76,I5)') 'Used number of processes per group:', mp2_env%mp2_num_proc
      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'Maximum allowed memory usage per MPI process:', &
         mp2_env%mp2_memory, ' MiB'

      IF ((mp2_env%method .NE. mp2_method_gpw) .AND. &
          (mp2_env%method .NE. ri_mp2_method_gpw) .AND. &
          (mp2_env%method .NE. ri_rpa_method_gpw) .AND. &
          (mp2_env%method .NE. ri_mp2_laplace)) THEN
         CALL m_memory(mem)
         mem_real = (mem + 1024*1024 - 1)/(1024*1024)
         CALL mp_max(mem_real, para_env%group)
         mp2_env%mp2_memory = mp2_env%mp2_memory - mem_real
         IF (mp2_env%mp2_memory < 0.0_dp) mp2_env%mp2_memory = 1.0_dp

         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'Available memory per MPI process for MP2:', &
            mp2_env%mp2_memory, ' MiB'
      END IF

      IF (unit_nr > 0) CALL m_flush(unit_nr)

      nspins = dft_control%nspins
      natom = SIZE(particle_set, 1)

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)
      nkind = SIZE(atomic_kind_set, 1)

      do_admm_rpa_exx = mp2_env%ri_rpa%do_admm
      CALL hfx_create_basis_types(basis_parameter, basis_info, qs_kind_set, do_admm_rpa_exx)

      dimen = 0
      max_nset = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         dimen = dimen + SUM(basis_parameter(ikind)%nsgf)
         max_nset = MAX(max_nset, basis_parameter(ikind)%nset)
      END DO

      CALL get_mo_set(mo_set=mos(1)%mo_set, nao=nao)
!      CPASSERT(dimen == nao)

      ! diagonalize the KS matrix in order to have the full set of MO's
      ! get S and KS matrices in fm_type (create also a working array)
      NULLIFY (fm_matrix_s)
      NULLIFY (fm_matrix_ks)
      NULLIFY (fm_matrix_work)
      NULLIFY (fm_struct)
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nfullrows_total, nfullcols_total=nfullcols_total)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nfullrows_total, &
                               ncol_global=nfullcols_total, para_env=para_env)
      CALL cp_fm_create(fm_matrix_s, fm_struct, name="fm_matrix_s")
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, fm_matrix_s)

      CALL cp_fm_create(fm_matrix_ks, fm_struct, name="fm_matrix_ks")

      CALL cp_fm_create(fm_matrix_work, fm_struct, name="fm_matrix_work")
      CALL cp_fm_set_all(matrix=fm_matrix_work, alpha=0.0_dp)

      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(matrix=fm_matrix_ks, nrow_block=nrow_block, ncol_block=ncol_block)

      IF (scf_env%cholesky_method == cholesky_off) THEN
         CALL cp_fm_power(fm_matrix_s, fm_matrix_work, -0.5_dp, scf_control%eps_eigval, ndep)
         cholesky_method = cholesky_off
      ELSE
         ! calculate S^(-1/2) (cholesky decomposition)
         CALL cp_fm_cholesky_decompose(fm_matrix_s)
         CALL cp_fm_triangular_invert(fm_matrix_s)
         cholesky_method = cholesky_inverse
      END IF

      NULLIFY (mos_mp2)
      ALLOCATE (mos_mp2(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, maxocc=maxocc, nelectron=nelectron)

         NULLIFY (mos_mp2(ispin)%mo_set)
         CALL allocate_mo_set(mo_set=mos_mp2(ispin)%mo_set, &
                              nao=nao, &
                              nmo=nao, &
                              nelectron=nelectron, &
                              n_el_f=REAL(nelectron, dp), &
                              maxocc=maxocc, &
                              flexible_electron_count=dft_control%relax_multiplicity)
      END DO

      NULLIFY (my_mpools)
      CALL mpools_create(mpools=my_mpools)
      CALL mpools_rebuild_fm_pools(mpools=my_mpools, &
                                   mos=mos_mp2, &
                                   blacs_env=blacs_env, &
                                   para_env=para_env)

      DO ispin = 1, nspins

         ! If ADMM we should make the ks matrix up-to-date
         IF (dft_control%do_admm) THEN
            CALL admm_correct_for_eigenvalues(ispin, admm_env, matrix_ks(ispin)%matrix)
         END IF

         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, fm_matrix_ks)

         IF (dft_control%do_admm) THEN
            CALL admm_uncorrect_for_eigenvalues(ispin, admm_env, matrix_ks(ispin)%matrix)
         END IF

         CALL init_mo_set(mos_mp2(ispin)%mo_set, &
                          fm_pool=my_mpools%ao_mo_fm_pools(ispin)%pool, &
                          name="mp2_mos")

         IF (cholesky_method == cholesky_inverse) THEN

            ! diagonalize KS matrix
            CALL eigensolver(matrix_ks_fm=fm_matrix_ks, &
                             mo_set=mos_mp2(ispin)%mo_set, &
                             ortho=fm_matrix_s, &
                             work=fm_matrix_work, &
                             cholesky_method=cholesky_method, &
                             do_level_shift=.FALSE., &
                             level_shift=0.0_dp, &
                             matrix_u_fm=null(), &
                             use_jacobi=.FALSE.)

         ELSE IF (cholesky_method == cholesky_off) THEN

            CALL eigensolver_symm(matrix_ks_fm=fm_matrix_ks, &
                                  mo_set=mos_mp2(ispin)%mo_set, &
                                  ortho=fm_matrix_s, &
                                  work=fm_matrix_work, &
                                  do_level_shift=.FALSE., &
                                  level_shift=0.0_dp, &
                                  matrix_u_fm=null(), &
                                  use_jacobi=.FALSE., &
                                  jacobi_threshold=0.0_dp)

         END IF

      END DO

      CALL cp_fm_release(fm_matrix_s)
      CALL cp_fm_release(fm_matrix_ks)
      CALL cp_fm_release(fm_matrix_work)
      CALL mpools_release(mpools=my_mpools)

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")

      !   build the table of index
      t1 = m_walltime()
      ALLOCATE (mp2_biel%index_table(natom, max_nset))

      CALL build_index_table(natom, max_nset, mp2_biel%index_table, basis_parameter, kind_of)

      ! free the hfx_container (for now if forces are required we don't release the HFX stuff)
      free_hfx_buffer = .FALSE.
      IF (ASSOCIATED(qs_env%x_data)) THEN
         free_hfx_buffer = .TRUE.
         IF (calc_forces .AND. (.NOT. mp2_env%ri_mp2%free_hfx_buffer)) free_hfx_buffer = .FALSE.
      END IF

      IF (.NOT. do_kpoints_cubic_RPA) THEN
      IF (virial%pv_numer) THEN
         ! in the case of numerical stress we don't have to free the HFX integrals
         free_hfx_buffer = .FALSE.
         mp2_env%ri_mp2%free_hfx_buffer = .FALSE.
      END IF
      END IF

      ! calculate the matrix sigma_x - vxc for G0W0
      t3 = 0
      IF (mp2_env%ri_rpa%do_ri_g0w0) THEN
         CALL compute_vec_Sigma_x_minus_vxc_gw(qs_env, mp2_env, mos_mp2, E_ex_from_GW, t3, unit_nr)
      END IF
      IF (free_hfx_buffer) THEN
         CALL timeset(routineN//"_free_hfx", handle2)
         CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)
         n_threads = 1
!$       n_threads = omp_get_max_threads()

         DO irep = 1, n_rep_hf
            DO i_thread = 0, n_threads - 1
               actual_x_data => qs_env%x_data(irep, i_thread + 1)

               do_dynamic_load_balancing = .TRUE.
               IF (n_threads == 1 .OR. actual_x_data%memory_parameter%do_disk_storage) do_dynamic_load_balancing = .FALSE.

               IF (do_dynamic_load_balancing) THEN
                  my_bin_size = SIZE(actual_x_data%distribution_energy)
               ELSE
                  my_bin_size = 1
               END IF

               IF (.NOT. actual_x_data%memory_parameter%do_all_on_the_fly) THEN
                  CALL dealloc_containers(actual_x_data%store_ints, actual_x_data%memory_parameter%actual_memory_usage)
                  !  CALL alloc_containers(actual_x_data, my_bin_size, hfx_do_eval_energy)
                  !
                  !  DO bin=1, my_bin_size
                  !    maxval_container => actual_x_data%maxval_container(bin)
                  !    integral_containers => actual_x_data%integral_containers(:,bin)
                  !    CALL hfx_init_container(maxval_container, actual_x_data%memory_parameter%actual_memory_usage, .FALSE.)
                  !    DO i=1,64
                  !      CALL hfx_init_container(integral_containers(i), &
                  !               actual_x_data%memory_parameter%actual_memory_usage, .FALSE.)
                  !    END DO
                  !  END DO
               END IF
            END DO
         END DO
         CALL timestop(handle2)
      END IF

      Emp2 = 0.D+00
      Emp2_Cou = 0.D+00
      Emp2_ex = 0.D+00
      calc_ex = .TRUE.

      t1 = m_walltime()
      SELECT CASE (mp2_env%method)
      CASE (mp2_method_direct)
         IF (unit_nr > 0) WRITE (unit_nr, *)

         IF (nspins == 2) THEN
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Unrestricted Canonical Direct Methods:'
            ! for now, require the mos to be always present

            ! get the alpha coeff and eigenvalues
            CALL get_mo_set(mo_set=mos_mp2(1)%mo_set, &
                            nelectron=nelec_alpha, &
                            eigenvalues=mo_eigenvalues, &
                            mo_coeff=mo_coeff)
            ALLOCATE (C_alpha(dimen, dimen))

            ALLOCATE (Auto_alpha(dimen))

            CALL cp_fm_get_submatrix(mo_coeff, C_alpha, 1, 1, dimen, dimen, .FALSE.)
            Auto_alpha(:) = mo_eigenvalues(:)

            ! get the beta coeff and eigenvalues
            CALL get_mo_set(mo_set=mos_mp2(2)%mo_set, &
                            nelectron=nelec_beta, &
                            eigenvalues=mo_eigenvalues, &
                            mo_coeff=mo_coeff)
            ALLOCATE (C_beta(dimen, dimen))

            ALLOCATE (Auto_beta(dimen))

            CALL cp_fm_get_submatrix(mo_coeff, C_beta, 1, 1, dimen, dimen, .FALSE.)
            Auto_beta(:) = mo_eigenvalues(:)

            ! calculate the alpha-alpha MP2
            Emp2_AA = 0.0_dp
            Emp2_AA_Cou = 0.0_dp
            Emp2_AA_ex = 0.0_dp
            CALL mp2_direct_energy(dimen, nelec_alpha, nelec_alpha, mp2_biel, &
                                   mp2_env, C_alpha, Auto_alpha, Emp2_AA, Emp2_AA_Cou, Emp2_AA_ex, &
                                   qs_env, para_env, unit_nr)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy Alpha-Alpha = ', Emp2_AA
            IF (unit_nr > 0) WRITE (unit_nr, *)

            Emp2_BB = 0.0_dp
            Emp2_BB_Cou = 0.0_dp
            Emp2_BB_ex = 0.0_dp
            CALL mp2_direct_energy(dimen, nelec_beta, nelec_beta, mp2_biel, mp2_env, &
                                   C_beta, Auto_beta, Emp2_BB, Emp2_BB_Cou, Emp2_BB_ex, &
                                   qs_env, para_env, unit_nr)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy Beta-Beta= ', Emp2_BB
            IF (unit_nr > 0) WRITE (unit_nr, *)

            Emp2_AB = 0.0_dp
            Emp2_AB_Cou = 0.0_dp
            Emp2_AB_ex = 0.0_dp
            CALL mp2_direct_energy(dimen, nelec_alpha, nelec_beta, mp2_biel, mp2_env, C_alpha, &
                                   Auto_alpha, Emp2_AB, Emp2_AB_Cou, Emp2_AB_ex, &
                                   qs_env, para_env, unit_nr, C_beta, Auto_beta)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy Alpha-Beta= ', Emp2_AB
            IF (unit_nr > 0) WRITE (unit_nr, *)

            Emp2 = Emp2_AA + Emp2_BB + Emp2_AB*2.0_dp !+Emp2_BA
            Emp2_Cou = Emp2_AA_Cou + Emp2_BB_Cou + Emp2_AB_Cou*2.0_dp !+Emp2_BA
            Emp2_ex = Emp2_AA_ex + Emp2_BB_ex + Emp2_AB_ex*2.0_dp !+Emp2_BA

            Emp2_S = Emp2_AB*2.0_dp
            Emp2_T = Emp2_AA + Emp2_BB

         ELSE

            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Canonical Direct Methods:'
            CALL get_mo_set(mo_set=mos_mp2(1)%mo_set, &
                            nelectron=nelectron, &
                            eigenvalues=mo_eigenvalues, &
                            mo_coeff=mo_coeff)
            ALLOCATE (C(dimen, dimen))

            ALLOCATE (Auto(dimen))

            CALL cp_fm_get_submatrix(mo_coeff, C, 1, 1, dimen, dimen, .FALSE.)
            Auto(:) = mo_eigenvalues(:)

            CALL mp2_direct_energy(dimen, nelectron/2, nelectron/2, mp2_biel, mp2_env, C, Auto, Emp2, Emp2_Cou, Emp2_ex, &
                                   qs_env, para_env, unit_nr)

         END IF

      CASE (mp2_ri_optimize_basis)
         ! optimize ri basis set or tests for RI-MP2 gradients
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T3,A)') 'Optimization of the auxiliary RI-MP2 basis'
            WRITE (unit_nr, *)
         END IF
         CALL get_mo_set(mo_set=mos_mp2(1)%mo_set, &
                         nelectron=nelectron, &
                         eigenvalues=mo_eigenvalues, &
                         mo_coeff=mo_coeff)
         ALLOCATE (C(dimen, dimen))

         ALLOCATE (Auto(dimen))

         CALL cp_fm_get_submatrix(mo_coeff, C, 1, 1, dimen, dimen, .FALSE.)
         Auto(:) = mo_eigenvalues(:)

         IF (nspins == 2) THEN
            ! get the beta coeff and eigenvalues open-shell case
            CALL get_mo_set(mo_set=mos_mp2(2)%mo_set, &
                            nelectron=nelec_beta, &
                            eigenvalues=mo_eigenvalues, &
                            mo_coeff=mo_coeff)
            ALLOCATE (C_beta(dimen, dimen))

            ALLOCATE (Auto_beta(dimen))

            CALL cp_fm_get_submatrix(mo_coeff, C_beta, 1, 1, dimen, dimen, .FALSE.)
            Auto_beta(:) = mo_eigenvalues(:)

            ! optimize basis
            CALL optimize_ri_basis_main(Emp2, Emp2_Cou, Emp2_ex, Emp2_S, Emp2_T, dimen, natom, nelectron, &
                                        mp2_biel, mp2_env, C, Auto, &
                                        kind_of, qs_env, para_env, unit_nr, &
                                        nelec_beta, C_beta, Auto_beta)

         ELSE
            ! optimize basis
            CALL optimize_ri_basis_main(Emp2, Emp2_Cou, Emp2_ex, Emp2_S, Emp2_T, dimen, natom, nelectron/2, &
                                        mp2_biel, mp2_env, C, Auto, &
                                        kind_of, qs_env, para_env, unit_nr)
         END IF

      CASE (mp2_method_gpw)
         ! check if calculate the exchange contribution
         IF (mp2_env%scale_T == 0.0_dp) calc_ex = .FALSE.

         ! go with mp2_gpw
         CALL mp2_gpw_main(qs_env, mp2_env, Emp2, Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T, &
                           mos_mp2, para_env, unit_nr, calc_forces, calc_ex)

      CASE (ri_mp2_method_gpw)
         ! check if calculate the exchange contribution
         IF (mp2_env%scale_T == 0.0_dp) calc_ex = .FALSE.

         ! go with mp2_gpw
         CALL mp2_gpw_main(qs_env, mp2_env, Emp2, Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T, &
                           mos_mp2, para_env, unit_nr, calc_forces, calc_ex, do_ri_mp2=.TRUE.)

      CASE (ri_rpa_method_gpw)
         ! perform RI-RPA energy calculation (since most part of the calculation
         ! is actually equal to the RI-MP2-GPW we decided to put RPA in the MP2
         ! section to avoid code replication)

         calc_ex = .FALSE.

         ! go with ri_rpa_gpw
         CALL mp2_gpw_main(qs_env, mp2_env, Emp2, Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T, &
                           mos_mp2, para_env, unit_nr, calc_forces, calc_ex, do_ri_rpa=.TRUE.)

         ! Scale energy contributions
         Emp2 = Emp2*mp2_env%ri_rpa%scale_rpa
         mp2_env%ri_rpa%ener_axk = mp2_env%ri_rpa%ener_axk*mp2_env%ri_rpa%scale_rpa

      CASE (ri_mp2_laplace)
         ! perform RI-SOS-Laplace-MP2 energy calculation, most part of the code in common
         ! with the RI-RPA part

         ! In SOS-MP2 only the coulomb-like contribution of the MP2 energy is computed
         calc_ex = .FALSE.

         ! go with sos_laplace_mp2_gpw
         CALL mp2_gpw_main(qs_env, mp2_env, Emp2, Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T, &
                           mos_mp2, para_env, unit_nr, calc_forces, calc_ex, do_ri_sos_laplace_mp2=.TRUE.)

      CASE DEFAULT
         CPABORT("")
      END SELECT
      t2 = m_walltime()
      IF (unit_nr > 0) WRITE (unit_nr, *)
      IF (mp2_env%method .NE. ri_rpa_method_gpw) THEN
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.6)') 'Total MP2 Time=', t2 - t1
         IF (mp2_env%method == ri_mp2_laplace) THEN
            Emp2_S = Emp2
            Emp2_T = 0.0_dp
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy SO component (singlet) = ', Emp2_S
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'Scaling factor SO                 = ', mp2_env%scale_S
         ELSE
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Coulomb Energy = ', Emp2_Cou/2.0_dp
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Exchange Energy = ', Emp2_ex
            IF (nspins == 1) THEN
               ! valid only in the closed shell case
               Emp2_S = Emp2_Cou/2.0_dp
               IF (calc_ex) THEN
                  Emp2_T = Emp2_ex + Emp2_Cou/2.0_dp
               ELSE
                  ! unknown if Emp2_ex is not computed
                  Emp2_T = 0.0_dp
               ENDIF
            END IF
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy SO component (singlet) = ', Emp2_S
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy SS component (triplet) = ', Emp2_T
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'Scaling factor SO                 = ', mp2_env%scale_S
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'Scaling factor SS                 = ', mp2_env%scale_T
         END IF
         Emp2_S = Emp2_S*mp2_env%scale_S
         Emp2_T = Emp2_T*mp2_env%scale_T
         Emp2 = Emp2_S + Emp2_T
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'Second order perturbation energy  =   ', Emp2
      ELSE
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.6)') 'Total RI-RPA Time=', t2 - t1
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'RI-RPA energy  =   ', Emp2
         IF (mp2_env%ri_rpa%do_ri_axk) THEN
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'RI-RPA-AXK energy=', mp2_env%ri_rpa%ener_axk
         ENDIF
         IF (mp2_env%ri_rpa%do_rse) THEN
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'Diagonal singles correction (dRSE) = ', &
               mp2_env%ri_rpa%rse_corr_diag
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'Full singles correction (RSE) =', &
               mp2_env%ri_rpa%rse_corr
         ENDIF
      END IF
      IF (unit_nr > 0) WRITE (unit_nr, *)

      ! we have it !!!!
      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         Emp2 = Emp2 + mp2_env%ri_rpa%ener_axk
      ENDIF
      IF (mp2_env%ri_rpa%do_rse) THEN
         Emp2 = Emp2 + mp2_env%ri_rpa%rse_corr
      ENDIF
      energy%mp2 = Emp2
      energy%total = energy%total + Emp2

      DO ispin = 1, nspins
         CALL deallocate_mo_set(mo_set=mos_mp2(ispin)%mo_set)
      END DO
      DEALLOCATE (mos_mp2)

      ! if necessary reallocate hfx buffer
      IF (free_hfx_buffer .AND. (.NOT. calc_forces)) THEN
         CALL timeset(routineN//"_alloc_hfx", handle2)
         DO irep = 1, n_rep_hf
            DO i_thread = 0, n_threads - 1
               actual_x_data => qs_env%x_data(irep, i_thread + 1)

               do_dynamic_load_balancing = .TRUE.
               IF (n_threads == 1 .OR. actual_x_data%memory_parameter%do_disk_storage) do_dynamic_load_balancing = .FALSE.

               IF (do_dynamic_load_balancing) THEN
                  my_bin_size = SIZE(actual_x_data%distribution_energy)
               ELSE
                  my_bin_size = 1
               END IF

               IF (.NOT. actual_x_data%memory_parameter%do_all_on_the_fly) THEN
                  ! CALL dealloc_containers(actual_x_data, hfx_do_eval_energy)
                  CALL alloc_containers(actual_x_data%store_ints, my_bin_size)

                  DO bin = 1, my_bin_size
                     maxval_container => actual_x_data%store_ints%maxval_container(bin)
                     integral_containers => actual_x_data%store_ints%integral_containers(:, bin)
                     CALL hfx_init_container(maxval_container, actual_x_data%memory_parameter%actual_memory_usage, .FALSE.)
                     DO i = 1, 64
                        CALL hfx_init_container(integral_containers(i), actual_x_data%memory_parameter%actual_memory_usage, .FALSE.)
                     END DO
                  END DO
               END IF
            END DO
         END DO
         CALL timestop(handle2)
      END IF

      CALL hfx_release_basis_types(basis_parameter)

      ! if required calculate the EXX contribution from the DFT density
      IF (mp2_env%method == ri_rpa_method_gpw) THEN
         do_exx = .FALSE.
         hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")
         CALL section_vals_get(hfx_sections, explicit=do_exx)
         IF (do_exx) THEN
            do_gw = mp2_env%ri_rpa%do_ri_g0w0
            do_admm = mp2_env%ri_rpa%do_admm
            CALL calculate_exx(qs_env, unit_nr, do_gw, do_admm, E_ex_from_GW, t3)
         END IF
      END IF

      CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                        "DFT%XC%WF_CORRELATION%MP2%MP2_INFO")

      CALL timestop(handle)

   END SUBROUTINE mp2_main

! **************************************************************************************************
!> \brief ...
!> \param natom ...
!> \param max_nset ...
!> \param index_table ...
!> \param basis_parameter ...
!> \param kind_of ...
! **************************************************************************************************
   PURE SUBROUTINE build_index_table(natom, max_nset, index_table, basis_parameter, kind_of)
      INTEGER, INTENT(IN)                                :: natom, max_nset
      INTEGER, DIMENSION(natom, max_nset), INTENT(OUT)   :: index_table
      TYPE(hfx_basis_type), DIMENSION(:), POINTER        :: basis_parameter
      INTEGER, DIMENSION(natom), INTENT(IN)              :: kind_of

      INTEGER                                            :: counter, iatom, ikind, iset, nset

      index_table = -HUGE(0)
      counter = 0
      DO iatom = 1, natom
         ikind = kind_of(iatom)
         nset = basis_parameter(ikind)%nset
         DO iset = 1, nset
            index_table(iatom, iset) = counter + 1
            counter = counter + basis_parameter(ikind)%nsgf(iset)
         END DO
      END DO

   END SUBROUTINE build_index_table

! **************************************************************************************************
!> \brief ...
!> \param dimen ...
!> \param occ_i ...
!> \param occ_j ...
!> \param mp2_biel ...
!> \param mp2_env ...
!> \param C_i ...
!> \param Auto_i ...
!> \param Emp2 ...
!> \param Emp2_Cou ...
!> \param Emp2_ex ...
!> \param qs_env ...
!> \param para_env ...
!> \param unit_nr ...
!> \param C_j ...
!> \param Auto_j ...
! **************************************************************************************************
   SUBROUTINE mp2_direct_energy(dimen, occ_i, occ_j, mp2_biel, mp2_env, C_i, Auto_i, Emp2, Emp2_Cou, Emp2_ex, &
                                qs_env, para_env, &
                                unit_nr, C_j, Auto_j)
      INTEGER, INTENT(IN)                                :: dimen, occ_i, occ_j
      TYPE(mp2_biel_type), INTENT(IN)                    :: mp2_biel
      TYPE(mp2_type), POINTER                            :: mp2_env
      REAL(KIND=dp), DIMENSION(dimen, dimen), INTENT(IN) :: C_i
      REAL(KIND=dp), DIMENSION(dimen), INTENT(IN)        :: Auto_i
      REAL(KIND=dp), INTENT(OUT)                         :: Emp2, Emp2_Cou, Emp2_ex
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(dimen, dimen), &
         INTENT(IN), OPTIONAL                            :: C_j
      REAL(KIND=dp), DIMENSION(dimen), INTENT(IN), &
         OPTIONAL                                        :: Auto_j

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_direct_energy', &
         routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: zero = 0.D+00

      INTEGER :: batch_number, color_sub, comm_sub, counter, elements_ij_proc, group_counter, &
         handle, i, i_batch, i_batch_start, i_group_counter, j, j_batch_start, j_group_counter, &
         last_batch, max_batch_number, max_batch_size, max_set, minimum_memory_needed, &
         my_batch_size, my_I_batch_size, my_I_occupied_end, my_I_occupied_start, my_J_batch_size, &
         my_J_occupied_end, my_J_occupied_start, natom, Ni_occupied, Nj_occupied, number_groups, &
         number_i_subset, number_j_subset, one, sqrt_number_groups, total_I_size_batch_group, &
         total_J_size_batch_group, virt_i, virt_j
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_sizes, batch_sizes_tmp, &
                                                            vector_batch_I_size_group, &
                                                            vector_batch_J_size_group
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: ij_list_proc, ij_list_proc_temp, &
                                                            ij_matrix
      LOGICAL                                            :: alpha_beta_case
      TYPE(cp_para_env_type), POINTER                    :: para_env_sub

      CALL timeset(routineN, handle)

      alpha_beta_case = .FALSE.
      IF (PRESENT(C_j) .AND. PRESENT(Auto_j)) alpha_beta_case = .TRUE.

      IF (unit_nr > 0 .AND. mp2_env%potential_parameter%potential_type == do_potential_TShPSC) THEN
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T64,F12.6,A5)') 'Truncated MP2 method, Rt=', &
            mp2_env%potential_parameter%cutoff_radius, ' Bohr'
      END IF

      ! create the local para env
      ! each para_env_sub corresponds to a group that is going to compute
      ! all the integrals. To each group a batch I is assigned and the
      ! communication takes place only inside the group
      number_groups = para_env%num_pe/mp2_env%mp2_num_proc
      IF (number_groups*mp2_env%mp2_num_proc /= para_env%num_pe) THEN
         CPABORT(" The number of processors needs to be a multiple of the processors per group. ")
      END IF
      IF (number_groups > occ_i*occ_j) THEN
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Number of groups greater then the number of IJ pairs!'
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Consider using more processors per group for improved efficiency'
      END IF

      color_sub = para_env%mepos/mp2_env%mp2_num_proc
      CALL mp_comm_split_direct(para_env%group, comm_sub, color_sub)
      NULLIFY (para_env_sub)
      CALL cp_para_env_create(para_env_sub, comm_sub)

      ! calculate the maximal size of the batch, according to the maximum RS size
      max_set = SIZE(mp2_biel%index_table, 2)
      minimum_memory_needed = (8*(max_set**4))/1024**2
      IF (minimum_memory_needed > mp2_env%mp2_memory) THEN
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T66,I12,A3)') 'Memory required below the minimum, new memory:', &
            minimum_memory_needed, 'MiB'
         mp2_env%mp2_memory = minimum_memory_needed
      END IF

      ! Distribute the batches over the groups in
      ! a rectangular fashion, bigger size for J index
      ! the sizes of the I batches should be as small as possible
      sqrt_number_groups = INT(SQRT(REAL(number_groups, KIND=dp)))
      DO i = 1, number_groups
         IF (MOD(number_groups, i) == 0) THEN
            IF (sqrt_number_groups/i <= 1) THEN
               number_j_subset = i
               EXIT
            END IF
         END IF
      END DO
      number_i_subset = number_groups/number_j_subset

      IF (number_i_subset < number_j_subset) THEN
         number_i_subset = number_j_subset
         number_j_subset = number_groups/number_i_subset
      END IF

      ! Distribute the I index and the J index over groups
      total_I_size_batch_group = occ_i/number_i_subset
      IF (total_I_size_batch_group < 1) total_I_size_batch_group = 1
      ALLOCATE (vector_batch_I_size_group(0:number_i_subset - 1))

      vector_batch_I_size_group = 0
      DO i = 0, number_i_subset - 1
         vector_batch_I_size_group(i) = total_I_size_batch_group
      END DO
      IF (SUM(vector_batch_I_size_group) /= occ_i) THEN
         one = 1
         IF (SUM(vector_batch_I_size_group) > occ_i) one = -1
         i = -1
         DO
            i = i + 1
            vector_batch_I_size_group(i) = vector_batch_I_size_group(i) + one
            IF (SUM(vector_batch_I_size_group) == occ_i) EXIT
            IF (i == number_i_subset - 1) i = -1
         END DO
      END IF

      total_J_size_batch_group = occ_j/number_j_subset
      IF (total_J_size_batch_group < 1) total_J_size_batch_group = 1
      ALLOCATE (vector_batch_J_size_group(0:number_j_subset - 1))

      vector_batch_J_size_group = 0
      DO i = 0, number_J_subset - 1
         vector_batch_J_size_group(i) = total_J_size_batch_group
      END DO
      IF (SUM(vector_batch_J_size_group) /= occ_j) THEN
         one = 1
         IF (SUM(vector_batch_J_size_group) > occ_j) one = -1
         i = -1
         DO
            i = i + 1
            vector_batch_J_size_group(i) = vector_batch_J_size_group(i) + one
            IF (SUM(vector_batch_J_size_group) == occ_j) EXIT
            IF (i == number_J_subset - 1) i = -1
         END DO
      END IF

      ! now the starting and ending I and J occupied orbitals are assigned to each group
      group_counter = 0
      i_group_counter = 0
      my_I_occupied_start = 1
      DO i = 0, number_i_subset - 1
         my_J_occupied_start = 1
         j_group_counter = 0
         DO j = 0, number_j_subset - 1
            group_counter = group_counter + 1
            IF (color_sub == group_counter - 1) EXIT
            my_J_occupied_start = my_J_occupied_start + vector_batch_J_size_group(j)
            j_group_counter = j_group_counter + 1
         END DO
         IF (color_sub == group_counter - 1) EXIT
         my_I_occupied_start = my_I_occupied_start + vector_batch_I_size_group(i)
         i_group_counter = i_group_counter + 1
      END DO
      my_I_occupied_end = my_I_occupied_start + vector_batch_I_size_group(i_group_counter) - 1
      my_I_batch_size = vector_batch_I_size_group(i_group_counter)
      my_J_occupied_end = my_J_occupied_start + vector_batch_J_size_group(j_group_counter) - 1
      my_J_batch_size = vector_batch_J_size_group(j_group_counter)

      DEALLOCATE (vector_batch_I_size_group)
      DEALLOCATE (vector_batch_J_size_group)

      max_batch_size = MIN( &
                       MAX(1, &
                           INT(mp2_env%mp2_memory*INT(1024, KIND=int_8)**2/ &
                               (8*(2*dimen - occ_i)*INT(dimen, KIND=int_8)*my_J_batch_size/para_env_sub%num_pe))) &
                       , my_I_batch_size)
      IF (max_batch_size < 1) THEN
         max_batch_size = INT((8*(occ_i + 1)*INT(dimen, KIND=int_8)**2/para_env%num_pe)/1024**2)
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T72,I6,A3)') 'More memory required, at least:', max_batch_size, 'MiB'
         max_batch_size = 1
      END IF

      ! create the size of the batches inside the group
      my_batch_size = my_I_batch_size
      ALLOCATE (batch_sizes(my_batch_size))

      batch_sizes = -HUGE(0)
      batch_number = 0
      DO i = 1, my_batch_size
         IF (i*max_batch_size > my_batch_size) EXIT
         batch_number = batch_number + 1
         batch_sizes(i) = max_batch_size
      END DO
      last_batch = my_batch_size - max_batch_size*batch_number
      IF (last_batch > 0) THEN
         batch_number = batch_number + 1
         batch_sizes(batch_number) = last_batch
      END IF

      ALLOCATE (batch_sizes_tmp(batch_number))
      batch_sizes_tmp(1:batch_number) = batch_sizes(1:batch_number)
      DEALLOCATE (batch_sizes)
      ALLOCATE (batch_sizes(batch_number))
      batch_sizes(:) = batch_sizes_tmp
      DEALLOCATE (batch_sizes_tmp)

      max_batch_size = MAXVAL(batch_sizes)
      CALL mp_max(max_batch_size, para_env%group)
      max_batch_number = batch_number
      CALL mp_max(max_batch_number, para_env%group)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A,T76,I5)') 'Maximum used batch size: ', max_batch_size
         WRITE (unit_nr, '(T3,A,T76,I5)') 'Number of integral recomputations: ', max_batch_number
         CALL m_flush(unit_nr)
      END IF

      ! Batches sizes exceed the occupied orbitals allocated for group
      CPASSERT(SUM(batch_sizes) <= my_batch_size)

      virt_i = dimen - occ_i
      virt_j = dimen - occ_j
      natom = SIZE(mp2_biel%index_table, 1)

      CALL mp_sync(para_env%group)
      Emp2 = zero
      Emp2_Cou = zero
      Emp2_ex = zero
      i_batch_start = my_I_occupied_start - 1
      j_batch_start = my_J_occupied_start - 1
      Nj_occupied = my_J_batch_size
      DO i_batch = 1, batch_number

         Ni_occupied = batch_sizes(i_batch)

         counter = -1
         ALLOCATE (ij_matrix(Ni_occupied, Nj_occupied))

         ij_matrix = 0
         DO i = 1, Ni_occupied
            DO j = 1, Nj_occupied
               counter = counter + 1
               IF (MOD(counter, para_env_sub%num_pe) == para_env_sub%mepos) THEN
                  ij_matrix(i, j) = ij_matrix(i, j) + 1
               END IF
            END DO
         END DO

         ALLOCATE (ij_list_proc_temp(Ni_occupied*occ_j, 2))

         elements_ij_proc = 0
         DO i = 1, Ni_occupied
            DO j = 1, Nj_occupied
               IF (ij_matrix(i, j) == 0) CYCLE
               elements_ij_proc = elements_ij_proc + 1
               ij_list_proc_temp(elements_ij_proc, 1) = i
               ij_list_proc_temp(elements_ij_proc, 2) = j
            END DO
         END DO
         DEALLOCATE (ij_matrix)

         ALLOCATE (ij_list_proc(elements_ij_proc, 2))
         DO i = 1, elements_ij_proc
            ij_list_proc(i, 1) = ij_list_proc_temp(i, 1)
            ij_list_proc(i, 2) = ij_list_proc_temp(i, 2)
         END DO
         DEALLOCATE (ij_list_proc_temp)

         IF (.NOT. alpha_beta_case) THEN
            CALL mp2_canonical_direct_single_batch(Emp2, Emp2_Cou, Emp2_ex, mp2_env, qs_env, para_env_sub, &
                                                   mp2_biel, dimen, C_i, Auto_i, i_batch_start, Ni_occupied, occ_i, &
                                                   elements_ij_proc, ij_list_proc, Nj_occupied, j_batch_start)
         ELSE
            CALL mp2_canonical_direct_single_batch(Emp2, Emp2_Cou, Emp2_ex, mp2_env, qs_env, para_env_sub, &
                                                   mp2_biel, dimen, C_i, Auto_i, i_batch_start, Ni_occupied, occ_i, &
                                                   elements_ij_proc, ij_list_proc, Nj_occupied, j_batch_start, &
                                                   occ_j, C_j, Auto_j)
         END IF

         i_batch_start = i_batch_start + Ni_occupied

         DEALLOCATE (ij_list_proc)

      END DO

      CALL mp_sum(Emp2_Cou, para_env%group)
      CALL mp_sum(Emp2_Ex, para_env%group)
      CALL mp_sum(Emp2, para_env%group)

      CALL cp_para_env_release(para_env_sub)

      CALL timestop(handle)

   END SUBROUTINE mp2_direct_energy

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param unit_nr ...
!> \param do_gw ...
!> \param do_admm ...
!> \param E_ex_from_GW ...
!> \param t3 ...
! **************************************************************************************************
   SUBROUTINE calculate_exx(qs_env, unit_nr, do_gw, do_admm, E_ex_from_GW, t3)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: do_gw, do_admm
      REAL(KIND=dp), INTENT(IN)                          :: E_ex_from_GW, t3

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_exx', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, irep, n_rep_hf, ns
      REAL(KIND=dp)                                      :: eh1, ehfx, t1, t2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux_fit, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_2d, rho_ao_2d
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      NULLIFY (hfx_sections, input, para_env, matrix_ks, matrix_ks_aux_fit, rho, rho_ao)

      ! if ADMM for exchange energy, get rho_aux_fit and matrix_ks_aux_fit
      IF (do_admm) THEN
         CALL get_qs_env(qs_env=qs_env, &
                         input=input, &
                         para_env=para_env, &
                         energy=energy, &
                         rho_aux_fit=rho, &
                         matrix_ks=matrix_ks, &
                         matrix_ks_aux_fit=matrix_ks_aux_fit)
      ELSE
         CALL get_qs_env(qs_env=qs_env, &
                         input=input, &
                         para_env=para_env, &
                         energy=energy, &
                         rho=rho, &
                         matrix_ks=matrix_ks, &
                         matrix_ks_aux_fit=matrix_ks_aux_fit)
      END IF

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")

      CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)

      ! put matrix_ks to zero
      DO i = 1, SIZE(matrix_ks)
         CALL dbcsr_set(matrix_ks(i)%matrix, 0.0_dp)
         IF (do_admm) THEN
            CALL dbcsr_set(matrix_ks_aux_fit(i)%matrix, 0.0_dp)
         END IF
      END DO

      ! Remove the Exchange-correlation energy contributions from the total energy
      energy%total = energy%total - (energy%exc + energy%exc1 + energy%ex + &
                                     energy%exc_aux_fit + energy%exc1_aux_fit)

      energy%exc = 0.0_dp
      energy%exc1 = 0.0_dp
      energy%exc_aux_fit = 0.0_dp
      energy%exc1_aux_fit = 0.0_dp
      energy%ex = 0.0_dp

      ! take the exact exchange energy from GW or calculate it
      IF (do_gw) THEN
         energy%total = energy%total + E_ex_from_GW
         energy%ex = E_ex_from_GW
      ELSE

         ehfx = 0.0_dp
         DO irep = 1, n_rep_hf
            ns = SIZE(rho_ao)
            rho_ao_2d(1:ns, 1:1) => rho_ao(1:ns)
            ns = SIZE(matrix_ks)
            IF (do_admm) THEN
               matrix_ks_2d(1:ns, 1:1) => matrix_ks_aux_fit(1:ns)
            ELSE
               matrix_ks_2d(1:ns, 1:1) => matrix_ks(1:ns)
            END IF
            CALL integrate_four_center(qs_env, qs_env%mp2_env%ri_rpa%x_data, matrix_ks_2d, eh1, &
                                       rho_ao_2d, hfx_sections, &
                                       para_env, .TRUE., irep, .TRUE., &
                                       ispin=1)
            ehfx = ehfx + eh1
         END DO

         ! include the EXX contribution to the total energy
         energy%ex = ehfx
         energy%total = energy%total + energy%ex
      END IF

      t2 = m_walltime()

      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.6)') 'Total EXX Time=', t2 - t1 + t3
      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'EXX energy  =   ', energy%ex

      ! reset to zero the Hartree-Fock energy
      energy%ex = 0.0_dp

      CALL timestop(handle)

   END SUBROUTINE calculate_exx

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mp2_env ...
!> \param mos_mp2 ...
!> \param energy_ex ...
!> \param t3 ...
!> \param unit_nr ...
!> \par History
!>      04.2015 created
!> \author Jan Wilhelm
! **************************************************************************************************
   SUBROUTINE compute_vec_Sigma_x_minus_vxc_gw(qs_env, mp2_env, mos_mp2, energy_ex, t3, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mp2_type), POINTER                            :: mp2_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos_mp2
      REAL(KIND=dp), INTENT(OUT)                         :: energy_ex, t3
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_vec_Sigma_x_minus_vxc_gw', &
         routineP = moduleN//':'//routineN

      CHARACTER(4)                                       :: occ_virt
      CHARACTER(LEN=40)                                  :: line
      INTEGER :: dimen, gw_corr_lev_occ, gw_corr_lev_virt, handle, homo, i_img, ikp, irep, ispin, &
         iunit, myfun, myfun_aux, myfun_prim, n_level_gw, n_level_gw_ref, n_rep_hf, nkp, nmo, ns, &
         nspins, print_exx
      LOGICAL                                            :: charge_constrain_tmp, do_admm_rpa, &
                                                            do_hfx, do_kpoints_cubic_RPA, &
                                                            do_ri_Sigma_x, really_read_line
      REAL(KIND=dp)                                      :: eh1, ehfx, exx_minus_vxc, hfx_fraction, &
                                                            t1, t2, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: vec_Sigma_x_minus_vxc_gw, &
                                                            vec_Sigma_x_minus_vxc_gw_im
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_ks_aux_fit, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER :: matrix_ks_2d, matrix_ks_kp_im, &
         matrix_ks_kp_re, matrix_ks_transl, matrix_sigma_x_minus_vxc, matrix_sigma_x_minus_vxc_im, &
         rho_ao_2d
      TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_tmp_2, mo_coeff_b
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: hfx_sections, input, xc_section, &
                                                            xc_section_admm_aux, &
                                                            xc_section_admm_prim

      NULLIFY (admm_env, matrix_ks, matrix_ks_aux_fit, rho_ao, matrix_sigma_x_minus_vxc, input, &
               xc_section, xc_section_admm_aux, xc_section_admm_prim, hfx_sections, rho, &
               dft_control, para_env, ks_env, mo_coeff, matrix_sigma_x_minus_vxc_im)

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      do_admm_rpa = mp2_env%ri_rpa%do_admm
      do_ri_Sigma_x = mp2_env%ri_g0w0%do_ri_Sigma_x
      do_kpoints_cubic_RPA = qs_env%mp2_env%ri_rpa_im_time%do_im_time_kpoints

      IF (do_kpoints_cubic_RPA) THEN
         CPASSERT(do_ri_Sigma_x)
      END IF

      ! In case we do admm for RPA/GW, get the aux. density matrix rho_aux_fit, otherwise
      ! just the normal density matrix rho
      IF (do_admm_rpa) THEN
         CALL get_qs_env(qs_env, &
                         admm_env=admm_env, &
                         matrix_ks=matrix_ks, &
                         matrix_ks_aux_fit=matrix_ks_aux_fit, &
                         rho_aux_fit=rho, &
                         input=input, &
                         dft_control=dft_control, &
                         para_env=para_env, &
                         ks_env=ks_env)

         ! RPA/GW with ADMM for EXX or the exchange self-energy only implemented for
         ! ADMM_PURIFICATION_METHOD  NONE
         ! METHOD                    BASIS_PROJECTION
         ! in the admm section
         CPASSERT(admm_env%purification_method == do_admm_purify_none)
         CPASSERT(dft_control%admm_control%method == do_admm_basis_projection)

      ELSE

         IF (do_kpoints_cubic_RPA) THEN

            CALL get_qs_env(qs_env, &
                            matrix_ks_kp=matrix_ks_transl, &
                            rho=rho, &
                            input=input, &
                            dft_control=dft_control, &
                            para_env=para_env, &
                            kpoints=kpoints, &
                            ks_env=ks_env)

            nkp = kpoints%nkp

         ELSE

            CALL get_qs_env(qs_env, &
                            matrix_ks=matrix_ks, &
                            rho=rho, &
                            input=input, &
                            dft_control=dft_control, &
                            para_env=para_env, &
                            ks_env=ks_env)

            nkp = 1

         END IF

      END IF

      nspins = dft_control%nspins

      IF (.NOT. do_kpoints_cubic_RPA) THEN
         CALL qs_rho_get(rho, rho_ao=rho_ao)
      END IF

      ! If ADMM we should make the ks matrix up-to-date
      IF (dft_control%do_admm) THEN
      DO ispin = 1, SIZE(matrix_ks)
         CALL admm_correct_for_eigenvalues(ispin, admm_env, matrix_ks(ispin)%matrix)
      END DO
      END IF

      IF (do_kpoints_cubic_RPA) THEN

         CALL allocate_matrix_ks_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)
         CALL transform_matrix_ks_to_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)

         DO ispin = 1, nspins
         DO i_img = 1, SIZE(matrix_ks_transl, 2)
            CALL dbcsr_set(matrix_ks_transl(ispin, i_img)%matrix, 0.0_dp)
         END DO
         END DO

      END IF

      ! initialize matrix_sigma_x_minus_vxc
      NULLIFY (matrix_sigma_x_minus_vxc)
      CALL dbcsr_allocate_matrix_set(matrix_sigma_x_minus_vxc, nspins, nkp)
      IF (do_kpoints_cubic_RPA) THEN
         NULLIFY (matrix_sigma_x_minus_vxc_im)
         CALL dbcsr_allocate_matrix_set(matrix_sigma_x_minus_vxc_im, nspins, nkp)
      END IF

      DO ispin = 1, nspins
         DO ikp = 1, nkp

            IF (do_kpoints_cubic_RPA) THEN

               ALLOCATE (matrix_sigma_x_minus_vxc(ispin, ikp)%matrix)
               CALL dbcsr_create(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, &
                                 template=matrix_ks_kp_re(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_symmetric)

               CALL dbcsr_copy(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, matrix_ks_kp_re(ispin, ikp)%matrix)
               CALL dbcsr_set(matrix_ks_kp_re(ispin, ikp)%matrix, 0.0_dp)

               ALLOCATE (matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix)
               CALL dbcsr_create(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix, &
                                 template=matrix_ks_kp_im(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)

               CALL dbcsr_copy(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix, matrix_ks_kp_im(ispin, ikp)%matrix)
               CALL dbcsr_set(matrix_ks_kp_im(ispin, ikp)%matrix, 0.0_dp)

            ELSE

               ALLOCATE (matrix_sigma_x_minus_vxc(ispin, ikp)%matrix)
               CALL dbcsr_create(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, &
                                 template=matrix_ks(1)%matrix)

               CALL dbcsr_copy(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, matrix_ks(ispin)%matrix)
               CALL dbcsr_set(matrix_ks(ispin)%matrix, 0.0_dp)

            END IF

         ENDDO
      ENDDO

      ! set DFT functional to none and hfx_fraction to zero
      hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%HF")
      CALL section_vals_get(hfx_sections, explicit=do_hfx)

      IF (do_hfx) THEN
         hfx_fraction = qs_env%x_data(1, 1)%general_parameter%fraction
         qs_env%x_data(:, :)%general_parameter%fraction = 0.0_dp
      END IF
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      CALL section_vals_val_get(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=xc_none)

      ! in ADMM, also set the XC functional for ADMM correction to none
      ! do not do this if we do ADMM for Sigma_x
      IF (dft_control%do_admm) THEN
         xc_section_admm_aux => section_vals_get_subs_vals(admm_env%xc_section_aux, &
                                                           "XC_FUNCTIONAL")
         CALL section_vals_val_get(xc_section_admm_aux, "_SECTION_PARAMETERS_", &
                                   i_val=myfun_aux)
         CALL section_vals_val_set(xc_section_admm_aux, "_SECTION_PARAMETERS_", &
                                   i_val=xc_none)

         ! the same for the primary basis
         xc_section_admm_prim => section_vals_get_subs_vals(admm_env%xc_section_primary, &
                                                            "XC_FUNCTIONAL")
         CALL section_vals_val_get(xc_section_admm_aux, "_SECTION_PARAMETERS_", &
                                   i_val=myfun_prim)
         CALL section_vals_val_set(xc_section_admm_prim, "_SECTION_PARAMETERS_", &
                                   i_val=xc_none)

         ! for ADMMQ/S, set the charge_constrain to false (otherwise wrong results)
         charge_constrain_tmp = .FALSE.
         IF (admm_env%charge_constrain) THEN
            admm_env%charge_constrain = .FALSE.
            charge_constrain_tmp = .TRUE.
         END IF

      END IF

      ! if we do ADMM for Sigma_x, we write the ADMM correction into matrix_ks_aux_fit
      ! and therefore we should set it to zero
      IF (do_admm_rpa) THEN
         DO ispin = 1, nspins
            CALL dbcsr_set(matrix_ks_aux_fit(ispin)%matrix, 0.0_dp)
         END DO
      END IF

      ! calculate KS-matrix without XC and without HF
      CALL qs_ks_build_kohn_sham_matrix(qs_env=qs_env, calculate_forces=.FALSE., &
                                        just_energy=.FALSE.)

      ! set the DFT functional and HF fraction back
      CALL section_vals_val_set(xc_section, "XC_FUNCTIONAL%_SECTION_PARAMETERS_", &
                                i_val=myfun)
      IF (do_hfx) THEN
         qs_env%x_data(:, :)%general_parameter%fraction = hfx_fraction
      END IF

      IF (dft_control%do_admm) THEN
         xc_section_admm_aux => section_vals_get_subs_vals(admm_env%xc_section_aux, &
                                                           "XC_FUNCTIONAL")
         xc_section_admm_prim => section_vals_get_subs_vals(admm_env%xc_section_primary, &
                                                            "XC_FUNCTIONAL")

         CALL section_vals_val_set(xc_section_admm_aux, "_SECTION_PARAMETERS_", &
                                   i_val=myfun_aux)
         CALL section_vals_val_set(xc_section_admm_prim, "_SECTION_PARAMETERS_", &
                                   i_val=myfun_prim)
         IF (charge_constrain_tmp) THEN
            admm_env%charge_constrain = .TRUE.
         END IF
      END IF

      IF (do_kpoints_cubic_RPA) THEN
         CALL transform_matrix_ks_to_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)
      END IF

      ! remove the single-particle part (kin. En + Hartree pot) and change the sign
      DO ispin = 1, nspins
         IF (do_kpoints_cubic_RPA) THEN
            DO ikp = 1, nkp
               CALL dbcsr_add(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, matrix_ks_kp_re(ispin, ikp)%matrix, -1.0_dp, 1.0_dp)
               CALL dbcsr_add(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix, matrix_ks_kp_im(ispin, ikp)%matrix, -1.0_dp, 1.0_dp)
            END DO
         ELSE
            CALL dbcsr_add(matrix_sigma_x_minus_vxc(ispin, 1)%matrix, matrix_ks(ispin)%matrix, -1.0_dp, 1.0_dp)
         END IF
      END DO

      IF (do_kpoints_cubic_RPA) THEN

         CALL transform_sigma_x_minus_vxc_to_MO_basis(kpoints, matrix_sigma_x_minus_vxc, &
                                                      matrix_sigma_x_minus_vxc_im, &
                                                      vec_Sigma_x_minus_vxc_gw, &
                                                      vec_Sigma_x_minus_vxc_gw_im, &
                                                      para_env, nmo, mp2_env)

      ELSE

         DO ispin = 1, nspins
            CALL dbcsr_set(matrix_ks(ispin)%matrix, 0.0_dp)
            IF (do_admm_rpa) THEN
               CALL dbcsr_set(matrix_ks_aux_fit(ispin)%matrix, 0.0_dp)
            END IF
         END DO

         hfx_sections => section_vals_get_subs_vals(input, "DFT%XC%WF_CORRELATION%RI_RPA%HF")

         CALL section_vals_get(hfx_sections, n_repetition=n_rep_hf)

         ! in most cases, we calculate the exchange self-energy here. But if we do only RI for
         ! the exchange self-energy, we do not calculate exchange here
         ehfx = 0.0_dp
         IF (.NOT. do_ri_Sigma_x) THEN

            ! add here HFX (=Sigma_exchange) to matrix_sigma_x_minus_vxc
            DO irep = 1, n_rep_hf
               ns = SIZE(rho_ao)
               rho_ao_2d(1:ns, 1:1) => rho_ao(1:ns)
               ns = SIZE(matrix_ks)
               IF (do_admm_rpa) THEN
                  matrix_ks_2d(1:ns, 1:1) => matrix_ks_aux_fit(1:ns)
               ELSE
                  matrix_ks_2d(1:ns, 1:1) => matrix_ks(1:ns)
               END IF

               CALL integrate_four_center(qs_env, qs_env%mp2_env%ri_rpa%x_data, matrix_ks_2d, eh1, &
                                          rho_ao_2d, hfx_sections, &
                                          para_env, .TRUE., irep, .TRUE., &
                                          ispin=1)
               ehfx = ehfx + eh1
            END DO
         END IF
         energy_ex = ehfx

         ! transform Fock-Matrix (calculated in integrate_four_center, written in matrix_ks_aux_fit in case
         ! of ADMM) from ADMM basis to primary basis
         IF (do_admm_rpa) THEN
            CALL admm_mo_merge_ks_matrix(qs_env)
         END IF

         DO ispin = 1, nspins
            CALL dbcsr_add(matrix_sigma_x_minus_vxc(ispin, 1)%matrix, matrix_ks(ispin)%matrix, 1.0_dp, 1.0_dp)
         END DO

         CALL dbcsr_desymmetrize(matrix_ks(1)%matrix, mo_coeff_b)
         CALL dbcsr_set(mo_coeff_b, 0.0_dp)

         ! Transform matrix_sigma_x_minus_vxc to MO basis
         DO ispin = 1, nspins

            CALL get_mo_set(mo_set=mos_mp2(ispin)%mo_set, &
                            mo_coeff=mo_coeff, &
                            nmo=nmo, &
                            homo=homo, &
                            nao=dimen)

            IF (ispin == 1) THEN

               ALLOCATE (vec_Sigma_x_minus_vxc_gw(nmo, nspins, nkp))
               vec_Sigma_x_minus_vxc_gw = 0.0_dp
            END IF

            CALL dbcsr_set(mo_coeff_b, 0.0_dp)
            CALL copy_fm_to_dbcsr(mo_coeff, mo_coeff_b, keep_sparsity=.FALSE.)

            ! initialize matrix_tmp and matrix_tmp2
            IF (ispin == 1) THEN
               CALL dbcsr_create(matrix_tmp, template=mo_coeff_b)
               CALL dbcsr_copy(matrix_tmp, mo_coeff_b)
               CALL dbcsr_set(matrix_tmp, 0.0_dp)

               CALL dbcsr_create(matrix_tmp_2, template=mo_coeff_b)
               CALL dbcsr_copy(matrix_tmp_2, mo_coeff_b)
               CALL dbcsr_set(matrix_tmp_2, 0.0_dp)
            END IF

            gw_corr_lev_occ = mp2_env%ri_g0w0%corr_mos_occ
            gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_virt
            ! if corrected occ/virt levels exceed the number of occ/virt levels,
            ! correct all occ/virt level energies
            IF (gw_corr_lev_occ > homo) gw_corr_lev_occ = homo
            IF (gw_corr_lev_virt > dimen - homo) gw_corr_lev_virt = dimen - homo
            IF (ispin == 1) THEN
               mp2_env%ri_g0w0%corr_mos_occ = gw_corr_lev_occ
               mp2_env%ri_g0w0%corr_mos_virt = gw_corr_lev_virt
            ELSE IF (ispin == 2) THEN
               ! ensure that the total number of corrected MOs is the same for alpha and beta, important
               ! for parallelization
               IF (mp2_env%ri_g0w0%corr_mos_occ + mp2_env%ri_g0w0%corr_mos_virt /= &
                   gw_corr_lev_occ + gw_corr_lev_virt) THEN
                  gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_occ + mp2_env%ri_g0w0%corr_mos_virt - gw_corr_lev_occ
               END IF
               mp2_env%ri_g0w0%corr_mos_occ_beta = gw_corr_lev_occ
               mp2_env%ri_g0w0%corr_mos_virt_beta = gw_corr_lev_virt

            END IF

            CALL dbcsr_multiply('N', 'N', 1.0_dp, matrix_sigma_x_minus_vxc(ispin, 1)%matrix, &
                                mo_coeff_b, 0.0_dp, matrix_tmp, first_column=homo + 1 - gw_corr_lev_occ, &
                                last_column=homo + gw_corr_lev_virt)

            CALL dbcsr_multiply('T', 'N', 1.0_dp, mo_coeff_b, &
                                matrix_tmp, 0.0_dp, matrix_tmp_2, first_row=homo + 1 - gw_corr_lev_occ, &
                                last_row=homo + gw_corr_lev_virt)

            CALL dbcsr_get_diag(matrix_tmp_2, vec_Sigma_x_minus_vxc_gw(:, ispin, 1))

            CALL dbcsr_set(matrix_tmp, 0.0_dp)
            CALL dbcsr_set(matrix_tmp_2, 0.0_dp)

         END DO

         CALL mp_sum(vec_Sigma_x_minus_vxc_gw, para_env%group)

      END IF

      CALL dbcsr_release(mo_coeff_b)
      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_tmp_2)
      IF (do_kpoints_cubic_RPA) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_ks_kp_re)
         CALL dbcsr_deallocate_matrix_set(matrix_ks_kp_im)
      END IF

      DO ispin = 1, nspins
         DO ikp = 1, nkp
            CALL dbcsr_release_p(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix)
            IF (do_kpoints_cubic_RPA) THEN
               CALL dbcsr_release_p(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix)
            END IF
         END DO
      END DO

      ALLOCATE (mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(nmo, nspins, nkp))

      print_exx = mp2_env%ri_g0w0%print_exx

      IF (print_exx == gw_print_exx) THEN

         IF (unit_nr > 0) THEN

            WRITE (unit_nr, '(T3,A)') 'Exchange energies'
            WRITE (unit_nr, '(T3,A)') '-----------------'
            WRITE (unit_nr, '(T3,A)') ''
            WRITE (unit_nr, '(T6,2A)') 'MO                    Sigma_x-vxc'
            DO n_level_gw = 1, gw_corr_lev_occ + gw_corr_lev_virt

               n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ
               IF (n_level_gw <= gw_corr_lev_occ) THEN
                  occ_virt = 'occ'
               ELSE
                  occ_virt = 'vir'
               END IF

               exx_minus_vxc = REAL(vec_Sigma_x_minus_vxc_gw(n_level_gw_ref, 1, 1)*evolt)

               WRITE (unit_nr, '(T4,I4,3A,3F21.3)') &
                  n_level_gw_ref, ' ( ', occ_virt, ')  ', exx_minus_vxc

            END DO
            WRITE (unit_nr, '(T3,A)') 'End of exchange energies'
            WRITE (unit_nr, '(T3,A)') '------------------------'
            WRITE (unit_nr, '(T3,A)') ''
         END IF

      END IF

      IF (print_exx == gw_read_exx) THEN

         CALL open_file(unit_number=iunit, file_name="exx.out")

         really_read_line = .FALSE.

         DO WHILE (.TRUE.)

            READ (iunit, '(A)') line

            IF (line == "  End of exchange energies              ") EXIT

            IF (really_read_line) THEN

               READ (line(1:7), *) n_level_gw_ref
               READ (line(17:40), *) tmp

               DO ikp = 1, SIZE(vec_Sigma_x_minus_vxc_gw, 3)
                  vec_Sigma_x_minus_vxc_gw(n_level_gw_ref, 1, ikp) = tmp/evolt
               END DO

            END IF

            IF (line == "     MO                    Sigma_x-vxc  ") really_read_line = .TRUE.

         END DO

         CALL close_file(iunit)

      END IF

      ! store vec_Sigma_x_minus_vxc_gw in the mp2_environment
      mp2_env%ri_g0w0%vec_Sigma_x_minus_vxc_gw(:, :, :) = vec_Sigma_x_minus_vxc_gw(:, :, :)

      ! clean up
      DEALLOCATE (matrix_sigma_x_minus_vxc, vec_Sigma_x_minus_vxc_gw)
      IF (do_kpoints_cubic_RPA) THEN
         DEALLOCATE (matrix_sigma_x_minus_vxc_im)
      END IF

      t2 = m_walltime()

      t3 = t2 - t1

      CALL timestop(handle)

   END SUBROUTINE compute_vec_Sigma_x_minus_vxc_gw

! **************************************************************************************************
!> \brief ...
!> \param kpoints ...
!> \param matrix_sigma_x_minus_vxc ...
!> \param matrix_sigma_x_minus_vxc_im ...
!> \param vec_Sigma_x_minus_vxc_gw ...
!> \param vec_Sigma_x_minus_vxc_gw_im ...
!> \param para_env ...
!> \param nmo ...
!> \param mp2_env ...
! **************************************************************************************************
   SUBROUTINE transform_sigma_x_minus_vxc_to_MO_basis(kpoints, matrix_sigma_x_minus_vxc, &
                                                      matrix_sigma_x_minus_vxc_im, vec_Sigma_x_minus_vxc_gw, &
                                                      vec_Sigma_x_minus_vxc_gw_im, para_env, nmo, mp2_env)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_sigma_x_minus_vxc, &
                                                            matrix_sigma_x_minus_vxc_im
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: vec_Sigma_x_minus_vxc_gw, &
                                                            vec_Sigma_x_minus_vxc_gw_im
      TYPE(cp_para_env_type), POINTER                    :: para_env
      INTEGER                                            :: nmo
      TYPE(mp2_type), POINTER                            :: mp2_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'transform_sigma_x_minus_vxc_to_MO_basis', &
         routineP = moduleN//':'//routineN

      INTEGER :: dimen, gw_corr_lev_occ, gw_corr_lev_virt, handle, homo, i_global, iiB, ikp, &
         ispin, j_global, jjB, ncol_local, nkp, nrow_local, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: imval, reval
      TYPE(cp_cfm_type), POINTER                         :: cfm_mos, cfm_sigma_x_minus_vxc, &
                                                            cfm_sigma_x_minus_vxc_mo_basis, cfm_tmp
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: fwork_im, fwork_re
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(mo_set_type), POINTER                         :: mo_set, mo_set_im, mo_set_re

      CALL timeset(routineN, handle)

      mo_set => kpoints%kp_env(1)%kpoint_env%mos(1, 1)%mo_set
      CALL get_mo_set(mo_set, nmo=nmo)

      nspins = SIZE(matrix_sigma_x_minus_vxc, 1)
      CALL get_kpoint_info(kpoints, nkp=nkp, kp_range=kp_range)

      ! if this CPASSERT is wrong, please make sure that the kpoint group size PARALLEL_GROUP_SIZE
      ! in the kpoint environment &DFT &KPOINTS is -1
      CPASSERT(kp_range(1) == 1 .AND. kp_range(2) == nkp)

      ALLOCATE (vec_Sigma_x_minus_vxc_gw(nmo, nspins, nkp))
      vec_Sigma_x_minus_vxc_gw = 0.0_dp

      ALLOCATE (vec_Sigma_x_minus_vxc_gw_im(nmo, nspins, nkp))
      vec_Sigma_x_minus_vxc_gw_im = 0.0_dp

      NULLIFY (fwork_re, fwork_im, cfm_mos, cfm_sigma_x_minus_vxc, &
               cfm_sigma_x_minus_vxc_mo_basis, cfm_tmp)
      CALL cp_fm_get_info(mo_set%mo_coeff, matrix_struct=matrix_struct)
      CALL cp_fm_create(fwork_re, matrix_struct)
      CALL cp_fm_create(fwork_im, matrix_struct)
      CALL cp_cfm_create(cfm_mos, matrix_struct)
      CALL cp_cfm_create(cfm_sigma_x_minus_vxc, matrix_struct)
      CALL cp_cfm_create(cfm_sigma_x_minus_vxc_mo_basis, matrix_struct)
      CALL cp_cfm_create(cfm_tmp, matrix_struct)

      CALL cp_cfm_get_info(matrix=cfm_sigma_x_minus_vxc_mo_basis, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      ! Transform matrix_sigma_x_minus_vxc to MO basis
      DO ikp = 1, nkp

         kp => kpoints%kp_env(ikp)%kpoint_env

         DO ispin = 1, nspins

            ! v_xc_n to fm matrix
            CALL copy_dbcsr_to_fm(matrix_sigma_x_minus_vxc(ispin, ikp)%matrix, fwork_re)
            CALL copy_dbcsr_to_fm(matrix_sigma_x_minus_vxc_im(ispin, ikp)%matrix, fwork_im)

            CALL cp_cfm_scale_and_add_fm(z_zero, cfm_sigma_x_minus_vxc, z_one, fwork_re)
            CALL cp_cfm_scale_and_add_fm(z_one, cfm_sigma_x_minus_vxc, gaussi, fwork_im)

            ! get real part (1) and imag. part (2) of the mo coeffs
            mo_set_re => kp%mos(1, ispin)%mo_set
            mo_set_im => kp%mos(2, ispin)%mo_set

            CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mos, z_one, mo_set_re%mo_coeff)
            CALL cp_cfm_scale_and_add_fm(z_one, cfm_mos, gaussi, mo_set_im%mo_coeff)

            ! tmp = V(k)*C(k)
            CALL cp_cfm_gemm('N', 'N', nmo, nmo, nmo, z_one, cfm_sigma_x_minus_vxc, &
                             cfm_mos, z_zero, cfm_tmp)

            ! V_n(k) = C^H(k)*tmp
            CALL cp_cfm_gemm('C', 'N', nmo, nmo, nmo, z_one, cfm_mos, cfm_tmp, &
                             z_zero, cfm_sigma_x_minus_vxc_mo_basis)

            DO jjB = 1, ncol_local

               j_global = col_indices(jjB)

               DO iiB = 1, nrow_local

                  i_global = row_indices(iiB)

                  IF (j_global == i_global .AND. i_global <= nmo) THEN

                     reval = REAL(cfm_sigma_x_minus_vxc_mo_basis%local_data(iiB, jjB))
                     imval = AIMAG(cfm_sigma_x_minus_vxc_mo_basis%local_data(iiB, jjB))

                     vec_Sigma_x_minus_vxc_gw(i_global, ispin, ikp) = reval
                     vec_Sigma_x_minus_vxc_gw_im(i_global, ispin, ikp) = imval

                  END IF

               END DO

            END DO

         END DO

      END DO

      CALL mp_sum(vec_Sigma_x_minus_vxc_gw, para_env%group)
      CALL mp_sum(vec_Sigma_x_minus_vxc_gw_im, para_env%group)

      ! also adjust in the case of kpoints too big gw_corr_lev_occ and gw_corr_lev_virt
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=kpoints%kp_env(1)%kpoint_env%mos(ispin, 1)%mo_set, &
                         homo=homo, nao=dimen)
         gw_corr_lev_occ = mp2_env%ri_g0w0%corr_mos_occ
         gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_virt
         ! if corrected occ/virt levels exceed the number of occ/virt levels,
         ! correct all occ/virt level energies
         IF (gw_corr_lev_occ > homo) gw_corr_lev_occ = homo
         IF (gw_corr_lev_virt > dimen - homo) gw_corr_lev_virt = dimen - homo
         IF (ispin == 1) THEN
            mp2_env%ri_g0w0%corr_mos_occ = gw_corr_lev_occ
            mp2_env%ri_g0w0%corr_mos_virt = gw_corr_lev_virt
         ELSE IF (ispin == 2) THEN
            ! ensure that the total number of corrected MOs is the same for alpha and beta, important
            ! for parallelization
            IF (mp2_env%ri_g0w0%corr_mos_occ + mp2_env%ri_g0w0%corr_mos_virt /= &
                gw_corr_lev_occ + gw_corr_lev_virt) THEN
               gw_corr_lev_virt = mp2_env%ri_g0w0%corr_mos_occ + mp2_env%ri_g0w0%corr_mos_virt - gw_corr_lev_occ
            END IF
            mp2_env%ri_g0w0%corr_mos_occ_beta = gw_corr_lev_occ
            mp2_env%ri_g0w0%corr_mos_virt_beta = gw_corr_lev_virt
         END IF
      END DO

      CALL cp_fm_release(fwork_re)
      CALL cp_fm_release(fwork_im)
      CALL cp_cfm_release(cfm_mos)
      CALL cp_cfm_release(cfm_sigma_x_minus_vxc)
      CALL cp_cfm_release(cfm_sigma_x_minus_vxc_mo_basis)
      CALL cp_cfm_release(cfm_tmp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_s ...
!> \param matrix_ks ...
!> \param mos ...
!> \param matrix_s_kp ...
!> \param matrix_ks_transl ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE get_gamma(matrix_s, matrix_ks, mos, matrix_s_kp, matrix_ks_transl, kpoints)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_ks
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, matrix_ks_transl
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: nspins

      nspins = SIZE(matrix_ks_transl, 1)

      matrix_ks(1:nspins) => matrix_ks_transl(1:nspins, 1)
      matrix_s(1:1) => matrix_s_kp(1:1, 1)
      mos(1:nspins) => kpoints%kp_env(1)%kpoint_env%mos(1:nspins, 1)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_ks_transl ...
!> \param matrix_ks_kp_re ...
!> \param matrix_ks_kp_im ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE transform_matrix_ks_to_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_transl, matrix_ks_kp_re, &
                                                            matrix_ks_kp_im
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'transform_matrix_ks_to_kp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikp, ispin, nkp, nspin
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      NULLIFY (sab_nl)
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, sab_nl=sab_nl, cell_to_index=cell_to_index)

      CPASSERT(ASSOCIATED(sab_nl))

      nspin = SIZE(matrix_ks_transl, 1)

      DO ikp = 1, nkp
         DO ispin = 1, nspin

            CALL dbcsr_set(matrix_ks_kp_re(ispin, ikp)%matrix, 0.0_dp)
            CALL dbcsr_set(matrix_ks_kp_im(ispin, ikp)%matrix, 0.0_dp)
            CALL rskp_transform(rmatrix=matrix_ks_kp_re(ispin, ikp)%matrix, &
                                cmatrix=matrix_ks_kp_im(ispin, ikp)%matrix, &
                                rsmat=matrix_ks_transl, ispin=ispin, &
                                xkp=xkp(1:3, ikp), cell_to_index=cell_to_index, sab_nl=sab_nl)

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_ks_transl ...
!> \param matrix_ks_kp_re ...
!> \param matrix_ks_kp_im ...
!> \param kpoints ...
! **************************************************************************************************
   SUBROUTINE allocate_matrix_ks_kp(matrix_ks_transl, matrix_ks_kp_re, matrix_ks_kp_im, kpoints)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_transl, matrix_ks_kp_re, &
                                                            matrix_ks_kp_im
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_matrix_ks_kp', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikp, ispin, nkp, nspin
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      NULLIFY (sab_nl)
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, sab_nl=sab_nl, cell_to_index=cell_to_index)

      CPASSERT(ASSOCIATED(sab_nl))

      nspin = SIZE(matrix_ks_transl, 1)

      NULLIFY (matrix_ks_kp_re, matrix_ks_kp_im)
      CALL dbcsr_allocate_matrix_set(matrix_ks_kp_re, nspin, nkp)
      CALL dbcsr_allocate_matrix_set(matrix_ks_kp_im, nspin, nkp)

      DO ikp = 1, nkp
      DO ispin = 1, nspin

         ALLOCATE (matrix_ks_kp_re(ispin, ikp)%matrix)
         ALLOCATE (matrix_ks_kp_im(ispin, ikp)%matrix)

         CALL dbcsr_create(matrix_ks_kp_re(ispin, ikp)%matrix, &
                           template=matrix_ks_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_create(matrix_ks_kp_im(ispin, ikp)%matrix, &
                           template=matrix_ks_transl(1, 1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)

         CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks_kp_re(ispin, ikp)%matrix, sab_nl)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks_kp_im(ispin, ikp)%matrix, sab_nl)

         CALL dbcsr_set(matrix_ks_kp_re(ispin, ikp)%matrix, 0.0_dp)
         CALL dbcsr_set(matrix_ks_kp_im(ispin, ikp)%matrix, 0.0_dp)

      END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE

END MODULE mp2

