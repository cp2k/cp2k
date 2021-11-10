!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate MP2 energy
!> \par History
!>      05.2011 created [Mauro Del Ben]
!> \author Mauro Del Ben
! **************************************************************************************************
MODULE mp2
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
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
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
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_get_info,&
                                              dbcsr_p_type
   USE hfx_types,                       ONLY: &
        alloc_containers, dealloc_containers, hfx_basis_info_type, hfx_basis_type, &
        hfx_container_type, hfx_create_basis_types, hfx_init_container, hfx_release_basis_types, &
        hfx_type
   USE input_constants,                 ONLY: cholesky_inverse,&
                                              cholesky_off,&
                                              do_eri_gpw,&
                                              do_eri_mme
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_flush,&
                                              m_memory,&
                                              m_walltime
   USE message_passing,                 ONLY: mp_max
   USE mp2_direct_method,               ONLY: mp2_direct_energy
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
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_matrix_pools,                 ONLY: mpools_create,&
                                              mpools_rebuild_fm_pools,&
                                              mpools_release,&
                                              qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type
   USE qs_scf_methods,                  ONLY: eigensolver,&
                                              eigensolver_symm
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE rpa_gw_sigma,                    ONLY: compute_vec_Sigma_x_minus_vxc_gw
   USE rpa_hfx,                         ONLY: calculate_exx
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

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mp2_main'

      INTEGER :: bin, cholesky_method, dimen, handle, handle2, i, i_thread, iatom, ikind, irep, &
         ispin, max_nset, my_bin_size, n_rep_hf, n_threads, nao, natom, ncol_block, ndep, &
         nfullcols_total, nfullrows_total, nkind, nrow_block, nspins, unit_nr
      INTEGER(KIND=int_8)                                :: mem
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, nelec
      LOGICAL :: calc_ex, do_admm, do_admm_rpa_exx, do_dynamic_load_balancing, do_exx, do_gw, &
         do_kpoints_cubic_RPA, free_hfx_buffer, update_xc_energy
      REAL(KIND=dp) :: E_ex_from_GW, Emp2, Emp2_AA, Emp2_AA_Cou, Emp2_AA_ex, Emp2_AB, Emp2_AB_Cou, &
         Emp2_AB_ex, Emp2_BB, Emp2_BB_Cou, Emp2_BB_ex, Emp2_Cou, Emp2_ex, Emp2_S, Emp2_T, maxocc, &
         mem_real, t1, t2, t3
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Auto
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: C
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
      END IF

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
      END IF

      IF (calc_forces) THEN
         CALL cite_reference(DelBen2015b)
         CALL cite_reference(Rybkin2016)
         IF (mp2_env%method /= ri_mp2_method_gpw) THEN
            CPABORT("No forces/gradients for the selected method.")
         END IF
      END IF

      IF (.NOT. do_kpoints_cubic_RPA) THEN
         IF (virial%pv_availability .AND. (.NOT. virial%pv_numer) .AND. mp2_env%eri_method == do_eri_mme) THEN
            CPABORT("analytical stress not implemented with ERI_METHOD = MME")
         END IF
      END IF

      IF (mp2_env%do_im_time .AND. mp2_env%eri_method .NE. do_eri_gpw) THEN
         mp2_env%mp2_num_proc = 1
      END IF

      IF (mp2_env%mp2_num_proc <= 0 .OR. mp2_env%mp2_num_proc > para_env%num_pe .OR. &
          MOD(para_env%num_pe, mp2_env%mp2_num_proc) .NE. 0) THEN
         IF (unit_nr > 0 .AND. mp2_env%mp2_num_proc .NE. -1) &
            WRITE (unit_nr, '(T3,A,T76,I5)') 'Requested number of processes per group:', mp2_env%mp2_num_proc
         mp2_env%mp2_num_proc = para_env%num_pe
      END IF
      IF (.NOT. mp2_env%do_im_time) THEN
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T76,I5)') 'Used number of processes per group:', mp2_env%mp2_num_proc
         IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T68,F9.2,A4)') 'Maximum allowed memory usage per MPI process:', &
            mp2_env%mp2_memory, ' MiB'
      END IF

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
      ALLOCATE (nelec(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, maxocc=maxocc, nelectron=nelec(ispin))

         NULLIFY (mos_mp2(ispin)%mo_set)
         CALL allocate_mo_set(mo_set=mos_mp2(ispin)%mo_set, &
                              nao=nao, &
                              nmo=nao, &
                              nelectron=nelec(ispin), &
                              n_el_f=REAL(nelec(ispin), dp), &
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
         IF (qs_env%x_data(1, 1)%do_hfx_ri) free_hfx_buffer = .FALSE.
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

         ALLOCATE (Auto(dimen, nspins))
         ALLOCATE (C(dimen, dimen, nspins))

         DO ispin = 1, nspins
            ! get the alpha coeff and eigenvalues
            CALL get_mo_set(mo_set=mos_mp2(ispin)%mo_set, &
                            eigenvalues=mo_eigenvalues, &
                            mo_coeff=mo_coeff)

            CALL cp_fm_get_submatrix(mo_coeff, C(:, :, ispin), 1, 1, dimen, dimen, .FALSE.)
            Auto(:, ispin) = mo_eigenvalues(:)
         END DO

         IF (nspins == 2) THEN
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Unrestricted Canonical Direct Methods:'
            ! for now, require the mos to be always present

            ! calculate the alpha-alpha MP2
            Emp2_AA = 0.0_dp
            Emp2_AA_Cou = 0.0_dp
            Emp2_AA_ex = 0.0_dp
            CALL mp2_direct_energy(dimen, nelec(1), nelec(1), mp2_biel, &
                                   mp2_env, C(:, :, 1), Auto(:, 1), Emp2_AA, Emp2_AA_Cou, Emp2_AA_ex, &
                                   qs_env, para_env, unit_nr)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy Alpha-Alpha = ', Emp2_AA
            IF (unit_nr > 0) WRITE (unit_nr, *)

            Emp2_BB = 0.0_dp
            Emp2_BB_Cou = 0.0_dp
            Emp2_BB_ex = 0.0_dp
            CALL mp2_direct_energy(dimen, nelec(2), nelec(2), mp2_biel, mp2_env, &
                                   C(:, :, 2), Auto(:, 2), Emp2_BB, Emp2_BB_Cou, Emp2_BB_ex, &
                                   qs_env, para_env, unit_nr)
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy Beta-Beta= ', Emp2_BB
            IF (unit_nr > 0) WRITE (unit_nr, *)

            Emp2_AB = 0.0_dp
            Emp2_AB_Cou = 0.0_dp
            Emp2_AB_ex = 0.0_dp
            CALL mp2_direct_energy(dimen, nelec(1), nelec(2), mp2_biel, mp2_env, C(:, :, 1), &
                                   Auto(:, 1), Emp2_AB, Emp2_AB_Cou, Emp2_AB_ex, &
                                   qs_env, para_env, unit_nr, C(:, :, 2), Auto(:, 2))
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'MP2 Energy Alpha-Beta= ', Emp2_AB
            IF (unit_nr > 0) WRITE (unit_nr, *)

            Emp2 = Emp2_AA + Emp2_BB + Emp2_AB*2.0_dp !+Emp2_BA
            Emp2_Cou = Emp2_AA_Cou + Emp2_BB_Cou + Emp2_AB_Cou*2.0_dp !+Emp2_BA
            Emp2_ex = Emp2_AA_ex + Emp2_BB_ex + Emp2_AB_ex*2.0_dp !+Emp2_BA

            Emp2_S = Emp2_AB*2.0_dp
            Emp2_T = Emp2_AA + Emp2_BB

         ELSE

            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A)') 'Canonical Direct Methods:'

            CALL mp2_direct_energy(dimen, nelec(1)/2, nelec(1)/2, mp2_biel, mp2_env, &
                                   C(:, :, 1), Auto(:, 1), Emp2, Emp2_Cou, Emp2_ex, &
                                   qs_env, para_env, unit_nr)

         END IF

         DEALLOCATE (C, Auto)

      CASE (mp2_ri_optimize_basis)
         ! optimize ri basis set or tests for RI-MP2 gradients
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, '(T3,A)') 'Optimization of the auxiliary RI-MP2 basis'
            WRITE (unit_nr, *)
         END IF

         ALLOCATE (Auto(dimen, nspins))
         ALLOCATE (C(dimen, dimen, nspins))

         DO ispin = 1, nspins
            ! get the alpha coeff and eigenvalues
            CALL get_mo_set(mo_set=mos_mp2(ispin)%mo_set, &
                            eigenvalues=mo_eigenvalues, &
                            mo_coeff=mo_coeff)

            CALL cp_fm_get_submatrix(mo_coeff, C(:, :, ispin), 1, 1, dimen, dimen, .FALSE.)
            Auto(:, ispin) = mo_eigenvalues(:)
         END DO

         ! optimize basis
         IF (nspins == 2) THEN
            CALL optimize_ri_basis_main(Emp2, Emp2_Cou, Emp2_ex, Emp2_S, Emp2_T, dimen, natom, nelec(1), &
                                        mp2_biel, mp2_env, C(:, :, 1), Auto(:, 1), &
                                        kind_of, qs_env, para_env, unit_nr, &
                                        nelec(2), C(:, :, 2), Auto(:, 2))

         ELSE
            CALL optimize_ri_basis_main(Emp2, Emp2_Cou, Emp2_ex, Emp2_S, Emp2_T, dimen, natom, nelec(1)/2, &
                                        mp2_biel, mp2_env, C(:, :, 1), Auto(:, 1), &
                                        kind_of, qs_env, para_env, unit_nr)
         END IF

         DEALLOCATE (Auto, C)

      CASE (mp2_method_gpw)
         ! check if calculate the exchange contribution
         IF (mp2_env%scale_T == 0.0_dp .AND. (nspins == 2)) calc_ex = .FALSE.

         ! go with mp2_gpw
         CALL mp2_gpw_main(qs_env, mp2_env, Emp2, Emp2_Cou, Emp2_EX, Emp2_S, Emp2_T, &
                           mos_mp2, para_env, unit_nr, calc_forces, calc_ex)

      CASE (ri_mp2_method_gpw)
         ! check if calculate the exchange contribution
         IF (mp2_env%scale_T == 0.0_dp .AND. (nspins == 2)) calc_ex = .FALSE.

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
               END IF
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

         update_xc_energy = .TRUE.
         IF (mp2_env%ri_rpa%do_ri_g0w0 .AND. .NOT. mp2_env%ri_g0w0%update_xc_energy) update_xc_energy = .FALSE.
         IF (.NOT. update_xc_energy) Emp2 = 0.0_dp

         IF (unit_nr > 0 .AND. update_xc_energy) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'RI-RPA energy  =   ', Emp2
         IF (mp2_env%ri_rpa%do_ri_axk) THEN
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'RI-RPA-AXK energy=', mp2_env%ri_rpa%ener_axk
         END IF
         IF (mp2_env%ri_rpa%do_rse) THEN
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'Diagonal singles correction (dRSE) = ', &
               mp2_env%ri_rpa%rse_corr_diag
            IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') 'Full singles correction (RSE) =', &
               mp2_env%ri_rpa%rse_corr
         END IF
      END IF
      IF (unit_nr > 0) WRITE (unit_nr, *)

      ! we have it !!!!
      IF (mp2_env%ri_rpa%do_ri_axk) THEN
         Emp2 = Emp2 + mp2_env%ri_rpa%ener_axk
      END IF
      IF (mp2_env%ri_rpa%do_rse) THEN
         Emp2 = Emp2 + mp2_env%ri_rpa%rse_corr
      END IF
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
!> \param matrix_s ...
!> \param matrix_ks ...
!> \param mos ...
!> \param matrix_s_kp ...
!> \param matrix_ks_transl ...
!> \param kpoints ...
! **************************************************************************************************
   PURE SUBROUTINE get_gamma(matrix_s, matrix_ks, mos, matrix_s_kp, matrix_ks_transl, kpoints)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_ks
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, matrix_ks_transl
      TYPE(kpoint_type), POINTER                         :: kpoints

      INTEGER                                            :: nspins

      nspins = SIZE(matrix_ks_transl, 1)

      matrix_ks(1:nspins) => matrix_ks_transl(1:nspins, 1)
      matrix_s(1:1) => matrix_s_kp(1:1, 1)
      mos(1:nspins) => kpoints%kp_env(1)%kpoint_env%mos(1:nspins, 1)

   END SUBROUTINE get_gamma

END MODULE mp2

