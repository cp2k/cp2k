!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to somehow generate an initial guess
!> \par History
!>       2006.03 Moved here from qs_scf.F [Joost VandeVondele]
! **************************************************************************************************
MODULE qs_initial_guess
   USE atom_kind_orbitals,              ONLY: calculate_atomic_orbitals
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              cp_fm_to_dbcsr_row_template
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_init_random, cp_fm_p_type, &
        cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_checksum, dbcsr_copy, dbcsr_dot, dbcsr_filter, dbcsr_get_diag, &
        dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_get_occupation, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_multiply, dbcsr_nfullrows_total, dbcsr_p_type, dbcsr_release, dbcsr_scale, &
        dbcsr_set, dbcsr_set_diag, dbcsr_type, dbcsr_verify_matrix
   USE external_potential_types,        ONLY: all_potential_type,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: atomic_guess,&
                                              core_guess,&
                                              history_guess,&
                                              mopac_guess,&
                                              no_guess,&
                                              random_guess,&
                                              restart_guess,&
                                              sparse_guess
   USE input_cp2k_hfx,                  ONLY: ri_mo
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE kpoint_io,                       ONLY: read_kpoints_restart
   USE kpoint_types,                    ONLY: kpoint_type
   USE message_passing,                 ONLY: mp_bcast,&
                                              mp_sum
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_io,                        ONLY: read_mo_set_from_restart,&
                                              wfn_restart_file_name
   USE qs_mo_methods,                   ONLY: make_basis_lowdin,&
                                              make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_restrict,&
                                              reassign_allocated_mos
   USE qs_mom_methods,                  ONLY: do_mom_guess
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: eigensolver,&
                                              eigensolver_simple
   USE qs_scf_types,                    ONLY: block_davidson_diag_method_nr,&
                                              block_krylov_diag_method_nr,&
                                              general_diag_method_nr,&
                                              ot_diag_method_nr,&
                                              qs_scf_env_type
   USE qs_wf_history_methods,           ONLY: wfi_update
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: sort
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_initial_guess'

   PUBLIC ::  calculate_first_density_matrix, calculate_atomic_block_dm, calculate_mopac_dm
   PUBLIC ::  calculate_atomic_fock_matrix

   TYPE atom_matrix_type
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER   :: mat
   END TYPE atom_matrix_type

CONTAINS

! **************************************************************************************************
!> \brief can use a variety of methods to come up with an initial
!>      density matrix and optionally an initial wavefunction
!> \param scf_env  SCF environment information
!> \param qs_env   QS environment
!> \par History
!>      03.2006 moved here from qs_scf [Joost VandeVondele]
!>      06.2007 allow to skip the initial guess [jgh]
!>      08.2014 kpoints [JGH]
!>      10.2019 tot_corr_zeff, switch_surf_dip [SGh]
!> \note
!>      badly needs to be split in subroutines each doing one of the possible
!>      schemes
! **************************************************************************************************
   SUBROUTINE calculate_first_density_matrix(scf_env, qs_env)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_first_density_matrix'

      CHARACTER(LEN=default_path_length)                 :: file_name, filename
      INTEGER :: atom_a, blk, density_guess, group, handle, homo, i, iatom, ic, icol, id_nr, &
         ikind, irow, iseed(4), ispin, istart_col, istart_row, j, last_read, n, n_cols, n_rows, &
         nao, natom, natoms, natoms_tmp, nblocks, nelectron, nmo, nmo_tmp, not_read, nsgf, nspin, &
         nvec, ounit, qs_env_id, safe_density_guess, size_atomic_kind_set, z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, kind_of, last_sgf
      INTEGER, DIMENSION(2)                              :: nelectron_spin
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, elec_conf, nelec_kind, &
                                                            sort_kind
      LOGICAL                                            :: did_guess, do_hfx_ri_mo, do_kpoints, &
                                                            do_std_diag, exist, has_unit_metric, &
                                                            natom_mismatch, need_mos, need_wm, &
                                                            ofgpw
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: buff, buff2
      REAL(dp), DIMENSION(:, :), POINTER                 :: pdata
      REAL(KIND=dp)                                      :: checksum, eps, length, maxocc, occ, &
                                                            rscale, tot_corr_zeff, trps1, zeff
      REAL(KIND=dp), DIMENSION(0:3)                      :: edftb
      TYPE(atom_matrix_type), DIMENSION(:), POINTER      :: pmat
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: work1
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_struct, ao_mo_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, moa, mob, ortho, sv, work2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: h_core_sparse, matrix_ks, p_rmpv, &
                                                            s_sparse
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h_kp, matrix_ks_kp, matrix_s_kp, &
                                                            rho_ao_kp
      TYPE(dbcsr_type)                                   :: mo_dbcsr, mo_tmp_dbcsr
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array, mos_last_converged
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, subsys_section

      logger => cp_get_default_logger()
      NULLIFY (atomic_kind, qs_kind, mo_coeff, sv, orb_basis_set, atomic_kind_set, &
               qs_kind_set, particle_set, ortho, work2, work1, mo_array, s_sparse, &
               scf_control, dft_control, p_rmpv, ortho, work2, work1, para_env, &
               s_sparse, scf_control, dft_control, h_core_sparse, matrix_ks, rho, &
               mos_last_converged)
      NULLIFY (dft_section, input, subsys_section)
      NULLIFY (matrix_s_kp, matrix_h_kp, matrix_ks_kp, rho_ao_kp)
      NULLIFY (moa, mob)
      NULLIFY (atom_list, elec_conf, kpoints)
      edftb = 0.0_dp
      tot_corr_zeff = 0.0_dp

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      mos=mo_array, &
                      matrix_s_kp=matrix_s_kp, &
                      matrix_h_kp=matrix_h_kp, &
                      matrix_ks_kp=matrix_ks_kp, &
                      input=input, &
                      scf_control=scf_control, &
                      id_nr=qs_env_id, &
                      dft_control=dft_control, &
                      has_unit_metric=has_unit_metric, &
                      do_kpoints=do_kpoints, &
                      kpoints=kpoints, &
                      rho=rho, &
                      nelectron_spin=nelectron_spin, &
                      para_env=para_env, &
                      x_data=x_data)

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      IF (dft_control%switch_surf_dip) THEN
         CALL get_qs_env(qs_env, mos_last_converged=mos_last_converged)
      ENDIF

      ! just initialize the first image, the other density are set to zero
      DO ispin = 1, dft_control%nspins
         DO ic = 1, SIZE(rho_ao_kp, 2)
            CALL dbcsr_set(rho_ao_kp(ispin, ic)%matrix, 0.0_dp)
         END DO
      END DO
      s_sparse => matrix_s_kp(:, 1)
      h_core_sparse => matrix_h_kp(:, 1)
      matrix_ks => matrix_ks_kp(:, 1)
      p_rmpv => rho_ao_kp(:, 1)

      work1 => scf_env%scf_work1
      work2 => scf_env%scf_work2
      ortho => scf_env%ortho

      dft_section => section_vals_get_subs_vals(input, "DFT")

      nspin = dft_control%nspins
      ofgpw = dft_control%qs_control%ofgpw
      density_guess = scf_control%density_guess
      do_std_diag = .FALSE.

      do_hfx_ri_mo = .FALSE.
      IF (ASSOCIATED(x_data)) THEN
         IF (x_data(1, 1)%do_hfx_ri) THEN
            IF (x_data(1, 1)%ri_data%flavor == ri_mo) do_hfx_ri_mo = .TRUE.
         ENDIF
      ENDIF

      IF (ASSOCIATED(scf_env%krylov_space)) do_std_diag = (scf_env%krylov_space%eps_std_diag > 0.0_dp)

      need_mos = scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr .OR. &
                 (scf_env%method == block_krylov_diag_method_nr .AND. .NOT. do_std_diag) &
                 .OR. dft_control%do_admm .OR. scf_env%method == block_davidson_diag_method_nr &
                 .OR. do_hfx_ri_mo

      safe_density_guess = atomic_guess
      IF (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb .OR. &
          dft_control%qs_control%xtb) THEN
         IF (density_guess == atomic_guess) density_guess = mopac_guess
         safe_density_guess = mopac_guess
      END IF

      IF (scf_control%use_ot .AND. &
          (.NOT. ((density_guess == random_guess) .OR. &
                  (density_guess == atomic_guess) .OR. &
                  (density_guess == core_guess) .OR. &
                  (density_guess == mopac_guess) .OR. &
                  (density_guess == sparse_guess) .OR. &
                  (((density_guess == restart_guess) .OR. &
                    (density_guess == history_guess)) .AND. &
                   (scf_control%level_shift == 0.0_dp))))) THEN
         CALL cp_abort(__LOCATION__, &
                       "OT needs GUESS ATOMIC / CORE / RANDOM / SPARSE / RESTART / HISTORY RESTART: other options NYI")
      END IF

      ! if a restart was requested, check that the file exists,
      ! if not we fall back to an atomic guess. No kidding, the file name should remain
      ! in sync with read_mo_set_from_restart
      id_nr = 0
      IF (density_guess == restart_guess) THEN
         ! only check existence on I/O node, otherwise if file exists there but
         ! not on compute nodes, everything goes crazy even though only I/O
         ! node actually reads the file
         IF (do_kpoints) THEN
            IF (para_env%ionode) THEN
               CALL wfn_restart_file_name(file_name, exist, dft_section, logger, kp=.TRUE.)
            END IF
         ELSE
            IF (para_env%ionode) THEN
               CALL wfn_restart_file_name(file_name, exist, dft_section, logger)
            END IF
         ENDIF
         CALL mp_bcast(exist, para_env%source, para_env%group)
         CALL mp_bcast(file_name, para_env%source, para_env%group)
         IF (.NOT. exist) THEN
            CALL cp_warn(__LOCATION__, &
                         "User requested to restart the wavefunction from the file named: "// &
                         TRIM(file_name)//". This file does not exist. Please check the existence of"// &
                         " the file or change properly the value of the keyword WFN_RESTART_FILE_NAME."// &
                         " Calculation continues using ATOMIC GUESS. ")
            density_guess = safe_density_guess
         END IF
      ELSE IF (density_guess == history_guess) THEN
         IF (do_kpoints) THEN
            CPABORT("calculate_first_density_matrix: history_guess not implemented for k-points")
         ENDIF
         IF (para_env%ionode) &
            CALL wfn_restart_file_name(file_name, exist, dft_section, logger)
         CALL mp_bcast(exist, para_env%source, para_env%group)
         CALL mp_bcast(file_name, para_env%source, para_env%group)
         nvec = qs_env%wf_history%memory_depth
         not_read = nvec + 1
         ! At this level we read the saved backup RESTART files..
         DO i = 1, nvec
            j = i - 1
            filename = TRIM(file_name)
            IF (j /= 0) filename = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(j))
            IF (para_env%ionode) &
               INQUIRE (FILE=filename, exist=exist)
            CALL mp_bcast(exist, para_env%source, para_env%group)
            IF ((.NOT. exist) .AND. (i < not_read)) THEN
               not_read = i
            END IF
         END DO
         IF (not_read == 1) THEN
            density_guess = restart_guess
            filename = TRIM(file_name)
            IF (para_env%ionode) INQUIRE (FILE=filename, exist=exist)
            CALL mp_bcast(exist, para_env%source, para_env%group)
            IF (.NOT. exist) THEN
               CALL cp_warn(__LOCATION__, &
                            "User requested to restart the wavefunction from a series of restart files named: "// &
                            TRIM(file_name)//" with extensions (.bak-n). These files do not exist."// &
                            " Even trying to switch to a plain restart wave-function failes because the"// &
                            " file named: "//TRIM(file_name)//" does not exist. Please check the existence of"// &
                            " the file or change properly the value of the keyword WFN_RESTART_FILE_NAME. "// &
                            " Calculation continues using ATOMIC GUESS. ")
               density_guess = safe_density_guess
            END IF
         END IF
         last_read = not_read - 1
      END IF

      did_guess = .FALSE.

      IF (dft_control%correct_el_density_dip) THEN
         tot_corr_zeff = qs_env%total_zeff_corr
         !WRITE(*,*) "tot_corr_zeff = ", tot_corr_zeff
         IF ((ABS(tot_corr_zeff) > 0.0_dp) .AND. (density_guess /= restart_guess)) THEN
            CALL cp_warn(__LOCATION__, &
                         "Use SCF_GUESS RESTART in conjunction with "// &
                         "CORE_CORRECTION /= 0.0 and SURFACE_DIPOLE_CORRECTION TRUE. "// &
                         "It is always advisable to perform SURFACE_DIPOLE_CORRECTION "// &
                         "after a simulation without the surface dipole correction "// &
                         "and using the ensuing wavefunction restart file. ")
         ENDIF
      ENDIF

      IF (density_guess == restart_guess) THEN

         IF (do_kpoints) THEN
            natoms = SIZE(particle_set)
            CALL read_kpoints_restart(rho_ao_kp, kpoints, work1, &
                                      natoms, para_env, id_nr, dft_section, natom_mismatch)
            IF (natom_mismatch) density_guess = safe_density_guess
         ELSE
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
            CALL read_mo_set_from_restart(mo_array, atomic_kind_set, qs_kind_set, particle_set, para_env, &
                                          id_nr=id_nr, multiplicity=dft_control%multiplicity, dft_section=dft_section, &
                                          natom_mismatch=natom_mismatch)

            IF (natom_mismatch) THEN
               density_guess = safe_density_guess
            ELSE
               DO ispin = 1, nspin
                  IF (scf_control%level_shift /= 0.0_dp) THEN
                     CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, mo_coeff=mo_coeff)
                     CALL cp_fm_to_fm(mo_coeff, ortho)
                  END IF

                  ! make all nmo vectors present orthonormal
                  CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                                  mo_coeff=mo_coeff, nmo=nmo, homo=homo)

                  IF (has_unit_metric) THEN
                     CALL make_basis_simple(mo_coeff, nmo)
                  ELSEIF (dft_control%smear) THEN
                     CALL make_basis_lowdin(vmatrix=mo_coeff, ncol=nmo, &
                                            matrix_s=s_sparse(1)%matrix)
                  ELSE
                     ! ortho so that one can restart for different positions (basis sets?)
                     CALL make_basis_sm(mo_coeff, homo, s_sparse(1)%matrix)
                  ENDIF
                  ! only alpha spin is kept for restricted
                  IF (dft_control%restricted) EXIT
               ENDDO
               IF (dft_control%restricted) CALL mo_set_restrict(mo_array)

               IF (.NOT. scf_control%diagonalization%mom) THEN
                  IF (dft_control%correct_surf_dip) THEN
                     IF (ABS(tot_corr_zeff) > 0.0_dp) THEN
                        CALL set_mo_occupation(mo_array, smear=qs_env%scf_control%smear, &
                                               tot_zeff_corr=tot_corr_zeff)
                     ELSE
                        CALL set_mo_occupation(mo_array, smear=qs_env%scf_control%smear)
                     ENDIF
                  ELSE
                     CALL set_mo_occupation(mo_array, smear=qs_env%scf_control%smear)
                  ENDIF
               ENDIF

               DO ispin = 1, nspin

                  IF (scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr) THEN !fm->dbcsr
                     CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff, &
                                           mo_array(ispin)%mo_set%mo_coeff_b) !fm->dbcsr
                  ENDIF !fm->dbcsr

                  CALL calculate_density_matrix(mo_array(ispin)%mo_set, &
                                                p_rmpv(ispin)%matrix)
               ENDDO
            ENDIF ! natom_mismatch

         END IF

         ! Maximum Overlap Method
         IF (scf_control%diagonalization%mom) THEN
            CALL do_mom_guess(nspin, mo_array, scf_control, p_rmpv)
         ENDIF

         did_guess = .TRUE.
      END IF

      IF (density_guess == history_guess) THEN
         IF (not_read > 1) THEN
            DO i = 1, last_read
               j = last_read - i
               CALL read_mo_set_from_restart(mo_array, atomic_kind_set, qs_kind_set, particle_set, para_env, &
                                             id_nr=j, multiplicity=dft_control%multiplicity, &
                                             dft_section=dft_section)

               DO ispin = 1, nspin
                  IF (scf_control%level_shift /= 0.0_dp) THEN
                     CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, mo_coeff=mo_coeff)
                     CALL cp_fm_to_fm(mo_coeff, ortho)
                  END IF

                  ! make all nmo vectors present orthonormal
                  CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, mo_coeff=mo_coeff, nmo=nmo, homo=homo)

                  IF (has_unit_metric) THEN
                     CALL make_basis_simple(mo_coeff, nmo)
                  ELSE
                     ! ortho so that one can restart for different positions (basis sets?)
                     CALL make_basis_sm(mo_coeff, homo, s_sparse(1)%matrix)
                  ENDIF
                  ! only alpha spin is kept for restricted
                  IF (dft_control%restricted) EXIT
               END DO
               IF (dft_control%restricted) CALL mo_set_restrict(mo_array)

               DO ispin = 1, nspin
                  CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set, &
                                         smear=qs_env%scf_control%smear)
               ENDDO

               DO ispin = 1, nspin
                  IF (scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr) THEN !fm->dbcsr
                     CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff, &
                                           mo_array(ispin)%mo_set%mo_coeff_b) !fm->dbcsr
                  END IF !fm->dbcsr
                  CALL calculate_density_matrix(mo_array(ispin)%mo_set, p_rmpv(ispin)%matrix)
               ENDDO

               ! Write to extrapolation pipeline
               CALL wfi_update(wf_history=qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
            END DO
         END IF

         did_guess = .TRUE.
      END IF

      IF (density_guess == random_guess) THEN

         DO ispin = 1, nspin
            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                            mo_coeff=mo_coeff, nmo=nmo)
            CALL cp_fm_init_random(mo_coeff, nmo)
            IF (has_unit_metric) THEN
               CALL make_basis_simple(mo_coeff, nmo)
            ELSE
               CALL make_basis_sm(mo_coeff, nmo, s_sparse(1)%matrix)
            ENDIF
            ! only alpha spin is kept for restricted
            IF (dft_control%restricted) EXIT
         ENDDO
         IF (dft_control%restricted) CALL mo_set_restrict(mo_array)

         DO ispin = 1, nspin
            CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set, &
                                   smear=qs_env%scf_control%smear)
         ENDDO

         DO ispin = 1, nspin

            IF (scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr) THEN !fm->dbcsr
               CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff, &
                                     mo_array(ispin)%mo_set%mo_coeff_b) !fm->dbcsr
            ENDIF !fm->dbcsr

            CALL calculate_density_matrix(mo_array(ispin)%mo_set, p_rmpv(ispin)%matrix)
         ENDDO

         did_guess = .TRUE.
      END IF

      IF (density_guess == core_guess) THEN

         IF (do_kpoints) THEN
            CPABORT("calculate_first_density_matrix: core_guess not implemented for k-points")
         ENDIF

         IF (.NOT. ASSOCIATED(work1)) THEN
            need_wm = .TRUE.
            CPASSERT(.NOT. ASSOCIATED(work2))
            CPASSERT(.NOT. ASSOCIATED(ortho))
         ELSE
            need_wm = .FALSE.
            CPASSERT(ASSOCIATED(work2))
         END IF

         IF (need_wm) THEN
            CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff=moa)
            CALL cp_fm_get_info(moa, matrix_struct=ao_mo_struct)
            CALL cp_fm_struct_get(ao_mo_struct, nrow_global=nao, nrow_block=nblocks)
            CALL cp_fm_struct_create(fmstruct=ao_ao_struct, &
                                     nrow_block=nblocks, &
                                     ncol_block=nblocks, &
                                     nrow_global=nao, &
                                     ncol_global=nao, &
                                     template_fmstruct=ao_mo_struct)
            ALLOCATE (work1(1))
            NULLIFY (work1(1)%matrix)
            CALL cp_fm_create(work1(1)%matrix, ao_ao_struct)
            CALL cp_fm_create(work2, ao_ao_struct)
            CALL cp_fm_create(ortho, ao_ao_struct)
            CALL copy_dbcsr_to_fm(matrix_s_kp(1, 1)%matrix, ortho)
            CALL cp_fm_cholesky_decompose(ortho)
            CALL cp_fm_struct_release(ao_ao_struct)
         END IF

         ispin = 1
         ! Load core Hamiltonian into work matrix
         CALL copy_dbcsr_to_fm(h_core_sparse(1)%matrix, work1(ispin)%matrix)

         ! Diagonalize the core Hamiltonian matrix and retrieve a first set of
         ! molecular orbitals (MOs)
         IF (has_unit_metric) THEN
            CALL eigensolver_simple(matrix_ks=work1(ispin)%matrix, &
                                    mo_set=mo_array(ispin)%mo_set, &
                                    work=work2, &
                                    do_level_shift=.FALSE., &
                                    level_shift=0.0_dp, &
                                    use_jacobi=.FALSE., jacobi_threshold=0._dp)
         ELSE
            CALL eigensolver(matrix_ks_fm=work1(ispin)%matrix, &
                             mo_set=mo_array(ispin)%mo_set, &
                             ortho=ortho, &
                             work=work2, &
                             cholesky_method=scf_env%cholesky_method, &
                             do_level_shift=.FALSE., &
                             level_shift=0.0_dp, &
                             matrix_u_fm=null(), &
                             use_jacobi=.FALSE.)
         END IF

         ! Open shell case: copy alpha MOs to beta MOs
         IF (nspin == 2) THEN
            CALL get_mo_set(mo_set=mo_array(1)%mo_set, mo_coeff=moa)
            CALL get_mo_set(mo_set=mo_array(2)%mo_set, mo_coeff=mob, nmo=nmo)
            CALL cp_fm_to_fm(moa, mob, nmo)
         END IF

         ! Build an initial density matrix (for each spin in the case of
         ! an open shell calculation) from the first MOs set
         DO ispin = 1, nspin
            CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set, smear=scf_control%smear)
            CALL calculate_density_matrix(mo_array(ispin)%mo_set, p_rmpv(ispin)%matrix)
         END DO

         ! release intermediate matrices
         IF (need_wm) THEN
            CALL cp_fm_release(ortho)
            CALL cp_fm_release(work2)
            CALL cp_fm_release(work1(1)%matrix)
            DEALLOCATE (work1)
            NULLIFY (work1, work2, ortho)
         END IF

         did_guess = .TRUE.
      END IF

      IF (density_guess == atomic_guess) THEN

         subsys_section => section_vals_get_subs_vals(input, "SUBSYS")
         ounit = cp_print_key_unit_nr(logger, subsys_section, "PRINT%KINDS", extension=".Log")
         IF (ounit > 0) THEN
            WRITE (UNIT=ounit, FMT="(/,(T2,A))") &
               "Atomic guess: The first density matrix is obtained in terms of atomic orbitals", &
               "              and electronic configurations assigned to each atomic kind"
         END IF

         CALL calculate_atomic_block_dm(p_rmpv, s_sparse(1)%matrix, &
                                        particle_set, atomic_kind_set, qs_kind_set, &
                                        nspin, nelectron_spin, ounit, para_env)

         DO ispin = 1, nspin

            ! The orbital transformation method (OT) requires not only an
            ! initial density matrix, but also an initial wavefunction (MO set)
            IF (ofgpw .AND. (scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr)) THEN
               ! get orbitals later
            ELSE
               IF (need_mos) THEN

                  IF (dft_control%restricted .AND. (ispin == 2)) THEN
                     CALL mo_set_restrict(mo_array)
                  ELSE
                     CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                                     mo_coeff=mo_coeff, &
                                     nmo=nmo, nao=nao, homo=homo)

                     CALL cp_fm_set_all(mo_coeff, 0.0_dp)
                     CALL cp_fm_init_random(mo_coeff, nmo)

                     CALL cp_fm_create(sv, mo_coeff%matrix_struct, "SV")
                     ! multiply times PS
                     IF (has_unit_metric) THEN
                        CALL cp_fm_to_fm(mo_coeff, sv)
                     ELSE
                        ! PS*C(:,1:nomo)+C(:,nomo+1:nmo) (nomo=NINT(nelectron/maxocc))
                        CALL cp_dbcsr_sm_fm_multiply(s_sparse(1)%matrix, mo_coeff, sv, nmo)
                     END IF
                     CALL cp_dbcsr_sm_fm_multiply(p_rmpv(ispin)%matrix, sv, mo_coeff, homo)

                     CALL cp_fm_release(sv)
                     ! and ortho the result
                     IF (has_unit_metric) THEN
                        CALL make_basis_simple(mo_coeff, nmo)
                     ELSE
                        CALL make_basis_sm(mo_coeff, nmo, s_sparse(1)%matrix)
                     END IF
                  END IF

                  CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set, &
                                         smear=qs_env%scf_control%smear)

                  CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff, &
                                        mo_array(ispin)%mo_set%mo_coeff_b) !fm->dbcsr

                  CALL calculate_density_matrix(mo_array(ispin)%mo_set, &
                                                p_rmpv(ispin)%matrix)
               END IF
               ! adjust el_density in case surface_dipole_correction is switched
               ! on and CORE_CORRECTION is non-zero
               IF (scf_env%method == general_diag_method_nr) THEN
                  IF (dft_control%correct_surf_dip) THEN
                     IF (ABS(tot_corr_zeff) > 0.0_dp) THEN
                        CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                                        mo_coeff=mo_coeff, &
                                        nmo=nmo, nao=nao, homo=homo)

                        CALL cp_fm_set_all(mo_coeff, 0.0_dp)
                        CALL cp_fm_init_random(mo_coeff, nmo)

                        CALL cp_fm_create(sv, mo_coeff%matrix_struct, "SV")
                        ! multiply times PS
                        IF (has_unit_metric) THEN
                           CALL cp_fm_to_fm(mo_coeff, sv)
                        ELSE
                           ! PS*C(:,1:nomo)+C(:,nomo+1:nmo) (nomo=NINT(nelectron/maxocc))
                           CALL cp_dbcsr_sm_fm_multiply(s_sparse(1)%matrix, mo_coeff, sv, nmo)
                        END IF
                        CALL cp_dbcsr_sm_fm_multiply(p_rmpv(ispin)%matrix, sv, mo_coeff, homo)

                        CALL cp_fm_release(sv)
                        ! and ortho the result
                        IF (has_unit_metric) THEN
                           CALL make_basis_simple(mo_coeff, nmo)
                        ELSE
                           CALL make_basis_sm(mo_coeff, nmo, s_sparse(1)%matrix)
                        END IF

                        CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set, smear=qs_env%scf_control%smear, &
                                               tot_zeff_corr=tot_corr_zeff)

                        CALL calculate_density_matrix(mo_array(ispin)%mo_set, &
                                                      p_rmpv(ispin)%matrix)
                     ENDIF
                  ENDIF
               ENDIF

            END IF

         END DO

         IF (ofgpw .AND. (scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr)) THEN
            ! We fit a function to the square root of the density
            CALL qs_rho_update_rho(rho, qs_env)
            CPASSERT(1 == 0)
!         CALL cp_fm_create(sv,mo_coeff%matrix_struct,"SV")
!         DO ispin=1,nspin
!           CALL integrate_ppl_rspace(qs%rho%rho_r(ispin),qs_env)
!           CALL cp_cfm_solve(overlap,mos)
!           CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,&
!                           mo_coeff=mo_coeff, nmo=nmo, nao=nao)
!           CALL cp_fm_init_random(mo_coeff,nmo)
!         END DO
!         CALL cp_fm_release(sv)
         END IF

         IF (scf_control%diagonalization%mom) THEN
            CALL do_mom_guess(nspin, mo_array, scf_control, p_rmpv)
         ENDIF

         CALL cp_print_key_finished_output(ounit, logger, subsys_section, &
                                           "PRINT%KINDS")

         did_guess = .TRUE.
      END IF

      IF (density_guess == sparse_guess) THEN

         IF (ofgpw) CPABORT("SPARSE_GUESS not implemented for OFGPW")
         IF (.NOT. scf_control%use_ot) CPABORT("OT needed!")
         IF (do_kpoints) THEN
            CPABORT("calculate_first_density_matrix: sparse_guess not implemented for k-points")
         ENDIF

         eps = 1.0E-5_dp

         ounit = cp_logger_get_default_io_unit(logger)
         group = para_env%group
         natoms = SIZE(particle_set)
         ALLOCATE (kind_of(natoms))
         ALLOCATE (first_sgf(natoms), last_sgf(natoms))

         checksum = dbcsr_checksum(s_sparse(1)%matrix)
         i = dbcsr_get_num_blocks(s_sparse(1)%matrix); CALL mp_sum(i, group)
         IF (ounit > 0) WRITE (ounit, *) 'S nblks', i, ' checksum', checksum
         CALL dbcsr_filter(s_sparse(1)%matrix, eps)
         checksum = dbcsr_checksum(s_sparse(1)%matrix)
         i = dbcsr_get_num_blocks(s_sparse(1)%matrix); CALL mp_sum(i, group)
         IF (ounit > 0) WRITE (ounit, *) 'S nblks', i, ' checksum', checksum

         CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, &
                               last_sgf=last_sgf)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

         ALLOCATE (pmat(SIZE(atomic_kind_set)))

         rscale = 1._dp
         IF (nspin == 2) rscale = 0.5_dp
         DO ikind = 1, SIZE(atomic_kind_set)
            atomic_kind => atomic_kind_set(ikind)
            qs_kind => qs_kind_set(ikind)
            NULLIFY (pmat(ikind)%mat)
            CALL calculate_atomic_orbitals(atomic_kind, qs_kind, pmat=pmat(ikind)%mat)
            NULLIFY (atomic_kind)
         END DO

         DO ispin = 1, nspin
            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                            maxocc=maxocc, &
                            nelectron=nelectron)
            !
            CALL dbcsr_iterator_start(iter, p_rmpv(ispin)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, pdata, blk)
               ikind = kind_of(irow)
               IF (icol .EQ. irow) THEN
                  IF (ispin == 1) THEN
                     pdata(:, :) = pmat(ikind)%mat(:, :, 1)*rscale + &
                                   pmat(ikind)%mat(:, :, 2)*rscale
                  ELSE
                     pdata(:, :) = pmat(ikind)%mat(:, :, 1)*rscale - &
                                   pmat(ikind)%mat(:, :, 2)*rscale
                  END IF
               END IF
            ENDDO
            CALL dbcsr_iterator_stop(iter)

            !CALL dbcsr_verify_matrix(p_rmpv(ispin)%matrix)
            checksum = dbcsr_checksum(p_rmpv(ispin)%matrix)
            occ = dbcsr_get_occupation(p_rmpv(ispin)%matrix)
            IF (ounit > 0) WRITE (ounit, *) 'P_init occ', occ, ' checksum', checksum
            ! so far p needs to have the same sparsity as S
            !CALL dbcsr_filter(p_rmpv(ispin)%matrix, eps)
            !CALL dbcsr_verify_matrix(p_rmpv(ispin)%matrix)
            checksum = dbcsr_checksum(p_rmpv(ispin)%matrix)
            occ = dbcsr_get_occupation(p_rmpv(ispin)%matrix)
            IF (ounit > 0) WRITE (ounit, *) 'P_init occ', occ, ' checksum', checksum

            CALL dbcsr_dot(p_rmpv(ispin)%matrix, s_sparse(1)%matrix, trps1)
            rscale = REAL(nelectron, dp)/trps1
            CALL dbcsr_scale(p_rmpv(ispin)%matrix, rscale)

            !CALL dbcsr_verify_matrix(p_rmpv(ispin)%matrix)
            checksum = dbcsr_checksum(p_rmpv(ispin)%matrix)
            occ = dbcsr_get_occupation(p_rmpv(ispin)%matrix)
            IF (ounit > 0) WRITE (ounit, *) 'P occ', occ, ' checksum', checksum
            !
            ! The orbital transformation method (OT) requires not only an
            ! initial density matrix, but also an initial wavefunction (MO set)
            IF (dft_control%restricted .AND. (ispin == 2)) THEN
               CALL mo_set_restrict(mo_array)
            ELSE
               CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                               mo_coeff=mo_coeff, &
                               nmo=nmo, nao=nao, homo=homo)
               CALL cp_fm_set_all(mo_coeff, 0.0_dp)

               n = MAXVAL(last_sgf - first_sgf) + 1
               size_atomic_kind_set = SIZE(atomic_kind_set)

               ALLOCATE (buff(n, n), sort_kind(size_atomic_kind_set), &
                         nelec_kind(size_atomic_kind_set))
               !
               ! sort kind vs nbr electron
               DO ikind = 1, size_atomic_kind_set
                  atomic_kind => atomic_kind_set(ikind)
                  qs_kind => qs_kind_set(ikind)
                  CALL get_atomic_kind(atomic_kind=atomic_kind, &
                                       natom=natom, &
                                       atom_list=atom_list, &
                                       z=z)
                  CALL get_qs_kind(qs_kind, nsgf=nsgf, elec_conf=elec_conf, &
                                   basis_set=orb_basis_set, zeff=zeff)
                  nelec_kind(ikind) = SUM(elec_conf)
               ENDDO
               CALL sort(nelec_kind, size_atomic_kind_set, sort_kind)
               !
               ! a -very- naive sparse guess
               nmo_tmp = nmo
               natoms_tmp = natoms
               istart_col = 1
               iseed(1) = 4; iseed(2) = 3; iseed(3) = 2; iseed(4) = 1 ! set the seed for dlarnv
               DO i = 1, size_atomic_kind_set
                  ikind = sort_kind(i)
                  atomic_kind => atomic_kind_set(ikind)
                  CALL get_atomic_kind(atomic_kind=atomic_kind, &
                                       natom=natom, atom_list=atom_list)
                  DO iatom = 1, natom
                     !
                     atom_a = atom_list(iatom)
                     istart_row = first_sgf(atom_a)
                     n_rows = last_sgf(atom_a) - first_sgf(atom_a) + 1
                     !
                     ! compute the "potential" nbr of states for this atom
                     n_cols = MAX(INT(REAL(nmo_tmp, dp)/REAL(natoms_tmp, dp)), 1)
                     IF (n_cols .GT. n_rows) n_cols = n_rows
                     !
                     nmo_tmp = nmo_tmp - n_cols
                     natoms_tmp = natoms_tmp - 1
                     IF (nmo_tmp .LT. 0 .OR. natoms_tmp .LT. 0) THEN
                        CPABORT("Wrong1!")
                     END IF
                     DO j = 1, n_cols
                        CALL dlarnv(1, iseed, n_rows, buff(1, j))
                     ENDDO
                     CALL cp_fm_set_submatrix(mo_coeff, buff, istart_row, istart_col, &
                                              n_rows, n_cols)
                     istart_col = istart_col + n_cols
                  ENDDO
               ENDDO

               IF (istart_col .LE. nmo) THEN
                  CPABORT("Wrong2!")
               END IF

               DEALLOCATE (buff, nelec_kind, sort_kind)

               IF (.FALSE.) THEN
                  ALLOCATE (buff(nao, 1), buff2(nao, 1))
                  DO i = 1, nmo
                     CALL cp_fm_get_submatrix(mo_coeff, buff, 1, i, nao, 1)
                     IF (SUM(buff**2) .LT. 1E-10_dp) THEN
                        WRITE (*, *) 'wrong', i, SUM(buff**2)
                     ENDIF
                     length = SQRT(DOT_PRODUCT(buff(:, 1), buff(:, 1)))
                     buff(:, :) = buff(:, :)/length
                     DO j = i + 1, nmo
                        CALL cp_fm_get_submatrix(mo_coeff, buff2, 1, j, nao, 1)
                        length = SQRT(DOT_PRODUCT(buff2(:, 1), buff2(:, 1)))
                        buff2(:, :) = buff2(:, :)/length
                        IF (ABS(DOT_PRODUCT(buff(:, 1), buff2(:, 1)) - 1.0_dp) .LT. 1E-10_dp) THEN
                           WRITE (*, *) 'wrong2', i, j, DOT_PRODUCT(buff(:, 1), buff2(:, 1))
                           DO ikind = 1, nao
                              IF (ABS(mo_coeff%local_data(ikind, i)) .GT. 1e-10_dp) THEN
                                 WRITE (*, *) 'c1', ikind, mo_coeff%local_data(ikind, i)
                              ENDIF
                              IF (ABS(mo_coeff%local_data(ikind, j)) .GT. 1e-10_dp) THEN
                                 WRITE (*, *) 'c2', ikind, mo_coeff%local_data(ikind, j)
                              ENDIF
                           ENDDO
                           CPABORT("")
                        ENDIF
                     ENDDO
                  ENDDO
                  DEALLOCATE (buff, buff2)

               ENDIF
               !
               CALL cp_fm_to_dbcsr_row_template(mo_dbcsr, mo_coeff, s_sparse(1)%matrix)
               !CALL dbcsr_verify_matrix(mo_dbcsr)
               checksum = dbcsr_checksum(mo_dbcsr)

               occ = dbcsr_get_occupation(mo_dbcsr)
               IF (ounit > 0) WRITE (ounit, *) 'C occ', occ, ' checksum', checksum
               CALL dbcsr_filter(mo_dbcsr, eps)
               !CALL dbcsr_verify_matrix(mo_dbcsr)
               occ = dbcsr_get_occupation(mo_dbcsr)
               checksum = dbcsr_checksum(mo_dbcsr)
               IF (ounit > 0) WRITE (ounit, *) 'C occ', occ, ' checksum', checksum
               !
               ! multiply times PS
               IF (has_unit_metric) THEN
                  CPABORT("has_unit_metric will be removed soon")
               END IF
               !
               ! S*C
               CALL dbcsr_copy(mo_tmp_dbcsr, mo_dbcsr, name="mo_tmp")
               CALL dbcsr_multiply("N", "N", 1.0_dp, s_sparse(1)%matrix, mo_dbcsr, &
                                   0.0_dp, mo_tmp_dbcsr, &
                                   retain_sparsity=.TRUE.)
               !CALL dbcsr_verify_matrix(mo_tmp_dbcsr)
               checksum = dbcsr_checksum(mo_tmp_dbcsr)
               occ = dbcsr_get_occupation(mo_tmp_dbcsr)
               IF (ounit > 0) WRITE (ounit, *) 'S*C occ', occ, ' checksum', checksum
               CALL dbcsr_filter(mo_tmp_dbcsr, eps)
               !CALL dbcsr_verify_matrix(mo_tmp_dbcsr)
               checksum = dbcsr_checksum(mo_tmp_dbcsr)
               occ = dbcsr_get_occupation(mo_tmp_dbcsr)
               IF (ounit > 0) WRITE (ounit, *) 'S*C occ', occ, ' checksum', checksum
               !
               ! P*SC
               ! the destroy is needed for the moment to avoid memory leaks !
               ! This one is not needed because _destroy takes care of zeroing.
               CALL dbcsr_multiply("N", "N", 1.0_dp, p_rmpv(ispin)%matrix, &
                                   mo_tmp_dbcsr, 0.0_dp, mo_dbcsr)
               IF (.FALSE.) CALL dbcsr_verify_matrix(mo_dbcsr)
               checksum = dbcsr_checksum(mo_dbcsr)
               occ = dbcsr_get_occupation(mo_dbcsr)
               IF (ounit > 0) WRITE (ounit, *) 'P*SC occ', occ, ' checksum', checksum
               CALL dbcsr_filter(mo_dbcsr, eps)
               !CALL dbcsr_verify_matrix(mo_dbcsr)
               checksum = dbcsr_checksum(mo_dbcsr)
               occ = dbcsr_get_occupation(mo_dbcsr)
               IF (ounit > 0) WRITE (ounit, *) 'P*SC occ', occ, ' checksum', checksum
               !
               CALL copy_dbcsr_to_fm(mo_dbcsr, mo_coeff)

               CALL dbcsr_release(mo_dbcsr)
               CALL dbcsr_release(mo_tmp_dbcsr)

               ! and ortho the result
               CALL make_basis_sm(mo_coeff, nmo, s_sparse(1)%matrix)
            END IF

            CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set, &
                                   smear=qs_env%scf_control%smear)

            CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff, &
                                  mo_array(ispin)%mo_set%mo_coeff_b) !fm->dbcsr

            CALL calculate_density_matrix(mo_array(ispin)%mo_set, &
                                          p_rmpv(ispin)%matrix)
            DO ikind = 1, SIZE(atomic_kind_set)
               IF (ASSOCIATED(pmat(ikind)%mat)) THEN
                  DEALLOCATE (pmat(ikind)%mat)
               END IF
            END DO
         END DO

         DEALLOCATE (pmat)

         DEALLOCATE (kind_of)

         DEALLOCATE (first_sgf, last_sgf)

         did_guess = .TRUE.
      END IF
      IF (density_guess == mopac_guess) THEN

         CALL calculate_mopac_dm(p_rmpv, s_sparse(1)%matrix, has_unit_metric, dft_control, &
                                 particle_set, atomic_kind_set, qs_kind_set, nspin, nelectron_spin, para_env)

         DO ispin = 1, nspin
            ! The orbital transformation method (OT) requires not only an
            ! initial density matrix, but also an initial wavefunction (MO set)
            IF (need_mos) THEN
               IF (dft_control%restricted .AND. (ispin == 2)) THEN
                  CALL mo_set_restrict(mo_array)
               ELSE
                  CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                                  mo_coeff=mo_coeff, &
                                  nmo=nmo, homo=homo)
                  CALL cp_fm_init_random(mo_coeff, nmo)
                  CALL cp_fm_create(sv, mo_coeff%matrix_struct, "SV")
                  ! multiply times PS
                  IF (has_unit_metric) THEN
                     CALL cp_fm_to_fm(mo_coeff, sv)
                  ELSE
                     CALL cp_dbcsr_sm_fm_multiply(s_sparse(1)%matrix, mo_coeff, sv, nmo)
                  END IF
                  ! here we could easily multiply with the diag that we actually have replicated already
                  CALL cp_dbcsr_sm_fm_multiply(p_rmpv(ispin)%matrix, sv, mo_coeff, homo)
                  CALL cp_fm_release(sv)
                  ! and ortho the result
                  IF (has_unit_metric) THEN
                     CALL make_basis_simple(mo_coeff, nmo)
                  ELSE
                     CALL make_basis_sm(mo_coeff, nmo, s_sparse(1)%matrix)
                  END IF
               END IF

               CALL set_mo_occupation(mo_set=mo_array(ispin)%mo_set, &
                                      smear=qs_env%scf_control%smear)
               CALL copy_fm_to_dbcsr(mo_array(ispin)%mo_set%mo_coeff, &
                                     mo_array(ispin)%mo_set%mo_coeff_b)

               CALL calculate_density_matrix(mo_array(ispin)%mo_set, &
                                             p_rmpv(ispin)%matrix)
            END IF
         END DO

         did_guess = .TRUE.
      END IF
! switch_surf_dip [SGh]
      IF (dft_control%switch_surf_dip) THEN
         DO ispin = 1, nspin
            CALL reassign_allocated_mos(mos_last_converged(ispin)%mo_set, &
                                        mo_array(ispin)%mo_set)
         ENDDO
      ENDIF

      IF (density_guess == no_guess) THEN
         did_guess = .TRUE.
      END IF

      IF (.NOT. did_guess) THEN
         CPABORT("An invalid keyword for the initial density guess was specified")
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_first_density_matrix

! **************************************************************************************************
!> \brief returns a block diagonal density matrix. Blocks correspond to the atomic densities.
!> \param pmatrix ...
!> \param matrix_s ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param nspin ...
!> \param nelectron_spin ...
!> \param ounit ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE calculate_atomic_block_dm(pmatrix, matrix_s, particle_set, atomic_kind_set, &
                                        qs_kind_set, nspin, nelectron_spin, ounit, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: pmatrix
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_s
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: nspin
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nelectron_spin
      INTEGER, INTENT(IN)                                :: ounit
      TYPE(cp_para_env_type)                             :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atomic_block_dm'

      INTEGER                                            :: blk, group, handle, icol, ikind, irow, &
                                                            ispin, natom, nc, nkind, nocc(2)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nok
      REAL(dp), DIMENSION(:, :), POINTER                 :: pdata
      REAL(KIND=dp)                                      :: rds, rscale, trps1
      TYPE(atom_matrix_type), ALLOCATABLE, DIMENSION(:)  :: pmat
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: matrix_p
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      natom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      ALLOCATE (pmat(nkind))
      ALLOCATE (nok(2, nkind))

      ! precompute the atomic blocks corresponding to spherical atoms
      DO ikind = 1, nkind
         atomic_kind => atomic_kind_set(ikind)
         qs_kind => qs_kind_set(ikind)
         NULLIFY (pmat(ikind)%mat)
         IF (ounit > 0) THEN
            WRITE (UNIT=ounit, FMT="(/,T2,A)") &
               "Guess for atomic kind: "//TRIM(atomic_kind%name)
         END IF
         CALL calculate_atomic_orbitals(atomic_kind, qs_kind, iunit=ounit, &
                                        pmat=pmat(ikind)%mat, nocc=nocc)
         nok(1:2, ikind) = nocc(1:2)
      END DO

      rscale = 1.0_dp
      IF (nspin == 2) rscale = 0.5_dp

      DO ispin = 1, nspin
         IF ((ounit > 0) .AND. (nspin > 1)) THEN
            WRITE (UNIT=ounit, FMT="(/,T2,A,I0)") "Spin ", ispin
         END IF

         matrix_p => pmatrix(ispin)%matrix
         CALL dbcsr_set(matrix_p, 0.0_dp)

         nocc(ispin) = 0
         CALL dbcsr_iterator_start(iter, matrix_p)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, pdata, blk)
            ikind = kind_of(irow)
            IF (icol .EQ. irow) THEN
               IF (ispin == 1) THEN
                  pdata(:, :) = pmat(ikind)%mat(:, :, 1)*rscale + &
                                pmat(ikind)%mat(:, :, 2)*rscale
               ELSE
                  pdata(:, :) = pmat(ikind)%mat(:, :, 1)*rscale - &
                                pmat(ikind)%mat(:, :, 2)*rscale
               END IF
               nocc(ispin) = nocc(ispin) + nok(ispin, ikind)
            END IF
         ENDDO
         CALL dbcsr_iterator_stop(iter)

         CALL dbcsr_dot(matrix_p, matrix_s, trps1)
         rds = 0.0_dp
         ! could be a ghost-atoms-only simulation
         IF (nelectron_spin(ispin) > 0) THEN
            rds = REAL(nelectron_spin(ispin), dp)/trps1
         END IF
         CALL dbcsr_scale(matrix_p, rds)

         IF (ounit > 0) THEN
            IF (nspin > 1) THEN
               WRITE (UNIT=ounit, FMT="(T2,A,I1)") &
                  "Re-scaling the density matrix to get the right number of electrons for spin ", ispin
            ELSE
               WRITE (UNIT=ounit, FMT="(T2,A)") &
                  "Re-scaling the density matrix to get the right number of electrons"
            END IF
            WRITE (ounit, '(T19,A,T44,A,T67,A)') "# Electrons", "Trace(P)", "Scaling factor"
            WRITE (ounit, '(T20,I10,T40,F12.3,T67,F14.3)') nelectron_spin(ispin), trps1, rds
         END IF

         IF (nspin > 1) THEN
            group = para_env%group
            CALL mp_sum(nocc, group)
            IF (nelectron_spin(ispin) > nocc(ispin)) THEN
               rds = 0.99_dp
               CALL dbcsr_scale(matrix_p, rds)
               rds = (1.0_dp - rds)*nelectron_spin(ispin)
               CALL dbcsr_get_info(matrix_p, nfullcols_total=nc)
               rds = rds/REAL(nc, KIND=dp)
               CALL dbcsr_add_on_diag(matrix_p, rds)
               IF (ounit > 0) THEN
                  WRITE (UNIT=ounit, FMT="(T4,A,/,T4,A,T59,F20.12)") &
                     "More MOs than initial guess orbitals detected", &
                     "Add constant to diagonal elements ", rds
               END IF
            END IF
         END IF

      END DO

      DO ikind = 1, nkind
         IF (ASSOCIATED(pmat(ikind)%mat)) THEN
            DEALLOCATE (pmat(ikind)%mat)
         END IF
      END DO
      DEALLOCATE (pmat)

      DEALLOCATE (kind_of, nok)

      CALL timestop(handle)

   END SUBROUTINE calculate_atomic_block_dm

! **************************************************************************************************
!> \brief returns a block diagonal fock matrix.
!> \param matrix_f ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param ounit ...
! **************************************************************************************************
   SUBROUTINE calculate_atomic_fock_matrix(matrix_f, particle_set, atomic_kind_set, &
                                           qs_kind_set, ounit)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_f
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: ounit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atomic_fock_matrix'

      INTEGER                                            :: handle, icol, ikind, irow
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(dp), DIMENSION(:, :), POINTER                 :: block
      TYPE(atom_matrix_type), ALLOCATABLE, DIMENSION(:)  :: fmat
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      ALLOCATE (kind_of(SIZE(particle_set)))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      ALLOCATE (fmat(SIZE(atomic_kind_set)))

      ! precompute the atomic blocks for each atomic-kind
      DO ikind = 1, SIZE(atomic_kind_set)
         atomic_kind => atomic_kind_set(ikind)
         qs_kind => qs_kind_set(ikind)
         NULLIFY (fmat(ikind)%mat)
         IF (ounit > 0) WRITE (UNIT=ounit, FMT="(/,T2,A)") &
            "Calculating atomic Fock matrix for atomic kind: "//TRIM(atomic_kind%name)

         !Currently only ispin=1 is supported
         CALL calculate_atomic_orbitals(atomic_kind, qs_kind, iunit=ounit, &
                                        fmat=fmat(ikind)%mat)
      END DO

      ! zero result matrix
      CALL dbcsr_set(matrix_f, 0.0_dp)

      ! copy precomputed blocks onto diagonal of result matrix
      CALL dbcsr_iterator_start(iter, matrix_f)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, block)
         ikind = kind_of(irow)
         IF (icol .EQ. irow) block(:, :) = fmat(ikind)%mat(:, :, 1)
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      ! cleanup
      DO ikind = 1, SIZE(atomic_kind_set)
         DEALLOCATE (fmat(ikind)%mat)
      END DO
      DEALLOCATE (fmat, kind_of)

      CALL timestop(handle)

   END SUBROUTINE calculate_atomic_fock_matrix

! **************************************************************************************************
!> \brief returns a block diagonal density matrix. Blocks correspond to the mopac initial guess.
!> \param pmat ...
!> \param matrix_s ...
!> \param has_unit_metric ...
!> \param dft_control ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param nspin ...
!> \param nelectron_spin ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE calculate_mopac_dm(pmat, matrix_s, has_unit_metric, &
                                 dft_control, particle_set, atomic_kind_set, qs_kind_set, &
                                 nspin, nelectron_spin, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: pmat
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_s
      LOGICAL                                            :: has_unit_metric
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: nspin
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nelectron_spin
      TYPE(cp_para_env_type)                             :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_mopac_dm'

      INTEGER                                            :: atom_a, group, handle, iatom, ikind, &
                                                            iset, isgf, isgfa, ishell, ispin, la, &
                                                            maxl, maxll, na, nao, natom, ncount, &
                                                            nset, nsgf, z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      INTEGER, DIMENSION(25)                             :: laox, naox
      INTEGER, DIMENSION(5)                              :: occupation
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, elec_conf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, l, last_sgfa
      LOGICAL                                            :: has_pot
      REAL(KIND=dp)                                      :: maxocc, my_sum, nelec, occ, paa, rscale, &
                                                            trps1, trps2, yy, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: econf, pdiag, sdiag
      REAL(KIND=dp), DIMENSION(0:3)                      :: edftb
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(dbcsr_type), POINTER                          :: matrix_p
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      DO ispin = 1, nspin
         matrix_p => pmat(ispin)%matrix
         CALL dbcsr_set(matrix_p, 0.0_dp)
      END DO

      group = para_env%group
      natom = SIZE(particle_set)
      nao = dbcsr_nfullrows_total(pmat(1)%matrix)
      IF (nspin == 1) THEN
         maxocc = 2.0_dp
      ELSE
         maxocc = 1.0_dp
      ENDIF

      ALLOCATE (first_sgf(natom))

      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)
      CALL get_qs_kind_set(qs_kind_set, maxlgto=maxl)

      ALLOCATE (econf(0:maxl))

      ALLOCATE (pdiag(nao))
      pdiag(:) = 0.0_dp

      ALLOCATE (sdiag(nao))

      sdiag(:) = 0.0_dp
      IF (has_unit_metric) THEN
         sdiag(:) = 1.0_dp
      ELSE
         CALL dbcsr_get_diag(matrix_s, sdiag)
         CALL mp_sum(sdiag, group)
      END IF

      ncount = 0
      trps1 = 0.0_dp
      trps2 = 0.0_dp
      pdiag(:) = 0.0_dp

      IF (SUM(nelectron_spin) /= 0) THEN
         DO ikind = 1, SIZE(atomic_kind_set)

            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, &
                             all_potential=all_potential, &
                             gth_potential=gth_potential, &
                             sgp_potential=sgp_potential)
            has_pot = ASSOCIATED(all_potential) .OR. ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)

            IF (dft_control%qs_control%dftb) THEN
               CALL get_dftb_atom_param(qs_kind_set(ikind)%dftb_parameter, &
                                        lmax=maxll, occupation=edftb)
               maxll = MIN(maxll, maxl)
               econf(0:maxl) = edftb(0:maxl)
            ELSEIF (dft_control%qs_control%xtb) THEN
               CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
               CALL get_xtb_atom_param(xtb_kind, z=z, natorb=nsgf, nao=naox, lao=laox, occupation=occupation)
            ELSEIF (has_pot) THEN
               CALL get_atomic_kind(atomic_kind_set(ikind), z=z)
               CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf, elec_conf=elec_conf, zeff=zeff)
               maxll = MIN(SIZE(elec_conf) - 1, maxl)
               econf(:) = 0.0_dp
               econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll), dp)
            ELSE
               CYCLE
            END IF

            ! MOPAC TYPE GUESS
            IF (dft_control%qs_control%dftb) THEN
               DO iatom = 1, natom
                  atom_a = atom_list(iatom)
                  isgfa = first_sgf(atom_a)
                  DO la = 0, maxll
                     SELECT CASE (la)
                     CASE (0)
                        pdiag(isgfa) = econf(0)
                     CASE (1)
                        pdiag(isgfa + 1) = econf(1)/3._dp
                        pdiag(isgfa + 2) = econf(1)/3._dp
                        pdiag(isgfa + 3) = econf(1)/3._dp
                     CASE (2)
                        pdiag(isgfa + 4) = econf(2)/5._dp
                        pdiag(isgfa + 5) = econf(2)/5._dp
                        pdiag(isgfa + 6) = econf(2)/5._dp
                        pdiag(isgfa + 7) = econf(2)/5._dp
                        pdiag(isgfa + 8) = econf(2)/5._dp
                     CASE (3)
                        pdiag(isgfa + 9) = econf(3)/7._dp
                        pdiag(isgfa + 10) = econf(3)/7._dp
                        pdiag(isgfa + 11) = econf(3)/7._dp
                        pdiag(isgfa + 12) = econf(3)/7._dp
                        pdiag(isgfa + 13) = econf(3)/7._dp
                        pdiag(isgfa + 14) = econf(3)/7._dp
                        pdiag(isgfa + 15) = econf(3)/7._dp
                     CASE DEFAULT
                        CPABORT("")
                     END SELECT
                  END DO
               END DO
            ELSEIF (dft_control%qs_control%xtb) THEN
               DO iatom = 1, natom
                  atom_a = atom_list(iatom)
                  isgfa = first_sgf(atom_a)
                  IF (z == 1 .AND. nsgf == 2) THEN
                     ! Hydrogen 2s basis
                     pdiag(isgfa) = 1.0_dp
                     pdiag(isgfa + 1) = 0.0_dp
                  ELSE
                     DO isgf = 1, nsgf
                        na = naox(isgf)
                        la = laox(isgf)
                        occ = REAL(occupation(la + 1), dp)/REAL(2*la + 1, dp)
                        pdiag(isgfa + isgf - 1) = occ
                     END DO
                  END IF
               END DO
            ELSEIF (dft_control%qs_control%semi_empirical) THEN
               yy = REAL(dft_control%charge, KIND=dp)/REAL(nao, KIND=dp)
               DO iatom = 1, natom
                  atom_a = atom_list(iatom)
                  isgfa = first_sgf(atom_a)
                  SELECT CASE (nsgf)
                  CASE (1) ! s-basis
                     pdiag(isgfa) = (zeff - yy)*0.5_dp*maxocc
                  CASE (4) ! sp-basis
                     IF (z == 1) THEN
                        ! special case: hydrogen with sp basis
                        pdiag(isgfa) = (zeff - yy)*0.5_dp*maxocc
                        pdiag(isgfa + 1) = 0._dp
                        pdiag(isgfa + 2) = 0._dp
                        pdiag(isgfa + 3) = 0._dp
                     ELSE
                        pdiag(isgfa) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa + 1) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa + 2) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa + 3) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                     END IF
                  CASE (9) ! spd-basis
                     IF (z < 21 .OR. z > 30 .AND. z < 39 .OR. z > 48 .AND. z < 57) THEN
                        !   Main Group Element:  The "d" shell is formally empty.
                        pdiag(isgfa) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa + 1) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa + 2) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa + 3) = (zeff*0.25_dp - yy)*0.5_dp*maxocc
                        pdiag(isgfa + 4) = (-yy)*0.5_dp*maxocc
                        pdiag(isgfa + 5) = (-yy)*0.5_dp*maxocc
                        pdiag(isgfa + 6) = (-yy)*0.5_dp*maxocc
                        pdiag(isgfa + 7) = (-yy)*0.5_dp*maxocc
                        pdiag(isgfa + 8) = (-yy)*0.5_dp*maxocc
                     ELSE IF (z < 99) THEN
                        my_sum = zeff - 9.0_dp*yy
                        !   First, put 2 electrons in the 's' shell
                        pdiag(isgfa) = (MAX(0.0_dp, MIN(my_sum, 2.0_dp)))*0.5_dp*maxocc
                        my_sum = my_sum - 2.0_dp
                        IF (my_sum > 0.0_dp) THEN
                           !   Now put as many electrons as possible into the 'd' shell
                           pdiag(isgfa + 4) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           pdiag(isgfa + 5) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           pdiag(isgfa + 6) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           pdiag(isgfa + 7) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           pdiag(isgfa + 8) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           my_sum = MAX(0.0_dp, my_sum - 10.0_dp)
                           !   Put the remaining electrons in the 'p' shell
                           pdiag(isgfa + 1) = (my_sum/3.0_dp)*0.5_dp*maxocc
                           pdiag(isgfa + 2) = (my_sum/3.0_dp)*0.5_dp*maxocc
                           pdiag(isgfa + 3) = (my_sum/3.0_dp)*0.5_dp*maxocc
                        END IF
                     END IF
                  CASE DEFAULT
                     CPABORT("")
                  END SELECT
               END DO
            ELSE
               CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                      nset=nset, &
                                      nshell=nshell, &
                                      l=l, &
                                      first_sgf=first_sgfa, &
                                      last_sgf=last_sgfa)

               DO iset = 1, nset
                  DO ishell = 1, nshell(iset)
                     la = l(ishell, iset)
                     nelec = maxocc*REAL(2*la + 1, dp)
                     IF (econf(la) > 0.0_dp) THEN
                        IF (econf(la) >= nelec) THEN
                           paa = maxocc
                           econf(la) = econf(la) - nelec
                        ELSE
                           paa = maxocc*econf(la)/nelec
                           econf(la) = 0.0_dp
                           ncount = ncount + NINT(nelec/maxocc)
                        END IF
                        DO isgfa = first_sgfa(ishell, iset), last_sgfa(ishell, iset)
                           DO iatom = 1, natom
                              atom_a = atom_list(iatom)
                              isgf = first_sgf(atom_a) + isgfa - 1
                              pdiag(isgf) = paa
                              IF (paa == maxocc) THEN
                                 trps1 = trps1 + paa*sdiag(isgf)
                              ELSE
                                 trps2 = trps2 + paa*sdiag(isgf)
                              END IF
                           END DO
                        END DO
                     END IF
                  END DO ! ishell
               END DO ! iset
            END IF
         END DO ! ikind

         IF (trps2 == 0.0_dp) THEN
            DO isgf = 1, nao
               IF (sdiag(isgf) > 0.0_dp) pdiag(isgf) = pdiag(isgf)/sdiag(isgf)
            END DO
            DO ispin = 1, nspin
               IF (nelectron_spin(ispin) /= 0) THEN
                  matrix_p => pmat(ispin)%matrix
                  CALL dbcsr_set_diag(matrix_p, pdiag)
               END IF
            END DO
         ELSE
            DO ispin = 1, nspin
               IF (nelectron_spin(ispin) /= 0) THEN
                  rscale = (REAL(nelectron_spin(ispin), dp) - trps1)/trps2
                  DO isgf = 1, nao
                     IF (pdiag(isgf) < maxocc) pdiag(isgf) = rscale*pdiag(isgf)
                  END DO
                  matrix_p => pmat(ispin)%matrix
                  CALL dbcsr_set_diag(matrix_p, pdiag)
                  DO isgf = 1, nao
                     IF (pdiag(isgf) < maxocc) pdiag(isgf) = pdiag(isgf)/rscale
                  END DO
               END IF
            END DO
         END IF
      END IF

      DEALLOCATE (econf)

      DEALLOCATE (first_sgf)

      DEALLOCATE (pdiag)

      DEALLOCATE (sdiag)

      CALL timestop(handle)

   END SUBROUTINE calculate_mopac_dm

END MODULE qs_initial_guess
