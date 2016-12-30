!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to somehow generate an intial guess
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
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_submatrix, cp_fm_init_random, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_checksum, dbcsr_copy, dbcsr_filter, dbcsr_get_diag, dbcsr_get_num_blocks, &
        dbcsr_get_occupation, dbcsr_init, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_nfullrows_total, dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_set, &
        dbcsr_set_diag, dbcsr_trace, dbcsr_type, dbcsr_verify_matrix
   USE external_potential_types,        ONLY: all_potential_type,&
                                              gth_potential_type
   USE input_constants,                 ONLY: atomic_guess,&
                                              core_guess,&
                                              history_guess,&
                                              mopac_guess,&
                                              no_guess,&
                                              random_guess,&
                                              restart_guess,&
                                              sparse_guess
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
   USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_io,                        ONLY: read_mo_set,&
                                              wfn_restart_file_name
   USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                              make_basis_lowdin,&
                                              make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_restrict
   USE qs_mom_methods,                  ONLY: do_mom_guess
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: eigensolver,&
                                              eigensolver_simple
   USE qs_scf_types,                    ONLY: block_davidson_diag_method_nr,&
                                              block_krylov_diag_method_nr,&
                                              ot_diag_method_nr,&
                                              qs_scf_env_type
   USE qs_wf_history_methods,           ONLY: wfi_update
   USE scf_control_types,               ONLY: scf_control_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_initial_guess'

   PUBLIC ::  calculate_first_density_matrix, calculate_atomic_block_dm, calculate_mopac_dm
   PUBLIC ::  calculate_atomic_fock_matrix

   TYPE atom_matrix_type
      REAL(KIND=dp), DIMENSION(:, :), POINTER   :: mat
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
!> \note
!>      badly needs to be split in subroutines each doing one of the possible
!>      schemes
! **************************************************************************************************
   SUBROUTINE calculate_first_density_matrix(scf_env, qs_env)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_first_density_matrix', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: file_name, filename
      INTEGER :: atom_a, blk, density_guess, group, handle, homo, i, iatom, ic, icol, id_nr, &
         ikind, irow, iseed(4), ispin, istart_col, istart_row, j, last_read, n, n_cols, n_rows, &
         nao, natom, natoms, natoms_tmp, nelectron, nmo, nmo_tmp, not_read, nsgf, nspin, nvec, &
         output_unit, qs_env_id, safe_density_guess, size_atomic_kind_set, z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, kind_of, last_sgf
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, elec_conf, nelec_kind, &
                                                            sort_kind
      LOGICAL                                            :: did_guess, do_kpoints, do_std_diag, &
                                                            exist, has_unit_metric, &
                                                            natom_mismatch, ofgpw
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: buff, buff2
      REAL(dp), DIMENSION(:, :), POINTER                 :: DATA
      REAL(KIND=dp)                                      :: checksum, eps, length, maxocc, occ, &
                                                            scale, trps1, zeff
      REAL(KIND=dp), DIMENSION(0:3)                      :: edftb
      TYPE(atom_matrix_type), DIMENSION(:), POINTER      :: pmat
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: work1
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
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
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
               s_sparse, scf_control, dft_control, h_core_sparse, matrix_ks, rho)
      NULLIFY (dft_section, input, subsys_section)
      NULLIFY (matrix_s_kp, matrix_h_kp, matrix_ks_kp, rho_ao_kp)
      NULLIFY (moa, mob)
      NULLIFY (atom_list, elec_conf, kpoints)
      edftb = 0.0_dp

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
                      para_env=para_env)

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

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

      safe_density_guess = atomic_guess
      IF (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb .OR. &
          dft_control%qs_control%scptb) THEN
         IF (density_guess == atomic_guess) density_guess = mopac_guess
         ! in case we need to bail to a safe restart type later on
         safe_density_guess = mopac_guess
      END IF

      IF (scf_control%use_ot .AND. &
          (.NOT. ((density_guess == random_guess) .OR. &
                  (density_guess == atomic_guess) .OR. &
                  (density_guess == mopac_guess) .OR. &
                  (density_guess == sparse_guess) .OR. &
                  (((density_guess == restart_guess) .OR. &
                    (density_guess == history_guess)) .AND. &
                   (scf_control%level_shift == 0.0_dp))))) THEN
         CALL cp_abort(__LOCATION__, &
                       "OT needs GUESS ATOMIC / RANDOM / SPARSE / RESTART / HISTORY RESTART: other options NYI")
      END IF

      ! if a restart was requested, check that the file exists,
      ! if not we fall back to an atomic guess. No kidding, the file name should remain
      ! in sync with read_mo_set
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
         not_read = nvec+1
         ! At this level we read the saved backup RESTART files..
         DO i = 1, nvec
            j = i-1
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
         last_read = not_read-1
      END IF

      did_guess = .FALSE.

      IF (density_guess == restart_guess) THEN

         IF (do_kpoints) THEN
            natoms = SIZE(particle_set)
            CALL read_kpoints_restart(rho_ao_kp, kpoints, work1, &
                                      natoms, para_env, id_nr, dft_section, natom_mismatch)
            IF (natom_mismatch) density_guess = safe_density_guess
         ELSE
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set)
            CALL read_mo_set(mo_array, atomic_kind_set, qs_kind_set, particle_set, para_env, &
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

               IF (.NOT. scf_control%diagonalization%mom) &
                  CALL set_mo_occupation(mo_array, smear=qs_env%scf_control%smear)

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
               j = last_read-i
               CALL read_mo_set(mo_array, atomic_kind_set, qs_kind_set, particle_set, para_env, &
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

         did_guess = .TRUE.
      END IF

      IF (density_guess == atomic_guess) THEN

         subsys_section => section_vals_get_subs_vals(input, "SUBSYS")
         output_unit = cp_print_key_unit_nr(logger, subsys_section, "PRINT%KINDS", extension=".Log")
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T2,A))") &
               "Atomic guess: The first density matrix is obtained in terms of atomic orbitals", &
               "              and electronic configurations assigned to each atomic kind"
         END IF

         DO ispin = 1, nspin

            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, nelectron=nelectron)

            CALL calculate_atomic_block_dm(p_rmpv(ispin)%matrix, s_sparse(1)%matrix, &
                                           particle_set, atomic_kind_set, qs_kind_set, &
                                           ispin, nspin, nelectron, output_unit)

            ! The orbital transformation method (OT) requires not only an
            ! initial density matrix, but also an initial wavefunction (MO set)
            IF (ofgpw .AND. (scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr)) THEN
               ! get orbitals later
            ELSE
               IF (ASSOCIATED(scf_env%krylov_space)) do_std_diag = (scf_env%krylov_space%eps_std_diag > 0.0_dp)
               IF (scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr .OR. &
                   (scf_env%method == block_krylov_diag_method_nr .AND. .NOT. do_std_diag) &
                   .OR. dft_control%do_admm .OR. scf_env%method == block_davidson_diag_method_nr) THEN
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

         CALL cp_print_key_finished_output(output_unit, logger, subsys_section, &
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

         output_unit = cp_logger_get_default_io_unit(logger)
         group = para_env%group
         natoms = SIZE(particle_set)
         ALLOCATE (kind_of(natoms))
         ALLOCATE (first_sgf(natoms), last_sgf(natoms))

         checksum = dbcsr_checksum(s_sparse(1)%matrix)
         i = dbcsr_get_num_blocks(s_sparse(1)%matrix); CALL mp_sum(i, group)
         IF (output_unit > 0) WRITE (output_unit, *) 'S nblks', i, ' checksum', checksum
         CALL dbcsr_filter(s_sparse(1)%matrix, eps)
         checksum = dbcsr_checksum(s_sparse(1)%matrix)
         i = dbcsr_get_num_blocks(s_sparse(1)%matrix); CALL mp_sum(i, group)
         IF (output_unit > 0) WRITE (output_unit, *) 'S nblks', i, ' checksum', checksum

         CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, &
                               last_sgf=last_sgf)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

         ALLOCATE (pmat(SIZE(atomic_kind_set)))

         DO ispin = 1, nspin
            scale = 1._dp
            IF (nspin == 2) scale = 0.5_dp
            DO ikind = 1, SIZE(atomic_kind_set)
               atomic_kind => atomic_kind_set(ikind)
               qs_kind => qs_kind_set(ikind)
               NULLIFY (pmat(ikind)%mat)
               CALL calculate_atomic_orbitals(atomic_kind, qs_kind, pmat=pmat(ikind)%mat, ispin=ispin)
               NULLIFY (atomic_kind)
            END DO
            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                            maxocc=maxocc, &
                            nelectron=nelectron)
            !
            CALL dbcsr_iterator_start(iter, p_rmpv(ispin)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, DATA, blk)
               ikind = kind_of(irow)
               IF (icol .EQ. irow) DATA(:, :) = pmat(ikind)%mat(:, :)*scale
            ENDDO
            CALL dbcsr_iterator_stop(iter)

            !CALL dbcsr_verify_matrix(p_rmpv(ispin)%matrix)
            checksum = dbcsr_checksum(p_rmpv(ispin)%matrix)
            occ = dbcsr_get_occupation(p_rmpv(ispin)%matrix)
            IF (output_unit > 0) WRITE (output_unit, *) 'P_init occ', occ, ' checksum', checksum
            ! so far p needs to have the same sparsity as S
            !CALL dbcsr_filter(p_rmpv(ispin)%matrix, eps)
            !CALL dbcsr_verify_matrix(p_rmpv(ispin)%matrix)
            checksum = dbcsr_checksum(p_rmpv(ispin)%matrix)
            occ = dbcsr_get_occupation(p_rmpv(ispin)%matrix)
            IF (output_unit > 0) WRITE (output_unit, *) 'P_init occ', occ, ' checksum', checksum

            CALL dbcsr_trace(p_rmpv(ispin)%matrix, s_sparse(1)%matrix, trps1)
            scale = REAL(nelectron, dp)/trps1
            CALL dbcsr_scale(p_rmpv(ispin)%matrix, scale)

            !CALL dbcsr_verify_matrix(p_rmpv(ispin)%matrix)
            checksum = dbcsr_checksum(p_rmpv(ispin)%matrix)
            occ = dbcsr_get_occupation(p_rmpv(ispin)%matrix)
            IF (output_unit > 0) WRITE (output_unit, *) 'P occ', occ, ' checksum', checksum
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

               n = MAXVAL(last_sgf-first_sgf)+1
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
                     n_rows = last_sgf(atom_a)-first_sgf(atom_a)+1
                     !
                     ! compute the "potential" nbr of states for this atom
                     n_cols = MAX(INT(REAL(nmo_tmp, dp)/REAL(natoms_tmp, dp)), 1)
                     IF (n_cols .GT. n_rows) n_cols = n_rows
                     !
                     nmo_tmp = nmo_tmp-n_cols
                     natoms_tmp = natoms_tmp-1
                     IF (nmo_tmp .LT. 0 .OR. natoms_tmp .LT. 0) THEN
                        CPABORT("Wrong1!")
                     END IF
                     DO j = 1, n_cols
                        CALL dlarnv(1, iseed, n_rows, buff(1, j))
                     ENDDO
                     CALL cp_fm_set_submatrix(mo_coeff, buff, istart_row, istart_col, &
                                              n_rows, n_cols)
                     istart_col = istart_col+n_cols
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
                     DO j = i+1, nmo
                        CALL cp_fm_get_submatrix(mo_coeff, buff2, 1, j, nao, 1)
                        length = SQRT(DOT_PRODUCT(buff2(:, 1), buff2(:, 1)))
                        buff2(:, :) = buff2(:, :)/length
                        IF (ABS(DOT_PRODUCT(buff(:, 1), buff2(:, 1))-1.0_dp) .LT. 1E-10_dp) THEN
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
               CALL dbcsr_init(mo_dbcsr)
               CALL cp_fm_to_dbcsr_row_template(mo_dbcsr, mo_coeff, s_sparse(1)%matrix)
               !CALL dbcsr_verify_matrix(mo_dbcsr)
               checksum = dbcsr_checksum(mo_dbcsr)

               occ = dbcsr_get_occupation(mo_dbcsr)
               IF (output_unit > 0) WRITE (output_unit, *) 'C occ', occ, ' checksum', checksum
               CALL dbcsr_filter(mo_dbcsr, eps)
               !CALL dbcsr_verify_matrix(mo_dbcsr)
               occ = dbcsr_get_occupation(mo_dbcsr)
               checksum = dbcsr_checksum(mo_dbcsr)
               IF (output_unit > 0) WRITE (output_unit, *) 'C occ', occ, ' checksum', checksum
               !
               ! multiply times PS
               IF (has_unit_metric) THEN
                  CPABORT("has_unit_metric will be removed soon")
               END IF
               !
               ! S*C
               CALL dbcsr_init(mo_tmp_dbcsr)
               CALL dbcsr_copy(mo_tmp_dbcsr, mo_dbcsr, name="mo_tmp")
               CALL dbcsr_multiply("N", "N", 1.0_dp, s_sparse(1)%matrix, mo_dbcsr, &
                                   0.0_dp, mo_tmp_dbcsr, &
                                   retain_sparsity=.TRUE.)
               !CALL dbcsr_verify_matrix(mo_tmp_dbcsr)
               checksum = dbcsr_checksum(mo_tmp_dbcsr)
               occ = dbcsr_get_occupation(mo_tmp_dbcsr)
               IF (output_unit > 0) WRITE (output_unit, *) 'S*C occ', occ, ' checksum', checksum
               CALL dbcsr_filter(mo_tmp_dbcsr, eps)
               !CALL dbcsr_verify_matrix(mo_tmp_dbcsr)
               checksum = dbcsr_checksum(mo_tmp_dbcsr)
               occ = dbcsr_get_occupation(mo_tmp_dbcsr)
               IF (output_unit > 0) WRITE (output_unit, *) 'S*C occ', occ, ' checksum', checksum
               !
               ! P*SC
               ! the destroy is needed for the moment to avoid memory leaks !
               ! This one is not needed because _destroy takes care of zeroing.
               CALL dbcsr_multiply("N", "N", 1.0_dp, p_rmpv(ispin)%matrix, &
                                   mo_tmp_dbcsr, 0.0_dp, mo_dbcsr)
               IF (.FALSE.) CALL dbcsr_verify_matrix(mo_dbcsr)
               checksum = dbcsr_checksum(mo_dbcsr)
               occ = dbcsr_get_occupation(mo_dbcsr)
               IF (output_unit > 0) WRITE (output_unit, *) 'P*SC occ', occ, ' checksum', checksum
               CALL dbcsr_filter(mo_dbcsr, eps)
               !CALL dbcsr_verify_matrix(mo_dbcsr)
               checksum = dbcsr_checksum(mo_dbcsr)
               occ = dbcsr_get_occupation(mo_dbcsr)
               IF (output_unit > 0) WRITE (output_unit, *) 'P*SC occ', occ, ' checksum', checksum
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

         DO ispin = 1, nspin

            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, nelectron=nelectron)

            CALL calculate_mopac_dm(p_rmpv(ispin)%matrix, s_sparse(1)%matrix, has_unit_metric, dft_control, &
                                    particle_set, atomic_kind_set, qs_kind_set, nspin, nelectron, para_env)

            ! The orbital transformation method (OT) requires not only an
            ! initial density matrix, but also an initial wavefunction (MO set)
            IF (ASSOCIATED(scf_env%krylov_space)) do_std_diag = (scf_env%krylov_space%eps_std_diag > 0.0_dp)
            IF (scf_control%use_ot .OR. scf_env%method == ot_diag_method_nr .OR. &
                (scf_env%method == block_krylov_diag_method_nr .AND. .NOT. do_std_diag)) THEN
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
!> \param matrix_p ...
!> \param matrix_s ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param ispin ...
!> \param nspin ...
!> \param nelectron_spin ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE calculate_atomic_block_dm(matrix_p, matrix_s, particle_set, atomic_kind_set, &
                                        qs_kind_set, ispin, nspin, nelectron_spin, output_unit)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_p, matrix_s
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: ispin, nspin, nelectron_spin, output_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atomic_block_dm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, handle, icol, ikind, irow, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(dp), DIMENSION(:, :), POINTER                 :: DATA
      REAL(KIND=dp)                                      :: scale, trps1
      TYPE(atom_matrix_type), ALLOCATABLE, DIMENSION(:)  :: pmat
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      IF ((output_unit > 0) .AND. (nspin > 1)) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A,I0)") "Spin ", ispin
      END IF

      CALL dbcsr_set(matrix_p, 0.0_dp)

      natom = SIZE(particle_set)

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      ALLOCATE (pmat(SIZE(atomic_kind_set)))

      ! precompute the atomic blocks corresponding to spherical atoms
      DO ikind = 1, SIZE(atomic_kind_set)
         atomic_kind => atomic_kind_set(ikind)
         qs_kind => qs_kind_set(ikind)
         NULLIFY (pmat(ikind)%mat)
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
               "Guess for atomic kind: "//TRIM(atomic_kind%name)
         END IF
         CALL calculate_atomic_orbitals(atomic_kind, qs_kind, iunit=output_unit, pmat=pmat(ikind)%mat, &
                                        ispin=ispin)
      END DO

      scale = 1.0_dp
      IF (nspin == 2) scale = 0.5_dp

      CALL dbcsr_iterator_start(iter, matrix_p)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, DATA, blk)
         ikind = kind_of(irow)
         IF (icol .EQ. irow) DATA(:, :) = pmat(ikind)%mat(:, :)*scale
      ENDDO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_trace(matrix_p, matrix_s, trps1)
      scale = 0.0_dp
      IF (nelectron_spin > 0) & ! could be a ghost-atoms-only simulation
         scale = REAL(nelectron_spin, dp)/trps1
      CALL dbcsr_scale(matrix_p, scale)

      IF (output_unit > 0) THEN
         IF (nspin > 1) THEN
            WRITE (UNIT=output_unit, FMT="(T2,A,I1)") &
               "Re-scaling the density matrix to get the right number of electrons for spin ", ispin
         ELSE
            WRITE (UNIT=output_unit, FMT="(T2,A)") &
               "Re-scaling the density matrix to get the right number of electrons"
         END IF
         WRITE (output_unit, '(T19,A,T44,A,T67,A)') "# Electrons", "Trace(P)", "Scaling factor"
         WRITE (output_unit, '(T20,I10,T40,F12.3,T67,F14.3)') nelectron_spin, trps1, scale
      END IF

      DO ikind = 1, SIZE(atomic_kind_set)
         IF (ASSOCIATED(pmat(ikind)%mat)) THEN
            DEALLOCATE (pmat(ikind)%mat)
         END IF
      END DO

      DEALLOCATE (pmat)

      DEALLOCATE (kind_of)

      CALL timestop(handle)

   END SUBROUTINE calculate_atomic_block_dm

! **************************************************************************************************
!> \brief returns a block diagonal fock matrix.
!> \param matrix_f ...
!> \param particle_set ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE calculate_atomic_fock_matrix(matrix_f, particle_set, atomic_kind_set, &
                                           qs_kind_set, output_unit)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_f
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atomic_fock_matrix', &
         routineP = moduleN//':'//routineN

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
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(/,T2,A)") &
            "Calculating atomic Fock matrix for atomic kind: "//TRIM(atomic_kind%name)

         !Currently only ispin=1 is supported
         CALL calculate_atomic_orbitals(atomic_kind, qs_kind, iunit=output_unit, &
                                        fmat=fmat(ikind)%mat, ispin=1)
      END DO

      ! zero result matrix
      CALL dbcsr_set(matrix_f, 0.0_dp)

      ! copy precomputed blocks onto diagonal of result matrix
      CALL dbcsr_iterator_start(iter, matrix_f)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, block)
         ikind = kind_of(irow)
         IF (icol .EQ. irow) block(:, :) = fmat(ikind)%mat(:, :)
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
!> \param matrix_p ...
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
   SUBROUTINE calculate_mopac_dm(matrix_p, matrix_s, has_unit_metric, &
                                 dft_control, particle_set, atomic_kind_set, qs_kind_set, &
                                 nspin, nelectron_spin, para_env)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_p, matrix_s
      LOGICAL                                            :: has_unit_metric
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: nspin, nelectron_spin
      TYPE(cp_para_env_type)                             :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_mopac_dm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, group, handle, iatom, ikind, &
                                                            iset, isgf, isgfa, ishell, la, maxl, &
                                                            maxll, nao, natom, ncount, nset, nsgf, &
                                                            z
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, elec_conf, nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, l, last_sgfa
      REAL(KIND=dp)                                      :: maxocc, my_sum, nelec, paa, scale, &
                                                            trps1, trps2, yy, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: econf, pdiag, sdiag
      REAL(KIND=dp), DIMENSION(0:3)                      :: edftb
      TYPE(all_potential_type), POINTER                  :: all_potential
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      CALL timeset(routineN, handle)

      group = para_env%group
      natom = SIZE(particle_set)
      nao = dbcsr_nfullrows_total(matrix_p)
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

      IF (nelectron_spin /= 0) THEN
         DO ikind = 1, SIZE(atomic_kind_set)

            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, &
                             all_potential=all_potential, &
                             gth_potential=gth_potential)

            IF (dft_control%qs_control%dftb) THEN
               CALL get_dftb_atom_param(qs_kind_set(ikind)%dftb_parameter, &
                                        lmax=maxll, occupation=edftb)
               maxll = MIN(maxll, maxl)
               econf(0:maxl) = edftb(0:maxl)
            ELSEIF (ASSOCIATED(all_potential)) THEN
               CALL get_atomic_kind(atomic_kind_set(ikind), z=z)
               CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf, elec_conf=elec_conf, zeff=zeff)
               maxll = MIN(SIZE(elec_conf)-1, maxl)
               econf(:) = 0.0_dp
               econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll), dp)
            ELSE IF (ASSOCIATED(gth_potential)) THEN
               CALL get_atomic_kind(atomic_kind_set(ikind), z=z)
               CALL get_qs_kind(qs_kind_set(ikind), nsgf=nsgf, elec_conf=elec_conf, zeff=zeff)
               maxll = MIN(SIZE(elec_conf)-1, maxl)
               econf(:) = 0.0_dp
               econf(0:maxll) = 0.5_dp*maxocc*REAL(elec_conf(0:maxll), dp)
            ELSE
               CYCLE
            END IF

            ! MOPAC TYEP GUESS
            IF (dft_control%qs_control%dftb) THEN
               DO iatom = 1, natom
                  atom_a = atom_list(iatom)
                  isgfa = first_sgf(atom_a)
                  DO la = 0, maxll
                     SELECT CASE (la)
                     CASE (0)
                        pdiag(isgfa) = econf(0)
                     CASE (1)
                        pdiag(isgfa+1) = econf(1)/3._dp
                        pdiag(isgfa+2) = econf(1)/3._dp
                        pdiag(isgfa+3) = econf(1)/3._dp
                     CASE (2)
                        pdiag(isgfa+4) = econf(2)/5._dp
                        pdiag(isgfa+5) = econf(2)/5._dp
                        pdiag(isgfa+6) = econf(2)/5._dp
                        pdiag(isgfa+7) = econf(2)/5._dp
                        pdiag(isgfa+8) = econf(2)/5._dp
                     CASE (3)
                        pdiag(isgfa+9) = econf(3)/7._dp
                        pdiag(isgfa+10) = econf(3)/7._dp
                        pdiag(isgfa+11) = econf(3)/7._dp
                        pdiag(isgfa+12) = econf(3)/7._dp
                        pdiag(isgfa+13) = econf(3)/7._dp
                        pdiag(isgfa+14) = econf(3)/7._dp
                        pdiag(isgfa+15) = econf(3)/7._dp
                     CASE DEFAULT
                        CPABORT("")
                     END SELECT
                  END DO
               END DO
            ELSEIF (dft_control%qs_control%semi_empirical) THEN
               yy = REAL(dft_control%charge, KIND=dp)/REAL(nao, KIND=dp)
               DO iatom = 1, natom
                  atom_a = atom_list(iatom)
                  isgfa = first_sgf(atom_a)
                  SELECT CASE (nsgf)
                  CASE (1) ! s-basis
                     pdiag(isgfa) = (zeff-yy)*0.5_dp*maxocc
                  CASE (4) ! sp-basis
                     IF (z == 1) THEN
                        ! special case: hydrogen with sp basis
                        pdiag(isgfa) = (zeff-yy)*0.5_dp*maxocc
                        pdiag(isgfa+1) = 0._dp
                        pdiag(isgfa+2) = 0._dp
                        pdiag(isgfa+3) = 0._dp
                     ELSE
                        pdiag(isgfa) = (zeff*0.25_dp-yy)*0.5_dp*maxocc
                        pdiag(isgfa+1) = (zeff*0.25_dp-yy)*0.5_dp*maxocc
                        pdiag(isgfa+2) = (zeff*0.25_dp-yy)*0.5_dp*maxocc
                        pdiag(isgfa+3) = (zeff*0.25_dp-yy)*0.5_dp*maxocc
                     END IF
                  CASE (9) ! spd-basis
                     IF (z < 21 .OR. z > 30 .AND. z < 39 .OR. z > 48 .AND. z < 57) THEN
                        !   Main Group Element:  The "d" shell is formally empty.
                        pdiag(isgfa) = (zeff*0.25_dp-yy)*0.5_dp*maxocc
                        pdiag(isgfa+1) = (zeff*0.25_dp-yy)*0.5_dp*maxocc
                        pdiag(isgfa+2) = (zeff*0.25_dp-yy)*0.5_dp*maxocc
                        pdiag(isgfa+3) = (zeff*0.25_dp-yy)*0.5_dp*maxocc
                        pdiag(isgfa+4) = (-yy)*0.5_dp*maxocc
                        pdiag(isgfa+5) = (-yy)*0.5_dp*maxocc
                        pdiag(isgfa+6) = (-yy)*0.5_dp*maxocc
                        pdiag(isgfa+7) = (-yy)*0.5_dp*maxocc
                        pdiag(isgfa+8) = (-yy)*0.5_dp*maxocc
                     ELSE IF (z < 99) THEN
                        my_sum = zeff-9.0_dp*yy
                        !   First, put 2 electrons in the 's' shell
                        pdiag(isgfa) = (MAX(0.0_dp, MIN(my_sum, 2.0_dp)))*0.5_dp*maxocc
                        my_sum = my_sum-2.0_dp
                        IF (my_sum > 0.0_dp) THEN
                           !   Now put as many electrons as possible into the 'd' shell
                           pdiag(isgfa+4) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           pdiag(isgfa+5) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           pdiag(isgfa+6) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           pdiag(isgfa+7) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           pdiag(isgfa+8) = (MAX(0.0_dp, MIN(my_sum*0.2_dp, 2.0_dp)))*0.5_dp*maxocc
                           my_sum = MAX(0.0_dp, my_sum-10.0_dp)
                           !   Put the remaining electrons in the 'p' shell
                           pdiag(isgfa+1) = (my_sum/3.0_dp)*0.5_dp*maxocc
                           pdiag(isgfa+2) = (my_sum/3.0_dp)*0.5_dp*maxocc
                           pdiag(isgfa+3) = (my_sum/3.0_dp)*0.5_dp*maxocc
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
                     nelec = maxocc*REAL(2*la+1, dp)
                     IF (econf(la) > 0.0_dp) THEN
                        IF (econf(la) >= nelec) THEN
                           paa = maxocc
                           econf(la) = econf(la)-nelec
                        ELSE
                           paa = maxocc*econf(la)/nelec
                           econf(la) = 0.0_dp
                           ncount = ncount+NINT(nelec/maxocc)
                        END IF
                        DO isgfa = first_sgfa(ishell, iset), last_sgfa(ishell, iset)
                           DO iatom = 1, natom
                              atom_a = atom_list(iatom)
                              isgf = first_sgf(atom_a)+isgfa-1
                              pdiag(isgf) = paa
                              IF (paa == maxocc) THEN
                                 trps1 = trps1+paa*sdiag(isgf)
                              ELSE
                                 trps2 = trps2+paa*sdiag(isgf)
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
         ELSE
            scale = (REAL(nelectron_spin, dp)-trps1)/trps2
            DO isgf = 1, nao
               IF (pdiag(isgf) < maxocc) pdiag(isgf) = scale*pdiag(isgf)
            END DO
         END IF
      END IF

      CALL dbcsr_set_diag(matrix_p, pdiag)

      DEALLOCATE (econf)

      DEALLOCATE (first_sgf)

      DEALLOCATE (pdiag)

      DEALLOCATE (sdiag)

      CALL timestop(handle)

   END SUBROUTINE calculate_mopac_dm

END MODULE qs_initial_guess
