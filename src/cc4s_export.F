!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Molecular and Gamma-point periodic CC4S export from restricted RI-HF.
!>        RI_FLAVOR MO stores a weighted AO-pair triangle. The full vertex is
!>        C^T (T + T^T) C, dressed with the symmetric Coulomb inverse square root.
!>        The SCF tensors and orbitals are read only. The dense MO vertex is assembled
!>        and written on the source rank.
! **************************************************************************************************
MODULE cc4s_export
   USE cell_types,                      ONLY: cell_type,&
                                              plane_distance,&
                                              use_perd_none,&
                                              use_perd_xyz
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_diag,                      ONLY: cp_fm_geeig_canon
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              cp_printkey_is_on
   USE dbt_api,                         ONLY: &
        dbt_get_block, dbt_get_info, dbt_iterator_blocks_left, dbt_iterator_next_block, &
        dbt_iterator_start, dbt_iterator_stop, dbt_iterator_type, dbt_nblks_total, dbt_type
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_truncated,&
                                              sic_none,&
                                              xc_none
   USE input_cp2k_hfx,                  ONLY: ri_mo
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp,&
                                              dp_size,&
                                              int_8
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_HOST,&
                                              local_gemm_ctxt_type
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              has_uniform_occupation,&
                                              mo_set_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: write_cc4s
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cc4s_export'

CONTAINS

! **************************************************************************************************
!> \brief Check the requested export, canonicalize private Fock/overlap copies and write tensors.
!> \param qs_env QS environment after SCF
!> \param converged Whether the SCF actually converged (including the outer SCF)
! **************************************************************************************************
   SUBROUTINE write_cc4s(qs_env, converged)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: converged

      CHARACTER(len=*), PARAMETER                        :: routineN = 'write_cc4s'
      REAL(dp), PARAMETER :: hf_fraction_tolerance = 1.0e-14_dp

      CHARACTER(len=8)                                   :: method
      CHARACTER(len=default_path_length)                 :: filename
      INTEGER :: frozen, functional, handle, homo, imo, inu, jmo, max_iter, method_id, nao, nmo, &
         nri, nstored, output_unit, poisson_periodic, stage_handle
      INTEGER, DIMENSION(2)                              :: dims2
      INTEGER, DIMENSION(3)                              :: dims3
      LOGICAL                                            :: do_kpoints, enabled, periodic
      REAL(dp) :: amp_conv, cutoff_radius, energy_conv, export_start, fermi_energy, max_cutoff, &
         max_memory_mb, memory_mb, write_start
      REAL(dp), ALLOCATABLE                              :: coeff(:, :), dressed(:, :), &
                                                            eigenvalues(:), gamma(:, :, :), &
                                                            half(:, :), metric(:, :), pair(:)
      REAL(dp), ALLOCATABLE, TARGET                      :: eri(:, :, :)
      REAL(dp), DIMENSION(2)                             :: elapsed
      REAL(dp), POINTER                                  :: eri_matrix(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_coeff, fm_fock, fm_overlap, fm_work
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(local_gemm_ctxt_type)                         :: gemm_context
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: input, print_section

      logger => cp_get_default_logger()
      CALL get_qs_env(qs_env, input=input)
      print_section => section_vals_get_subs_vals(input, "DFT%PRINT%CC4S")
      enabled = cp_printkey_is_on(logger%iter_info, print_section)
      IF (.NOT. enabled) RETURN
      IF (.NOT. qs_env%single_point_run) THEN
         CALL cp_abort(__LOCATION__, "CC4S: export requires a single-point RUN_TYPE ENERGY or ENERGY_FORCE")
      END IF
      IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, print_section), cp_p_file)) RETURN
      CALL section_vals_val_get(print_section, "FILENAME", c_val=filename)
      IF (filename(1:1) == '=' .OR. filename == '__STD_OUT__') THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: FILENAME must be a file prefix, not an exact name (=name) or __STD_OUT__")
      END IF
      CALL timeset(routineN, handle)
      export_start = m_walltime()
      elapsed = 0.0_dp
      output_unit = cp_logger_get_default_io_unit(logger)

      CALL section_vals_val_get(print_section, "METHOD", i_val=method_id)
      CALL section_vals_val_get(print_section, "N_FROZEN_CORE", i_val=frozen)
      CALL section_vals_val_get(print_section, "MAX_ITER", i_val=max_iter)
      CALL section_vals_val_get(print_section, "ENERGY_CONVERGENCE", r_val=energy_conv)
      CALL section_vals_val_get(print_section, "AMPLITUDES_CONVERGENCE", r_val=amp_conv)
      CALL section_vals_val_get(print_section, "MAX_MEMORY", r_val=max_memory_mb)
      SELECT CASE (method_id)
      CASE (1)
         method = 'MP2'
      CASE (2)
         method = 'CCSD'
      CASE (3)
         method = 'CCSD_T'
      CASE DEFAULT
         CPABORT("CC4S: unsupported method")
      END SELECT
      IF (max_iter < 1 .OR. energy_conv <= 0.0_dp .OR. amp_conv <= 0.0_dp .OR. max_memory_mb <= 0.0_dp) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: iteration, convergence and memory limits must be positive")
      END IF
      CALL get_qs_env(qs_env, x_data=x_data, mos=mos, para_env=para_env, dft_control=dft_control, &
                      cell=cell, do_kpoints=do_kpoints, scf_control=scf_control, &
                      matrix_ks=matrix_ks, matrix_s=matrix_s)
      IF (dft_control%qs_control%cdft) THEN
         CPABORT("CC4S: constrained DFT is unsupported")
      END IF
      IF (.NOT. converged) THEN
         CPABORT("CC4S: export requires a converged SCF")
      END IF
      CALL section_vals_val_get(input, "DFT%POISSON%PERIODIC", i_val=poisson_periodic)
      IF (do_kpoints) THEN
         CPABORT("CC4S: k-point sampling is unsupported")
      END IF
      periodic = ALL(cell%perd == 1) .AND. poisson_periodic == use_perd_xyz
      IF (.NOT. periodic .AND. &
          (ANY(cell%perd /= 0) .OR. poisson_periodic /= use_perd_none)) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: CELL and POISSON must both use PERIODIC NONE or both PERIODIC XYZ")
      END IF
      IF (periodic .AND. dft_control%charge /= 0) THEN
         CPABORT("CC4S: periodic cells must be neutral")
      END IF
      IF (dft_control%nspins /= 1 .OR. dft_control%uks .OR. dft_control%roks .OR. dft_control%smear) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: a restricted closed-shell reference without smearing is required")
      END IF
      IF (dft_control%do_admm .OR. dft_control%dft_plus_u .OR. dft_control%sic_method_id /= sic_none) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: ADMM, DFT+U and self-interaction corrections are unsupported")
      END IF
      CALL section_vals_val_get(input, "DFT%XC%XC_FUNCTIONAL%_SECTION_PARAMETERS_", i_val=functional)
      IF (functional /= xc_none) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: canonical HF orbitals required; use XC_FUNCTIONAL NONE")
      END IF
      IF (.NOT. ASSOCIATED(x_data)) THEN
         CPABORT("CC4S: enable HF with RI_FLAVOR MO")
      END IF
      IF (SIZE(x_data, 1) /= 1) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: exactly one HF section is allowed")
      END IF
      IF (ABS(x_data(1, 1)%general_parameter%fraction - 1.0_dp) > hf_fraction_tolerance) THEN
         CPABORT("CC4S: HF FRACTION must be 1")
      END IF
      IF (.NOT. ASSOCIATED(x_data(1, 1)%ri_data)) THEN
         CPABORT("CC4S: enable HF/RI with RI_FLAVOR MO")
      END IF
      IF (x_data(1, 1)%ri_data%flavor /= ri_mo) THEN
         CPABORT("CC4S: HF/RI requires RI_FLAVOR MO")
      END IF
      cutoff_radius = 0.0_dp
      IF (periodic) THEN
         IF (.NOT. x_data(1, 1)%ri_data%same_op .OR. &
             x_data(1, 1)%ri_data%hfx_pot%potential_type /= do_potential_truncated) THEN
            CALL cp_abort(__LOCATION__, &
                          "CC4S: periodic export requires TRUNCATED HF interaction and the same RI metric; use RI_METRIC HFX")
         END IF
         cutoff_radius = x_data(1, 1)%ri_data%hfx_pot%cutoff_radius
         max_cutoff = 0.5_dp*MIN(plane_distance(1, 0, 0, cell), &
                                 plane_distance(0, 1, 0, cell), plane_distance(0, 0, 1, cell))
         IF (cutoff_radius <= 0.0_dp .OR. cutoff_radius >= max_cutoff) THEN
            CALL cp_abort(__LOCATION__, &
                          "CC4S: positive HF CUTOFF_RADIUS must be smaller than half the minimal cell plane spacing")
         END IF
      ELSE
         IF (.NOT. x_data(1, 1)%ri_data%same_op .OR. &
             x_data(1, 1)%ri_data%hfx_pot%potential_type /= do_potential_coulomb) THEN
            CALL cp_abort(__LOCATION__, &
                          "CC4S: molecular HF interaction and RI metric must both be COULOMB")
         END IF
      END IF

      CALL get_mo_set(mos(1), mo_coeff=mo_coeff, homo=homo, nmo=nstored, nao=nao)
      IF (homo < 1 .OR. homo > nstored) THEN
         CPABORT("CC4S: no occupied orbitals")
      END IF
      IF (.NOT. has_uniform_occupation(mos(1), first_mo=1, last_mo=homo, occupation=2.0_dp) .OR. &
          .NOT. has_uniform_occupation(mos(1), first_mo=homo + 1, last_mo=nstored, occupation=0.0_dp)) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: occupations must be 2 below the HOMO and 0 above")
      END IF
      IF (frozen < 0 .OR. frozen >= homo) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: N_FROZEN_CORE must leave at least one occupied orbital")
      END IF
      CALL dbt_get_info(x_data(1, 1)%ri_data%t_3c_int_ctr_1(1, 1), nfull_total=dims3)
      nri = dims3(2)
      CPASSERT(dims3(1) == nao .AND. dims3(3) == nao .AND. nri > 0)
      IF (REAL(nri, dp)*REAL(nao, dp)**2 > REAL(HUGE(1), dp)) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: dense tensor exceeds the supported integer indexing range")
      END IF
      ! Conservative dense workspace estimate; excludes CP2K's existing SCF storage.
      memory_mb = 8.0_dp*(4.0_dp*REAL(nri, dp)*REAL(nao, dp)**2 + &
                          REAL(nri, dp)**2 + 6.0_dp*REAL(nao, dp)**2)/1.0e6_dp
      IF (output_unit > 0) WRITE (output_unit, '(/,T3,A,T69,F12.3)') &
         "CC4S| Estimated memory required for export per MPI rank [MB]:", memory_mb
      IF (memory_mb > max_memory_mb) THEN
         CPABORT("CC4S: dense export exceeds MAX_MEMORY [MB]")
      END IF

      ! Always diagonalize private copies, including for OT or a complete stored MO set.
      ! Drop overlap null modes using the same cutoff as SCF; do not export placeholder energies.
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, &
                               template_fmstruct=mo_coeff%matrix_struct)
      CALL cp_fm_create(fm_fock, fm_struct)
      CALL cp_fm_create(fm_overlap, fm_struct)
      CALL cp_fm_create(fm_coeff, fm_struct)
      CALL cp_fm_create(fm_work, fm_struct)
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(matrix_ks(1)%matrix, fm_fock)
      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, fm_overlap)
      ALLOCATE (eigenvalues(nao))
      CALL cp_fm_geeig_canon(fm_fock, fm_overlap, fm_coeff, eigenvalues, fm_work, &
                             scf_control%eps_eigval, nmo_retained=nmo)
      IF (nmo <= homo) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: no virtual orbitals after removing overlap null modes")
      END IF
      IF (eigenvalues(homo + 1) <= eigenvalues(homo)) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: a positive HOMO-LUMO gap is required")
      END IF
      fermi_energy = 0.5_dp*(eigenvalues(homo) + eigenvalues(homo + 1))
      nmo = nmo - frozen
      ALLOCATE (coeff(nao, nmo))
      CALL cp_fm_get_submatrix(fm_coeff, coeff, start_col=frozen + 1)
      CALL cp_fm_release(fm_fock)
      CALL cp_fm_release(fm_overlap)
      CALL cp_fm_release(fm_coeff)
      CALL cp_fm_release(fm_work)

      ALLOCATE (eri(nri, nao, nao), metric(nri, nri))
      CALL timeset("cc4s_gather", stage_handle)
      CALL gather_3c_pao(x_data(1, 1)%ri_data%t_3c_int_ctr_1(1, 1), para_env, nri, nao, eri)
      CALL dbt_get_info(x_data(1, 1)%ri_data%t_2c_int(1, 1), nfull_total=dims2)
      CPASSERT(ALL(dims2 == nri))
      CALL gather_2c_dense(x_data(1, 1)%ri_data%t_2c_int(1, 1), para_env, nri, metric)
      CALL timestop(stage_handle)

      IF (para_env%is_source()) THEN
         ALLOCATE (dressed(nri, nao*nao), half(nri, nmo), gamma(nri, nmo, nmo), pair(nri))
         CALL gemm_context%create(LOCAL_GEMM_PU_HOST)
         ! Identical HF and RI operators give a symmetric inverse square root.
         ! Remap the contiguous AO-pair dimensions without copying the dense tensor.
         CALL timeset("cc4s_metric", stage_handle)
         eri_matrix(1:nri, 1:nao*nao) => eri
         CALL gemm_context%gemm('T', 'N', nri, nao*nao, nri, 1.0_dp, metric, nri, eri_matrix, nri, &
                                0.0_dp, dressed, nri)
         NULLIFY (eri_matrix)
         CALL timestop(stage_handle)
         CALL timeset("cc4s_transform", stage_handle)
         gamma(:, :, :) = 0.0_dp
         DO inu = 1, nao
            CALL gemm_context%gemm('N', 'N', nri, nmo, nao, 1.0_dp, dressed(:, (inu - 1)*nao + 1:inu*nao), nri, &
                                   coeff, nao, 0.0_dp, half, nri)
            DO imo = 1, nmo
               gamma(:, :, imo) = gamma(:, :, imo) + coeff(inu, imo)*half(:, :)
            END DO
         END DO
         CALL gemm_context%destroy()
         ! Restore the weighted AO-pair triangle AFTER the linear MO transform.
         ! Sum, not average: diagonal atom blocks in the stored tensor have weight 1/2.
         DO imo = 1, nmo
            DO jmo = imo, nmo
               pair(:) = gamma(:, imo, jmo) + gamma(:, jmo, imo)
               gamma(:, imo, jmo) = pair
               gamma(:, jmo, imo) = pair
            END DO
         END DO
         CALL timestop(stage_handle)
         CALL timeset("cc4s_write_files", stage_handle)
         write_start = m_walltime()
         CALL write_eigen_energies_files(logger, print_section, nmo, homo - frozen, &
                                         eigenvalues(frozen + 1:frozen + nmo), fermi_energy)
         ! A finite-range periodic interaction has no Ewald/Madelung correction.
         ! Preserve the SCF eigenenergies and its live truncated-interaction vertex.
         CALL write_coulomb_vertex_files(logger, print_section, nri, nmo, gamma, periodic, cutoff_radius, cell%hmat)
         CALL write_cc4s_in(logger, print_section, method, max_iter, energy_conv, amp_conv)
         elapsed(2) = m_walltime() - write_start
         CALL timestop(stage_handle)
         IF (output_unit > 0) THEN
            IF (method == 'CCSD_T') THEN
               WRITE (output_unit, '(T3,A,T73,A8)') "CC4S| Requested method:", 'CCSD(T)'
            ELSE
               WRITE (output_unit, '(T3,A,T73,A8)') "CC4S| Requested method:", TRIM(method)
            END IF
            IF (method /= 'MP2') THEN
               WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Maximum CCSD iterations:", max_iter
               WRITE (output_unit, '(T3,A,T73,ES8.1)') "CC4S| CCSD energy convergence [Ha]:", energy_conv
               WRITE (output_unit, '(T3,A,T73,ES8.1)') "CC4S| CCSD amplitudes convergence:", amp_conv
            END IF
            WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Number of auxiliary basis functions:", nri
            WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Number of exported orbitals:", nmo
            WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Number of exported occupied orbitals:", homo - frozen
            WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Number of exported virtual orbitals:", nmo - homo + frozen
            WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Number of frozen occupied orbitals:", frozen
            WRITE (output_unit, '(T3,A,T56,ES25.16)') "CC4S| Vertex norm squared [Ha]:", SUM(gamma**2)
            WRITE (output_unit, '(T3,A,T56,ES25.16)') &
               "CC4S| Sum of exported orbital energies [Ha]:", SUM(eigenvalues(frozen + 1:frozen + nmo))
            IF (periodic) THEN
               WRITE (output_unit, '(T3,A)') "CC4S| Gamma-point TRUNCATED interaction"
               WRITE (output_unit, '(T3,A,T69,F12.6)') "CC4S| Coulomb cutoff radius [bohr]:", cutoff_radius
            END IF
            filename = cp_print_key_generate_filename(logger, print_section, extension='.in', my_local=.FALSE.)
            WRITE (output_unit, '(T3,A,T41,A)') "CC4S| Generated input file:", &
               ADJUSTR(filename(:MAX(40, LEN_TRIM(filename))))
         END IF
         DEALLOCATE (dressed, half, gamma, pair)
      END IF
      DEALLOCATE (eri, metric, coeff, eigenvalues)
      CALL para_env%sync()
      elapsed(1) = m_walltime() - export_start
      ! Report the longest rank-local elapsed times; only the source rank writes files.
      CALL para_env%max(elapsed)
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T3,A,T69,F12.3)') "CC4S| Total export time [s]:", elapsed(1)
         WRITE (output_unit, '(T3,A,T69,F12.3)') "CC4S| File writing time [s]:", elapsed(2)
      END IF
      CALL timestop(handle)
   END SUBROUTINE write_cc4s

! **************************************************************************************************
!> \brief  Gather the (AO, RI, AO) = (mu, Q, nu) block tensor into a dense array
!>         with index order (Q, mu, nu) so Q varies fastest in Fortran storage.
!>         Contributions are reduced to the source rank.
!> \param tensor Distributed RI tensor
!> \param para_env MPI parallel environment
!> \param nRI Number of auxiliary functions
!> \param nAO Number of AO basis functions
!> \param eri_pao Dense auxiliary-AO-AO tensor on the source rank
! **************************************************************************************************
   SUBROUTINE gather_3c_pao(tensor, para_env, nRI, nAO, eri_pao)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: nRI, nAO
      REAL(KIND=dp), DIMENSION(nRI, nAO, nAO), &
         INTENT(OUT)                                     :: eri_pao

      INTEGER                                            :: imu, inu, nb1, nb2, nb3
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: off1, off2, off3
      INTEGER, DIMENSION(3)                              :: ind3
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: blk3
      TYPE(dbt_iterator_type)                            :: iter

      eri_pao(:, :, :) = 0.0_dp

      ! Allocate offset arrays sized by number of blocks in each dimension
      nb1 = dbt_nblks_total(tensor, 1)
      nb2 = dbt_nblks_total(tensor, 2)
      nb3 = dbt_nblks_total(tensor, 3)
      ALLOCATE (off1(nb1), off2(nb2), off3(nb3))
      CALL dbt_get_info(tensor, blk_offset_1=off1, blk_offset_2=off2, blk_offset_3=off3)

      ! DBT iterator requires an OpenMP parallel region; each thread processes its own shard
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor, off1, off2, off3, eri_pao) &
!$OMP PRIVATE(iter, ind3, blk3, found, imu, inu)
      CALL dbt_iterator_start(iter, tensor)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind3)
         CALL dbt_get_block(tensor, ind3, blk3, found)
         IF (.NOT. found) CYCLE

         ! Source: blk3(mu_loc, Q_loc, nu_loc), tensor index order = (AO, RI, AO)
         ! Target: eri_pao(Q, mu, nu)  — transpose first two indices
         DO inu = 1, SIZE(blk3, 3)
            DO imu = 1, SIZE(blk3, 1)
               eri_pao(off2(ind3(2)):off2(ind3(2)) + SIZE(blk3, 2) - 1, &
                       off1(ind3(1)) + imu - 1, &
                       off3(ind3(3)) + inu - 1) = blk3(imu, :, inu)
            END DO
         END DO

         DEALLOCATE (blk3)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      DEALLOCATE (off1, off2, off3)

      ! Reduce: collect contributions on the source rank
      CALL para_env%sum(eri_pao, para_env%source)

   END SUBROUTINE gather_3c_pao

! **************************************************************************************************
!> \brief  Gather a (RI, RI) block tensor to a dense matrix.  Reduces to the source rank.
!> \param tensor Distributed RI tensor
!> \param para_env MPI parallel environment
!> \param nRI Number of auxiliary functions
!> \param dense Dense auxiliary matrix on the source rank
! **************************************************************************************************
   SUBROUTINE gather_2c_dense(tensor, para_env, nRI, dense)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: nRI
      REAL(KIND=dp), DIMENSION(nRI, nRI), INTENT(OUT)    :: dense

      INTEGER                                            :: nb1, nb2
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: off1, off2
      INTEGER, DIMENSION(2)                              :: ind2
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: blk2
      TYPE(dbt_iterator_type)                            :: iter

      dense(:, :) = 0.0_dp

      nb1 = dbt_nblks_total(tensor, 1)
      nb2 = dbt_nblks_total(tensor, 2)
      ALLOCATE (off1(nb1), off2(nb2))
      CALL dbt_get_info(tensor, blk_offset_1=off1, blk_offset_2=off2)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(tensor, off1, off2, dense) &
!$OMP PRIVATE(iter, ind2, blk2, found)
      CALL dbt_iterator_start(iter, tensor)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind2)
         CALL dbt_get_block(tensor, ind2, blk2, found)
         IF (.NOT. found) CYCLE
         dense(off1(ind2(1)):off1(ind2(1)) + SIZE(blk2, 1) - 1, &
               off2(ind2(2)):off2(ind2(2)) + SIZE(blk2, 2) - 1) = blk2
         DEALLOCATE (blk2)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      DEALLOCATE (off1, off2)

      CALL para_env%sum(dense, para_env%source)

   END SUBROUTINE gather_2c_dense

! **************************************************************************************************
!> \brief  Write EigenEnergies.yaml and EigenEnergies.elements (CC4S format, Hartree units).
!> \param logger Logger for output selection and naming
!> \param print_section CC4S print key
!> \param nMO Number of active canonical orbitals
!> \param homo Number of active occupied orbitals
!> \param eigenvalues Active orbital energies in Hartree
!> \param fermi_energy Energy inside the HOMO-LUMO gap
! **************************************************************************************************
   SUBROUTINE write_eigen_energies_files(logger, print_section, nMO, homo, eigenvalues, fermi_energy)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_section
      INTEGER, INTENT(IN)                                :: nMO, homo
      REAL(KIND=dp), INTENT(IN)                          :: eigenvalues(nMO), fermi_energy

      INTEGER                                            :: i, iunit

      iunit = cp_print_key_unit_nr(logger, print_section, extension='.yaml', middle_name='EigenEnergies', &
                                   file_status='REPLACE', file_position='REWIND')
      WRITE (iunit, '(A)') 'version: 100'
      WRITE (iunit, '(A)') 'type: Tensor'
      WRITE (iunit, '(A)') 'scalarType: Real64'
      WRITE (iunit, '(A)') 'dimensions:'
      WRITE (iunit, '(A,I0)') '- length: ', nMO
      WRITE (iunit, '(A)') '  type: State'
      WRITE (iunit, '(A)') 'elements:'
      WRITE (iunit, '(A)') '  type: TextFile'
      WRITE (iunit, '(A)') 'unit: 1.0'
      WRITE (iunit, '(A)') 'metaData:'
      WRITE (iunit, '(A,ES25.16)') '  fermiEnergy: ', fermi_energy
      WRITE (iunit, '(A,I0)') '  homoIndex: ', homo
      WRITE (iunit, '(A)') '  energies:'
      DO i = 1, nMO
         WRITE (iunit, '(A,ES25.16)') '  - ', eigenvalues(i)
      END DO
      CALL cp_print_key_finished_output(iunit, logger, print_section)

      iunit = cp_print_key_unit_nr(logger, print_section, extension='.elements', middle_name='EigenEnergies', &
                                   file_status='REPLACE', file_position='REWIND')
      DO i = 1, nMO
         WRITE (iunit, '(ES25.16)') eigenvalues(i)
      END DO
      CALL cp_print_key_finished_output(iunit, logger, print_section)

   END SUBROUTINE write_eigen_energies_files

! **************************************************************************************************
!> \brief Stream the dense vertex in blocks of at most 8 MiB without copying it.
!>        Preserve Fortran element order with the auxiliary index varying fastest.
!> \param logger Logger for output selection and naming
!> \param print_section CC4S print key
!> \param nRI Number of auxiliary functions
!> \param nMO Number of active canonical orbitals
!> \param gamma Coulomb vertex in auxiliary-orbital-orbital order
!> \param periodic Whether the cell is three-dimensionally periodic at Gamma
!> \param cutoff_radius Truncated Coulomb cutoff in bohr (zero for molecular Coulomb)
!> \param hmat Cell vectors as columns, in bohr
! **************************************************************************************************
   SUBROUTINE write_coulomb_vertex_files(logger, print_section, nRI, nMO, gamma, periodic, cutoff_radius, hmat)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_section
      INTEGER, INTENT(IN)                                :: nRI, nMO
      REAL(KIND=dp), INTENT(IN), TARGET                  :: gamma(nRI, nMO, nMO)
      LOGICAL, INTENT(IN)                                :: periodic
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_radius, hmat(3, 3)

      INTEGER(KIND=int_8), PARAMETER :: max_write_elements = 1048576_int_8

      INTEGER                                            :: iunit
      INTEGER(KIND=int_8)                                :: count, nvalues, offset
      REAL(KIND=dp), CONTIGUOUS, POINTER                 :: values(:)

      CALL open_coulomb_vertex_file(logger, print_section, nRI, nMO, periodic, cutoff_radius, hmat, iunit)
      nvalues = SIZE(gamma, KIND=int_8)
      values(1:nvalues) => gamma
      offset = 0_int_8
      DO WHILE (offset < nvalues)
         count = MIN(max_write_elements, nvalues - offset)
         CALL write_coulomb_vertex_block(iunit, offset, values(offset + 1:offset + count))
         offset = offset + count
      END DO
      NULLIFY (values)
      CALL cp_print_key_finished_output(iunit, logger, print_section)
   END SUBROUTINE write_coulomb_vertex_files

! **************************************************************************************************
!> \brief Write the Coulomb vertex header and open its raw binary elements file.
!> \param logger Logger for output selection and naming
!> \param print_section CC4S print key
!> \param nRI Number of auxiliary functions
!> \param nMO Number of active canonical orbitals
!> \param periodic Whether the cell is three-dimensionally periodic at Gamma
!> \param cutoff_radius Truncated Coulomb cutoff in bohr (zero for molecular Coulomb)
!> \param hmat Cell vectors as columns, in bohr
!> \param iunit Open stream unit for vertex elements
! **************************************************************************************************
   SUBROUTINE open_coulomb_vertex_file(logger, print_section, nRI, nMO, periodic, cutoff_radius, hmat, iunit)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_section
      INTEGER, INTENT(IN)                                :: nRI, nMO
      LOGICAL, INTENT(IN)                                :: periodic
      REAL(KIND=dp), INTENT(IN)                          :: cutoff_radius, hmat(3, 3)
      INTEGER, INTENT(OUT)                               :: iunit

      INTEGER                                            :: i

      iunit = cp_print_key_unit_nr(logger, print_section, extension='.yaml', middle_name='CoulombVertex', &
                                   file_status='REPLACE', file_position='REWIND')
      WRITE (iunit, '(A)') 'version: 100'
      WRITE (iunit, '(A)') 'type: Tensor'
      WRITE (iunit, '(A)') 'scalarType: Real64'
      WRITE (iunit, '(A)') 'dimensions:'
      WRITE (iunit, '(A,I0)') '- length: ', nRI
      WRITE (iunit, '(A)') '  type: AuxiliaryField'
      WRITE (iunit, '(A,I0)') '- length: ', nMO
      WRITE (iunit, '(A)') '  type: State'
      WRITE (iunit, '(A,I0)') '- length: ', nMO
      WRITE (iunit, '(A)') '  type: State'
      WRITE (iunit, '(A)') 'elements:'
      WRITE (iunit, '(A)') '  type: IeeeBinaryFile'
      WRITE (iunit, '(A)') 'unit: 1.0'
      WRITE (iunit, '(A)') 'metaData:'
      IF (periodic) THEN
         WRITE (iunit, '(A)') '  boundaryConditions: "periodic-3d-gamma"'
         WRITE (iunit, '(A)') '  interactionPotential: "TRUNCATED"'
         WRITE (iunit, '(A,ES25.16)') '  cutoffRadiusBohr: ', cutoff_radius
         WRITE (iunit, '(A)') '  madelungCorrection: 0.0'
         WRITE (iunit, '(A)') '  latticeVectorsBohr:'
         DO i = 1, 3
            WRITE (iunit, '(A,ES25.16,A,ES25.16,A,ES25.16,A)') &
               '  - [', hmat(1, i), ', ', hmat(2, i), ', ', hmat(3, i), ']'
         END DO
      ELSE
         WRITE (iunit, '(A)') '  boundaryConditions: "molecular"'
         WRITE (iunit, '(A)') '  interactionPotential: "COULOMB"'
      END IF
      CALL cp_print_key_finished_output(iunit, logger, print_section)

      iunit = cp_print_key_unit_nr(logger, print_section, extension='.elements', middle_name='CoulombVertex', &
                                   file_status='REPLACE', file_position='REWIND', &
                                   file_form='UNFORMATTED', file_access='STREAM')
   END SUBROUTINE open_coulomb_vertex_file

! **************************************************************************************************
!> \brief Write a contiguous vertex block at an explicit position in the raw elements file.
!> \param iunit Open stream unit for vertex elements
!> \param offset Zero-based element offset in auxiliary-orbital-orbital storage order
!> \param values Contiguous block of vertex elements
! **************************************************************************************************
   SUBROUTINE write_coulomb_vertex_block(iunit, offset, values)
      INTEGER, INTENT(IN)                                :: iunit
      INTEGER(KIND=int_8), INTENT(IN)                    :: offset
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: values(:)

      INTEGER                                            :: io_status
      INTEGER(KIND=int_8)                                :: position

      ! Use CP2K's Real64 byte size, independent of direct-access record-length units.
      CPASSERT(offset >= 0_int_8)
      CPASSERT(SIZE(values, KIND=int_8) > 0_int_8)
      CPASSERT(SIZE(values, KIND=int_8) <= HUGE(position)/INT(dp_size, int_8) - offset)
      position = offset*INT(dp_size, int_8) + 1_int_8
      WRITE (iunit, POS=position, IOSTAT=io_status) values
      IF (io_status /= 0) THEN
         CPABORT("CC4S: failed writing CoulombVertex.elements")
      END IF
   END SUBROUTINE write_coulomb_vertex_block

! **************************************************************************************************
!> \brief  Write the CC4S algorithm sequence and solver settings to cc4s.in.
!> \param logger Logger for output selection and naming
!> \param print_section CC4S print key
!> \param method Requested correlation method: MP2, CCSD or CCSD_T
!> \param max_iter Maximum coupled-cluster iterations
!> \param energy_conv CC4S energy convergence threshold
!> \param amp_conv CC4S amplitude convergence threshold
! **************************************************************************************************
   SUBROUTINE write_cc4s_in(logger, print_section, method, max_iter, energy_conv, amp_conv)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_section
      CHARACTER(LEN=*), INTENT(IN)                       :: method
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(dp), INTENT(IN)                               :: energy_conv, amp_conv

      CHARACTER(len=default_path_length)                 :: filename
      INTEGER                                            :: iunit

      iunit = cp_print_key_unit_nr(logger, print_section, extension='.in', &
                                   file_status='REPLACE', file_position='REWIND')

      ! --- common header: read tensors, slice, build integrals ---
      filename = cp_print_key_generate_filename(logger, print_section, 'CoulombVertex', '.yaml', .FALSE.)
      CALL write_tensor_read(iunit, filename, 'CoulombVertex')
      filename = cp_print_key_generate_filename(logger, print_section, 'EigenEnergies', '.yaml', .FALSE.)
      CALL write_tensor_read(iunit, filename, 'EigenEnergies')
      WRITE (iunit, '(A)') '- name: DefineHolesAndParticles'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    eigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '- name: SliceOperator'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '    operator: CoulombVertex'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    slicedOperator: CoulombVertex'
      WRITE (iunit, '(A)') '- name: VertexCoulombIntegrals'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    slicedCoulombVertex: CoulombVertex'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'

      ! --- MP2 (always included) ---
      WRITE (iunit, '(A)') '- name: SecondOrderPerturbationTheory'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    energy: Mp2Energy'

      ! --- CoupledCluster (CCSD, CCSD_T) ---
      SELECT CASE (TRIM(method))
      CASE ('MP2')
         ! The MP2 driver is complete above.
      CASE ('CCSD', 'CCSD_T')
         CALL write_coupled_cluster_step(iunit, 'Ccsd', 'Amplitudes', max_iter, energy_conv, amp_conv)
      CASE DEFAULT
         CPABORT("CC4S: unsupported method")
      END SELECT

      ! --- PerturbativeTriples (CCSD_T only) ---
      IF (TRIM(method) == 'CCSD_T') THEN
         WRITE (iunit, '(A)') '- name: PerturbativeTriples'
         WRITE (iunit, '(A)') '  in:'
         WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
         WRITE (iunit, '(A)') '    amplitudes: Amplitudes'
         WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'
         WRITE (iunit, '(A)') '  out: {}'
      END IF

      CALL cp_print_key_finished_output(iunit, logger, print_section)

   END SUBROUTINE write_cc4s_in

! **************************************************************************************************
!> \brief  Write a single CoupledCluster step to the cc4s.in file.
!> \param iunit Output file unit
!> \param cc_method CC4S coupled-cluster method
!> \param amp_key CC4S output key for amplitudes
!> \param max_iter Maximum coupled-cluster iterations
!> \param energy_conv CC4S energy convergence threshold
!> \param amp_conv CC4S amplitude convergence threshold
! **************************************************************************************************
   SUBROUTINE write_coupled_cluster_step(iunit, cc_method, amp_key, max_iter, energy_conv, amp_conv)
      INTEGER, INTENT(IN)                                :: iunit
      CHARACTER(LEN=*), INTENT(IN)                       :: cc_method, amp_key
      INTEGER, INTENT(IN)                                :: max_iter
      REAL(dp), INTENT(IN)                               :: energy_conv, amp_conv

      WRITE (iunit, '(A)') '- name: CoupledCluster'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    method: '//TRIM(cc_method)
      WRITE (iunit, '(A)') '    integralsSliceSize: 100'
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'
      WRITE (iunit, '(A)') '    slicedCoulombVertex: CoulombVertex'
      WRITE (iunit, '(A,I0)') '    maxIterations: ', max_iter
      WRITE (iunit, '(A,ES25.16)') '    energyConvergence: ', energy_conv
      WRITE (iunit, '(A,ES25.16)') '    amplitudesConvergence: ', amp_conv
      WRITE (iunit, '(A)') '    mixer:'
      WRITE (iunit, '(A)') '      type: DiisMixer'
      WRITE (iunit, '(A)') '      maxResidua: 5'
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    amplitudes: '//TRIM(amp_key)
      WRITE (iunit, '(A)') '    energy: CoupledClusterEnergy'

   END SUBROUTINE write_coupled_cluster_step

! **************************************************************************************************
!> \brief Write a tensor reference relative to the generated CC4S input's directory.
!> \param iunit CC4S input unit
!> \param filename Generated tensor header filename
!> \param destination CC4S tensor name
! **************************************************************************************************
   SUBROUTINE write_tensor_read(iunit, filename, destination)
      INTEGER, INTENT(IN)                                :: iunit
      CHARACTER(len=*), INTENT(IN)                       :: filename, destination

      INTEGER                                            :: i

      WRITE (iunit, '(A)') '- name: Read'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)', ADVANCE='NO') "    fileName: '"
      ! A quoted YAML scalar preserves spaces and punctuation; apostrophes are doubled.
      DO i = INDEX(filename, '/', BACK=.TRUE.) + 1, LEN_TRIM(filename)
         WRITE (iunit, '(A)', ADVANCE='NO') filename(i:i)
         IF (filename(i:i) == "'") WRITE (iunit, '(A)', ADVANCE='NO') "'"
      END DO
      WRITE (iunit, '(A)') "'"
      WRITE (iunit, '(A)') '  out:'
      WRITE (iunit, '(A)') '    destination: '//destination
   END SUBROUTINE write_tensor_read

END MODULE cc4s_export
