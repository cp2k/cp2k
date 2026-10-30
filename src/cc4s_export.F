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
!>        The SCF tensors and orbitals are read only. MO transformations are distributed
!>        and bounded output pieces are gathered to the writing rank.
! **************************************************************************************************
MODULE cc4s_export
   USE cell_types,                      ONLY: cell_type,&
                                              plane_distance,&
                                              use_perd_none,&
                                              use_perd_xyz
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_fm_to_dbcsr_row_template
   USE cp_fm_diag,                      ONLY: cp_fm_geeig_canon
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
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
        dbt_contract, dbt_copy, dbt_copy_matrix_to_tensor, dbt_create, dbt_destroy, dbt_get_block, &
        dbt_get_info, dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, &
        dbt_iterator_stop, dbt_iterator_type, dbt_nblks_total, dbt_pgrid_create, &
        dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE hfx_types,                       ONLY: hfx_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
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
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              has_uniform_occupation,&
                                              mo_set_type
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              split_block_sizes
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: write_cc4s
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cc4s_export'
   REAL(dp), PARAMETER :: elements_per_mib = 1024.0_dp**2/REAL(dp_size, dp)

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
      INTEGER :: frozen, functional, handle, homo, max_iter, method_id, nao, nmo, nri, nstored, &
         output_unit, poisson_periodic, stage_handle
      INTEGER(KIND=int_8)                                :: buffer_elements
      INTEGER, DIMENSION(2)                              :: dims2
      INTEGER, DIMENSION(3)                              :: dims3
      LOGICAL                                            :: do_kpoints, enabled, periodic
      REAL(dp) :: amp_conv, cutoff_radius, energy_conv, export_start, fermi_energy, max_cutoff, &
         max_memory_mb, memory_mb, output_buffer_mib, vertex_norm, write_start
      REAL(dp), ALLOCATABLE                              :: eigenvalues(:)
      REAL(dp), DIMENSION(5)                             :: elapsed
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_active, fm_coeff, fm_fock, &
                                                            fm_overlap, fm_work
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dbcsr_type)                                   :: coefficient_matrix
      TYPE(dbt_type)                                     :: coefficients
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
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
      CALL section_vals_val_get(print_section, "OUTPUT_BUFFER_SIZE", r_val=output_buffer_mib)
      IF (.NOT. ieee_is_finite(output_buffer_mib)) THEN
         CPABORT("CC4S: OUTPUT_BUFFER_SIZE must be finite")
      END IF
      IF (output_buffer_mib < 1.0_dp/elements_per_mib .OR. &
          output_buffer_mib > REAL(HUGE(1), dp)/elements_per_mib) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: OUTPUT_BUFFER_SIZE must hold at least one Real64 value and fit the MPI integer count")
      END IF
      buffer_elements = FLOOR(output_buffer_mib*elements_per_mib, KIND=int_8)
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
      IF (.NOT. ieee_is_finite(max_memory_mb)) CPABORT("CC4S: MAX_MEMORY must be finite")
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
      ! Canonicalization precedes batch selection; reject an insufficient budget before allocation.
      memory_mb = REAL(dp_size, dp)*6.0_dp*REAL(nao, dp)**2/1.0e6_dp
      IF (memory_mb > max_memory_mb) THEN
         CPABORT("CC4S: canonicalization exceeds MAX_MEMORY [MB]")
      END IF
      elapsed = 0.0_dp

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
      buffer_elements = MIN(buffer_elements, INT(nri, int_8)*INT(nmo, int_8)**2)
      IF (output_unit > 0) WRITE (output_unit, '(/,T3,A,T69,F12.6)') &
         "CC4S| Actual output buffer size per MPI rank [MiB]:", REAL(buffer_elements, dp)/elements_per_mib
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nmo, template_fmstruct=mo_coeff%matrix_struct)
      CALL cp_fm_create(fm_active, fm_struct, name="CC4S active orbitals")
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_to_fm(fm_coeff, fm_active, ncol=nmo, source_start=frozen + 1)
      CALL cp_fm_to_dbcsr_row_template(coefficient_matrix, fm_active, matrix_s(1)%matrix)
      CALL dbt_create(coefficient_matrix, coefficients, name="CC4S coefficients")
      CALL dbt_copy_matrix_to_tensor(coefficient_matrix, coefficients)
      CALL dbcsr_release(coefficient_matrix)
      CALL cp_fm_release(fm_active)
      CALL cp_fm_release(fm_fock)
      CALL cp_fm_release(fm_overlap)
      CALL cp_fm_release(fm_coeff)
      CALL cp_fm_release(fm_work)

      CALL dbt_get_info(x_data(1, 1)%ri_data%t_2c_int(1, 1), nfull_total=dims2)
      CPASSERT(ALL(dims2 == nri))
      ! Transform and export one orbital range at a time; never retain the complete vertex.
      CALL export_vertex_batches(x_data(1, 1)%ri_data%t_3c_int_ctr_1(1, 1), &
                                 x_data(1, 1)%ri_data%t_2c_int(1, 1), coefficients, para_env, &
                                 logger, print_section, nmo, max_memory_mb, buffer_elements, &
                                 periodic, cutoff_radius, cell%hmat, output_unit, vertex_norm, elapsed, &
                                 homo - frozen, eigenvalues(frozen + 1:frozen + nmo), fermi_energy)
      CALL dbt_destroy(coefficients)
      CALL timeset("cc4s_write_files", stage_handle)
      write_start = m_walltime()
      IF (para_env%is_source()) THEN
         CALL write_cc4s_in(logger, print_section, method, homo - frozen, max_iter, energy_conv, amp_conv)
      END IF
      elapsed(5) = elapsed(5) + m_walltime() - write_start
      CALL timestop(stage_handle)
      IF (para_env%is_source()) THEN
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
            WRITE (output_unit, '(T3,A,T56,ES25.16)') "CC4S| Vertex norm squared [Ha]:", vertex_norm
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
      END IF
      DEALLOCATE (eigenvalues)
      CALL para_env%sync()
      elapsed(1) = m_walltime() - export_start
      ! Report the longest rank-local elapsed times; only the source rank writes files.
      CALL para_env%max(elapsed)
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T3,A,T69,F12.3)') "CC4S| Total export time [s]:", elapsed(1)
         WRITE (output_unit, '(T3,A,T69,F12.3)') "CC4S| Orbital transformation time [s]:", elapsed(2)
         WRITE (output_unit, '(T3,A,T69,F12.3)') "CC4S| Metric contraction time [s]:", elapsed(3)
         WRITE (output_unit, '(T3,A,T69,F12.3)') "CC4S| Output gathering time [s]:", elapsed(4)
         WRITE (output_unit, '(T3,A,T69,F12.3)') "CC4S| File writing time [s]:", elapsed(5)
      END IF
      CALL timestop(handle)
   END SUBROUTINE write_cc4s

! **************************************************************************************************
!> \brief Transform and stream orbital batches, using private auxiliary slices of the SCF triangle.
!> \param integrals Live SCF tensor, ordered (AO, RI, AO)
!> \param metric Live SCF inverse square root of the RI metric
!> \param coefficients Distributed active canonical orbitals
!> \param para_env MPI parallel environment
!> \param logger Logger for output selection and naming
!> \param print_section CC4S print key
!> \param nmo Number of active orbitals
!> \param max_memory_mb Estimated additional workspace limit per rank, in decimal MB
!> \param buffer_elements Maximum number of Real64 values in each output piece
!> \param periodic Whether the cell is three-dimensionally periodic at Gamma
!> \param cutoff_radius Truncated Coulomb cutoff in bohr
!> \param hmat Cell vectors in bohr
!> \param output_unit Log output unit
!> \param vertex_norm Squared vertex norm on the writing rank
!> \param elapsed Accumulated wall times: total, transformation, metric, gathering, writing
!> \param homo Number of active occupied orbitals
!> \param eigenvalues Active canonical orbital energies
!> \param fermi_energy Energy inside the HOMO-LUMO gap
! **************************************************************************************************
   SUBROUTINE export_vertex_batches(integrals, metric, coefficients, para_env, logger, print_section, &
                                    nmo, max_memory_mb, buffer_elements, periodic, cutoff_radius, hmat, &
                                    output_unit, vertex_norm, elapsed, homo, eigenvalues, fermi_energy)
      TYPE(dbt_type), INTENT(INOUT)                      :: integrals, metric, coefficients
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_section
      INTEGER, INTENT(IN)                                :: nmo
      REAL(dp), INTENT(IN)                               :: max_memory_mb
      INTEGER(KIND=int_8), INTENT(IN)                    :: buffer_elements
      LOGICAL, INTENT(IN)                                :: periodic
      REAL(dp), INTENT(IN)                               :: cutoff_radius, hmat(3, 3)
      INTEGER, INTENT(IN)                                :: output_unit
      REAL(dp), INTENT(OUT)                              :: vertex_norm
      REAL(dp), INTENT(INOUT)                            :: elapsed(5)
      INTEGER, INTENT(IN)                                :: homo
      REAL(dp), INTENT(IN)                               :: eigenvalues(nmo), fermi_energy

      INTEGER                                            :: bounds2(2, 2), bounds3(2, 3), handle, &
                                                            im, ir, iunit, mb, mfirst, mlast, nao, &
                                                            nri, pdims(3), pdims2(2), rb, rfirst, &
                                                            rlast
      INTEGER, ALLOCATABLE                               :: ao(:), d1(:), d2(:), d3(:), mo(:), ri(:)
      REAL(dp)                                           :: memory_mb, start, trial_m, trial_r
      REAL(dp), ALLOCATABLE                              :: values(:)
      TYPE(dbt_pgrid_type)                               :: grid, grid2
      TYPE(dbt_type)                                     :: ao_triangle, coeff, coeff_batch, gamma, &
                                                            half, metric_copy, triangle

      ALLOCATE (ao(dbt_nblks_total(integrals, 1)), ri(dbt_nblks_total(integrals, 2)))
      CALL dbt_get_info(integrals, blk_size_1=ao, blk_size_2=ri)
      nao = SUM(ao)
      nri = SUM(ri)
      CALL split_block_sizes([nmo], mo, max_size=16)
      mb = SIZE(mo)
      rb = SIZE(ri)
      ! Bound the busiest rank without assuming balanced distribution. Include private operands,
      ! simultaneously live results, and conservative redistribution/communication copies.
      memory_mb = batch_workspace(nao, nri, nmo, mo, ri, mb, rb, buffer_elements)
      DO WHILE (memory_mb > max_memory_mb .AND. (mb > 1 .OR. rb > 1))
         trial_m = HUGE(0.0_dp)
         trial_r = HUGE(0.0_dp)
         IF (mb > 1) trial_m = batch_workspace(nao, nri, nmo, mo, ri, (mb + 1)/2, rb, buffer_elements)
         IF (rb > 1) trial_r = batch_workspace(nao, nri, nmo, mo, ri, mb, (rb + 1)/2, buffer_elements)
         IF (trial_m <= trial_r) THEN
            mb = (mb + 1)/2
         ELSE
            rb = (rb + 1)/2
         END IF
         memory_mb = MIN(trial_m, trial_r)
      END DO
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T3,A,T69,F12.3)') &
            "CC4S| Estimated memory required for export per MPI rank [MB]:", memory_mb
         WRITE (output_unit, '(T3,A,T69,F12.3)') &
            "CC4S| Maximum memory allowed for export per MPI rank [MB]:", max_memory_mb
      END IF
      IF (memory_mb > max_memory_mb) THEN
         CALL cp_abort(__LOCATION__, &
                       "CC4S: minimum transformation batch exceeds MAX_MEMORY; increase the limit or reduce OUTPUT_BUFFER_SIZE")
      END IF
      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Maximum orbitals per transformation batch:", batch_width(mo, mb)
         WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Maximum auxiliary functions per batch:", batch_width(ri, rb)
         WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Number of orbital batches:", (SIZE(mo) + mb - 1)/mb
         WRITE (output_unit, '(T3,A,T69,I12)') "CC4S| Number of auxiliary batches:", (SIZE(ri) + rb - 1)/rb
      END IF
      pdims = 0
      pdims2 = 0
      CALL dbt_pgrid_create(para_env, pdims, grid, tensor_dims=[SIZE(mo), SIZE(ri), SIZE(ao)])
      CALL dbt_pgrid_create(para_env, pdims2, grid2, tensor_dims=[SIZE(ao), SIZE(mo)])
      CALL create_2c_tensor(coeff, d1, d2, grid2, ao, mo, name="CC4S active coefficients")
      DEALLOCATE (d1, d2)
      CALL dbt_copy(coefficients, coeff)
      ! The metric and integrals can have different auxiliary block partitions.
      CALL create_2c_tensor(metric_copy, d1, d2, grid2, ri, ri, name="CC4S inverse metric")
      DEALLOCATE (d1, d2)
      CALL dbt_copy(metric, metric_copy)
      ALLOCATE (values(buffer_elements))
      vertex_norm = 0.0_dp
      iunit = -1
      start = m_walltime()
      IF (para_env%is_source()) THEN
         CALL write_eigen_energies_files(logger, print_section, nmo, homo, eigenvalues, fermi_energy)
         CALL open_coulomb_vertex_file(logger, print_section, nri, nmo, periodic, cutoff_radius, hmat, iunit)
      END IF
      elapsed(5) = elapsed(5) + m_walltime() - start
      DO im = 1, SIZE(mo), mb
         CALL dbt_create(coeff, coeff_batch, name="CC4S coefficient batch")
         CALL create_3c_tensor(gamma, d1, d2, d3, grid, ri, mo, mo, [1, 2], [3], &
                               name="CC4S Coulomb vertex batch")
         DEALLOCATE (d1, d2, d3)
         mfirst = SUM(mo(:im - 1)) + 1
         mlast = SUM(mo(:MIN(SIZE(mo), im + mb - 1)))
         bounds2(:, 1) = [1, nao]
         bounds2(:, 2) = [mfirst, mlast]
         CALL dbt_copy(coeff, coeff_batch, bounds=bounds2)
         DO ir = 1, SIZE(ri), rb
            rfirst = SUM(ri(:ir - 1)) + 1
            rlast = SUM(ri(:MIN(SIZE(ri), ir + rb - 1)))
            start = m_walltime()
            CALL timeset("cc4s_transform", handle)
            bounds3(:, 1) = [1, nao]
            bounds3(:, 2) = [rfirst, rlast]
            bounds3(:, 3) = [1, nao]
            CALL dbt_create(integrals, ao_triangle)
            CALL create_3c_tensor(half, d1, d2, d3, grid, mo, ri, ao, [1], [2, 3], &
                                  name="CC4S half transform")
            DEALLOCATE (d1, d2, d3)
            ! Copy only this auxiliary slice. Contractions never receive live SCF tensors.
            CALL dbt_copy(integrals, ao_triangle, bounds=bounds3)
            CALL dbt_contract(1.0_dp, coeff_batch, ao_triangle, 0.0_dp, half, &
                              contract_1=[1], notcontract_1=[2], contract_2=[1], notcontract_2=[2, 3], &
                              map_1=[1], map_2=[2, 3], filter_eps=0.0_dp)
            ! Restore the weighted AO triangle before transforming the second orbital index.
            ! Sum, not average: diagonal atom blocks carry weight 1/2.
            CALL dbt_contract(1.0_dp, coeff_batch, ao_triangle, 1.0_dp, half, &
                              contract_1=[1], notcontract_1=[2], contract_2=[3], notcontract_2=[2, 1], &
                              map_1=[1], map_2=[2, 3], filter_eps=0.0_dp)
            CALL dbt_destroy(ao_triangle)
            CALL create_3c_tensor(triangle, d1, d2, d3, grid, ri, mo, mo, [1, 2], [3], &
                                  name="CC4S MO batch")
            DEALLOCATE (d1, d2, d3)
            CALL dbt_contract(1.0_dp, half, coeff, 0.0_dp, triangle, &
                              contract_1=[3], notcontract_1=[1, 2], contract_2=[1], notcontract_2=[2], &
                              map_1=[3, 1], map_2=[2], filter_eps=0.0_dp)
            CALL dbt_destroy(half)
            CALL timestop(handle)
            elapsed(2) = elapsed(2) + m_walltime() - start
            start = m_walltime()
            CALL timeset("cc4s_metric", handle)
            CALL dbt_contract(1.0_dp, metric_copy, triangle, 1.0_dp, gamma, &
                              contract_1=[1], notcontract_1=[2], contract_2=[1], notcontract_2=[2, 3], &
                              map_1=[1], map_2=[2, 3], filter_eps=0.0_dp)
            CALL dbt_destroy(triangle)
            CALL timestop(handle)
            elapsed(3) = elapsed(3) + m_walltime() - start
         END DO
         CALL write_coulomb_vertex_batch(para_env, nri, nmo, gamma, iunit, mfirst, mlast, values, vertex_norm, elapsed)
         CALL dbt_destroy(gamma)
         CALL dbt_destroy(coeff_batch)
      END DO
      start = m_walltime()
      IF (para_env%is_source()) CALL cp_print_key_finished_output(iunit, logger, print_section)
      elapsed(5) = elapsed(5) + m_walltime() - start
      DEALLOCATE (values)
      CALL dbt_destroy(metric_copy)
      CALL dbt_destroy(coeff)
      CALL dbt_pgrid_destroy(grid)
      CALL dbt_pgrid_destroy(grid2)
      DEALLOCATE (ao, ri, mo)
   END SUBROUTINE export_vertex_batches

! **************************************************************************************************
!> \brief Largest number of functions in a group of consecutive tensor blocks.
!> \param sizes Block sizes
!> \param blocks Number of blocks per batch
!> \return Maximum batch width
! **************************************************************************************************
   FUNCTION batch_width(sizes, blocks) RESULT(width)
      INTEGER, INTENT(IN)                                :: sizes(:), blocks
      INTEGER                                            :: width

      INTEGER                                            :: i

      width = 0
      DO i = 1, SIZE(sizes), blocks
         width = MAX(width, SUM(sizes(i:MIN(SIZE(sizes), i + blocks - 1))))
      END DO
   END FUNCTION batch_width

! **************************************************************************************************
!> \brief Conservative busiest-rank workspace estimate, without a balanced-distribution assumption.
!>        Four copies allow for private operands, contraction results and redistribution buffers.
!>        Excludes existing SCF storage, allocator overhead and CC4S. This is not an RSS limit.
!> \param nao Number of AO functions
!> \param nri Number of auxiliary functions
!> \param nmo Number of active orbitals
!> \param mo Orbital block sizes
!> \param ri Auxiliary block sizes
!> \param mb Orbital blocks per batch
!> \param rb Auxiliary blocks per batch
!> \param buffer_elements Reusable output buffer length
!> \return Estimated additional memory per rank in decimal MB
! **************************************************************************************************
   FUNCTION batch_workspace(nao, nri, nmo, mo, ri, mb, rb, buffer_elements) RESULT(memory_mb)
      INTEGER, INTENT(IN)                                :: nao, nri, nmo, mo(:), ri(:), mb, rb
      INTEGER(KIND=int_8), INTENT(IN)                    :: buffer_elements
      REAL(dp)                                           :: memory_mb

      REAL(dp)                                           :: a, b, fixed, m, q, r, work

      a = REAL(nao, dp)
      r = REAL(nri, dp)
      m = REAL(nmo, dp)
      b = REAL(batch_width(mo, mb), dp)
      q = REAL(batch_width(ri, rb), dp)
      fixed = 4.0_dp*(r*r + a*m) + 2.0_dp*REAL(buffer_elements, dp)
      work = 4.0_dp*(a*a*q + b*q*a + q*m*b + r*m*b + a*b)
      memory_mb = REAL(dp_size, dp)*MAX(6.0_dp*a*a, fixed + work)/1.0e6_dp
   END FUNCTION batch_workspace

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
!> \brief Gather and write one orbital range in bounded pieces, preserving Fortran element order.
!> \param para_env MPI parallel environment
!> \param nRI Number of auxiliary functions
!> \param nMO Number of active canonical orbitals
!> \param gamma Distributed vertex, with data only in the selected third-index range
!> \param iunit Open stream unit on the writing rank
!> \param first_mo First orbital in this batch
!> \param last_mo Last orbital in this batch
!> \param values Reusable output buffer
!> \param vertex_norm Accumulated squared norm on the writing rank
!> \param elapsed Accumulated wall times
! **************************************************************************************************
   SUBROUTINE write_coulomb_vertex_batch(para_env, nRI, nMO, gamma, iunit, first_mo, last_mo, &
                                         values, vertex_norm, elapsed)
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: nRI, nMO
      TYPE(dbt_type), INTENT(INOUT)                      :: gamma
      INTEGER, INTENT(IN)                                :: iunit, first_mo, last_mo
      REAL(dp), INTENT(INOUT)                            :: values(:), vertex_norm, elapsed(5)

      INTEGER                                            :: column, handle, hi, ind(3), j, k, lo, &
                                                            sizes(3)
      INTEGER(KIND=int_8)                                :: count, first, last, mirror, nvalues, &
                                                            offset, row, segment, start, used, &
                                                            within
      INTEGER, ALLOCATABLE                               :: off1(:), off2(:), off3(:)
      LOGICAL                                            :: found
      REAL(dp)                                           :: start_time
      REAL(dp), ALLOCATABLE                              :: block(:, :, :)
      TYPE(dbt_iterator_type)                            :: iter

      start_time = m_walltime()
      CALL timeset("cc4s_gather", handle)
      ALLOCATE (off1(dbt_nblks_total(gamma, 1)), off2(dbt_nblks_total(gamma, 2)), &
                off3(dbt_nblks_total(gamma, 3)))
      CALL dbt_get_info(gamma, blk_offset_1=off1, blk_offset_2=off2, blk_offset_3=off3)
      DO column = first_mo, last_mo
         offset = INT(nRI, int_8)*INT(nMO, int_8)*INT(column - 1, int_8)
         nvalues = offset + INT(nRI, int_8)*INT(column, int_8)
         ! Gather one triangle only. Write its transpose from identical values for exact symmetry.
         DO WHILE (offset < nvalues)
            count = MIN(SIZE(values, KIND=int_8), nvalues - offset)
            values = 0.0_dp
            ! Blocks own disjoint elements; threads fill disjoint parts of the bounded buffer.
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(gamma, off1, off2, off3, nRI, nMO, offset, count, values) &
!$OMP PRIVATE(iter, ind, sizes, block, found, j, k, first, last, start, lo, hi)
            CALL dbt_iterator_start(iter, gamma)
            DO WHILE (dbt_iterator_blocks_left(iter))
               CALL dbt_iterator_next_block(iter, ind, blk_size=sizes)
               first = INT(off1(ind(1)) - 1, int_8) + INT(nRI, int_8)* &
                       (INT(off2(ind(2)) - 1, int_8) + INT(nMO, int_8)*INT(off3(ind(3)) - 1, int_8))
               last = first + INT(sizes(1), int_8) + INT(nRI, int_8)* &
                      (INT(sizes(2) - 1, int_8) + INT(nMO, int_8)*INT(sizes(3) - 1, int_8))
               IF (first >= offset + count .OR. last <= offset) CYCLE
               CALL dbt_get_block(gamma, ind, block, found)
               CPASSERT(found)
               DO k = 1, sizes(3)
                  DO j = 1, sizes(2)
                     start = first + INT(nRI, int_8)*(INT(j - 1, int_8) + INT(nMO, int_8)*INT(k - 1, int_8))
                     lo = INT(MAX(0_int_8, offset - start)) + 1
                     hi = INT(MIN(INT(sizes(1), int_8), offset + count - start))
                     IF (lo > hi) CYCLE
                     values(start + lo - offset:start + hi - offset) = block(lo:hi, j, k)
                  END DO
               END DO
               DEALLOCATE (block)
            END DO
            CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
            CALL para_env%sum(values(:count), para_env%source)
            CALL timestop(handle)
            elapsed(4) = elapsed(4) + m_walltime() - start_time
            start_time = m_walltime()
            CALL timeset("cc4s_write_files", handle)
            IF (para_env%is_source()) THEN
               CALL write_coulomb_vertex_block(iunit, offset, values(:count))
               vertex_norm = vertex_norm + SUM(values(:count)**2)
               used = 0_int_8
               DO WHILE (used < count)
                  within = MODULO(offset + used, INT(nRI, int_8)*INT(nMO, int_8))
                  row = within/INT(nRI, int_8) + 1_int_8
                  within = MODULO(within, INT(nRI, int_8))
                  segment = MIN(count - used, INT(nRI, int_8) - within)
                  IF (row < column) THEN
                     mirror = within + INT(nRI, int_8)*(INT(column - 1, int_8) + INT(nMO, int_8)*(row - 1_int_8))
                     CALL write_coulomb_vertex_block(iunit, mirror, values(used + 1:used + segment))
                     vertex_norm = vertex_norm + SUM(values(used + 1:used + segment)**2)
                  END IF
                  used = used + segment
               END DO
            END IF
            CALL timestop(handle)
            elapsed(5) = elapsed(5) + m_walltime() - start_time
            start_time = m_walltime()
            CALL timeset("cc4s_gather", handle)
            offset = offset + count
         END DO
      END DO
      DEALLOCATE (off1, off2, off3)
      CALL timestop(handle)
      elapsed(4) = elapsed(4) + m_walltime() - start_time
   END SUBROUTINE write_coulomb_vertex_batch

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
!> \param nocc Number of exported occupied spatial orbitals
!> \param max_iter Maximum coupled-cluster iterations
!> \param energy_conv CC4S energy convergence threshold
!> \param amp_conv CC4S amplitude convergence threshold
! **************************************************************************************************
   SUBROUTINE write_cc4s_in(logger, print_section, method, nocc, max_iter, energy_conv, amp_conv)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: print_section
      CHARACTER(LEN=*), INTENT(IN)                       :: method
      INTEGER, INTENT(IN)                                :: nocc, max_iter
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
         CALL write_coupled_cluster_step(iunit, 'Ccsd', 'Amplitudes', nocc, max_iter, energy_conv, amp_conv)
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
!> \param nocc Number of exported occupied spatial orbitals
!> \param max_iter Maximum coupled-cluster iterations
!> \param energy_conv CC4S energy convergence threshold
!> \param amp_conv CC4S amplitude convergence threshold
! **************************************************************************************************
   SUBROUTINE write_coupled_cluster_step(iunit, cc_method, amp_key, nocc, max_iter, energy_conv, amp_conv)
      INTEGER, INTENT(IN)                                :: iunit
      CHARACTER(LEN=*), INTENT(IN)                       :: cc_method, amp_key
      INTEGER, INTENT(IN)                                :: nocc, max_iter
      REAL(dp), INTENT(IN)                               :: energy_conv, amp_conv

      WRITE (iunit, '(A)') '- name: CoupledCluster'
      WRITE (iunit, '(A)') '  in:'
      WRITE (iunit, '(A)') '    method: '//TRIM(cc_method)
      ! Match CC4S v1.1 defaults using the occupied space after frozen-core removal.
      WRITE (iunit, '(A,I0)') '    integralsSliceSize: ', nocc
      WRITE (iunit, '(A)') '    slicedEigenEnergies: EigenEnergies'
      WRITE (iunit, '(A)') '    coulombIntegrals: CoulombIntegrals'
      WRITE (iunit, '(A)') '    slicedCoulombVertex: CoulombVertex'
      WRITE (iunit, '(A,I0)') '    maxIterations: ', max_iter
      WRITE (iunit, '(A,ES25.16)') '    energyConvergence: ', energy_conv
      WRITE (iunit, '(A,ES25.16)') '    amplitudesConvergence: ', amp_conv
      WRITE (iunit, '(A)') '    mixer:'
      WRITE (iunit, '(A)') '      type: DiisMixer'
      WRITE (iunit, '(A)') '      maxResidua: 4'
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
