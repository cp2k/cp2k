!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Does all kind of post scf calculations for DFTB
!> \par History
!>      Started as a copy from the GPW file
!> \author JHU (03.2013)
! **************************************************************************************************
MODULE qs_scf_post_tb
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_reduce,&
                                              cp_fm_cholesky_restore
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_init_random,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE input_constants,                 ONLY: ot_precond_full_all
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_ivals,&
                                              section_get_lval,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: twopi
   USE memory_utilities,                ONLY: reallocate
   USE molden_utils,                    ONLY: write_mos_molden
   USE moments_utils,                   ONLY: get_reference_point
   USE mulliken,                        ONLY: mulliken_charges
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: debye
   USE population_analyses,             ONLY: lowdin_population_analysis
   USE preconditioner_types,            ONLY: preconditioner_type
   USE pw_env_methods,                  ONLY: pw_env_create,&
                                              pw_env_rebuild
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_release,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_types,                ONLY: do_ewald_none,&
                                              greens_fn_type,&
                                              pw_green_create,&
                                              pw_green_release,&
                                              pw_poisson_analytic,&
                                              pw_poisson_parameter_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_rho_core,&
                                              calculate_rho_elec,&
                                              calculate_wavefunction
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
   USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
   USE qs_dos,                          ONLY: calculate_dos,&
                                              calculate_dos_kp
   USE qs_elf_methods,                  ONLY: qs_elf_calc
   USE qs_energy_window,                ONLY: energy_windows
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_io,                        ONLY: write_mo_set
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues,&
                                              make_mo_eig
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_pdos,                         ONLY: calculate_projected_dos
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_csr_write,                ONLY: write_ks_matrix_csr,&
                                              write_s_matrix_csr
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
   USE qs_scf_wfn_mix,                  ONLY: wfn_mix
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE scf_control_types,               ONLY: scf_control_type
   USE stm_images,                      ONLY: th_stm_image
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              task_list_type
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_tb'
   PUBLIC :: scf_post_calculation_tb, make_lumo_tb

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!> \param qs_env ...
!> \param tb_type ...
!> \param no_mos ...
!> \par History
!>      03.2013 copy of scf_post_gpw
!> \author JHU
!> \note
! **************************************************************************************************
   SUBROUTINE scf_post_calculation_tb(qs_env, tb_type, no_mos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*)                                   :: tb_type
      LOGICAL, INTENT(IN)                                :: no_mos

      CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_tb', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=6)                                   :: ana
      CHARACTER(LEN=default_string_length)               :: aname
      INTEGER                                            :: handle, homo, iat, iatom, ikind, ispin, &
                                                            nat, natom, nkind, nlumo_stm, nlumos, &
                                                            nspins, print_level, unit_nr
      LOGICAL                                            :: do_cube, do_kpoints, explicit, has_homo, &
                                                            print_it, rebuild
      REAL(KIND=dp)                                      :: zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mcharge
      REAL(KIND=dp), DIMENSION(2, 2)                     :: homo_lumo
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: charges
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: unoccupied_evals_stm
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: unoccupied_orbs_stm
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_deriv
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, moments_section, print_key, &
                                                            print_section, wfn_mix_section
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      NULLIFY (dft_control, rho, para_env, matrix_s, matrix_p)
      CALL get_qs_env(qs_env, scf_env=scf_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      dft_control=dft_control, rho=rho, natom=natom, para_env=para_env, &
                      particle_set=particle_set, do_kpoints=do_kpoints, matrix_s_kp=matrix_s)
      nspins = dft_control%nspins
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ! Mulliken charges
      ALLOCATE (charges(natom, nspins), mcharge(natom))
      !
      CALL mulliken_charges(matrix_p, matrix_s, para_env, charges)
      !
      nkind = SIZE(atomic_kind_set)
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
         SELECT CASE (TRIM(tb_type))
         CASE ("DFTB")
            CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
            CALL get_dftb_atom_param(dftb_kind, zeff=zeff)
         CASE ("xTB")
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, zeff=zeff)
         CASE DEFAULT
            CPABORT("unknown TB type")
         END SELECT
         DO iatom = 1, nat
            iat = atomic_kind_set(ikind)%atom_list(iatom)
            mcharge(iat) = zeff - SUM(charges(iat, 1:nspins))
         END DO
      END DO

      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")
      print_section => section_vals_get_subs_vals(dft_section, "PRINT")

      ! Mulliken
      print_key => section_vals_get_subs_vals(print_section, "MULLIKEN")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_section, "MULLIKEN", &
                                        extension=".mulliken", log_filename=.FALSE.)
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(/,/,T2,A)") "MULLIKEN POPULATION ANALYSIS"
            IF (nspins == 1) THEN
               WRITE (UNIT=unit_nr, FMT="(/,T2,A,T70,A)") &
                  " # Atom   Element   Kind        Atomic population", " Net charge"
               DO ikind = 1, nkind
                  CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
                  aname = ""
                  SELECT CASE (tb_type)
                  CASE ("DFTB")
                     CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
                     CALL get_dftb_atom_param(dftb_kind, name=aname)
                  CASE ("xTB")
                     CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
                     CALL get_xtb_atom_param(xtb_kind, symbol=aname)
                  CASE DEFAULT
                     CPABORT("unknown TB type")
                  END SELECT
                  ana = ADJUSTR(TRIM(ADJUSTL(aname)))
                  DO iatom = 1, nat
                     iat = atomic_kind_set(ikind)%atom_list(iatom)
                     WRITE (UNIT=unit_nr, &
                            FMT="(T2,I7,5X,A6,I6,T39,F12.6,T69,F12.6)") &
                        iat, ADJUSTL(ana), ikind, charges(iat, 1), mcharge(iat)
                  END DO
               END DO
               WRITE (UNIT=unit_nr, &
                      FMT="(T2,A,T39,F12.6,T69,F12.6,/)") &
                  "# Total charge", SUM(charges(:, 1)), SUM(mcharge(:))
            ELSE
               WRITE (UNIT=unit_nr, FMT="(/,T2,A)") &
                  "# Atom  Element  Kind  Atomic population (alpha,beta)   Net charge  Spin moment"
               DO ikind = 1, nkind
                  CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
                  aname = ""
                  SELECT CASE (tb_type)
                  CASE ("DFTB")
                     CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
                     CALL get_dftb_atom_param(dftb_kind, name=aname)
                  CASE ("xTB")
                     CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
                     CALL get_xtb_atom_param(xtb_kind, symbol=aname)
                  CASE DEFAULT
                     CPABORT("unknown TB type")
                  END SELECT
                  ana = ADJUSTR(TRIM(ADJUSTL(aname)))
                  DO iatom = 1, nat
                     iat = atomic_kind_set(ikind)%atom_list(iatom)
                     WRITE (UNIT=unit_nr, &
                            FMT="(T2,I6,3X,A6,I6,T29,4(1X,F12.6))") &
                        iat, ADJUSTL(ana), ikind, charges(iat, 1:2), mcharge(iat), &
                        charges(iat, 1) - charges(iat, 2)
                  END DO
               END DO
               WRITE (UNIT=unit_nr, &
                      FMT="(T2,A,T29,4(1X,F12.6),/)") &
                  "# Total charge and spin", SUM(charges(:, 1)), SUM(charges(:, 2)), SUM(mcharge(:))
            END IF
            CALL m_flush(unit_nr)
         END IF
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      ! Compute the Lowdin charges
      print_key => section_vals_get_subs_vals(print_section, "LOWDIN")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         SELECT CASE (tb_type)
         CASE ("DFTB")
            CPWARN("Lowdin population analysis not implemented for DFTB method.")
         CASE ("xTB")
            unit_nr = cp_print_key_unit_nr(logger, print_section, "LOWDIN", extension=".lowdin", &
                                           log_filename=.FALSE.)
            print_level = 1
            CALL section_vals_val_get(print_key, "PRINT_GOP", l_val=print_it)
            IF (print_it) print_level = 2
            CALL section_vals_val_get(print_key, "PRINT_ALL", l_val=print_it)
            IF (print_it) print_level = 3
            IF (do_kpoints) THEN
               CPWARN("Lowdin charges not implemented for k-point calculations!")
            ELSE
               CALL lowdin_population_analysis(qs_env, unit_nr, print_level)
            END IF
            CALL cp_print_key_finished_output(unit_nr, logger, print_section, "LOWDIN")
         CASE DEFAULT
            CPABORT("unknown TB type")
         END SELECT
      END IF

      ! Hirshfeld
      print_key => section_vals_get_subs_vals(print_section, "HIRSHFELD")
      CALL section_vals_get(print_key, explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CPWARN("Hirshfeld charges not available for TB methods.")
         END IF
      END IF

      ! MAO
      print_key => section_vals_get_subs_vals(print_section, "MAO_ANALYSIS")
      CALL section_vals_get(print_key, explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CPWARN("MAO analysis not available for TB methods.")
         END IF
      END IF

      ! Dipole Moments
      print_key => section_vals_get_subs_vals(print_section, "MOMENTS")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_section, "MOMENTS", &
                                        extension=".data", middle_name="tb_dipole", log_filename=.FALSE.)
         moments_section => section_vals_get_subs_vals(print_section, "MOMENTS")
         CALL tb_dipole(qs_env, moments_section, unit_nr, mcharge)
         CALL cp_print_key_finished_output(unit_nr, logger, print_key)
      END IF

      DEALLOCATE (charges, mcharge)

      ! MO
      IF (.NOT. no_mos) THEN
         print_key => section_vals_get_subs_vals(print_section, "MO")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CALL print_mos(qs_env, dft_section, tb_type)
         END IF
      END IF

      ! Wavefunction mixing
      IF (.NOT. no_mos) THEN
         wfn_mix_section => section_vals_get_subs_vals(dft_section, "PRINT%WFN_MIX")
         CALL section_vals_get(wfn_mix_section, explicit=explicit)
         IF (explicit .AND. .NOT. qs_env%run_rtp) CALL wfn_mix_tb(qs_env, dft_section, scf_env)
      END IF

      IF (.NOT. no_mos) THEN
         print_key => section_vals_get_subs_vals(print_section, "DOS")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            IF (do_kpoints) THEN
               CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
               CALL calculate_dos_kp(kpoints, qs_env, dft_section)
            ELSE
               CALL get_qs_env(qs_env, mos=mos)
               CALL calculate_dos(mos, qs_env, dft_section)
            ENDIF
         ENDIF
      ENDIF

      ! PDOS
      IF (.NOT. no_mos) THEN
         print_key => section_vals_get_subs_vals(print_section, "PDOS")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            IF (do_kpoints) THEN
               CPWARN("Projected density of states not implemented for k-points.")
            ELSE
               CALL get_qs_env(qs_env, mos=mos, matrix_ks=ks_rmpv)
               DO ispin = 1, dft_control%nspins
                  IF (scf_env%method == ot_method_nr) THEN
                     CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                                     eigenvalues=mo_eigenvalues)
                     IF (ASSOCIATED(qs_env%mo_derivs)) THEN
                        mo_coeff_deriv => qs_env%mo_derivs(ispin)%matrix
                     ELSE
                        mo_coeff_deriv => NULL()
                     ENDIF
                     CALL calculate_subspace_eigenvalues(mo_coeff, ks_rmpv(ispin)%matrix, mo_eigenvalues, &
                                                         do_rotation=.TRUE., &
                                                         co_rotate_dbcsr=mo_coeff_deriv)
                     CALL set_mo_occupation(mo_set=mos(ispin)%mo_set)
                  END IF
                  IF (dft_control%nspins == 2) THEN
                     CALL calculate_projected_dos(mos(ispin)%mo_set, atomic_kind_set, &
                                                  qs_kind_set, particle_set, qs_env, dft_section, ispin=ispin)
                  ELSE
                     CALL calculate_projected_dos(mos(ispin)%mo_set, atomic_kind_set, &
                                                  qs_kind_set, particle_set, qs_env, dft_section)
                  END IF
               ENDDO
            ENDIF
         END IF
      END IF

      ! can we do CUBE files?
      SELECT CASE (tb_type)
      CASE ("DFTB")
         do_cube = .FALSE.
         rebuild = .FALSE.
      CASE ("xTB")
         do_cube = .TRUE.
         rebuild = .TRUE.
      CASE DEFAULT
         CPABORT("unknown TB type")
      END SELECT

      ! Energy Windows for LS code
      print_key => section_vals_get_subs_vals(print_section, "ENERGY_WINDOWS")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (do_cube) THEN
            IF (do_kpoints) THEN
               CPWARN("Energy Windows not implemented for k-points.")
            ELSE
               IF (rebuild) THEN
                  CALL rebuild_pw_env(qs_env)
                  rebuild = .FALSE.
               END IF
               CALL energy_windows(qs_env)
            END IF
         ELSE
            CPWARN("Energy Windows not implemented for TB methods.")
         END IF
      ENDIF

      ! DENSITY CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "E_DENSITY_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (do_cube) THEN
            IF (rebuild) THEN
               CALL rebuild_pw_env(qs_env)
               rebuild = .FALSE.
            END IF
            CALL print_e_density(qs_env, print_key)
         ELSE
            CPWARN("Electronic density cube file not implemented for TB methods.")
         END IF
      END IF

      ! TOTAL DENSITY CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "TOT_DENSITY_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (do_cube) THEN
            IF (rebuild) THEN
               CALL rebuild_pw_env(qs_env)
               rebuild = .FALSE.
            END IF
            CALL print_density_cubes(qs_env, print_key, total_density=.TRUE.)
         ELSE
            CPWARN("Total density cube file not implemented for TB methods.")
         END IF
      END IF

      ! V_Hartree CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "V_HARTREE_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (do_cube) THEN
            IF (rebuild) THEN
               CALL rebuild_pw_env(qs_env)
               rebuild = .FALSE.
            END IF
            CALL print_density_cubes(qs_env, print_key, v_hartree=.TRUE.)
         ELSE
            CPWARN("Hartree potential cube file not implemented for TB methods.")
         END IF
      END IF

      ! EFIELD CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "EFIELD_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (do_cube) THEN
            IF (rebuild) THEN
               CALL rebuild_pw_env(qs_env)
               rebuild = .FALSE.
            END IF
            CALL print_density_cubes(qs_env, print_key, efield=.TRUE.)
         ELSE
            CPWARN("Efield cube file not implemented for TB methods.")
         END IF
      END IF

      ! ELF
      print_key => section_vals_get_subs_vals(print_section, "ELF_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (do_cube) THEN
            IF (rebuild) THEN
               CALL rebuild_pw_env(qs_env)
               rebuild = .FALSE.
            END IF
            CALL print_elf(qs_env, print_key)
         ELSE
            CPWARN("ELF not implemented for TB methods.")
         END IF
      END IF

      ! MO CUBES
      IF (.NOT. no_mos) THEN
         print_key => section_vals_get_subs_vals(print_section, "MO_CUBES")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            IF (do_cube) THEN
               IF (rebuild) THEN
                  CALL rebuild_pw_env(qs_env)
                  rebuild = .FALSE.
               END IF
               CALL print_mo_cubes(qs_env, print_key)
            ELSE
               CPWARN("Printing of MO cube files not implemented for TB methods.")
            END IF
         END IF
      END IF

      ! STM
      IF (.NOT. no_mos) THEN
         print_key => section_vals_get_subs_vals(print_section, "STM")
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            IF (do_cube) THEN
               IF (rebuild) THEN
                  CALL rebuild_pw_env(qs_env)
                  rebuild = .FALSE.
               END IF
               IF (do_kpoints) THEN
                  CPWARN("STM not implemented for k-point calculations!")
               ELSE
                  nlumo_stm = section_get_ival(print_key, "NLUMO")
                  CPASSERT(.NOT. dft_control%restricted)
                  CALL get_qs_env(qs_env, mos=mos, mo_derivs=mo_derivs, &
                                  scf_control=scf_control, matrix_ks=ks_rmpv)
                  CALL make_mo_eig(mos, dft_control%nspins, ks_rmpv, scf_control, mo_derivs)
                  DO ispin = 1, dft_control%nspins
                     CALL get_mo_set(mo_set=mos(ispin)%mo_set, eigenvalues=mo_eigenvalues, homo=homo)
                     homo_lumo(ispin, 1) = mo_eigenvalues(homo)
                  END DO
                  has_homo = .TRUE.
                  IF (nlumo_stm > 0) THEN
                     ALLOCATE (unoccupied_orbs_stm(dft_control%nspins))
                     ALLOCATE (unoccupied_evals_stm(dft_control%nspins))
                     CALL make_lumo_tb(qs_env, scf_env, unoccupied_orbs_stm, unoccupied_evals_stm, &
                                       nlumo_stm, nlumos)
                  END IF

                  CALL get_qs_env(qs_env, subsys=subsys)
                  CALL qs_subsys_get(subsys, particles=particles)
                  CALL th_stm_image(qs_env, print_key, particles, unoccupied_orbs_stm, &
                                    unoccupied_evals_stm)

                  IF (nlumo_stm > 0) THEN
                     DO ispin = 1, dft_control%nspins
                        DEALLOCATE (unoccupied_evals_stm(ispin)%array)
                        CALL cp_fm_release(unoccupied_orbs_stm(ispin)%matrix)
                     ENDDO
                     DEALLOCATE (unoccupied_evals_stm)
                     DEALLOCATE (unoccupied_orbs_stm)
                  END IF
               END IF
            END IF
         END IF
      END IF

      ! these print keys are not supported in TB

      ! V_XC CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "V_XC_CUBE")
      CALL section_vals_get(print_key, explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CPWARN("XC potential cube file not available for TB methods.")
         END IF
      END IF

      ! Electric field gradients
      print_key => section_vals_get_subs_vals(print_section, "ELECTRIC_FIELD_GRADIENT")
      CALL section_vals_get(print_key, explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CPWARN("Electric field gradient not implemented for TB methods.")
         END IF
      END IF

      ! KINETIC ENERGY
      print_key => section_vals_get_subs_vals(print_section, "KINETIC_ENERGY")
      CALL section_vals_get(print_key, explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CPWARN("Kinetic energy not available for TB methods.")
         END IF
      END IF

      ! Xray diffraction spectrum
      print_key => section_vals_get_subs_vals(print_section, "XRAY_DIFFRACTION_SPECTRUM")
      CALL section_vals_get(print_key, explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CPWARN("Xray diffraction spectrum not implemented for TB methods.")
         END IF
      END IF

      ! EPR Hyperfine Coupling
      print_key => section_vals_get_subs_vals(print_section, "HYPERFINE_COUPLING_TENSOR")
      CALL section_vals_get(print_key, explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CPWARN("Hyperfine Coupling not implemented for TB methods.")
         END IF
      END IF

      ! PLUS_U
      print_key => section_vals_get_subs_vals(print_section, "PLUS_U")
      CALL section_vals_get(print_key, explicit=explicit)
      IF (explicit) THEN
         IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
            CPWARN("DFT+U method not implemented for TB methods.")
         END IF
      END IF

      CALL write_ks_matrix_csr(qs_env, qs_env%input)
      CALL write_s_matrix_csr(qs_env, qs_env%input)

      CALL timestop(handle)

   END SUBROUTINE scf_post_calculation_tb

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input ...
!> \param unit_nr ...
!> \param charges ...
! **************************************************************************************************
   SUBROUTINE tb_dipole(qs_env, input, unit_nr, charges)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input
      INTEGER, INTENT(in)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(:), INTENT(in)            :: charges

      CHARACTER(len=*), PARAMETER :: routineN = 'tb_dipole', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: description, dipole_type
      COMPLEX(KIND=dp)                                   :: dzeta, dzphase(3), zeta, zphase(3)
      COMPLEX(KIND=dp), DIMENSION(3)                     :: dggamma, ggamma
      INTEGER                                            :: i, iat, ikind, j, nat, reference
      LOGICAL                                            :: do_berry
      REAL(KIND=dp) :: charge_tot, ci(3), dci(3), dipole(3), dipole_deriv(3), drcc(3), dria(3), &
         dtheta, gvec(3), q, rcc(3), ria(3), theta, tmp(3), via(3)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (atomic_kind_set, cell, results)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, cell=cell, results=results)

      ! Reference point
      reference = section_get_ival(input, keyword_name="REFERENCE")
      NULLIFY (ref_point)
      description = '[DIPOLE]'
      CALL section_vals_val_get(input, "REF_POINT", r_vals=ref_point)
      CALL section_vals_val_get(input, "PERIODIC", l_val=do_berry)

      CALL get_reference_point(rcc, drcc, qs_env=qs_env, reference=reference, ref_point=ref_point)

      ! Dipole deriv will be the derivative of the Dipole(dM/dt=\sum e_j v_j)
      dipole_deriv = 0.0_dp
      dipole = 0.0_dp
      IF (do_berry) THEN
         dipole_type = "[BERRY PHASE]"
         rcc = pbc(rcc, cell)
         charge_tot = 0._dp
         charge_tot = SUM(charges)
         ria = twopi*MATMUL(cell%h_inv, rcc)
         zphase = CMPLX(COS(ria), SIN(ria), dp)**charge_tot

         dria = twopi*MATMUL(cell%h_inv, drcc)
         dzphase = charge_tot*CMPLX(-SIN(ria), COS(ria), dp)**(charge_tot - 1.0_dp)*dria

         ggamma = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
         dggamma = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
            DO i = 1, nat
               iat = atomic_kind_set(ikind)%atom_list(i)
               ria = particle_set(iat)%r(:)
               ria = pbc(ria, cell)
               via = particle_set(iat)%v(:)
               q = charges(iat)
               DO j = 1, 3
                  gvec = twopi*cell%h_inv(j, :)
                  theta = SUM(ria(:)*gvec(:))
                  dtheta = SUM(via(:)*gvec(:))
                  zeta = CMPLX(COS(theta), SIN(theta), KIND=dp)**(-q)
                  dzeta = -q*CMPLX(-SIN(theta), COS(theta), KIND=dp)**(-q - 1.0_dp)*dtheta
                  dggamma(j) = dggamma(j)*zeta + ggamma(j)*dzeta
                  ggamma(j) = ggamma(j)*zeta
               END DO
            ENDDO
         END DO
         dggamma = dggamma*zphase + ggamma*dzphase
         ggamma = ggamma*zphase
         IF (ALL(REAL(ggamma, KIND=dp) /= 0.0_dp)) THEN
            tmp = AIMAG(ggamma)/REAL(ggamma, KIND=dp)
            ci = ATAN(tmp)
            dci = (1.0_dp/(1.0_dp + tmp**2))* &
                  (AIMAG(dggamma)*REAL(ggamma, KIND=dp) - AIMAG(ggamma)*REAL(dggamma, KIND=dp))/(REAL(ggamma, KIND=dp))**2
            dipole = MATMUL(cell%hmat, ci)/twopi
            dipole_deriv = MATMUL(cell%hmat, dci)/twopi
         END IF
      ELSE
         dipole_type = "[Non Periodic]"
         DO i = 1, SIZE(particle_set)
            ! no pbc(particle_set(i)%r(:),cell) so that the total dipole is the sum of the molecular dipoles
            ria = particle_set(i)%r(:)
            q = charges(i)
            dipole = dipole - q*(ria - rcc)
            dipole_deriv(:) = dipole_deriv(:) - q*(particle_set(i)%v(:) - drcc)
         END DO
      END IF
      CALL cp_results_erase(results=results, description=description)
      CALL put_results(results=results, description=description, &
                       values=dipole(1:3))
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, '(1X,A,T48,3F11.6)') "DIPOLE "//TRIM(dipole_type)//"(A.U.)|", dipole
         WRITE (unit_nr, '(1X,A,T48,3F11.6)') "DIPOLE "//TRIM(dipole_type)//"(Debye)|", dipole*debye
         WRITE (unit_nr, '(1X,A,T48,3F11.6)') "DIPOLE "//TRIM(dipole_type)//" DERIVATIVE(A.U.)|", dipole_deriv
      END IF

   END SUBROUTINE tb_dipole

! **************************************************************************************************
!> \brief computes the MOs and calls the wavefunction mixing routine.
!> \param qs_env ...
!> \param dft_section ...
!> \param scf_env ...
!> \author Florian Schiffmann
!> \note
! **************************************************************************************************

   SUBROUTINE wfn_mix_tb(qs_env, dft_section, scf_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'wfn_mix_tb', routineP = moduleN//':'//routineN

      INTEGER                                            :: iounit, ispin, nao, nmo
      REAL(dp), DIMENSION(:), POINTER                    :: mo_eigenvalues
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: lumos
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fmstruct, ao_lumo_struct
      TYPE(cp_fm_type), POINTER                          :: KS_tmp, mo_coeff, MO_tmp, S_tmp, work
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: wfn_mix_section

      logger => cp_get_default_logger()
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s, matrix_ks=matrix_ks, &
                      particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, mos=mos)

      wfn_mix_section => section_vals_get_subs_vals(dft_section, "PRINT%WFN_MIX")

      CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)

      CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct, nrow_global=nao, ncol_global=nao, &
                               template_fmstruct=mo_coeff%matrix_struct)
      CALL cp_fm_create(S_tmp, matrix_struct=ao_ao_fmstruct)
      CALL cp_fm_create(KS_tmp, matrix_struct=ao_ao_fmstruct)
      CALL cp_fm_create(MO_tmp, matrix_struct=ao_ao_fmstruct)
      CALL cp_fm_create(work, matrix_struct=ao_ao_fmstruct)
      ALLOCATE (lumos(SIZE(mos)))

      CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, S_tmp)
      CALL cp_fm_cholesky_decompose(S_tmp)

      DO ispin = 1, SIZE(mos)
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, eigenvalues=mo_eigenvalues, nmo=nmo)
         CALL cp_fm_struct_create(fmstruct=ao_lumo_struct, nrow_global=nao, ncol_global=nao - nmo, &
                                  template_fmstruct=mo_coeff%matrix_struct)

         CALL cp_fm_create(lumos(ispin)%matrix, matrix_struct=ao_lumo_struct)
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, KS_tmp)
         CALL cp_fm_cholesky_reduce(KS_tmp, S_tmp)
         CALL choose_eigv_solver(KS_tmp, work, mo_eigenvalues)
         CALL cp_fm_cholesky_restore(work, nao, S_tmp, MO_tmp, "SOLVE")
         CALL cp_fm_to_fm_submat(MO_tmp, mo_coeff, nao, nmo, 1, 1, 1, 1)
         CALL cp_fm_to_fm_submat(MO_tmp, lumos(ispin)%matrix, nao, nao - nmo, 1, nmo + 1, 1, 1)

         CALL cp_fm_struct_release(ao_lumo_struct)
      END DO

      iounit = cp_logger_get_default_io_unit(logger)
      CALL wfn_mix(mos, particle_set, dft_section, qs_kind_set, &
                   lumos, scf_env, matrix_s, iounit)

      DO ispin = 1, SIZE(mos)
         CALL cp_fm_release(lumos(ispin)%matrix)
      END DO
      DEALLOCATE (lumos)
      CALL cp_fm_release(S_tmp)
      CALL cp_fm_release(MO_tmp)
      CALL cp_fm_release(KS_tmp)
      CALL cp_fm_release(work)
      CALL cp_fm_struct_release(ao_ao_fmstruct)

   END SUBROUTINE wfn_mix_tb

! **************************************************************************************************
!> \brief MO printing
!> \param qs_env ...
!> \param dft_section ...
!> \param tb_type ...
!> \par History
!> \author JHU
!> \note
! **************************************************************************************************
   SUBROUTINE print_mos(qs_env, dft_section, tb_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: dft_section
      CHARACTER(LEN=*)                                   :: tb_type

      CHARACTER(len=*), PARAMETER :: routineN = 'print_mos', routineP = moduleN//':'//routineN

      INTEGER                                            :: ik, iounit
      LOGICAL                                            :: do_kpoints
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(mo_set_p_type), DIMENSION(:, :), POINTER      :: mos_k
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: sprint_section

      CPASSERT(ASSOCIATED(qs_env))
      IF (qs_env%run_rtp) RETURN

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set)

      IF (.NOT. do_kpoints) THEN
         CALL get_qs_env(qs_env, mos=mos)
         IF (dft_control%nspins == 2) THEN
            CALL write_mo_set(mos(1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                              dft_section, spin="ALPHA", last=.TRUE.)
            CALL write_mo_set(mos(2)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                              dft_section, spin="BETA", last=.TRUE.)
         ELSE
            CALL write_mo_set(mos(1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                              dft_section, last=.TRUE.)
         END IF
         !
         SELECT CASE (tb_type)
         CASE ("DFTB")
            !
         CASE ("xTB")
            sprint_section => section_vals_get_subs_vals(dft_section, "PRINT%MO_MOLDEN")
            CALL write_mos_molden(mos, qs_kind_set, particle_set, sprint_section)
         CASE DEFAULT
            CPABORT("unknown TB type")
         END SELECT
      ELSE
         CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
         IF (kpoints%nkp /= 0) THEN
            DO ik = 1, SIZE(kpoints%kp_env)
               mos_k => kpoints%kp_env(ik)%kpoint_env%mos
               CPASSERT(ASSOCIATED(mos_k))
               IF (dft_control%nspins == 2) THEN
                  CALL write_mo_set(mos_k(1, 1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                    dft_section, spin="ALPHA", last=.TRUE., kpt=ik)
                  CALL write_mo_set(mos_k(1, 2)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                    dft_section, spin="BETA", last=.TRUE., kpt=ik)
               ELSE
                  CALL write_mo_set(mos_k(1, 1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                    dft_section, last=.TRUE., kpt=ik)
               END IF
            END DO
         END IF
      END IF

   END SUBROUTINE print_mos

! **************************************************************************************************
!> \brief Gets the lumos, and eigenvalues for the lumos
!> \param qs_env ...
!> \param scf_env ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param nlumo ...
!> \param nlumos ...
! **************************************************************************************************
   SUBROUTINE make_lumo_tb(qs_env, scf_env, unoccupied_orbs, unoccupied_evals, nlumo, nlumos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: unoccupied_orbs
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: unoccupied_evals
      INTEGER                                            :: nlumo
      INTEGER, INTENT(OUT)                               :: nlumos

      CHARACTER(len=*), PARAMETER :: routineN = 'make_lumo_tb', routineP = moduleN//':'//routineN

      INTEGER                                            :: homo, iounit, ispin, n, nao, nmo
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(preconditioner_type), POINTER                 :: local_preconditioner
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (mos, ks_rmpv, scf_control, dft_control, para_env, blacs_env)
      CALL get_qs_env(qs_env, &
                      mos=mos, &
                      matrix_ks=ks_rmpv, &
                      scf_control=scf_control, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      DO ispin = 1, dft_control%nspins
         NULLIFY (unoccupied_orbs(ispin)%matrix)
         NULLIFY (unoccupied_evals(ispin)%array)
         ! Always write eigenvalues
         IF (iounit > 0) WRITE (iounit, *) " "
         IF (iounit > 0) WRITE (iounit, *) " Lowest Eigenvalues of the unoccupied subspace spin ", ispin
         IF (iounit > 0) WRITE (iounit, FMT='(1X,A)') "-----------------------------------------------------"
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, homo=homo, nao=nao, nmo=nmo)
         CALL cp_fm_get_info(mo_coeff, nrow_global=n)
         nlumos = MAX(1, MIN(nlumo, nao - nmo))
         IF (nlumo == -1) nlumos = nao - nmo
         ALLOCATE (unoccupied_evals(ispin)%array(nlumos))
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=n, ncol_global=nlumos)
         CALL cp_fm_create(unoccupied_orbs(ispin)%matrix, fm_struct_tmp, name="lumos")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_init_random(unoccupied_orbs(ispin)%matrix, nlumos)

         ! the full_all preconditioner makes not much sense for lumos search
         NULLIFY (local_preconditioner)
         IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
            local_preconditioner => scf_env%ot_preconditioner(1)%preconditioner
            ! this one can for sure not be right (as it has to match a given C0)
            IF (local_preconditioner%in_use == ot_precond_full_all) THEN
               NULLIFY (local_preconditioner)
            ENDIF
         ENDIF

         CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix, matrix_s=matrix_s(1)%matrix, &
                             matrix_c_fm=unoccupied_orbs(ispin)%matrix, &
                             matrix_orthogonal_space_fm=mo_coeff, &
                             eps_gradient=scf_control%eps_lumos, &
                             preconditioner=local_preconditioner, &
                             iter_max=scf_control%max_iter_lumos, &
                             size_ortho_space=nmo)

         CALL calculate_subspace_eigenvalues(unoccupied_orbs(ispin)%matrix, ks_rmpv(ispin)%matrix, &
                                             unoccupied_evals(ispin)%array, scr=iounit, &
                                             ionode=iounit > 0)

      END DO

   END SUBROUTINE make_lumo_tb

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE rebuild_pw_env(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'rebuild_pw_env', routineP = moduleN//':'//routineN

      LOGICAL                                            :: skip_load_balance_distributed
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: new_pw_env
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(task_list_type), POINTER                      :: task_list

      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control, pw_env=new_pw_env)
      IF (.NOT. ASSOCIATED(new_pw_env)) THEN
         CALL pw_env_create(new_pw_env)
         CALL set_ks_env(ks_env, pw_env=new_pw_env)
         CALL pw_env_release(new_pw_env)
      ENDIF
      CALL get_qs_env(qs_env, pw_env=new_pw_env, dft_control=dft_control, cell=cell)

      new_pw_env%cell_hmat = cell%hmat
      CALL pw_env_rebuild(new_pw_env, qs_env=qs_env)

      NULLIFY (task_list)
      CALL get_ks_env(ks_env, task_list=task_list)
      IF (.NOT. ASSOCIATED(task_list)) THEN
         CALL allocate_task_list(task_list)
         CALL set_ks_env(ks_env, task_list=task_list)
      ENDIF
      skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
      CALL generate_qs_task_list(ks_env, task_list, &
                                 reorder_rs_grid_ranks=.TRUE., soft_valid=.FALSE., &
                                 skip_load_balance_distributed=skip_load_balance_distributed)
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_rebuild(rho, qs_env=qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)

   END SUBROUTINE rebuild_pw_env

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cube_section ...
! **************************************************************************************************
   SUBROUTINE print_e_density(qs_env, cube_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: cube_section

      CHARACTER(len=*), PARAMETER :: routineN = 'print_e_density', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, mpi_filename, my_pos_cube
      INTEGER                                            :: iounit, ispin, unit_nr
      LOGICAL                                            :: append_cube, mpi_io
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_elec_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL get_qs_env(qs_env, dft_control=dft_control)

      append_cube = section_get_lval(cube_section, "APPEND")
      my_pos_cube = "REWIND"
      IF (append_cube) my_pos_cube = "APPEND"

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      ! we need to construct the density on a realspace grid
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho)
      NULLIFY (rho_r, rho_g, tot_rho_r)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp, &
                      rho_r=rho_r, rho_g=rho_g, tot_rho_r=tot_rho_r)
      DO ispin = 1, dft_control%nspins
         rho_ao => rho_ao_kp(ispin, :)
         CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                 rho=rho_r(ispin), &
                                 rho_gspace=rho_g(ispin), &
                                 total_rho=tot_rho_r(ispin), &
                                 ks_env=ks_env)
      END DO
      CALL qs_rho_set(rho, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      IF (dft_control%nspins > 1) THEN
         IF (iounit > 0) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A,T51,2F15.6)") &
               "Integrated alpha and beta electronic density:", tot_rho_r(1:2)
         END IF
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
         CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
         CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                                pw=rho_elec_rspace%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_copy(rho_r(1)%pw, rho_elec_rspace%pw)
         CALL pw_axpy(rho_r(2)%pw, rho_elec_rspace%pw)
         filename = "ELECTRON_DENSITY"
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, cube_section, '', &
                                        extension=".cube", middle_name=TRIM(filename), &
                                        file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                        fout=mpi_filename)
         IF (iounit > 0) THEN
            IF (.NOT. mpi_io) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
            ELSE
               filename = mpi_filename
            END IF
            WRITE (UNIT=iounit, FMT="(T2,A,/,T2,A79)") &
               "The sum of alpha and beta density is written in cube file format to the file:", ADJUSTR(TRIM(filename))
         END IF
         CALL cp_pw_to_cube(rho_elec_rspace%pw, unit_nr, "SUM OF ALPHA AND BETA DENSITY", &
                            particles=particles, stride=section_get_ivals(cube_section, "STRIDE"), &
                            mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, cube_section, '', mpi_io=mpi_io)
         CALL pw_copy(rho_r(1)%pw, rho_elec_rspace%pw)
         CALL pw_axpy(rho_r(2)%pw, rho_elec_rspace%pw, alpha=-1.0_dp)
         filename = "SPIN_DENSITY"
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, cube_section, '', &
                                        extension=".cube", middle_name=TRIM(filename), &
                                        file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                        fout=mpi_filename)
         IF (iounit > 0) THEN
            IF (.NOT. mpi_io) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
            ELSE
               filename = mpi_filename
            END IF
            WRITE (UNIT=iounit, FMT="(T2,A,/,T2,A79)") &
               "The spin density is written in cube file format to the file:", ADJUSTR(TRIM(filename))
         END IF
         CALL cp_pw_to_cube(rho_elec_rspace%pw, unit_nr, "SPIN DENSITY", &
                            particles=particles, &
                            stride=section_get_ivals(cube_section, "STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, cube_section, '', mpi_io=mpi_io)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_elec_rspace%pw)
      ELSE
         IF (iounit > 0) THEN
            WRITE (UNIT=iounit, FMT="(/,T2,A,T66,F15.6)") &
               "Integrated electronic density:", tot_rho_r(1)
         END IF
         filename = "ELECTRON_DENSITY"
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, cube_section, '', &
                                        extension=".cube", middle_name=TRIM(filename), &
                                        file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                        fout=mpi_filename)
         IF (iounit > 0) THEN
            IF (.NOT. mpi_io) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
            ELSE
               filename = mpi_filename
            END IF
            WRITE (UNIT=iounit, FMT="(T2,A,/,T2,A79)") &
               "The electron density is written in cube file format to the file:", ADJUSTR(TRIM(filename))
         END IF
         CALL cp_pw_to_cube(rho_r(1)%pw, unit_nr, "ELECTRON DENSITY", &
                            particles=particles, &
                            stride=section_get_ivals(cube_section, "STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, cube_section, '', mpi_io=mpi_io)
      END IF ! nspins

   END SUBROUTINE print_e_density
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cube_section ...
!> \param total_density ...
!> \param v_hartree ...
!> \param efield ...
! **************************************************************************************************
   SUBROUTINE print_density_cubes(qs_env, cube_section, total_density, v_hartree, efield)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: cube_section
      LOGICAL, INTENT(IN), OPTIONAL                      :: total_density, v_hartree, efield

      CHARACTER(len=*), PARAMETER :: routineN = 'print_density_cubes', &
         routineP = moduleN//':'//routineN
      CHARACTER(len=1), DIMENSION(3), PARAMETER          :: cdir = (/"x", "y", "z"/)

      CHARACTER(LEN=default_path_length)                 :: filename, mpi_filename, my_pos_cube
      INTEGER                                            :: id, iounit, ispin, nd(3), unit_nr
      LOGICAL                                            :: append_cube, mpi_io, my_efield, &
                                                            my_total_density, my_v_hartree
      REAL(KIND=dp)                                      :: total_rho_core_rspace, udvol
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(greens_fn_type), POINTER                      :: green_fft
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: pwdummy
      TYPE(pw_p_type)                                    :: rho_core, rho_tot_gspace, &
                                                            rho_tot_rspace, vhartree
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r
      TYPE(pw_poisson_parameter_type)                    :: poisson_params
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control)

      append_cube = section_get_lval(cube_section, "APPEND")
      my_pos_cube = "REWIND"
      IF (append_cube) my_pos_cube = "APPEND"

      IF (PRESENT(total_density)) THEN
         my_total_density = total_density
      ELSE
         my_total_density = .FALSE.
      END IF
      IF (PRESENT(v_hartree)) THEN
         my_v_hartree = v_hartree
      ELSE
         my_v_hartree = .FALSE.
      END IF
      IF (PRESENT(efield)) THEN
         my_efield = efield
      ELSE
         my_efield = .FALSE.
      END IF

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      ! we need to construct the density on a realspace grid
      CALL get_qs_env(qs_env, ks_env=ks_env, rho=rho)
      NULLIFY (rho_r, rho_g, tot_rho_r)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp, &
                      rho_r=rho_r, rho_g=rho_g, tot_rho_r=tot_rho_r)
      DO ispin = 1, dft_control%nspins
         rho_ao => rho_ao_kp(ispin, :)
         CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                 rho=rho_r(ispin), &
                                 rho_gspace=rho_g(ispin), &
                                 total_rho=tot_rho_r(ispin), &
                                 ks_env=ks_env)
      END DO
      CALL qs_rho_set(rho, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
      CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                             pw=rho_core%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL calculate_rho_core(rho_core, total_rho_core_rspace, qs_env)

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A,T66,F15.6)") &
            "Integrated electronic density:", SUM(tot_rho_r(:))
         WRITE (UNIT=iounit, FMT="(T2,A,T66,F15.6)") &
            "Integrated core density:", total_rho_core_rspace
      END IF

      CALL pw_pool_create_pw(pool=auxbas_pw_pool, pw=rho_tot_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_transfer(rho_core%pw, rho_tot_rspace%pw)
      DO ispin = 1, dft_control%nspins
         CALL pw_axpy(rho_r(ispin)%pw, rho_tot_rspace%pw)
      END DO

      IF (my_total_density) THEN
         filename = "TOTAL_DENSITY"
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, cube_section, '', &
                                        extension=".cube", middle_name=TRIM(filename), file_position=my_pos_cube, &
                                        log_filename=.FALSE., mpi_io=mpi_io, fout=mpi_filename)
         IF (iounit > 0) THEN
            IF (.NOT. mpi_io) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
            ELSE
               filename = mpi_filename
            END IF
            WRITE (UNIT=iounit, FMT="(T2,A,/,T2,A79)") &
               "The total density is written in cube file format to the file:", ADJUSTR(TRIM(filename))
         END IF
         CALL cp_pw_to_cube(rho_tot_rspace%pw, unit_nr, "TOTAL DENSITY", &
                            particles=particles, &
                            stride=section_get_ivals(cube_section, "STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, cube_section, '', mpi_io=mpi_io)
      END IF
      IF (my_v_hartree .OR. my_efield) THEN
         CALL pw_pool_create_pw(pool=auxbas_pw_pool, pw=rho_tot_gspace%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_transfer(rho_tot_rspace%pw, rho_tot_gspace%pw)
         poisson_params%solver = pw_poisson_analytic
         poisson_params%periodic = cell%perd
         poisson_params%ewald_type = do_ewald_none
         NULLIFY (green_fft, pwdummy)
         CALL pw_green_create(green_fft, poisson_params, cell%hmat, auxbas_pw_pool, pwdummy, pwdummy)
         rho_tot_gspace%pw%cc(:) = rho_tot_gspace%pw%cc(:)*green_fft%influence_fn%cc(:)
         CALL pw_green_release(green_fft, auxbas_pw_pool)
         IF (my_v_hartree) THEN
            CALL pw_pool_create_pw(pool=auxbas_pw_pool, pw=vhartree%pw, &
                                   use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_transfer(rho_tot_gspace%pw, vhartree%pw)
            filename = "V_HARTREE"
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, cube_section, '', &
                                           extension=".cube", middle_name=TRIM(filename), file_position=my_pos_cube, &
                                           log_filename=.FALSE., mpi_io=mpi_io, fout=mpi_filename)
            IF (iounit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=unit_nr, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=iounit, FMT="(T2,A,/,T2,A79)") &
                  "The Hartree potential is written in cube file format to the file:", ADJUSTR(TRIM(filename))
            END IF
            CALL cp_pw_to_cube(vhartree%pw, unit_nr, "Hartree Potential", &
                               particles=particles, &
                               stride=section_get_ivals(cube_section, "STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, cube_section, '', mpi_io=mpi_io)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, vhartree%pw)
         END IF
         IF (my_efield) THEN
            CALL pw_pool_create_pw(pool=auxbas_pw_pool, pw=vhartree%pw, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            udvol = 1.0_dp/rho_tot_rspace%pw%pw_grid%dvol
            DO id = 1, 3
               CALL pw_transfer(rho_tot_gspace%pw, vhartree%pw)
               nd = 0
               nd(id) = 1
               CALL pw_derive(vhartree%pw, nd)
               CALL pw_transfer(vhartree%pw, rho_tot_rspace%pw)
               CALL pw_scale(rho_tot_rspace%pw, udvol)

               filename = "EFIELD_"//cdir(id)
               mpi_io = .TRUE.
               unit_nr = cp_print_key_unit_nr(logger, cube_section, '', &
                                              extension=".cube", middle_name=TRIM(filename), file_position=my_pos_cube, &
                                              log_filename=.FALSE., mpi_io=mpi_io, fout=mpi_filename)
               IF (iounit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=unit_nr, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=iounit, FMT="(T2,A,/,T2,A79)") &
                     "The Efield is written in cube file format to the file:", ADJUSTR(TRIM(filename))
               END IF
               CALL cp_pw_to_cube(rho_tot_rspace%pw, unit_nr, "EFIELD "//cdir(id), &
                                  particles=particles, &
                                  stride=section_get_ivals(cube_section, "STRIDE"), mpi_io=mpi_io)
               CALL cp_print_key_finished_output(unit_nr, logger, cube_section, '', mpi_io=mpi_io)
            END DO
            CALL pw_pool_give_back_pw(auxbas_pw_pool, vhartree%pw)
         END IF
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_core%pw)

   END SUBROUTINE print_density_cubes

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param elf_section ...
! **************************************************************************************************
   SUBROUTINE print_elf(qs_env, elf_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: elf_section

      CHARACTER(len=*), PARAMETER :: routineN = 'print_elf', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, mpi_filename, my_pos_cube, &
                                                            title
      INTEGER                                            :: iounit, ispin, unit_nr
      LOGICAL                                            :: append_cube, mpi_io
      REAL(KIND=dp)                                      :: rho_cutoff
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: elf_r, rho_g, rho_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      ! we need to construct the density on a realspace grid
      CALL get_qs_env(qs_env, dft_control=dft_control, ks_env=ks_env, rho=rho)
      NULLIFY (rho_r, rho_g, tot_rho_r)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp, &
                      rho_r=rho_r, rho_g=rho_g, tot_rho_r=tot_rho_r)
      DO ispin = 1, dft_control%nspins
         rho_ao => rho_ao_kp(ispin, :)
         CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                 rho=rho_r(ispin), &
                                 rho_gspace=rho_g(ispin), &
                                 total_rho=tot_rho_r(ispin), &
                                 ks_env=ks_env)
      END DO
      CALL qs_rho_set(rho, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      ALLOCATE (elf_r(dft_control%nspins))
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
      DO ispin = 1, dft_control%nspins
         CALL pw_pool_create_pw(auxbas_pw_pool, elf_r(ispin)%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_zero(elf_r(ispin)%pw)
      END DO

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "ELF is computed on the real space grid -----"
      END IF
      rho_cutoff = section_get_rval(elf_section, "density_cutoff")
      CALL qs_elf_calc(qs_env, elf_r, rho_cutoff)

      ! write ELF into cube file
      append_cube = section_get_lval(elf_section, "APPEND")
      my_pos_cube = "REWIND"
      IF (append_cube) my_pos_cube = "APPEND"
      DO ispin = 1, dft_control%nspins
         WRITE (filename, '(a5,I1.1)') "ELF_S", ispin
         WRITE (title, *) "ELF spin ", ispin
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, elf_section, '', extension=".cube", &
                                        middle_name=TRIM(filename), file_position=my_pos_cube, &
                                        log_filename=.FALSE., mpi_io=mpi_io, fout=mpi_filename)
         IF (iounit > 0) THEN
            IF (.NOT. mpi_io) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
            ELSE
               filename = mpi_filename
            END IF
            WRITE (UNIT=iounit, FMT="(T2,A,/,T2,A79)") &
               "ELF is written in cube file format to the file:", ADJUSTR(TRIM(filename))
         END IF

         CALL cp_pw_to_cube(elf_r(ispin)%pw, unit_nr, title, particles=particles, &
                            stride=section_get_ivals(elf_section, "STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, elf_section, '', mpi_io=mpi_io)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, elf_r(ispin)%pw)
      END DO

      DEALLOCATE (elf_r)

   END SUBROUTINE print_elf
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cube_section ...
! **************************************************************************************************
   SUBROUTINE print_mo_cubes(qs_env, cube_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: cube_section

      CHARACTER(len=*), PARAMETER :: routineN = 'print_mo_cubes', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube, title
      INTEGER                                            :: homo, i, ifirst, ilast, iounit, ir, &
                                                            ispin, ivector, n_rep, nhomo, nlist, &
                                                            nlumo, nmo, shomo, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: list, list_index
      LOGICAL                                            :: append_cube, mpi_io, write_cube
      REAL(KIND=dp)                                      :: homo_lumo(2, 2)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, mo_derivs
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env, mos=mos, matrix_ks=ks_rmpv, scf_control=scf_control)
      CALL get_qs_env(qs_env, dft_control=dft_control, mo_derivs=mo_derivs)
      CALL make_mo_eig(mos, dft_control%nspins, ks_rmpv, scf_control, mo_derivs)
      NULLIFY (mo_eigenvalues)
      homo = 0
      DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, eigenvalues=mo_eigenvalues, homo=shomo)
         homo_lumo(ispin, 1) = mo_eigenvalues(shomo)
         homo = MAX(homo, shomo)
      END DO
      write_cube = section_get_lval(cube_section, "WRITE_CUBE")
      nlumo = section_get_ival(cube_section, "NLUMO")
      nhomo = section_get_ival(cube_section, "NHOMO")
      NULLIFY (list_index)
      CALL section_vals_val_get(cube_section, "HOMO_LIST", n_rep_val=n_rep)
      IF (n_rep > 0) THEN
         nlist = 0
         DO ir = 1, n_rep
            NULLIFY (list)
            CALL section_vals_val_get(cube_section, "HOMO_LIST", i_rep_val=ir, i_vals=list)
            IF (ASSOCIATED(list)) THEN
               CALL reallocate(list_index, 1, nlist + SIZE(list))
               DO i = 1, SIZE(list)
                  list_index(i + nlist) = list(i)
               END DO
               nlist = nlist + SIZE(list)
            END IF
         END DO
         nhomo = MAXVAL(list_index)
      ELSE
         IF (nhomo == -1) nhomo = homo
         nlist = homo - MAX(1, homo - nhomo + 1) + 1
         ALLOCATE (list_index(nlist))
         DO i = 1, nlist
            list_index(i) = MAX(1, homo - nhomo + 1) + i - 1
         END DO
      END IF

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)

      CALL get_qs_env(qs_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)

      append_cube = section_get_lval(cube_section, "APPEND")
      my_pos_cube = "REWIND"
      IF (append_cube) THEN
         my_pos_cube = "APPEND"
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      particle_set=particle_set)

      IF (nhomo >= 0) THEN
         DO ispin = 1, dft_control%nspins
            ! Prints the cube files of OCCUPIED ORBITALS
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                            eigenvalues=mo_eigenvalues, homo=homo, nmo=nmo)
            IF (write_cube) THEN
               DO i = 1, nlist
                  ivector = list_index(i)
                  IF (ivector > homo) CYCLE
                  CALL calculate_wavefunction(mo_coeff, ivector, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                              cell, dft_control, particle_set, pw_env)
                  WRITE (filename, '(a4,I5.5,a1,I1.1)') "WFN_", ivector, "_", ispin
                  mpi_io = .TRUE.
                  unit_nr = cp_print_key_unit_nr(logger, cube_section, '', extension=".cube", &
                                                 middle_name=TRIM(filename), file_position=my_pos_cube, &
                                                 log_filename=.FALSE., mpi_io=mpi_io)
                  WRITE (title, *) "WAVEFUNCTION ", ivector, " spin ", ispin, " i.e. HOMO - ", ivector - homo
                  CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, &
                                     stride=section_get_ivals(cube_section, "STRIDE"), mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(unit_nr, logger, cube_section, '', mpi_io=mpi_io)
               ENDDO
            ENDIF
         END DO
      ENDIF

      IF (nlumo /= 0) THEN
         DO ispin = 1, dft_control%nspins
            ! Prints the cube files of UNOCCUPIED ORBITALS
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                            eigenvalues=mo_eigenvalues, homo=homo, nmo=nmo)
            IF (write_cube) THEN
               ifirst = homo + 1
               IF (nlumo == -1) THEN
                  ilast = nmo
               ELSE
                  ilast = ifirst + nlumo - 1
                  ilast = MIN(nmo, ilast)
               END IF
               DO ivector = ifirst, ilast
                  CALL calculate_wavefunction(mo_coeff, ivector, wf_r, wf_g, atomic_kind_set, &
                                              qs_kind_set, cell, dft_control, particle_set, pw_env)
                  WRITE (filename, '(a4,I5.5,a1,I1.1)') "WFN_", ivector, "_", ispin
                  mpi_io = .TRUE.
                  unit_nr = cp_print_key_unit_nr(logger, cube_section, '', extension=".cube", &
                                                 middle_name=TRIM(filename), file_position=my_pos_cube, &
                                                 log_filename=.FALSE., mpi_io=mpi_io)
                  WRITE (title, *) "WAVEFUNCTION ", ivector, " spin ", ispin, " i.e. LUMO + ", ivector - ifirst
                  CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, &
                                     stride=section_get_ivals(cube_section, "STRIDE"), mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(unit_nr, logger, cube_section, '', mpi_io=mpi_io)
               ENDDO
            END IF
         END DO
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
      IF (ASSOCIATED(list_index)) DEALLOCATE (list_index)

   END SUBROUTINE print_mo_cubes

! **************************************************************************************************

END MODULE qs_scf_post_tb
