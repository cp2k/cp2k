!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE input_constants,                 ONLY: ot_precond_full_all
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: twopi
   USE molden_utils,                    ONLY: write_mos_molden
   USE moments_utils,                   ONLY: get_reference_point
   USE mulliken,                        ONLY: mulliken_charges
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: debye
   USE population_analyses,             ONLY: lowdin_population_analysis
   USE preconditioner_types,            ONLY: preconditioner_type
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
   USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_mo_io,                        ONLY: write_mo_set
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues,&
                                              make_mo_eig
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_pdos,                         ONLY: calculate_projected_dos
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
   USE qs_scf_wfn_mix,                  ONLY: wfn_mix
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE scf_control_types,               ONLY: scf_control_type
   USE stm_images,                      ONLY: th_stm_image
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_tb'
   PUBLIC :: scf_post_calculation_tb

CONTAINS

! **************************************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!> \param qs_env ...
!> \param scf_env ...
!> \param dft_section ...
!> \param tb_type ...
!> \par History
!>      03.2013 copy of scf_post_gpw
!> \author JHU
!> \note
! **************************************************************************************************
   SUBROUTINE scf_post_calculation_tb(qs_env, scf_env, dft_section, tb_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env
      TYPE(section_vals_type), POINTER                   :: dft_section
      CHARACTER(LEN=*)                                   :: tb_type

      CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_tb', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=2)                                   :: ana
      CHARACTER(LEN=default_string_length)               :: aname
      INTEGER                                            :: handle, homo, iat, iatom, ikind, ispin, &
                                                            nat, natom, nkind, nlumo_stm, nlumos, &
                                                            nspins, print_level, unit_nr
      LOGICAL                                            :: do_kpoints, explicit, has_homo, print_it, do_cube, rebuild
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
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: moments_section, print_key, &
                                                            print_section, wfn_mix_section
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(qs_env))
      NULLIFY (dft_control, rho, para_env, matrix_s, matrix_p)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
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
         SELECT CASE (TRIM (tb_type))
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
            mcharge(iat) = zeff-SUM(charges(iat, 1:nspins))
         END DO
      END DO

      print_section => section_vals_get_subs_vals(dft_section, "PRINT")

      ! Mulliken
      print_key => section_vals_get_subs_vals(print_section, "MULLIKEN")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_section, "MULLIKEN", extension=".mulliken", log_filename=.FALSE.)
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(/,/,T2,A)") "MULLIKEN POPULATION ANALYSIS"
            IF (nspins == 1) THEN
               WRITE (UNIT=unit_nr, FMT="(/,T2,A,T70,A)") &
                  " # Atom   Element   Kind        Atomic population", " Net charge"
               DO ikind = 1, nkind
                  CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
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
                            FMT="(T2,I7,6X,A2,3X,I6,T39,F12.6,T69,F12.6)") &
                        iat, ana, ikind, charges(iat, 1), mcharge(iat)
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
                            FMT="(T2,I6,5X,A2,2X,I6,T29,4(1X,F12.6))") &
                        iat, ana, ikind, charges(iat, 1:2), mcharge(iat), charges(iat, 1)-charges(iat, 2)
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
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Hirshfeld charges not available for TB methods.")
      END IF

      ! MAO
      print_key => section_vals_get_subs_vals(print_section, "MAO_ANALYSIS")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("MAO analysis not available for TB methods.")
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
      print_key => section_vals_get_subs_vals(print_section, "MO")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CALL print_mos(qs_env, scf_env, dft_section, tb_type)
      END IF

      ! Wavefunction mixing
      wfn_mix_section => section_vals_get_subs_vals(dft_section, "PRINT%WFN_MIX")
      CALL section_vals_get(wfn_mix_section, explicit=explicit)
      IF (explicit .AND. .NOT. qs_env%run_rtp) CALL wfn_mix_tb(qs_env, dft_section, scf_env)

      ! PDOS
      print_key => section_vals_get_subs_vals(print_section, "PDOS")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         IF (do_kpoints) THEN
            CPWARN("Projected density of states not implemented for k-points.")
         ELSE
            CALL get_qs_env(qs_env, mos=mos, matrix_ks=ks_rmpv)
            DO ispin = 1, dft_control%nspins
               IF (PRESENT(scf_env)) THEN
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

      ! V_XC CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "V_XC_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("XC potential cube file not available for TB methods.")
      END IF

      ! EFIELD CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "EFIELD_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Efield cube file not implemented for TB methods.")
      END IF

      ! DENSITY CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "E_DENSITY_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Electronic density cube file not implemented for TB methods.")
      END IF

      ! TOTAL DENSITY CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "TOT_DENSITY_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Total density cube file not implemented for TB methods.")
      END IF

      ! V_Hartree CUBE FILE
      print_key => section_vals_get_subs_vals(print_section, "V_HARTREE_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Hartree potential cube file not implemented for TB methods.")
      END IF

      ! ELF
      print_key => section_vals_get_subs_vals(print_section, "ELF_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("ELF not implemented for TB methods.")
      END IF

      ! MO CUBES
      print_key => section_vals_get_subs_vals(print_section, "MO_CUBES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Printing of MO cube files not implemented for TB methods.")
      END IF

      ! STM
      print_key => section_vals_get_subs_vals(print_section, "STM")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
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
               CALL make_lumo(qs_env, scf_env, unoccupied_orbs_stm, unoccupied_evals_stm, &
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

      ! these print keys are not supported in TB

      ! Electric field gradients
      print_key => section_vals_get_subs_vals(print_section, "ELECTRIC_FIELD_GRADIENT")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Electric field gradient not implemented for TB methods.")
      END IF

      ! KINETIC ENERGY
      print_key => section_vals_get_subs_vals(print_section, "KINETIC_ENERGY")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Kinetic energy not available for TB methods.")
      END IF

      ! Xray diffraction spectrum
      print_key => section_vals_get_subs_vals(print_section, "XRAY_DIFFRACTION_SPECTRUM")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Xray diffraction spectrum not implemented for TB methods.")
      END IF

      ! EPR Hyperfine Coupling
      print_key => section_vals_get_subs_vals(print_section, "HYPERFINE_COUPLING_TENSOR")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("Hyperfine Coupling not implemented for TB methods.")
      END IF

      ! PLUS_U
      print_key => section_vals_get_subs_vals(print_section, "PLUS_U")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         CPWARN("DFT+U method not implemented for TB methods.")
      END IF

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
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (atomic_kind_set, cell)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, cell=cell)

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
         dzphase = charge_tot*CMPLX(-SIN(ria), COS(ria), dp)**(charge_tot-1.0_dp)*dria

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
                  dzeta = -q*CMPLX(-SIN(theta), COS(theta), KIND=dp)**(-q-1.0_dp)*dtheta
                  dggamma(j) = dggamma(j)*zeta+ggamma(j)*dzeta
                  ggamma(j) = ggamma(j)*zeta
               END DO
            ENDDO
         END DO
         dggamma = dggamma*zphase+ggamma*dzphase
         ggamma = ggamma*zphase
         IF (ALL(REAL(ggamma, KIND=dp) /= 0.0_dp)) THEN
            tmp = AIMAG(ggamma)/REAL(ggamma, KIND=dp)
            ci = ATAN(tmp)
            dci = (1.0_dp/(1.0_dp+tmp**2))* &
                  (AIMAG(dggamma)*REAL(ggamma, KIND=dp)-AIMAG(ggamma)*REAL(dggamma, KIND=dp))/(REAL(ggamma, KIND=dp))**2
            dipole = MATMUL(cell%hmat, ci)/twopi
            dipole_deriv = MATMUL(cell%hmat, dci)/twopi
         END IF
      ELSE
         dipole_type = "[Non Periodic]"
         DO i = 1, SIZE(particle_set)
            ! no pbc(particle_set(i)%r(:),cell) so that the total dipole is the sum of the molecular dipoles
            ria = particle_set(i)%r(:)
            q = charges(i)
            dipole = dipole-q*(ria-rcc)
            dipole_deriv(:) = dipole_deriv(:)-q*(particle_set(i)%v(:)-drcc)
         END DO
      END IF
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
         CALL cp_fm_struct_create(fmstruct=ao_lumo_struct, nrow_global=nao, ncol_global=nao-nmo, &
                                  template_fmstruct=mo_coeff%matrix_struct)

         CALL cp_fm_create(lumos(ispin)%matrix, matrix_struct=ao_lumo_struct)
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, KS_tmp)
         CALL cp_fm_cholesky_reduce(KS_tmp, S_tmp)
         CALL choose_eigv_solver(KS_tmp, work, mo_eigenvalues)
         CALL cp_fm_cholesky_restore(work, nao, S_tmp, MO_tmp, "SOLVE")
         CALL cp_fm_to_fm_submat(MO_tmp, mo_coeff, nao, nmo, 1, 1, 1, 1)
         CALL cp_fm_to_fm_submat(MO_tmp, lumos(ispin)%matrix, nao, nao-nmo, 1, nmo+1, 1, 1)

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
!> \param scf_env ...
!> \param dft_section ...
!> \param tb_type ...
!> \par History
!> \author JHU
!> \note
! **************************************************************************************************
   SUBROUTINE print_mos(qs_env, scf_env, dft_section, tb_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env
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
            sprint_section => section_vals_get_subs_vals(dft_section, "SCF%PRINT")
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
   SUBROUTINE make_lumo(qs_env, scf_env, unoccupied_orbs, unoccupied_evals, nlumo, nlumos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: unoccupied_orbs
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: unoccupied_evals
      INTEGER                                            :: nlumo
      INTEGER, INTENT(OUT)                               :: nlumos

      CHARACTER(len=*), PARAMETER :: routineN = 'make_lumo', routineP = moduleN//':'//routineN

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
         nlumos = MAX(1, MIN(nlumo, nao-nmo))
         IF (nlumo == -1) nlumos = nao-nmo
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

   END SUBROUTINE make_lumo
! **************************************************************************************************

END MODULE qs_scf_post_tb
