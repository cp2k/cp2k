!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief NEGF based quantum transport calculations
! **************************************************************************************************

MODULE negf_methods
   USE bibliography,                    ONLY: Bailey2006,&
                                              Papior2017,&
                                              cite_reference
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_gemm,&
                                              cp_cfm_scale,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_trace
   USE cp_cfm_types,                    ONLY: &
        copy_cfm_info_type, cp_cfm_cleanup_copy_general, cp_cfm_create, &
        cp_cfm_finish_copy_general, cp_cfm_get_info, cp_cfm_get_submatrix, cp_cfm_p_type, &
        cp_cfm_release, cp_cfm_retain, cp_cfm_set_submatrix, cp_cfm_start_copy_general, &
        cp_cfm_to_fm, cp_cfm_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_info, cp_fm_p_type, cp_fm_release, &
        cp_fm_retain, cp_fm_set_all, cp_fm_to_fm, cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_dot,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_p_type,&
                                              force_env_type
   USE global_types,                    ONLY: global_environment_type
   USE input_constants,                 ONLY: negfint_method_cc,&
                                              negfint_method_simpson
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: pi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_sum
   USE negf_control_types,              ONLY: negf_control_create,&
                                              negf_control_release,&
                                              negf_control_type,&
                                              read_negf_control
   USE negf_env_types,                  ONLY: negf_env_create,&
                                              negf_env_release,&
                                              negf_env_type
   USE negf_green_cache,                ONLY: green_functions_cache_expand,&
                                              green_functions_cache_release,&
                                              green_functions_cache_reorder,&
                                              green_functions_cache_type
   USE negf_green_methods,              ONLY: do_sancho,&
                                              negf_contact_broadening_matrix,&
                                              negf_contact_self_energy,&
                                              negf_retarded_green_function,&
                                              sancho_work_matrices_create,&
                                              sancho_work_matrices_release,&
                                              sancho_work_matrices_type
   USE negf_integr_cc,                  ONLY: &
        cc_interval_full, cc_interval_half, cc_shape_arc, cc_shape_linear, &
        ccquad_double_number_of_points, ccquad_init, ccquad_reduce_and_append_zdata, &
        ccquad_refine_integral, ccquad_release, ccquad_type
   USE negf_integr_simpson,             ONLY: simpsonrule_get_next_nodes,&
                                              simpsonrule_init,&
                                              simpsonrule_refine_integral,&
                                              simpsonrule_release,&
                                              simpsonrule_type,&
                                              sr_shape_arc,&
                                              sr_shape_linear
   USE negf_matrix_utils,               ONLY: invert_cell_to_index,&
                                              negf_copy_fm_submat_to_dbcsr,&
                                              negf_copy_sym_dbcsr_to_fm_submat
   USE negf_subgroup_types,             ONLY: negf_sub_env_create,&
                                              negf_sub_env_release,&
                                              negf_subgroup_env_type
   USE physcon,                         ONLY: e_charge,&
                                              evolt,&
                                              kelvin,&
                                              seconds
   USE qs_density_mixing_types,         ONLY: direct_mixing_nr,&
                                              gspace_mixing_nr
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_gspace_mixing,                ONLY: gspace_mixing
   USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
   USE qs_mixing_utils,                 ONLY: mixing_allocate,&
                                              mixing_init
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_methods,                  ONLY: scf_env_density_mixing
   USE qs_subsys_types,                 ONLY: qs_subsys_type
   USE string_utilities,                ONLY: integer_to_string
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_methods'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: do_negf

! **************************************************************************************************
!> \brief Type to accumulate the total number of points used in integration as well as
!>        the final error estimate
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE integration_status_type
      INTEGER                                            :: npoints
      REAL(kind=dp)                                      :: error
   END TYPE integration_status_type

CONTAINS

! **************************************************************************************************
!> \brief Perform NEGF calculation.
!> \param force_env  Force environment
!> \par History
!>    * 01.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE do_negf(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'do_negf', routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: contact_id_str
      INTEGER                                            :: handle, icontact, ispin, log_unit, &
                                                            ncontacts, npoints, nspins, print_unit
      LOGICAL                                            :: should_output
      REAL(kind=dp)                                      :: energy_max, energy_min
      REAL(kind=dp), DIMENSION(2)                        :: current
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(force_env_p_type), DIMENSION(:), POINTER      :: sub_force_env
      TYPE(global_environment_type), POINTER             :: global_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_env_type)                                :: negf_env
      TYPE(negf_subgroup_env_type)                       :: sub_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: negf_contact_section, &
                                                            negf_mixing_section, negf_section, &
                                                            print_section, root_section

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      log_unit = cp_logger_get_default_io_unit()

      CALL cite_reference(Bailey2006)
      CALL cite_reference(Papior2017)

      NULLIFY (blacs_env, cp_subsys, global_env, qs_env, root_section, sub_force_env)
      CALL force_env_get(force_env, globenv=global_env, qs_env=qs_env, root_section=root_section, &
                         sub_force_env=sub_force_env, subsys=cp_subsys)

      CALL get_qs_env(qs_env, blacs_env=blacs_env)

      negf_section => section_vals_get_subs_vals(root_section, "NEGF")
      negf_contact_section => section_vals_get_subs_vals(negf_section, "CONTACT")
      negf_mixing_section => section_vals_get_subs_vals(negf_section, "MIXING")

      NULLIFY (negf_control)
      CALL negf_control_create(negf_control)
      CALL read_negf_control(negf_control, root_section, cp_subsys)

      CALL negf_sub_env_create(sub_env, negf_control, blacs_env, global_env%blacs_grid_layout, global_env%blacs_repeatable)
      CALL negf_env_create(negf_env, sub_env, negf_control, force_env, negf_mixing_section, log_unit)

      ! compute contact Fermi level as well as requested properties
      ncontacts = SIZE(negf_control%contacts)
      DO icontact = 1, ncontacts
         NULLIFY (qs_env)
         IF (negf_control%contacts(icontact)%force_env_index > 0) THEN
            CALL force_env_get(sub_force_env(negf_control%contacts(icontact)%force_env_index)%force_env, qs_env=qs_env)
         ELSE
            CALL force_env_get(force_env, qs_env=qs_env)
         END IF
         CALL guess_fermi_level(icontact, negf_env, negf_control, sub_env, qs_env, log_unit)

         print_section => section_vals_get_subs_vals(negf_contact_section, "PRINT", i_rep_section=icontact)
         should_output = BTEST(cp_print_key_should_output(logger%iter_info, print_section, "DOS"), cp_p_file)

         IF (should_output) THEN
            CALL section_vals_val_get(print_section, "DOS%FROM_ENERGY", r_val=energy_min)
            CALL section_vals_val_get(print_section, "DOS%TILL_ENERGY", r_val=energy_max)
            CALL section_vals_val_get(print_section, "DOS%N_GRIDPOINTS", i_val=npoints)

            CALL integer_to_string(icontact, contact_id_str)
            print_unit = cp_print_key_unit_nr(logger, print_section, "DOS", &
                                              extension=".dos", &
                                              middle_name=TRIM(ADJUSTL(contact_id_str)), &
                                              file_status="REPLACE")

            CALL negf_print_dos(print_unit, energy_min, energy_max, npoints, &
                                v_shift=0.0_dp, negf_env=negf_env, negf_control=negf_control, &
                                sub_env=sub_env, base_contact=icontact, just_contact=icontact)

            CALL cp_print_key_finished_output(print_unit, logger, print_section, "DOS")
         END IF
      END DO

      IF (ncontacts > 1) THEN
         NULLIFY (qs_env)
         CALL force_env_get(force_env, qs_env=qs_env)
         CALL shift_potential(negf_env, negf_control, sub_env, qs_env, base_contact=1, log_unit=log_unit)

         CALL converge_density(negf_env, negf_control, sub_env, qs_env, negf_control%v_shift, base_contact=1, log_unit=log_unit)

         ! current
         CALL get_qs_env(qs_env, dft_control=dft_control)

         nspins = dft_control%nspins

         CPASSERT(nspins <= 2)
         DO ispin = 1, nspins
            ! compute the electric current flown through a pair of electrodes
            ! contact_id1 -> extended molecule -> contact_id2.
            ! Only extended systems with two electrodes are supported at the moment,
            ! so for the time being the contacts' indices are hardcoded.
            current(ispin) = negf_compute_current(contact_id1=1, contact_id2=2, &
                                                  v_shift=negf_control%v_shift, &
                                                  negf_env=negf_env, &
                                                  negf_control=negf_control, &
                                                  sub_env=sub_env, &
                                                  ispin=ispin, &
                                                  blacs_env_global=blacs_env)
         END DO

         IF (log_unit > 0) THEN
            IF (nspins > 1) THEN
               WRITE (log_unit, '(/,T2,A,T60,ES20.7E2)') "NEGF| Alpha-spin electric current (A)", current(1)
               WRITE (log_unit, '(T2,A,T60,ES20.7E2)') "NEGF|  Beta-spin electric current (A)", current(2)
            ELSE
               WRITE (log_unit, '(/,T2,A,T60,ES20.7E2)') "NEGF|  Electric current (A)", 2.0_dp*current(1)
            END IF
         END IF

         ! density of states
         print_section => section_vals_get_subs_vals(negf_section, "PRINT")
         should_output = BTEST(cp_print_key_should_output(logger%iter_info, print_section, "DOS"), cp_p_file)

         IF (should_output) THEN
            CALL section_vals_val_get(print_section, "DOS%FROM_ENERGY", r_val=energy_min)
            CALL section_vals_val_get(print_section, "DOS%TILL_ENERGY", r_val=energy_max)
            CALL section_vals_val_get(print_section, "DOS%N_GRIDPOINTS", i_val=npoints)

            CALL integer_to_string(0, contact_id_str)
            print_unit = cp_print_key_unit_nr(logger, print_section, "DOS", &
                                              extension=".dos", &
                                              middle_name=TRIM(ADJUSTL(contact_id_str)), &
                                              file_status="REPLACE")

            CALL negf_print_dos(print_unit, energy_min, energy_max, npoints, negf_control%v_shift, &
                                negf_env=negf_env, negf_control=negf_control, &
                                sub_env=sub_env, base_contact=1)

            CALL cp_print_key_finished_output(print_unit, logger, print_section, "DOS")
         END IF

         ! transmission coefficient
         should_output = BTEST(cp_print_key_should_output(logger%iter_info, print_section, "TRANSMISSION"), cp_p_file)

         IF (should_output) THEN
            CALL section_vals_val_get(print_section, "TRANSMISSION%FROM_ENERGY", r_val=energy_min)
            CALL section_vals_val_get(print_section, "TRANSMISSION%TILL_ENERGY", r_val=energy_max)
            CALL section_vals_val_get(print_section, "TRANSMISSION%N_GRIDPOINTS", i_val=npoints)

            CALL integer_to_string(0, contact_id_str)
            print_unit = cp_print_key_unit_nr(logger, print_section, "TRANSMISSION", &
                                              extension=".transm", &
                                              middle_name=TRIM(ADJUSTL(contact_id_str)), &
                                              file_status="REPLACE")

            CALL negf_print_transmission(print_unit, energy_min, energy_max, npoints, negf_control%v_shift, &
                                         negf_env=negf_env, negf_control=negf_control, &
                                         sub_env=sub_env, contact_id1=1, contact_id2=2)

            CALL cp_print_key_finished_output(print_unit, logger, print_section, "TRANSMISSION")
         END IF

      END IF

      CALL negf_env_release(negf_env)
      CALL negf_sub_env_release(sub_env)

      CALL negf_control_release(negf_control)
      CALL timestop(handle)
   END SUBROUTINE do_negf

! **************************************************************************************************
!> \brief Compute the contact's Fermi level.
!> \param contact_id    index of the contact
!> \param negf_env      NEGF environment
!> \param negf_control  NEGF control
!> \param sub_env       NEGF parallel (sub)group environment
!> \param qs_env        QuickStep environment
!> \param log_unit      output unit
!> \par History
!>    * 10.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE guess_fermi_level(contact_id, negf_env, negf_control, sub_env, qs_env, log_unit)
      INTEGER, INTENT(in)                                :: contact_id
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'guess_fermi_level', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: temperature_str
      COMPLEX(kind=dp)                                   :: lbound_cpath, lbound_lpath, ubound_lpath
      INTEGER                                            :: direction_axis_abs, handle, image, &
                                                            ispin, nao, nimages, nspins, step
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: do_kpoints
      REAL(kind=dp) :: delta_au, energy_ubound_minus_fermi, fermi_level_guess, fermi_level_max, &
         fermi_level_min, nelectrons_guess, nelectrons_max, nelectrons_min, nelectrons_qs_cell0, &
         nelectrons_qs_cell1, offset_au, rscale, t1, t2, trace
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: matrix_s_fm, rho_ao_fm, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env_global
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, rho_ao_qs_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(green_functions_cache_type)                   :: g_surf_cache
      TYPE(integration_status_type)                      :: stats
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL timeset(routineN, handle)

      IF (negf_control%contacts(contact_id)%compute_fermi_level) THEN
         CALL get_qs_env(qs_env, &
                         blacs_env=blacs_env_global, &
                         dft_control=dft_control, &
                         do_kpoints=do_kpoints, &
                         kpoints=kpoints, &
                         matrix_s_kp=matrix_s_kp, &
                         para_env=para_env_global, &
                         rho=rho_struct, subsys=subsys)
         CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_qs_kp)

         nimages = dft_control%nimages
         nspins = dft_control%nspins
         direction_axis_abs = ABS(negf_env%contacts(contact_id)%direction_axis)

         CPASSERT(SIZE(negf_env%contacts(contact_id)%h_00) == nspins)

         IF (sub_env%ngroups > 1) THEN
            NULLIFY (matrix_s_fm, rho_ao_fm, fm_struct)

            CALL cp_fm_get_info(negf_env%contacts(contact_id)%s_00, nrow_global=nao)
            CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env_global)
            CALL cp_fm_create(rho_ao_fm, fm_struct)

            CALL cp_fm_create(matrix_s_fm, fm_struct)
            CALL cp_fm_struct_release(fm_struct)

            IF (sub_env%group_distribution(sub_env%mepos_global) == 0) THEN
               work_fm => negf_env%contacts(contact_id)%s_00
            ELSE
               NULLIFY (work_fm)
            END IF

            CALL cp_fm_copy_general(work_fm, matrix_s_fm, para_env_global)
         ELSE
            matrix_s_fm => negf_env%contacts(contact_id)%s_00
            CALL cp_fm_retain(matrix_s_fm)
            CALL cp_fm_get_info(matrix_s_fm, matrix_struct=fm_struct)
            CALL cp_fm_create(rho_ao_fm, fm_struct)
         END IF

         IF (do_kpoints) THEN
            CALL get_kpoint_info(kpoints, cell_to_index=cell_to_index)
         ELSE
            ALLOCATE (cell_to_index(0:0, 0:0, 0:0))
            cell_to_index(0, 0, 0) = 1
         END IF

         ALLOCATE (index_to_cell(3, nimages))
         CALL invert_cell_to_index(cell_to_index, nimages, index_to_cell)
         IF (.NOT. do_kpoints) DEALLOCATE (cell_to_index)

         IF (nspins == 1) THEN
            ! spin-restricted calculation: number of electrons must be doubled
            rscale = 2.0_dp
         ELSE
            rscale = 1.0_dp
         END IF

         ! compute the refence number of electrons using the electron density
         nelectrons_qs_cell0 = 0.0_dp
         nelectrons_qs_cell1 = 0.0_dp
         DO image = 1, nimages
            IF (index_to_cell(direction_axis_abs, image) == 0) THEN
               DO ispin = 1, nspins
                  CALL dbcsr_dot(rho_ao_qs_kp(ispin, image)%matrix, matrix_s_kp(1, image)%matrix, trace)
                  nelectrons_qs_cell0 = nelectrons_qs_cell0 + trace
               END DO
            ELSE IF (ABS(index_to_cell(direction_axis_abs, image)) == 1) THEN
               DO ispin = 1, nspins
                  CALL dbcsr_dot(rho_ao_qs_kp(ispin, image)%matrix, matrix_s_kp(1, image)%matrix, trace)
                  nelectrons_qs_cell1 = nelectrons_qs_cell1 + trace
               END DO
            END IF
         END DO
         DEALLOCATE (index_to_cell)

         IF (log_unit > 0) THEN
            WRITE (temperature_str, '(F11.3)') negf_control%contacts(contact_id)%temperature*kelvin
            WRITE (log_unit, '(/,T2,A,I0,A)') "COMPUTE FERMI LEVEL OF CONTACT ", &
               contact_id, " AT "//TRIM(ADJUSTL(temperature_str))//" KELVIN"
            WRITE (log_unit, '(/,T2,A,T60,F20.10,/)') "Electronic density of the isolated contact unit cell:", &
               -1.0_dp*(nelectrons_qs_cell0 + nelectrons_qs_cell1)
            WRITE (log_unit, '(T3,A)') "Step     Integration method      Time      Fermi level   Convergence (density)"
            WRITE (log_unit, '(T3,78("-"))')
         END IF

         nelectrons_qs_cell0 = 0.0_dp
         DO ispin = 1, nspins
            CALL cp_fm_trace(negf_env%contacts(contact_id)%rho_00(ispin)%matrix, &
                             negf_env%contacts(contact_id)%s_00, trace)
            nelectrons_qs_cell0 = nelectrons_qs_cell0 + trace
         END DO

         ! Use orbital energies of HOMO and LUMO as reference points and then
         ! refine the Fermi level by using a simple linear interpolation technique
         IF (negf_control%homo_lumo_gap > 0.0_dp) THEN
            IF (negf_control%contacts(contact_id)%refine_fermi_level) THEN
               fermi_level_min = negf_control%contacts(contact_id)%fermi_level
            ELSE
               fermi_level_min = negf_env%contacts(contact_id)%homo_energy
            END IF
            fermi_level_max = fermi_level_min + negf_control%homo_lumo_gap
         ELSE
            IF (negf_control%contacts(contact_id)%refine_fermi_level) THEN
               fermi_level_max = negf_control%contacts(contact_id)%fermi_level
            ELSE
               fermi_level_max = negf_env%contacts(contact_id)%homo_energy
            END IF
            fermi_level_min = fermi_level_max + negf_control%homo_lumo_gap
         END IF

         step = 0
         lbound_cpath = CMPLX(negf_control%energy_lbound, negf_control%eta, kind=dp)
         delta_au = REAL(negf_control%delta_npoles, kind=dp)*twopi*negf_control%contacts(contact_id)%temperature
         offset_au = REAL(negf_control%gamma_kT, kind=dp)*negf_control%contacts(contact_id)%temperature
         energy_ubound_minus_fermi = -2.0_dp*LOG(negf_control%conv_density)*negf_control%contacts(contact_id)%temperature
         t1 = m_walltime()

         DO
            step = step + 1

            SELECT CASE (step)
            CASE (1)
               fermi_level_guess = fermi_level_min
            CASE (2)
               fermi_level_guess = fermi_level_max
            CASE DEFAULT
               fermi_level_guess = fermi_level_min - (nelectrons_min - nelectrons_qs_cell0)* &
                                   (fermi_level_max - fermi_level_min)/(nelectrons_max - nelectrons_min)
            END SELECT

            negf_control%contacts(contact_id)%fermi_level = fermi_level_guess
            nelectrons_guess = 0.0_dp

            lbound_lpath = CMPLX(fermi_level_guess - offset_au, delta_au, kind=dp)
            ubound_lpath = CMPLX(fermi_level_guess + energy_ubound_minus_fermi, delta_au, kind=dp)

            CALL integration_status_reset(stats)

            DO ispin = 1, nspins
               CALL negf_init_rho_equiv_residuals(rho_ao_fm=rho_ao_fm, &
                                                  v_shift=0.0_dp, &
                                                  ignore_bias=.TRUE., &
                                                  negf_env=negf_env, &
                                                  negf_control=negf_control, &
                                                  sub_env=sub_env, &
                                                  ispin=ispin, &
                                                  base_contact=contact_id, &
                                                  just_contact=contact_id)

               CALL negf_add_rho_equiv_low(rho_ao_fm=rho_ao_fm, &
                                           stats=stats, &
                                           v_shift=0.0_dp, &
                                           ignore_bias=.TRUE., &
                                           negf_env=negf_env, &
                                           negf_control=negf_control, &
                                           sub_env=sub_env, &
                                           ispin=ispin, &
                                           base_contact=contact_id, &
                                           integr_lbound=lbound_cpath, &
                                           integr_ubound=lbound_lpath, &
                                           matrix_s_global=matrix_s_fm, &
                                           is_circular=.TRUE., &
                                           g_surf_cache=g_surf_cache, &
                                           just_contact=contact_id)
               CALL green_functions_cache_release(g_surf_cache)

               CALL negf_add_rho_equiv_low(rho_ao_fm=rho_ao_fm, &
                                           stats=stats, &
                                           v_shift=0.0_dp, &
                                           ignore_bias=.TRUE., &
                                           negf_env=negf_env, &
                                           negf_control=negf_control, &
                                           sub_env=sub_env, &
                                           ispin=ispin, &
                                           base_contact=contact_id, &
                                           integr_lbound=lbound_lpath, &
                                           integr_ubound=ubound_lpath, &
                                           matrix_s_global=matrix_s_fm, &
                                           is_circular=.FALSE., &
                                           g_surf_cache=g_surf_cache, &
                                           just_contact=contact_id)
               CALL green_functions_cache_release(g_surf_cache)

               CALL cp_fm_trace(rho_ao_fm, matrix_s_fm, trace)
               nelectrons_guess = nelectrons_guess + trace
            END DO
            nelectrons_guess = nelectrons_guess*rscale

            t2 = m_walltime()

            IF (log_unit > 0) THEN
               WRITE (log_unit, '(T2,I5,T12,A,T32,F8.1,T42,F15.8,T60,ES20.5E2)') &
                  step, get_method_description_string(stats, negf_control%integr_method), &
                  t2 - t1, fermi_level_guess, nelectrons_guess - nelectrons_qs_cell0
            END IF

            IF (ABS(nelectrons_qs_cell0 - nelectrons_guess) < negf_control%conv_density) EXIT

            SELECT CASE (step)
            CASE (1)
               nelectrons_min = nelectrons_guess
            CASE (2)
               nelectrons_max = nelectrons_guess
            CASE DEFAULT
               IF (fermi_level_guess < fermi_level_min) THEN
                  fermi_level_max = fermi_level_min
                  nelectrons_max = nelectrons_min
                  fermi_level_min = fermi_level_guess
                  nelectrons_min = nelectrons_guess
               ELSE IF (fermi_level_guess > fermi_level_max) THEN
                  fermi_level_min = fermi_level_max
                  nelectrons_min = nelectrons_max
                  fermi_level_max = fermi_level_guess
                  nelectrons_max = nelectrons_guess
               ELSE IF (fermi_level_max - fermi_level_guess < fermi_level_guess - fermi_level_min) THEN
                  fermi_level_max = fermi_level_guess
                  nelectrons_max = nelectrons_guess
               ELSE
                  fermi_level_min = fermi_level_guess
                  nelectrons_min = nelectrons_guess
               END IF
            END SELECT

            t1 = t2
         END DO

         negf_control%contacts(contact_id)%fermi_level = fermi_level_guess

         CALL cp_fm_release(matrix_s_fm)
         CALL cp_fm_release(rho_ao_fm)
      END IF

      IF (log_unit > 0) THEN
         WRITE (temperature_str, '(F11.3)') negf_control%contacts(contact_id)%temperature*kelvin
         WRITE (log_unit, '(/,T2,A,I0)') "NEGF| Contact No. ", contact_id
         WRITE (log_unit, '(T2,A,T62,F18.8)') "NEGF|    Fermi level at "//TRIM(ADJUSTL(temperature_str))// &
            " Kelvin (a.u.):", negf_control%contacts(contact_id)%fermi_level
         WRITE (log_unit, '(T2,A,T62,F18.8)') "NEGF|    Electric potential (V):", &
            negf_control%contacts(contact_id)%v_external*evolt
      END IF

      CALL timestop(handle)
   END SUBROUTINE guess_fermi_level

! **************************************************************************************************
!> \brief Compute shift in Hartree potential
!> \param negf_env      NEGF environment
!> \param negf_control  NEGF control
!> \param sub_env       NEGF parallel (sub)group environment
!> \param qs_env        QuickStep environment
!> \param base_contact  index of the reference contact
!> \param log_unit      output unit
! **************************************************************************************************
   SUBROUTINE shift_potential(negf_env, negf_control, sub_env, qs_env, base_contact, log_unit)
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(in)                                :: base_contact, log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'shift_potential', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp)                                   :: lbound_cpath, ubound_cpath, ubound_lpath
      INTEGER                                            :: handle, ispin, iter_count, nao, &
                                                            ncontacts, nspins
      LOGICAL                                            :: do_kpoints
      REAL(kind=dp) :: mu_base, nelectrons_guess, nelectrons_max, nelectrons_min, nelectrons_ref, &
         t1, t2, temperature, trace, v_shift_guess, v_shift_max, v_shift_min
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: rho_ao_fm
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: matrix_s_fm, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_qs_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(green_functions_cache_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: g_surf_circular, g_surf_linear
      TYPE(integration_status_type)                      :: stats
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(qs_subsys_type), POINTER                      :: subsys

      ncontacts = SIZE(negf_control%contacts)
      ! nothing to do
      IF (.NOT. (ALLOCATED(negf_env%h_s) .AND. ALLOCATED(negf_env%h_sc) .AND. &
                 ASSOCIATED(negf_env%s_s) .AND. ALLOCATED(negf_env%s_sc))) RETURN
      IF (ncontacts < 2) RETURN

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, blacs_env=blacs_env, do_kpoints=do_kpoints, dft_control=dft_control, &
                      para_env=para_env, rho=rho_struct, subsys=subsys)
      CPASSERT(.NOT. do_kpoints)

      ! apply external NEGF potential
      t1 = m_walltime()

      ! need a globally distributed overlap matrix in order to compute integration errors
      IF (sub_env%ngroups > 1) THEN
         NULLIFY (matrix_s_fm, fm_struct)

         CALL cp_fm_get_info(negf_env%s_s, nrow_global=nao)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)

         CALL cp_fm_create(matrix_s_fm, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         IF (sub_env%group_distribution(sub_env%mepos_global) == 0) THEN
            work_fm => negf_env%s_s
         ELSE
            NULLIFY (work_fm)
         END IF

         CALL cp_fm_copy_general(work_fm, matrix_s_fm, para_env)
      ELSE
         matrix_s_fm => negf_env%s_s
         CALL cp_fm_retain(matrix_s_fm)
         CALL cp_fm_get_info(matrix_s_fm, matrix_struct=fm_struct)
      END IF

      CALL cp_fm_get_info(matrix_s_fm, matrix_struct=fm_struct)

      nspins = SIZE(negf_env%h_s)

      mu_base = negf_control%contacts(base_contact)%fermi_level

      ! keep the initial charge density matrix and Kohn-Sham matrix
      CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_qs_kp)

      ! extract the reference density matrix blocks
      nelectrons_ref = 0.0_dp
      ALLOCATE (rho_ao_fm(nspins))
      DO ispin = 1, nspins
         NULLIFY (rho_ao_fm(ispin)%matrix)
         CALL cp_fm_create(rho_ao_fm(ispin)%matrix, fm_struct)

         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=rho_ao_qs_kp(ispin, 1)%matrix, &
                                               fm=rho_ao_fm(ispin)%matrix, &
                                               atomlist_row=negf_control%atomlist_S_screening, &
                                               atomlist_col=negf_control%atomlist_S_screening, &
                                               subsys=subsys, mpi_comm_global=para_env%group, &
                                               do_upper_diag=.TRUE., do_lower=.TRUE.)

         CALL cp_fm_trace(rho_ao_fm(ispin)%matrix, matrix_s_fm, trace)
         nelectrons_ref = nelectrons_ref + trace
      END DO

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A)') "COMPUTE SHIFT IN HARTREE POTENTIAL"
         WRITE (log_unit, '(/,T2,A,T55,F25.14,/)') "Initial electronic density of the scattering region:", -1.0_dp*nelectrons_ref
         WRITE (log_unit, '(T3,A)') "Step     Integration method      Time        V shift     Convergence (density)"
         WRITE (log_unit, '(T3,78("-"))')
      END IF

      temperature = negf_control%contacts(base_contact)%temperature

      ! integration limits: C-path (arch)
      lbound_cpath = CMPLX(negf_control%energy_lbound, negf_control%eta, kind=dp)
      ubound_cpath = CMPLX(mu_base - REAL(negf_control%gamma_kT, kind=dp)*temperature, &
                           REAL(negf_control%delta_npoles, kind=dp)*twopi*temperature, kind=dp)

      ! integration limits: L-path (linear)
      ubound_lpath = CMPLX(mu_base - LOG(negf_control%conv_density)*temperature, &
                           REAL(negf_control%delta_npoles, kind=dp)*twopi*temperature, kind=dp)

      v_shift_min = negf_control%v_shift
      v_shift_max = negf_control%v_shift + negf_control%v_shift_offset

      ALLOCATE (g_surf_circular(nspins), g_surf_linear(nspins))

      DO iter_count = 1, negf_control%v_shift_maxiters
         SELECT CASE (iter_count)
         CASE (1)
            v_shift_guess = v_shift_min
         CASE (2)
            v_shift_guess = v_shift_max
         CASE DEFAULT
            v_shift_guess = v_shift_min - (nelectrons_min - nelectrons_ref)* &
                            (v_shift_max - v_shift_min)/(nelectrons_max - nelectrons_min)
         END SELECT

         ! compute an updated density matrix
         CALL integration_status_reset(stats)

         DO ispin = 1, nspins
            ! closed contour: residuals
            CALL negf_init_rho_equiv_residuals(rho_ao_fm=rho_ao_fm(ispin)%matrix, &
                                               v_shift=v_shift_guess, &
                                               ignore_bias=.TRUE., &
                                               negf_env=negf_env, &
                                               negf_control=negf_control, &
                                               sub_env=sub_env, &
                                               ispin=ispin, &
                                               base_contact=base_contact)

            ! closed contour: C-path
            CALL negf_add_rho_equiv_low(rho_ao_fm=rho_ao_fm(ispin)%matrix, &
                                        stats=stats, &
                                        v_shift=v_shift_guess, &
                                        ignore_bias=.TRUE., &
                                        negf_env=negf_env, &
                                        negf_control=negf_control, &
                                        sub_env=sub_env, &
                                        ispin=ispin, &
                                        base_contact=base_contact, &
                                        integr_lbound=lbound_cpath, &
                                        integr_ubound=ubound_cpath, &
                                        matrix_s_global=matrix_s_fm, &
                                        is_circular=.TRUE., &
                                        g_surf_cache=g_surf_circular(ispin))
            IF (negf_control%disable_cache) &
               CALL green_functions_cache_release(g_surf_circular(ispin))

            ! closed contour: L-path
            CALL negf_add_rho_equiv_low(rho_ao_fm=rho_ao_fm(ispin)%matrix, &
                                        stats=stats, &
                                        v_shift=v_shift_guess, &
                                        ignore_bias=.TRUE., &
                                        negf_env=negf_env, &
                                        negf_control=negf_control, &
                                        sub_env=sub_env, &
                                        ispin=ispin, &
                                        base_contact=base_contact, &
                                        integr_lbound=ubound_cpath, &
                                        integr_ubound=ubound_lpath, &
                                        matrix_s_global=matrix_s_fm, &
                                        is_circular=.FALSE., &
                                        g_surf_cache=g_surf_linear(ispin))
            IF (negf_control%disable_cache) &
               CALL green_functions_cache_release(g_surf_linear(ispin))
         END DO

         IF (nspins > 1) THEN
            DO ispin = 2, nspins
               CALL cp_fm_scale_and_add(1.0_dp, rho_ao_fm(1)%matrix, 1.0_dp, rho_ao_fm(ispin)%matrix)
            END DO
         ELSE
            CALL cp_fm_scale(2.0_dp, rho_ao_fm(1)%matrix)
         END IF

         CALL cp_fm_trace(rho_ao_fm(1)%matrix, matrix_s_fm, nelectrons_guess)

         t2 = m_walltime()

         IF (log_unit > 0) THEN
            WRITE (log_unit, '(T2,I5,T12,A,T32,F8.1,T42,F15.8,T60,ES20.5E2)') &
               iter_count, get_method_description_string(stats, negf_control%integr_method), &
               t2 - t1, v_shift_guess, nelectrons_guess - nelectrons_ref
         END IF

         IF (ABS(nelectrons_guess - nelectrons_ref) < negf_control%conv_scf) EXIT

         ! compute correction
         SELECT CASE (iter_count)
         CASE (1)
            nelectrons_min = nelectrons_guess
         CASE (2)
            nelectrons_max = nelectrons_guess
         CASE DEFAULT
            IF (v_shift_guess < v_shift_min) THEN
               v_shift_max = v_shift_min
               nelectrons_max = nelectrons_min
               v_shift_min = v_shift_guess
               nelectrons_min = nelectrons_guess
            ELSE IF (v_shift_guess > v_shift_max) THEN
               v_shift_min = v_shift_max
               nelectrons_min = nelectrons_max
               v_shift_max = v_shift_guess
               nelectrons_max = nelectrons_guess
            ELSE IF (v_shift_max - v_shift_guess < v_shift_guess - v_shift_min) THEN
               v_shift_max = v_shift_guess
               nelectrons_max = nelectrons_guess
            ELSE
               v_shift_min = v_shift_guess
               nelectrons_min = nelectrons_guess
            END IF
         END SELECT

         t1 = t2
      END DO

      negf_control%v_shift = v_shift_guess

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(T2,A,T62,F18.8)') "NEGF|    Shift in Hartree potential", negf_control%v_shift
      END IF

      DO ispin = nspins, 1, -1
         CALL green_functions_cache_release(g_surf_circular(ispin))
         CALL green_functions_cache_release(g_surf_linear(ispin))
      END DO
      DEALLOCATE (g_surf_circular, g_surf_linear)

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(rho_ao_fm(ispin)%matrix)
      END DO
      DEALLOCATE (rho_ao_fm)

      IF (ASSOCIATED(matrix_s_fm)) &
         CALL cp_fm_release(matrix_s_fm)

      CALL timestop(handle)
   END SUBROUTINE shift_potential

! **************************************************************************************************
!> \brief Converge electronic density of the scattering region.
!> \param negf_env      NEGF environment
!> \param negf_control  NEGF control
!> \param sub_env       NEGF parallel (sub)group environment
!> \param qs_env        QuickStep environment
!> \param v_shift       shift in Hartree potential
!> \param base_contact  index of the reference contact
!> \param log_unit      output unit
!> \par History
!>    * 06.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE converge_density(negf_env, negf_control, sub_env, qs_env, v_shift, base_contact, log_unit)
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(kind=dp), INTENT(in)                          :: v_shift
      INTEGER, INTENT(in)                                :: base_contact, log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'converge_density', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER :: threshold = 16.0_dp*EPSILON(0.0_dp)

      COMPLEX(kind=dp)                                   :: lbound_cpath, ubound_cpath, ubound_lpath
      INTEGER                                            :: handle, icontact, image, ispin, &
                                                            iter_count, nao, ncontacts, nimages, &
                                                            nspins
      LOGICAL                                            :: do_kpoints
      REAL(kind=dp)                                      :: iter_delta, mu_base, nelectrons, &
                                                            nelectrons_diff, t1, t2, temperature, &
                                                            trace, v_base, v_contact
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:)      :: rho_ao_delta_fm, rho_ao_new_fm
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: ao_ao_fm_global, matrix_s_fm, work_fm
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_initial_kp, matrix_ks_qs_kp, &
                                                            rho_ao_initial_kp, rho_ao_new_kp, &
                                                            rho_ao_qs_kp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(green_functions_cache_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: g_surf_circular, g_surf_linear, &
                                                            g_surf_nonequiv
      TYPE(integration_status_type)                      :: stats
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(qs_subsys_type), POINTER                      :: subsys

      ncontacts = SIZE(negf_control%contacts)
      ! nothing to do
      IF (.NOT. (ALLOCATED(negf_env%h_s) .AND. ALLOCATED(negf_env%h_sc) .AND. &
                 ASSOCIATED(negf_env%s_s) .AND. ALLOCATED(negf_env%s_sc))) RETURN
      IF (ncontacts < 2) RETURN

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, blacs_env=blacs_env, do_kpoints=do_kpoints, dft_control=dft_control, &
                      matrix_ks_kp=matrix_ks_qs_kp, para_env=para_env, rho=rho_struct, subsys=subsys)
      CPASSERT(.NOT. do_kpoints)

      ! apply external NEGF potential
      t1 = m_walltime()

      ! need a globally distributed overlap matrix in order to compute integration errors
      IF (sub_env%ngroups > 1) THEN
         NULLIFY (matrix_s_fm, fm_struct)

         CALL cp_fm_get_info(negf_env%s_s, nrow_global=nao)
         CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, context=blacs_env)

         CALL cp_fm_create(matrix_s_fm, fm_struct)
         CALL cp_fm_struct_release(fm_struct)

         IF (sub_env%group_distribution(sub_env%mepos_global) == 0) THEN
            work_fm => negf_env%s_s
         ELSE
            NULLIFY (work_fm)
         END IF

         CALL cp_fm_copy_general(work_fm, matrix_s_fm, para_env)
      ELSE
         matrix_s_fm => negf_env%s_s
         CALL cp_fm_retain(matrix_s_fm)
         CALL cp_fm_get_info(matrix_s_fm, matrix_struct=fm_struct)
      END IF

      CALL cp_fm_get_info(matrix_s_fm, matrix_struct=fm_struct)

      nspins = SIZE(negf_env%h_s)
      nimages = dft_control%nimages

      v_base = negf_control%contacts(base_contact)%v_external
      mu_base = negf_control%contacts(base_contact)%fermi_level + v_base

      ! the current subroutine works for the general case as well, but the Poisson solver does not
      IF (ncontacts > 2) THEN
         CPABORT("Poisson solver does not support the general NEGF setup (>2 contacts).")
      END IF

      ! keep the initial charge density matrix and Kohn-Sham matrix
      CALL qs_rho_get(rho_struct, rho_ao_kp=rho_ao_qs_kp)

      NULLIFY (matrix_ks_initial_kp, rho_ao_initial_kp, rho_ao_new_kp)
      CALL dbcsr_allocate_matrix_set(matrix_ks_initial_kp, nspins, nimages)
      CALL dbcsr_allocate_matrix_set(rho_ao_initial_kp, nspins, nimages)
      CALL dbcsr_allocate_matrix_set(rho_ao_new_kp, nspins, nimages)

      DO image = 1, nimages
         DO ispin = 1, nspins
            CALL dbcsr_init_p(matrix_ks_initial_kp(ispin, image)%matrix)
            CALL dbcsr_copy(matrix_b=matrix_ks_initial_kp(ispin, image)%matrix, matrix_a=matrix_ks_qs_kp(ispin, image)%matrix)

            CALL dbcsr_init_p(rho_ao_initial_kp(ispin, image)%matrix)
            CALL dbcsr_copy(matrix_b=rho_ao_initial_kp(ispin, image)%matrix, matrix_a=rho_ao_qs_kp(ispin, image)%matrix)

            CALL dbcsr_init_p(rho_ao_new_kp(ispin, image)%matrix)
            CALL dbcsr_copy(matrix_b=rho_ao_new_kp(ispin, image)%matrix, matrix_a=rho_ao_qs_kp(ispin, image)%matrix)
         END DO
      END DO

      ! extract the reference density matrix blocks
      nelectrons = 0.0_dp
      ALLOCATE (rho_ao_delta_fm(nspins), rho_ao_new_fm(nspins))
      DO ispin = 1, nspins
         NULLIFY (rho_ao_delta_fm(ispin)%matrix, rho_ao_new_fm(ispin)%matrix)
         CALL cp_fm_create(rho_ao_delta_fm(ispin)%matrix, fm_struct)
         CALL cp_fm_create(rho_ao_new_fm(ispin)%matrix, fm_struct)

         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=rho_ao_qs_kp(ispin, 1)%matrix, &
                                               fm=rho_ao_delta_fm(ispin)%matrix, &
                                               atomlist_row=negf_control%atomlist_S_screening, &
                                               atomlist_col=negf_control%atomlist_S_screening, &
                                               subsys=subsys, mpi_comm_global=para_env%group, &
                                               do_upper_diag=.TRUE., do_lower=.TRUE.)

         CALL cp_fm_trace(rho_ao_delta_fm(ispin)%matrix, matrix_s_fm, trace)
         nelectrons = nelectrons + trace
      END DO

      NULLIFY (ao_ao_fm_global)
      IF (sub_env%ngroups > 1) &
         CALL cp_fm_create(ao_ao_fm_global, fm_struct)

      ! mixing storage allocation
      IF (negf_env%mixing_method >= gspace_mixing_nr) THEN
         CALL mixing_allocate(qs_env, negf_env%mixing_method, nspins=nspins, mixing_store=negf_env%mixing_storage)
         IF (dft_control%qs_control%dftb .OR. dft_control%qs_control%xtb) THEN
            CPABORT('TB Code not available')
         ELSE IF (dft_control%qs_control%semi_empirical) THEN
            CPABORT('SE Code not possible')
         ELSE
            CALL mixing_init(negf_env%mixing_method, rho_struct, negf_env%mixing_storage, para_env)
         END IF
      END IF

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A)') "NEGF SELF-CONSISTENT PROCEDURE"
         WRITE (log_unit, '(/,T2,A,T55,F25.14,/)') "Initial electronic density of the scattering region:", -1.0_dp*nelectrons
         WRITE (log_unit, '(T3,A)') "Step     Integration method      Time     Electronic density      Convergence"
         WRITE (log_unit, '(T3,78("-"))')
      END IF

      temperature = negf_control%contacts(base_contact)%temperature

      DO icontact = 1, ncontacts
         IF (icontact /= base_contact) THEN
            v_contact = negf_control%contacts(icontact)%v_external

            ! integration limits: C-path (arch)
            lbound_cpath = CMPLX(negf_control%energy_lbound, negf_control%eta, kind=dp)
            ubound_cpath = CMPLX(mu_base - REAL(negf_control%gamma_kT, kind=dp)*temperature, &
                                 REAL(negf_control%delta_npoles, kind=dp)*twopi*temperature, kind=dp)

            ! integration limits: L-path (linear)
            ubound_lpath = CMPLX(mu_base - LOG(negf_control%conv_density)*temperature, &
                                 REAL(negf_control%delta_npoles, kind=dp)*twopi*temperature, kind=dp)

            ALLOCATE (g_surf_circular(nspins), g_surf_linear(nspins), g_surf_nonequiv(nspins))

            DO iter_count = 1, negf_control%max_scf
               ! compute an updated density matrix
               CALL integration_status_reset(stats)

               DO ispin = 1, nspins
                  ! closed contour: residuals
                  CALL negf_init_rho_equiv_residuals(rho_ao_fm=rho_ao_new_fm(ispin)%matrix, &
                                                     v_shift=v_shift, &
                                                     ignore_bias=.FALSE., &
                                                     negf_env=negf_env, &
                                                     negf_control=negf_control, &
                                                     sub_env=sub_env, &
                                                     ispin=ispin, &
                                                     base_contact=base_contact)

                  ! closed contour: C-path
                  CALL negf_add_rho_equiv_low(rho_ao_fm=rho_ao_new_fm(ispin)%matrix, &
                                              stats=stats, &
                                              v_shift=v_shift, &
                                              ignore_bias=.FALSE., &
                                              negf_env=negf_env, &
                                              negf_control=negf_control, &
                                              sub_env=sub_env, &
                                              ispin=ispin, &
                                              base_contact=base_contact, &
                                              integr_lbound=lbound_cpath, &
                                              integr_ubound=ubound_cpath, &
                                              matrix_s_global=matrix_s_fm, &
                                              is_circular=.TRUE., &
                                              g_surf_cache=g_surf_circular(ispin))
                  IF (negf_control%disable_cache) &
                     CALL green_functions_cache_release(g_surf_circular(ispin))

                  ! closed contour: L-path
                  CALL negf_add_rho_equiv_low(rho_ao_fm=rho_ao_new_fm(ispin)%matrix, &
                                              stats=stats, &
                                              v_shift=v_shift, &
                                              ignore_bias=.FALSE., &
                                              negf_env=negf_env, &
                                              negf_control=negf_control, &
                                              sub_env=sub_env, &
                                              ispin=ispin, &
                                              base_contact=base_contact, &
                                              integr_lbound=ubound_cpath, &
                                              integr_ubound=ubound_lpath, &
                                              matrix_s_global=matrix_s_fm, &
                                              is_circular=.FALSE., &
                                              g_surf_cache=g_surf_linear(ispin))
                  IF (negf_control%disable_cache) &
                     CALL green_functions_cache_release(g_surf_linear(ispin))

                  ! non-equilibrium part
                  IF (ABS(negf_control%contacts(icontact)%v_external - &
                          negf_control%contacts(base_contact)%v_external) >= threshold) THEN
                     CALL negf_add_rho_nonequiv(rho_ao_fm=rho_ao_new_fm(ispin)%matrix, &
                                                stats=stats, &
                                                v_shift=v_shift, &
                                                negf_env=negf_env, &
                                                negf_control=negf_control, &
                                                sub_env=sub_env, &
                                                ispin=ispin, &
                                                base_contact=base_contact, &
                                                matrix_s_global=matrix_s_fm, &
                                                g_surf_cache=g_surf_nonequiv(ispin))
                     IF (negf_control%disable_cache) &
                        CALL green_functions_cache_release(g_surf_nonequiv(ispin))
                  END IF
               END DO

               IF (nspins == 1) CALL cp_fm_scale(2.0_dp, rho_ao_new_fm(1)%matrix)

               nelectrons = 0.0_dp
               nelectrons_diff = 0.0_dp
               DO ispin = 1, nspins
                  CALL cp_fm_trace(rho_ao_new_fm(ispin)%matrix, matrix_s_fm, trace)
                  nelectrons = nelectrons + trace

                  ! rho_ao_delta_fm contains the original (non-mixed) density matrix from the previous iteration
                  CALL cp_fm_scale_and_add(1.0_dp, rho_ao_delta_fm(ispin)%matrix, -1.0_dp, rho_ao_new_fm(ispin)%matrix)
                  CALL cp_fm_trace(rho_ao_delta_fm(ispin)%matrix, matrix_s_fm, trace)
                  nelectrons_diff = nelectrons_diff + trace

                  ! rho_ao_new_fm -> rho_ao_delta_fm
                  CALL cp_fm_to_fm(rho_ao_new_fm(ispin)%matrix, rho_ao_delta_fm(ispin)%matrix)
               END DO

               t2 = m_walltime()

               IF (log_unit > 0) THEN
                  WRITE (log_unit, '(T2,I5,T12,A,T32,F8.1,T43,F20.8,T65,ES15.5E2)') &
                     iter_count, get_method_description_string(stats, negf_control%integr_method), &
                     t2 - t1, -1.0_dp*nelectrons, nelectrons_diff
               END IF

               IF (ABS(nelectrons_diff) < negf_control%conv_scf) EXIT

               t1 = t2

               ! mix density matrices
               IF (negf_env%mixing_method == direct_mixing_nr) THEN
                  DO image = 1, nimages
                     DO ispin = 1, nspins
                        CALL dbcsr_copy(matrix_b=rho_ao_new_kp(ispin, image)%matrix, &
                                        matrix_a=rho_ao_initial_kp(ispin, image)%matrix)
                     END DO
                  END DO

                  DO ispin = 1, nspins
                     CALL negf_copy_fm_submat_to_dbcsr(fm=rho_ao_new_fm(ispin)%matrix, &
                                                       matrix=rho_ao_new_kp(ispin, 1)%matrix, &
                                                       atomlist_row=negf_control%atomlist_S_screening, &
                                                       atomlist_col=negf_control%atomlist_S_screening, &
                                                       subsys=subsys)
                  END DO

                  CALL scf_env_density_mixing(rho_ao_new_kp, negf_env%mixing_storage, rho_ao_qs_kp, &
                                              para_env, iter_delta, iter_count)

                  DO image = 1, nimages
                     DO ispin = 1, nspins
                        CALL dbcsr_copy(rho_ao_qs_kp(ispin, image)%matrix, rho_ao_new_kp(ispin, image)%matrix)
                     END DO
                  END DO
               ELSE
                  ! store the updated density matrix directly into the variable 'rho_ao_qs_kp'
                  ! (which is qs_env%rho%rho_ao_kp); density mixing will be done on an inverse-space grid
                  DO image = 1, nimages
                     DO ispin = 1, nspins
                        CALL dbcsr_copy(matrix_b=rho_ao_qs_kp(ispin, image)%matrix, &
                                        matrix_a=rho_ao_initial_kp(ispin, image)%matrix)
                     END DO
                  END DO

                  DO ispin = 1, nspins
                     CALL negf_copy_fm_submat_to_dbcsr(fm=rho_ao_new_fm(ispin)%matrix, &
                                                       matrix=rho_ao_qs_kp(ispin, 1)%matrix, &
                                                       atomlist_row=negf_control%atomlist_S_screening, &
                                                       atomlist_col=negf_control%atomlist_S_screening, &
                                                       subsys=subsys)
                  END DO
               END IF

               CALL qs_rho_update_rho(rho_struct, qs_env=qs_env)

               IF (negf_env%mixing_method >= gspace_mixing_nr) THEN
                  CALL gspace_mixing(qs_env, negf_env%mixing_method, negf_env%mixing_storage, &
                                     rho_struct, para_env, iter_count)
               END IF

               ! update KS-matrix
               CALL qs_ks_build_kohn_sham_matrix(qs_env, calculate_forces=.FALSE., just_energy=.FALSE.)

               ! extract blocks from the updated Kohn-Sham matrix
               DO ispin = 1, nspins
                  CALL negf_copy_sym_dbcsr_to_fm_submat(matrix=matrix_ks_qs_kp(ispin, 1)%matrix, &
                                                        fm=negf_env%h_s(ispin)%matrix, &
                                                        atomlist_row=negf_control%atomlist_S_screening, &
                                                        atomlist_col=negf_control%atomlist_S_screening, &
                                                        subsys=subsys, mpi_comm_global=para_env%group, &
                                                        do_upper_diag=.TRUE., do_lower=.TRUE.)

               END DO
            END DO

            IF (log_unit > 0) THEN
               IF (iter_count <= negf_control%max_scf) THEN
                  WRITE (log_unit, '(/,T11,1X,A,I0,A)') "*** NEGF run converged in ", iter_count, " iteration(s) ***"
               ELSE
                  WRITE (log_unit, '(/,T11,1X,A,I0,A)') "*** NEGF run did NOT converge after ", iter_count - 1, " iteration(s) ***"
               END IF
            END IF

            DO ispin = nspins, 1, -1
               CALL green_functions_cache_release(g_surf_circular(ispin))
               CALL green_functions_cache_release(g_surf_linear(ispin))
               CALL green_functions_cache_release(g_surf_nonequiv(ispin))
            END DO
            DEALLOCATE (g_surf_circular, g_surf_linear, g_surf_nonequiv)
         END IF
      END DO

      CALL cp_fm_release(ao_ao_fm_global)

      DO ispin = nspins, 1, -1
         CALL cp_fm_release(rho_ao_new_fm(ispin)%matrix)
         CALL cp_fm_release(rho_ao_delta_fm(ispin)%matrix)
      END DO
      DEALLOCATE (rho_ao_delta_fm, rho_ao_new_fm)

      DO image = 1, nimages
         DO ispin = 1, nspins
            CALL dbcsr_copy(matrix_b=matrix_ks_qs_kp(ispin, image)%matrix, matrix_a=matrix_ks_initial_kp(ispin, image)%matrix)
            CALL dbcsr_copy(matrix_b=rho_ao_qs_kp(ispin, image)%matrix, matrix_a=rho_ao_initial_kp(ispin, image)%matrix)

            CALL dbcsr_deallocate_matrix(matrix_ks_initial_kp(ispin, image)%matrix)
            CALL dbcsr_deallocate_matrix(rho_ao_initial_kp(ispin, image)%matrix)
            CALL dbcsr_deallocate_matrix(rho_ao_new_kp(ispin, image)%matrix)
         END DO
      END DO
      DEALLOCATE (matrix_ks_initial_kp, rho_ao_new_kp, rho_ao_initial_kp)

      IF (ASSOCIATED(matrix_s_fm)) &
         CALL cp_fm_release(matrix_s_fm)

      CALL timestop(handle)
   END SUBROUTINE converge_density

! **************************************************************************************************
!> \brief Compute the surface retarded Green's function at a set of points in parallel.
!> \param g_surf      set of surface Green's functions computed within the given parallel group
!> \param omega       list of energy points where the surface Green's function need to be computed
!> \param h0          diagonal block of the Kohn-Sham matrix (must be Hermitian)
!> \param s0          diagonal block of the overlap matrix (must be Hermitian)
!> \param h1          off-fiagonal block of the Kohn-Sham matrix
!> \param s1          off-fiagonal block of the overlap matrix
!> \param sub_env     NEGF parallel (sub)group environment
!> \param v_external  applied electric potential
!> \param conv        convergence threshold
!> \param transp      flag which indicates that the matrices h1 and s1 should be transposed
!> \par History
!>    * 07.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_surface_green_function_batch(g_surf, omega, h0, s0, h1, s1, sub_env, v_external, conv, transp)
      TYPE(cp_cfm_p_type), DIMENSION(:), INTENT(inout)   :: g_surf
      COMPLEX(kind=dp), DIMENSION(:), INTENT(in)         :: omega
      TYPE(cp_fm_type), POINTER                          :: h0, s0, h1, s1
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      REAL(kind=dp), INTENT(in)                          :: v_external, conv
      LOGICAL, INTENT(in)                                :: transp

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_surface_green_function_batch', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, igroup, ipoint, npoints
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(sancho_work_matrices_type)                    :: work

      CALL timeset(routineN, handle)
      npoints = SIZE(omega)

      NULLIFY (fm_struct)
      CALL cp_fm_get_info(s0, matrix_struct=fm_struct)
      CALL sancho_work_matrices_create(work, fm_struct)

      igroup = sub_env%group_distribution(sub_env%mepos_global)

      DO ipoint = 1, npoints
         NULLIFY (g_surf(ipoint)%matrix)
      END DO

      DO ipoint = igroup + 1, npoints, sub_env%ngroups
         IF (debug_this_module) THEN
            CPASSERT(.NOT. ASSOCIATED(g_surf(ipoint)%matrix))
         END IF
         CALL cp_cfm_create(g_surf(ipoint)%matrix, fm_struct)

         CALL do_sancho(g_surf(ipoint)%matrix, omega(ipoint) - v_external, &
                        h0, s0, h1, s1, conv, transp, work)
      END DO

      CALL sancho_work_matrices_release(work)
      CALL timestop(handle)
   END SUBROUTINE negf_surface_green_function_batch

! **************************************************************************************************
!> \brief Compute the retarded Green's function and related properties at a set of points in parallel.
!> \param omega              list of energy points
!> \param v_shift            shift in Hartree potential
!> \param ignore_bias        ignore v_external from negf_control
!> \param negf_env           NEGF environment
!> \param negf_control       NEGF control
!> \param sub_env            (sub)group environment
!> \param ispin              spin component to compute
!> \param g_surf_contacts    set of surface Green's functions for every contact that computed
!>                           within the given parallel group
!> \param g_ret_s            globally distributed matrices to store retarded Green's functions
!> \param g_ret_scale        scale factor for retarded Green's functions
!> \param gamma_contacts     2-D array of globally distributed matrices to store broadening matrices
!>                           for every contact ([n_contacts, npoints])
!> \param gret_gamma_gadv    2-D array of globally distributed matrices to store the spectral function:
!>                           g_ret_s * gamma * g_ret_s^C for every contact ([n_contacts, n_points])
!> \param dos                density of states at 'omega' ([n_points])
!> \param transm_coeff       transmission coefficients between two contacts 'transm_contact1'
!>                           and 'transm_contact2' computed at points 'omega' ([n_points])
!> \param transm_contact1    index of the first contact
!> \param transm_contact2    index of the second contact
!> \param just_contact       if present, compute the retarded Green's function of the system
!>                           lead1 -- device -- lead2. All 3 regions have the same Kohn-Sham
!>                           matrices which are taken from 'negf_env%contacts(just_contact)%h'.
!>                           Useful to apply NEGF procedure a single contact in order to compute
!>                           its Fermi level
!> \par History
!>    * 07.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_retarded_green_function_batch(omega, v_shift, ignore_bias, negf_env, negf_control, sub_env, ispin, &
                                                 g_surf_contacts, &
                                                 g_ret_s, g_ret_scale, gamma_contacts, gret_gamma_gadv, dos, &
                                                 transm_coeff, transm_contact1, transm_contact2, just_contact)
      COMPLEX(kind=dp), DIMENSION(:), INTENT(in)         :: omega
      REAL(kind=dp), INTENT(in)                          :: v_shift
      LOGICAL, INTENT(in)                                :: ignore_bias
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      INTEGER, INTENT(in)                                :: ispin
      TYPE(cp_cfm_p_type), DIMENSION(:, :), INTENT(in)   :: g_surf_contacts
      TYPE(cp_cfm_p_type), DIMENSION(:), INTENT(in), &
         OPTIONAL                                        :: g_ret_s
      COMPLEX(kind=dp), DIMENSION(:), INTENT(in), &
         OPTIONAL                                        :: g_ret_scale
      TYPE(cp_cfm_p_type), DIMENSION(:, :), INTENT(in), &
         OPTIONAL                                        :: gamma_contacts, gret_gamma_gadv
      REAL(kind=dp), DIMENSION(:), INTENT(out), OPTIONAL :: dos
      COMPLEX(kind=dp), DIMENSION(:), INTENT(out), &
         OPTIONAL                                        :: transm_coeff
      INTEGER, INTENT(in), OPTIONAL                      :: transm_contact1, transm_contact2, &
                                                            just_contact

      CHARACTER(len=*), PARAMETER :: routineN = 'negf_retarded_green_function_batch', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icontact, igroup, ipoint, &
                                                            ncontacts, npoints, nrows
      REAL(kind=dp)                                      :: v_external
      TYPE(copy_cfm_info_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: info1
      TYPE(copy_cfm_info_type), ALLOCATABLE, &
         DIMENSION(:, :)                                 :: info2
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:)     :: g_ret_s_group, self_energy_contacts, &
                                                            zwork1_contacts, zwork2_contacts
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:, :)  :: gamma_contacts_group, &
                                                            gret_gamma_gadv_group
      TYPE(cp_cfm_type), POINTER                         :: matrix_s_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: g_ret_imag, matrix_s
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      npoints = SIZE(omega)
      ncontacts = SIZE(negf_env%contacts)
      CPASSERT(SIZE(negf_control%contacts) == ncontacts)

      NULLIFY (matrix_s_cfm)

      IF (PRESENT(just_contact)) THEN
         CPASSERT(just_contact <= ncontacts)
         ncontacts = 2
      END IF

      CPASSERT(ncontacts >= 2)

      IF (ignore_bias) v_external = 0.0_dp

      IF (PRESENT(transm_coeff) .OR. PRESENT(transm_contact1) .OR. PRESENT(transm_contact2)) THEN
         CPASSERT(PRESENT(transm_coeff))
         CPASSERT(PRESENT(transm_contact1))
         CPASSERT(PRESENT(transm_contact2))
         CPASSERT(.NOT. PRESENT(just_contact))
      END IF

      ALLOCATE (self_energy_contacts(ncontacts), zwork1_contacts(ncontacts), zwork2_contacts(ncontacts))

      IF (PRESENT(just_contact)) THEN
         CALL cp_fm_get_info(negf_env%contacts(just_contact)%s_01, matrix_struct=fm_struct)
         DO icontact = 1, ncontacts
            NULLIFY (zwork1_contacts(icontact)%matrix, zwork2_contacts(icontact)%matrix)
            CALL cp_cfm_create(zwork1_contacts(icontact)%matrix, fm_struct)
            CALL cp_cfm_create(zwork2_contacts(icontact)%matrix, fm_struct)
         END DO

         CALL cp_fm_get_info(negf_env%contacts(just_contact)%s_00, nrow_global=nrows, matrix_struct=fm_struct)
         DO icontact = 1, ncontacts
            NULLIFY (self_energy_contacts(icontact)%matrix)
            CALL cp_cfm_create(self_energy_contacts(icontact)%matrix, fm_struct)
         END DO
      ELSE
         DO icontact = 1, ncontacts
            CALL cp_fm_get_info(negf_env%s_sc(icontact)%matrix, matrix_struct=fm_struct)
            NULLIFY (zwork1_contacts(icontact)%matrix, zwork2_contacts(icontact)%matrix)
            CALL cp_cfm_create(zwork1_contacts(icontact)%matrix, fm_struct)
            CALL cp_cfm_create(zwork2_contacts(icontact)%matrix, fm_struct)
         END DO

         CALL cp_fm_get_info(negf_env%s_s, nrow_global=nrows, matrix_struct=fm_struct)
         DO icontact = 1, ncontacts
            NULLIFY (self_energy_contacts(icontact)%matrix)
            CALL cp_cfm_create(self_energy_contacts(icontact)%matrix, fm_struct)
         END DO
      END IF

      IF (PRESENT(g_ret_s) .OR. PRESENT(gret_gamma_gadv) .OR. &
          PRESENT(dos) .OR. PRESENT(transm_coeff)) THEN
         ALLOCATE (g_ret_s_group(npoints))

         IF (sub_env%ngroups > 1 .OR. .NOT. PRESENT(g_ret_s)) THEN
            DO ipoint = 1, npoints
               NULLIFY (g_ret_s_group(ipoint)%matrix)
            END DO
         ELSE
            DO ipoint = 1, npoints
               g_ret_s_group(ipoint)%matrix => g_ret_s(ipoint)%matrix
               CALL cp_cfm_retain(g_ret_s_group(ipoint)%matrix)
            END DO
         END IF
      END IF

      IF (PRESENT(gamma_contacts) .OR. PRESENT(gret_gamma_gadv) .OR. PRESENT(transm_coeff)) THEN
         IF (debug_this_module .AND. PRESENT(gamma_contacts)) THEN
            CPASSERT(SIZE(gamma_contacts, 1) == ncontacts)
         END IF

         ALLOCATE (gamma_contacts_group(ncontacts, npoints))
         IF (sub_env%ngroups > 1 .OR. .NOT. PRESENT(gamma_contacts)) THEN
            DO ipoint = 1, npoints
               DO icontact = 1, ncontacts
                  NULLIFY (gamma_contacts_group(icontact, ipoint)%matrix)
               END DO
            END DO
         ELSE
            DO ipoint = 1, npoints
               DO icontact = 1, ncontacts
                  gamma_contacts_group(icontact, ipoint)%matrix => gamma_contacts(icontact, ipoint)%matrix
                  CALL cp_cfm_retain(gamma_contacts_group(icontact, ipoint)%matrix)
               END DO
            END DO
         END IF
      END IF

      IF (PRESENT(gret_gamma_gadv)) THEN
         IF (debug_this_module .AND. PRESENT(gret_gamma_gadv)) THEN
            CPASSERT(SIZE(gret_gamma_gadv, 1) == ncontacts)
         END IF

         ALLOCATE (gret_gamma_gadv_group(ncontacts, npoints))
         IF (sub_env%ngroups > 1) THEN
            DO ipoint = 1, npoints
               DO icontact = 1, ncontacts
                  NULLIFY (gret_gamma_gadv_group(icontact, ipoint)%matrix)
               END DO
            END DO
         ELSE
            DO ipoint = 1, npoints
               DO icontact = 1, ncontacts
                  gret_gamma_gadv_group(icontact, ipoint)%matrix => gret_gamma_gadv(icontact, ipoint)%matrix

                  IF (ASSOCIATED(gret_gamma_gadv_group(icontact, ipoint)%matrix)) &
                     CALL cp_cfm_retain(gret_gamma_gadv_group(icontact, ipoint)%matrix)
               END DO
            END DO
         END IF
      END IF

      igroup = sub_env%group_distribution(sub_env%mepos_global)

      DO ipoint = 1, npoints
         IF (ASSOCIATED(g_surf_contacts(1, ipoint)%matrix)) THEN
            IF (sub_env%ngroups > 1 .OR. .NOT. PRESENT(g_ret_s)) THEN
               ! create a group-specific matrix to store retarded Green's function if there are
               ! at least two parallel groups; otherwise pointers to group-specific matrices have
               ! already been initialised and they point to globally distributed matrices
               IF (ALLOCATED(g_ret_s_group)) THEN
                  CALL cp_cfm_create(g_ret_s_group(ipoint)%matrix, fm_struct)
               END IF
            END IF

            IF (sub_env%ngroups > 1 .OR. .NOT. PRESENT(gamma_contacts)) THEN
               IF (ALLOCATED(gamma_contacts_group)) THEN
                  DO icontact = 1, ncontacts
                     CALL cp_cfm_create(gamma_contacts_group(icontact, ipoint)%matrix, fm_struct)
                  END DO
               END IF
            END IF

            IF (sub_env%ngroups > 1) THEN
               IF (ALLOCATED(gret_gamma_gadv_group)) THEN
                  DO icontact = 1, ncontacts
                     IF (ASSOCIATED(gret_gamma_gadv(icontact, ipoint)%matrix)) &
                        CALL cp_cfm_create(gret_gamma_gadv_group(icontact, ipoint)%matrix, fm_struct)
                  END DO
               END IF
            END IF

            IF (PRESENT(just_contact)) THEN
               ! self energy of the "left" (1) and "right" contacts
               DO icontact = 1, ncontacts
                  CALL negf_contact_self_energy(self_energy_c=self_energy_contacts(icontact)%matrix, &
                                                omega=omega(ipoint), &
                                                g_surf_c=g_surf_contacts(icontact, ipoint)%matrix, &
                                                h_sc0=negf_env%contacts(just_contact)%h_01(ispin)%matrix, &
                                                s_sc0=negf_env%contacts(just_contact)%s_01, &
                                                zwork1=zwork1_contacts(icontact)%matrix, &
                                                zwork2=zwork2_contacts(icontact)%matrix, &
                                                transp=(icontact == 1))
               END DO
            ELSE
               ! contact self energies
               DO icontact = 1, ncontacts
                  IF (.NOT. ignore_bias) v_external = negf_control%contacts(icontact)%v_external

                  CALL negf_contact_self_energy(self_energy_c=self_energy_contacts(icontact)%matrix, &
                                                omega=omega(ipoint) - v_external, &
                                                g_surf_c=g_surf_contacts(icontact, ipoint)%matrix, &
                                                h_sc0=negf_env%h_sc(ispin, icontact)%matrix, &
                                                s_sc0=negf_env%s_sc(icontact)%matrix, &
                                                zwork1=zwork1_contacts(icontact)%matrix, &
                                                zwork2=zwork2_contacts(icontact)%matrix, &
                                                transp=.FALSE.)
               END DO
            END IF

            ! broadening matrices
            IF (ALLOCATED(gamma_contacts_group)) THEN
               DO icontact = 1, ncontacts
                  CALL negf_contact_broadening_matrix(gamma_c=gamma_contacts_group(icontact, ipoint)%matrix, &
                                                      self_energy_c=self_energy_contacts(icontact)%matrix)
               END DO
            END IF

            IF (ALLOCATED(g_ret_s_group)) THEN
               ! sum up self energies for all contacts
               DO icontact = 2, ncontacts
                  CALL cp_cfm_scale_and_add(z_one, self_energy_contacts(1)%matrix, z_one, self_energy_contacts(icontact)%matrix)
               END DO

               ! retarded Green's function for the scattering region
               IF (PRESENT(just_contact)) THEN
                  CALL negf_retarded_green_function(g_ret_s=g_ret_s_group(ipoint)%matrix, &
                                                    omega=omega(ipoint) - v_shift, &
                                                    self_energy_ret_sum=self_energy_contacts(1)%matrix, &
                                                    h_s=negf_env%contacts(just_contact)%h_00(ispin)%matrix, &
                                                    s_s=negf_env%contacts(just_contact)%s_00, &
                                                    v_hartree_s=null())
               ELSE IF (ignore_bias) THEN
                  CALL negf_retarded_green_function(g_ret_s=g_ret_s_group(ipoint)%matrix, &
                                                    omega=omega(ipoint) - v_shift, &
                                                    self_energy_ret_sum=self_energy_contacts(1)%matrix, &
                                                    h_s=negf_env%h_s(ispin)%matrix, &
                                                    s_s=negf_env%s_s, &
                                                    v_hartree_s=null())
               ELSE
                  CALL negf_retarded_green_function(g_ret_s=g_ret_s_group(ipoint)%matrix, &
                                                    omega=omega(ipoint) - v_shift, &
                                                    self_energy_ret_sum=self_energy_contacts(1)%matrix, &
                                                    h_s=negf_env%h_s(ispin)%matrix, &
                                                    s_s=negf_env%s_s, &
                                                    v_hartree_s=negf_env%v_hartree_s)
               END IF

               IF (PRESENT(g_ret_scale)) THEN
                  IF (g_ret_scale(ipoint) /= z_one) CALL cp_cfm_scale(g_ret_scale(ipoint), g_ret_s_group(ipoint)%matrix)
               END IF
            END IF

            IF (ALLOCATED(gret_gamma_gadv_group)) THEN
               ! we do not need contact self energies any longer, so we can use
               ! the array 'self_energy_contacts' as a set of work matrices
               DO icontact = 1, ncontacts
                  IF (ASSOCIATED(gret_gamma_gadv_group(icontact, ipoint)%matrix)) THEN
                     CALL cp_cfm_gemm('N', 'C', nrows, nrows, nrows, &
                                      z_one, gamma_contacts_group(icontact, ipoint)%matrix, &
                                      g_ret_s_group(ipoint)%matrix, &
                                      z_zero, self_energy_contacts(icontact)%matrix)
                     CALL cp_cfm_gemm('N', 'N', nrows, nrows, nrows, &
                                      z_one, g_ret_s_group(ipoint)%matrix, &
                                      self_energy_contacts(icontact)%matrix, &
                                      z_zero, gret_gamma_gadv_group(icontact, ipoint)%matrix)
                  END IF
               END DO
            END IF
         END IF
      END DO

      ! redistribute locally stored matrices
      IF (PRESENT(g_ret_s)) THEN
         IF (sub_env%ngroups > 1) THEN
            NULLIFY (para_env)
            DO ipoint = 1, npoints
               IF (ASSOCIATED(g_ret_s(ipoint)%matrix)) THEN
                  CALL cp_cfm_get_info(g_ret_s(ipoint)%matrix, para_env=para_env)
                  EXIT
               END IF
            END DO

            IF (ASSOCIATED(para_env)) THEN
               ALLOCATE (info1(npoints))

               DO ipoint = 1, npoints
                  IF (ASSOCIATED(g_ret_s(ipoint)%matrix)) THEN
                     CALL cp_cfm_start_copy_general(g_ret_s_group(ipoint)%matrix, &
                                                    g_ret_s(ipoint)%matrix, &
                                                    para_env, info1(ipoint))
                  END IF
               END DO

               DO ipoint = 1, npoints
                  IF (ASSOCIATED(g_ret_s(ipoint)%matrix)) THEN
                     CALL cp_cfm_finish_copy_general(g_ret_s(ipoint)%matrix, info1(ipoint))
                     IF (ASSOCIATED(g_ret_s_group(ipoint)%matrix)) &
                        CALL cp_cfm_cleanup_copy_general(g_ret_s_group(ipoint)%matrix, info1(ipoint))
                  END IF
               END DO

               DEALLOCATE (info1)
            END IF
         END IF
      END IF

      IF (PRESENT(gamma_contacts)) THEN
         IF (sub_env%ngroups > 1) THEN
            NULLIFY (para_env)
            pnt1: DO ipoint = 1, npoints
               DO icontact = 1, ncontacts
                  IF (ASSOCIATED(gamma_contacts(icontact, ipoint)%matrix)) THEN
                     CALL cp_cfm_get_info(gamma_contacts(icontact, ipoint)%matrix, para_env=para_env)
                     EXIT pnt1
                  END IF
               END DO
            END DO pnt1

            IF (ASSOCIATED(para_env)) THEN
               ALLOCATE (info2(ncontacts, npoints))

               DO ipoint = 1, npoints
                  DO icontact = 1, ncontacts
                     IF (ASSOCIATED(gamma_contacts(icontact, ipoint)%matrix)) THEN
                        CALL cp_cfm_start_copy_general(gamma_contacts_group(icontact, ipoint)%matrix, &
                                                       gamma_contacts(icontact, ipoint)%matrix, &
                                                       para_env, info2(icontact, ipoint))
                     END IF
                  END DO
               END DO

               DO ipoint = 1, npoints
                  DO icontact = 1, ncontacts
                     IF (ASSOCIATED(gamma_contacts(icontact, ipoint)%matrix)) THEN
                        CALL cp_cfm_finish_copy_general(gamma_contacts(icontact, ipoint)%matrix, info2(icontact, ipoint))
                        IF (ASSOCIATED(gamma_contacts_group(icontact, ipoint)%matrix)) THEN
                           CALL cp_cfm_cleanup_copy_general(gamma_contacts_group(icontact, ipoint)%matrix, &
                                                            info2(icontact, ipoint))
                        END IF
                     END IF
                  END DO
               END DO

               DEALLOCATE (info2)
            END IF
         END IF
      END IF

      IF (PRESENT(gret_gamma_gadv)) THEN
         IF (sub_env%ngroups > 1) THEN
            NULLIFY (para_env)

            pnt2: DO ipoint = 1, npoints
               DO icontact = 1, ncontacts
                  IF (ASSOCIATED(gret_gamma_gadv(icontact, ipoint)%matrix)) THEN
                     CALL cp_cfm_get_info(gret_gamma_gadv(icontact, ipoint)%matrix, para_env=para_env)
                     EXIT pnt2
                  END IF
               END DO
            END DO pnt2

            IF (ASSOCIATED(para_env)) THEN
               ALLOCATE (info2(ncontacts, npoints))

               DO ipoint = 1, npoints
                  DO icontact = 1, ncontacts
                     IF (ASSOCIATED(gret_gamma_gadv(icontact, ipoint)%matrix)) THEN
                        CALL cp_cfm_start_copy_general(gret_gamma_gadv_group(icontact, ipoint)%matrix, &
                                                       gret_gamma_gadv(icontact, ipoint)%matrix, &
                                                       para_env, info2(icontact, ipoint))
                     END IF
                  END DO
               END DO

               DO ipoint = 1, npoints
                  DO icontact = 1, ncontacts
                     IF (ASSOCIATED(gret_gamma_gadv(icontact, ipoint)%matrix)) THEN
                        CALL cp_cfm_finish_copy_general(gret_gamma_gadv(icontact, ipoint)%matrix, info2(icontact, ipoint))
                        IF (ASSOCIATED(gret_gamma_gadv_group(icontact, ipoint)%matrix)) THEN
                           CALL cp_cfm_cleanup_copy_general(gret_gamma_gadv_group(icontact, ipoint)%matrix, &
                                                            info2(icontact, ipoint))
                        END IF
                     END IF
                  END DO
               END DO

               DEALLOCATE (info2)
            END IF
         END IF
      END IF

      IF (PRESENT(dos)) THEN
         dos(:) = 0.0_dp

         IF (PRESENT(just_contact)) THEN
            matrix_s => negf_env%contacts(just_contact)%s_00
         ELSE
            matrix_s => negf_env%s_s
         END IF

         CALL cp_fm_get_info(matrix_s, matrix_struct=fm_struct)
         NULLIFY (g_ret_imag)
         CALL cp_fm_create(g_ret_imag, fm_struct)

         DO ipoint = 1, npoints
            IF (ASSOCIATED(g_ret_s_group(ipoint)%matrix)) THEN
               CALL cp_cfm_to_fm(g_ret_s_group(ipoint)%matrix, mtargeti=g_ret_imag)
               CALL cp_fm_trace(g_ret_imag, matrix_s, dos(ipoint))
               IF (sub_env%para_env%mepos /= 0) dos(ipoint) = 0.0_dp
            END IF
         END DO

         CALL cp_fm_release(g_ret_imag)

         CALL mp_sum(dos, sub_env%mpi_comm_global)
         dos(:) = -1.0_dp/pi*dos(:)
      END IF

      IF (PRESENT(transm_coeff)) THEN
         transm_coeff(:) = z_zero

         DO ipoint = 1, npoints
            IF (ASSOCIATED(g_ret_s_group(ipoint)%matrix)) THEN
               ! gamma_1 * g_adv_s * gamma_2
               CALL cp_cfm_gemm('N', 'C', nrows, nrows, nrows, &
                                z_one, gamma_contacts_group(transm_contact1, ipoint)%matrix, &
                                g_ret_s_group(ipoint)%matrix, &
                                z_zero, self_energy_contacts(transm_contact1)%matrix)
               CALL cp_cfm_gemm('N', 'N', nrows, nrows, nrows, &
                                z_one, self_energy_contacts(transm_contact1)%matrix, &
                                gamma_contacts_group(transm_contact2, ipoint)%matrix, &
                                z_zero, self_energy_contacts(transm_contact2)%matrix)

               !  Trace[ g_ret_s * gamma_1 * g_adv_s * gamma_2 ]
               CALL cp_cfm_trace(g_ret_s_group(ipoint)%matrix, &
                                 self_energy_contacts(transm_contact2)%matrix, &
                                 transm_coeff(ipoint))
               IF (sub_env%para_env%mepos /= 0) transm_coeff(ipoint) = 0.0_dp
            END IF
         END DO

         ! transmission coefficients are scaled by 2/pi
         CALL mp_sum(transm_coeff, sub_env%mpi_comm_global)
         !transm_coeff(:) = 0.5_dp/pi*transm_coeff(:)
      END IF

      ! -- deallocate temporary matrices
      IF (ALLOCATED(g_ret_s_group)) THEN
         DO ipoint = npoints, 1, -1
            IF (ASSOCIATED(g_ret_s_group(ipoint)%matrix)) CALL cp_cfm_release(g_ret_s_group(ipoint)%matrix)
         END DO
         DEALLOCATE (g_ret_s_group)
      END IF

      IF (ASSOCIATED(matrix_s_cfm)) CALL cp_cfm_release(matrix_s_cfm)

      IF (ALLOCATED(gamma_contacts_group)) THEN
         DO ipoint = npoints, 1, -1
            DO icontact = ncontacts, 1, -1
               IF (ASSOCIATED(gamma_contacts_group(icontact, ipoint)%matrix)) &
                  CALL cp_cfm_release(gamma_contacts_group(icontact, ipoint)%matrix)
            END DO
         END DO
         DEALLOCATE (gamma_contacts_group)
      END IF

      IF (ALLOCATED(gret_gamma_gadv_group)) THEN
         DO ipoint = npoints, 1, -1
            DO icontact = ncontacts, 1, -1
               IF (ASSOCIATED(gret_gamma_gadv_group(icontact, ipoint)%matrix)) &
                  CALL cp_cfm_release(gret_gamma_gadv_group(icontact, ipoint)%matrix)
            END DO
         END DO
         DEALLOCATE (gret_gamma_gadv_group)
      END IF

      IF (ALLOCATED(self_energy_contacts)) THEN
         DO icontact = ncontacts, 1, -1
            IF (ASSOCIATED(self_energy_contacts(icontact)%matrix)) &
               CALL cp_cfm_release(self_energy_contacts(icontact)%matrix)
         END DO
         DEALLOCATE (self_energy_contacts)
      END IF

      IF (ALLOCATED(zwork1_contacts)) THEN
         DO icontact = ncontacts, 1, -1
            IF (ASSOCIATED(zwork1_contacts(icontact)%matrix)) &
               CALL cp_cfm_release(zwork1_contacts(icontact)%matrix)
         END DO
         DEALLOCATE (zwork1_contacts)
      END IF

      IF (ALLOCATED(zwork2_contacts)) THEN
         DO icontact = ncontacts, 1, -1
            IF (ASSOCIATED(zwork2_contacts(icontact)%matrix)) &
               CALL cp_cfm_release(zwork2_contacts(icontact)%matrix)
         END DO
         DEALLOCATE (zwork2_contacts)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_retarded_green_function_batch

! **************************************************************************************************
!> \brief Fermi function (exp(E/(kT)) + 1) ^ {-1} .
!> \param omega       'energy' point on the complex plane
!> \param temperature temperature in atomic units
!> \return value
!> \par History
!>    * 05.2017 created [Sergey Chulkov]
! **************************************************************************************************
   PURE FUNCTION fermi_function(omega, temperature) RESULT(val)
      COMPLEX(kind=dp), INTENT(in)                       :: omega
      REAL(kind=dp), INTENT(in)                          :: temperature
      COMPLEX(kind=dp)                                   :: val

      REAL(kind=dp), PARAMETER :: max_ln_omega_over_T = LOG(HUGE(0.0_dp))/16.0_dp

      IF (REAL(omega, kind=dp) <= temperature*max_ln_omega_over_T) THEN
         ! exp(omega / T) < huge(0), so EXP() should not return infinity
         val = z_one/(EXP(omega/temperature) + z_one)
      ELSE
         val = z_zero
      END IF
   END FUNCTION fermi_function

! **************************************************************************************************
!> \brief Compute contribution to the density matrix from the poles of the Fermi function.
!> \param rho_ao_fm     density matrix (initialised on exit)
!> \param v_shift       shift in Hartree potential
!> \param ignore_bias   ignore v_external from negf_control
!> \param negf_env      NEGF environment
!> \param negf_control  NEGF control
!> \param sub_env       NEGF parallel (sub)group environment
!> \param ispin         spin conponent to proceed
!> \param base_contact  index of the reference contact
!> \param just_contact  ...
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_init_rho_equiv_residuals(rho_ao_fm, v_shift, ignore_bias, negf_env, &
                                            negf_control, sub_env, ispin, base_contact, just_contact)
      TYPE(cp_fm_type), POINTER                          :: rho_ao_fm
      REAL(kind=dp), INTENT(in)                          :: v_shift
      LOGICAL, INTENT(in)                                :: ignore_bias
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      INTEGER, INTENT(in)                                :: ispin, base_contact
      INTEGER, INTENT(in), OPTIONAL                      :: just_contact

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_init_rho_equiv_residuals', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: omega
      INTEGER                                            :: handle, icontact, ipole, ncontacts, &
                                                            npoles
      REAL(kind=dp)                                      :: mu_base, pi_temperature, temperature, &
                                                            v_external
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:)     :: g_ret_s
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(green_functions_cache_type)                   :: g_surf_cache

      CALL timeset(routineN, handle)

      temperature = negf_control%contacts(base_contact)%temperature
      IF (ignore_bias) THEN
         mu_base = negf_control%contacts(base_contact)%fermi_level
         v_external = 0.0_dp
      ELSE
         mu_base = negf_control%contacts(base_contact)%fermi_level + negf_control%contacts(base_contact)%v_external
      END IF

      pi_temperature = pi*temperature
      npoles = negf_control%delta_npoles

      ncontacts = SIZE(negf_env%contacts)
      CPASSERT(base_contact <= ncontacts)
      IF (PRESENT(just_contact)) THEN
         ncontacts = 2
         CPASSERT(just_contact == base_contact)
      END IF

      IF (npoles > 0) THEN
         CALL cp_fm_get_info(rho_ao_fm, para_env=para_env, matrix_struct=fm_struct)

         ALLOCATE (omega(npoles), g_ret_s(npoles))

         DO ipole = 1, npoles
            NULLIFY (g_ret_s(ipole)%matrix)
            CALL cp_cfm_create(g_ret_s(ipole)%matrix, fm_struct)

            omega(ipole) = CMPLX(mu_base, REAL(2*ipole - 1, kind=dp)*pi_temperature, kind=dp)
         END DO

         CALL green_functions_cache_expand(g_surf_cache, ncontacts, npoles)

         IF (PRESENT(just_contact)) THEN
            ! do not apply the external potential when computing the Fermi level of a bulk contact.
            ! We are using a fictitious electronic device, which identical to the bulk contact in question;
            ! icontact == 1 corresponds to the "left" contact, so the matrices h_01 and s_01 needs to be transposed,
            ! while icontact == 2 correspond to the "right" contact and we should use the matrices h_01 and s_01 as is.
            DO icontact = 1, ncontacts
               CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, :), &
                                                      omega=omega(:), &
                                                      h0=negf_env%contacts(just_contact)%h_00(ispin)%matrix, &
                                                      s0=negf_env%contacts(just_contact)%s_00, &
                                                      h1=negf_env%contacts(just_contact)%h_01(ispin)%matrix, &
                                                      s1=negf_env%contacts(just_contact)%s_01, &
                                                      sub_env=sub_env, v_external=0.0_dp, &
                                                      conv=negf_control%conv_green, transp=(icontact == 1))
            END DO
         ELSE
            DO icontact = 1, ncontacts
               IF (.NOT. ignore_bias) v_external = negf_control%contacts(icontact)%v_external

               CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, :), &
                                                      omega=omega(:), &
                                                      h0=negf_env%contacts(icontact)%h_00(ispin)%matrix, &
                                                      s0=negf_env%contacts(icontact)%s_00, &
                                                      h1=negf_env%contacts(icontact)%h_01(ispin)%matrix, &
                                                      s1=negf_env%contacts(icontact)%s_01, &
                                                      sub_env=sub_env, &
                                                      v_external=v_external, &
                                                      conv=negf_control%conv_green, transp=.FALSE.)
            END DO
         END IF

         CALL negf_retarded_green_function_batch(omega=omega(:), &
                                                 v_shift=v_shift, &
                                                 ignore_bias=ignore_bias, &
                                                 negf_env=negf_env, &
                                                 negf_control=negf_control, &
                                                 sub_env=sub_env, &
                                                 ispin=ispin, &
                                                 g_surf_contacts=g_surf_cache%g_surf_contacts, &
                                                 g_ret_s=g_ret_s, &
                                                 just_contact=just_contact)

         CALL green_functions_cache_release(g_surf_cache)

         DO ipole = 2, npoles
            CALL cp_cfm_scale_and_add(z_one, g_ret_s(1)%matrix, z_one, g_ret_s(ipole)%matrix)
         END DO

         !Re(-i * (-2*pi*i*kB*T/(-pi) * [Re(G)+i*Im(G)]) == 2*kB*T * Re(G)
         CALL cp_cfm_to_fm(g_ret_s(1)%matrix, mtargetr=rho_ao_fm)
         CALL cp_fm_scale(2.0_dp*temperature, rho_ao_fm)

         DO ipole = npoles, 1, -1
            CALL cp_cfm_release(g_ret_s(ipole)%matrix)
         END DO
         DEALLOCATE (g_ret_s, omega)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_init_rho_equiv_residuals

! **************************************************************************************************
!> \brief Compute equilibrium contribution to the density matrix.
!> \param rho_ao_fm       density matrix (initialised on exit)
!> \param stats           integration statistics (updated on exit)
!> \param v_shift         shift in Hartree potential
!> \param ignore_bias     ignore v_external from negf_control
!> \param negf_env        NEGF environment
!> \param negf_control    NEGF control
!> \param sub_env         NEGF parallel (sub)group environment
!> \param ispin           spin conponent to proceed
!> \param base_contact    index of the reference contact
!> \param integr_lbound   integration lower bound
!> \param integr_ubound   integration upper bound
!> \param matrix_s_global globally distributed overlap matrix
!> \param is_circular     compute the integral along the circular path
!> \param g_surf_cache    set of precomputed surface Green's functions (updated on exit)
!> \param just_contact    ...
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_add_rho_equiv_low(rho_ao_fm, stats, v_shift, ignore_bias, negf_env, negf_control, sub_env, &
                                     ispin, base_contact, integr_lbound, integr_ubound, matrix_s_global, &
                                     is_circular, g_surf_cache, just_contact)
      TYPE(cp_fm_type), POINTER                          :: rho_ao_fm
      TYPE(integration_status_type), INTENT(inout)       :: stats
      REAL(kind=dp), INTENT(in)                          :: v_shift
      LOGICAL, INTENT(in)                                :: ignore_bias
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      INTEGER, INTENT(in)                                :: ispin, base_contact
      COMPLEX(kind=dp), INTENT(in)                       :: integr_lbound, integr_ubound
      TYPE(cp_fm_type), POINTER                          :: matrix_s_global
      LOGICAL, INTENT(in)                                :: is_circular
      TYPE(green_functions_cache_type), INTENT(inout)    :: g_surf_cache
      INTEGER, INTENT(in), OPTIONAL                      :: just_contact

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_add_rho_equiv_low', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: xnodes, zscale
      INTEGER :: handle, icontact, interval_id, ipoint, max_points, min_points, ncontacts, &
         npoints, npoints_exist, npoints_tmp, npoints_total, shape_id
      LOGICAL                                            :: do_surface_green
      REAL(kind=dp)                                      :: conv_integr, mu_base, temperature, &
                                                            v_external
      TYPE(ccquad_type)                                  :: cc_env
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:)     :: zdata, zdata_tmp
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: integral_imag
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(simpsonrule_type)                             :: sr_env

      CALL timeset(routineN, handle)

      ! convergence criteria for the integral of the retarded Green's function. This integral needs to be
      ! computed for both spin-components and needs to be scaled by -1/pi to obtain the electron density.
      conv_integr = 0.5_dp*negf_control%conv_density*pi

      IF (ignore_bias) THEN
         mu_base = negf_control%contacts(base_contact)%fermi_level
         v_external = 0.0_dp
      ELSE
         mu_base = negf_control%contacts(base_contact)%fermi_level + negf_control%contacts(base_contact)%v_external
      END IF

      min_points = negf_control%integr_min_points
      max_points = negf_control%integr_max_points
      temperature = negf_control%contacts(base_contact)%temperature

      ncontacts = SIZE(negf_env%contacts)
      CPASSERT(base_contact <= ncontacts)
      IF (PRESENT(just_contact)) THEN
         ncontacts = 2
         CPASSERT(just_contact == base_contact)
      END IF

      do_surface_green = .NOT. ALLOCATED(g_surf_cache%tnodes)

      IF (do_surface_green) THEN
         npoints = min_points
      ELSE
         npoints = SIZE(g_surf_cache%tnodes)
      END IF
      npoints_total = 0

      NULLIFY (integral_imag)
      CALL cp_fm_get_info(rho_ao_fm, para_env=para_env, matrix_struct=fm_struct)
      CALL cp_fm_create(integral_imag, fm_struct)

      SELECT CASE (negf_control%integr_method)
      CASE (negfint_method_cc)
         ! Adaptive Clenshaw-Curtis method
         ALLOCATE (xnodes(npoints))

         IF (is_circular) THEN
            shape_id = cc_shape_arc
            interval_id = cc_interval_full
         ELSE
            shape_id = cc_shape_linear
            interval_id = cc_interval_half
         END IF

         IF (do_surface_green) THEN
            CALL ccquad_init(cc_env, xnodes, npoints, integr_lbound, integr_ubound, &
                             interval_id, shape_id, matrix_s_global)
         ELSE
            CALL ccquad_init(cc_env, xnodes, npoints, integr_lbound, integr_ubound, &
                             interval_id, shape_id, matrix_s_global, tnodes_restart=g_surf_cache%tnodes)
         END IF

         ALLOCATE (zdata(npoints))
         DO ipoint = 1, npoints
            NULLIFY (zdata(ipoint)%matrix)
            CALL cp_cfm_create(zdata(ipoint)%matrix, fm_struct)
         END DO

         DO
            IF (do_surface_green) THEN
               CALL green_functions_cache_expand(g_surf_cache, ncontacts, npoints)

               IF (PRESENT(just_contact)) THEN
                  ! do not apply the external potential when computing the Fermi level of a bulk contact.
                  DO icontact = 1, ncontacts
                     CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, npoints_total + 1:), &
                                                            omega=xnodes(1:npoints), &
                                                            h0=negf_env%contacts(just_contact)%h_00(ispin)%matrix, &
                                                            s0=negf_env%contacts(just_contact)%s_00, &
                                                            h1=negf_env%contacts(just_contact)%h_01(ispin)%matrix, &
                                                            s1=negf_env%contacts(just_contact)%s_01, &
                                                            sub_env=sub_env, v_external=0.0_dp, &
                                                            conv=negf_control%conv_green, transp=(icontact == 1))
                  END DO
               ELSE
                  DO icontact = 1, ncontacts
                     IF (.NOT. ignore_bias) v_external = negf_control%contacts(icontact)%v_external

                     CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, npoints_total + 1:), &
                                                            omega=xnodes(1:npoints), &
                                                            h0=negf_env%contacts(icontact)%h_00(ispin)%matrix, &
                                                            s0=negf_env%contacts(icontact)%s_00, &
                                                            h1=negf_env%contacts(icontact)%h_01(ispin)%matrix, &
                                                            s1=negf_env%contacts(icontact)%s_01, &
                                                            sub_env=sub_env, &
                                                            v_external=v_external, &
                                                            conv=negf_control%conv_green, transp=.FALSE.)
                  END DO
               END IF
            END IF

            ALLOCATE (zscale(npoints))

            IF (temperature >= 0.0_dp) THEN
               DO ipoint = 1, npoints
                  zscale(ipoint) = fermi_function(xnodes(ipoint) - mu_base, temperature)
               END DO
            ELSE
               zscale(:) = z_one
            END IF

            CALL negf_retarded_green_function_batch(omega=xnodes(1:npoints), &
                                                    v_shift=v_shift, &
                                                    ignore_bias=ignore_bias, &
                                                    negf_env=negf_env, &
                                                    negf_control=negf_control, &
                                                    sub_env=sub_env, &
                                                    ispin=ispin, &
                                                    g_surf_contacts=g_surf_cache%g_surf_contacts(:, npoints_total + 1:), &
                                                    g_ret_s=zdata(1:npoints), &
                                                    g_ret_scale=zscale(1:npoints), &
                                                    just_contact=just_contact)

            DEALLOCATE (xnodes, zscale)
            npoints_total = npoints_total + npoints

            CALL ccquad_reduce_and_append_zdata(cc_env, zdata)
            CALL MOVE_ALLOC(zdata, zdata_tmp)

            CALL ccquad_refine_integral(cc_env)

            IF (cc_env%error <= conv_integr) EXIT
            IF (2*(npoints_total - 1) + 1 > max_points) EXIT

            ! all cached points have been reused at the first iteration;
            ! we need to compute surface Green's function at extra points if the integral has not been converged
            do_surface_green = .TRUE.

            npoints_tmp = npoints
            CALL ccquad_double_number_of_points(cc_env, xnodes)
            npoints = SIZE(xnodes)

            ALLOCATE (zdata(npoints))

            npoints_exist = 0
            DO ipoint = 1, npoints_tmp
               IF (ASSOCIATED(zdata_tmp(ipoint)%matrix)) THEN
                  npoints_exist = npoints_exist + 1
                  zdata(npoints_exist)%matrix => zdata_tmp(ipoint)%matrix
               END IF
            END DO
            DEALLOCATE (zdata_tmp)

            DO ipoint = npoints_exist + 1, npoints
               NULLIFY (zdata(ipoint)%matrix)
               CALL cp_cfm_create(zdata(ipoint)%matrix, fm_struct)
            END DO
         END DO

         ! the obtained integral will be scaled by -1/pi, so scale the error extimate as well
         stats%error = stats%error + cc_env%error/pi

         DO ipoint = SIZE(zdata_tmp), 1, -1
            IF (ASSOCIATED(zdata_tmp(ipoint)%matrix)) CALL cp_cfm_release(zdata_tmp(ipoint)%matrix)
         END DO
         DEALLOCATE (zdata_tmp)

         CALL cp_cfm_to_fm(cc_env%integral, mtargeti=integral_imag)

         ! keep the cache
         IF (do_surface_green) THEN
            CALL green_functions_cache_reorder(g_surf_cache, cc_env%tnodes)
         END IF
         CALL ccquad_release(cc_env)

      CASE (negfint_method_simpson)
         ! Adaptive Simpson's rule method
         ALLOCATE (xnodes(npoints), zdata(npoints), zscale(npoints))

         IF (is_circular) THEN
            shape_id = sr_shape_arc
         ELSE
            shape_id = sr_shape_linear
         END IF

         IF (do_surface_green) THEN
            CALL simpsonrule_init(sr_env, xnodes, npoints, integr_lbound, integr_ubound, &
                                  shape_id, conv_integr, matrix_s_global)
         ELSE
            CALL simpsonrule_init(sr_env, xnodes, npoints, integr_lbound, integr_ubound, &
                                  shape_id, conv_integr, matrix_s_global, tnodes_restart=g_surf_cache%tnodes)
         END IF

         DO WHILE (npoints > 0 .AND. npoints_total < max_points)
            DO ipoint = 1, npoints
               NULLIFY (zdata(ipoint)%matrix)
               CALL cp_cfm_create(zdata(ipoint)%matrix, fm_struct)
            END DO

            IF (do_surface_green) THEN
               CALL green_functions_cache_expand(g_surf_cache, ncontacts, npoints)

               IF (PRESENT(just_contact)) THEN
                  ! do not apply the external potential when computing the Fermi level of a bulk contact.
                  DO icontact = 1, ncontacts
                     CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, npoints_total + 1:), &
                                                            omega=xnodes(1:npoints), &
                                                            h0=negf_env%contacts(just_contact)%h_00(ispin)%matrix, &
                                                            s0=negf_env%contacts(just_contact)%s_00, &
                                                            h1=negf_env%contacts(just_contact)%h_01(ispin)%matrix, &
                                                            s1=negf_env%contacts(just_contact)%s_01, &
                                                            sub_env=sub_env, v_external=0.0_dp, &
                                                            conv=negf_control%conv_green, transp=(icontact == 1))
                  END DO
               ELSE
                  DO icontact = 1, ncontacts
                     IF (.NOT. ignore_bias) v_external = negf_control%contacts(icontact)%v_external

                     CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, npoints_total + 1:), &
                                                            omega=xnodes(1:npoints), &
                                                            h0=negf_env%contacts(icontact)%h_00(ispin)%matrix, &
                                                            s0=negf_env%contacts(icontact)%s_00, &
                                                            h1=negf_env%contacts(icontact)%h_01(ispin)%matrix, &
                                                            s1=negf_env%contacts(icontact)%s_01, &
                                                            sub_env=sub_env, &
                                                            v_external=v_external, &
                                                            conv=negf_control%conv_green, transp=.FALSE.)
                  END DO
               END IF
            END IF

            IF (temperature >= 0.0_dp) THEN
               DO ipoint = 1, npoints
                  zscale(ipoint) = fermi_function(xnodes(ipoint) - mu_base, temperature)
               END DO
            ELSE
               zscale(:) = z_one
            END IF

            CALL negf_retarded_green_function_batch(omega=xnodes(1:npoints), &
                                                    v_shift=v_shift, &
                                                    ignore_bias=ignore_bias, &
                                                    negf_env=negf_env, &
                                                    negf_control=negf_control, &
                                                    sub_env=sub_env, &
                                                    ispin=ispin, &
                                                    g_surf_contacts=g_surf_cache%g_surf_contacts(:, npoints_total + 1:), &
                                                    g_ret_s=zdata(1:npoints), &
                                                    g_ret_scale=zscale(1:npoints), &
                                                    just_contact=just_contact)

            npoints_total = npoints_total + npoints

            CALL simpsonrule_refine_integral(sr_env, zdata(1:npoints))

            IF (sr_env%error <= conv_integr) EXIT

            ! all cached points have been reused at the first iteration;
            ! if the integral has not been converged, turn on the 'do_surface_green' flag
            ! in order to add more points
            do_surface_green = .TRUE.

            npoints = max_points - npoints_total
            IF (npoints <= 0) EXIT
            IF (npoints > SIZE(xnodes)) npoints = SIZE(xnodes)

            CALL simpsonrule_get_next_nodes(sr_env, xnodes, npoints)
         END DO

         ! the obtained integral will be scaled by -1/pi, so scale the error extimate as well
         stats%error = stats%error + sr_env%error/pi

         CALL cp_cfm_to_fm(sr_env%integral, mtargeti=integral_imag)

         ! keep the cache
         IF (do_surface_green) THEN
            CALL green_functions_cache_reorder(g_surf_cache, sr_env%tnodes)
         END IF

         CALL simpsonrule_release(sr_env)
         DEALLOCATE (xnodes, zdata, zscale)

      CASE DEFAULT
         CPABORT("Unimplemented integration method")
      END SELECT

      stats%npoints = stats%npoints + npoints_total

      CALL cp_fm_scale_and_add(1.0_dp, rho_ao_fm, -1.0_dp/pi, integral_imag)
      CALL cp_fm_release(integral_imag)

      CALL timestop(handle)
   END SUBROUTINE negf_add_rho_equiv_low

! **************************************************************************************************
!> \brief Compute non-equilibrium contribution to the density matrix.
!> \param rho_ao_fm       density matrix (initialised on exit)
!> \param stats           integration statistics (updated on exit)
!> \param v_shift         shift in Hartree potential
!> \param negf_env        NEGF environment
!> \param negf_control    NEGF control
!> \param sub_env         NEGF parallel (sub)group environment
!> \param ispin           spin conponent to proceed
!> \param base_contact    index of the reference contact
!> \param matrix_s_global globally distributed overlap matrix
!> \param g_surf_cache    set of precomputed surface Green's functions (updated on exit)
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_add_rho_nonequiv(rho_ao_fm, stats, v_shift, negf_env, negf_control, sub_env, &
                                    ispin, base_contact, matrix_s_global, g_surf_cache)
      TYPE(cp_fm_type), POINTER                          :: rho_ao_fm
      TYPE(integration_status_type), INTENT(inout)       :: stats
      REAL(kind=dp), INTENT(in)                          :: v_shift
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      INTEGER, INTENT(in)                                :: ispin, base_contact
      TYPE(cp_fm_type), POINTER                          :: matrix_s_global
      TYPE(green_functions_cache_type), INTENT(inout)    :: g_surf_cache

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_add_rho_nonequiv', &
         routineP = moduleN//':'//routineN

      COMPLEX(kind=dp)                                   :: fermi_base, fermi_contact, &
                                                            integr_lbound, integr_ubound
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: xnodes
      INTEGER                                            :: handle, icontact, ipoint, jcontact, &
                                                            max_points, min_points, ncontacts, &
                                                            npoints, npoints_total
      LOGICAL                                            :: do_surface_green
      REAL(kind=dp)                                      :: conv_density, conv_integr, eta, &
                                                            ln_conv_density, mu_base, mu_contact, &
                                                            temperature_base, temperature_contact
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:, :)  :: zdata
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: integral_real
      TYPE(simpsonrule_type)                             :: sr_env

      CALL timeset(routineN, handle)

      ncontacts = SIZE(negf_env%contacts)
      CPASSERT(base_contact <= ncontacts)

      ! the current subroutine works for the general case as well, but the Poisson solver does not
      IF (ncontacts > 2) THEN
         CPABORT("Poisson solver does not support the general NEGF setup (>2 contacts).")
      END IF

      mu_base = negf_control%contacts(base_contact)%fermi_level + negf_control%contacts(base_contact)%v_external
      min_points = negf_control%integr_min_points
      max_points = negf_control%integr_max_points
      temperature_base = negf_control%contacts(base_contact)%temperature
      eta = negf_control%eta
      conv_density = negf_control%conv_density
      ln_conv_density = LOG(conv_density)

      ! convergence criteria for the integral. This integral needs to be computed for both
      ! spin-components and needs to be scaled by -1/pi to obtain the electron density.
      conv_integr = 0.5_dp*conv_density*pi

      DO icontact = 1, ncontacts
         IF (icontact /= base_contact) THEN
            mu_contact = negf_control%contacts(icontact)%fermi_level + negf_control%contacts(icontact)%v_external
            temperature_contact = negf_control%contacts(icontact)%temperature

            integr_lbound = CMPLX(MIN(mu_base + ln_conv_density*temperature_base, &
                                      mu_contact + ln_conv_density*temperature_contact), eta, kind=dp)
            integr_ubound = CMPLX(MAX(mu_base - ln_conv_density*temperature_base, &
                                      mu_contact - ln_conv_density*temperature_contact), eta, kind=dp)

            do_surface_green = .NOT. ALLOCATED(g_surf_cache%tnodes)

            IF (do_surface_green) THEN
               npoints = min_points
            ELSE
               npoints = SIZE(g_surf_cache%tnodes)
            END IF
            npoints_total = 0

            CALL cp_fm_get_info(rho_ao_fm, matrix_struct=fm_struct)

            ALLOCATE (xnodes(npoints), zdata(ncontacts, npoints))
            DO ipoint = 1, npoints
               DO jcontact = 1, ncontacts
                  NULLIFY (zdata(jcontact, ipoint)%matrix)
               END DO
            END DO

            IF (do_surface_green) THEN
               CALL simpsonrule_init(sr_env, xnodes, npoints, integr_lbound, integr_ubound, &
                                     sr_shape_linear, conv_integr, matrix_s_global)
            ELSE
               CALL simpsonrule_init(sr_env, xnodes, npoints, integr_lbound, integr_ubound, &
                                     sr_shape_linear, conv_integr, matrix_s_global, tnodes_restart=g_surf_cache%tnodes)
            END IF

            DO WHILE (npoints > 0 .AND. npoints_total < max_points)
               IF (do_surface_green) THEN
                  CALL green_functions_cache_expand(g_surf_cache, ncontacts, npoints)

                  DO jcontact = 1, ncontacts
                     CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(jcontact, npoints_total + 1:), &
                                                            omega=xnodes(1:npoints), &
                                                            h0=negf_env%contacts(jcontact)%h_00(ispin)%matrix, &
                                                            s0=negf_env%contacts(jcontact)%s_00, &
                                                            h1=negf_env%contacts(jcontact)%h_01(ispin)%matrix, &
                                                            s1=negf_env%contacts(jcontact)%s_01, &
                                                            sub_env=sub_env, &
                                                            v_external=negf_control%contacts(jcontact)%v_external, &
                                                            conv=negf_control%conv_green, transp=.FALSE.)
                  END DO
               END IF

               DO ipoint = 1, npoints
                  CALL cp_cfm_create(zdata(icontact, ipoint)%matrix, fm_struct)
               END DO

               CALL negf_retarded_green_function_batch(omega=xnodes(1:npoints), &
                                                       v_shift=v_shift, &
                                                       ignore_bias=.FALSE., &
                                                       negf_env=negf_env, &
                                                       negf_control=negf_control, &
                                                       sub_env=sub_env, &
                                                       ispin=ispin, &
                                                       g_surf_contacts=g_surf_cache%g_surf_contacts(:, npoints_total + 1:), &
                                                       gret_gamma_gadv=zdata(:, 1:npoints))

               DO ipoint = 1, npoints
                  fermi_base = fermi_function(CMPLX(REAL(xnodes(ipoint), kind=dp) - mu_base, 0.0_dp, kind=dp), &
                                              temperature_base)
                  fermi_contact = fermi_function(CMPLX(REAL(xnodes(ipoint), kind=dp) - mu_contact, 0.0_dp, kind=dp), &
                                                 temperature_contact)
                  CALL cp_cfm_scale(fermi_contact - fermi_base, zdata(icontact, ipoint)%matrix)
               END DO

               npoints_total = npoints_total + npoints

               CALL simpsonrule_refine_integral(sr_env, zdata(icontact, 1:npoints))

               IF (sr_env%error <= conv_integr) EXIT

               ! not enought cached points to achieve target accuracy
               do_surface_green = .TRUE.

               npoints = max_points - npoints_total
               IF (npoints <= 0) EXIT
               IF (npoints > SIZE(xnodes)) npoints = SIZE(xnodes)

               CALL simpsonrule_get_next_nodes(sr_env, xnodes, npoints)
            END DO

            NULLIFY (integral_real)
            CALL cp_fm_create(integral_real, fm_struct)

            CALL cp_cfm_to_fm(sr_env%integral, mtargetr=integral_real)
            CALL cp_fm_scale_and_add(1.0_dp, rho_ao_fm, 0.5_dp/pi, integral_real)

            CALL cp_fm_release(integral_real)

            DEALLOCATE (xnodes, zdata)

            stats%error = stats%error + sr_env%error*0.5_dp/pi
            stats%npoints = stats%npoints + npoints_total

            ! keep the cache
            IF (do_surface_green) THEN
               CALL green_functions_cache_reorder(g_surf_cache, sr_env%tnodes)
            END IF

            CALL simpsonrule_release(sr_env)
         END IF
      END DO

      CALL timestop(handle)
   END SUBROUTINE negf_add_rho_nonequiv

! **************************************************************************************************
!> \brief Reset integration statistics.
!> \param stats integration statistics
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE integration_status_reset(stats)
      TYPE(integration_status_type), INTENT(out)         :: stats

      stats%npoints = 0
      stats%error = 0.0_dp
   END SUBROUTINE integration_status_reset

! **************************************************************************************************
!> \brief Generate an integration method description string.
!> \param stats              integration statistics
!> \param integration_method integration method used
!> \return description string
!> \author Sergey Chulkov
! **************************************************************************************************
   PURE FUNCTION get_method_description_string(stats, integration_method) RESULT(method_descr)
      TYPE(integration_status_type), INTENT(in)          :: stats
      INTEGER, INTENT(in)                                :: integration_method
      CHARACTER(len=18)                                  :: method_descr

      CHARACTER(len=2)                                   :: method_abbr
      CHARACTER(len=6)                                   :: npoints_str

      SELECT CASE (integration_method)
      CASE (negfint_method_cc)
         ! Adaptive Clenshaw-Curtis method
         method_abbr = "CC"
      CASE (negfint_method_simpson)
         ! Adaptive Simpson's rule method
         method_abbr = "SR"
      CASE DEFAULT
         method_abbr = "??"
      END SELECT

      WRITE (npoints_str, '(I6)') stats%npoints
      WRITE (method_descr, '(A2,T4,A,T11,ES8.2E2)') method_abbr, TRIM(ADJUSTL(npoints_str)), stats%error
   END FUNCTION get_method_description_string

! **************************************************************************************************
!> \brief Compute electric current for one spin-channel through the scattering region.
!> \param contact_id1       reference contact
!> \param contact_id2       another contact
!> \param v_shift           shift in Hartree potential
!> \param negf_env          NEFG environment
!> \param negf_control      NEGF control
!> \param sub_env           NEGF parallel (sub)group environment
!> \param ispin             spin conponent to proceed
!> \param blacs_env_global  global BLACS environment
!> \return electric current in Amper
!> \author Sergey Chulkov
! **************************************************************************************************
   FUNCTION negf_compute_current(contact_id1, contact_id2, v_shift, negf_env, negf_control, sub_env, ispin, &
                                 blacs_env_global) RESULT(current)
      INTEGER, INTENT(in)                                :: contact_id1, contact_id2
      REAL(kind=dp), INTENT(in)                          :: v_shift
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      INTEGER, INTENT(in)                                :: ispin
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      REAL(kind=dp)                                      :: current

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_compute_current', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER :: threshold = 16.0_dp*EPSILON(0.0_dp)

      COMPLEX(kind=dp)                                   :: fermi_contact1, fermi_contact2, &
                                                            integr_lbound, integr_ubound
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: transm_coeff, xnodes
      COMPLEX(kind=dp), DIMENSION(1, 1)                  :: transmission
      INTEGER                                            :: handle, icontact, ipoint, max_points, &
                                                            min_points, ncontacts, npoints, &
                                                            npoints_total
      REAL(kind=dp) :: conv_density, energy, eta, ln_conv_density, mu_contact1, mu_contact2, &
         temperature_contact1, temperature_contact2, v_contact1, v_contact2
      TYPE(cp_cfm_p_type), ALLOCATABLE, DIMENSION(:)     :: zdata
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_single
      TYPE(cp_fm_type), POINTER                          :: weights
      TYPE(green_functions_cache_type)                   :: g_surf_cache
      TYPE(simpsonrule_type)                             :: sr_env

      current = 0.0_dp
      ! nothing to do
      IF (.NOT. ASSOCIATED(negf_env%s_s)) RETURN

      CALL timeset(routineN, handle)

      ncontacts = SIZE(negf_env%contacts)
      CPASSERT(contact_id1 <= ncontacts)
      CPASSERT(contact_id2 <= ncontacts)
      CPASSERT(contact_id1 /= contact_id2)

      v_contact1 = negf_control%contacts(contact_id1)%v_external
      mu_contact1 = negf_control%contacts(contact_id1)%fermi_level + v_contact1
      v_contact2 = negf_control%contacts(contact_id2)%v_external
      mu_contact2 = negf_control%contacts(contact_id2)%fermi_level + v_contact2

      IF (ABS(mu_contact1 - mu_contact2) < threshold) THEN
         CALL timestop(handle)
         RETURN
      END IF

      min_points = negf_control%integr_min_points
      max_points = negf_control%integr_max_points
      temperature_contact1 = negf_control%contacts(contact_id1)%temperature
      temperature_contact2 = negf_control%contacts(contact_id2)%temperature
      eta = negf_control%eta
      conv_density = negf_control%conv_density
      ln_conv_density = LOG(conv_density)

      integr_lbound = CMPLX(MIN(mu_contact1 + ln_conv_density*temperature_contact1, &
                                mu_contact2 + ln_conv_density*temperature_contact2), eta, kind=dp)
      integr_ubound = CMPLX(MAX(mu_contact1 - ln_conv_density*temperature_contact1, &
                                mu_contact2 - ln_conv_density*temperature_contact2), eta, kind=dp)

      npoints_total = 0
      npoints = min_points

      NULLIFY (fm_struct_single, weights)
      CALL cp_fm_struct_create(fm_struct_single, nrow_global=1, ncol_global=1, context=blacs_env_global)
      CALL cp_fm_create(weights, fm_struct_single)
      CALL cp_fm_set_all(weights, 1.0_dp)

      ALLOCATE (transm_coeff(npoints), xnodes(npoints), zdata(npoints))

      CALL simpsonrule_init(sr_env, xnodes, npoints, integr_lbound, integr_ubound, &
                            sr_shape_linear, negf_control%conv_density, weights)

      DO WHILE (npoints > 0 .AND. npoints_total < max_points)
         CALL green_functions_cache_expand(g_surf_cache, ncontacts, npoints)

         DO icontact = 1, ncontacts
            CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, 1:npoints), &
                                                   omega=xnodes(1:npoints), &
                                                   h0=negf_env%contacts(icontact)%h_00(ispin)%matrix, &
                                                   s0=negf_env%contacts(icontact)%s_00, &
                                                   h1=negf_env%contacts(icontact)%h_01(ispin)%matrix, &
                                                   s1=negf_env%contacts(icontact)%s_01, &
                                                   sub_env=sub_env, &
                                                   v_external=negf_control%contacts(icontact)%v_external, &
                                                   conv=negf_control%conv_green, transp=.FALSE.)
         END DO

         CALL negf_retarded_green_function_batch(omega=xnodes(1:npoints), &
                                                 v_shift=v_shift, &
                                                 ignore_bias=.FALSE., &
                                                 negf_env=negf_env, &
                                                 negf_control=negf_control, &
                                                 sub_env=sub_env, &
                                                 ispin=ispin, &
                                                 g_surf_contacts=g_surf_cache%g_surf_contacts(:, 1:npoints), &
                                                 transm_coeff=transm_coeff(1:npoints), &
                                                 transm_contact1=contact_id1, &
                                                 transm_contact2=contact_id2)

         DO ipoint = 1, npoints
            CALL cp_cfm_create(zdata(ipoint)%matrix, fm_struct_single)

            energy = REAL(xnodes(ipoint), kind=dp)
            fermi_contact1 = fermi_function(CMPLX(energy - mu_contact1, 0.0_dp, kind=dp), temperature_contact1)
            fermi_contact2 = fermi_function(CMPLX(energy - mu_contact2, 0.0_dp, kind=dp), temperature_contact2)

            transmission(1, 1) = transm_coeff(ipoint)*(fermi_contact1 - fermi_contact2)
            CALL cp_cfm_set_submatrix(zdata(ipoint)%matrix, transmission)
         END DO

         CALL green_functions_cache_release(g_surf_cache)

         npoints_total = npoints_total + npoints

         CALL simpsonrule_refine_integral(sr_env, zdata(1:npoints))

         IF (sr_env%error <= negf_control%conv_density) EXIT

         npoints = max_points - npoints_total
         IF (npoints <= 0) EXIT
         IF (npoints > SIZE(xnodes)) npoints = SIZE(xnodes)

         CALL simpsonrule_get_next_nodes(sr_env, xnodes, npoints)
      END DO

      CALL cp_cfm_get_submatrix(sr_env%integral, transmission)

      current = 0.5_dp/pi*REAL(transmission(1, 1), kind=dp)*e_charge/seconds

      CALL cp_fm_release(weights)
      CALL cp_fm_struct_release(fm_struct_single)

      CALL simpsonrule_release(sr_env)
      DEALLOCATE (transm_coeff, xnodes, zdata)

      CALL timestop(handle)
   END FUNCTION negf_compute_current

! **************************************************************************************************
!> \brief Print the Density of States.
!> \param log_unit     output unit
!> \param energy_min   energy point to start with
!> \param energy_max   energy point to end with
!> \param npoints      number of points to compute
!> \param v_shift      shift in Hartree potential
!> \param negf_env     NEFG environment
!> \param negf_control NEGF control
!> \param sub_env      NEGF parallel (sub)group environment
!> \param base_contact index of the reference contact
!> \param just_contact compute DOS for the given contact rather than for a scattering region
!> \param volume       unit cell volume
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_print_dos(log_unit, energy_min, energy_max, npoints, v_shift, negf_env, &
                             negf_control, sub_env, base_contact, just_contact, volume)
      INTEGER, INTENT(in)                                :: log_unit
      REAL(kind=dp), INTENT(in)                          :: energy_min, energy_max
      INTEGER, INTENT(in)                                :: npoints
      REAL(kind=dp), INTENT(in)                          :: v_shift
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      INTEGER, INTENT(in)                                :: base_contact
      INTEGER, INTENT(in), OPTIONAL                      :: just_contact
      REAL(kind=dp), INTENT(in), OPTIONAL                :: volume

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_print_dos', routineP = moduleN//':'//routineN

      CHARACTER                                          :: uks_str
      CHARACTER(len=15)                                  :: units_str
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: xnodes
      INTEGER                                            :: handle, icontact, ipoint, ispin, &
                                                            ncontacts, npoints_bundle, &
                                                            npoints_remain, nspins
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: dos
      TYPE(green_functions_cache_type)                   :: g_surf_cache

      CALL timeset(routineN, handle)

      IF (PRESENT(just_contact)) THEN
         nspins = SIZE(negf_env%contacts(just_contact)%h_00)
      ELSE
         nspins = SIZE(negf_env%h_s)
      END IF

      IF (log_unit > 0) THEN
         IF (PRESENT(volume)) THEN
            units_str = ' (angstroms^-3)'
         ELSE
            units_str = ''
         END IF

         IF (nspins > 1) THEN
            ! [alpha , beta]
            uks_str = ','
         ELSE
            ! [alpha + beta]
            uks_str = '+'
         END IF

         IF (PRESENT(just_contact)) THEN
            WRITE (log_unit, '(3A,T70,I11)') "# Density of states", TRIM(units_str), " for the contact No. ", just_contact
         ELSE
            WRITE (log_unit, '(3A)') "# Density of states", TRIM(units_str), " for the scattering region"
         END IF

         WRITE (log_unit, '(A,T10,A,T43,3A)') "#", "Energy (a.u.)", "Number of states [alpha ", uks_str, " beta]"

         WRITE (log_unit, '("#", T3,78("-"))')
      END IF

      ncontacts = SIZE(negf_env%contacts)
      CPASSERT(base_contact <= ncontacts)
      IF (PRESENT(just_contact)) THEN
         ncontacts = 2
         CPASSERT(just_contact == base_contact)
      END IF
      MARK_USED(base_contact)

      npoints_bundle = 4*sub_env%ngroups
      IF (npoints_bundle > npoints) npoints_bundle = npoints

      ALLOCATE (dos(npoints_bundle, nspins), xnodes(npoints_bundle))

      npoints_remain = npoints
      DO WHILE (npoints_remain > 0)
         IF (npoints_bundle > npoints_remain) npoints_bundle = npoints_remain

         IF (npoints > 1) THEN
            DO ipoint = 1, npoints_bundle
               xnodes(ipoint) = CMPLX(energy_min + REAL(npoints - npoints_remain + ipoint - 1, kind=dp)/ &
                                      REAL(npoints - 1, kind=dp)*(energy_max - energy_min), negf_control%eta, kind=dp)
            END DO
         ELSE
            xnodes(ipoint) = CMPLX(energy_min, negf_control%eta, kind=dp)
         END IF

         DO ispin = 1, nspins
            CALL green_functions_cache_expand(g_surf_cache, ncontacts, npoints_bundle)

            IF (PRESENT(just_contact)) THEN
               DO icontact = 1, ncontacts
                  CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, :), &
                                                         omega=xnodes(1:npoints_bundle), &
                                                         h0=negf_env%contacts(just_contact)%h_00(ispin)%matrix, &
                                                         s0=negf_env%contacts(just_contact)%s_00, &
                                                         h1=negf_env%contacts(just_contact)%h_01(ispin)%matrix, &
                                                         s1=negf_env%contacts(just_contact)%s_01, &
                                                         sub_env=sub_env, v_external=0.0_dp, &
                                                         conv=negf_control%conv_green, transp=(icontact == 1))
               END DO
            ELSE
               DO icontact = 1, ncontacts
                  CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, :), &
                                                         omega=xnodes(1:npoints_bundle), &
                                                         h0=negf_env%contacts(icontact)%h_00(ispin)%matrix, &
                                                         s0=negf_env%contacts(icontact)%s_00, &
                                                         h1=negf_env%contacts(icontact)%h_01(ispin)%matrix, &
                                                         s1=negf_env%contacts(icontact)%s_01, &
                                                         sub_env=sub_env, &
                                                         v_external=negf_control%contacts(icontact)%v_external, &
                                                         conv=negf_control%conv_green, transp=.FALSE.)
               END DO
            END IF

            CALL negf_retarded_green_function_batch(omega=xnodes(1:npoints_bundle), &
                                                    v_shift=v_shift, &
                                                    ignore_bias=.FALSE., &
                                                    negf_env=negf_env, &
                                                    negf_control=negf_control, &
                                                    sub_env=sub_env, &
                                                    ispin=ispin, &
                                                    g_surf_contacts=g_surf_cache%g_surf_contacts, &
                                                    dos=dos(1:npoints_bundle, ispin), &
                                                    just_contact=just_contact)

            CALL green_functions_cache_release(g_surf_cache)
         END DO

         IF (log_unit > 0) THEN
            DO ipoint = 1, npoints_bundle
               IF (nspins > 1) THEN
                  ! spin-polarised calculations: print alpha- and beta-spin components separately
                  WRITE (log_unit, '(T2,F20.8,T30,2ES25.11E3)') REAL(xnodes(ipoint), kind=dp), dos(ipoint, 1), dos(ipoint, 2)
               ELSE
                  ! spin-restricted calculations: print alpha- and beta-spin components together
                  WRITE (log_unit, '(T2,F20.8,T43,ES25.11E3)') REAL(xnodes(ipoint), kind=dp), 2.0_dp*dos(ipoint, 1)
               END IF
            END DO
         END IF

         npoints_remain = npoints_remain - npoints_bundle
      END DO

      DEALLOCATE (dos, xnodes)
      CALL timestop(handle)
   END SUBROUTINE negf_print_dos

! **************************************************************************************************
!> \brief Print the transmission coefficient.
!> \param log_unit     output unit
!> \param energy_min   energy point to start with
!> \param energy_max   energy point to end with
!> \param npoints      number of points to compute
!> \param v_shift      shift in Hartree potential
!> \param negf_env     NEFG environment
!> \param negf_control NEGF control
!> \param sub_env      NEGF parallel (sub)group environment
!> \param contact_id1  index of a reference contact
!> \param contact_id2  index of another contact
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE negf_print_transmission(log_unit, energy_min, energy_max, npoints, v_shift, negf_env, &
                                      negf_control, sub_env, contact_id1, contact_id2)
      INTEGER, INTENT(in)                                :: log_unit
      REAL(kind=dp), INTENT(in)                          :: energy_min, energy_max
      INTEGER, INTENT(in)                                :: npoints
      REAL(kind=dp), INTENT(in)                          :: v_shift
      TYPE(negf_env_type), INTENT(in)                    :: negf_env
      TYPE(negf_control_type), POINTER                   :: negf_control
      TYPE(negf_subgroup_env_type), INTENT(in)           :: sub_env
      INTEGER, INTENT(in)                                :: contact_id1, contact_id2

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_print_transmission', &
         routineP = moduleN//':'//routineN

      CHARACTER                                          :: uks_str
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:)        :: xnodes
      COMPLEX(kind=dp), ALLOCATABLE, DIMENSION(:, :)     :: transm_coeff
      INTEGER                                            :: handle, icontact, ipoint, ispin, &
                                                            ncontacts, npoints_bundle, &
                                                            npoints_remain, nspins
      REAL(kind=dp)                                      :: rscale
      TYPE(green_functions_cache_type)                   :: g_surf_cache

      CALL timeset(routineN, handle)

      nspins = SIZE(negf_env%h_s)

      IF (nspins > 1) THEN
         ! [alpha , beta]
         uks_str = ','
      ELSE
         ! [alpha + beta]
         uks_str = '+'
      END IF

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(A)') "# Transmission coefficient (G0 = 2 e^2/h) for the scattering region"

         WRITE (log_unit, '(A,T10,A,T39,3A)') "#", "Energy (a.u.)", "Transmission coefficient [alpha ", uks_str, " beta]"
         WRITE (log_unit, '("#", T3,78("-"))')
      END IF

      ncontacts = SIZE(negf_env%contacts)
      CPASSERT(contact_id1 <= ncontacts)
      CPASSERT(contact_id2 <= ncontacts)

      IF (nspins == 1) THEN
         rscale = 2.0_dp
      ELSE
         rscale = 1.0_dp
      END IF

      ! print transmission coefficients in terms of G0 = 2 * e^2 / h = 1 / pi ;
      ! transmission coefficients returned by negf_retarded_green_function_batch() are already multiplied by 2 / pi
      rscale = 0.5_dp*rscale

      npoints_bundle = 4*sub_env%ngroups
      IF (npoints_bundle > npoints) npoints_bundle = npoints

      ALLOCATE (transm_coeff(npoints_bundle, nspins), xnodes(npoints_bundle))

      npoints_remain = npoints
      DO WHILE (npoints_remain > 0)
         IF (npoints_bundle > npoints_remain) npoints_bundle = npoints_remain

         IF (npoints > 1) THEN
            DO ipoint = 1, npoints_bundle
               xnodes(ipoint) = CMPLX(energy_min + REAL(npoints - npoints_remain + ipoint - 1, kind=dp)/ &
                                      REAL(npoints - 1, kind=dp)*(energy_max - energy_min), negf_control%eta, kind=dp)
            END DO
         ELSE
            xnodes(ipoint) = CMPLX(energy_min, negf_control%eta, kind=dp)
         END IF

         DO ispin = 1, nspins
            CALL green_functions_cache_expand(g_surf_cache, ncontacts, npoints_bundle)

            DO icontact = 1, ncontacts
               CALL negf_surface_green_function_batch(g_surf=g_surf_cache%g_surf_contacts(icontact, :), &
                                                      omega=xnodes(1:npoints_bundle), &
                                                      h0=negf_env%contacts(icontact)%h_00(ispin)%matrix, &
                                                      s0=negf_env%contacts(icontact)%s_00, &
                                                      h1=negf_env%contacts(icontact)%h_01(ispin)%matrix, &
                                                      s1=negf_env%contacts(icontact)%s_01, &
                                                      sub_env=sub_env, &
                                                      v_external=negf_control%contacts(icontact)%v_external, &
                                                      conv=negf_control%conv_green, transp=.FALSE.)
            END DO

            CALL negf_retarded_green_function_batch(omega=xnodes(1:npoints_bundle), &
                                                    v_shift=v_shift, &
                                                    ignore_bias=.FALSE., &
                                                    negf_env=negf_env, &
                                                    negf_control=negf_control, &
                                                    sub_env=sub_env, &
                                                    ispin=ispin, &
                                                    g_surf_contacts=g_surf_cache%g_surf_contacts, &
                                                    transm_coeff=transm_coeff(1:npoints_bundle, ispin), &
                                                    transm_contact1=contact_id1, &
                                                    transm_contact2=contact_id2)

            CALL green_functions_cache_release(g_surf_cache)
         END DO

         IF (log_unit > 0) THEN
            DO ipoint = 1, npoints_bundle
               IF (nspins > 1) THEN
                  ! spin-polarised calculations: print alpha- and beta-spin components separately
                  WRITE (log_unit, '(T2,F20.8,T30,2ES25.11E3)') &
                     REAL(xnodes(ipoint), kind=dp), rscale*REAL(transm_coeff(ipoint, 1:2), kind=dp)
               ELSE
                  ! spin-restricted calculations: print alpha- and beta-spin components together
                  WRITE (log_unit, '(T2,F20.8,T43,ES25.11E3)') &
                     REAL(xnodes(ipoint), kind=dp), rscale*REAL(transm_coeff(ipoint, 1), kind=dp)
               END IF
            END DO
         END IF

         npoints_remain = npoints_remain - npoints_bundle
      END DO

      DEALLOCATE (transm_coeff, xnodes)
      CALL timestop(handle)
   END SUBROUTINE negf_print_transmission
END MODULE negf_methods
