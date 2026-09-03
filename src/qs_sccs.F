!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Self-consistent continuum solvation (SCCS) model implementation
!> \author  Matthias Krack (MK)
!> \version 1.0
!> \par Literature:
!>      - J.-L. Fattebert and F. Gygi,
!>        Density functional theory for efficient ab initio molecular dynamics
!>        simulations in solution, J. Comput. Chem. 23, 662-666 (2002)
!>      - O. Andreussi, I. Dabo, and N. Marzari,
!>        Revised self-consistent continuum solvation in electronic-structure
!>        calculations, J. Chem. Phys. 136, 064102-20 (2012)
!> \par History:
!>      - Creation (10.10.2013,MK)
!>      - Derivatives using finite differences (26.11.2013,MK)
!>      - Cube file dump of the dielectric function (19.12.2013,MK)
!>      - Cube file dump of the polarisation potential (20.12.2013,MK)
!>      - Calculation of volume and surface of the cavity (21.12.2013,MK)
!>      - Functional derivative of the cavitation energy (28.12.2013,MK)
!>      - Update printout (11.11.2022,MK)
! **************************************************************************************************

MODULE qs_sccs

   USE aics_interface,                  ONLY: &
        clear_aics_field_registrations, cp2k_aics_cancel_initialize_from_fortran, &
        cp2k_aics_get_last_error_from_fortran, cp2k_aics_initialize_from_fortran, &
        cp2k_aics_initialize_prepared_from_fortran, cp2k_aics_set_electrostatic_potential, &
        cp2k_aics_set_fd_field_squared_sensitivity_parallel, &
        cp2k_aics_set_fd_field_squared_sensitivity_perpendicular, cp2k_aics_set_grid_metadata, &
        cp2k_aics_set_relative_permittivity_parallel, &
        cp2k_aics_set_relative_permittivity_perpendicular, cp2k_aics_set_solute_charge_density, &
        cp2k_aics_solve_from_fortran
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type,&
                                              sccs_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              low_print_level
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_realspace_grid_init,          ONLY: init_input_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_constants,                 ONLY: &
        aics_interface_plane_ab, aics_interface_plane_bc, aics_interface_plane_ca, &
        aics_poisson_solver_finite_difference, aics_poisson_solver_finite_element, sccs_andreussi, &
        sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7, sccs_derivative_fft, &
        sccs_fattebert_gygi, sccs_saa_andreussi
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE iso_c_binding,                   ONLY: c_double,&
                                              c_int,&
                                              c_loc,&
                                              c_null_ptr,&
                                              c_size_t
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_integral_ab,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_func_u_convolution,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_analytic,&
                                              pw_poisson_mt,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_type,&
                                              realspace_grid_input_type,&
                                              realspace_grid_type,&
                                              rs_grid_create,&
                                              rs_grid_create_descriptor,&
                                              rs_grid_release,&
                                              rs_grid_release_descriptor
   USE rs_methods,                      ONLY: derive_fdm_cd3,&
                                              derive_fdm_cd5,&
                                              derive_fdm_cd7
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_sccs'

   TYPE, PRIVATE :: aics_interface_geometry_type
      REAL(KIND=dp) :: upper_interface_coordinate
      REAL(KIND=dp) :: lower_interface_coordinate
      REAL(KIND=dp) :: unit_normal(3)
      INTEGER :: normal_axis
      INTEGER :: normal_sign
      INTEGER :: n_upper_surface_atoms
      INTEGER :: n_lower_surface_atoms
      LOGICAL, ALLOCATABLE :: upper_surface_atom(:)
      LOGICAL, ALLOCATABLE :: lower_surface_atom(:)
   END TYPE aics_interface_geometry_type

   PUBLIC :: get_aics_interface_normal, print_sccs_results, sccs, sccs_aics

CONTAINS

! **************************************************************************************************
!> \brief   Return the Cartesian unit normal for the selected AICS interface plane.
!> \param   cell              CP2K simulation cell.
!> \param   interface_plane   Selected AB, BC, or CA interface plane.
!> \param   unit_normal       Cartesian unit normal.
! **************************************************************************************************
   SUBROUTINE get_aics_interface_normal(cell, interface_plane, unit_normal)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      INTEGER, INTENT(IN)                                :: interface_plane
      REAL(KIND=dp), DIMENSION(3), INTENT(OUT)           :: unit_normal

      TYPE(aics_interface_geometry_type)                 :: interface_geometry

      CALL determine_aics_orientation(cell, interface_plane, interface_geometry)
      unit_normal = interface_geometry%unit_normal
   END SUBROUTINE get_aics_interface_normal

! **************************************************************************************************
!> \brief   Self-consistent continuum solvation (SCCS) model implementation
!> \param qs_env ...
!> \param rho_tot_gspace ...
!> \param v_hartree_gspace ...
!> \param v_sccs ...
!> \param h_stress ...
!> \par History:
!>      - Creation (10.10.2013,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! **************************************************************************************************

   SUBROUTINE sccs(qs_env, rho_tot_gspace, v_hartree_gspace, v_sccs, h_stress)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_tot_gspace, v_hartree_gspace
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_sccs
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: h_stress

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sccs'
      REAL(KIND=dp), PARAMETER                           :: epstol = 1.0E-8_dp

      CHARACTER(LEN=4*default_string_length)             :: message, my_pos_cube
      CHARACTER(LEN=default_path_length)                 :: cube_path, filename, mpi_filename, &
                                                            print_path
      INTEGER                                            :: cube_unit, handle, i, ispin, iter, j, k, &
                                                            nspin, output_unit, print_level
      INTEGER(KIND=int_8)                                :: ngpts
      INTEGER, DIMENSION(3)                              :: lb, ub
      LOGICAL                                            :: append_cube, calculate_stress_tensor, &
                                                            do_kpoints, mpi_io, should_output
      REAL(KIND=dp) :: alpha_zeta, cavity_surface, cavity_volume, cell_volume, delta_eta, &
         delta_zeta, dphi2, dvol, e_tot, epsilon_solvent, f, f0, polarisation_charge, r, R_solv, &
         rho_delta, rho_delta_avg, rho_delta_max, rho_iter_new, tot_rho_elec, tot_rho_solute
      REAL(KIND=dp), DIMENSION(3)                        :: abc, lxyz, uxyz
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                               :: d_s_rhoel, deps_elec, dtf, eps_elec, ff, &
                                                            p_e_interface, s, t, tem_convolution, &
                                                            tem_func, u
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: dln_eps_elec, dphi_tot
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_pw_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_pw_r_sccs
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool)
      NULLIFY (cp_subsys)
      NULLIFY (dft_control)
      NULLIFY (energy)
      NULLIFY (input)
      NULLIFY (logger)
      NULLIFY (para_env)
      NULLIFY (particles)
      NULLIFY (poisson_env)
      NULLIFY (pw_env)
      NULLIFY (pw_pools)
      NULLIFY (rho)
      NULLIFY (sccs_control)
      NULLIFY (scf_env)

      ! Load data from Quickstep environment
      CALL get_qs_env(qs_env=qs_env, &
                      cp_subsys=cp_subsys, &
                      do_kpoints=do_kpoints, &
                      dft_control=dft_control, &
                      energy=energy, &
                      input=input, &
                      para_env=para_env, &
                      pw_env=pw_env, &
                      rho=rho, &
                      scf_env=scf_env)
      CALL cp_subsys_get(cp_subsys, particles=particles)

      sccs_control => dft_control%sccs_control

      CPASSERT(ASSOCIATED(qs_env))

      IF (do_kpoints) THEN
         CPWARN("SCCS with k-points has not yet been fully validated")
      END IF

      IF (PRESENT(h_stress)) THEN
         calculate_stress_tensor = .TRUE.
         h_stress(:, :) = 0.0_dp
         CPWARN("The stress tensor for SCCS has not yet been fully validated")
      ELSE
         calculate_stress_tensor = .FALSE.
      END IF

      ! Get access to the PW grid pool
      CALL pw_env_get(pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools, &
                      poisson_env=poisson_env)

      CALL pw_zero(v_sccs)

      ! Calculate no SCCS contribution, if the requested SCF convergence threshold is not reached yet
      IF (.NOT. sccs_control%sccs_activated) THEN
         IF (sccs_control%eps_scf > 0.0_dp) THEN
            IF ((scf_env%iter_delta > sccs_control%eps_scf) .OR. &
                ((qs_env%scf_env%outer_scf%iter_count == 0) .AND. &
                 (qs_env%scf_env%iter_count <= 1))) THEN
               IF (calculate_stress_tensor) THEN
                  ! Request also the calculation of the stress tensor contribution
                  CALL pw_poisson_solve(poisson_env=poisson_env, &
                                        density=rho_tot_gspace, &
                                        ehartree=energy%hartree, &
                                        vhartree=v_hartree_gspace, &
                                        h_stress=h_stress)
               ELSE
                  CALL pw_poisson_solve(poisson_env=poisson_env, &
                                        density=rho_tot_gspace, &
                                        ehartree=energy%hartree, &
                                        vhartree=v_hartree_gspace)
               END IF
               energy%sccs_pol = 0.0_dp
               energy%sccs_cav = 0.0_dp
               energy%sccs_dis = 0.0_dp
               energy%sccs_rep = 0.0_dp
               energy%sccs_sol = 0.0_dp
               energy%sccs_hartree = energy%hartree
               CALL timestop(handle)
               RETURN
            END IF
         END IF
         sccs_control%sccs_activated = .TRUE.
      END IF

      nspin = dft_control%nspins

      ! Manage print output control
      logger => cp_get_default_logger()
      print_level = logger%iter_info%print_level
      print_path = "DFT%PRINT%SCCS"
      should_output = (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                        TRIM(print_path)), cp_p_file))
      output_unit = cp_print_key_unit_nr(logger, input, TRIM(print_path), &
                                         extension=".sccs", &
                                         ignore_should_output=should_output, &
                                         log_filename=.FALSE.)

      ! Get rho
      CALL qs_rho_get(rho_struct=rho, &
                      rho_r=rho_pw_r, &
                      rho_r_sccs=rho_pw_r_sccs)

      ! Retrieve the last rho_iter from the previous SCCS cycle if available
      CPASSERT(ASSOCIATED(rho_pw_r_sccs))

      ! Retrieve the total electronic density in r-space
      BLOCK
         TYPE(pw_r3d_rs_type) :: rho_elec
         CALL auxbas_pw_pool%create_pw(rho_elec)

         ! Retrieve grid parameters
         ngpts = rho_elec%pw_grid%ngpts
         dvol = rho_elec%pw_grid%dvol
         cell_volume = rho_elec%pw_grid%vol
         abc(1:3) = REAL(rho_elec%pw_grid%npts(1:3), KIND=dp)*rho_elec%pw_grid%dr(1:3)
         lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
         ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

         CALL pw_copy(rho_pw_r(1), rho_elec)
         DO ispin = 2, nspin
            CALL pw_axpy(rho_pw_r(ispin), rho_elec)
         END DO
         tot_rho_elec = pw_integrate_function(rho_elec)

         ! Calculate the dielectric (smoothed) function of rho_elec in r-space
         CALL auxbas_pw_pool%create_pw(eps_elec)
         CALL auxbas_pw_pool%create_pw(deps_elec)
         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(s)
            CALL auxbas_pw_pool%create_pw(d_s_rhoel)
            CALL auxbas_pw_pool%create_pw(u)
            CALL auxbas_pw_pool%create_pw(ff)
            CALL auxbas_pw_pool%create_pw(t)
            CALL auxbas_pw_pool%create_pw(dtf)
            f0 = sccs_control%f0
            delta_eta = sccs_control%delta_eta
            alpha_zeta = sccs_control%alpha_zeta
            R_solv = sccs_control%R_solv
            delta_zeta = sccs_control%delta_zeta
         END IF

         ! Relative permittivity or dielectric constant of the solvent (medium)
         epsilon_solvent = sccs_control%epsilon_solvent
         SELECT CASE (sccs_control%method_id)
         CASE (sccs_andreussi)
            CALL andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%rho_max, &
                           sccs_control%rho_min)
         CASE (sccs_fattebert_gygi)
            CALL fattebert_gygi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%beta, &
                                sccs_control%rho_zero)
         CASE (sccs_saa_andreussi)
            CALL sa_andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%rho_max, &
                              sccs_control%rho_min, s, d_s_rhoel)
         CASE DEFAULT
            CPABORT("Invalid method specified for SCCS model")
         END SELECT

         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN
            ! Optional output of the dielectric function in cube file format
            filename = "DIELECTRIC_FUNCTION"
            cube_path = TRIM(print_path)//"%"//TRIM(filename)
            IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                      cp_p_file)) THEN
               append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
               my_pos_cube = "REWIND"
               IF (append_cube) my_pos_cube = "APPEND"
               mpi_io = .TRUE.
               cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                extension=".cube", middle_name=TRIM(filename), &
                                                file_position=my_pos_cube, log_filename=.FALSE., &
                                                mpi_io=mpi_io, fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=cube_unit, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(T3,A)") &
                     "SCCS| The dielectric function is written in cube file format to the file:", &
                     "SCCS| "//TRIM(filename)
               END IF
               CALL cp_pw_to_cube(eps_elec, cube_unit, TRIM(filename), particles=particles, &
                                  stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                  mpi_io=mpi_io)
               CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
            END IF
         END IF

         ! Calculate the (quantum) volume and surface of the solute cavity
         cavity_surface = 0.0_dp
         cavity_volume = 0.0_dp

         IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN

            BLOCK
               TYPE(pw_r3d_rs_type) :: theta, norm_drho_elec
               TYPE(pw_r3d_rs_type), DIMENSION(3)                        :: drho_elec
               CALL auxbas_pw_pool%create_pw(theta)
               CALL pw_zero(theta)

               ! Calculate the (quantum) volume of the solute cavity
               f = 1.0_dp/(epsilon_solvent - 1.0_dp)
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(epsilon_solvent,eps_elec,f,lb,theta,ub)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        theta%array(i, j, k) = f*(epsilon_solvent - eps_elec%array(i, j, k))
                     END DO
                  END DO
               END DO
!$OMP    END PARALLEL DO
               cavity_volume = pw_integrate_function(theta)

               ! Calculate the derivative of the electronic density in r-space
               ! TODO: Could be retrieved from the QS environment
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_elec(i))
               END DO
               CALL derive(rho_elec, drho_elec, sccs_derivative_fft, pw_env, input)

               CALL auxbas_pw_pool%create_pw(norm_drho_elec)

               ! Calculate the norm of the gradient of the electronic density in r-space
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(drho_elec,lb,norm_drho_elec,ub)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        norm_drho_elec%array(i, j, k) = SQRT(drho_elec(1)%array(i, j, k)* &
                                                             drho_elec(1)%array(i, j, k) + &
                                                             drho_elec(2)%array(i, j, k)* &
                                                             drho_elec(2)%array(i, j, k) + &
                                                             drho_elec(3)%array(i, j, k)* &
                                                             drho_elec(3)%array(i, j, k))
                     END DO
                  END DO
               END DO
!$OMP    END PARALLEL DO

               ! Optional output of the norm of the density gradient in cube file format
               filename = "DENSITY_GRADIENT"
               cube_path = TRIM(print_path)//"%"//TRIM(filename)
               IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                         cp_p_file)) THEN
                  append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
                  my_pos_cube = "REWIND"
                  IF (append_cube) my_pos_cube = "APPEND"
                  mpi_io = .TRUE.
                  cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                   extension=".cube", middle_name=TRIM(filename), &
                                                   file_position=my_pos_cube, log_filename=.FALSE., &
                                                   mpi_io=mpi_io, fout=mpi_filename)
                  IF (output_unit > 0) THEN
                     IF (.NOT. mpi_io) THEN
                        INQUIRE (UNIT=cube_unit, NAME=filename)
                     ELSE
                        filename = mpi_filename
                     END IF
                     WRITE (UNIT=output_unit, FMT="(T3,A)") &
                        "SCCS| The norm of the density gradient is written in cube file format to the file:", &
                        "SCCS| "//TRIM(filename)
                  END IF
                  CALL cp_pw_to_cube(norm_drho_elec, cube_unit, TRIM(filename), particles=particles, &
                                     stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                     mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
               END IF

               ! Calculate the (quantum) surface of the solute cavity
               SELECT CASE (sccs_control%method_id)
               CASE (sccs_andreussi)
                  CALL surface_andreussi(rho_elec, norm_drho_elec, theta, epsilon_solvent, &
                                         sccs_control%rho_max, sccs_control%rho_min, &
                                         sccs_control%delta_rho)
               CASE (sccs_fattebert_gygi)
                  CALL surface_fattebert_gygi(rho_elec, norm_drho_elec, theta, epsilon_solvent, &
                                              sccs_control%beta, sccs_control%rho_zero, &
                                              sccs_control%delta_rho)
               CASE (sccs_saa_andreussi)
               CASE DEFAULT
                  CPABORT("Invalid method specified for SCCS model")
               END SELECT
               cavity_surface = pw_integrate_function(theta)

               ! Release storage
               CALL auxbas_pw_pool%give_back_pw(theta)
               CALL auxbas_pw_pool%give_back_pw(norm_drho_elec)
               DO i = 1, 3
                  CALL auxbas_pw_pool%give_back_pw(drho_elec(i))
               END DO
            END BLOCK

         END IF ! epsilon_solvent > 1

         CALL auxbas_pw_pool%give_back_pw(rho_elec)
      END BLOCK

      BLOCK
         TYPE(pw_r3d_rs_type) :: rho_tot, phi_tot, rho_solute, rho_tot_zero
         ! Retrieve the total charge density (core + elec) of the solute in r-space
         CALL auxbas_pw_pool%create_pw(rho_solute)
         CALL pw_zero(rho_solute)
         CALL pw_transfer(rho_tot_gspace, rho_solute)
         tot_rho_solute = pw_integrate_function(rho_solute)

         ! Check total charge
         IF (ABS(tot_rho_solute) >= 1.0E-6_dp) THEN
            IF ((poisson_env%parameters%solver /= pw_poisson_analytic) .AND. &
                (poisson_env%parameters%solver /= pw_poisson_mt)) THEN
               WRITE (UNIT=message, FMT="(A,SP,F0.6,A)") &
                  "The system (solute) has a non-negligible charge of ", -tot_rho_solute, &
                  ". It is recommended to use non-periodic boundary conditions (PERIODIC none) "// &
                  "combined with an appropriate Poisson solver (POISSON_SOLVER MT or analytic)"
               CPWARN(message)
            END IF
         END IF

         ! Reassign work storage to rho_tot_zero, because rho_elec is no longer needed
         CALL auxbas_pw_pool%create_pw(rho_tot_zero)

         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN

            ! Build the initial (rho_iter = 0) total charge density (solute plus polarisation) in r-space
            ! eps_elec <- ln(eps_elec)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(eps_elec,lb,message,output_unit,para_env,ub) &
!$OMP             SHARED(rho_solute,rho_tot_zero)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (eps_elec%array(i, j, k) < 1.0_dp) THEN
                        WRITE (UNIT=message, FMT="(A,ES12.3,A,3(I0,A))") &
                           "SCCS| Invalid dielectric function value ", eps_elec%array(i, j, k), &
                           " encountered at grid point (", i, ",", j, ",", k, ")"
                        CPABORT(message)
                     END IF
                     rho_tot_zero%array(i, j, k) = rho_solute%array(i, j, k)/eps_elec%array(i, j, k)
                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         ELSE

            lxyz(:) = REAL(u%pw_grid%bounds(1, :), dp)
            uxyz(:) = REAL(u%pw_grid%bounds(2, :), dp)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,r) &
!$OMP             SHARED(lb,message,output_unit,para_env,ub,dtf,uxyz,lxyz) &
!$OMP             SHARED(u, alpha_zeta, R_solv, delta_zeta, eps_elec)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (eps_elec%array(i, j, k) < 1.0_dp) THEN
                        WRITE (UNIT=message, FMT="(A,ES12.3,A,3(I0,A))") &
                           "SCCS| Invalid dielectric function value ", eps_elec%array(i, j, k), &
                           " encountered at grid point (", i, ",", j, ",", k, ")"
                        CPABORT(message)
                     END IF
                     IF (i > -1) r = ((REAL(uxyz(1), dp) - REAL(i, dp) + 0.5_dp)*u%pw_grid%dr(1))**2.0_dp
                     IF (i < 0) r = ((REAL(i, dp) - REAL(lxyz(1), dp) + 0.5_dp)*u%pw_grid%dr(1))**2.0_dp
                     IF (j > -1) r = ((REAL(uxyz(2), dp) - REAL(j, dp) + 0.5_dp)*u%pw_grid%dr(2))**2.0_dp + r
                     IF (j < 0) r = ((REAL(j, dp) - REAL(lxyz(2), dp) + 0.5_dp)*u%pw_grid%dr(2))**2.0_dp + r
                     IF (k > -1) r = ((REAL(uxyz(3), dp) - REAL(k, dp) + 0.5_dp)*u%pw_grid%dr(3))**2.0_dp + r
                     IF (k < 0) r = ((REAL(k, dp) - REAL(lxyz(3), dp) + 0.5_dp)*u%pw_grid%dr(3))**2.0_dp + r
                     r = SQRT(r)
                     u%array(i, j, k) = (1.0_dp/2.0_dp)*erfc((r - alpha_zeta*R_solv)/delta_zeta)
                     dtf%array(i, j, k) = 1.0_dp
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL pw_func_u_convolution(poisson_env=poisson_env, func=dtf, convolution=ff, u=u)
            CALL pw_scale(u, 1.0_dp/ff%array(0, 0, 0))
            CALL pw_func_u_convolution(poisson_env=poisson_env, func=s, convolution=ff, u=u)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(lb,para_env,ub, delta_eta, f0, epsilon_solvent) &
!$OMP             SHARED(ff, t, dtf, eps_elec, rho_tot_zero, rho_solute, s)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     t%array(i, j, k) = 0.5_dp*(1.0_dp + erf((ff%array(i, j, k) - f0)/delta_eta))
                     dtf%array(i, j, k) = (1.0_dp/SQRT(pi))*EXP(-(ff%array(i, j, k) - f0)**2.0/(delta_eta**2.0))
                     eps_elec%array(i, j, k) = EXP(LOG(epsilon_solvent)*(1 - (s%array(i, j, k) &
                                                                              + (1 - s%array(i, j, k))*t%array(i, j, k))))
                     rho_tot_zero%array(i, j, k) = rho_solute%array(i, j, k)/eps_elec%array(i, j, k)
!                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            ! Optional output of the dielectric function in cube file format (solvent aware algorithm case)
            filename = "DIELECTRIC_FUNCTION"
            cube_path = TRIM(print_path)//"%"//TRIM(filename)
            IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                      cp_p_file)) THEN
               append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
               my_pos_cube = "REWIND"
               IF (append_cube) my_pos_cube = "APPEND"
               mpi_io = .TRUE.
               cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                extension=".cube", middle_name=TRIM(filename), &
                                                file_position=my_pos_cube, log_filename=.FALSE., &
                                                mpi_io=mpi_io, fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=cube_unit, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(T3,A)") &
                     "SCCS| The dielectric function is written in cube file format to the file:", &
                     "SCCS| "//TRIM(filename)
               END IF
               CALL cp_pw_to_cube(eps_elec, cube_unit, TRIM(filename), particles=particles, &
                                  stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                  mpi_io=mpi_io)
               CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
            END IF

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(lb,para_env,ub, delta_eta, f0, epsilon_solvent) &
!$OMP             SHARED(ff, t, dtf, eps_elec, rho_tot_zero, rho_solute, s)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         END IF

         ! Build the derivative of LOG(eps_elec)
         DO i = 1, 3
            CALL auxbas_pw_pool%create_pw(dln_eps_elec(i))
            CALL pw_zero(dln_eps_elec(i))
         END DO
         CALL derive(eps_elec, dln_eps_elec, sccs_control%derivative_method, pw_env, input)
         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%give_back_pw(eps_elec)
         END IF

         ! Print header for the SCCS cycle
         IF (should_output .AND. (output_unit > 0)) THEN
            IF (print_level > low_print_level) THEN
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
                  "SCCS| Total electronic charge density ", -tot_rho_elec, &
                  "SCCS| Total charge density (solute)   ", -tot_rho_solute
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.3)") &
                  "SCCS| Volume of the cell           [bohr^3]", cell_volume, &
                  "SCCS|                              [angstrom^3]", &
                  cp_unit_from_cp2k(cell_volume, "angstrom^3")
               IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.3)") &
                     "SCCS| Volume of the solute cavity  [bohr^3]", cavity_volume, &
                     "SCCS|                              [angstrom^3]", &
                     cp_unit_from_cp2k(cavity_volume, "angstrom^3"), &
                     "SCCS| Surface of the solute cavity [bohr^2]", cavity_surface, &
                     "SCCS|                              [angstrom^2]", &
                     cp_unit_from_cp2k(cavity_surface, "angstrom^2")
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS|", &
                  "SCCS|   Step    Average residual    Maximum residual         E(Hartree) [a.u.]"
            END IF
         END IF

         ! Get storage for the derivative of the total potential (field) in r-space
         DO i = 1, 3
            CALL auxbas_pw_pool%create_pw(dphi_tot(i))
         END DO

         ! Initialise the total charge density in r-space rho_tot with rho_tot_zero + rho_iter_zero
         CALL auxbas_pw_pool%create_pw(rho_tot)
         CALL pw_copy(rho_tot_zero, rho_tot)
         CALL pw_axpy(rho_pw_r_sccs, rho_tot)

         CALL auxbas_pw_pool%create_pw(phi_tot)
         CALL pw_zero(phi_tot)

         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(p_e_interface)
            CALL pw_zero(p_e_interface)
         END IF

         ! Main SCCS iteration loop
         iter = 0

         iter_loop: DO

            ! Increment iteration counter
            iter = iter + 1

            ! Check if the requested maximum number of SCCS iterations is reached
            IF (iter > sccs_control%max_iter) THEN
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,/,T3,A,I0,A)") &
                     "SCCS| Maximum number of SCCS iterations reached", &
                     "SCCS| Iteration cycle did not converge in ", sccs_control%max_iter, " steps"
               ELSE
                  WRITE (UNIT=message, FMT="(A,I0,A)") &
                     "The SCCS iteration cycle did not converge in ", sccs_control%max_iter, " steps"
                  CPWARN(message)
               END IF
               EXIT iter_loop
            END IF

            ! Calculate derivative of the current total potential in r-space
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_tot, &
                                  vhartree=phi_tot, &
                                  dvhartree=dphi_tot)
            energy%sccs_hartree = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)

            ! Update total charge density (solute plus polarisation) in r-space
            ! based on the iterated polarisation charge density
            f = 1.0_dp/fourpi
            rho_delta_avg = 0.0_dp
            rho_delta_max = 0.0_dp
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k,rho_delta,rho_iter_new) &
!$OMP                SHARED(dln_eps_elec,dphi_tot,f,lb,rho_pw_r_sccs,ub) &
!$OMP                SHARED(rho_tot,rho_tot_zero,sccs_control) &
!$OMP                REDUCTION(+:rho_delta_avg) &
!$OMP                REDUCTION(MAX:rho_delta_max)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     rho_iter_new = (dln_eps_elec(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                                     dln_eps_elec(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                                     dln_eps_elec(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k))*f
                     rho_iter_new = rho_pw_r_sccs%array(i, j, k) + &
                                    sccs_control%mixing*(rho_iter_new - rho_pw_r_sccs%array(i, j, k))
                     rho_delta = ABS(rho_iter_new - rho_pw_r_sccs%array(i, j, k))
                     rho_delta_max = MAX(rho_delta, rho_delta_max)
                     rho_delta_avg = rho_delta_avg + rho_delta
                     rho_tot%array(i, j, k) = rho_tot_zero%array(i, j, k) + rho_iter_new
                     rho_pw_r_sccs%array(i, j, k) = rho_iter_new
                  END DO
               END DO
            END DO
!$OMP    END PARALLEL DO

            CALL para_env%sum(rho_delta_avg)
            rho_delta_avg = rho_delta_avg/REAL(ngpts, KIND=dp)
            CALL para_env%max(rho_delta_max)

            IF (should_output .AND. (output_unit > 0)) THEN
               IF (print_level > low_print_level) THEN
                  IF ((ABS(rho_delta_avg) < 1.0E-8_dp) .OR. &
                      (ABS(rho_delta_avg) >= 1.0E5_dp)) THEN
                     WRITE (UNIT=output_unit, FMT="(T3,A,I6,4X,ES16.4,4X,ES16.4,1X,F25.12)") &
                        "SCCS| ", iter, rho_delta_avg, rho_delta_max, energy%sccs_hartree
                  ELSE
                     WRITE (UNIT=output_unit, FMT="(T3,A,I6,4X,F16.8,4X,F16.8,1X,F25.12)") &
                        "SCCS| ", iter, rho_delta_avg, rho_delta_max, energy%sccs_hartree
                  END IF
               END IF
            END IF

            ! Check if the SCCS iteration cycle is converged to the requested tolerance
            IF (rho_delta_max <= sccs_control%eps_sccs) THEN
               IF (should_output .AND. (output_unit > 0)) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,I0,A)") &
                     "SCCS| Iteration cycle converged in ", iter, " steps"
               END IF
               EXIT iter_loop
            END IF

         END DO iter_loop

         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(tem_convolution)
            CALL auxbas_pw_pool%create_pw(tem_func)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k,dphi2) &
!$OMP                SHARED(lb,ub, tem_func, s, dtf) &
!$OMP                SHARED(dphi_tot, p_e_interface, eps_elec, epsilon_solvent)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     p_e_interface%array(i, j, k) = 1.0_dp/(4.0_dp*twopi)*(EXP(eps_elec%array(i, j, k))*LOG(epsilon_solvent))
                     dphi2 = dphi_tot(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                             dphi_tot(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                             dphi_tot(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k)
                     p_e_interface%array(i, j, k) = p_e_interface%array(i, j, k)*dphi2
                     tem_func%array(i, j, k) = (1.0_dp - s%array(i, j, k))*p_e_interface%array(i, j, k)*dtf%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL pw_func_u_convolution(poisson_env=poisson_env, func=tem_func, convolution=tem_convolution, u=u)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(lb,ub) &
!$OMP                SHARED(p_e_interface, t, tem_convolution, phi_tot, d_s_rhoel)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     p_e_interface%array(i, j, k) = p_e_interface%array(i, j, k)*(1.0_dp - t%array(i, j, k)) &
                                                    + tem_convolution%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         END IF

         ! Release work storage which is no longer needed
         CALL auxbas_pw_pool%give_back_pw(rho_tot_zero)
         DO i = 1, 3
            CALL auxbas_pw_pool%give_back_pw(dln_eps_elec(i))
         END DO

         ! Optional output of the total charge density in cube file format
         filename = "TOTAL_CHARGE_DENSITY"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The total SCCS charge density is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            CALL cp_pw_to_cube(rho_tot, cube_unit, TRIM(filename), particles=particles, &
                               stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         ! Calculate the total SCCS Hartree energy, potential, and its
         ! derivatives of the solute and the implicit solvent
         CALL pw_transfer(rho_tot, rho_tot_gspace)
         IF (calculate_stress_tensor) THEN
            ! Request also the calculation of the stress tensor contribution
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_tot_gspace, &
                                  ehartree=e_tot, &
                                  vhartree=v_hartree_gspace, &
                                  dvhartree=dphi_tot, &
                                  h_stress=h_stress)
         ELSE
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_tot_gspace, &
                                  ehartree=e_tot, &
                                  vhartree=v_hartree_gspace, &
                                  dvhartree=dphi_tot)
         END IF
         CALL pw_transfer(v_hartree_gspace, phi_tot)
         energy%sccs_hartree = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)

         ! Calculate the Hartree energy and potential of the solute only
         BLOCK
            TYPE(pw_r3d_rs_type) :: phi_solute
            CALL auxbas_pw_pool%create_pw(phi_solute)
            CALL pw_zero(phi_solute)
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_solute, &
                                  ehartree=energy%hartree, &
                                  vhartree=phi_solute)

            ! Calculate the polarisation potential (store it in phi_tot)
            ! phi_pol = phi_tot - phi_solute
            CALL pw_axpy(phi_solute, phi_tot, alpha=-1.0_dp)
            CALL auxbas_pw_pool%give_back_pw(phi_solute)
         END BLOCK

         ! Optional output of the SCCS polarisation potential in cube file format
         filename = "POLARISATION_POTENTIAL"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                   cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The SCCS polarisation potential is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            CALL cp_pw_to_cube(phi_tot, cube_unit, TRIM(filename), particles=particles, &
                               stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         ! Calculate the polarisation charge (store it in rho_tot)
         ! rho_pol = rho_tot - rho_solute
         CALL pw_axpy(rho_solute, rho_tot, alpha=-1.0_dp)
         polarisation_charge = pw_integrate_function(rho_tot)

         ! Optional output of the SCCS polarisation charge in cube file format
         filename = "POLARISATION_CHARGE_DENSITY"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                   cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The SCCS polarisation charge density is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            CALL cp_pw_to_cube(rho_tot, cube_unit, TRIM(filename), particles=particles, &
                               stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         ! Calculate SCCS polarisation energy
         energy%sccs_pol = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)
         CALL auxbas_pw_pool%give_back_pw(rho_solute)
         CALL auxbas_pw_pool%give_back_pw(phi_tot)
         CALL auxbas_pw_pool%give_back_pw(rho_tot)
         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%give_back_pw(s)
            CALL auxbas_pw_pool%give_back_pw(u)
            CALL auxbas_pw_pool%give_back_pw(ff)
            CALL auxbas_pw_pool%give_back_pw(t)
            CALL auxbas_pw_pool%give_back_pw(dtf)
            CALL auxbas_pw_pool%give_back_pw(tem_convolution)
            CALL auxbas_pw_pool%give_back_pw(tem_func)
            CALL auxbas_pw_pool%give_back_pw(eps_elec)
         END IF
      END BLOCK

      ! Calculate additional solvation terms
      energy%sccs_cav = sccs_control%gamma_solvent*cavity_surface
      energy%sccs_dis = sccs_control%beta_solvent*cavity_volume
      energy%sccs_rep = sccs_control%alpha_solvent*cavity_surface
      ! Calculate solvation free energy: \delta G^el + (alpha + gamma)*S + beta*V
      energy%sccs_sol = energy%sccs_pol + energy%sccs_rep + energy%sccs_cav + energy%sccs_dis

      IF (should_output .AND. (output_unit > 0)) THEN
         WRITE (UNIT=output_unit, FMT="(T3,A)") &
            "SCCS|"
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
            "SCCS| Polarisation charge", polarisation_charge
         !MK   "SCCS| Total interaction energy [a.u.]", e_tot
         WRITE (UNIT=output_unit, FMT="(T3,A)") &
            "SCCS|"
         CALL print_sccs_results(energy, sccs_control, output_unit)
      END IF

      ! Calculate SCCS contribution to the Kohn-Sham potential
      IF (sccs_control%method_id /= sccs_saa_andreussi) THEN

         f = -0.5_dp*dvol/fourpi
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dphi2,i,j,k) &
!$OMP             SHARED(f,deps_elec,dphi_tot) &
!$OMP             SHARED(lb,ub,v_sccs)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  dphi2 = dphi_tot(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                          dphi_tot(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                          dphi_tot(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k)
                  v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + f*deps_elec%array(i, j, k)*dphi2
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

      ELSE

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(lb,ub, dvol) &
!$OMP                SHARED(p_E_interface, v_sccs, d_s_rhoel)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  v_sccs%array(i, j, k) = dvol*d_s_rhoel%array(i, j, k)*p_E_interface%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

         CALL auxbas_pw_pool%give_back_pw(d_s_rhoel)
         CALL auxbas_pw_pool%give_back_pw(p_E_interface)

      END IF

      CALL auxbas_pw_pool%give_back_pw(deps_elec)
      DO i = 1, 3
         CALL auxbas_pw_pool%give_back_pw(dphi_tot(i))
      END DO

      ! Release the SCCS printout environment
      CALL cp_print_key_finished_output(output_unit, logger, input, TRIM(print_path), &
                                        ignore_should_output=should_output)

      CALL timestop(handle)

   END SUBROUTINE sccs

! **************************************************************************************************
!> \brief   Compute the AICS anisotropic dielectric electrostatic contribution.
!> \param   qs_env ...
!> \param   rho_tot_gspace ...
!> \param   v_hartree_gspace ...
!> \param   v_sccs ...
!> \param   calculate_forces ...
!> \param   h_stress Optional AICS stress contribution; its presence is rejected because it is
!>                   unsupported.
!> \author  Ziwei Chai
!> \note    The supported scope is the ANDREUSSI electrostatic/polarization contribution. The
!>          selected AB, BC, or CA cell-vector plane is mapped to an oriented Cartesian normal.
!>          Its normal axis uses the perpendicular parameters and the two in-plane Cartesian axes
!>          use the parallel parameters. Current discretizations require an axis-aligned
!>          orthorhombic cell.
!> \note    The explicitly selected surface coordinates must describe one contiguous, unwrapped
!>          slab along the selected normal, with solvent at both normal cell boundaries. A slab
!>          wrapping across the periodic normal boundary is unsupported. AICS stress is not
!>          implemented.
!> \note    FEM Cartesian potential gradients use CP2K's fixed seven-point central-difference
!>          CD7 derivative. The FD backend instead uses its discrete dielectric sensitivities.
!>          The parent SCCS DERIVATIVE_METHOD keyword does not control AICS.
! **************************************************************************************************
   SUBROUTINE sccs_aics(qs_env, rho_tot_gspace, v_hartree_gspace, v_sccs, calculate_forces, h_stress)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_tot_gspace, v_hartree_gspace
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_sccs
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: h_stress

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sccs_aics'

      CHARACTER(LEN=1023)                                :: bridge_cause, bridge_error
      CHARACTER(LEN=default_path_length)                 :: filename, mpi_filename
      CHARACTER(LEN=default_string_length)               :: component_filename, cube_file_position, &
                                                            cube_path, print_path
      INTEGER :: bridge_failed, cube_unit, failing_rank, handle, i, iatom, ispin, j, k, &
         mpi_comm_handle, natom, nspin, output_unit, parallel_axis
      INTEGER(KIND=c_int)                                :: bridge_status
      INTEGER, DIMENSION(3)                              :: aics_parameter_direction, glb, lb, &
                                                            local_end, local_npts, local_start, ub
      LOGICAL                                            :: append_cube, mpi_io, should_output, &
                                                            use_fd_dielectric_sensitivity
      REAL(KIND=dp) :: density_weighted_gradient_squared, dielectric_potential_prefactor, dvol, &
         epsilon_solvent, grid_point_coordinate, lower_interface_force_integral, &
         upper_interface_force_integral
      REAL(KIND=dp), DIMENSION(3)                        :: alpha, beta, d0, epsilon_dist_near, sigma
      TYPE(aics_interface_geometry_type)                 :: interface_geometry
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(2) :: dielectric_field_squared_sensitivity
      TYPE(pw_r3d_rs_type), DIMENSION(3) :: aics_permittivity, d_permittivity_d_density, &
         interface_permittivity_derivative, potential_gradient
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_pw_r
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(qs_env))
      IF (PRESENT(h_stress)) THEN
         CALL cp_abort(__LOCATION__, &
                       "The stress tensor is not implemented for AICS. AICS cannot currently be used "// &
                       "when FORCE_EVAL requests stress.")
      END IF

      NULLIFY (auxbas_pw_pool)
      NULLIFY (cp_subsys)
      NULLIFY (cell)
      NULLIFY (dft_control)
      NULLIFY (energy)
      NULLIFY (input)
      NULLIFY (logger)
      NULLIFY (para_env)
      NULLIFY (particles)
      NULLIFY (particle_set)
      NULLIFY (poisson_env)
      NULLIFY (pw_env)
      NULLIFY (rho)
      NULLIFY (rho_pw_r)
      NULLIFY (sccs_control)
      NULLIFY (scf_env)

      ! Load data from the Quickstep environment.
      CALL get_qs_env(qs_env=qs_env, &
                      cp_subsys=cp_subsys, &
                      cell=cell, &
                      dft_control=dft_control, &
                      energy=energy, &
                      input=input, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      pw_env=pw_env, &
                      rho=rho, &
                      scf_env=scf_env)
      CALL cp_subsys_get(cp_subsys, particles=particles)

      CPASSERT(ASSOCIATED(particle_set))

      sccs_control => dft_control%sccs_control
      logger => cp_get_default_logger()
      CALL determine_aics_orientation(cell, sccs_control%aics%interface_plane, interface_geometry)
      CALL determine_aics_interface_geometry(particle_set, &
                                             sccs_control%aics%upper_surface_atoms, &
                                             sccs_control%aics%lower_surface_atoms, &
                                             interface_geometry)
      aics_parameter_direction(:) = 1
      aics_parameter_direction(interface_geometry%normal_axis) = 2
      parallel_axis = MERGE(2, 1, interface_geometry%normal_axis == 1)
      use_fd_dielectric_sensitivity = &
         sccs_control%aics%poisson_solver == aics_poisson_solver_finite_difference

      CALL pw_env_get(pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      IF (.NOT. ALL(cell%perd == 1)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Both current AICS backends support only three-dimensional XYZ periodic boundary "// &
                       "conditions; set SUBSYS%CELL%PERIODIC XYZ.")
      END IF
      IF (.NOT. ALL(poisson_env%parameters%periodic == 1)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Both current AICS backends support only three-dimensional XYZ periodic boundary "// &
                       "conditions; set DFT%POISSON%PERIODIC XYZ.")
      END IF

      CALL pw_zero(v_sccs)

      ! Delay activation of the implicit-solvent contribution until the requested
      ! electronic SCF threshold has been reached.
      IF (.NOT. sccs_control%sccs_activated) THEN
         IF (sccs_control%eps_scf > 0.0_dp) THEN
            IF ((scf_env%iter_delta > sccs_control%eps_scf) .OR. &
                ((qs_env%scf_env%outer_scf%iter_count == 0) .AND. &
                 (qs_env%scf_env%iter_count <= 1))) THEN
               CALL pw_poisson_solve(poisson_env=poisson_env, &
                                     density=rho_tot_gspace, &
                                     ehartree=energy%hartree, &
                                     vhartree=v_hartree_gspace)
               energy%sccs_pol = 0.0_dp
               energy%sccs_cav = 0.0_dp
               energy%sccs_dis = 0.0_dp
               energy%sccs_rep = 0.0_dp
               energy%sccs_sol = 0.0_dp
               energy%sccs_hartree = energy%hartree
               CALL timestop(handle)
               RETURN
            END IF
         END IF

         ! Initialise the external anisotropic Poisson solver once SCCS/AICS is activated.
         mpi_comm_handle = para_env%get_handle()
         glb(1:3) = v_sccs%pw_grid%bounds(1, 1:3)
         lb(1:3) = v_sccs%pw_grid%bounds_local(1, 1:3)
         ub(1:3) = v_sccs%pw_grid%bounds_local(2, 1:3)
         local_npts(1:3) = MAX(ub(1:3) - lb(1:3) + 1, 0)
         local_start(1:3) = lb(1:3) - glb(1:3)
         local_end(1:3) = local_start(1:3) + local_npts(1:3)

         CALL clear_aics_field_registrations()
         CALL cp2k_aics_set_grid_metadata(INT(v_sccs%pw_grid%npts(1), KIND=c_size_t), &
                                          INT(v_sccs%pw_grid%npts(2), KIND=c_size_t), &
                                          INT(v_sccs%pw_grid%npts(3), KIND=c_size_t), &
                                          INT(local_npts(1), KIND=c_size_t), &
                                          INT(local_npts(2), KIND=c_size_t), &
                                          INT(local_npts(3), KIND=c_size_t), &
                                          INT(local_start(1), KIND=c_size_t), &
                                          INT(local_start(2), KIND=c_size_t), &
                                          INT(local_start(3), KIND=c_size_t), &
                                          INT(local_end(1), KIND=c_size_t), &
                                          INT(local_end(2), KIND=c_size_t), &
                                          INT(local_end(3), KIND=c_size_t), &
                                          REAL(v_sccs%pw_grid%dr(1), KIND=c_double), &
                                          REAL(v_sccs%pw_grid%dr(2), KIND=c_double), &
                                          REAL(v_sccs%pw_grid%dr(3), KIND=c_double))
         CALL cp2k_aics_initialize_from_fortran(mpi_comm_handle, &
                                                MERGE(1, 0, para_env%is_source()), &
                                                sccs_control%aics%poisson_solver, &
                                                interface_geometry%normal_axis - 1, &
                                                sccs_control%aics%relative_tolerance, &
                                                sccs_control%aics%absolute_tolerance, &
                                                sccs_control%aics%max_iterations, bridge_status)
         bridge_failed = MERGE(1, 0, bridge_status /= 0_c_int)
         CALL para_env%max(bridge_failed)
         IF (bridge_failed /= 0) THEN
            CALL cp2k_aics_cancel_initialize_from_fortran()
            CALL clear_aics_field_registrations()
            failing_rank = MERGE(para_env%mepos, para_env%num_pe, bridge_status /= 0_c_int)
            CALL para_env%min(failing_rank)
            bridge_error = ""
            IF (para_env%mepos == failing_rank) THEN
               CALL cp2k_aics_get_last_error_from_fortran(bridge_error)
            END IF
            CALL para_env%bcast(bridge_error, failing_rank)
            IF (LEN_TRIM(bridge_error) == 0) THEN
               bridge_error = "CP2K-AICS Python initialization failed."
            END IF
            bridge_cause = bridge_error
            WRITE (bridge_error, '(A,I0,A,A)') "CP2K-AICS pre-Dup initialization failed on rank ", &
               failing_rank, ": ", TRIM(bridge_cause)
            CPABORT(TRIM(bridge_error))
         END IF
         CALL cp2k_aics_initialize_prepared_from_fortran(bridge_status)
         CALL clear_aics_field_registrations()
         bridge_failed = MERGE(1, 0, bridge_status /= 0_c_int)
         CALL para_env%max(bridge_failed)
         IF (bridge_failed /= 0) THEN
            failing_rank = MERGE(para_env%mepos, para_env%num_pe, bridge_status /= 0_c_int)
            CALL para_env%min(failing_rank)
            bridge_error = ""
            IF (para_env%mepos == failing_rank) THEN
               CALL cp2k_aics_get_last_error_from_fortran(bridge_error)
            END IF
            CALL para_env%bcast(bridge_error, failing_rank)
            IF (LEN_TRIM(bridge_error) == 0) THEN
               bridge_error = "CP2K-AICS Python context initialization failed."
            END IF
            CPABORT(TRIM(bridge_error))
         END IF
         sccs_control%sccs_activated = .TRUE.
      END IF

      nspin = dft_control%nspins

      ! Manage SCCS print output using the standard CP2K print-key machinery.
      print_path = "DFT%PRINT%SCCS"
      should_output = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                       TRIM(print_path)), cp_p_file)
      output_unit = cp_print_key_unit_nr(logger, input, TRIM(print_path), &
                                         extension=".sccs", &
                                         ignore_should_output=should_output, &
                                         log_filename=.FALSE.)

      CALL qs_rho_get(rho_struct=rho, &
                      rho_r=rho_pw_r)

      ! Construct the anisotropic dielectric tensor and its density- and
      ! interface-position derivatives.
      BLOCK
         TYPE(pw_r3d_rs_type) :: rho_elec

         CALL auxbas_pw_pool%create_pw(rho_elec)
         dvol = rho_elec%pw_grid%dvol
         lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
         ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

         CALL pw_copy(rho_pw_r(1), rho_elec)
         DO ispin = 2, nspin
            CALL pw_axpy(rho_pw_r(ispin), rho_elec)
         END DO

         DO i = 1, 3
            CALL auxbas_pw_pool%create_pw(aics_permittivity(i))
            CALL auxbas_pw_pool%create_pw(d_permittivity_d_density(i))
            CALL auxbas_pw_pool%create_pw(interface_permittivity_derivative(i))
            CALL pw_zero(aics_permittivity(i))
            CALL pw_zero(d_permittivity_d_density(i))
            CALL pw_zero(interface_permittivity_derivative(i))
         END DO

         epsilon_solvent = sccs_control%epsilon_solvent
         ! Expand the physical (parallel, perpendicular) controls into the
         ! Cartesian (x, y, z) ordering consumed by the existing evaluator.
         epsilon_dist_near = sccs_control%aics%epsilon_dist_near(aics_parameter_direction)
         sigma = sccs_control%aics%sigma(aics_parameter_direction)
         d0 = sccs_control%aics%d0(aics_parameter_direction)
         alpha = sccs_control%aics%alpha(aics_parameter_direction)
         beta = sccs_control%aics%beta(aics_parameter_direction)

         SELECT CASE (sccs_control%method_id)
         CASE (sccs_andreussi)
            CALL build_aics_dielectric_fields(rho_elec, aics_permittivity, d_permittivity_d_density, &
                                              interface_permittivity_derivative, epsilon_solvent, &
                                              epsilon_dist_near, sigma, d0, alpha, beta, &
                                              sccs_control%rho_max, sccs_control%rho_min, &
                                              particle_set, cell, interface_geometry, auxbas_pw_pool)
         CASE DEFAULT
            CPABORT("AICS currently supports only the Andreussi density-switching method")
         END SELECT

         ! Use the first in-plane construction as the parallel dielectric. Keep the other
         ! in-plane component identical for dielectric cube output and Fortran bookkeeping.
         DO i = 2, 3
            IF (aics_parameter_direction(i) == aics_parameter_direction(1)) THEN
               CALL pw_copy(aics_permittivity(1), aics_permittivity(i))
            END IF
         END DO

         ! Optional output of the three diagonal dielectric-tensor components.
         filename = "DIELECTRIC_FUNCTION"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                   cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            cube_file_position = "REWIND"
            IF (append_cube) cube_file_position = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=cube_file_position, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The dielectric tensor is written in cube format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            DO i = 1, 3
               SELECT CASE (i)
               CASE (1)
                  component_filename = "DIELECTRIC_FUNCTION_X"
               CASE (2)
                  component_filename = "DIELECTRIC_FUNCTION_Y"
               CASE (3)
                  component_filename = "DIELECTRIC_FUNCTION_Z"
               END SELECT
               CALL cp_pw_to_cube(aics_permittivity(i), cube_unit, TRIM(component_filename), particles=particles, &
                                  stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                  mpi_io=mpi_io)
            END DO
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         CALL auxbas_pw_pool%give_back_pw(rho_elec)
      END BLOCK

      ! Solve the anisotropic Poisson equation and construct the electrostatic energy.
      BLOCK
         TYPE(pw_r3d_rs_type) :: phi_tot, rho_solute

         CALL auxbas_pw_pool%create_pw(rho_solute)
         CALL pw_zero(rho_solute)
         CALL pw_transfer(rho_tot_gspace, rho_solute)

         IF (sccs_control%aics%poisson_solver == aics_poisson_solver_finite_element) THEN
            DO i = 1, 3
               CALL auxbas_pw_pool%create_pw(potential_gradient(i))
               CALL pw_zero(potential_gradient(i))
            END DO
         END IF
         IF (use_fd_dielectric_sensitivity) THEN
            DO i = 1, 2
               CALL auxbas_pw_pool%create_pw(dielectric_field_squared_sensitivity(i))
               CALL pw_zero(dielectric_field_squared_sensitivity(i))
            END DO
         END IF

         CALL auxbas_pw_pool%create_pw(phi_tot)
         CALL pw_zero(phi_tot)

         IF (use_fd_dielectric_sensitivity) THEN
            CALL solve_aics_poisson_zero_copy(rho_solute, aics_permittivity, phi_tot, para_env, &
                                              interface_geometry%normal_axis, &
                                              dielectric_field_squared_sensitivity)
         ELSE
            CALL solve_aics_poisson_zero_copy(rho_solute, aics_permittivity, phi_tot, para_env, &
                                              interface_geometry%normal_axis)
         END IF

         IF (sccs_control%aics%poisson_solver == aics_poisson_solver_finite_element) THEN
            CALL derive(phi_tot, potential_gradient, sccs_derivative_cd7, pw_env, input)
         END IF
         CALL pw_transfer(phi_tot, v_hartree_gspace)

         energy%sccs_hartree = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)

         ! Vacuum Hartree energy of the same solute charge density.
         CALL pw_poisson_solve(poisson_env=poisson_env, &
                               density=rho_solute, &
                               ehartree=energy%hartree)

         energy%sccs_pol = energy%sccs_hartree - energy%hartree

         CALL auxbas_pw_pool%give_back_pw(rho_solute)
         CALL auxbas_pw_pool%give_back_pw(phi_tot)
      END BLOCK

      ! The present AICS implementation contains only the electrostatic term.
      IF ((ABS(sccs_control%gamma_solvent) > 0.0_dp) .OR. &
          (ABS(sccs_control%beta_solvent) > 0.0_dp) .OR. &
          (ABS(sccs_control%alpha_solvent) > 0.0_dp)) THEN
         CALL cp_abort(__LOCATION__, &
                       "AICS supports only the ANDREUSSI electrostatic/polarization contribution. "// &
                       "Top-level SCCS ALPHA, BETA, and GAMMA must be zero; these controls are distinct "// &
                       "from SCCS%AICS%ALPHA and SCCS%AICS%BETA.")
      END IF
      energy%sccs_cav = 0.0_dp
      energy%sccs_dis = 0.0_dp
      energy%sccs_rep = 0.0_dp
      energy%sccs_sol = energy%sccs_pol

      IF (should_output .AND. (output_unit > 0)) THEN
         WRITE (UNIT=output_unit, FMT="(T3,A)") "SCCS|"
         CALL print_sccs_results(energy, sccs_control, output_unit)
      END IF

      natom = SIZE(particle_set)

      glb(1:3) = d_permittivity_d_density(1)%pw_grid%bounds(1, 1:3)

      ! Density-dependent contribution to the Kohn-Sham potential and the
      ! force associated with motion of the two dielectric interfaces.
      ! The FD backend returns the unprefactored g_parallel and
      ! g_perpendicular quantities, so the -dV/(8*pi) prefactor is applied
      ! explicitly here. In contrast, the interface-permittivity derivative
      ! used below already includes (1/(8*pi))*d epsilon/dd, consistent with
      ! the established FEM formulation.
      dielectric_potential_prefactor = -0.5_dp*dvol/fourpi
      upper_interface_force_integral = 0.0_dp
      lower_interface_force_integral = 0.0_dp

      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP             PRIVATE(density_weighted_gradient_squared,grid_point_coordinate,i,j,k) &
      !$OMP             SHARED(calculate_forces,d_permittivity_d_density,dielectric_potential_prefactor,dvol,glb, &
      !$OMP                    dielectric_field_squared_sensitivity,interface_geometry, &
      !$OMP                    interface_permittivity_derivative,lb,parallel_axis,potential_gradient,ub, &
      !$OMP                    use_fd_dielectric_sensitivity,v_sccs) &
      !$OMP             REDUCTION(+:lower_interface_force_integral,upper_interface_force_integral)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               IF (use_fd_dielectric_sensitivity) THEN
                  density_weighted_gradient_squared = &
                     d_permittivity_d_density(parallel_axis)%array(i, j, k)* &
                     dielectric_field_squared_sensitivity(1)%array(i, j, k) + &
                     d_permittivity_d_density(interface_geometry%normal_axis)%array(i, j, k)* &
                     dielectric_field_squared_sensitivity(2)%array(i, j, k)
               ELSE
                  density_weighted_gradient_squared = &
                     d_permittivity_d_density(1)%array(i, j, k)*potential_gradient(1)%array(i, j, k)**2 + &
                     d_permittivity_d_density(2)%array(i, j, k)*potential_gradient(2)%array(i, j, k)**2 + &
                     d_permittivity_d_density(3)%array(i, j, k)*potential_gradient(3)%array(i, j, k)**2
               END IF

               v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + &
                                       dielectric_potential_prefactor*density_weighted_gradient_squared

               IF (calculate_forces) THEN
                  SELECT CASE (interface_geometry%normal_axis)
                  CASE (1)
                     grid_point_coordinate = interface_geometry%normal_sign* &
                                             (0.5_dp + REAL(i - glb(1), KIND=dp))*d_permittivity_d_density(1)%pw_grid%dr(1)
                  CASE (2)
                     grid_point_coordinate = interface_geometry%normal_sign* &
                                             (0.5_dp + REAL(j - glb(2), KIND=dp))*d_permittivity_d_density(1)%pw_grid%dr(2)
                  CASE (3)
                     grid_point_coordinate = interface_geometry%normal_sign* &
                                             (0.5_dp + REAL(k - glb(3), KIND=dp))*d_permittivity_d_density(1)%pw_grid%dr(3)
                  END SELECT
                  IF (grid_point_coordinate > interface_geometry%upper_interface_coordinate) THEN
                     IF (use_fd_dielectric_sensitivity) THEN
                        upper_interface_force_integral = upper_interface_force_integral + &
                                                         (interface_permittivity_derivative(parallel_axis)%array(i, j, k)* &
                                                          dielectric_field_squared_sensitivity(1)%array(i, j, k) + &
                                                 interface_permittivity_derivative(interface_geometry%normal_axis)%array(i, j, k)* &
                                                          dielectric_field_squared_sensitivity(2)%array(i, j, k))*dvol
                     ELSE
                        upper_interface_force_integral = upper_interface_force_integral + &
                                                         (interface_permittivity_derivative(1)%array(i, j, k)* &
                                                          potential_gradient(1)%array(i, j, k)**2 + &
                                                          interface_permittivity_derivative(2)%array(i, j, k)* &
                                                          potential_gradient(2)%array(i, j, k)**2 + &
                                                          interface_permittivity_derivative(3)%array(i, j, k)* &
                                                          potential_gradient(3)%array(i, j, k)**2)*dvol
                     END IF

                  ELSE IF (grid_point_coordinate < interface_geometry%lower_interface_coordinate) THEN
                     IF (use_fd_dielectric_sensitivity) THEN
                        lower_interface_force_integral = lower_interface_force_integral + &
                                                         (interface_permittivity_derivative(parallel_axis)%array(i, j, k)* &
                                                          dielectric_field_squared_sensitivity(1)%array(i, j, k) + &
                                                 interface_permittivity_derivative(interface_geometry%normal_axis)%array(i, j, k)* &
                                                          dielectric_field_squared_sensitivity(2)%array(i, j, k))*dvol
                     ELSE
                        lower_interface_force_integral = lower_interface_force_integral + &
                                                         (interface_permittivity_derivative(1)%array(i, j, k)* &
                                                          potential_gradient(1)%array(i, j, k)**2 + &
                                                          interface_permittivity_derivative(2)%array(i, j, k)* &
                                                          potential_gradient(2)%array(i, j, k)**2 + &
                                                          interface_permittivity_derivative(3)%array(i, j, k)* &
                                                          potential_gradient(3)%array(i, j, k)**2)*dvol
                     END IF
                  END IF
               END IF

            END DO
         END DO
      END DO
      !$OMP END PARALLEL DO

      IF (calculate_forces) THEN
         ! Reduce the direct interface-force integrals and apply each integral
         ! to atoms in its corresponding explicit surface list.
         CALL para_env%sum(upper_interface_force_integral)
         CALL para_env%sum(lower_interface_force_integral)

         DO iatom = 1, natom
            IF (interface_geometry%upper_surface_atom(iatom)) THEN
               particle_set(iatom)%f(interface_geometry%normal_axis) = &
                  particle_set(iatom)%f(interface_geometry%normal_axis) - &
                  interface_geometry%normal_sign*upper_interface_force_integral
            END IF

            IF (interface_geometry%lower_surface_atom(iatom)) THEN
               particle_set(iatom)%f(interface_geometry%normal_axis) = &
                  particle_set(iatom)%f(interface_geometry%normal_axis) - &
                  interface_geometry%normal_sign*lower_interface_force_integral
            END IF
         END DO
      END IF

      DEALLOCATE (interface_geometry%upper_surface_atom)
      DEALLOCATE (interface_geometry%lower_surface_atom)

      DO i = 1, 3
         CALL auxbas_pw_pool%give_back_pw(aics_permittivity(i))
         CALL auxbas_pw_pool%give_back_pw(d_permittivity_d_density(i))
         CALL auxbas_pw_pool%give_back_pw(interface_permittivity_derivative(i))
         IF (sccs_control%aics%poisson_solver == aics_poisson_solver_finite_element) THEN
            CALL auxbas_pw_pool%give_back_pw(potential_gradient(i))
         END IF
      END DO
      IF (use_fd_dielectric_sensitivity) THEN
         DO i = 1, 2
            CALL auxbas_pw_pool%give_back_pw(dielectric_field_squared_sensitivity(i))
         END DO
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, input, TRIM(print_path), &
                                        ignore_should_output=should_output)

      CALL timestop(handle)

   END SUBROUTINE sccs_aics

! **************************************************************************************************
!> \brief   Perform one synchronous collective zero-copy AICS Poisson solve.
!> \param   rho_solute                              Local real-space solute charge-density field
!>                                                  borrowed by the bridge.
!> \param   aics_permittivity                       Three Cartesian diagonal dielectric fields on
!>                                                  the local PW domain.
!> \param   phi_tot                                 Local real-space electrostatic-potential output
!>                                                  field borrowed by the bridge.
!> \param   para_env ...
!> \param   normal_axis                             Cartesian component index of the oriented
!>                                                  interface normal.
!> \param   dielectric_field_squared_sensitivity    Optional CP2K-owned parallel/perpendicular FD
!>                                                  outputs.
!> \author  Ziwei Chai
!> \note    Each MPI rank owns only its local portion of the global real-space grid. The bridge
!>          passes the local array pointers together with their global-grid location and size to
!>          the Python backend for this solve only. A rank with no local grid points passes
!>          C_NULL_PTR and zero extents. The pointers are borrowed from CP2K and must not be
!>          retained after the synchronous solve returns.
! **************************************************************************************************
   SUBROUTINE solve_aics_poisson_zero_copy(rho_solute, aics_permittivity, phi_tot, para_env, normal_axis, &
                                           dielectric_field_squared_sensitivity)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_solute
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: aics_permittivity
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: phi_tot
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: normal_axis
      TYPE(pw_r3d_rs_type), DIMENSION(2), INTENT(INOUT), OPTIONAL :: &
         dielectric_field_squared_sensitivity

      CHARACTER(LEN=1023)                                :: bridge_error
      INTEGER                                            :: bridge_failed, failing_rank, &
                                                            parallel_axis
      INTEGER(KIND=c_int)                                :: bridge_status
      INTEGER, DIMENSION(3)                              :: glb, lb, local_end, local_npts, &
                                                            local_start, ub

      glb(1:3) = aics_permittivity(1)%pw_grid%bounds(1, 1:3)
      lb(1:3) = aics_permittivity(1)%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = aics_permittivity(1)%pw_grid%bounds_local(2, 1:3)
      local_npts(1:3) = MAX(ub(1:3) - lb(1:3) + 1, 0)
      local_start(1:3) = lb(1:3) - glb(1:3)
      local_end(1:3) = local_start(1:3) + local_npts(1:3)

      parallel_axis = MERGE(2, 1, normal_axis == 1)
      CALL clear_aics_field_registrations()
      IF (ALL(local_npts > 0)) THEN
         CALL cp2k_aics_set_solute_charge_density(c_loc(rho_solute%array(lb(1), lb(2), lb(3))))
         CALL cp2k_aics_set_relative_permittivity_parallel( &
            c_loc(aics_permittivity(parallel_axis)%array(lb(1), lb(2), lb(3))))
         CALL cp2k_aics_set_relative_permittivity_perpendicular( &
            c_loc(aics_permittivity(normal_axis)%array(lb(1), lb(2), lb(3))))
         CALL cp2k_aics_set_electrostatic_potential(c_loc(phi_tot%array(lb(1), lb(2), lb(3))))
         IF (PRESENT(dielectric_field_squared_sensitivity)) THEN
            CALL cp2k_aics_set_fd_field_squared_sensitivity_parallel( &
               c_loc(dielectric_field_squared_sensitivity(1)%array(lb(1), lb(2), lb(3))))
            CALL cp2k_aics_set_fd_field_squared_sensitivity_perpendicular( &
               c_loc(dielectric_field_squared_sensitivity(2)%array(lb(1), lb(2), lb(3))))
         END IF
      ELSE
         CALL cp2k_aics_set_solute_charge_density(c_null_ptr)
         CALL cp2k_aics_set_relative_permittivity_parallel(c_null_ptr)
         CALL cp2k_aics_set_relative_permittivity_perpendicular(c_null_ptr)
         CALL cp2k_aics_set_electrostatic_potential(c_null_ptr)
         IF (PRESENT(dielectric_field_squared_sensitivity)) THEN
            CALL cp2k_aics_set_fd_field_squared_sensitivity_parallel(c_null_ptr)
            CALL cp2k_aics_set_fd_field_squared_sensitivity_perpendicular(c_null_ptr)
         END IF
      END IF

      CALL cp2k_aics_set_grid_metadata(INT(aics_permittivity(1)%pw_grid%npts(1), KIND=c_size_t), &
                                       INT(aics_permittivity(1)%pw_grid%npts(2), KIND=c_size_t), &
                                       INT(aics_permittivity(1)%pw_grid%npts(3), KIND=c_size_t), &
                                       INT(local_npts(1), KIND=c_size_t), &
                                       INT(local_npts(2), KIND=c_size_t), &
                                       INT(local_npts(3), KIND=c_size_t), &
                                       INT(local_start(1), KIND=c_size_t), &
                                       INT(local_start(2), KIND=c_size_t), &
                                       INT(local_start(3), KIND=c_size_t), &
                                       INT(local_end(1), KIND=c_size_t), &
                                       INT(local_end(2), KIND=c_size_t), &
                                       INT(local_end(3), KIND=c_size_t), &
                                       REAL(aics_permittivity(1)%pw_grid%dr(1), KIND=c_double), &
                                       REAL(aics_permittivity(1)%pw_grid%dr(2), KIND=c_double), &
                                       REAL(aics_permittivity(1)%pw_grid%dr(3), KIND=c_double))

      CALL cp2k_aics_solve_from_fortran(bridge_status)
      CALL clear_aics_field_registrations()
      bridge_failed = MERGE(1, 0, bridge_status /= 0_c_int)
      CALL para_env%max(bridge_failed)
      IF (bridge_failed /= 0) THEN
         failing_rank = MERGE(para_env%mepos, para_env%num_pe, bridge_status /= 0_c_int)
         CALL para_env%min(failing_rank)
         bridge_error = ""
         IF (para_env%mepos == failing_rank) THEN
            CALL cp2k_aics_get_last_error_from_fortran(bridge_error)
         END IF
         CALL para_env%bcast(bridge_error, failing_rank)
         IF (LEN_TRIM(bridge_error) == 0) THEN
            bridge_error = "CP2K-AICS Python solve failed."
         END IF
         CPABORT(TRIM(bridge_error))
      END IF

   END SUBROUTINE solve_aics_poisson_zero_copy

! **************************************************************************************************
!> \brief   Construct explicit upper and lower surface masks and their mean coordinates along the
!>          oriented interface normal.
!> \param   particle_set ...
!> \param   upper_surface_atoms     Global 1-based atom indices defining the upper reference
!>                                  plane.
!> \param   lower_surface_atoms     Global 1-based atom indices defining the lower reference
!>                                  plane.
!> \param   interface_geometry      Interface orientation together with the constructed surface
!>                                  masks, atom counts, and mean coordinates.
!> \author  Ziwei Chai
!> \note    Each reference coordinate is the arithmetic mean of the selected atomic coordinates
!>          projected onto unit_normal. The selected atoms must represent one contiguous,
!>          unwrapped slab, with solvent at both normal cell boundaries. A slab crossing the
!>          periodic normal boundary is not supported.
! **************************************************************************************************
   SUBROUTINE determine_aics_interface_geometry(particle_set, upper_surface_atoms, &
                                                lower_surface_atoms, interface_geometry)

      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, DIMENSION(:), INTENT(IN)                  :: upper_surface_atoms, lower_surface_atoms
      TYPE(aics_interface_geometry_type), INTENT(INOUT)  :: interface_geometry

      INTEGER                                            :: iatom, natom
      REAL(KIND=dp)                                      :: atom_coordinate

      CPASSERT(ASSOCIATED(particle_set))

      natom = SIZE(particle_set)
      ALLOCATE (interface_geometry%upper_surface_atom(natom))
      ALLOCATE (interface_geometry%lower_surface_atom(natom))
      interface_geometry%upper_surface_atom(:) = .FALSE.
      interface_geometry%lower_surface_atom(:) = .FALSE.

      IF (ANY(upper_surface_atoms > natom)) THEN
         CALL cp_abort(__LOCATION__, &
                       "DFT%SCCS%AICS%UPPER_SURFACE_ATOMS contains an atom index greater than NATOM")
      END IF
      IF (ANY(lower_surface_atoms > natom)) THEN
         CALL cp_abort(__LOCATION__, &
                       "DFT%SCCS%AICS%LOWER_SURFACE_ATOMS contains an atom index greater than NATOM")
      END IF

      interface_geometry%upper_interface_coordinate = 0.0_dp
      interface_geometry%lower_interface_coordinate = 0.0_dp
      interface_geometry%n_upper_surface_atoms = SIZE(upper_surface_atoms)
      interface_geometry%n_lower_surface_atoms = SIZE(lower_surface_atoms)
      DO iatom = 1, SIZE(upper_surface_atoms)
         interface_geometry%upper_surface_atom(upper_surface_atoms(iatom)) = .TRUE.
         atom_coordinate = DOT_PRODUCT(particle_set(upper_surface_atoms(iatom))%r, &
                                       interface_geometry%unit_normal)
         interface_geometry%upper_interface_coordinate = &
            interface_geometry%upper_interface_coordinate + atom_coordinate
      END DO
      DO iatom = 1, SIZE(lower_surface_atoms)
         interface_geometry%lower_surface_atom(lower_surface_atoms(iatom)) = .TRUE.
         atom_coordinate = DOT_PRODUCT(particle_set(lower_surface_atoms(iatom))%r, &
                                       interface_geometry%unit_normal)
         interface_geometry%lower_interface_coordinate = &
            interface_geometry%lower_interface_coordinate + atom_coordinate
      END DO

      interface_geometry%upper_interface_coordinate = interface_geometry%upper_interface_coordinate/ &
                                                      REAL(interface_geometry%n_upper_surface_atoms, KIND=dp)
      interface_geometry%lower_interface_coordinate = interface_geometry%lower_interface_coordinate/ &
                                                      REAL(interface_geometry%n_lower_surface_atoms, KIND=dp)

      IF (interface_geometry%upper_interface_coordinate == &
          interface_geometry%lower_interface_coordinate) THEN
         CALL cp_abort(__LOCATION__, &
                       "DFT%SCCS%AICS%UPPER_SURFACE_ATOMS and "// &
                       "DFT%SCCS%AICS%LOWER_SURFACE_ATOMS define coincident planes")
      END IF
      IF (interface_geometry%upper_interface_coordinate < &
          interface_geometry%lower_interface_coordinate) THEN
         CALL cp_abort(__LOCATION__, &
                       "DFT%SCCS%AICS%UPPER_SURFACE_ATOMS must define a plane above "// &
                       "DFT%SCCS%AICS%LOWER_SURFACE_ATOMS along the oriented interface normal")
      END IF

   END SUBROUTINE determine_aics_interface_geometry

! **************************************************************************************************
!> \brief   Map an AB, BC, or CA cell-vector plane to its oriented Cartesian normal.
!> \param   cell                 CP2K cell whose vector columns define the A, B, and C lattice
!>                               vectors.
!> \param   interface_plane      Selected AB, BC, or CA interface plane.
!> \param   interface_geometry   Stores the unit normal, aligned Cartesian axis, and orientation
!>                               sign of the selected interface plane.
!> \author  Ziwei Chai
!> \note    The oriented normals are defined by the cross products A x B, B x C, and C x A.
!>          unit_normal is the corresponding normalized Cartesian vector, normal_axis identifies
!>          its aligned Cartesian component, and normal_sign gives the sign of that component.
!>          The current implementation requires an axis-aligned orthorhombic cell, for which AB
!>          maps to z, BC to x, and CA to y.
! **************************************************************************************************
   SUBROUTINE determine_aics_orientation(cell, interface_plane, interface_geometry)
      TYPE(cell_type), INTENT(IN), POINTER               :: cell
      INTEGER, INTENT(IN)                                :: interface_plane
      TYPE(aics_interface_geometry_type), INTENT(INOUT)  :: interface_geometry

      REAL(KIND=dp), PARAMETER :: alignment_tolerance = 1.0E-12_dp

      INTEGER                                            :: i, j, vector_1, vector_2
      REAL(KIND=dp)                                      :: cell_scale, normal_norm
      REAL(KIND=dp), DIMENSION(3)                        :: normal

      SELECT CASE (interface_plane)
      CASE (aics_interface_plane_ab)
         vector_1 = 1; vector_2 = 2
      CASE (aics_interface_plane_bc)
         vector_1 = 2; vector_2 = 3
      CASE (aics_interface_plane_ca)
         vector_1 = 3; vector_2 = 1
      CASE DEFAULT
         CPABORT("Invalid AICS INTERFACE_PLANE")
      END SELECT
      normal = [cell%hmat(2, vector_1)*cell%hmat(3, vector_2) - &
                cell%hmat(3, vector_1)*cell%hmat(2, vector_2), &
                cell%hmat(3, vector_1)*cell%hmat(1, vector_2) - &
                cell%hmat(1, vector_1)*cell%hmat(3, vector_2), &
                cell%hmat(1, vector_1)*cell%hmat(2, vector_2) - &
                cell%hmat(2, vector_1)*cell%hmat(1, vector_2)]
      normal_norm = SQRT(DOT_PRODUCT(normal, normal))
      cell_scale = MAXVAL(ABS(cell%hmat))
      IF (normal_norm <= alignment_tolerance*cell_scale*cell_scale) THEN
         CPABORT("AICS INTERFACE_PLANE is degenerate in the current CP2K cell")
      END IF
      interface_geometry%unit_normal = normal/normal_norm
      interface_geometry%normal_axis = MAXLOC(ABS(interface_geometry%unit_normal), DIM=1)
      IF (MAXVAL(ABS(interface_geometry%unit_normal), &
                 MASK=[(i /= interface_geometry%normal_axis, i=1, 3)]) > alignment_tolerance) THEN
         CALL cp_abort(__LOCATION__, &
                       "The current AICS implementation requires the normal of INTERFACE_PLANE "// &
                       "to be aligned with a Cartesian real-space grid axis.")
      END IF
      interface_geometry%normal_sign = MERGE(1, -1, &
                                             interface_geometry%unit_normal(interface_geometry%normal_axis) > 0.0_dp)

      ! The AICS Cartesian dr-based coordinates, seven-point stencil, and box FE mesh
      ! do not represent metric cross terms of a non-orthorhombic working cell.
      DO j = 1, 3
         DO i = 1, 3
            IF (i /= j .AND. ABS(cell%hmat(i, j)) > alignment_tolerance*cell_scale) THEN
               CALL cp_abort(__LOCATION__, &
                             "The current AICS Cartesian real-space discretizations require an "// &
                             "axis-aligned orthorhombic CP2K cell.")
            END IF
         END DO
      END DO
   END SUBROUTINE determine_aics_orientation

! **************************************************************************************************
!> \brief   Build the Cartesian diagonal AICS dielectric tensor and its derivative fields.
!> \param   rho_elec ...
!> \param   aics_permittivity                   Output Cartesian diagonal dielectric components
!>                                              epsilon_xx, epsilon_yy, and epsilon_zz.
!> \param   d_permittivity_d_density            Output derivatives of each dielectric component
!>                                              with respect to the electronic density.
!> \param   interface_permittivity_derivative   Output derivatives associated with displacement
!>                                              of the dielectric interfaces.
!> \param   bulk_relative_permittivity          Bulk-solvent relative permittivity.
!> \param   epsilon_dist_near                   Three Cartesian near-interface relative
!>                                              permittivities.
!> \param   sigma                               Three Cartesian interface-transition widths in
!>                                              CP2K internal units.
!> \param   d0                                  Three Cartesian interface midpoint distances in
!>                                              CP2K internal units.
!> \param   alpha                               Three Cartesian dimensionless mapping-asymmetry
!>                                              parameters.
!> \param   beta                                Three Cartesian mapping parameters in inverse
!>                                              CP2K length units.
!> \param   rho_max                             Upper electronic-density switching threshold.
!> \param   rho_min                             Lower electronic-density switching threshold.
!> \param   particle_set                        Atomic positions used to construct the exclusion
!>                                              mask.
!> \param   cell                                Periodic simulation cell used for minimum-image
!>                                              atom-grid displacements.
!> \param   interface_geometry                  Oriented interface normal and explicit upper and
!>                                              lower reference surfaces.
!> \param   auxbas_pw_pool                      PW pool supplying the temporary atom-exclusion
!>                                              mask.
!> \author  Ziwei Chai
!> \note    All Cartesian dielectric components use the same scalar distance coordinate relative
!>          to the explicit interfaces. The Cartesian axis aligned with the oriented normal uses
!>          the perpendicular parameter set, while the two in-plane axes use the parallel set.
!>          Negative exclusion-mask values identify atom-excluded grid points, where the
!>          dielectric components are set to unity and their derivatives are set to zero.
! **************************************************************************************************
   SUBROUTINE build_aics_dielectric_fields(rho_elec, aics_permittivity, d_permittivity_d_density, &
                                           interface_permittivity_derivative, bulk_relative_permittivity, &
                                           epsilon_dist_near, sigma, d0, alpha, beta, rho_max, &
                                           rho_min, particle_set, cell, interface_geometry, auxbas_pw_pool)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT) :: aics_permittivity, &
         d_permittivity_d_density, interface_permittivity_derivative
      REAL(KIND=dp), INTENT(IN)                          :: bulk_relative_permittivity, &
                                                            epsilon_dist_near(3), sigma(3), d0(3), &
                                                            alpha(3), beta(3), rho_max, rho_min
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(aics_interface_geometry_type), INTENT(IN)     :: interface_geometry
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_aics_dielectric_fields'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, icomponent, interface_side, &
                                                            j, k
      INTEGER, DIMENSION(3)                              :: glb, lb, ub
      REAL(KIND=dp) :: density_log, density_switch_derivative_factor, density_switch_exponent, &
         density_switch_phase, density_switch_phase_derivative, diff, exclusion_mask_value, &
         grid_point_coordinate, interface_distance, ln_rho_max, ln_rho_min, &
         mapped_interface_distance, near_interface_log_factor, q, rho
      TYPE(pw_r3d_rs_type)                               :: atom_exclusion_mask

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(ASSOCIATED(auxbas_pw_pool))

      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) THEN
         CPABORT("DFT%SCCS%ANDREUSSI%RHO_MAX - DFT%SCCS%ANDREUSSI%RHO_MIN must be at least 1.0E-6 for AICS")
      END IF
      IF (rho_min < rhotol) THEN
         CPABORT("DFT%SCCS%ANDREUSSI%RHO_MIN must be at least 1.0E-12 for AICS")
      END IF

      ln_rho_max = LOG(rho_max)
      ln_rho_min = LOG(rho_min)
      q = twopi/(ln_rho_max - ln_rho_min)

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)
      glb(1:3) = rho_elec%pw_grid%bounds(1, 1:3)

      CALL auxbas_pw_pool%create_pw(atom_exclusion_mask)
      CALL build_aics_atom_exclusion_mask(atom_exclusion_mask, particle_set, cell, glb, lb, ub)

      ! Construct the diagonal dielectric and its density and interface-position derivatives.
      DO icomponent = 1, 3
         IF (ABS(alpha(icomponent)) >= 1.0_dp) THEN
            CPABORT("DFT%SCCS%AICS%ALPHA values must satisfy ABS(ALPHA) < 1")
         END IF
         IF (sigma(icomponent) <= 0.0_dp) THEN
            CPABORT("DFT%SCCS%AICS%SIGMA values must be strictly positive")
         END IF

         near_interface_log_factor = LOG(epsilon_dist_near(icomponent))/twopi
         density_switch_derivative_factor = -near_interface_log_factor*q

         ! Calculate the dielectric function and its derivatives.
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(density_log,density_switch_exponent,density_switch_phase, &
!$OMP                     density_switch_phase_derivative,exclusion_mask_value,grid_point_coordinate,i, &
!$OMP                     interface_distance,interface_side,j,k,mapped_interface_distance,rho) &
!$OMP             SHARED(aics_permittivity,alpha,atom_exclusion_mask,beta,bulk_relative_permittivity, &
!$OMP                    d_permittivity_d_density,density_switch_derivative_factor,glb, &
!$OMP                    icomponent,interface_permittivity_derivative,lb,ln_rho_max, &
!$OMP                    epsilon_dist_near,interface_geometry,near_interface_log_factor, &
!$OMP                    q,rho_elec,rho_max,rho_min,d0,sigma,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  rho = rho_elec%array(i, j, k)
                  exclusion_mask_value = atom_exclusion_mask%array(i, j, k)

                  IF ((rho < rho_min) .AND. (exclusion_mask_value > 0.0_dp)) THEN
                     ! Pure solvent. Solvent is permitted only outside the two
                     ! reference surfaces; low-density points inside the slab are
                     ! treated as solute to avoid an undefined interface direction.
                     SELECT CASE (interface_geometry%normal_axis)
                     CASE (1)
                        grid_point_coordinate = interface_geometry%normal_sign* &
                                                (0.5_dp + REAL(i - glb(1), KIND=dp))*aics_permittivity(icomponent)%pw_grid%dr(1)
                     CASE (2)
                        grid_point_coordinate = interface_geometry%normal_sign* &
                                                (0.5_dp + REAL(j - glb(2), KIND=dp))*aics_permittivity(icomponent)%pw_grid%dr(2)
                     CASE (3)
                        grid_point_coordinate = interface_geometry%normal_sign* &
                                                (0.5_dp + REAL(k - glb(3), KIND=dp))*aics_permittivity(icomponent)%pw_grid%dr(3)
                     END SELECT
                     IF (grid_point_coordinate > interface_geometry%upper_interface_coordinate) THEN
                        interface_distance = &
                           grid_point_coordinate - interface_geometry%upper_interface_coordinate
                        interface_side = 1
                     ELSE IF (grid_point_coordinate < interface_geometry%lower_interface_coordinate) THEN
                        interface_distance = &
                           interface_geometry%lower_interface_coordinate - grid_point_coordinate
                        interface_side = -1
                     ELSE
                        aics_permittivity(icomponent)%array(i, j, k) = 1.0_dp
                        d_permittivity_d_density(icomponent)%array(i, j, k) = 0.0_dp
                        interface_permittivity_derivative(icomponent)%array(i, j, k) = 0.0_dp
                        CYCLE
                     END IF

                     mapped_interface_distance = interface_distance + &
                                                 (interface_distance - d0(icomponent))/ &
                                                 (1.0_dp + alpha(icomponent)*TANH(beta(icomponent)* &
                                                                                  (interface_distance - d0(icomponent))))

                     aics_permittivity(icomponent)%array(i, j, k) = epsilon_dist_near(icomponent) + &
                                                                    (bulk_relative_permittivity - &
                                                                     epsilon_dist_near(icomponent))* &
                                                                    (1.0_dp + ERF((mapped_interface_distance - &
                                                                                   d0(icomponent))/sigma(icomponent)))*0.5_dp
                     d_permittivity_d_density(icomponent)%array(i, j, k) = 0.0_dp

                     ! Analytical permittivity derivative for displacement of one
                     ! contributing interface atom.  The lower-surface convention
                     ! is negative and the upper-surface convention is positive;
                     ! each is divided by its corresponding surface-atom count.
                     interface_permittivity_derivative(icomponent)%array(i, j, k) = &
                        (0.5_dp/fourpi)*(bulk_relative_permittivity - &
                                         epsilon_dist_near(icomponent))/(SQRT(twopi/2.0_dp)*sigma(icomponent))* &
                        EXP(-((mapped_interface_distance - d0(icomponent))/sigma(icomponent))**2)* &
                        (1.0_dp + &
                         ((1.0_dp + alpha(icomponent)*TANH(beta(icomponent)*(interface_distance - &
                                                                           d0(icomponent)))) - alpha(icomponent)*beta(icomponent)* &
                          (interface_distance - d0(icomponent))* &
                          (1.0_dp - TANH(beta(icomponent)*(interface_distance - &
                                                           d0(icomponent)))**2))/ &
                         (1.0_dp + alpha(icomponent)*TANH(beta(icomponent)*(interface_distance - &
                                                                            d0(icomponent))))**2)

                     IF (interface_side < 0) THEN
                        interface_permittivity_derivative(icomponent)%array(i, j, k) = &
                           -interface_permittivity_derivative(icomponent)%array(i, j, k)/ &
                           REAL(interface_geometry%n_lower_surface_atoms, KIND=dp)
                     ELSE
                        interface_permittivity_derivative(icomponent)%array(i, j, k) = &
                           interface_permittivity_derivative(icomponent)%array(i, j, k)/ &
                           REAL(interface_geometry%n_upper_surface_atoms, KIND=dp)
                     END IF

                  ELSE IF ((rho <= rho_max) .AND. (exclusion_mask_value > 0.0_dp)) THEN
                     ! Density-defined transition region.
                     density_log = LOG(rho)
                     density_switch_phase = q*(ln_rho_max - density_log)
                     density_switch_exponent = near_interface_log_factor* &
                                               (density_switch_phase - SIN(density_switch_phase))
                     aics_permittivity(icomponent)%array(i, j, k) = EXP(density_switch_exponent)
                     density_switch_phase_derivative = density_switch_derivative_factor* &
                                                       (1.0_dp - COS(density_switch_phase))
                     d_permittivity_d_density(icomponent)%array(i, j, k) = &
                        aics_permittivity(icomponent)%array(i, j, k)*density_switch_phase_derivative/rho
                     interface_permittivity_derivative(icomponent)%array(i, j, k) = 0.0_dp

                  ELSE
                     ! Pure solute or an atom-centred exclusion point.
                     aics_permittivity(icomponent)%array(i, j, k) = 1.0_dp
                     d_permittivity_d_density(icomponent)%array(i, j, k) = 0.0_dp
                     interface_permittivity_derivative(icomponent)%array(i, j, k) = 0.0_dp
                  END IF
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      END DO

      CALL auxbas_pw_pool%give_back_pw(atom_exclusion_mask)
      CALL timestop(handle)

   END SUBROUTINE build_aics_dielectric_fields

! **************************************************************************************************
!> \brief   Build the hard all-atom exclusion mask for the local AICS grid domain.
!> \param   atom_exclusion_mask   Output mask: +1 permits dielectric mapping and -1 excludes it.
!> \param   particle_set ...
!> \param   cell ...
!> \param   glb                   Global PW lower bounds used to construct Cartesian grid-point
!>                                coordinates.
!> \param   lb                    Inclusive local lower PW bounds.
!> \param   ub                    Inclusive local upper PW bounds.
!> \author  Ziwei Chai
!> \note    Every atom defines a hard exclusion sphere of radius 0.25 angstrom on the cell-centred
!>          real-space grid, using a strict distance test and periodic minimum images. No analytical
!>          force is included for motion of the discontinuous mask boundary.
! **************************************************************************************************
   SUBROUTINE build_aics_atom_exclusion_mask(atom_exclusion_mask, particle_set, cell, glb, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: atom_exclusion_mask
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, DIMENSION(3), INTENT(IN)                  :: glb, lb, ub

      REAL(KIND=dp), PARAMETER :: exclusion_radius_angstrom = 0.25_dp

      INTEGER                                            :: i, iatom, j, k, natom
      REAL(KIND=dp)                                      :: atom_distance, exclusion_radius
      REAL(KIND=dp), DIMENSION(3)                        :: atom_displacement, grid_point

      natom = SIZE(particle_set)

      exclusion_radius = exclusion_radius_angstrom/angstrom

      ! Prevent solvent dielectric occupation within the atom-centred exclusion radius.
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(atom_displacement,atom_distance,grid_point,i,iatom,j,k) &
!$OMP             SHARED(atom_exclusion_mask,cell,exclusion_radius,glb,lb,natom,particle_set,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               atom_exclusion_mask%array(i, j, k) = 1.0_dp
               grid_point = [(0.5_dp + REAL(i - glb(1), KIND=dp))*atom_exclusion_mask%pw_grid%dr(1), &
                             (0.5_dp + REAL(j - glb(2), KIND=dp))*atom_exclusion_mask%pw_grid%dr(2), &
                             (0.5_dp + REAL(k - glb(3), KIND=dp))*atom_exclusion_mask%pw_grid%dr(3)]
               DO iatom = 1, natom
                  atom_displacement = pbc(grid_point - particle_set(iatom)%r, cell)
                  atom_distance = SQRT(DOT_PRODUCT(atom_displacement, atom_displacement))
                  IF (atom_distance < exclusion_radius) THEN
                     atom_exclusion_mask%array(i, j, k) = -1.0_dp
                     EXIT
                  END IF
               END DO
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE build_aics_atom_exclusion_mask

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of Andreussi et al.
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \par History:
!>      - Creation (16.10.2013,MK)
!>      - Finite difference of isosurfaces implemented (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.1
! **************************************************************************************************
   SUBROUTINE andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, rho_max, &
                        rho_min)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, dq, dt, f, ln_rho_max, ln_rho_min, &
                                                            q, rho, t, x, y

      CALL timeset(routineN, handle)

      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
         dq = -f*q
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dt,i,j,k,rho,t,x,y) &
!$OMP             SHARED(deps_elec,dq,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(ln_rho_max,rho_elec,q,rho_max,rho_min)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rho_min) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  deps_elec%array(i, j, k) = 0.0_dp
               ELSE IF (rho <= rho_max) THEN
                  x = LOG(rho)
                  y = q*(ln_rho_max - x)
                  t = f*(y - SIN(y))
                  eps_elec%array(i, j, k) = EXP(t)
                  dt = dq*(1.0_dp - COS(y))
                  deps_elec%array(i, j, k) = eps_elec%array(i, j, k)*dt/rho
               ELSE
                  eps_elec%array(i, j, k) = 1.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE andreussi

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of Fattebert and Gygi
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param beta ...
!> \param rho_zero ...
!> \par History:
!>      - Creation (15.10.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE fattebert_gygi(rho_elec, eps_elec, deps_elec, epsilon_solvent, beta, &
                             rho_zero)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, beta, rho_zero

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fattebert_gygi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: df, f, p, q, rho, s, t, twobeta

      CALL timeset(routineN, handle)

      df = (1.0_dp - epsilon_solvent)/rho_zero
      f = 0.5_dp*(epsilon_solvent - 1.0_dp)
      q = 1.0_dp/rho_zero
      twobeta = 2.0_dp*beta

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the smoothed dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,p,rho,s,t) &
!$OMP             SHARED(df,deps_elec,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(q,rho_elec,twobeta)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rhotol) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  deps_elec%array(i, j, k) = 0.0_dp
               ELSE
                  s = rho*q
                  p = s**twobeta
                  t = 1.0_dp/(1.0_dp + p)
                  eps_elec%array(i, j, k) = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
                  deps_elec%array(i, j, k) = df*twobeta*t*t*p/s
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE fattebert_gygi

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of solvent aware algorithm
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \param s ...
!> \param d_s_rhoel ...
!> \author     Ziwei Chai
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE sa_andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, rho_max, &
                           rho_min, s, d_s_rhoel)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: s, d_s_rhoel

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sa_andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, dq, dt, f, ln_rho_max, ln_rho_min, &
                                                            q, rho, t, x, y

      CALL timeset(routineN, handle)

      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
         dq = -f*q
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dt,i,j,k,rho,t,x,y) &
!$OMP             SHARED(deps_elec,dq,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(ln_rho_max,rho_elec,q,rho_max,rho_min, s, d_s_rhoel)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rho_min) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  s%array(i, j, k) = 0.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
                  d_s_rhoel%array(i, j, k) = 0.0_dp
               ELSE IF (rho <= rho_max) THEN
                  x = LOG(rho)
                  y = q*(ln_rho_max - x)
                  t = f*(y - SIN(y))
                  s%array(i, j, k) = 1 - t/(f*twopi)
                  eps_elec%array(i, j, k) = EXP(t)
                  dt = dq*(1.0_dp - COS(y))
                  deps_elec%array(i, j, k) = eps_elec%array(i, j, k)*dt/rho
                  d_s_rhoel%array(i, j, k) = deps_elec%array(i, j, k)/(-eps_elec%array(i, j, k)*LOG(epsilon_solvent))
               ELSE
                  eps_elec%array(i, j, k) = 1.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
                  s%array(i, j, k) = 1.0_dp
                  d_s_rhoel%array(i, j, k) = 0.0_dp
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE sa_andreussi

! **************************************************************************************************
!> \brief      Build the numerical derivative of a function on realspace grid
!> \param f ...
!> \param df ...
!> \param method ...
!> \param pw_env ...
!> \param input ...
!> \par History:
!>      - Creation (15.11.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE derive(f, df, method, pw_env, input)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: f
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: df
      INTEGER, INTENT(IN)                                :: method
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'derive'

      INTEGER                                            :: border_points, handle, i
      INTEGER, DIMENSION(3)                              :: lb, n, ub
      TYPE(pw_c1d_gs_type), DIMENSION(2)                 :: work_g1d
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_desc_type), POINTER            :: rs_desc
      TYPE(realspace_grid_input_type)                    :: input_settings
      TYPE(realspace_grid_type), POINTER                 :: rs_grid
      TYPE(section_vals_type), POINTER                   :: rs_grid_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(pw_env))

      ! Perform method specific setup
      SELECT CASE (method)
      CASE (sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7)
         NULLIFY (rs_desc)
         rs_grid_section => section_vals_get_subs_vals(input, "DFT%MGRID%RS_GRID")
         SELECT CASE (method)
         CASE (sccs_derivative_cd3)
            border_points = 1
         CASE (sccs_derivative_cd5)
            border_points = 2
         CASE (sccs_derivative_cd7)
            border_points = 3
         END SELECT
         CALL init_input_type(input_settings, 2*border_points + 1, rs_grid_section, &
                              1, [-1, -1, -1])
         CALL rs_grid_create_descriptor(rs_desc, f%pw_grid, input_settings, &
                                        border_points=border_points)
         ALLOCATE (rs_grid)
         CALL rs_grid_create(rs_grid, rs_desc)
!MK      CALL rs_grid_print(rs_grid, 6)
      CASE (sccs_derivative_fft)
         lb(1:3) = f%pw_grid%bounds_local(1, 1:3)
         ub(1:3) = f%pw_grid%bounds_local(2, 1:3)
         NULLIFY (auxbas_pw_pool)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         ! Get work storage for the 1d grids in g-space (derivative calculation)
         DO i = 1, SIZE(work_g1d)
            CALL auxbas_pw_pool%create_pw(work_g1d(i))
         END DO
      END SELECT

      ! Calculate the derivatives
      SELECT CASE (method)
      CASE (sccs_derivative_cd3)
         CALL derive_fdm_cd3(f, df, rs_grid)
      CASE (sccs_derivative_cd5)
         CALL derive_fdm_cd5(f, df, rs_grid)
      CASE (sccs_derivative_cd7)
         CALL derive_fdm_cd7(f, df, rs_grid)
      CASE (sccs_derivative_fft)
         ! FFT
         CALL pw_transfer(f, work_g1d(1))
         DO i = 1, 3
            n(:) = 0
            n(i) = 1
            CALL pw_copy(work_g1d(1), work_g1d(2))
            CALL pw_derive(work_g1d(2), n(:))
            CALL pw_transfer(work_g1d(2), df(i))
         END DO
      CASE DEFAULT
         CPABORT("Invalid derivative method for SCCS specified")
      END SELECT

      ! Perform method specific cleanup
      SELECT CASE (method)
      CASE (sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7)
         CALL rs_grid_release(rs_grid)
         DEALLOCATE (rs_grid)
         CALL rs_grid_release_descriptor(rs_desc)
      CASE (sccs_derivative_fft)
         DO i = 1, SIZE(work_g1d)
            CALL auxbas_pw_pool%give_back_pw(work_g1d(i))
         END DO
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE derive

! **************************************************************************************************
!> \brief      Calculate the finite difference between two isosurfaces of the
!>             electronic density. The smoothed dielectric function of
!>             Andreussi et al. is used as switching function eventually
!>             defining the quantum volume and surface of the cavity.
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param dtheta ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \param delta_rho ...
!> \par History:
!>      - Creation (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE surface_andreussi(rho_elec, norm_drho_elec, dtheta, &
                                epsilon_solvent, rho_max, rho_min, delta_rho)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, norm_drho_elec, dtheta
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min, &
                                                            delta_rho

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'surface_andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k, l
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, e, eps_elec, f, ln_rho_max, &
                                                            ln_rho_min, q, rho, t, x, y
      REAL(KIND=dp), DIMENSION(2)                        :: theta

      CALL timeset(routineN, handle)

      e = epsilon_solvent - 1.0_dp
      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate finite difference between two isosurfaces
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(eps_elec,i,j,k,l,rho,t,theta,x,y) &
!$OMP             SHARED(delta_rho,dtheta,e,epsilon_solvent,f,lb) &
!$OMP             SHARED(ln_rho_max,norm_drho_elec,rho_elec,q,rho_max,rho_min,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               DO l = 1, 2
                  rho = rho_elec%array(i, j, k) + (REAL(l, KIND=dp) - 1.5_dp)*delta_rho
                  IF (rho < rho_min) THEN
                     eps_elec = epsilon_solvent
                  ELSE IF (rho <= rho_max) THEN
                     x = LOG(rho)
                     y = q*(ln_rho_max - x)
                     t = f*(y - SIN(y))
                     eps_elec = EXP(t)
                  ELSE
                     eps_elec = 1.0_dp
                  END IF
                  theta(l) = (epsilon_solvent - eps_elec)/e
               END DO
               dtheta%array(i, j, k) = (theta(2) - theta(1))*norm_drho_elec%array(i, j, k)/delta_rho
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE surface_andreussi

! **************************************************************************************************
!> \brief      Calculate the finite difference between two isosurfaces of the
!>             the electronic density. The smoothed dielectric function of
!>             Fattebert and Gygi is used as switching function eventually
!>             defining the quantum volume and surface of the cavity.
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param dtheta ...
!> \param epsilon_solvent ...
!> \param beta ...
!> \param rho_zero ...
!> \param delta_rho ...
!> \par History:
!>      - Creation (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE surface_fattebert_gygi(rho_elec, norm_drho_elec, dtheta, &
                                     epsilon_solvent, beta, rho_zero, delta_rho)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, norm_drho_elec, dtheta
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, beta, rho_zero, &
                                                            delta_rho

      CHARACTER(LEN=*), PARAMETER :: routineN = 'surface_fattebert_gygi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k, l
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: e, eps_elec, f, p, q, rho, s, t, twobeta
      REAL(KIND=dp), DIMENSION(2)                        :: theta

      CALL timeset(routineN, handle)

      e = epsilon_solvent - 1.0_dp
      f = 0.5_dp*e
      q = 1.0_dp/rho_zero
      twobeta = 2.0_dp*beta

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate finite difference between two isosurfaces
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(eps_elec,i,j,k,l,p,rho,s,t,theta) &
!$OMP             SHARED(delta_rho,dtheta,e,epsilon_solvent,f,lb) &
!$OMP             SHARED(norm_drho_elec,q,rho_elec,twobeta,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               DO l = 1, 2
                  rho = rho_elec%array(i, j, k) + (REAL(l, KIND=dp) - 1.5_dp)*delta_rho
                  IF (rho < rhotol) THEN
                     eps_elec = epsilon_solvent
                  ELSE
                     s = rho*q
                     p = s**twobeta
                     t = 1.0_dp/(1.0_dp + p)
                     eps_elec = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
                  END IF
                  theta(l) = (epsilon_solvent - eps_elec)/e
               END DO
               dtheta%array(i, j, k) = (theta(2) - theta(1))*norm_drho_elec%array(i, j, k)/delta_rho
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE surface_fattebert_gygi

! **************************************************************************************************
!> \brief      Print SCCS results
!> \param energy ...
!> \param sccs_control ...
!> \param output_unit ...
!> \par History:
!>      - Creation (11.11.2022,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE print_sccs_results(energy, sccs_control, output_unit)

      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      INTEGER, INTENT(IN)                                :: output_unit

      IF (output_unit > 0) THEN
         CPASSERT(ASSOCIATED(energy))
         CPASSERT(ASSOCIATED(sccs_control))
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
            "SCCS| Hartree energy of solute and solvent [Hartree]", energy%sccs_hartree, &
            "SCCS| Hartree energy of the solute only    [Hartree]", energy%hartree
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14,/,T3,A,T61,F20.3)") &
            "SCCS| Polarisation energy                  [Hartree]", energy%sccs_pol, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_pol, "kcalmol"), &
            "SCCS| Cavitation energy                    [Hartree]", energy%sccs_cav, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_cav, "kcalmol"), &
            "SCCS| Dispersion free energy               [Hartree]", energy%sccs_dis, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_dis, "kcalmol"), &
            "SCCS| Repulsion free energy                [Hartree]", energy%sccs_rep, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_rep, "kcalmol"), &
            "SCCS| Solvation free energy                [Hartree]", energy%sccs_sol, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_sol, "kcalmol")
      END IF

   END SUBROUTINE print_sccs_results

END MODULE qs_sccs
