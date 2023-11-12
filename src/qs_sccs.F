!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
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
   USE input_constants,                 ONLY: sccs_andreussi,&
                                              sccs_derivative_cd3,&
                                              sccs_derivative_cd5,&
                                              sccs_derivative_cd7,&
                                              sccs_derivative_fft,&
                                              sccs_fattebert_gygi
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE mathconstants,                   ONLY: fourpi,&
                                              twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_list_types,             ONLY: particle_list_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_integral_ab,&
                                              pw_integrate_function,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_analytic,&
                                              pw_poisson_mt,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
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

   PUBLIC :: print_sccs_results, sccs

CONTAINS

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
      TYPE(pw_type), INTENT(INOUT)                       :: rho_tot_gspace, v_hartree_gspace, v_sccs
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: h_stress

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sccs'
      REAL(KIND=dp), PARAMETER                           :: epstol = 1.0E-8_dp

      CHARACTER(LEN=4*default_string_length)             :: message
      CHARACTER(LEN=default_path_length)                 :: mpi_filename
      CHARACTER(LEN=default_string_length)               :: cube_path, filename, my_pos_cube, &
                                                            print_path
      INTEGER                                            :: cube_unit, handle, i, ispin, iter, j, k, &
                                                            nspin, output_unit, print_level
      INTEGER(KIND=int_8)                                :: ngpts
      INTEGER, DIMENSION(3)                              :: lb, ub
      LOGICAL                                            :: append_cube, calculate_stress_tensor, &
                                                            mpi_io, should_output
      REAL(KIND=dp) :: cavity_surface, cavity_volume, cell_volume, dphi2, dvol, e_tot, &
         epsilon_solvent, f, polarisation_charge, rho_delta, rho_delta_avg, rho_delta_max, &
         rho_iter_new, tot_rho_elec, tot_rho_solute
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type)                                      :: deps_elec, eps_elec
      TYPE(pw_type), DIMENSION(3)                        :: dln_eps_elec, dphi_tot
      TYPE(pw_type), DIMENSION(:), POINTER               :: rho_pw_r
      TYPE(pw_type), POINTER                             :: rho_pw_r_sccs
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
         TYPE(pw_type) :: rho_elec
         CALL auxbas_pw_pool%create_pw(rho_elec, &
                                       use_data=REALDATA3D, &
                                       in_space=REALSPACE)

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
         CALL auxbas_pw_pool%create_pw(eps_elec, &
                                       use_data=REALDATA3D, &
                                       in_space=REALSPACE)
         CALL auxbas_pw_pool%create_pw(deps_elec, &
                                       use_data=REALDATA3D, &
                                       in_space=REALSPACE)
         ! Relative permittivity or dielectric constant of the solvent (medium)
         epsilon_solvent = sccs_control%epsilon_solvent
         SELECT CASE (sccs_control%method_id)
         CASE (sccs_andreussi)
            CALL andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%rho_max, &
                           sccs_control%rho_min)
         CASE (sccs_fattebert_gygi)
            CALL fattebert_gygi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%beta, &
                                sccs_control%rho_zero)
         CASE DEFAULT
            CPABORT("Invalid method specified for SCCS model")
         END SELECT

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

         ! Calculate the (quantum) volume and surface of the solute cavity
         cavity_surface = 0.0_dp
         cavity_volume = 0.0_dp

         IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN

            BLOCK
               TYPE(pw_type) :: theta, norm_drho_elec
               TYPE(pw_type), DIMENSION(3)                        :: drho_elec
               CALL auxbas_pw_pool%create_pw(theta, &
                                             use_data=REALDATA3D, &
                                             in_space=REALSPACE)
               CALL pw_zero(theta)

               ! Calculate the (quantum) volume of the solute cavity
               f = 1.0_dp/(epsilon_solvent - 1.0_dp)
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(epsilon_solvent,eps_elec,f,lb,theta,ub)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        theta%cr3d(i, j, k) = f*(epsilon_solvent - eps_elec%cr3d(i, j, k))
                     END DO
                  END DO
               END DO
!$OMP    END PARALLEL DO
               cavity_volume = pw_integrate_function(theta)

               ! Calculate the derivative of the electronic density in r-space
               ! TODO: Could be retrieved from the QS environment
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_elec(i), &
                                                use_data=REALDATA3D, &
                                                in_space=REALSPACE)
               END DO
               CALL derive(rho_elec, drho_elec, sccs_derivative_fft, pw_env, input)

               CALL auxbas_pw_pool%create_pw(norm_drho_elec, &
                                             use_data=REALDATA3D, &
                                             in_space=REALSPACE)

               ! Calculate the norm of the gradient of the electronic density in r-space
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(drho_elec,lb,norm_drho_elec,ub)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        norm_drho_elec%cr3d(i, j, k) = SQRT(drho_elec(1)%cr3d(i, j, k)* &
                                                            drho_elec(1)%cr3d(i, j, k) + &
                                                            drho_elec(2)%cr3d(i, j, k)* &
                                                            drho_elec(2)%cr3d(i, j, k) + &
                                                            drho_elec(3)%cr3d(i, j, k)* &
                                                            drho_elec(3)%cr3d(i, j, k))
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
         TYPE(pw_type) :: rho_tot, phi_tot, rho_solute, rho_tot_zero
         ! Retrieve the total charge density (core + elec) of the solute in r-space
         CALL auxbas_pw_pool%create_pw(rho_solute, &
                                       use_data=REALDATA3D, &
                                       in_space=REALSPACE)
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
         CALL auxbas_pw_pool%create_pw(rho_tot_zero, &
                                       use_data=REALDATA3D, &
                                       in_space=REALSPACE)

         ! Build the initial (rho_iter = 0) total charge density (solute plus polarisation) in r-space
         ! eps_elec <- ln(eps_elec)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(eps_elec,lb,message,output_unit,para_env,ub) &
!$OMP             SHARED(rho_solute,rho_tot_zero)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  IF (eps_elec%cr3d(i, j, k) < 1.0_dp) THEN
                     WRITE (UNIT=message, FMT="(A,ES12.3,A,3(I0,A))") &
                        "SCCS| Invalid dielectric function value ", eps_elec%cr3d(i, j, k), &
                        " encountered at grid point (", i, ",", j, ",", k, ")"
                     CPABORT(message)
                  END IF
                  rho_tot_zero%cr3d(i, j, k) = rho_solute%cr3d(i, j, k)/eps_elec%cr3d(i, j, k)
                  eps_elec%cr3d(i, j, k) = LOG(eps_elec%cr3d(i, j, k))
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

         ! Build the derivative of LOG(eps_elec)
         DO i = 1, 3
            CALL auxbas_pw_pool%create_pw(dln_eps_elec(i), &
                                          use_data=REALDATA3D, &
                                          in_space=REALSPACE)
            CALL pw_zero(dln_eps_elec(i))
         END DO
         CALL derive(eps_elec, dln_eps_elec, sccs_control%derivative_method, pw_env, input)
         CALL auxbas_pw_pool%give_back_pw(eps_elec)

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
            CALL auxbas_pw_pool%create_pw(dphi_tot(i), &
                                          use_data=REALDATA3D, &
                                          in_space=REALSPACE)
         END DO

         ! Initialise the total charge density in r-space rho_tot with rho_tot_zero + rho_iter_zero
         CALL auxbas_pw_pool%create_pw(rho_tot, &
                                       use_data=REALDATA3D, &
                                       in_space=REALSPACE)
         CALL pw_copy(rho_tot_zero, rho_tot)
         CALL pw_axpy(rho_pw_r_sccs, rho_tot)

         CALL auxbas_pw_pool%create_pw(phi_tot, &
                                       use_data=REALDATA3D, &
                                       in_space=REALSPACE)
         CALL pw_zero(phi_tot)

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
                     rho_iter_new = (dln_eps_elec(1)%cr3d(i, j, k)*dphi_tot(1)%cr3d(i, j, k) + &
                                     dln_eps_elec(2)%cr3d(i, j, k)*dphi_tot(2)%cr3d(i, j, k) + &
                                     dln_eps_elec(3)%cr3d(i, j, k)*dphi_tot(3)%cr3d(i, j, k))*f
                     rho_iter_new = rho_pw_r_sccs%cr3d(i, j, k) + &
                                    sccs_control%mixing*(rho_iter_new - rho_pw_r_sccs%cr3d(i, j, k))
                     rho_delta = ABS(rho_iter_new - rho_pw_r_sccs%cr3d(i, j, k))
                     rho_delta_max = MAX(rho_delta, rho_delta_max)
                     rho_delta_avg = rho_delta_avg + rho_delta
                     rho_tot%cr3d(i, j, k) = rho_tot_zero%cr3d(i, j, k) + rho_iter_new
                     rho_pw_r_sccs%cr3d(i, j, k) = rho_iter_new
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
            TYPE(pw_type) :: phi_solute
            CALL auxbas_pw_pool%create_pw(phi_solute, &
                                          use_data=REALDATA3D, &
                                          in_space=REALSPACE)
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
      f = -0.5_dp*dvol/fourpi
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dphi2,i,j,k) &
!$OMP             SHARED(f,deps_elec,dphi_tot) &
!$OMP             SHARED(lb,ub,v_sccs)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               dphi2 = dphi_tot(1)%cr3d(i, j, k)*dphi_tot(1)%cr3d(i, j, k) + &
                       dphi_tot(2)%cr3d(i, j, k)*dphi_tot(2)%cr3d(i, j, k) + &
                       dphi_tot(3)%cr3d(i, j, k)*dphi_tot(3)%cr3d(i, j, k)
               v_sccs%cr3d(i, j, k) = v_sccs%cr3d(i, j, k) + f*deps_elec%cr3d(i, j, k)*dphi2
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

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

      TYPE(pw_type), INTENT(IN)                          :: rho_elec, eps_elec, deps_elec
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
               rho = rho_elec%cr3d(i, j, k)
               IF (rho < rho_min) THEN
                  eps_elec%cr3d(i, j, k) = epsilon_solvent
                  deps_elec%cr3d(i, j, k) = 0.0_dp
               ELSE IF (rho <= rho_max) THEN
                  x = LOG(rho)
                  y = q*(ln_rho_max - x)
                  t = f*(y - SIN(y))
                  eps_elec%cr3d(i, j, k) = EXP(t)
                  dt = dq*(1.0_dp - COS(y))
                  deps_elec%cr3d(i, j, k) = eps_elec%cr3d(i, j, k)*dt/rho
               ELSE
                  eps_elec%cr3d(i, j, k) = 1.0_dp
                  deps_elec%cr3d(i, j, k) = 0.0_dp
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

      TYPE(pw_type), INTENT(IN)                          :: rho_elec, eps_elec, deps_elec
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
               rho = rho_elec%cr3d(i, j, k)
               IF (rho < rhotol) THEN
                  eps_elec%cr3d(i, j, k) = epsilon_solvent
                  deps_elec%cr3d(i, j, k) = 0.0_dp
               ELSE
                  s = rho*q
                  p = s**twobeta
                  t = 1.0_dp/(1.0_dp + p)
                  eps_elec%cr3d(i, j, k) = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
                  deps_elec%cr3d(i, j, k) = df*twobeta*t*t*p/s
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE fattebert_gygi

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

      TYPE(pw_type), INTENT(IN)                          :: f
      TYPE(pw_type), DIMENSION(3), INTENT(INOUT)         :: df
      INTEGER, INTENT(IN)                                :: method
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'derive'

      INTEGER                                            :: border_points, handle, i
      INTEGER, DIMENSION(3)                              :: lb, n, ub
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), DIMENSION(2)                        :: work_g1d
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
                              1, (/-1, -1, -1/))
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
            CALL auxbas_pw_pool%create_pw(work_g1d(i), &
                                          use_data=COMPLEXDATA1D, &
                                          in_space=RECIPROCALSPACE)
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

      TYPE(pw_type), INTENT(IN)                          :: rho_elec, norm_drho_elec, dtheta
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
                  rho = rho_elec%cr3d(i, j, k) + (REAL(l, KIND=dp) - 1.5_dp)*delta_rho
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
               dtheta%cr3d(i, j, k) = (theta(2) - theta(1))*norm_drho_elec%cr3d(i, j, k)/delta_rho
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

      TYPE(pw_type), INTENT(IN)                          :: rho_elec, norm_drho_elec, dtheta
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
                  rho = rho_elec%cr3d(i, j, k) + (REAL(l, KIND=dp) - 1.5_dp)*delta_rho
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
               dtheta%cr3d(i, j, k) = (theta(2) - theta(1))*norm_drho_elec%cr3d(i, j, k)/delta_rho
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
