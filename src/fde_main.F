!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE fde_main
   USE bibliography,                    ONLY: cite_reference,&
                                              Schreder2024_3
   USE cp_control_types,                ONLY: dft_control_type,&
                                              fde_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE fde_methods,                     ONLY: fde_build_vemb_A,&
                                              fde_read_pw,&
                                              fde_read_results,&
                                              fde_spawn,&
                                              fde_v_es,&
                                              fde_write_pw
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integral_ab,&
                                              pw_scale,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_external_potential,           ONLY: interpolate_external_potential
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fde_main'

   PUBLIC :: fde_run

CONTAINS

! **************************************************************************************************
!> \brief Join a working directory and a file name into a path.
!> \param workdir ...
!> \param filename ...
!> \return ...
! **************************************************************************************************
   FUNCTION fde_join_path(workdir, filename) RESULT(path)
      CHARACTER(len=*), INTENT(IN)       :: workdir, filename
      CHARACTER(len=default_path_length) :: path

      IF (LEN_TRIM(workdir) > 0) THEN
         path = TRIM(workdir) // "/" // TRIM(filename)
      ELSE
         path = TRIM(filename)
      END IF

   END FUNCTION fde_join_path

! **************************************************************************************************
!> \brief Integrate a real-space pw
!> \param pw ...
!> \param para_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION fde_integrate_pw(pw, para_env) RESULT(res)
      TYPE(pw_r3d_rs_type), INTENT(IN)   :: pw
      TYPE(mp_para_env_type), INTENT(IN) :: para_env
      REAL(KIND=dp)                      :: res

      res = SUM(pw%array) * pw%pw_grid%dvol
      CALL para_env%sum(res)

   END FUNCTION fde_integrate_pw

! **************************************************************************************************
!> \brief Run the FDE self-consistency loop
!> \param qs_env ..
! **************************************************************************************************
   SUBROUTINE fde_run(qs_env)
      TYPE(qs_environment_type), POINTER          :: qs_env

      CHARACTER(len=*), PARAMETER                 :: routineN = 'fde_scf'

      CHARACTER(len=default_path_length)          :: grid_path, input, v_path, res_path, rho_path
      INTEGER                                     :: handle, iter, iounit, max_iter
      INTEGER                                     :: i_a, n_a
      LOGICAL                                     :: converged
      REAL(KIND=dp)                               :: drho, e_sub, e_sub_bare, e_emb, e_es_e, &
                                                     e_es_int, e_es_n, e_env, e_kin_na, e_xc_na, &
                                                     e_n, e_tot, n_el, n_el_pw, phi_at_atom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)    :: z
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: r
      TYPE(cp_logger_type), POINTER               :: logger
      TYPE(dft_control_type), POINTER             :: dft_control
      TYPE(fde_control_type), POINTER             :: fde_control
      TYPE(mp_para_env_type), POINTER             :: para_env
      TYPE(pw_env_type), POINTER                  :: pw_env
      TYPE(pw_pool_type), POINTER                 :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type)                        :: rho, rho_new, v_emb
      TYPE(pw_r3d_rs_type), TARGET                :: phi_es
      TYPE(pw_r3d_rs_type), POINTER               :: phi_es_ptr
      TYPE(qs_energy_type), POINTER               :: qs_energy

      CALL timeset(routineN, handle)

      CALL cite_reference(Schreder2024_3)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      NULLIFY (dft_control, fde_control, para_env, pw_env, auxbas_pw_pool, qs_energy)
      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env, pw_env=pw_env, &
                      energy=qs_energy)
      fde_control => dft_control%fde_control
      CPASSERT(ASSOCIATED(fde_control))

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      ! socket files
      v_path = fde_join_path(fde_control%workdir, fde_control%pot_file)
      grid_path = fde_join_path(fde_control%workdir, fde_control%grid_file)
      rho_path = fde_join_path(fde_control%workdir, fde_control%rho_file)
      res_path = fde_join_path(fde_control%workdir, fde_control%result_file)
      input = fde_control%molcas_input_template

      ! bare environment DFT energy
      e_env = qs_energy%total

      ! initial guess = zero
      CALL auxbas_pw_pool%create_pw(rho)
      CALL pw_zero(rho)
      CALL auxbas_pw_pool%create_pw(rho_new)

      max_iter = 1
      IF (fde_control%freeze_and_thaw) max_iter = fde_control%max_ft_iter

      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T2,A)') REPEAT("=", 78)
         WRITE (iounit, '(T2,A)') "FDE:  Frozen-density embedding"
         WRITE (iounit, '(T2,A)') REPEAT("=", 78)
      END IF

      converged = .FALSE.
      DO iter = 1, max_iter
         CALL fde_build_vemb_A(qs_env, rho, v_emb, e_xc_na, e_kin_na)
         CALL fde_write_pw(v_emb, TRIM(v_path), write_value=.TRUE.)
         CALL fde_write_pw(v_emb, TRIM(grid_path), write_value=.FALSE.)

         CALL fde_spawn(fde_control, TRIM(input), para_env)

         CALL fde_read_pw(rho_new, TRIM(rho_path))
         CALL fde_read_results(TRIM(res_path), para_env, e_sub, e_emb, e_n, n_el, n_a, z, r)

         n_el_pw = fde_integrate_pw(rho_new, para_env)
         IF (n_el_pw > 1.0e-10_dp .AND. n_el > 0.0_dp) &
            CALL pw_scale(rho_new, n_el/n_el_pw)

         drho = fde_drho(rho, rho_new, para_env)

         CALL fde_v_es(qs_env, phi_es)
         phi_es_ptr => phi_es
         e_es_e = pw_integral_ab(rho_new, phi_es)

         e_es_n = 0.0_dp
         DO i_a = 1, n_a
            CALL interpolate_external_potential(r(:, i_a), phi_es_ptr, func=phi_at_atom)
            e_es_n = e_es_n - z(i_a)*phi_at_atom
         END DO
         e_es_int = e_es_e + e_es_n
         CALL auxbas_pw_pool%give_back_pw(phi_es)
         DEALLOCATE (z, r)

         e_sub_bare = e_sub - e_emb
         e_tot = e_sub_bare + e_env + e_es_int + e_xc_na + e_kin_na

         IF (iounit > 0) THEN
            WRITE (iounit, '(/,T2,A,I0)') "FDE freeze-and-thaw iteration ", iter
            WRITE (iounit, '(T4,A,F20.10)') "E(sub, embedded)       = ", e_sub
            WRITE (iounit, '(T4,A,F20.10)') "E(emb)                 = ", e_emb
            WRITE (iounit, '(T4,A,F20.10)') "E(sub)                 = ", e_sub_bare
            WRITE (iounit, '(T4,A,F20.10)') "E(env)                 = ", e_env
            WRITE (iounit, '(T4,A,F20.10)') "E(es)                  = ", e_es_int
            WRITE (iounit, '(T6,A,F20.10)') "E(elec)                = ", e_es_e
            WRITE (iounit, '(T6,A,F20.10)') "E(nuc)                 = ", e_es_n
            WRITE (iounit, '(T4,A,F20.10)') "E(XC, nonadditive)     = ", e_xc_na
            WRITE (iounit, '(T4,A,F20.10)') "E(kin, nonadditive)    = ", e_kin_na
            WRITE (iounit, '(T4,A,F20.10)') "E(FDE)                 = ", e_tot
            WRITE (iounit, '(T4,A,F14.6,A,F14.6)') "N(sub, e): grid = ", n_el_pw, &
                                                   " reported = ", n_el
            WRITE (iounit, '(T4,A,ES14.6)') "int|rho(sub, new) - rho(sub, old)|dr = ", drho
         END IF

         CALL auxbas_pw_pool%give_back_pw(v_emb)

         IF (fde_control%freeze_and_thaw) THEN
            IF (drho < fde_control%eps_ft_density) THEN
               converged = .TRUE.
            ELSE
               CALL pw_scale(rho, 1.0_dp - fde_control%density_mixing)
               CALL pw_axpy(rho_new, rho, fde_control%density_mixing)
            END IF
         ELSE
            converged = .TRUE.
         END IF

         IF (converged) EXIT

      END DO

      IF (iounit > 0) THEN
         WRITE (iounit, '(/,T2,A)') REPEAT("=", 78)
         IF (fde_control%freeze_and_thaw .AND. .NOT. converged) &
            WRITE (iounit, '(T2,A)') "FDE  WARNING: freeze-and-thaw did not converge"
         WRITE (iounit, '(T2,A,F20.10)') "FDE final total energy = ", e_tot
         WRITE (iounit, '(T2,A)') REPEAT("=", 78)
      END IF

      qs_energy%total = e_tot

      CALL auxbas_pw_pool%give_back_pw(rho)
      CALL auxbas_pw_pool%give_back_pw(rho_new)

      CALL timestop(handle)

   END SUBROUTINE fde_run

! **************************************************************************************************
!> \brief Integrated density difference: int |rho_new - rho_old| dr
!> \param rho_a ...
!> \param rho_b ...
!> \param para_env ...
!> \return ...
! **************************************************************************************************
   FUNCTION fde_drho(rho_a, rho_b, para_env) RESULT(drho)
      TYPE(pw_r3d_rs_type), INTENT(IN)   :: rho_a, rho_b
      TYPE(mp_para_env_type), INTENT(IN) :: para_env
      REAL(KIND=dp)                      :: drho

      drho = SUM(ABS(rho_b%array - rho_a%array)) * rho_b%pw_grid%dvol
      CALL para_env%sum(drho)

   END FUNCTION fde_drho

END MODULE fde_main
