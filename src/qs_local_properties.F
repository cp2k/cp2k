!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for calculating local energy and stress tensor
!> \author JGH
!> \par History
!>      - 07.2019 created
! **************************************************************************************************
MODULE qs_local_properties
   USE bibliography,                    ONLY: Cohen2000,&
                                              Filippetti2000,&
                                              Rogers2002,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: det_3x3
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_integrate_function,&
                                              pw_multiply,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_core_energies,                ONLY: calculate_ptrace
   USE qs_energy_matrix_w,              ONLY: qs_energies_compute_w
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_xc_density
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_local_properties'

   PUBLIC :: qs_local_energy, qs_local_stress

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Routine to calculate the local energy
!> \param qs_env the qs_env to update
!> \param energy_density ...
!> \par History
!>      07.2019 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_local_energy(qs_env, energy_density)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), INTENT(INOUT)                     :: energy_density

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_local_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, img, iounit, ispin, nimages, &
                                                            nkind, nspins
      LOGICAL                                            :: gapw, gapw_xc
      REAL(KIND=dp)                                      :: eban, eband, eh, exc, ovol
      REAL(KIND=dp), DIMENSION(2)                        :: band_energy
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_w, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: band_density, edens_g, edens_r, &
                                                            hartree_density, rho_tot_gspace, &
                                                            rho_tot_rspace, v_hartree_rspace, &
                                                            xc_density
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_p_type), POINTER                           :: rho_core
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_struct
      TYPE(section_vals_type), POINTER                   :: input, xc_section

      CALL timeset(routineN, handle)

      CALL cite_reference(Cohen2000)

      CPASSERT(ASSOCIATED(qs_env))
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit()

      ! Check for GAPW method : additional terms for local densities
      CALL get_qs_env(qs_env, nkind=nkind, dft_control=dft_control)
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      nimages = dft_control%nimages
      nspins = dft_control%nspins

      ! get working arrays
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool, band_density%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, hartree_density%pw, use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, xc_density%pw, use_data=REALDATA3D, in_space=REALSPACE)

      ! w matrix
      CALL get_qs_env(qs_env, matrix_w_kp=matrix_w)
      IF (.NOT. ASSOCIATED(matrix_w)) THEN
         CALL get_qs_env(qs_env, &
                         ks_env=ks_env, &
                         matrix_s_kp=matrix_s)
         matrix => matrix_s(1, 1)%matrix
         CALL dbcsr_allocate_matrix_set(matrix_w, nspins, nimages)
         DO ispin = 1, nspins
            DO img = 1, nimages
               ALLOCATE (matrix_w(ispin, img)%matrix)
               CALL dbcsr_copy(matrix_w(ispin, img)%matrix, matrix, name="W MATRIX")
               CALL dbcsr_set(matrix_w(ispin, img)%matrix, 0.0_dp)
            END DO
         END DO
         CALL set_ks_env(ks_env, matrix_w_kp=matrix_w)
      END IF
      ! band structure energy density
      CALL qs_energies_compute_w(qs_env, .TRUE.)
      band_energy = 0.0_dp
      CALL get_qs_env(qs_env, ks_env=ks_env, matrix_w_kp=matrix_w)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             edens_r%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             edens_g%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_zero(band_density%pw)
      DO ispin = 1, nspins
         rho_ao => matrix_w(ispin, :)
         CALL calculate_rho_elec(matrix_p_kp=rho_ao, &
                                 rho=edens_r, &
                                 rho_gspace=edens_g, &
                                 total_rho=band_energy(ispin), &
                                 ks_env=ks_env, soft_valid=(gapw .OR. gapw_xc))
         CALL pw_axpy(edens_r%pw, band_density%pw)
      END DO
      CALL pw_pool_give_back_pw(auxbas_pw_pool, edens_r%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, edens_g%pw)

      ! Hartree energy density correction = -0.5 * V_H(r) * [rho(r) - rho_core(r)]
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             rho_tot_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      NULLIFY(rho_core)
      CALL get_qs_env(qs_env, &
                      v_hartree_rspace=v_hartree_rspace%pw, &
                      rho_core=rho_core, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r)
      IF(ASSOCIATED(rho_core)) THEN
         CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)
         CALL pw_transfer(rho_core%pw, rho_tot_rspace%pw)
      ELSE
         CALL pw_zero(rho_tot_rspace%pw)
      END IF
      DO ispin = 1, nspins
         CALL pw_axpy(rho_r(ispin)%pw, rho_tot_rspace%pw, alpha=-1.0_dp)
      END DO
      CALL pw_zero(hartree_density%pw)
      ovol = 0.5_dp/hartree_density%pw%pw_grid%dvol
      CALL pw_multiply(hartree_density%pw, v_hartree_rspace%pw, rho_tot_rspace%pw, alpha=ovol)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_rspace%pw)

      IF (dft_control%do_admm) THEN
         CALL cp_warn(__LOCATION__, "ADMM not supported for local energy calculation")
      END IF
      IF (gapw_xc .OR. gapw) THEN
         CALL cp_warn(__LOCATION__, "GAPW/GAPW_XC not supported for local energy calculation")
      END IF
      ! XC energy density correction = E_xc(r) - V_xc(r)*rho(r)
      CALL get_qs_env(qs_env, input=input)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      CALL get_qs_env(qs_env=qs_env, rho=rho_struct)
      !
      CALL qs_xc_density(ks_env, rho_struct, xc_section, xc_density)
      !
      ! energies
      CALL get_qs_env(qs_env=qs_env, energy=energy)
      eban = pw_integrate_function(band_density%pw)
      eh = pw_integrate_function(hartree_density%pw)
      exc = pw_integrate_function(xc_density%pw)

      ! band energy
      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
      CALL calculate_ptrace(matrix_ks, rho_ao_kp, eband, nspins)

      ! get full density
      CALL pw_copy(band_density%pw, energy_density%pw)
      CALL pw_axpy(hartree_density%pw, energy_density%pw)
      CALL pw_axpy(xc_density%pw, energy_density%pw)

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T3,A)") REPEAT("=", 78)
         WRITE (UNIT=iounit, FMT="(T4,A,T52,A,T75,A)") "Local Energy Calculation", "GPW/GAPW", "Local"
         WRITE (UNIT=iounit, FMT="(T4,A,T45,F15.8,T65,F15.8)") "Band Energy", eband, eban
         WRITE (UNIT=iounit, FMT="(T4,A,T65,F15.8)") "Hartree Energy Correction", eh
         WRITE (UNIT=iounit, FMT="(T4,A,T65,F15.8)") "XC Energy Correction", exc
         WRITE (UNIT=iounit, FMT="(T4,A,T45,F15.8,T65,F15.8)") "Total Energy", &
            energy%total, eban + eh + exc + energy%core_overlap + energy%core_self + energy%dispersion
         WRITE (UNIT=iounit, FMT="(T3,A)") REPEAT("=", 78)
      END IF

      ! return temp arrays
      CALL pw_pool_give_back_pw(auxbas_pw_pool, band_density%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, hartree_density%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, xc_density%pw)

      CALL timestop(handle)

   END SUBROUTINE qs_local_energy

! **************************************************************************************************
!> \brief Routine to calculate the local stress
!> \param qs_env the qs_env to update
!> \param stress_tensor ...
!> \param pressure ...
!> \param beta ...
!> \par History
!>      07.2019 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_local_stress(qs_env, stress_tensor, pressure, beta)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: stress_tensor
      TYPE(pw_p_type), INTENT(INOUT), OPTIONAL           :: pressure
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_local_stress', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iounit, j, nimages, nkind, &
                                                            nspins
      LOGICAL                                            :: do_pressure, do_stress, gapw, gapw_xc, &
                                                            use_virial
      REAL(KIND=dp)                                      :: my_beta
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_loc
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: v_hartree_gspace, v_hartree_rspace, &
                                                            xc_density
      TYPE(pw_p_type), DIMENSION(3)                      :: efield
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(section_vals_type), POINTER                   :: input, xc_section
      TYPE(virial_type), POINTER                         :: virial

      CALL cp_warn(__LOCATION__, "Local Stress Tensor code is not working, skipping")
      RETURN

      CALL timeset(routineN, handle)

      CALL cite_reference(Filippetti2000)
      CALL cite_reference(Rogers2002)

      CPASSERT(ASSOCIATED(qs_env))

      IF (PRESENT(stress_tensor)) THEN
         do_stress = .TRUE.
      ELSE
         do_stress = .FALSE.
      END IF
      IF (PRESENT(pressure)) THEN
         do_pressure = .TRUE.
      ELSE
         do_pressure = .FALSE.
      END IF
      IF (PRESENT(beta)) THEN
         my_beta = beta
      ELSE
         my_beta = 0.0_dp
      END IF

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit()

      !!!!!!
      CALL cp_warn(__LOCATION__, "Local Stress Tensor code is not tested")
      !!!!!!

      ! Check for GAPW method : additional terms for local densities
      CALL get_qs_env(qs_env, nkind=nkind, dft_control=dft_control)
      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc

      nimages = dft_control%nimages
      nspins = dft_control%nspins

      ! get working arrays
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool, xc_density%pw, use_data=REALDATA3D, in_space=REALSPACE)

      ! init local stress tensor
      IF (do_stress) THEN
         DO i = 1, 3
            DO j = 1, 3
               CALL pw_zero(stress_tensor(i, j)%pw)
            END DO
         END DO
      END IF

      IF (dft_control%do_admm) THEN
         CALL cp_warn(__LOCATION__, "ADMM not supported for local energy calculation")
      END IF
      IF (gapw_xc .OR. gapw) THEN
         CALL cp_warn(__LOCATION__, "GAPW/GAPW_XC not supported for local energy calculation")
      END IF
      ! XC energy density correction = E_xc(r) - V_xc(r)*rho(r)
      CALL get_qs_env(qs_env, ks_env=ks_env, input=input, rho=rho_struct)
      xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      !
      CALL qs_xc_density(ks_env, rho_struct, xc_section, xc_density)

      ! Electrical field terms
      CALL get_qs_env(qs_env, v_hartree_rspace=v_hartree_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_transfer(v_hartree_rspace%pw, v_hartree_gspace%pw)
      DO i = 1, 3
         NULLIFY (efield(i)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, efield(i)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         CALL pw_copy(v_hartree_gspace%pw, efield(i)%pw)
      END DO
      CALL pw_derive(efield(1)%pw, (/1, 0, 0/))
      CALL pw_derive(efield(2)%pw, (/0, 1, 0/))
      CALL pw_derive(efield(3)%pw, (/0, 0, 1/))
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      DO i = 1, 3
         CALL pw_pool_give_back_pw(auxbas_pw_pool, efield(i)%pw)
      END DO

      pv_loc = 0.0_dp

      CALL get_qs_env(qs_env=qs_env, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      IF (.NOT. use_virial) THEN
         CALL cp_warn(__LOCATION__, "Local stress should be used with standard stress calculation.")
      END IF
      IF (iounit > 0 .AND. use_virial) THEN
         WRITE (UNIT=iounit, FMT="(/,T3,A)") REPEAT("=", 78)
         WRITE (UNIT=iounit, FMT="(T4,A)") "Local Stress Calculation"
         WRITE (UNIT=iounit, FMT="(T42,A,T64,A)") "       1/3 Trace", "     Determinant"
         WRITE (UNIT=iounit, FMT="(T4,A,T42,F16.8,T64,F16.8)") "Total Stress", &
            (pv_loc(1, 1) + pv_loc(2, 2) + pv_loc(3, 3))/3.0_dp, det_3x3(pv_loc)
         WRITE (UNIT=iounit, FMT="(T3,A)") REPEAT("=", 78)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_local_stress

! **************************************************************************************************

END MODULE qs_local_properties
