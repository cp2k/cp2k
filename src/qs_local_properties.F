!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for calculating local energy and stress tensor
!> \author JGH
!> \par History
!>      - 07.2019 created
! **************************************************************************************************
MODULE qs_local_properties
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE bibliography,                    ONLY: Cohen2000,&
                                              Rogers2002,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
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
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_atom_methods,             ONLY: allocate_rho_atom_internals,&
                                              calculate_rho_atom,&
                                              calculate_rho_atom_coeff
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_xc_density
   USE qs_vxc_atom,                     ONLY: calculate_vxc_atom
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_local_properties'

   PUBLIC :: qs_local_energy

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Routine to calculate the local energy
!> \param qs_env the qs_env to update
!> \param band_density ...
!> \param hartree_density ...
!> \param xc_density ...
!> \param eden_atom_set ...
!> \par History
!>      07.2019 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_local_energy(qs_env, band_density, hartree_density, xc_density, eden_atom_set)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), INTENT(INOUT)                     :: band_density, hartree_density, xc_density
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: eden_atom_set

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_local_energy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, img, iounit, ispin, na, &
                                                            nimages, nkind, nspins
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: gapw, gapw_xc, paw_atom
      REAL(KIND=dp)                                      :: eban, eband, eh, exc, ovol
      REAL(KIND=dp), DIMENSION(2)                        :: band_energy, eden_h, eden_s
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_w, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: edens_g, edens_r, rho_tot_gspace, &
                                                            rho_tot_rspace, v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_p_type), POINTER                           :: rho_core
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_struct, rho_xc
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
      CALL get_qs_env(qs_env, ks_env=ks_env, matrix_w_kp=matrix_w, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
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

      ! atomic energy densities
      IF (gapw .OR. gapw_xc) THEN
         CPASSERT(ASSOCIATED(eden_atom_set))
         CALL get_qs_env(qs_env, &
                         qs_kind_set=qs_kind_set, &
                         atomic_kind_set=atomic_kind_set, &
                         para_env=para_env)
         CALL allocate_rho_atom_internals(qs_env, eden_atom_set)
         CALL calculate_rho_atom_coeff(qs_env, matrix_w, eden_atom_set)
         eden_h = 0.0_dp
         eden_s = 0.0_dp
         DO ikind = 1, nkind
            NULLIFY (grid_atom, harmonics)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=na)
            CALL get_qs_kind(qs_kind_set(ikind), grid_atom=grid_atom, &
                             paw_atom=paw_atom, harmonics=harmonics)
            IF (paw_atom) THEN
               CALL calculate_rho_atom(para_env, eden_atom_set, qs_kind_set(ikind), &
                                       atom_list, grid_atom, na, nspins, eden_h, eden_s)
            ENDIF
         ENDDO
      END IF

      ! Hartree energy density correction = -0.5 * V_H(r) * [rho(r) - rho_core(r)]
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             rho_tot_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
      CALL get_qs_env(qs_env, &
                      v_hartree_rspace=v_hartree_rspace%pw, &
                      rho_core=rho_core, rho=rho)
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)
      CALL qs_rho_get(rho, rho_r=rho_r)
      CALL pw_transfer(rho_core%pw, rho_tot_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_axpy(rho_r(ispin)%pw, rho_tot_rspace%pw, alpha=-1.0_dp)
      END DO
      CALL pw_zero(hartree_density%pw)
      ovol = 0.5_dp/hartree_density%pw%pw_grid%dvol
      CALL pw_multiply(hartree_density%pw, v_hartree_rspace%pw, rho_tot_rspace%pw, alpha=ovol)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_rspace%pw)

      ! XC energy density correction = E_xc(r) - V_xc(r)*rho(r)
      IF (dft_control%do_admm) THEN
         CALL cp_warn(__LOCATION__, "ADMM not supported for local energy calculation")
      ELSE
         CALL get_qs_env(qs_env, input=input)
         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
      END IF
      IF (gapw_xc) THEN
         CALL get_qs_env(qs_env=qs_env, rho_xc=rho_struct)
         CPASSERT(ASSOCIATED(rho_xc))
      ELSE
         CALL get_qs_env(qs_env=qs_env, rho=rho_struct)
      END IF
      !
      CALL qs_xc_density(ks_env, rho_struct, xc_section, xc_density)
      !
      IF (gapw .OR. gapw_xc) THEN
         CALL calculate_vxc_atom(qs_env, .FALSE.)
      ENDIF

      ! energies
      CALL get_qs_env(qs_env=qs_env, energy=energy)
      eban = pw_integrate_function(band_density%pw)
      eh = pw_integrate_function(hartree_density%pw)
      exc = pw_integrate_function(xc_density%pw)

      ! band energy
      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
      CALL calculate_ptrace(matrix_ks, rho_ao_kp, eband, nspins)

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T3,A)") REPEAT("=", 78)
         WRITE (UNIT=iounit, FMT="(T4,A,T52,A,T75,A)") "Local Energy Calculation", "GPW/GAPW", "Local"
         WRITE (UNIT=iounit, FMT="(T4,A,T45,F15.8,T65,F15.8)") "Band Energy", eband, eban
         WRITE (UNIT=iounit, FMT="(T4,A,T65,F15.8)") "Hartree Energy Correction", eh
         WRITE (UNIT=iounit, FMT="(T4,A,T65,F15.8)") "XC Energy Correction", exc
         WRITE (UNIT=iounit, FMT="(T4,A,T45,F15.8,T65,F15.8)") "Total Energy", &
            energy%total, eban+eh+exc+energy%core_overlap+energy%core_self
         WRITE (UNIT=iounit, FMT="(T3,A)") REPEAT("=", 78)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_local_energy

! **************************************************************************************************
!> \brief Routine to calculate the local stress
!> \param qs_env the qs_env to update
!> \par History
!>      07.2019 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_local_stress(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_local_stress', &
         routineP = moduleN//':'//routineN

      CALL cite_reference(Rogers2002)

      CPASSERT(ASSOCIATED(qs_env))

   END SUBROUTINE qs_local_stress

! **************************************************************************************************

END MODULE qs_local_properties
