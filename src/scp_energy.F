!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/scp_energy [1.0] *
!!
!!   NAME
!!     scp_energy
!!
!!   FUNCTION
!!     scp_energy minim routine
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!   SOURCE
!******************************************************************************

MODULE scp_energy
!  USE input_cp2k_restarts,             ONLY: write_restart
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             get_aux_coeff
  USE scp_density_methods,             ONLY: integrate_rhoscp_vrspace
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE scp_force_types,                 ONLY: scp_force_type
  USE scp_hartree_1center,             ONLY: integrate_vhscp_gcore,&
                                             integrate_vhscp_gscp
  USE scp_rho_types,                   ONLY: scp_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: scp_qs_energies

!******************************************************************************

CONTAINS

!!****f** scp_energy/scp_qs_energies [1.1] *
!!
!!   NAME
!!    scp_qs_energies
!!
!!   FUNCTION
!!     Calculates the pieces of the SCP energy involved
!!     in the optimization of the polarization coefficients [AUX_COEF_SET]
!!     within a GPW calculation
!!
!!   AUTHOR
!!    CJM
!!
!!   MODIFICATION HISTORY
!!
!!*** *************************************************************************
SUBROUTINE scp_qs_energies ( qs_env, calculate_forces, just_energy, error )

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces, just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE scp_qs_energies (MODULE scp_energy)"

    INTEGER                                  :: handle, i, icoef_atom, ikind, &
                                                n_els, ncoef_atom, nkind, &
                                                output_unit, stat
    LOGICAL                                  :: energy_only, failure, ionode, &
                                                my_calculate_forces
    REAL(dp)                                 :: alpha_pol, ecoeff
    REAL(dp), DIMENSION(:, :), POINTER       :: coeff, force
    REAL(KIND=dp)                            :: t1
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(scp_energy_type), POINTER           :: scp_energy
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(scp_force_type), POINTER            :: scp_force( : )
    TYPE(scp_rho_type), POINTER              :: rho_scp

!    TYPE(section_vals_type), POINTER         :: root_section

  CALL timeset ( 'scp_qs_energies','I',' ',handle )

  failure = .FALSE.
  NULLIFY(atomic_kind_set, aux_coeff_set, scp_env, rho_scp, scp_energy)
  NULLIFY(para_env,logger,force, coeff, local_coeffs, local_particles, scp_force)
  para_env=>qs_env%para_env
  logger => cp_error_get_logger(error)
  ionode=para_env%ionode
  output_unit=cp_logger_get_default_unit_nr(logger)
  stat = 0
  t1 = m_walltime()
  my_calculate_forces = .FALSE.
  IF ( PRESENT ( calculate_forces ) ) my_calculate_forces = calculate_forces
  energy_only = .FALSE.
  IF ( PRESENT ( just_energy ) ) energy_only = just_energy

  CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, &
                    atomic_kind_set=atomic_kind_set,    &
                    local_particles=local_particles, &
                    particle_set = particle_set,  &
                    energy=qs_energy, pw_env = pw_env, &
                    dft_control = dft_control, cell = cell, &
                    error = error )

! Getting the coeffs and the density structures
  CALL get_scp_env ( scp_env=scp_env,&
                     aux_coeff_set=aux_coeff_set,&
                     rho_scp=rho_scp, &
                     energy=scp_energy, error=error)

  IF ( my_calculate_forces ) THEN 
    CALL get_scp_env ( scp_env=scp_env,&
                       force=scp_force, error=error)
! Zeroing forces if needed
    DO ikind = 1, SIZE ( atomic_kind_set )
      scp_force ( ikind ) % f_scp ( :, : ) = 0.0
      scp_force ( ikind ) % f_total ( :, : ) = 0.0
    END DO
  END IF

! Coeff forces have been zeroed by a call to aux_coeff_zero_fc in integrate_a_vhscp_b

! 1 center contributions using V_HSCP using the basis-set
  CALL integrate_vhscp_gscp ( scp_env, atomic_kind_set, local_particles, &
                              just_energy = energy_only, error=error )
  CALL integrate_vhscp_gcore ( scp_env, atomic_kind_set, local_particles, &
                               just_energy = energy_only, error=error )
! Coefficient forces due to SCP density and V_H+HSCP (e.g. d/dc [E_H [\rho_KS + \rho_SCP]}
  CALL integrate_rhoscp_vrspace ( scp_env, qs_env%ks_env%v_hartree_rspace, &
                                  pw_env, atomic_kind_set, &
                                  particle_set, local_particles, cell, &
                                  dft_control%qs_control%eps_gvg_rspace, &
                                  calculate_forces = my_calculate_forces, &
                                  just_energy = energy_only, error=error )
! Compute the contribution of the kernel (c^2/2a) to energy and coefficient forces
  ecoeff = 0._dp
  nkind = SIZE ( atomic_kind_set )
  DO ikind = 1, nkind 
    atomic_kind => atomic_kind_set ( ikind )
    CALL get_atomic_kind ( atomic_kind = atomic_kind, alpha_pol = alpha_pol )
    local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
    IF ( ASSOCIATED ( local_coeffs ) ) THEN
      CALL get_aux_coeff ( coeffs = local_coeffs, fc = force, c = coeff, &
                           n_els = n_els, ncoef_atom = ncoef_atom, &
                           error = error )
      DO i = 1, n_els
        DO icoef_atom = 1, ncoef_atom
          ecoeff = ecoeff +  coeff ( i, icoef_atom ) * coeff ( i, icoef_atom ) / 2._dp / alpha_pol
          IF ( .NOT. energy_only ) THEN
            force ( i, icoef_atom ) = force ( i, icoef_atom ) - coeff ( i, icoef_atom ) / alpha_pol
          END IF 
        END DO
      END DO
! DBG
!  force ( 1, 1 ) = 0._dp
!  force ( 1, 2 ) = 0._dp
!  force ( 1, 3 ) = 0._dp
! DBG
    ENDIF
  END DO
  CALL mp_sum ( ecoeff, para_env % group )
  scp_energy % e_scp_kernel = ecoeff
  t1 = m_walltime()

  scp_energy%e_scp_total = scp_energy%e_scp_self + scp_energy%e_scp_core &
                  + scp_energy%e_scp_hartree + scp_energy%e_scp_kernel   &
                  + scp_energy % e_scp_ks + scp_energy % e_scp_ks_self   &
                  + scp_energy % e_scp_qmmm
! DBG
!  WRITE ( 55, * ) coeff ( 1, 1 ), scp_energy % e_scp_total
! DBG
! print
  IF ( qs_env % qmmm )  WRITE (UNIT=output_unit,FMT="((T3,A,T55,F25.14))")&
       "SCP Hartree (SCP density, QMMM potential)       ",scp_energy % e_scp_qmmm
  WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
       "SCP Hartree (SCP,SCP)                          ",scp_energy % e_scp_hartree,&
       "SCP Hartree Self (SCP,SCP)                     ",scp_energy%e_scp_self,&
       "SCP Hartree (KS,SCP)                           ",scp_energy % e_scp_ks,&
       "SCP Hartree Self (KS,SCP)                      ",scp_energy%e_scp_ks_self,&
       "SCP Hartree Self Core (KS core,SCP)            ",scp_energy%e_scp_core,&
       "Polarization Kernel                            ",scp_energy % e_scp_kernel, &
       "KS TOTAL                                       ",qs_energy % total, &
       "SCP TOTAL                                      ",scp_energy % e_scp_total

! Add the scp energies to qs_energy%total for consistent reporting in the SCF
!  qs_energy%total = qs_energy%total + scp_energy%e_scp_total
         
  CALL m_flush(output_unit)

!  CALL write_restart ( scp_env=scp_env, root_section=root_section,
!                       force_env_section=scp_env%input, error=error)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE scp_qs_energies

!******************************************************************************
END MODULE scp_energy
