!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/scp_scf [1.0] *
!!
!!   NAME
!!     scp_scf
!!
!!   FUNCTION
!!     SCP minimization routine
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!   SOURCE
!******************************************************************************

MODULE scp_scf
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_integral_ab,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_zero
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_scp_methods,               ONLY: qs_ks_scp_did_change,&
                                             qs_ks_scp_update
  USE qs_ks_scp_types,                 ONLY: qs_ks_scp_env_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             get_aux_coeff
  USE scp_density_methods,             ONLY: integrate_rhoscp_vrspace,&
                                             update_rhoscp
  USE scp_diis_utils,                  ONLY: diis_info_create,&
                                             diis_info_destroy,&
                                             diis_info_type,&
                                             do_diis_step
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE scp_hartree_1center,             ONLY: integrate_a_vhscp_b,&
                                             integrate_vhscp_gcore,&
                                             integrate_vhscp_gscp
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: scp_qs_scf

!******************************************************************************

CONTAINS

!******************************************************************************
!!****f** scp_scf/scp_qs_scf [1.1] *
!!
!!   NAME
!!    scp_do_scf
!!
!!   FUNCTION
!!   Performs SCF for SCP variables only within qs_environment
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!*** *************************************************************************
SUBROUTINE scp_qs_scf ( qs_env, error )

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE scp_qs_scf (MODULE scp_scf)"

    INTEGER                                  :: handle, handle2, i, ikind, &
                                                iter, max_diis, nkind, &
                                                output_unit, stat
    LOGICAL                                  :: ionode, should_stop
    REAL(dp), POINTER                        :: fc( :, : )
    REAL(KIND=dp)                            :: maxerr, maxgrad, t1, t2
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(diis_info_type), POINTER            :: diis_info
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), POINTER                 :: vtotal
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_environment_type), POINTER      :: scp_env

  CALL timeset ( 'scp_qs_scf','I',' ',handle )
  max_diis = 7
  NULLIFY ( scp_env, local_particles,vtotal )
  NULLIFY ( aux_coeff_set, local_coeffs, para_env )
  CALL get_qs_env(qs_env, &
                  scp_env = scp_env, &
                  local_particles = local_particles, &
                  atomic_kind_set = atomic_kind_set, &
                  pw_env = pw_env, &
                  para_env=para_env,error=error)

  CALL pw_env_get ( pw_env, auxbas_pw_pool = auxbas_pw_pool, error = error )

  CALL get_scp_env ( scp_env, aux_coeff_set = aux_coeff_set, &
                     energy = energy, error = error )

  ALLOCATE ( vtotal, stat = stat )
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
  CALL pw_pool_create_pw(auxbas_pw_pool,vtotal%pw,&
       use_data=REALDATA3D,in_space=REALSPACE,error=error)

  CALL pw_zero(vtotal%pw) 

  logger => cp_error_get_logger(error)
  ionode=para_env%ionode
  output_unit=cp_logger_get_default_unit_nr(logger)
  stat = 0
  iter = 0

  t1 = m_walltime()
  
  CALL diis_info_create ( diis_info, max_diis, aux_coeff_set )

  CALL get_scp_energy ( qs_env, vtotal, error )

  IF (ionode) THEN
     WRITE (UNIT=output_unit,&
        FMT="(/,T3,A,T16,A,T49,A,T68,A,/,T3,A)")&
        "Step","maxerr","Time","Convergence","Total energy",&
         REPEAT("-",77)
  END IF

scf_loop: DO

    CALL timeset('scp_scf_iter','I','',handle2)
    IF (ionode) CALL m_flush(output_unit)
    iter = iter + 1

    CALL do_diis_step ( scp_env, atomic_kind_set, local_particles,  &
                        diis_info, maxerr, error=error )

    CALL get_scp_energy ( qs_env, vtotal, error )
! print
  IF ( qs_env % qmmm )  WRITE (UNIT=output_unit,FMT="((T3,A,T55,F25.14))")&
       "SCP Hartree (SCP density, QMMM potential)       ",energy % e_scp_qmmm
  WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
       "SCP Hartree (SCP density, SCP potential)       ",energy % e_scp_hartree-energy%e_scp_self,&
       "SCP Hartree (KS density, SCP potential)        ",energy % e_scp_ks-energy%e_scp_ks_self-energy%e_scp_core,&
       "Polarization Kernel                            ",energy % e_scp_kernel, &
       "SCP TOTAL                                      ",energy % e_scp_total

! get max_gradient
    nkind = SIZE ( aux_coeff_set % coeffs_of_kind )
    maxgrad=0._dp
    DO ikind= 1, nkind
      local_coeffs => aux_coeff_set%coeffs_of_kind(ikind)%coeffs
      IF(ASSOCIATED(local_coeffs)) THEN
        CALL get_aux_coeff(coeffs=local_coeffs,fc=fc,error=error)
        maxgrad=MAX(maxgrad,MAXVAL(ABS(fc)))
      END IF
    END DO
    CALL mp_max(maxgrad,para_env%group)
    t2 = m_walltime()

    IF (ionode) THEN
       WRITE (UNIT=output_unit,&
          FMT="(T2,I5,2X,E10.2,T32,F8.2,T40,2F20.10)")&
          iter,maxerr,t2 - t1,maxgrad,energy%e_scp_total
    END IF

! ** convergence check
    IF (para_env%source==para_env%mepos) &
            INQUIRE(file="EXIT_SCP",exist=should_stop)
    CALL mp_bcast(should_stop,para_env%source,para_env%group)
    IF (maxgrad < 1.e-8_dp) THEN
      IF (ionode) THEN
         WRITE(UNIT=output_unit,FMT="(/,T3,A,I5,A/)")&
               "*** SCF run converged in", iter, "steps ***"
      END IF
      CALL timestop(0.0_dp,handle2)
      EXIT scf_loop
    ELSE IF (should_stop.OR.&
        iter == 100) THEN
      IF (ionode) THEN
             WRITE(UNIT=output_unit,FMT="(/,T3,A,/)")&
                  "*** SCF run NOT converged ***"
      END IF
      CALL timestop(0.0_dp,handle2)
      EXIT
    END IF

    t1 = m_walltime()

    CALL timestop(0.0_dp,handle2)

  END DO scf_loop

  CALL pw_release ( vtotal%pw,error=error)

  DEALLOCATE ( vtotal, stat = stat )
  CPPostconditionNoFail(stat==0,cp_warning_level,routine,error)
 
  CALL diis_info_destroy(diis_info)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE scp_qs_scf
!******************************************************************************
SUBROUTINE get_scp_energy ( qs_env, v, error )
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type), POINTER                 :: v
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: i, icoef_atom, ikind, ispin, &
                                                n_els, ncoef_atom, nkind
    REAL(dp)                                 :: alpha_pol, ecoeff, energy
    REAL(dp), POINTER                        :: c( :, : ), fc( :, : )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind, &
                                                atomic_kind_set( : )
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(scp_energy_type), POINTER           :: scp_energy
    TYPE(scp_environment_type), POINTER      :: scp_env

  NULLIFY ( local_coeffs, scp_energy, aux_coeff_set, particle_set )
  NULLIFY ( atomic_kind_set, local_particles, scp_env )
  NULLIFY ( rho, cell, matrix_ks, ks_scp_env, ks_env )
  NULLIFY ( dft_control, pw_env, para_env, scp_energy )

  CALL get_qs_env(qs_env,matrix_ks=matrix_ks,rho=rho,&
                  ks_env = ks_env, &
                  scp_env = scp_env, &
                  dft_control = dft_control, &
                  atomic_kind_set = atomic_kind_set, &
                  local_particles = local_particles, &
                  particle_set = particle_set, &
                  pw_env = pw_env, &
                  cell = cell, &
                  para_env=para_env,error=error)

  CALL get_scp_env ( scp_env, energy = scp_energy, &
                     ks_scp_env = ks_scp_env, &
                     aux_coeff_set = aux_coeff_set, &
                     error = error )

  CALL integrate_a_vhscp_b ( qs_env, matrix_ks, rho%rho_ao, .FALSE., error=error ) 

  CALL update_rhoscp ( qs_env, error )
  CALL qs_ks_scp_did_change ( qs_env, full_reset = .TRUE., error = error )
  CALL qs_ks_scp_update ( qs_env, just_energy = .FALSE., error = error )
  CALL integrate_a_vhscp_b ( qs_env, matrix_ks, rho%rho_ao, .FALSE., error=error ) 

  ! 1 center contributions using V_HSCP using the basis-set
  CALL integrate_vhscp_gscp ( scp_env, atomic_kind_set, local_particles, &
                              just_energy = .FALSE., error=error )
  CALL integrate_vhscp_gcore ( scp_env, atomic_kind_set, local_particles, &
                               just_energy = .FALSE., error=error )

  energy = 0.0_dp
  DO ispin=1, dft_control%nspins
    energy = energy + pw_integral_ab ( rho%rho_r(ispin)%pw, ks_scp_env%v_scp_rspace%pw )
  END DO
  scp_energy % e_scp_ks = energy

! Coefficient forces due to SCP density and V_H+HSCP (e.g. d/dc [E_H [\rho_KS + \rho_SCP]}

  v%pw%cr3d(:,:,:) = ks_env%v_hartree_rspace%pw%cr3d(:,:,:) + &
                          ks_scp_env%v_scp_rspace%pw%cr3d(:,:,:)* &
                          ks_env%v_hartree_rspace%pw%pw_grid%dvol 
  CALL integrate_rhoscp_vrspace ( scp_env, v, &
                                  pw_env, atomic_kind_set, &
                                  particle_set, local_particles, cell, &
                                  dft_control%qs_control%eps_gvg_rspace, &
                                  .FALSE., &
                                  .FALSE., error=error )
! Compute the contribution of the kernel (c^2/2a) to energy and coefficient forces
  ecoeff = 0._dp
  nkind = SIZE ( atomic_kind_set )
  DO ikind = 1, nkind 
    atomic_kind => atomic_kind_set ( ikind )
    CALL get_atomic_kind ( atomic_kind = atomic_kind, alpha_pol = alpha_pol )
    local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
    IF ( ASSOCIATED ( local_coeffs ) ) THEN
      CALL get_aux_coeff ( coeffs = local_coeffs, fc = fc, c = c, &
                           n_els = n_els, ncoef_atom = ncoef_atom, &
                           error = error )
      DO i = 1, n_els
        DO icoef_atom = 1, ncoef_atom
          ecoeff = ecoeff +  c ( i, icoef_atom ) * c ( i, icoef_atom ) / 2._dp / alpha_pol
          fc ( i, icoef_atom ) = fc ( i, icoef_atom ) - c ( i, icoef_atom ) / alpha_pol
        END DO
      END DO
    ENDIF
  END DO
  CALL mp_sum ( ecoeff, para_env % group )
  scp_energy % e_scp_kernel = ecoeff

  scp_energy % e_scp_total = scp_energy % e_scp_kernel + scp_energy % e_scp_hartree + &
                             scp_energy % e_scp_ks + scp_energy % e_scp_ks_self + &
                             scp_energy % e_scp_core + scp_energy % e_scp_qmmm + &
                             scp_energy % e_scp_self

END SUBROUTINE get_scp_energy

!******************************************************************************
END MODULE scp_scf
