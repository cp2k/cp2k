!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_rho_methods [1.0] *
!!
!!   NAME
!!     qs_rho_methods
!!
!!   FUNCTION
!!     methods of the rho structure (defined in qs_rho_types)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_rho_methods
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dft_types,                       ONLY: dft_control_type
  USE kinds,                           ONLY: wp=>dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_core,&
                                             calculate_rho_elec,&
                                             calculate_total_rho
  USE qs_core_energies,                ONLY: calculate_ecore,&
                                             calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_create,&
                                             qs_rho_did_change,&
                                             qs_rho_get,&
                                             qs_rho_release,&
                                             qs_rho_retain,&
                                             qs_rho_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type,&
                                             real_matrix_type,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_rho_methods'

  PUBLIC :: qs_rho_update_rho

!***
!****************************************************************************
CONTAINS

!!****f* qs_rho_methods/qs_rho_update_rho [1.0] *
!!
!!   NAME
!!     qs_rho_update_rho
!!
!!   FUNCTION
!!     updates rho_r and rho_g to the new (given) ao.
!!
!!   NOTES
!!     make the default qs_chages qs_env%qs_charges?
!!
!!   INPUTS
!!     - rho_struct: the rho structure that should be updated
!!     - rho_ao_new: the new value for rho_ao
!!     - qs_env: the qs_env rho_struct refers to
!!     - qs_charges: the object where to store
!!       the integrated charge in r space
!!     - qs_energy: the object where to store the energies
!!     - rebuild: if false builds rho_r by just updating the parts
!!       of rho that changed with respect to the old density,
!!       (by default looks at how may updates were done without rebuild).
!!     - update_rho_ao: if rho_ao should also be updated (defaults
!!       to true)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_rho_update_rho(rho_struct, rho_ao_new,qs_env,&
     qs_charges,energy,rebuild,update_rho_ao, error)
  TYPE(qs_rho_type), POINTER :: rho_struct
  TYPE(qs_environment_type), INTENT(in) :: qs_env
  TYPE(real_matrix_p_type),DIMENSION(:), INTENT(in) :: rho_ao_new
  LOGICAL, INTENT(in),OPTIONAL :: rebuild, update_rho_ao
  TYPE(qs_charges_type), POINTER, OPTIONAL :: qs_charges
  TYPE(qs_energy_type), POINTER, OPTIONAL :: energy
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, rebuild_local, update_rho_ao_local
  INTEGER :: nspins, ispin, stat
  CHARACTER(len=*), PARAMETER :: routineN='qs_rho_update_rho',&
        routineP=moduleN//':'//routineN
  REAL(kind=wp) :: ecore, charge_tot
  TYPE(dft_control_type), POINTER :: dft_control
  TYPE(coeff_type) :: tmp_gspace
  TYPE(coeff_type), DIMENSION(:), POINTER :: mgrid_rspace, mgrid_gspace
  TYPE(pw_env_type), POINTER :: pw_env
  TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
  TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h, s
  failure=.FALSE.
  update_rho_ao_local=.TRUE.
  IF (PRESENT(update_rho_ao)) update_rho_ao_local=update_rho_ao
  NULLIFY(dft_control, mgrid_rspace, mgrid_gspace, pw_env, auxbas_pw_pool,&
       pw_pools, para_env, h, s)

  CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
  CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL get_qs_env(qs_env, dft_control=dft_control, h=h,&
          s=s,pw_env=pw_env, para_env=para_env)
     nspins=dft_control%nspins
     ! move update of ecore to some other place ?
     energy%core=0.0_wp
     DO ispin=1,nspins
        CALL calculate_ecore(h=h(1)%matrix,&
             p=rho_ao_new(ispin)%matrix,&
             s=s(1)%matrix,&
             ecore=ecore,&
             para_env=para_env,&
             error=error)
        energy%core=energy%core+ecore
     END DO
     
     ! decides if to rebuild or just update
     rebuild_local=rho_struct%updates_since_rebuild > rho_struct%rebuild_each
     IF (PRESENT(rebuild)) rebuild_local=rebuild
     IF (.NOT.(rho_struct%rho_ao_valid.AND.rho_struct%rho_elec_r_valid)) THEN
        rebuild_local=.TRUE.
     END IF
     IF (rebuild_local) THEN
        rho_struct%updates_since_rebuild=0
     ELSE
        rho_struct%updates_since_rebuild=rho_struct%updates_since_rebuild+1
     END IF
     
     ! alloc multigrids
     CALL pw_env_get(pw_env, pw_pools=pw_pools, auxbas_pw_pool=auxbas_pw_pool,&
          error=error)
     ALLOCATE(mgrid_gspace(SIZE(pw_pools)), mgrid_rspace(SIZE(pw_pools)),&
          stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     ! get tmp grids
     CALL pw_pool_init_coeff(auxbas_pw_pool, tmp_gspace, &
          use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)
     CALL pw_pools_init_coeffs(pw_pools, mgrid_rspace,&
          use_data=REALDATA3D, in_space=REALSPACE, error=error)
     CALL pw_pools_init_coeffs(pw_pools, mgrid_gspace,&
          use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE, error=error)
     
     DO ispin=1,nspins
        CALL calculate_rho_elec(p=rho_ao_new(ispin)%matrix,&
             p_old=rho_struct%rho_ao(ispin)%matrix,&
             rebuild=rebuild_local,&
             rho=rho_struct%rho_elec_r(ispin),&
             rho_gspace=tmp_gspace,&
             rho_gspace_old=rho_struct%rho_g(ispin),&
             mgrid_rspace=mgrid_rspace,&
             mgrid_gspace=mgrid_gspace,&
             gridlevel_info=pw_env%gridlevel_info,&
             cube_info=pw_env%cube_info,&
             l_info=pw_env%l_info,&
             total_rho=qs_charges%total_rho_elec_rspace(ispin),&
             qs_env=qs_env,error=error)
     END DO

     ! put back grids
     CALL pw_pool_give_back_coeff(auxbas_pw_pool, tmp_gspace, error=error)
     CALL pw_pools_give_back_coeffs(pw_pools, mgrid_rspace, error=error)
     CALL pw_pools_give_back_coeffs(pw_pools, mgrid_gspace, error=error)
     DEALLOCATE(mgrid_gspace, mgrid_rspace,stat=stat)
     
     rho_struct%rho_ao_valid=.FALSE.
     rho_struct%rho_elec_r_valid=.TRUE.
     rho_struct%rho_g_valid=.TRUE.
     IF (update_rho_ao_local) THEN
        DO ispin=1,nspins
           CALL transfer_matrix(rho_ao_new(ispin)%matrix,&
                rho_struct%rho_ao(ispin)%matrix)
        END DO
        rho_struct%rho_ao_valid=.TRUE.
     END IF
  END IF
END SUBROUTINE qs_rho_update_rho
!***************************************************************************

END MODULE qs_rho_methods
