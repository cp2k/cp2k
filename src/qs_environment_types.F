!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****** cp2k/qs_environment_types [1.0] *
!!
!!   NAME
!!     qs_environment_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MK (23.01.2002)
!!
!!   MODIFICATION HISTORY
!!     - mo_set_p_type added to qs_env (23.04.02,MK)
!!     - qs_force_type added to qs_env (05.06.02,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_environment_types
  !USE cp_matrix_utils,        ONLY: cp_sm_set_dist_sparsity_from_sm
  !USE pw_grid_types,          ONLY: HALFSPACE,&
  !                                  pw_grid_type
  !USE pw_grids,               ONLY: pw_grid_construct,&
  !                                  pw_grid_setup
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type,&
                                             cp_bmstruct_create,&
                                             cp_bmstruct_release
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
  USE cp_f_matrix_struct,              ONLY: cp_f_matrix_struct_type,&
                                             cp_fmstruct_create,&
                                             cp_fmstruct_release,&
                                             cp_fmstruct_retain
  USE cp_fm_pool,                      ONLY: cp_fm_pool_type,&
                                             fm_pool_create,&
                                             fm_pool_release,&
                                             fm_pool_retain
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dft_types,                       ONLY: dft_control_type
  USE kinds,                           ONLY: wp => dp
  USE pao_types,                       ONLY: pao_env_type
  USE particle_types,                  ONLY: particle_type
  USE pw_pool,                         ONLY: pw_pool_create,&
                                             pw_pool_p_type,&
                                             pw_pool_release,&
                                             pw_pool_retain,&
                                             pw_pool_type,&
                                             pw_pools_copy,&
                                             pw_pools_dealloc
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_mo_types,                     ONLY: mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_retain,&
                                             qs_rho_type
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_environment_types'

  TYPE qs_environment_type
!   PRIVATE
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(cell_type), POINTER                              :: cell
    TYPE(dft_control_type), POINTER                       :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER            :: c
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_erfc,&
                                                             sac_ppl,sca_ppnl,&
                                                             scb_ppnl
    TYPE(pao_env_type), POINTER                           :: pao_env
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(qs_energy_type), POINTER                         :: energy
    TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: h,k,p,p_old,s,kinetic
    TYPE(cp_b_matrix_struct_type), POINTER                :: s_mstruct
    TYPE(qs_rho_type), POINTER                            :: rho
    TYPE(scf_control_type), POINTER                       :: scf_control

    TYPE(cp_fm_pool_type), POINTER                        :: ao_mo_fm_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER           :: pw_pools

    TYPE(cp_para_env_type), POINTER                       :: para_env
  END TYPE qs_environment_type

! *** Public data types ***

  PUBLIC :: qs_environment_type

! *** Public subroutines ***

  PUBLIC :: get_qs_env,&
            init_qs_env,&
            set_qs_env

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,p_old,s, s_mstruct,&
                        rho, ao_mo_fm_pool, pw_pools, para_env, scf_control,&
                        kinetic,error)

!   Purpose: Get the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(IN)       :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER          :: cell
    TYPE(dft_control_type), OPTIONAL, POINTER   :: dft_control
    TYPE(mo_set_p_type),  DIMENSION(:),&
      OPTIONAL, POINTER                         :: c
    TYPE(neighbor_list_set_p_type),&
      DIMENSION(:), OPTIONAL, POINTER           :: sab_orb,sac_erfc,sac_ppl,&
                                                   sca_ppnl,scb_ppnl
    TYPE(pao_env_type), OPTIONAL, POINTER       :: pao_env
    TYPE(particle_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER     :: energy
    TYPE(qs_force_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: force
    TYPE(real_matrix_p_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: h,k,p,p_old,s,kinetic
    TYPE(cp_b_matrix_struct_type), OPTIONAL,&
         POINTER                                :: s_mstruct
    TYPE(qs_rho_type), OPTIONAL, POINTER        :: rho
    TYPE(cp_fm_pool_type), OPTIONAL, POINTER    :: ao_mo_fm_pool
    TYPE(pw_pool_p_type), OPTIONAL, &
         DIMENSION(:), POINTER                  :: pw_pools
    TYPE(cp_para_env_type), OPTIONAL, POINTER   :: para_env
    TYPE(scf_control_type), OPTIONAL, POINTER   :: scf_control
    TYPE(cp_error_type), OPTIONAL,INTENT(inout) :: error
!   ---------------------------------------------------------------------------

    IF (PRESENT(atomic_kind_set)) atomic_kind_set => qs_env%atomic_kind_set
    IF (PRESENT(cell)) cell => qs_env%cell
    IF (PRESENT(dft_control)) dft_control => qs_env%dft_control
    IF (PRESENT(c)) c => qs_env%c
    IF (PRESENT(sab_orb)) sab_orb => qs_env%sab_orb
    IF (PRESENT(sac_erfc)) sac_erfc => qs_env%sac_erfc
    IF (PRESENT(sac_ppl)) sac_ppl => qs_env%sac_ppl
    IF (PRESENT(sca_ppnl)) sca_ppnl => qs_env%sca_ppnl
    IF (PRESENT(scb_ppnl)) scb_ppnl => qs_env%scb_ppnl
    IF (PRESENT(pao_env)) pao_env => qs_env%pao_env
    IF (PRESENT(particle_set)) particle_set => qs_env%particle_set
    IF (PRESENT(energy)) energy => qs_env%energy
    IF (PRESENT(force)) force => qs_env%force
    IF (PRESENT(h)) h => qs_env%h
    IF (PRESENT(k)) k => qs_env%k
    IF (PRESENT(p)) p => qs_env%p
    IF (PRESENT(kinetic)) kinetic => qs_env%kinetic
    IF (PRESENT(p_old)) p_old => qs_env%p_old
    IF (PRESENT(s)) s => qs_env%s
    IF (PRESENT(s_mstruct)) s_mstruct => qs_env%s_mstruct
    IF (PRESENT(rho)) rho => qs_env%rho
    IF (PRESENT(pw_pools)) pw_pools => qs_env%pw_pools
    IF (PRESENT(ao_mo_fm_pool)) ao_mo_fm_pool => qs_env%ao_mo_fm_pool
    IF (PRESENT(para_env)) para_env => qs_env%para_env
    IF (PRESENT(scf_control)) scf_control => qs_env%scf_control

  END SUBROUTINE get_qs_env

! *****************************************************************************

  SUBROUTINE init_qs_env(qs_env,para_env,error)

!   Purpose: Initialise the QUICKSTEP environment.

!   History: - Creation (25.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)       :: qs_env
    TYPE(cp_para_env_type), POINTER              :: para_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   ---------------------------------------------------------------------------

    NULLIFY (qs_env%atomic_kind_set)
    NULLIFY (qs_env%cell)
    NULLIFY (qs_env%dft_control)
    NULLIFY (qs_env%c)
    NULLIFY (qs_env%sab_orb)
    NULLIFY (qs_env%sac_erfc)
    NULLIFY (qs_env%sac_ppl)
    NULLIFY (qs_env%sca_ppnl)
    NULLIFY (qs_env%scb_ppnl)
    NULLIFY (qs_env%pao_env)
    NULLIFY (qs_env%particle_set)
    NULLIFY (qs_env%energy)
    NULLIFY (qs_env%force)
    NULLIFY (qs_env%h)
    NULLIFY (qs_env%k)
    NULLIFY (qs_env%p)
    NULLIFY (qs_env%kinetic)
    NULLIFY (qs_env%p_old)
    NULLIFY (qs_env%s)
    NULLIFY (qs_env%s_mstruct)
    NULLIFY (qs_env%rho)
    NULLIFY (qs_env%ao_mo_fm_pool)
    NULLIFY (qs_env%pw_pools)
    NULLIFY (qs_env%scf_control)

    CALL cp_para_env_retain(para_env,error=error)
    qs_env%para_env => para_env

  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,p_old,s,&
                        pw_pools,ao_mo_fm_pool,rho,kinetic,&
                        scf_control,error)

!   Purpose: Set the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)    :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:),&
      OPTIONAL, POINTER                       :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER        :: cell
    TYPE(dft_control_type), OPTIONAL, POINTER :: dft_control
    TYPE(mo_set_p_type),  DIMENSION(:),&
      OPTIONAL, POINTER                       :: c
    TYPE(neighbor_list_set_p_type),&
      DIMENSION(:), OPTIONAL, POINTER         :: sab_orb,sac_erfc,sac_ppl,&
                                                 sca_ppnl,scb_ppnl
    TYPE(pao_env_type), OPTIONAL, POINTER     :: pao_env
    TYPE(particle_type), DIMENSION(:),&
      OPTIONAL, POINTER                       :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER   :: energy
    TYPE(qs_force_type), DIMENSION(:),&
      OPTIONAL, POINTER                       :: force
    TYPE(real_matrix_p_type), DIMENSION(:),&
         OPTIONAL, POINTER                    :: h,k,p,p_old,s,kinetic
    TYPE(qs_rho_type), POINTER,OPTIONAL       :: rho
    TYPE(cp_fm_pool_type), OPTIONAL, POINTER  :: ao_mo_fm_pool
    TYPE(pw_pool_p_type), OPTIONAL, &
         DIMENSION(:), POINTER                :: pw_pools
    TYPE(scf_control_type), OPTIONAL, POINTER :: scf_control
    TYPE(cp_error_type),INTENT(inout),OPTIONAL:: error

!   ---------------------------------------------------------------------------
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: new_pools
    NULLIFY(new_pools)

    IF (PRESENT(atomic_kind_set)) qs_env%atomic_kind_set => atomic_kind_set
    IF (PRESENT(cell)) qs_env%cell => cell
    IF (PRESENT(dft_control)) qs_env%dft_control => dft_control
    IF (PRESENT(c)) qs_env%c => c
    IF (PRESENT(sab_orb)) qs_env%sab_orb => sab_orb
    IF (PRESENT(sac_erfc)) qs_env%sac_erfc => sac_erfc
    IF (PRESENT(sac_ppl)) qs_env%sac_ppl => sac_ppl
    IF (PRESENT(sca_ppnl)) qs_env%sca_ppnl => sca_ppnl
    IF (PRESENT(scb_ppnl)) qs_env%scb_ppnl => scb_ppnl
    IF (PRESENT(pao_env)) qs_env%pao_env => pao_env
    IF (PRESENT(particle_set)) qs_env%particle_set => particle_set
    IF (PRESENT(energy)) qs_env%energy => energy
    IF (PRESENT(force)) qs_env%force => force
    IF (PRESENT(h)) qs_env%h => h
    IF (PRESENT(k)) qs_env%k => k
    IF (PRESENT(p)) qs_env%p => p
    IF (PRESENT(p_old)) qs_env%p_old => p_old
    IF (PRESENT(s)) qs_env%s => s
    IF (PRESENT(rho)) THEN ! accept also null pointers?
      CALL qs_rho_retain(rho,error=error)
      CALL qs_rho_release(qs_env%rho,error=error)
      qs_env%rho => rho
    END IF
    IF (PRESENT(pw_pools)) THEN
      CALL pw_pools_copy(pw_pools,new_pools,error=error)
      CALL pw_pools_dealloc(qs_env%pw_pools,error=error)
      qs_env%pw_pools => new_pools
    END IF
    IF (PRESENT(ao_mo_fm_pool)) THEN ! accept also null pointers?
      CALL fm_pool_retain(ao_mo_fm_pool,error=error)
      CALL fm_pool_release(qs_env%ao_mo_fm_pool,error=error)
      qs_env%ao_mo_fm_pool => ao_mo_fm_pool
    END IF
    IF (PRESENT(scf_control)) THEN
       CALL scf_c_retain(scf_control,error=error)
       CALL scf_c_release(qs_env%scf_control,error=error)
       qs_env%scf_control => scf_control
    END IF
    IF (PRESENT(kinetic)) qs_env%kinetic => kinetic

  END SUBROUTINE set_qs_env

! *****************************************************************************

END MODULE qs_environment_types
