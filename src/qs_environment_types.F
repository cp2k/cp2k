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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             deallocate_atomic_kind_set,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type,&
                                             cp_bmstruct_create,&
                                             cp_bmstruct_release
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error,&
                                             cp_unimplemented_error_nr
  USE cp_fm_basic_linalg,              ONLY: cp_sm_fm_multiply
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create_fm,&
                                             fm_pool_give_back_fm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_full_matrix_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_matrix_utils,                 ONLY: sm_struct_add_sm_local_blocks
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             cp_sm_pool_type,&
                                             sm_pool_create,&
                                             sm_pool_get_mstruct,&
                                             sm_pool_release,&
                                             sm_pool_retain,&
                                             sm_pools_copy,&
                                             sm_pools_dealloc,&
                                             sm_pools_flush_cache
  USE cp_sm_struct,                    ONLY: cp_sm_struct_type,&
                                             sm_struct_clear_local_blocks,&
                                             sm_struct_create,&
                                             sm_struct_release,&
                                             sm_struct_retain
  USE dft_types,                       ONLY: dft_control_type
  USE kinds,                           ONLY: wp => dp
  USE pao_types,                       ONLY: pao_env_type
  USE particle_types,                  ONLY: deallocate_particle_set,&
                                             particle_type
  USE pw_env_types,                    ONLY: pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE qs_blacs,                        ONLY: blacs_make_basis
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_retain,&
                                             qs_charges_type
  USE qs_energy_types,                 ONLY: deallocate_qs_energy,&
                                             qs_energy_type
  USE qs_force_types,                  ONLY: deallocate_qs_force,&
                                             qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             mpools_release,&
                                             mpools_retain,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_retain,&
                                             qs_rho_type
  USE qs_scf_types,                    ONLY: ot_method_nr,&
                                             qs_scf_env_type,&
                                             scf_env_release,&
                                             scf_env_retain
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_release,&
                                             wfi_retain
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             real_matrix_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_environment_types'
  INTEGER, PRIVATE, SAVE :: last_qs_env_id_nr=0

! *** Public data types ***

  PUBLIC :: qs_environment_type, qs_environment_p_type

! *** Public subroutines ***

  PUBLIC :: get_qs_env,&
            init_qs_env,&
            set_qs_env,&
            qs_env_create, qs_env_retain, qs_env_release,&
            qs_env_reorthogonalize_vectors

!!***
!****************************************************************************

!!****t* qs_environment_types/qs_environment_type [1.0] *
!!
!!   NAME
!!     qs_environment_type
!!
!!   FUNCTION
!!     contains all the info nedeed by quickstep to do dft scf
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number, guarateed to be unique
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - atomic_kind_set: array with infos about the species (atomic_kinds)
!!       present in the system
!!     - cell: the simulation cell (box)
!!     - dft_control: object that contains the values of various parameters
!!       relevant to a dft calculation
!!     - c: the molecular orbitals, in LSD array with the different spins
!!     - sab_orb,sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl,sca_ppnl,scb_ppnl:
!!       neighbor lists for the calculation of the core Hamiltonian matrix
!!     - pao_env: things nedeed by the pao method
!!     - particle_set: info on the atoms you simulate, pos,...
!!     - energy: place where the energies calculated are stored
!!     - h: core hamiltonian (h(1)), if needed by the calculation also
!!       its derivatives wrt. x,y, and z (h(2:4))
!!     - k: Khon-Sham matrix (filtred by the structure of S).
!!       with LSD an array with the differents orbitals
!!     - s: the overlap matrix (s(1)), if needed by the calculation also
!!       its derivatives wrt. x,y, and z (s(2:4))
!!     - kinetic: kinetic part of h
!!     - rho: the (old) density, in various representations (ao+grid)
!!     - scf_control: various parameters relevant to the scf cycle
!!     - pw_env: various things to create and get info on the grids used
!!       by qs
!!     - para_env: parallel environement that should be  used to do the qs
!!       related calculations
!!     - blacs_env: parallel environement needed by scalapack (full matrixes)
!!     - rho_core: a grid with the charges of the cores of the atoms in the
!!       reciprocal (g) space
!!     - qs_charges: an object where to put the charges that are present on 
!!       the various grids
!!     - ks_env: environement for the calculation of the Kohn Sham matrix
!!     - wf_history: place where the snapshots of the previous states are
!!       stored (used during the md to extrapolate the new starting point)
!!     - scf_env: place where temporary information used in the scf is stored
!!     - mpools: all the pools of matrixes
!!     
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     11.2002 added doc and attribute description [fawzi]
!!
!!*** *********************************************************************
  TYPE qs_environment_type
!   PRIVATE
    INTEGER :: id_nr, ref_count
    TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
    TYPE(cell_type), POINTER                              :: cell
    TYPE(dft_control_type), POINTER                       :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER            :: c
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_erfc,&
                                                             sac_ppl,sac_ppnl,&
                                                             sbc_ppnl,sca_ppnl,&
                                                             scb_ppnl
    TYPE(pao_env_type), POINTER                           :: pao_env
    TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
    TYPE(qs_energy_type), POINTER                         :: energy
    TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: h,k,s,kinetic
    TYPE(qs_rho_type), POINTER                            :: rho
    TYPE(scf_control_type), POINTER                       :: scf_control
    TYPE(pw_env_type), POINTER                            :: pw_env
    TYPE(cp_para_env_type), POINTER                       :: para_env
    TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
    TYPE(coeff_type), POINTER                             :: rho_core
    TYPE(qs_charges_type), POINTER                        :: qs_charges
    TYPE(qs_ks_env_type), POINTER                         :: ks_env
    TYPE(qs_wf_history_type), POINTER                     :: wf_history
    TYPE(qs_scf_env_type), POINTER                        :: scf_env
    TYPE(qs_matrix_pools_type), POINTER                   :: mpools
  END TYPE qs_environment_type
! *****************************************************************************

!!****s* qs_environment_types/qs_environment_p_type [1.0] *
!!
!!   NAME
!!     qs_environment_p_type
!!
!!   FUNCTION
!!     to build arrays of pointers
!!
!!   ATTRIBUTES
!!     - qs_env: the pointer to the qs_env
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE qs_environment_p_type
     TYPE(qs_environment_type), POINTER :: qs_env
  END TYPE qs_environment_p_type
!!***
  !****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE get_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl,sca_ppnl,scb_ppnl,&
                        pao_env,particle_set,energy,force,h,k,s,&
                        rho, pw_env, mpools,&
                        para_env, blacs_env, scf_control,&
                        kinetic,qs_charges, rho_core, ks_env,&
                        wf_history,scf_env,&
                        id_nr,error)

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
                                                   sac_ppnl,sbc_ppnl,sca_ppnl,&
                                                   scb_ppnl
    TYPE(pao_env_type), OPTIONAL, POINTER       :: pao_env
    TYPE(particle_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER     :: energy
    TYPE(qs_force_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: force
    TYPE(real_matrix_p_type), DIMENSION(:),&
      OPTIONAL, POINTER                         :: h,k,s,kinetic
    TYPE(qs_matrix_pools_type), OPTIONAL, POINTER :: mpools
    TYPE(qs_rho_type), OPTIONAL, POINTER        :: rho
    TYPE(pw_env_type), OPTIONAL, POINTER        :: pw_env
    TYPE(cp_para_env_type), OPTIONAL, POINTER   :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, POINTER  :: blacs_env
    TYPE(scf_control_type), OPTIONAL, POINTER   :: scf_control
    TYPE(qs_charges_type), OPTIONAL, POINTER    :: qs_charges
    TYPE(coeff_type), OPTIONAL, POINTER         :: rho_core
    TYPE(qs_ks_env_type), OPTIONAL, POINTER     :: ks_env
    TYPE(qs_wf_history_type), OPTIONAL, POINTER :: wf_history
    TYPE(qs_scf_env_type), OPTIONAL, POINTER    :: scf_env
    INTEGER, INTENT(out), OPTIONAL              :: id_nr
    TYPE(cp_error_type), OPTIONAL,INTENT(inout) :: error
!   ---------------------------------------------------------------------------

    IF (PRESENT(atomic_kind_set)) atomic_kind_set => qs_env%atomic_kind_set
    IF (PRESENT(cell)) cell => qs_env%cell
    IF (PRESENT(dft_control)) dft_control => qs_env%dft_control
    IF (PRESENT(c)) c => qs_env%c
    IF (PRESENT(sab_orb)) sab_orb => qs_env%sab_orb
    IF (PRESENT(sac_erfc)) sac_erfc => qs_env%sac_erfc
    IF (PRESENT(sac_ppl)) sac_ppl => qs_env%sac_ppl
    IF (PRESENT(sac_ppnl)) sac_ppnl => qs_env%sac_ppnl
    IF (PRESENT(sbc_ppnl)) sbc_ppnl => qs_env%sbc_ppnl
    IF (PRESENT(sca_ppnl)) sca_ppnl => qs_env%sca_ppnl
    IF (PRESENT(scb_ppnl)) scb_ppnl => qs_env%scb_ppnl
    IF (PRESENT(pao_env)) pao_env => qs_env%pao_env
    IF (PRESENT(particle_set)) particle_set => qs_env%particle_set
    IF (PRESENT(energy)) energy => qs_env%energy
    IF (PRESENT(force)) force => qs_env%force
    IF (PRESENT(h)) h => qs_env%h
    IF (PRESENT(k)) k => qs_env%k
    IF (PRESENT(kinetic)) kinetic => qs_env%kinetic
    IF (PRESENT(s)) s => qs_env%s
    IF (PRESENT(rho)) rho => qs_env%rho
    IF (PRESENT(pw_env)) pw_env => qs_env%pw_env
    IF (PRESENT(mpools)) mpools => qs_env%mpools
    IF (PRESENT(para_env)) para_env => qs_env%para_env
    IF (PRESENT(blacs_env)) blacs_env => qs_env%blacs_env
    IF (PRESENT(scf_control)) scf_control => qs_env%scf_control
    IF (PRESENT(rho_core)) rho_core => qs_env%rho_core
    IF (PRESENT(qs_charges)) qs_charges => qs_env%qs_charges
    IF (PRESENT(ks_env)) ks_env => qs_env%ks_env
    IF (PRESENT(wf_history)) wf_history => qs_env%wf_history
    IF (PRESENT(scf_env)) scf_env => qs_env%scf_env
    IF (PRESENT(id_nr)) id_nr=qs_env%id_nr

  END SUBROUTINE get_qs_env

! *****************************************************************************

  SUBROUTINE init_qs_env(qs_env,para_env, blacs_env, error)

!   Purpose: Initialise the QUICKSTEP environment.

!   History: - Creation (25.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)       :: qs_env
    TYPE(cp_para_env_type), POINTER              :: para_env
    TYPE(cp_blacs_env_type), POINTER             :: blacs_env
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   ---------------------------------------------------------------------------

    NULLIFY (qs_env%atomic_kind_set)
    NULLIFY (qs_env%cell)
    NULLIFY (qs_env%dft_control)
    NULLIFY (qs_env%c)
    NULLIFY (qs_env%sab_orb)
    NULLIFY (qs_env%sac_erfc)
    NULLIFY (qs_env%sac_ppl)
    NULLIFY (qs_env%sac_ppnl)
    NULLIFY (qs_env%sbc_ppnl)
    NULLIFY (qs_env%sca_ppnl)
    NULLIFY (qs_env%scb_ppnl)
    NULLIFY (qs_env%pao_env)
    NULLIFY (qs_env%particle_set)
    NULLIFY (qs_env%energy)
    NULLIFY (qs_env%force)
    NULLIFY (qs_env%h)
    NULLIFY (qs_env%k)
    NULLIFY (qs_env%kinetic)
    NULLIFY (qs_env%s)
    NULLIFY (qs_env%mpools)
    NULLIFY (qs_env%rho)
    NULLIFY (qs_env%pw_env)
    NULLIFY (qs_env%scf_control)
    NULLIFY (qs_env%qs_charges)
    NULLIFY (qs_env%rho_core)
    NULLIFY (qs_env%ks_env)
    NULLIFY (qs_env%wf_history)
    NULLIFY (qs_env%scf_env)

    CALL cp_para_env_retain(para_env,error=error)
    qs_env%para_env => para_env
    CALL cp_blacs_env_retain(blacs_env,error=error)
    qs_env%blacs_env => blacs_env
    qs_env%ref_count=1
    last_qs_env_id_nr=last_qs_env_id_nr+1
    qs_env%id_nr=last_qs_env_id_nr
  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl,sca_ppnl,scb_ppnl,&
                        pao_env,particle_set,energy,force,h,k,s,&
                        pw_env,mpools,rho,kinetic,&
                        scf_control,qs_charges,rho_core, ks_env,&
                        wf_history, scf_env, error)

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
                                                 sac_ppnl,sbc_ppnl,sca_ppnl,&
                                                 scb_ppnl
    TYPE(pao_env_type), OPTIONAL, POINTER     :: pao_env
    TYPE(particle_type), DIMENSION(:),&
      OPTIONAL, POINTER                       :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER   :: energy
    TYPE(qs_force_type), DIMENSION(:),&
      OPTIONAL, POINTER                       :: force
    TYPE(real_matrix_p_type), DIMENSION(:),&
         OPTIONAL, POINTER                    :: h,k,s,kinetic
    TYPE(qs_rho_type), POINTER,OPTIONAL       :: rho
    TYPE(pw_env_type), OPTIONAL, POINTER      :: pw_env
    TYPE(scf_control_type), OPTIONAL, POINTER :: scf_control
    TYPE(qs_charges_type), OPTIONAL, POINTER  :: qs_charges
    TYPE(coeff_type), POINTER, OPTIONAL       :: rho_core
    TYPE(qs_ks_env_type), POINTER, OPTIONAL   :: ks_env
    TYPE(qs_wf_history_type), POINTER, OPTIONAL :: wf_history
    TYPE(qs_scf_env_type), POINTER, OPTIONAL  :: scf_env
    TYPE(qs_matrix_pools_type), POINTER, OPTIONAL  :: mpools
    TYPE(cp_error_type),INTENT(inout),OPTIONAL:: error

!   ---------------------------------------------------------------------------

    IF (PRESENT(atomic_kind_set)) qs_env%atomic_kind_set => atomic_kind_set
    IF (PRESENT(cell)) qs_env%cell => cell
    IF (PRESENT(dft_control)) qs_env%dft_control => dft_control
    IF (PRESENT(c)) qs_env%c => c
    IF (PRESENT(sab_orb)) qs_env%sab_orb => sab_orb
    IF (PRESENT(sac_erfc)) qs_env%sac_erfc => sac_erfc
    IF (PRESENT(sac_ppl)) qs_env%sac_ppl => sac_ppl
    IF (PRESENT(sac_ppnl)) qs_env%sac_ppnl => sac_ppnl
    IF (PRESENT(sbc_ppnl)) qs_env%sbc_ppnl => sbc_ppnl
    IF (PRESENT(sca_ppnl)) qs_env%sca_ppnl => sca_ppnl
    IF (PRESENT(scb_ppnl)) qs_env%scb_ppnl => scb_ppnl
    IF (PRESENT(pao_env)) qs_env%pao_env => pao_env
    IF (PRESENT(particle_set)) qs_env%particle_set => particle_set
    IF (PRESENT(energy)) qs_env%energy => energy
    IF (PRESENT(force)) qs_env%force => force
    IF (PRESENT(h)) qs_env%h => h
    IF (PRESENT(k)) qs_env%k => k
    IF (PRESENT(s)) qs_env%s => s
    IF (PRESENT(rho)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho)) CALL qs_rho_retain(rho,error=error)
      CALL qs_rho_release(qs_env%rho,error=error)
      qs_env%rho => rho
    END IF
    IF (PRESENT(pw_env)) THEN ! accept also null pointers?
      CALL pw_env_retain(pw_env,error=error)
      CALL pw_env_release(qs_env%pw_env,error=error)
      qs_env%pw_env => pw_env
    END IF
    IF (PRESENT(scf_control)) THEN ! accept also null pointers?
       CALL scf_c_retain(scf_control,error=error)
       CALL scf_c_release(qs_env%scf_control,error=error)
       qs_env%scf_control => scf_control
    END IF
    IF (PRESENT(kinetic)) qs_env%kinetic => kinetic
    IF (PRESENT(rho_core)) qs_env%rho_core => rho_core
    IF (PRESENT(qs_charges)) THEN
       CALL qs_charges_retain(qs_charges,error=error)
       CALL qs_charges_release(qs_env%qs_charges,error=error)
       qs_env%qs_charges => qs_charges
    END IF
    IF (PRESENT(ks_env)) THEN ! accept also null pointers?
       CALL qs_ks_retain(ks_env, error=error)
       CALL qs_ks_release(qs_env%ks_env, error=error)
       qs_env%ks_env => ks_env
    END IF
    IF (PRESENT(wf_history)) THEN ! accept also null pointers ?
       CALL wfi_retain(wf_history,error=error)
       CALL wfi_release(qs_env%wf_history,error=error)
       qs_env%wf_history => wf_history
    END IF
    IF (PRESENT(scf_env)) THEN ! accept also null pointers ?
       CALL scf_env_retain(scf_env,error=error)
       CALL scf_env_release(qs_env%scf_env, error=error)
       qs_env%scf_env => scf_env
    END IF
    IF (PRESENT(mpools)) THEN
       CALL mpools_retain(mpools,error=error)
       CALL mpools_release(mpools, error=error)
       qs_env%mpools => mpools
    END IF
  END SUBROUTINE set_qs_env

! *****************************************************************************

!!****f* qs_environment_types/qs_env_create [1.0] *
!!
!!   NAME
!!     qs_env_create
!!
!!   SYNOPSIS
!!     Subroutine qs_env_create(qs_env, para_env, blacs_env, error)
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Type(cp_para_env_type), Pointer:: para_env
!!       Type(cp_blacs_env_type), Pointer:: blacs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_create
!!
!!   FUNCTION
!!     allocates and intitializes a qs_env
!!
!!   ARGUMENTS
!!     - qs_env: the object to create
!!     - para_env: the parallel environement for the qs_env
!!     - blacs_env: the blacs env for the qs_env (for full matrixes,
!!       it must be coerent with the para_env)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_create(qs_env,para_env,blacs_env,error)
  TYPE(qs_environment_type), POINTER :: qs_env
  TYPE(cp_para_env_type), POINTER :: para_env
  TYPE(cp_blacs_env_type), POINTER :: blacs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_create',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  ALLOCATE(qs_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL init_qs_env(qs_env,para_env=para_env, blacs_env=blacs_env,&
          error=error)
  END IF
END SUBROUTINE qs_env_create
!***************************************************************************

!!****f* qs_environment_types/qs_env_retain [1.0] *
!!
!!   NAME
!!     qs_env_retain
!!
!!   SYNOPSIS
!!     Subroutine qs_env_retain(qs_env, error)
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_retain
!!
!!   FUNCTION
!!     retains the given qs_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - qs_env: the object to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_retain(qs_env,error)
  TYPE(qs_environment_type), POINTER :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_retain',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     qs_env%ref_count=qs_env%ref_count+1
  END IF
END SUBROUTINE qs_env_retain
!***************************************************************************

!!****f* qs_environment_types/qs_env_release [1.0] *
!!
!!   NAME
!!     qs_env_release
!!
!!   SYNOPSIS
!!     Subroutine qs_env_release(qs_env, error)
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_release
!!
!!   FUNCTION
!!     releases the given qs_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - qs_env: the object to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     12.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_release(qs_env,error)
  TYPE(qs_environment_type), POINTER :: qs_env
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_release',&
       routineP=moduleN//':'//routineN
  INTEGER :: i, stat
  failure=.FALSE.

  IF (ASSOCIATED(qs_env)) THEN
     CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     qs_env%ref_count=qs_env%ref_count-1
     IF (qs_env%ref_count<1) THEN
        IF (ASSOCIATED(qs_env%atomic_kind_set)) &
             CALL deallocate_atomic_kind_set(qs_env%atomic_kind_set)
        IF (ASSOCIATED(qs_env%c)) THEN
           DO i=1,SIZE(qs_env%c)
              CALL deallocate_mo_set(qs_env%c(i)%mo_set)
           END DO
           DEALLOCATE(qs_env%c, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(qs_env%sab_orb)) THEN
           DO i=1,SIZE(qs_env%sab_orb)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sab_orb(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sab_orb,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sac_erfc)) THEN
           DO i=1,SIZE(qs_env%sac_erfc)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sac_erfc(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sac_erfc,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sac_ppl)) THEN
           DO i=1,SIZE(qs_env%sac_ppl)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sac_ppl(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sac_ppl,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sac_ppnl)) THEN
           DO i=1,SIZE(qs_env%sac_ppnl)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sac_ppnl(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sac_ppnl,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sbc_ppnl)) THEN
           DO i=1,SIZE(qs_env%sbc_ppnl)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sbc_ppnl(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sbc_ppnl,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sca_ppnl)) THEN
           DO i=1,SIZE(qs_env%sca_ppnl)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sca_ppnl(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sca_ppnl,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%scb_ppnl)) THEN
           DO i=1,SIZE(qs_env%scb_ppnl)
              CALL deallocate_neighbor_list_set(&
                   qs_env%scb_ppnl(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%scb_ppnl,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

        IF (ASSOCIATED(qs_env%particle_set)) &
             CALL deallocate_particle_set(qs_env%particle_set)
        IF (ASSOCIATED(qs_env%energy)) &
             CALL deallocate_qs_energy(qs_env%energy)
        IF (ASSOCIATED(qs_env%force))&
             CALL deallocate_qs_force(qs_env%force)
        IF (ASSOCIATED(qs_env%h))&
             CALL deallocate_matrix_set(qs_env%h)
        IF (ASSOCIATED(qs_env%k))&
             CALL deallocate_matrix_set(qs_env%k)
        IF (ASSOCIATED(qs_env%s))&
             CALL deallocate_matrix_set(qs_env%s)
        IF (ASSOCIATED(qs_env%kinetic))&
             CALL deallocate_matrix_set(qs_env%kinetic)
        CALL qs_rho_release(qs_env%rho, error=error)
        CALL scf_c_release(qs_env%scf_control, error=error)
        CALL pw_env_release(qs_env%pw_env, error=error)
        CALL cp_para_env_release(qs_env%para_env, error=error)
        CALL cp_blacs_env_release(qs_env%blacs_env, error=error)
        IF (ASSOCIATED(qs_env%rho_core)) &
             CALL coeff_deallocate(qs_env%rho_core)
        CALL qs_charges_release(qs_env%qs_charges, error=error)
        CALL qs_ks_release(qs_env%ks_env, error=error)
        CALL mpools_release(qs_env%mpools,error=error)

!FM        if (associated(qs_env%pao_env)) then ! circular ref
!FM           call pao_env_dealloc_ref(pao_env_dealloc_ref,error=error)
!FM           deallocate(qs_env%pao_env,stat=stat)
!FM           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
!FM        end if

        CALL cp_assert(.FALSE.,cp_warning_level, cp_unimplemented_error_nr,&
             routineP, &
             "qs_env%cell, qs_env%pao_env and qs_env%dft_control "//&
             "not fully deallocated"//&
             CPSourceFileRef,&
             error)
        NULLIFY(qs_env%cell, qs_env%dft_control,qs_env%pao_env)
        CALL wfi_release(qs_env%wf_history,error=error)
        CALL scf_env_release(qs_env%scf_env, error=error)

        DEALLOCATE(qs_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(qs_env)
END SUBROUTINE qs_env_release
!***************************************************************************

!!****f* qs_environment_types/qs_env_reorthogonalize_vectors [1.0] *
!!
!!   NAME
!!     qs_env_reorthogonalize_vectors
!!
!!   SYNOPSIS
!!     Subroutine qs_env_reorthogonalize_vectors(qs_env, v_matrix, ispin,&
!!         n_col, error)
!!       Type(qs_environment_type), Pointer:: qs_env
!!       Type(cp_full_matrix_type), Pointer:: v_matrix
!!       Integer, Intent (IN):: ispin
!!       Integer, Intent (IN), Optional:: n_col
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine qs_env_reorthogonalize_vectors
!!
!!   FUNCTION
!!     reorthogonalizes the mos
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - qs_env: the qs_env in which to orthogonalize
!!     - v: the vectors to orthogonalize
!!     - ispin: the spin of the vector to orthogonalize (unused at the 
!!       moment)
!!     - n_col: number of column of v to orthogonalize
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_env_reorthogonalize_vectors(qs_env, v_matrix, ispin, n_col,&
     error)
  TYPE(qs_environment_type), POINTER :: qs_env
  TYPE(cp_full_matrix_type), POINTER :: v_matrix
  INTEGER, INTENT(in) :: ispin
  INTEGER, INTENT(in), OPTIONAL :: n_col
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, ortho_contains_cholesky
  CHARACTER(len=*), PARAMETER :: routineN='qs_env_reorthogonalize_vectors',&
        routineP=moduleN//':'//routineN
  INTEGER :: my_n_col, nmo
  TYPE(qs_scf_env_type), POINTER :: scf_env
  TYPE(scf_control_type), POINTER :: scf_control
  TYPE(cp_full_matrix_type), POINTER :: sv
  TYPE(cp_fm_pool_type), POINTER :: maxao_maxmo_fm_pool
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s
  TYPE(qs_matrix_pools_type), POINTER :: mpools
  
  failure=.FALSE.
  NULLIFY(scf_env,scf_control,sv,maxao_maxmo_fm_pool,s,mpools)
  
  CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(v_matrix),cp_failure_level,routineP,error,failure)
  
  IF (.NOT. failure) THEN

     CALL cp_fm_get_info(v_matrix,ncol_global=my_n_col,error=error)
     IF (PRESENT(n_col)) my_n_col=n_col
     CALL get_qs_env(qs_env,mpools=mpools,&
          scf_env=scf_env,&
          scf_control=scf_control,&
          s=s,&
          error=error)
     CALL mpools_get(mpools,maxao_maxmo_fm_pool=maxao_maxmo_fm_pool,&
          error=error)
     IF (ASSOCIATED(scf_env)) THEN
        ortho_contains_cholesky=(scf_env%method /= ot_method_nr).AND.&
             (scf_control%use_cholesky .OR. scf_control%use_arpack).AND.&
             ASSOCIATED(scf_env%ortho)
     ELSE
        ortho_contains_cholesky=.FALSE.
     END IF

     IF (ortho_contains_cholesky) THEN
        CALL blacs_make_basis(vmatrix=v_matrix,ncol=my_n_col,&
             ortho=scf_env%ortho,otype='CHOLESKY')
     ELSE
        CALL fm_pool_create_fm(maxao_maxmo_fm_pool,sv,error=error)
        CALL cp_sm_fm_multiply(sparse_matrix=s(1)%matrix,v_in=v_matrix,&
             v_out=sv,ncol=my_n_col,para_env=qs_env%para_env, &
             error=error)
        CALL blacs_make_basis(vmatrix=v_matrix,ncol=my_n_col,&
             ortho=sv,otype='SV')
        CALL fm_pool_give_back_fm(maxao_maxmo_fm_pool,sv,error=error)
     END IF
  END IF
END SUBROUTINE qs_env_reorthogonalize_vectors
!***************************************************************************

END MODULE qs_environment_types
