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
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE coefficient_types,               ONLY: coeff_type
  USE cp_b_matrix_structure,           ONLY: cp_b_matrix_struct_type,&
                                             cp_bmstruct_create,&
                                             cp_bmstruct_release
  USE cp_blacs_env,                    ONLY: cp_blacs_env_retain
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
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                             cp_fm_pool_type,&
                                             fm_pool_create,&
                                             fm_pool_release,&
                                             fm_pool_retain,&
                                             fm_pools_copy,&
                                             fm_pools_dealloc
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_retain,&
                                             cp_fm_struct_type
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
  USE pw_env_types,                    ONLY: pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_retain,&
                                             qs_charges_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
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

! *** Public data types ***

  PUBLIC :: qs_environment_type, qs_environment_p_type

! *** Public subroutines ***

  PUBLIC :: get_qs_env,&
            init_qs_env,&
            set_qs_env

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
!!     - atomic_kind_set: array with infos about the species (atomic_kinds)
!!       present in the system
!!     - cell: the simulation cell (box)
!!     - dft_control: object that contains the values of various parameters
!!       relevant to a dft calculation
!!     - c: the molecular orbitals, in LSD array with the different spins
!!     - sab_orb,sac_erfc,sac_ppl,sca_ppnl,scb_ppnl : various neighbor
!!       lists (Matthias expand this??)
!!     - pao_env: things nedeed by the pao method
!!     - particle_set: info on the atoms you simulate, pos,...
!!     - energy: place where the energies calculated are stored
!!     - h: core hamiltonian (h(1)), if needed by the calculation also
!!       its derivatives wrt. x,y, and z (h(2:4))
!!     - k: Khon-Sham matrix (filtred by the structure of S).
!!       with LSD an array with the differents orbitals
!!     - p: the density matrix (filtred by the structure of S).
!!       with LSD an array with the differents orbitals
!!     - s: the overlap matrix (s(1)), if needed by the calculation also
!!       its derivatives wrt. x,y, and z (s(2:4))
!!     - kinetic: kinetic part of h
!!     - s_mstruct: structure of s
!!     - rho: the (old) density, in various representations (ao+grid)
!!     - scf_control: various parameters relevant to the scf cycle
!!     - ao_mo_fm_pools: pools with (ao x mo) full matrixes (same order as
!!       c).
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
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: h,k,p,s,kinetic
    TYPE(cp_b_matrix_struct_type), POINTER                :: s_mstruct
    TYPE(qs_rho_type), POINTER                            :: rho
    TYPE(scf_control_type), POINTER                       :: scf_control
    TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER        :: ao_mo_fm_pools
    TYPE(pw_env_type), POINTER                            :: pw_env
    TYPE(cp_para_env_type), POINTER                       :: para_env
    TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
    TYPE(coeff_type), POINTER                             :: rho_core
    TYPE(qs_charges_type), POINTER                        :: qs_charges
    TYPE(qs_ks_env_type), POINTER                         :: ks_env
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
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,s, s_mstruct,&
                        rho, ao_mo_fm_pools,ao_maxmo_fm_pool, pw_env, &
                        para_env, blacs_env, scf_control,&
                        kinetic,qs_charges, rho_core, ks_env,error)

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
      OPTIONAL, POINTER                         :: h,k,p,s,kinetic
    TYPE(cp_b_matrix_struct_type), OPTIONAL,&
         POINTER                                :: s_mstruct
    TYPE(qs_rho_type), OPTIONAL, POINTER        :: rho
    TYPE(cp_fm_pool_p_type), DIMENSION(:),&
         OPTIONAL, POINTER                      :: ao_mo_fm_pools
    TYPE(pw_env_type), OPTIONAL, POINTER        :: pw_env
    TYPE(cp_fm_pool_type), OPTIONAL, POINTER    :: ao_maxmo_fm_pool
    TYPE(cp_para_env_type), OPTIONAL, POINTER   :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, POINTER  :: blacs_env
    TYPE(scf_control_type), OPTIONAL, POINTER   :: scf_control
    TYPE(qs_charges_type), OPTIONAL, POINTER    :: qs_charges
    TYPE(coeff_type), OPTIONAL, POINTER         :: rho_core
    TYPE(qs_ks_env_type), OPTIONAL, POINTER     :: ks_env
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
    IF (PRESENT(s)) s => qs_env%s
    IF (PRESENT(s_mstruct)) s_mstruct => qs_env%s_mstruct
    IF (PRESENT(rho)) rho => qs_env%rho
    IF (PRESENT(pw_env)) pw_env => qs_env%pw_env
    IF (PRESENT(ao_mo_fm_pools)) ao_mo_fm_pools => qs_env%ao_mo_fm_pools
    IF (PRESENT(ao_maxmo_fm_pool)) THEN
       IF (ASSOCIATED(qs_env%ao_mo_fm_pools)) THEN
          !TC what is this hardcoded 1??
          ao_maxmo_fm_pool => qs_env%ao_mo_fm_pools(1)%pool
       ELSE
          NULLIFY(ao_maxmo_fm_pool) ! raise an error?
       END IF
    END IF
    IF (PRESENT(para_env)) para_env => qs_env%para_env
    IF (PRESENT(blacs_env)) blacs_env => qs_env%blacs_env
    IF (PRESENT(scf_control)) scf_control => qs_env%scf_control
    IF (PRESENT(rho_core)) rho_core => qs_env%rho_core
    IF (PRESENT(qs_charges)) qs_charges => qs_env%qs_charges
    IF (PRESENT(ks_env)) ks_env => qs_env%ks_env

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
    NULLIFY (qs_env%s)
    NULLIFY (qs_env%s_mstruct)
    NULLIFY (qs_env%rho)
    NULLIFY (qs_env%ao_mo_fm_pools)
    NULLIFY (qs_env%pw_env)
    NULLIFY (qs_env%scf_control)
    NULLIFY (qs_env%qs_charges)
    NULLIFY (qs_env%rho_core)
    NULLIFY (qs_env%ks_env)

    CALL cp_para_env_retain(para_env,error=error)
    qs_env%para_env => para_env
    CALL cp_blacs_env_retain(blacs_env,error=error)
    qs_env%blacs_env => blacs_env
  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,sab_orb,&
                        sac_erfc,sac_ppl,sca_ppnl,scb_ppnl,pao_env,&
                        particle_set,energy,force,h,k,p,s,&
                        pw_env,ao_mo_fm_pools,rho,kinetic,&
                        scf_control,qs_charges,rho_core, ks_env, error)

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
         OPTIONAL, POINTER                    :: h,k,p,s,kinetic
    TYPE(qs_rho_type), POINTER,OPTIONAL       :: rho
    TYPE(cp_fm_pool_p_type), DIMENSION(:),&
         OPTIONAL, POINTER                    :: ao_mo_fm_pools
    TYPE(pw_env_type), OPTIONAL, POINTER      :: pw_env
    TYPE(scf_control_type), OPTIONAL, POINTER :: scf_control
    TYPE(qs_charges_type), OPTIONAL, POINTER  :: qs_charges
    TYPE(coeff_type), POINTER, OPTIONAL       :: rho_core
    TYPE(qs_ks_env_type), POINTER, OPTIONAL   :: ks_env
    TYPE(cp_error_type),INTENT(inout),OPTIONAL:: error

    TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER :: new_ao_mo_fm_pools

!   ---------------------------------------------------------------------------

    NULLIFY(new_ao_mo_fm_pools)

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
    IF (PRESENT(ao_mo_fm_pools)) THEN
      IF (ASSOCIATED(ao_mo_fm_pools)) THEN
        CALL fm_pools_copy(ao_mo_fm_pools,new_ao_mo_fm_pools,error=error)
      END IF
      CALL fm_pools_dealloc(qs_env%ao_mo_fm_pools,error=error)
      qs_env%ao_mo_fm_pools => new_ao_mo_fm_pools
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
  END SUBROUTINE set_qs_env

! *****************************************************************************

END MODULE qs_environment_types
