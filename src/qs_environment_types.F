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
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_create
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_retain,&
                                             dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_sm_pool_types,                ONLY: cp_sm_pool_p_type,&
                                             sm_pools_give_back_matrix_vect
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsys_set,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE kinds,                           ONLY: dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE pao_types,                       ONLY: pao_env_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE qs_charges_types,                ONLY: qs_charges_release,&
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
  USE qs_mo_methods,                   ONLY: make_basis
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
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE sparse_matrix_types,             ONLY: deallocate_matrix_set,&
                                             real_matrix_p_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='qs_environment_types'
  INTEGER, PRIVATE, SAVE :: last_qs_env_id_nr=0

! *** Public data types ***

  PUBLIC :: qs_environment_type, qs_environment_p_type

! *** Public subroutines ***

  PUBLIC :: get_qs_env,&
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
!!     contains all the info needed by quickstep to do dft scf
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - id_nr: identification number, guarateed to be unique
!!     - ref_count: reference count (see doc/ReferenceCounting.html)
!!     - cell: the simulation cell (box)
!!     - dft_control: object that contains the values of various parameters
!!       relevant to a dft calculation
!!     - c: the molecular orbitals, in LSD array with the different spins
!!     - sab_orb,sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl:
!!       neighbor lists for the calculation of the core Hamiltonian matrix
!!     - pao_env: things needed by the pao method
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
!!     - distribution_2d: distribution of the atom pairs between the processors
!!     - subsys: the particles, molecules,... of this environment
!!
!!     compatibility get (things that you should get from the subsys):
!!     - atomic_kind_set: array with infos about the species (atomic_kinds)
!!       present in the system
!!     - particle_set: info on the atoms you simulate, pos,...
!!     - local_particles: which particles ar local to this processor
!!     new:
!!     - local_molecules: which molecules are local to this processor
!!     - molecule_kind_set: description of the molecule kinds
!!     - molecule_set: all the molecule description
!!     
!!   AUTHOR
!!     Matthias Krack & fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2002 added doc and attribute description [fawzi]
!!
!!*** *********************************************************************

  TYPE qs_environment_type
!   PRIVATE
    INTEGER :: id_nr, ref_count
    TYPE(cell_type), POINTER                              :: cell
    TYPE(dft_control_type), POINTER                       :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER            :: c
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_erfc,&
                                                             sac_ppl,sac_ppnl,&
                                                             sbc_ppnl
    TYPE(pao_env_type), POINTER                           :: pao_env
    TYPE(qs_energy_type), POINTER                         :: energy
    TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: h,k,s,kinetic,w
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
    TYPE(distribution_2d_type), POINTER                   :: distribution_2d
    TYPE(cp_subsystem_type), POINTER                      :: subsys
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
                        sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl,&
                        pao_env,particle_set,energy,force,h,k,s,w,&
                        rho, pw_env, mpools,&
                        para_env, blacs_env, scf_control,&
                        kinetic,qs_charges, rho_core, ks_env,&
                        wf_history,scf_env,&
                        id_nr,local_particles,local_molecules,&
                        distribution_2d,molecule_kind_set,molecule_set,&
                        subsys,error)

!   Purpose: Get the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(in), &
      TARGET                                 :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: c
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_orb, sac_erfc, sac_ppl, &
                                                sac_ppnl, sbc_ppnl
    TYPE(pao_env_type), OPTIONAL, POINTER    :: pao_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: h, k, s, w
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env
    TYPE(qs_matrix_pools_type), OPTIONAL, &
      POINTER                                :: mpools
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                :: scf_control
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: kinetic
    TYPE(qs_charges_type), OPTIONAL, POINTER :: qs_charges
    TYPE(coeff_type), OPTIONAL, POINTER      :: rho_core
    TYPE(qs_ks_env_type), OPTIONAL, POINTER  :: ks_env
    TYPE(qs_wf_history_type), OPTIONAL, &
      POINTER                                :: wf_history
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    INTEGER, INTENT(out), OPTIONAL           :: id_nr
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_qs_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(atomic_kinds, particles, molecules_new,&
         molecule_kinds_new)
!    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env%subsys),cp_failure_level,routineP,error,failure)

    IF (PRESENT(cell)) cell => qs_env%cell
    IF (PRESENT(dft_control)) dft_control => qs_env%dft_control
    IF (PRESENT(c)) c => qs_env%c
    IF (PRESENT(sab_orb)) sab_orb => qs_env%sab_orb
    IF (PRESENT(sac_erfc)) sac_erfc => qs_env%sac_erfc
    IF (PRESENT(sac_ppl)) sac_ppl => qs_env%sac_ppl
    IF (PRESENT(sac_ppnl)) sac_ppnl => qs_env%sac_ppnl
    IF (PRESENT(sbc_ppnl)) sbc_ppnl => qs_env%sbc_ppnl
    IF (PRESENT(pao_env)) pao_env => qs_env%pao_env
    IF (PRESENT(energy)) energy => qs_env%energy
    IF (PRESENT(force)) force => qs_env%force
    IF (PRESENT(h)) h => qs_env%h
    IF (PRESENT(k)) k => qs_env%k
    IF (PRESENT(kinetic)) kinetic => qs_env%kinetic
    IF (PRESENT(s)) s => qs_env%s
    IF (PRESENT(w)) w => qs_env%w
    IF (PRESENT(rho)) rho => qs_env%rho
    IF (PRESENT(pw_env)) pw_env => qs_env%pw_env
    IF (PRESENT(mpools)) mpools => qs_env%mpools
    IF (PRESENT(blacs_env)) blacs_env => qs_env%blacs_env
    IF (PRESENT(scf_control)) scf_control => qs_env%scf_control
    IF (PRESENT(rho_core)) rho_core => qs_env%rho_core
    IF (PRESENT(qs_charges)) qs_charges => qs_env%qs_charges
    IF (PRESENT(ks_env)) ks_env => qs_env%ks_env
    IF (PRESENT(wf_history)) wf_history => qs_env%wf_history
    IF (PRESENT(scf_env)) scf_env => qs_env%scf_env
    IF (PRESENT(id_nr)) id_nr=qs_env%id_nr
    IF (PRESENT(subsys)) subsys => qs_env%subsys
    IF (PRESENT(para_env)) para_env => qs_env%para_env
    IF (PRESENT(distribution_2d)) distribution_2d => qs_env%distribution_2d

    CALL cp_subsys_get(qs_env%subsys,&
                      atomic_kinds=atomic_kinds,&
                      local_molecules_new=local_molecules,&
                      local_particles=local_particles,&
                      particles=particles,&
                      molecule_kinds_new=molecule_kinds_new,&
                      molecules_new=molecules_new,&
                      error=error)
    IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
    IF (PRESENT(particle_set)) particle_set => particles%els
    IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds_new%els
    IF (PRESENT(molecule_set)) molecule_set => molecules_new%els

  END SUBROUTINE get_qs_env

! *****************************************************************************

  SUBROUTINE init_qs_env(qs_env,para_env, error)

!   Purpose: Initialise the QUICKSTEP environment.

!   History: - Creation (25.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)   :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

!   ---------------------------------------------------------------------------

    NULLIFY (qs_env%cell)
    NULLIFY (qs_env%dft_control)
    NULLIFY (qs_env%c)
    NULLIFY (qs_env%sab_orb)
    NULLIFY (qs_env%sac_erfc)
    NULLIFY (qs_env%sac_ppl)
    NULLIFY (qs_env%sac_ppnl)
    NULLIFY (qs_env%sbc_ppnl)
    NULLIFY (qs_env%pao_env)
    NULLIFY (qs_env%energy)
    NULLIFY (qs_env%force)
    NULLIFY (qs_env%h)
    NULLIFY (qs_env%k)
    NULLIFY (qs_env%kinetic)
    NULLIFY (qs_env%s)
    NULLIFY (qs_env%w)
    NULLIFY (qs_env%mpools)
    NULLIFY (qs_env%rho)
    NULLIFY (qs_env%pw_env)
    NULLIFY (qs_env%scf_control)
    NULLIFY (qs_env%qs_charges)
    NULLIFY (qs_env%rho_core)
    NULLIFY (qs_env%ks_env)
    NULLIFY (qs_env%wf_history)
    NULLIFY (qs_env%scf_env)
    NULLIFY (qs_env%distribution_2d)
    NULLIFY (qs_env%subsys)

    CALL cp_subsys_create(qs_env%subsys, para_env=para_env, error=error)
    CALL cp_para_env_retain(para_env,error=error)
    qs_env%para_env => para_env
    CALL cp_blacs_env_create(qs_env%blacs_env,qs_env%para_env,error=error)

    qs_env%ref_count=1
    last_qs_env_id_nr=last_qs_env_id_nr+1
    qs_env%id_nr=last_qs_env_id_nr
  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,atomic_kind_set,cell,dft_control,c,&
                        sab_orb,sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl,&
                        pao_env,particle_set,energy,force,h,k,s,w,&
                        pw_env,mpools,rho,kinetic,&
                        scf_control,qs_charges,rho_core, ks_env,&
                        wf_history,scf_env,local_particles,local_molecules,&
                        distribution_2d,molecule_kind_set,molecule_set,&
                        subsys, error)

!   Purpose: Set the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), &
      INTENT(inout), TARGET                  :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: c
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_orb, sac_erfc, sac_ppl, &
                                                sac_ppnl, sbc_ppnl
    TYPE(pao_env_type), OPTIONAL, POINTER    :: pao_env
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: h, k, s, w
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env
    TYPE(qs_matrix_pools_type), OPTIONAL, &
      POINTER                                :: mpools
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: kinetic
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                :: scf_control
    TYPE(qs_charges_type), OPTIONAL, POINTER :: qs_charges
    TYPE(coeff_type), OPTIONAL, POINTER      :: rho_core
    TYPE(qs_ks_env_type), OPTIONAL, POINTER  :: ks_env
    TYPE(qs_wf_history_type), OPTIONAL, &
      POINTER                                :: wf_history
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "set_qs_env", &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles

!   ---------------------------------------------------------------------------

    failure=.FALSE.
!    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(qs_env%cell,error=error)
       qs_env%cell => cell
    END IF
    IF (PRESENT(dft_control)) THEN
       CALL dft_control_retain(dft_control, error=error)
       CALL dft_control_release(qs_env%dft_control, error=error)
       qs_env%dft_control => dft_control
    END IF
    IF (PRESENT(c)) qs_env%c => c
    IF (PRESENT(sab_orb)) qs_env%sab_orb => sab_orb
    IF (PRESENT(sac_erfc)) qs_env%sac_erfc => sac_erfc
    IF (PRESENT(sac_ppl)) qs_env%sac_ppl => sac_ppl
    IF (PRESENT(sac_ppnl)) qs_env%sac_ppnl => sac_ppnl
    IF (PRESENT(sbc_ppnl)) qs_env%sbc_ppnl => sbc_ppnl
    IF (PRESENT(pao_env)) qs_env%pao_env => pao_env
    IF (PRESENT(energy)) qs_env%energy => energy
    IF (PRESENT(force)) qs_env%force => force
    IF (PRESENT(h)) qs_env%h => h
    IF (PRESENT(k)) qs_env%k => k
    IF (PRESENT(s)) qs_env%s => s
    IF (PRESENT(w)) qs_env%w => w
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
       CALL mpools_release(qs_env%mpools, error=error) 
       qs_env%mpools => mpools
    END IF
    IF (PRESENT(subsys)) THEN
       CALL cp_subsys_retain(subsys,error=error)
       CALL cp_subsys_release(qs_env%subsys,error=error)
       qs_env%subsys => subsys
    END IF
    IF (PRESENT(distribution_2d)) THEN
       CALL distribution_2d_retain(distribution_2d)
       CALL distribution_2d_release(qs_env%distribution_2d)
       qs_env%distribution_2d => distribution_2d
    ENDIF

    IF (PRESENT(atomic_kind_set)) THEN
      CALL atomic_kind_list_create(atomic_kinds,&
                                   els_ptr=atomic_kind_set,&
                                   error=error)
      CALL cp_subsys_set(qs_env%subsys,&
                        atomic_kinds=atomic_kinds,&
                        error=error)
      CALL atomic_kind_list_release(atomic_kinds,error=error)
    END IF
    IF (PRESENT(particle_set)) THEN
      CALL particle_list_create(particles,&
                                els_ptr=particle_set,&
                                error=error)
      CALL cp_subsys_set(qs_env%subsys,&
                        particles=particles,&
                        error=error)
      CALL particle_list_release(particles,error=error)
    END IF
    IF (PRESENT(local_particles)) THEN
      CALL cp_subsys_set(qs_env%subsys,&
                        local_particles=local_particles,&
                        error=error)
    END IF
    IF (PRESENT(local_molecules)) THEN
      CALL cp_subsys_set(qs_env%subsys,&
                        local_molecules_new=local_molecules,&
                        error=error)
    END IF
    IF (PRESENT(molecule_kind_set)) THEN
      CALL mol_kind_new_list_create(molecule_kinds_new,&
                                    els_ptr=molecule_kind_set,&
                                    error=error)
      CALL cp_subsys_set(qs_env%subsys,&
                        molecule_kinds_new=molecule_kinds_new,&
                        error=error)
      CALL mol_kind_new_list_release(molecule_kinds_new,error=error)
    END IF
    IF (PRESENT(molecule_set)) THEN
      CALL mol_new_list_create(molecules_new,&
                               els_ptr=molecule_set,&
                               error=error)
      CALL cp_subsys_set(qs_env%subsys,&
                        molecules_new=molecules_new,&
                        error=error)
      CALL mol_new_list_release(molecules_new,error=error)
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
SUBROUTINE qs_env_create(qs_env,para_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(qs_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL init_qs_env(qs_env,para_env=para_env, error=error)
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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(cp_sm_pool_p_type), DIMENSION(:), &
      POINTER                                :: S_sm_pools

  NULLIFY(S_sm_pools)
  failure=.FALSE.

  IF (ASSOCIATED(qs_env)) THEN
     CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
     qs_env%ref_count=qs_env%ref_count-1
     IF (qs_env%ref_count<1) THEN
        CALL cell_release(qs_env%cell)
        CALL dft_control_release(qs_env%dft_control, error=error)
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
        IF (ASSOCIATED(qs_env%energy)) &
             CALL deallocate_qs_energy(qs_env%energy)
        IF (ASSOCIATED(qs_env%force))&
             CALL deallocate_qs_force(qs_env%force)
        IF (ASSOCIATED(qs_env%h)) CALL deallocate_matrix_set(qs_env%h)

        IF (ASSOCIATED(qs_env%k)) THEN
            CALL mpools_get(qs_env%mpools, S_sm_pools=S_sm_pools, error=error)
            CALL sm_pools_give_back_matrix_vect(S_sm_pools,qs_env % k,error=error)
        ENDIF 
        IF (ASSOCIATED(qs_env%s)) CALL deallocate_matrix_set(qs_env%s)
        IF (ASSOCIATED(qs_env%w)) CALL deallocate_matrix_set(qs_env%w)
        IF (ASSOCIATED(qs_env%kinetic))&
             CALL deallocate_matrix_set(qs_env%kinetic)
        CALL qs_rho_release(qs_env%rho, error=error)
        CALL scf_c_release(qs_env%scf_control, error=error)
        CALL pw_env_release(qs_env%pw_env, error=error)
        CALL cp_para_env_release(qs_env%para_env,error=error)
        CALL cp_blacs_env_release(qs_env%blacs_env, error=error)
        IF (ASSOCIATED(qs_env%rho_core)) THEN
             CALL coeff_deallocate(qs_env%rho_core)
             DEALLOCATE(qs_env%rho_core)
        ENDIF
        CALL qs_charges_release(qs_env%qs_charges, error=error)
        CALL qs_ks_release(qs_env%ks_env, error=error)
        CALL wfi_release(qs_env%wf_history,error=error)
        CALL scf_env_release(qs_env%scf_env, error=error)
        CALL mpools_release(qs_env%mpools,error=error)
        IF (ASSOCIATED(qs_env%distribution_2d)) THEN
           CALL distribution_2d_release(qs_env%distribution_2d)
        ENDIF
        CALL cp_subsys_release(qs_env%subsys,error=error)
        IF (ASSOCIATED(qs_env%pao_env)) THEN ! circular ref
           NULLIFY(qs_env%pao_env)
!FM           call pao_env_dealloc_ref(pao_env_dealloc_ref,error=error)
!FM           deallocate(qs_env%pao_env,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF

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
!!       Type(cp_fm_type), Pointer:: v_matrix
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
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: v_matrix
    INTEGER, INTENT(in)                      :: ispin
    INTEGER, INTENT(in), OPTIONAL            :: n_col
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'qs_env_reorthogonalize_vectors', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, my_n_col
    LOGICAL                                  :: failure, &
                                                ortho_contains_cholesky
    TYPE(cp_fm_pool_type), POINTER           :: maxao_maxmo_fm_pool
    TYPE(cp_fm_type), POINTER                :: sv
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: s
    TYPE(scf_control_type), POINTER          :: scf_control

  failure=.FALSE.
  NULLIFY(scf_env,scf_control,sv,maxao_maxmo_fm_pool,s,mpools)
  CALL timeset(routineN,'I',' ',handle)
  
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
             (scf_control%use_cholesky ).AND.&
             ASSOCIATED(scf_env%ortho)
     ELSE
        ortho_contains_cholesky=.FALSE.
     END IF

     IF (ortho_contains_cholesky) THEN
        CALL make_basis(vmatrix=v_matrix,ncol=my_n_col,&
             ortho=scf_env%ortho,otype='CHOLESKY')
     ELSE
        ! CALL fm_pool_create_fm(maxao_maxmo_fm_pool,sv,error=error)
        CALL cp_fm_create(sv,v_matrix%matrix_struct)
        CALL cp_sm_fm_multiply(sparse_matrix=s(1)%matrix,v_in=v_matrix,&
             v_out=sv,ncol=my_n_col,para_env=qs_env%para_env, &
             error=error)
        CALL make_basis(vmatrix=v_matrix,ncol=my_n_col,&
             ortho=sv,otype='SV')
        ! CALL fm_pool_give_back_fm(maxao_maxmo_fm_pool,sv,error=error)
        CALL cp_fm_release(sv)
     END IF
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE qs_env_reorthogonalize_vectors
!***************************************************************************

END MODULE qs_environment_types
