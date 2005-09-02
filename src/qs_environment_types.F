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
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE coefficient_types,               ONLY: coeff_deallocate,&
                                             coeff_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_retain,&
                                             dft_control_type,&
                                             tddfpt_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_type, &
                                             cp_fm_p_type
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
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE harris_env_types,                ONLY: harris_env_type, &
                                             harris_env_release
  USE hartree_local_types,             ONLY: ecoul_1center_type,&
                                             get_hartree_local,&
                                             hartree_local_create,&
                                             hartree_local_release,&
                                             hartree_local_type,&
                                             set_hartree_local
  USE input_section_types,             ONLY: section_vals_type, section_vals_retain,&
                                             section_vals_release
  USE kg_gpw_pw_env_types,             ONLY: kg_sub_pw_env_release,&
                                             kg_sub_pw_env_retain,&
                                             kg_sub_pw_env_type
  USE kinds,                           ONLY: dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_type
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
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type,&
                                             qs_ks_qmmm_release,&
                                             qs_ks_qmmm_retain
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
  USE qs_linres_types,                 ONLY: linres_control_release,&
                                             linres_control_retain,&
                                             linres_control_type
  USE qs_local_rho_types,              ONLY: local_rho_set_create,&
                                             local_rho_set_release,&
                                             get_local_rho,&
                                             local_rho_type,&
                                             rhoz_type,&
                                             set_local_rho
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             mpools_release,&
                                             mpools_retain,&
                                             qs_matrix_pools_type
  USE qs_mo_methods,                   ONLY: make_basis
  USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE qs_oce_types,                    ONLY: deallocate_oce_set,&
                                             oce_matrix_type,&
                                             vtriple_type
  USE qs_rho_atom_types,               ONLY: deallocate_rho_atom_set, &
                                             rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_retain,&
                                             qs_rho_type
  USE qs_rho0_types,                   ONLY: deallocate_rho0_atom,&
                                             deallocate_rho0_mpole,&
                                             rho0_atom_type,&
                                             rho0_mpole_type
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
  USE green_types,                     ONLY: greens_function_type,&
                                             pw_green_release,&
                                             pw_green_retain,&
                                             pw_green_create
  USE xas_env_types,                   ONLY: xas_env_release,&
                                             xas_env_retain,& 
                                             xas_environment_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: moduleN='qs_environment_types'
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
!!     - mos: the molecular orbitals, in LSD array with the different spins
!!     - sab_orb,sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl:
!!       neighbor lists for the calculation of the core Hamiltonian matrix
!!     - sac_oce,sbc_oce: neighbor lists for the calculation of the expansion
!!       coefficients of the local atomic densities rho1_hard and rho1_soft
!!     - sab_2c: neighbor lists for the calculation of the 2 centers
!!       hartree term Vh(rho0_hard)rho0_hard-Vh(rho0_soft)*rho0_soft
!!     - sac_3c, sbc_3c: neighbor lists for the calculation of the 3 centers
!!       hartree term Vh(rho0_hard-rho0_soft)*rho_soft
!!     - energy: place where the energies calculated are stored
!!     - matrix_h: core hamiltonian (matrix_h(1)), if needed by the calculation also
!!       its derivatives wrt. x,y, and z (matrix_h(2:4))
!!     - matrix_ks: Khon-Sham matrix (filtred by the structure of S).
!!       with LSD an array with the differents orbitals
!!     - matrix_s: the overlap matrix (matrix_s(1)), if needed by the calculation also
!!       its derivatives wrt. x,y, and z (matrix_s(2:4))
!!     - kinetic: kinetic part of h
!!     - rho: the (old) density, in various representations (ao+grid)
!!     - rho_xc: the (old) soft density, used for the GAPW_XC method only
!!     - scf_control: various parameters relevant to the scf cycle
!!     - pw_env: various things to create and get info on the grids used
!!       by qs
!!     - kg_sub_pw_env : it contains the set of pw_env for the different molecule kind
!!       to be used in the KG_GPW method in order to construct the molecular desities
!!       on different grids. The grids are built in the small box containing the single molecule. 
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
!!     - local_rho_set: contains the atomic, compensations and core densities
!!                      and the local parts of the xc terms
!!     - hartree_local: contains the 1, 2 and 3 centers coulomb terms 
!!     - requires_mo_derivs: logical, true if dE/dC is required (e.g. OT) 
!!     - mo_derivs: the actual derivatives of the total energy wrt to MO coeffs
!!     - neighbor_list_id: the current version of the neighbor_list
!!     - xas_env: temporary information for xas calculation
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
!!     - harris_env: The harris environment, which consist of the harris energy
!!                   and force type and the rho structure to pass the density of
!!                   the harris energy functional to the force calculation
!!     
!!   AUTHOR
!!     Matthias Krack & fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2002 added doc and attribute description [fawzi]
!!     08.2004 renamed some of the very short names (s,c,k,h) for easier grepping
!!
!!*** *********************************************************************

  TYPE qs_environment_type
!   PRIVATE
    INTEGER :: id_nr, ref_count
    LOGICAL :: use_ref_cell
    LOGICAL :: qmmm
    LOGICAL :: requires_mo_derivs
    TYPE(cell_type), POINTER                              :: cell, super_cell 
    TYPE(cell_type), POINTER                              :: cell_ref
    TYPE(greens_function_type), POINTER                   :: green
    TYPE(dft_control_type), POINTER                       :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER            :: mos
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER             :: mo_derivs
    INTEGER                                               :: neighbor_list_id
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_erfc,&
                                                             sac_ppl,sac_ppnl,&
                                                             sbc_ppnl,sac_oce,&
                                                             sbc_oce,sab_2c,&
                                                             sac_3c,sbc_3c
    TYPE(qs_energy_type), POINTER                         :: energy
    TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
    TYPE(harris_env_type), POINTER                        :: harris_env
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER       :: matrix_h,&
                                                             matrix_ks,&
                                                             matrix_s,&
                                                             kinetic,&
                                                             matrix_w
    TYPE(qs_rho_type), POINTER                            :: rho
    TYPE(qs_rho_type), POINTER                            :: rho_XC
    TYPE(scf_control_type), POINTER                       :: scf_control
    TYPE(pw_env_type), POINTER                            :: pw_env
    TYPE(cp_para_env_type), POINTER                       :: para_env
    TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
    TYPE(coeff_type), POINTER                             :: rho_core
    TYPE(qs_charges_type), POINTER                        :: qs_charges
    TYPE(qs_ks_env_type), POINTER                         :: ks_env
    TYPE(qs_ks_qmmm_env_type), POINTER                    :: ks_qmmm_env
    TYPE(qs_wf_history_type), POINTER                     :: wf_history
    TYPE(qs_scf_env_type), POINTER                        :: scf_env
    TYPE(qs_matrix_pools_type), POINTER                   :: mpools
    TYPE(distribution_2d_type), POINTER                   :: distribution_2d
    TYPE(cp_subsystem_type), POINTER                      :: subsys
    TYPE(oce_matrix_type), POINTER                        :: oce
    TYPE(local_rho_type), POINTER                         :: local_rho_set
    TYPE(hartree_local_type),  POINTER                    :: hartree_local
    TYPE(section_vals_type), POINTER                      :: input
    TYPE(kg_sub_pw_env_type), POINTER                     :: kg_sub_pw_env
    TYPE(linres_control_type), POINTER                    :: linres_control
    TYPE(xas_environment_type), POINTER                   :: xas_env
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

  SUBROUTINE get_qs_env(qs_env,atomic_kind_set,cell,super_cell,cell_ref,use_ref_cell,&
                        dft_control,mos,sab_orb,qmmm,green,&
                        sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl,&
                        sac_oce,sbc_oce,sab_2c,sac_3c,sbc_3c,&
                        particle_set,energy,force,harris_env,&
                        matrix_h,matrix_ks,matrix_s,matrix_w,&
                        rho, rho_xc, pw_env, mpools,input,&
                        para_env, blacs_env, scf_control,&
                        kinetic,qs_charges, rho_core, ks_env, ks_qmmm_env,&
                        wf_history,scf_env,&
                        id_nr,local_particles,local_molecules,&
                        distribution_2d,molecule_kind_set,molecule_set,&
                        subsys,error,oce,rho_atom_set,&
                        rho0_atom_set,rho0_mpole,nddo_mpole,&
                        rhoz_set,ecoul_1c,VH_3c_list,rho0_s_rs,rho0_s_gs, &
                        requires_mo_derivs,mo_derivs,tddfpt_control, &
                        neighbor_list_id,kg_sub_pw_env,linres_control,xas_env)

!   Purpose: Get the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(greens_function_type), POINTER, OPTIONAL  :: green
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref, super_cell
    LOGICAL, OPTIONAL                        :: use_ref_cell
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(section_vals_type), POINTER, optional :: input
    TYPE(mo_set_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_orb, sac_erfc, sac_ppl, &
                                                sac_ppnl, sbc_ppnl,&
                                                sac_oce,sbc_oce,sab_2c,&
                                                sac_3c,sbc_3c
    LOGICAL, OPTIONAL                        :: qmmm
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(harris_env_type), OPTIONAL, &
       POINTER                               :: harris_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h, matrix_ks, matrix_s, matrix_w
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho_xc
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
    TYPE(qs_ks_qmmm_env_type), OPTIONAL, &
      POINTER                                :: ks_qmmm_env
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
    TYPE(oce_matrix_type), OPTIONAL, POINTER :: oce
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_atom_set
    TYPE(rho0_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho0_atom_set
    TYPE(rho0_mpole_type), OPTIONAL, POINTER :: rho0_mpole
    TYPE(rho0_mpole_type), OPTIONAL, POINTER :: nddo_mpole
    TYPE(rhoz_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rhoz_set
    TYPE(ecoul_1center_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ecoul_1c
    TYPE(vtriple_type), DIMENSION(:,:), &
      OPTIONAL, POINTER                      :: VH_3c_list
    TYPE(coeff_type), OPTIONAL, POINTER      :: rho0_s_rs,&
                                                rho0_s_gs
    LOGICAL, OPTIONAL                        :: requires_mo_derivs
    TYPE(cp_fm_p_type), OPTIONAL, POINTER,&
      DIMENSION(:)                           :: mo_derivs 
    TYPE(tddfpt_control_type), &
         OPTIONAL, POINTER                   :: tddfpt_control
    INTEGER, OPTIONAL                        :: neighbor_list_id
    TYPE(kg_sub_pw_env_type), OPTIONAL,&
      POINTER                                :: kg_sub_pw_env
    TYPE(linres_control_type), OPTIONAL,&
      POINTER                                :: linres_control
    TYPE(xas_environment_type), OPTIONAL,&
      POINTER                                :: xas_env


    CHARACTER(len=*), PARAMETER :: routineN = 'get_qs_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(rho0_mpole_type), POINTER           :: rho0_m

!   ---------------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(atomic_kinds, particles, molecules_new,&
         molecule_kinds_new)
!    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env%subsys),cp_failure_level,routineP,error,failure)

    IF (PRESENT(cell)) cell => qs_env%cell
    IF (PRESENT(super_cell)) super_cell => qs_env%super_cell
    IF (PRESENT(cell_ref)) cell_ref => qs_env%cell_ref
    IF (PRESENT(green)) green => qs_env%green
    IF (PRESENT(use_ref_cell)) use_ref_cell = qs_env%use_ref_cell
    IF (PRESENT(qmmm)) qmmm = qs_env%qmmm
    IF (PRESENT(dft_control)) dft_control => qs_env%dft_control
    IF (PRESENT(mos)) mos => qs_env%mos
    IF (PRESENT(neighbor_list_id)) neighbor_list_id=qs_env%neighbor_list_id
    IF (PRESENT(sab_orb)) sab_orb => qs_env%sab_orb
    IF (PRESENT(sac_erfc)) sac_erfc => qs_env%sac_erfc
    IF (PRESENT(sac_ppl)) sac_ppl => qs_env%sac_ppl
    IF (PRESENT(sac_ppnl)) sac_ppnl => qs_env%sac_ppnl
    IF (PRESENT(sbc_ppnl)) sbc_ppnl => qs_env%sbc_ppnl
    IF (PRESENT(sac_oce)) sac_oce =>  qs_env%sac_oce
    IF (PRESENT(sbc_oce)) sbc_oce =>  qs_env%sbc_oce
    IF (PRESENT(sab_2c)) sab_2c =>  qs_env%sab_2c
    IF (PRESENT(sac_3c)) sac_3c =>  qs_env%sac_3c
    IF (PRESENT(sbc_3c)) sbc_3c =>  qs_env%sbc_3c
    IF (PRESENT(energy)) energy => qs_env%energy
    IF (PRESENT(force)) force => qs_env%force
    IF (PRESENT(harris_env)) harris_env => qs_env%harris_env
    IF (PRESENT(matrix_h)) matrix_h => qs_env%matrix_h
    IF (PRESENT(matrix_ks)) matrix_ks => qs_env%matrix_ks
    IF (PRESENT(kinetic)) kinetic => qs_env%kinetic
    IF (PRESENT(matrix_s)) matrix_s => qs_env%matrix_s
    IF (PRESENT(matrix_w)) matrix_w => qs_env%matrix_w
    IF (PRESENT(rho)) rho => qs_env%rho
    IF (PRESENT(rho_xc)) rho_xc => qs_env%rho_xc
    IF (PRESENT(pw_env)) pw_env => qs_env%pw_env
    IF (PRESENT(mpools)) mpools => qs_env%mpools
    IF (PRESENT(blacs_env)) blacs_env => qs_env%blacs_env
    IF (PRESENT(scf_control)) scf_control => qs_env%scf_control
    IF (PRESENT(rho_core)) rho_core => qs_env%rho_core
    IF (PRESENT(qs_charges)) qs_charges => qs_env%qs_charges
    IF (PRESENT(ks_env)) ks_env => qs_env%ks_env
    IF (PRESENT(ks_qmmm_env)) ks_qmmm_env => qs_env%ks_qmmm_env
    IF (PRESENT(wf_history)) wf_history => qs_env%wf_history
    IF (PRESENT(scf_env)) scf_env => qs_env%scf_env
    IF (PRESENT(id_nr)) id_nr=qs_env%id_nr
    IF (PRESENT(subsys)) subsys => qs_env%subsys
    IF (PRESENT(para_env)) para_env => qs_env%para_env
    IF (PRESENT(distribution_2d)) distribution_2d => qs_env%distribution_2d
    IF (PRESENT(oce)) oce =>  qs_env%oce
    IF (PRESENT(requires_mo_derivs)) requires_mo_derivs = qs_env%requires_mo_derivs
    IF (PRESENT(mo_derivs)) mo_derivs => qs_env%mo_derivs
    IF (PRESENT(tddfpt_control)) tddfpt_control => qs_env%dft_control%tddfpt_control
    IF (PRESENT(linres_control)) linres_control => qs_env%linres_control

    IF (PRESENT(rho_atom_set)) &
       CALL get_local_rho(qs_env%local_rho_set, rho_atom_set=rho_atom_set)
    IF (PRESENT(rho0_atom_set)) &
       CALL get_local_rho(qs_env%local_rho_set, rho0_atom_set=rho0_atom_set)
    IF (PRESENT(rho0_mpole)) &
         CALL get_local_rho(qs_env%local_rho_set,  rho0_mpole=rho0_mpole)
    IF (PRESENT(nddo_mpole)) &
         CALL get_local_rho(qs_env%local_rho_set,  nddo_mpole=nddo_mpole)
    IF (PRESENT(rhoz_set)) & 
         CALL get_local_rho(qs_env%local_rho_set, rhoz_set=rhoz_set)
    IF (PRESENT(ecoul_1c)) &
         CALL get_hartree_local(qs_env%hartree_local,ecoul_1c=ecoul_1c)
    IF (PRESENT(VH_3c_list)) &
         CALL get_hartree_local(qs_env%hartree_local,vtriple=VH_3c_list)
    IF (PRESENT(rho0_s_rs)) THEN
       CALL get_local_rho(qs_env%local_rho_set,  rho0_mpole=rho0_m)
       rho0_s_rs => rho0_m%rho0_s_rs
    END IF
    IF (PRESENT(rho0_s_gs)) THEN
       CALL get_local_rho(qs_env%local_rho_set,  rho0_mpole=rho0_m)
       IF(ASSOCIATED(rho0_m)) THEN      
         rho0_s_gs => rho0_m%rho0_s_gs
       END IF
    END IF

    IF(PRESENT(kg_sub_pw_env)) kg_sub_pw_env => qs_env%kg_sub_pw_env
    IF(PRESENT(xas_env)) xas_env => qs_env%xas_env

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
    IF (PRESENT(input)) input => qs_env%input

  END SUBROUTINE get_qs_env

! *****************************************************************************

  SUBROUTINE init_qs_env(qs_env,para_env, error)

!   Purpose: Initialise the QUICKSTEP environment.

!   History: - Creation (25.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

!   ---------------------------------------------------------------------------

    NULLIFY (qs_env%blacs_env)
    NULLIFY (qs_env%green)
    NULLIFY (qs_env%cell)
    NULLIFY (qs_env%super_cell)
    NULLIFY (qs_env%cell_ref)
    NULLIFY (qs_env%dft_control)
    NULLIFY (qs_env%mos)
    NULLIFY (qs_env%sab_orb)
    NULLIFY (qs_env%sac_erfc)
    NULLIFY (qs_env%sac_ppl)
    NULLIFY (qs_env%sac_ppnl)
    NULLIFY (qs_env%sbc_ppnl)
    NULLIFY (qs_env%sac_oce)
    NULLIFY (qs_env%sbc_oce)
    NULLIFY (qs_env%sab_2c)
    NULLIFY (qs_env%sac_3c)
    NULLIFY (qs_env%sbc_3c)
    NULLIFY (qs_env%energy)
    NULLIFY (qs_env%force)
    NULLIFY (qs_env%harris_env)
    NULLIFY (qs_env%matrix_h)
    NULLIFY (qs_env%matrix_ks)
    NULLIFY (qs_env%kinetic)
    NULLIFY (qs_env%matrix_s)
    NULLIFY (qs_env%matrix_w)
    NULLIFY (qs_env%mpools)
    NULLIFY (qs_env%rho)
    NULLIFY (qs_env%rho_xc)
    NULLIFY (qs_env%pw_env)
    NULLIFY (qs_env%scf_control)
    NULLIFY (qs_env%qs_charges)
    NULLIFY (qs_env%rho_core)
    NULLIFY (qs_env%ks_env)
    NULLIFY (qs_env%ks_qmmm_env)
    NULLIFY (qs_env%wf_history)
    NULLIFY (qs_env%scf_env)
    NULLIFY (qs_env%distribution_2d)
    NULLIFY (qs_env%subsys)
    NULLIFY (qs_env%blacs_env)
    NULLIFY (qs_env%oce)
    NULLIFY (qs_env%local_rho_set)
    NULLIFY (qs_env%hartree_local)
    NULLIFY (qs_env%input)
    NULLIFY (qs_env%kg_sub_pw_env)
    NULLIFY (qs_env%linres_control)
    NULLIFY (qs_env%xas_env)
    qs_env%neighbor_list_id=-1

    CALL local_rho_set_create(qs_env%local_rho_set, error=error)
    CALL hartree_local_create(qs_env%hartree_local, error=error)

    CALL cp_subsys_create(qs_env%subsys, para_env=para_env, error=error)
    CALL cp_para_env_retain(para_env,error=error)
    qs_env%para_env => para_env
    CALL cp_blacs_env_create(qs_env%blacs_env,qs_env%para_env,error=error)
    qs_env%ref_count=1
    last_qs_env_id_nr=last_qs_env_id_nr+1
    qs_env%id_nr=last_qs_env_id_nr
    qs_env%use_ref_cell=.FALSE.
    qs_env%qmmm=.FALSE.
    qs_env%requires_mo_derivs=.FALSE.
    NULLIFY(qs_env%mo_derivs)
  END SUBROUTINE init_qs_env

! *****************************************************************************

  SUBROUTINE set_qs_env(qs_env,cell,super_cell,cell_ref,use_ref_cell, &
                        dft_control,mos,qmmm,green,&
                        sab_orb,sac_erfc,sac_ppl,sac_ppnl,sbc_ppnl,&
                        sac_oce,sbc_oce,sab_2c,sac_3c,sbc_3c,&
                        energy,force,harris_env,&
                        matrix_h,matrix_ks,matrix_s,matrix_w,&
                        pw_env,mpools,rho,rho_xc,kinetic,&
                        scf_control,qs_charges,rho_core,ks_env,ks_qmmm_env,&
                        wf_history,scf_env,&
                        distribution_2d,dft_section,input,&
                        subsys, oce,rho_atom_set,rho0_atom_set,rho0_mpole,&
                        nddo_mpole,rhoz_set,rhoz_tot,ecoul_1c,VH_3c_list,error, &
                        requires_mo_derivs,mo_derivs,neighbor_list_id,&
                        kg_sub_pw_env,linres_control,xas_env)

!   Purpose: Set the QUICKSTEP environment.

!   History: - Creation (23.01.2002,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(greens_function_type), POINTER,OPTIONAL  :: green
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, super_cell, cell_ref
    LOGICAL, OPTIONAL                        :: use_ref_cell
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mos
    LOGICAL, OPTIONAL                        :: qmmm
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_orb, sac_erfc, sac_ppl, &
                                                sac_ppnl, sbc_ppnl,&
                                                sac_oce, sbc_oce,&
                                                sab_2c,sac_3c,sbc_3c
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(harris_env_type), POINTER, &
       OPTIONAL                              :: harris_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h, matrix_ks, matrix_s, matrix_w
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env
    TYPE(qs_matrix_pools_type), OPTIONAL, &
      POINTER                                :: mpools
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho_xc
    TYPE(real_matrix_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: kinetic
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                :: scf_control
    TYPE(qs_charges_type), OPTIONAL, POINTER :: qs_charges
    TYPE(coeff_type), OPTIONAL, POINTER      :: rho_core
    TYPE(qs_ks_env_type), OPTIONAL, POINTER  :: ks_env
    TYPE(qs_ks_qmmm_env_type), OPTIONAL, &
      POINTER                                :: ks_qmmm_env
    TYPE(qs_wf_history_type), OPTIONAL, &
      POINTER                                :: wf_history
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(cp_subsystem_type), OPTIONAL, &
      POINTER                                :: subsys
    TYPE(oce_matrix_type), OPTIONAL, POINTER :: oce
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_atom_set
    TYPE(rho0_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho0_atom_set
    TYPE(rho0_mpole_type), OPTIONAL, POINTER :: rho0_mpole
    TYPE(rho0_mpole_type), OPTIONAL, POINTER :: nddo_mpole
    TYPE(rhoz_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rhoz_set
    REAL(dp), OPTIONAL                       :: rhoz_tot
    TYPE(ecoul_1center_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ecoul_1c
    TYPE(vtriple_type), DIMENSION(:,:),&
      OPTIONAL, POINTER                      :: VH_3c_list
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    TYPE(section_vals_type), POINTER, optional :: input
    LOGICAL, OPTIONAL                        :: requires_mo_derivs
    TYPE(cp_fm_p_type), DIMENSION(:),&
       POINTER, OPTIONAL                     :: mo_derivs
    TYPE(section_vals_type), POINTER, optional :: dft_section
    INTEGER, OPTIONAL                        :: neighbor_list_id
    TYPE(kg_sub_pw_env_type), OPTIONAL,&
      POINTER                                :: kg_sub_pw_env
    TYPE(linres_control_type), OPTIONAL,&
      POINTER                                :: linres_control
    TYPE(xas_environment_type), OPTIONAL,&
      POINTER                                :: xas_env

    CHARACTER(len=*), PARAMETER :: routineN = "set_qs_env", &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure=.FALSE.
!    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(qs_env%cell,error=error)
       qs_env%cell => cell
    END IF
    IF (PRESENT(super_cell)) THEN
       CALL cell_retain(super_cell, error=error)
       CALL cell_release(qs_env%super_cell,error=error)
       qs_env%super_cell => super_cell
    END IF
    IF (PRESENT(cell_ref)) THEN
       CALL cell_retain(cell_ref, error=error)
       CALL cell_release(qs_env%cell_ref,error=error)
       qs_env%cell_ref => cell_ref
    END IF
    IF (PRESENT(green)) THEN
       CALL pw_green_retain(green, error=error)
       CALL pw_green_release(qs_env%green,error=error)
       qs_env%green => green
    END IF
    IF (PRESENT(dft_control)) THEN
       CALL dft_control_retain(dft_control, error=error)
       CALL dft_control_release(qs_env%dft_control, error=error)
       qs_env%dft_control => dft_control
    END IF
    IF ( PRESENT ( use_ref_cell ) ) qs_env%use_ref_cell = use_ref_cell
    IF ( PRESENT ( qmmm ) ) qs_env%qmmm = qmmm
    IF (PRESENT(mos)) qs_env%mos => mos
    IF (PRESENT(neighbor_list_id)) qs_env%neighbor_list_id=neighbor_list_id
    IF (PRESENT(sab_orb)) qs_env%sab_orb => sab_orb
    IF (PRESENT(sac_erfc)) qs_env%sac_erfc => sac_erfc
    IF (PRESENT(sac_ppl)) qs_env%sac_ppl => sac_ppl
    IF (PRESENT(sac_ppnl)) qs_env%sac_ppnl => sac_ppnl
    IF (PRESENT(sbc_ppnl)) qs_env%sbc_ppnl => sbc_ppnl
    IF (PRESENT(sac_oce)) qs_env%sac_oce => sac_oce
    IF (PRESENT(sbc_oce)) qs_env%sbc_oce => sbc_oce
    IF (PRESENT(sab_2c)) qs_env%sab_2c => sab_2c
    IF (PRESENT(sac_3c)) qs_env%sac_3c => sac_3c
    IF (PRESENT(sbc_3c)) qs_env%sbc_3c => sbc_3c
    IF (PRESENT(energy)) qs_env%energy => energy
    IF (PRESENT(force)) qs_env%force => force
    IF (PRESENT(harris_env)) qs_env%harris_env => harris_env
    IF (PRESENT(matrix_h)) qs_env%matrix_h => matrix_h
    IF (PRESENT(matrix_ks)) qs_env%matrix_ks => matrix_ks
    IF (PRESENT(matrix_s)) qs_env%matrix_s => matrix_s
    IF (PRESENT(matrix_w)) qs_env%matrix_w => matrix_w
    IF (PRESENT(oce)) qs_env%oce => oce
    IF (PRESENT(requires_mo_derivs)) qs_env%requires_mo_derivs=requires_mo_derivs
    IF (PRESENT(mo_derivs)) qs_env%mo_derivs=>mo_derivs
    IF (PRESENT(rho)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho)) CALL qs_rho_retain(rho,error=error)
      CALL qs_rho_release(qs_env%rho,error=error)
      qs_env%rho => rho
    END IF
    IF (PRESENT(rho_xc)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho_xc)) CALL qs_rho_retain(rho_xc,error=error)
      CALL qs_rho_release(qs_env%rho_xc,error=error)
      qs_env%rho_xc => rho_xc
    END IF

    IF(PRESENT(kg_sub_pw_env)) THEN
      CALL kg_sub_pw_env_retain(kg_sub_pw_env,error=error)
      CALL kg_sub_pw_env_release(qs_env%kg_sub_pw_env, error=error)
      qs_env%kg_sub_pw_env => kg_sub_pw_env
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
    IF (PRESENT(linres_control)) THEN ! accept also null pointers?
       CALL linres_control_retain(linres_control,error=error)
       CALL linres_control_release(qs_env%linres_control,error=error)
       qs_env%linres_control => linres_control
    END IF
    IF (PRESENT(kinetic)) qs_env%kinetic => kinetic
    IF (PRESENT(rho_core)) qs_env%rho_core => rho_core
    IF (PRESENT(qs_charges)) THEN
       CALL qs_charges_retain(qs_charges,error=error)
       CALL qs_charges_release(qs_env%qs_charges,error=error)
       qs_env%qs_charges => qs_charges
    END IF
    IF (PRESENT(ks_qmmm_env)) THEN
       CALL qs_ks_qmmm_retain(ks_qmmm_env, error=error)
       CALL qs_ks_qmmm_release(qs_env%ks_qmmm_env, error=error)
       qs_env%ks_qmmm_env => ks_qmmm_env
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
    IF (PRESENT(xas_env)) THEN ! accept also null pointers?
       CALL xas_env_retain(xas_env, error=error)
       CALL xas_env_release(qs_env%xas_env, error=error)
       qs_env%xas_env => xas_env
    END IF
    IF (PRESENT(mpools)) THEN
       CALL mpools_retain(mpools,error=error)
       CALL mpools_release(qs_env%mpools, error=error) 
       qs_env%mpools => mpools
    END IF
    IF (PRESENT(rho_atom_set)) THEN
        CALL set_local_rho(qs_env%local_rho_set,rho_atom_set=rho_atom_set)
    END IF
    IF (PRESENT(rho0_atom_set)) THEN
        CALL set_local_rho(qs_env%local_rho_set,rho0_atom_set=rho0_atom_set)
    END IF
    IF (PRESENT(rho0_mpole)) THEN
        CALL set_local_rho(qs_env%local_rho_set,rho0_mpole=rho0_mpole)
    END IF
    IF (PRESENT(nddo_mpole)) THEN
        CALL set_local_rho(qs_env%local_rho_set,nddo_mpole=nddo_mpole)
    END IF
    IF (PRESENT(rhoz_set)) THEN
        CALL set_local_rho(qs_env%local_rho_set,rhoz_set=rhoz_set)
    END IF
    IF (PRESENT(rhoz_tot)) qs_env%local_rho_set%rhoz_tot = rhoz_tot
    IF (PRESENT(ecoul_1c)) THEN
        CALL set_hartree_local(qs_env%hartree_local,ecoul_1c=ecoul_1c)
    END IF
    IF (PRESENT(VH_3c_list)) THEN
        CALL set_hartree_local(qs_env%hartree_local,vtriple=VH_3c_list)
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
    IF (PRESENT(input)) THEN
       CALL section_vals_retain(input,error=error)
       CALL section_vals_release(qs_env%input,error=error)
       qs_env%input => input
    END IF

  END SUBROUTINE set_qs_env

! *****************************************************************************

!!****f* qs_environment_types/qs_env_create [1.0] *
!!
!!   NAME
!!     qs_env_create
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
     CALL init_qs_env(qs_env,para_env=para_env,&
          error=error)
  END IF
END SUBROUTINE qs_env_create
!***************************************************************************

!!****f* qs_environment_types/qs_env_retain [1.0] *
!!
!!   NAME
!!     qs_env_retain
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
        CALL cell_release(qs_env%super_cell)
        CALL cell_release(qs_env%cell_ref)
        CALL pw_green_release(qs_env%green)
        CALL dft_control_release(qs_env%dft_control, error=error)
        IF (ASSOCIATED(qs_env%mos)) THEN
           DO i=1,SIZE(qs_env%mos)
              CALL deallocate_mo_set(qs_env%mos(i)%mo_set)
           END DO
           DEALLOCATE(qs_env%mos, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%mo_derivs)) THEN
           DO I=1,SIZE(qs_env%mo_derivs)
              CALL cp_fm_release(qs_env%mo_derivs(I)%matrix)
           ENDDO
           DEALLOCATE(qs_env%mo_derivs, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
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
        IF (ASSOCIATED(qs_env%sac_oce)) THEN
           DO i=1,SIZE(qs_env%sac_oce)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sac_oce(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sac_oce,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sbc_oce)) THEN
           DO i=1,SIZE(qs_env%sbc_oce)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sbc_oce(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sbc_oce,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sab_2c)) THEN
           DO i=1,SIZE(qs_env%sab_2c)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sab_2c(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sab_2c,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sac_3c)) THEN
           DO i=1,SIZE(qs_env%sac_3c)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sac_3c(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sac_3c,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%sbc_3c)) THEN
           DO i=1,SIZE(qs_env%sbc_3c)
              CALL deallocate_neighbor_list_set(&
                   qs_env%sbc_3c(i)%neighbor_list_set)
           END DO
           DEALLOCATE(qs_env%sbc_3c,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        IF (ASSOCIATED(qs_env%energy)) &
             CALL deallocate_qs_energy(qs_env%energy)
        IF (ASSOCIATED(qs_env%force))&
             CALL deallocate_qs_force(qs_env%force)
        IF (ASSOCIATED(qs_env%harris_env)) THEN
            CALL harris_env_release(qs_env%harris_env, error=error)
        END IF
        IF (ASSOCIATED(qs_env%matrix_h)) CALL deallocate_matrix_set(qs_env%matrix_h)

        IF (ASSOCIATED(qs_env%matrix_ks)) THEN
            CALL mpools_get(qs_env%mpools, S_sm_pools=S_sm_pools, error=error)
            CALL sm_pools_give_back_matrix_vect(S_sm_pools,qs_env % matrix_ks,error=error)
        ENDIF 
        IF (ASSOCIATED(qs_env%matrix_s)) CALL deallocate_matrix_set(qs_env%matrix_s)
        IF (ASSOCIATED(qs_env%matrix_w)) CALL deallocate_matrix_set(qs_env%matrix_w)
        IF (ASSOCIATED(qs_env%oce)) CALL deallocate_oce_set(qs_env%oce)

        IF (ASSOCIATED(qs_env%local_rho_set)) THEN 
           CALL local_rho_set_release(qs_env%local_rho_set)
        END IF
        IF (ASSOCIATED(qs_env%hartree_local)) THEN 
           CALL hartree_local_release(qs_env%hartree_local)
        END IF
        IF (ASSOCIATED(qs_env%kinetic))&
             CALL deallocate_matrix_set(qs_env%kinetic)

        IF(ASSOCIATED(qs_env%rho_xc)) THEN
          NULLIFY(qs_env%rho_xc%rho_ao)
          CALL qs_rho_release(qs_env%rho_xc, error=error)
        END IF
        CALL qs_rho_release(qs_env%rho, error=error)
        CALL scf_c_release(qs_env%scf_control, error=error)

        IF(ASSOCIATED(qs_env%linres_control)) THEN
          CALL linres_control_release(qs_env%linres_control, error=error)
        END IF

        IF(ASSOCIATED(qs_env%kg_sub_pw_env)) THEN
          CALL kg_sub_pw_env_release(qs_env%kg_sub_pw_env,error=error)
        END IF

        !Only if do_xas_calculation
        IF(ASSOCIATED(qs_env%xas_env)) THEN
          CALL xas_env_release(qs_env%xas_env,error=error)
        END IF

        CALL pw_env_release(qs_env%pw_env, error=error)
        CALL cp_para_env_release(qs_env%para_env,error=error)
        CALL cp_blacs_env_release(qs_env%blacs_env, error=error)
        IF (ASSOCIATED(qs_env%rho_core)) THEN
             CALL coeff_deallocate(qs_env%rho_core)
             DEALLOCATE(qs_env%rho_core)
        ENDIF
        CALL qs_charges_release(qs_env%qs_charges, error=error)
        CALL qs_ks_release(qs_env%ks_env, error=error)
        CALL qs_ks_qmmm_release(qs_env%ks_qmmm_env, error=error)
        CALL wfi_release(qs_env%wf_history,error=error)
        CALL scf_env_release(qs_env%scf_env, error=error)
        CALL mpools_release(qs_env%mpools,error=error)
        IF (ASSOCIATED(qs_env%distribution_2d)) THEN
           CALL distribution_2d_release(qs_env%distribution_2d)
        ENDIF
        CALL cp_subsys_release(qs_env%subsys,error=error)
        CALL section_vals_release(qs_env%input,error=error)
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
      POINTER                                :: matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control

  failure=.FALSE.
  NULLIFY(scf_env,scf_control,sv,maxao_maxmo_fm_pool,matrix_s,mpools)
  CALL timeset(routineN,'I',' ',handle)
  
  CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(v_matrix),cp_failure_level,routineP,error,failure)
  
  IF (.NOT. failure) THEN

     CALL cp_fm_get_info(v_matrix,ncol_global=my_n_col,error=error)
     IF (PRESENT(n_col)) my_n_col=n_col
     CALL get_qs_env(qs_env,mpools=mpools,&
          scf_env=scf_env,&
          scf_control=scf_control,&
          matrix_s=matrix_s,&
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
        CALL cp_sm_fm_multiply(sparse_matrix=matrix_s(1)%matrix,v_in=v_matrix,&
             v_out=sv,ncol=my_n_col, error=error)
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
