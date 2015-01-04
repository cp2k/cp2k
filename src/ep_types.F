!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief the types needed for the calculation of the enrgy with perturbation
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE ep_types
  USE cp_array_r_utils,                ONLY: cp_2d_r_p_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE iso_c_binding
  USE kinds,                           ONLY: dp
  USE preconditioner_types,            ONLY: preconditioner_p_type
  USE qs_environment_types,            ONLY: qs_environment_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_p_sparse_psi,                 ONLY: qs_p_projection_p_type
  USE replica_types,                   ONLY: replica_env_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ep_types'

  PUBLIC :: ep_env_type, ep_env_p_type, ep_energy_type, ep_force_type


! *****************************************************************************
!> \brief stores various energies calculated during an ep calculation
!> \param e_no_int energy of the single fragments (no interaction)
!> \param e 0: the energy of the orthogonalized frozen densities
!>        (energy of psi0, E(psi0))
!> \param e 1: the correction energy (E(psi0+psi1)-E(psi0))
!> \param e_tot total energy
!> \author fawzi
! *****************************************************************************
  TYPE ep_energy_type
     REAL(dp) :: e_no_int,e0,e1,e_tot
  END TYPE ep_energy_type

! *****************************************************************************
!> \brief stores various forces calculated during an ep calculation
!>      (only those not stored in main_qs_env)
!> \param f 0_internal: intermolecular forces due to e0
!> \author fawzi
! *****************************************************************************
  TYPE ep_force_type
     INTEGER :: id_nr, ref_count
     REAL(dp), DIMENSION(:,:), POINTER :: f0_internal
  END TYPE ep_force_type

! *****************************************************************************
!> \brief environement for an ep calculation
!> \param id_nr identification number (unique)
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param mol_envs qs_envs that represent the pieces of the full system
!>        handled by a replica_env
!> \param sub_proj projection from the full system to the pieces
!> \param main_qs_env qs_env of the full system
!> \param main_p_env perturbation environment for the full system
!> \param sub_p_env perturbation environement for the molecular subsystem
!> \param m_pi_Hrho_psi 0d: the rhs of the linear system to solve for
!>        the perturbation ( - pi Hrho psi0d )
!> \param nspins number of spins
!> \param nmol number of molecules replicas in mol_envs
!> \param nat number of atoms (in main_qs_env)
!> \param nat_per_mol number of atoms in each mol
!> \param sub_nao number of atom orbitals in the mol_env
!> \param sub_nmo number of molecular orbitals in the mol_env
!> \param full_nao number of atom orbitals in main_qs_env
!> \param full_nmo number of molecular orbitals in main_qs_env
!> \param at 2sub: which molecule owns each atom
!> \param para_env the global parallel environment for the ep calculation
!> \param input the input for the ep calculation
!> \param globenv global environment for the main_qs
!> \param energy the energy of the last configuration evaluated
!> \param local_H_full diagonal elements of the KohnSham matrix
!> \param local_matrix_w matrix_w for the intra forces
!> \par History
!>      11.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE ep_env_type
     INTEGER :: id_nr, ref_count, f_env_id,nspins,nat,nmol,nat_per_mol
     TYPE(global_environment_type), POINTER :: globenv
     TYPE(section_vals_type), POINTER         :: root_section
     TYPE(cp_para_env_type), POINTER :: para_env
     TYPE(section_vals_type), POINTER :: input
     TYPE(replica_env_type), POINTER :: mol_envs
     TYPE(qs_p_projection_p_type), &
          DIMENSION(:), POINTER        :: sub_proj
     TYPE(qs_environment_type), POINTER :: main_qs_env
     TYPE(qs_p_env_type), POINTER :: main_p_env, sub_p_env
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: m_pi_Hrho_psi0d
     TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: psi1
     TYPE(preconditioner_p_type), DIMENSION(:), POINTER :: precond
     INTEGER, DIMENSION(:), POINTER :: sub_nmo,sub_nao,full_nmo,full_nao,at2sub
     TYPE(ep_energy_type) :: energy
     TYPE(ep_force_type), POINTER :: force
     TYPE(cp_2d_r_p_type), DIMENSION(:), POINTER :: base_C0
  END TYPE ep_env_type

! *****************************************************************************
!> \brief type to have array of pointers
!> \param ep_env the pointer to the ep env
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE ep_env_p_type
     TYPE(ep_env_type), POINTER :: ep_env
  END TYPE ep_env_p_type


END MODULE ep_types
