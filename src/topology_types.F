!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH (30.11.2001) : new entries in setup_parameters_type
!>                         change name from input_file_name to coord_...
!>                         added topology file
!>                         added atom_names
!>      Teodoro Laino [tlaino] 12.2008 - Preparing for VIRTUAL SITE constraints
!>                                       (patch by Marcel Baer)
!> \author CJM & JGH
! **************************************************************************************************
MODULE topology_types
   USE cell_types,                      ONLY: cell_release,&
                                              cell_type
   USE colvar_types,                    ONLY: colvar_p_type,&
                                              colvar_release
   USE input_constants,                 ONLY: do_bondparm_covalent,&
                                              do_conn_generate,&
                                              do_constr_none,&
                                              do_skip_13
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

! **************************************************************************************************
   TYPE atom_info_type
      INTEGER, DIMENSION(:), POINTER         :: id_molname => NULL()
      INTEGER, DIMENSION(:), POINTER         :: id_resname => NULL()
      INTEGER, DIMENSION(:), POINTER         :: id_atmname => NULL()
      INTEGER, DIMENSION(:), POINTER         :: id_atom_names => NULL()
      INTEGER, DIMENSION(:), POINTER         :: id_element => NULL()
      INTEGER, POINTER                       :: resid(:) => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER :: r => NULL()
      INTEGER, POINTER                       :: map_mol_typ(:) => NULL()
      INTEGER, POINTER                       :: map_mol_num(:) => NULL()
      INTEGER, POINTER                       :: map_mol_res(:) => NULL()
      REAL(KIND=dp), POINTER                 :: atm_charge(:) => NULL()
      REAL(KIND=dp), POINTER                 :: atm_mass(:) => NULL()
      REAL(KIND=dp), POINTER                 :: occup(:) => NULL()
      REAL(KIND=dp), POINTER                 :: beta(:) => NULL()
   END TYPE atom_info_type

! **************************************************************************************************
   TYPE connectivity_info_type
      INTEGER, POINTER :: bond_a(:) => NULL(), bond_b(:) => NULL(), bond_type(:) => NULL()
      INTEGER, POINTER :: ub_a(:) => NULL(), ub_b(:) => NULL(), ub_c(:) => NULL()
      INTEGER, POINTER :: theta_a(:) => NULL(), theta_b(:) => NULL(), theta_c(:) => NULL(), theta_type(:) => NULL()
      INTEGER, POINTER :: phi_a(:) => NULL(), phi_b(:) => NULL(), phi_c(:) => NULL(), phi_d(:) => NULL(), phi_type(:) => NULL()
      INTEGER, POINTER :: impr_a(:) => NULL(), impr_b(:) => NULL(), impr_c(:) => NULL(), &
                          impr_d(:) => NULL(), impr_type(:) => NULL()
      INTEGER, POINTER :: onfo_a(:) => NULL(), onfo_b(:) => NULL()
      INTEGER, POINTER :: c_bond_a(:) => NULL(), c_bond_b(:) => NULL(), c_bond_type(:) => NULL()
   END TYPE connectivity_info_type

! **************************************************************************************************
   TYPE constraint_info_type
      ! Bonds involving Hydrogens
      LOGICAL                                       :: hbonds_restraint = .FALSE. ! Restraints control
      REAL(KIND=dp)                                 :: hbonds_k0 = -1.0_dp ! Restraints control
      ! Fixed Atoms
      INTEGER                                       :: nfixed_atoms = -1
     INTEGER, POINTER                              :: fixed_atoms(:) => NULL(), fixed_type(:) => NULL(), fixed_mol_type(:) => NULL()
      LOGICAL, POINTER                              :: fixed_restraint(:) => NULL() ! Restraints control
      REAL(KIND=dp), POINTER                        :: fixed_k0(:) => NULL() ! Restraints control
      ! Freeze QM or MM
      INTEGER                                       :: freeze_qm = -1, freeze_mm = -1, freeze_qm_type = -1, freeze_mm_type = -1
      LOGICAL                                       :: fixed_mm_restraint = .FALSE., fixed_qm_restraint = .FALSE. ! Restraints control
      REAL(KIND=dp)                                 :: fixed_mm_k0 = -1.0_dp, fixed_qm_k0 = -1.0_dp ! Restraints control
      ! Freeze with molnames
      LOGICAL, POINTER                              :: fixed_mol_restraint(:) => NULL() ! Restraints control
      REAL(KIND=dp), POINTER                        :: fixed_mol_k0(:) => NULL() ! Restraints control
      CHARACTER(LEN=default_string_length), POINTER :: fixed_molnames(:) => NULL()
      LOGICAL, POINTER, DIMENSION(:)                :: fixed_exclude_qm => NULL(), fixed_exclude_mm => NULL()
      ! Collective constraints
      INTEGER                                       :: nconst_colv = -1
      INTEGER, POINTER                              :: const_colv_mol(:) => NULL()
      CHARACTER(LEN=default_string_length), POINTER :: const_colv_molname(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_colv_target(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_colv_target_growth(:) => NULL()
      TYPE(colvar_p_type), POINTER, DIMENSION(:)    :: colvar_set => NULL()
      LOGICAL, POINTER                              :: colv_intermolecular(:) => NULL()
      LOGICAL, POINTER                              :: colv_restraint(:) => NULL() ! Restraints control
      REAL(KIND=dp), POINTER                        :: colv_k0(:) => NULL() ! Restraints control
      LOGICAL, POINTER, DIMENSION(:)                :: colv_exclude_qm => NULL(), colv_exclude_mm => NULL()
      ! G3x3
      INTEGER                                       :: nconst_g33 = -1
      INTEGER, POINTER                              :: const_g33_mol(:) => NULL()
      CHARACTER(LEN=default_string_length), POINTER :: const_g33_molname(:) => NULL()
      INTEGER, POINTER                              :: const_g33_a(:) => NULL()
      INTEGER, POINTER                              :: const_g33_b(:) => NULL()
      INTEGER, POINTER                              :: const_g33_c(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g33_dab(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g33_dac(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g33_dbc(:) => NULL()
      LOGICAL, POINTER                              :: g33_intermolecular(:) => NULL()
      LOGICAL, POINTER                              :: g33_restraint(:) => NULL() ! Restraints control
      REAL(KIND=dp), POINTER                        :: g33_k0(:) => NULL() ! Restraints control
      LOGICAL, POINTER, DIMENSION(:)                :: g33_exclude_qm => NULL(), g33_exclude_mm => NULL()
      ! G4x6
      INTEGER                                       :: nconst_g46 = -1
      INTEGER, POINTER                              :: const_g46_mol(:) => NULL()
      CHARACTER(LEN=default_string_length), POINTER :: const_g46_molname(:) => NULL()
      INTEGER, POINTER                              :: const_g46_a(:) => NULL()
      INTEGER, POINTER                              :: const_g46_b(:) => NULL()
      INTEGER, POINTER                              :: const_g46_c(:) => NULL()
      INTEGER, POINTER                              :: const_g46_d(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g46_dab(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g46_dac(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g46_dbc(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g46_dad(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g46_dbd(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_g46_dcd(:) => NULL()
      LOGICAL, POINTER                              :: g46_intermolecular(:) => NULL()
      LOGICAL, POINTER                              :: g46_restraint(:) => NULL() ! Restraints control
      REAL(KIND=dp), POINTER                        :: g46_k0(:) => NULL() ! Restraints control
      LOGICAL, POINTER, DIMENSION(:)                :: g46_exclude_qm => NULL(), g46_exclude_mm => NULL()
      ! virtual_site
      INTEGER                                       :: nconst_vsite = -1
      INTEGER, POINTER                              :: const_vsite_mol(:) => NULL()
      CHARACTER(LEN=default_string_length), POINTER :: const_vsite_molname(:) => NULL()
      INTEGER, POINTER                              :: const_vsite_a(:) => NULL()
      INTEGER, POINTER                              :: const_vsite_b(:) => NULL()
      INTEGER, POINTER                              :: const_vsite_c(:) => NULL()
      INTEGER, POINTER                              :: const_vsite_d(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_vsite_wbc(:) => NULL()
      REAL(KIND=dp), POINTER                        :: const_vsite_wdc(:) => NULL()
      LOGICAL, POINTER                              :: vsite_intermolecular(:) => NULL()
      LOGICAL, POINTER                              :: vsite_restraint(:) => NULL() ! Restraints control
      REAL(KIND=dp), POINTER                        :: vsite_k0(:) => NULL() ! Restraints control
      LOGICAL, POINTER, DIMENSION(:)                :: vsite_exclude_qm => NULL(), vsite_exclude_mm => NULL()
   END TYPE constraint_info_type

! **************************************************************************************************
   TYPE topology_parameters_type
      TYPE(atom_info_type), POINTER             :: atom_info => NULL()
      TYPE(connectivity_info_type), POINTER     :: conn_info => NULL()
      TYPE(constraint_info_type), POINTER       :: cons_info => NULL()
      TYPE(cell_type), POINTER                 :: cell => NULL(), cell_ref => NULL(), cell_muc => NULL()
      INTEGER                                   :: conn_type = -1
      INTEGER                                   :: coord_type = -1
      INTEGER                                   :: exclude_vdw = -1
      INTEGER                                   :: exclude_ei = -1
      INTEGER                                   :: bondparm_type = -1
      !TRY TO REMOVE THIS FIVE VARIABLE IN THE FUTURE
      INTEGER                                   :: natoms = -1, natom_type = -1
      INTEGER                                   :: nmol = -1, nmol_type = -1, nmol_conn = -1
      !TRY TO REMOVE THIS FIVE VARIABLE IN THE FUTURE
      LOGICAL                                   :: aa_element = .FALSE.
      LOGICAL                                   :: molname_generated = .FALSE.
      REAL(KIND=dp)                             :: bondparm_factor = -1.0_dp
      LOGICAL                                   :: create_molecules = .FALSE.
      LOGICAL                                   :: reorder_atom = .FALSE.
      LOGICAL                                   :: molecules_check = .FALSE.
      LOGICAL                                   :: coordinate = .FALSE.
      LOGICAL                                   :: use_g96_velocity = .FALSE.
      CHARACTER(LEN=default_path_length)       :: coord_file_name = ""
      CHARACTER(LEN=default_path_length)       :: conn_file_name = ""
      LOGICAL                                   :: const_atom = .FALSE.
      LOGICAL                                   :: const_hydr = .FALSE.
      LOGICAL                                   :: const_colv = .FALSE.
      LOGICAL                                   :: const_33 = .FALSE.
      LOGICAL                                   :: const_46 = .FALSE.
      LOGICAL                                   :: const_vsite = .FALSE.
      LOGICAL                                   :: charge_occup = .FALSE.
      LOGICAL                                   :: charge_beta = .FALSE.
      LOGICAL                                   :: charge_extended = .FALSE.
      LOGICAL                                   :: para_res = .FALSE.
   END TYPE topology_parameters_type

! **************************************************************************************************
   TYPE constr_list_type
      INTEGER, DIMENSION(:), POINTER :: constr => NULL()
   END TYPE constr_list_type

   PUBLIC :: atom_info_type, &
             connectivity_info_type, &
             constraint_info_type, &
             topology_parameters_type, &
             constr_list_type

   PUBLIC :: init_topology, &
             deallocate_topology, &
             pre_read_topology

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_types'
   PRIVATE

CONTAINS

! **************************************************************************************************
!> \brief 1. Just NULLIFY and zero all the stuff
!> \param topology ...
!> \par History
!>      none
! **************************************************************************************************
   SUBROUTINE init_topology(topology)
      TYPE(topology_parameters_type), INTENT(INOUT)      :: topology

!-----------------------------------------------------------------------------
! 1. Nullify and allocate things in topology
!-----------------------------------------------------------------------------

      ALLOCATE (topology%atom_info)
      ALLOCATE (topology%conn_info)
      ALLOCATE (topology%cons_info)
      !-----------------------------------------------------------------------------
      ! 2. Initialize and Nullify things in topology
      !-----------------------------------------------------------------------------
      NULLIFY (topology%cell, topology%cell_ref, topology%cell_muc)
      topology%natoms = 0
      topology%natom_type = 0
      topology%nmol = 0
      topology%nmol_type = 0
      topology%nmol_conn = 0
      topology%bondparm_type = do_bondparm_covalent
      topology%reorder_atom = .FALSE.
      topology%create_molecules = .FALSE.
      topology%molecules_check = .FALSE.
      topology%coordinate = .FALSE.
      topology%use_g96_velocity = .FALSE.
      topology%coord_type = -1
      topology%coord_file_name = ''
      topology%conn_type = do_conn_generate
      topology%conn_file_name = 'OFF'
      topology%const_atom = .FALSE.
      topology%const_hydr = .FALSE.
      topology%const_colv = .FALSE.
      topology%const_33 = .FALSE.
      topology%const_46 = .FALSE.
      topology%const_vsite = .FALSE.
      topology%charge_occup = .FALSE.
      topology%charge_beta = .FALSE.
      topology%charge_extended = .FALSE.
      topology%para_res = .FALSE.
      topology%molname_generated = .FALSE.
      topology%aa_element = .FALSE.
      topology%exclude_vdw = do_skip_13
      topology%exclude_ei = do_skip_13
      !-----------------------------------------------------------------------------
      ! 3. Initialize and Nullify things in topology%atom_info
      !-----------------------------------------------------------------------------
      NULLIFY (topology%atom_info%id_molname)
      NULLIFY (topology%atom_info%id_resname)
      NULLIFY (topology%atom_info%resid)
      NULLIFY (topology%atom_info%id_atmname)
      NULLIFY (topology%atom_info%id_atom_names)
      NULLIFY (topology%atom_info%r)
      NULLIFY (topology%atom_info%map_mol_typ)
      NULLIFY (topology%atom_info%map_mol_num)
      NULLIFY (topology%atom_info%map_mol_res)
      NULLIFY (topology%atom_info%atm_charge)
      NULLIFY (topology%atom_info%atm_mass)
      NULLIFY (topology%atom_info%occup)
      NULLIFY (topology%atom_info%beta)
      NULLIFY (topology%atom_info%id_element)
      !-----------------------------------------------------------------------------
      ! 4. Initialize and Nullify things in topology%conn_info
      !-----------------------------------------------------------------------------
      NULLIFY (topology%conn_info%bond_a)
      NULLIFY (topology%conn_info%bond_b)
      NULLIFY (topology%conn_info%bond_type)
      NULLIFY (topology%conn_info%ub_a)
      NULLIFY (topology%conn_info%ub_b)
      NULLIFY (topology%conn_info%ub_c)
      NULLIFY (topology%conn_info%theta_a)
      NULLIFY (topology%conn_info%theta_b)
      NULLIFY (topology%conn_info%theta_c)
      NULLIFY (topology%conn_info%theta_type)
      NULLIFY (topology%conn_info%phi_a)
      NULLIFY (topology%conn_info%phi_b)
      NULLIFY (topology%conn_info%phi_c)
      NULLIFY (topology%conn_info%phi_d)
      NULLIFY (topology%conn_info%phi_type)
      NULLIFY (topology%conn_info%impr_a)
      NULLIFY (topology%conn_info%impr_b)
      NULLIFY (topology%conn_info%impr_c)
      NULLIFY (topology%conn_info%impr_d)
      NULLIFY (topology%conn_info%impr_type)
      NULLIFY (topology%conn_info%onfo_a)
      NULLIFY (topology%conn_info%onfo_b)
      NULLIFY (topology%conn_info%c_bond_a)
      NULLIFY (topology%conn_info%c_bond_b)
      NULLIFY (topology%conn_info%c_bond_type)
      !-----------------------------------------------------------------------------
      ! 5. Initialize and Nullify things in topology%cons_info
      !-----------------------------------------------------------------------------
      CALL init_constraint(topology%cons_info)
   END SUBROUTINE init_topology

! **************************************************************************************************
!> \brief 1. Just NULLIFY and zero all the stuff
!> \param constraint_info ...
!> \par History
!>      none
! **************************************************************************************************
   SUBROUTINE init_constraint(constraint_info)
      TYPE(constraint_info_type), INTENT(INOUT)          :: constraint_info

! Bonds involving Hydrogens

      constraint_info%hbonds_restraint = .FALSE.
      ! Fixed Atoms
      constraint_info%nfixed_atoms = 0
      constraint_info%freeze_mm = do_constr_none
      constraint_info%freeze_qm = do_constr_none
      NULLIFY (constraint_info%fixed_atoms)
      NULLIFY (constraint_info%fixed_type)
      NULLIFY (constraint_info%fixed_mol_type)
      NULLIFY (constraint_info%fixed_molnames)
      NULLIFY (constraint_info%fixed_restraint)
      NULLIFY (constraint_info%fixed_k0)
      NULLIFY (constraint_info%fixed_mol_restraint)
      NULLIFY (constraint_info%fixed_mol_k0)
      NULLIFY (constraint_info%fixed_exclude_qm, constraint_info%fixed_exclude_mm)
      ! Collective Constraints
      constraint_info%nconst_colv = 0
      NULLIFY (constraint_info%colvar_set)
      NULLIFY (constraint_info%const_colv_mol)
      NULLIFY (constraint_info%const_colv_molname)
      NULLIFY (constraint_info%const_colv_target)
      NULLIFY (constraint_info%const_colv_target_growth)
      NULLIFY (constraint_info%colv_intermolecular)
      NULLIFY (constraint_info%colv_restraint)
      NULLIFY (constraint_info%colv_k0)
      NULLIFY (constraint_info%colv_exclude_qm, constraint_info%colv_exclude_mm)
      ! G3x3
      constraint_info%nconst_g33 = 0
      NULLIFY (constraint_info%const_g33_mol)
      NULLIFY (constraint_info%const_g33_molname)
      NULLIFY (constraint_info%const_g33_a)
      NULLIFY (constraint_info%const_g33_b)
      NULLIFY (constraint_info%const_g33_c)
      NULLIFY (constraint_info%const_g33_dab)
      NULLIFY (constraint_info%const_g33_dac)
      NULLIFY (constraint_info%const_g33_dbc)
      NULLIFY (constraint_info%g33_intermolecular)
      NULLIFY (constraint_info%g33_restraint)
      NULLIFY (constraint_info%g33_k0)
      NULLIFY (constraint_info%g33_exclude_qm, constraint_info%g33_exclude_mm)
      ! G4x6
      constraint_info%nconst_g46 = 0
      NULLIFY (constraint_info%const_g46_mol)
      NULLIFY (constraint_info%const_g46_molname)
      NULLIFY (constraint_info%const_g46_a)
      NULLIFY (constraint_info%const_g46_b)
      NULLIFY (constraint_info%const_g46_c)
      NULLIFY (constraint_info%const_g46_d)
      NULLIFY (constraint_info%const_g46_dab)
      NULLIFY (constraint_info%const_g46_dac)
      NULLIFY (constraint_info%const_g46_dbc)
      NULLIFY (constraint_info%const_g46_dad)
      NULLIFY (constraint_info%const_g46_dbd)
      NULLIFY (constraint_info%const_g46_dcd)
      NULLIFY (constraint_info%g46_intermolecular)
      NULLIFY (constraint_info%g46_restraint)
      NULLIFY (constraint_info%g46_k0)
      NULLIFY (constraint_info%g46_exclude_qm, constraint_info%g46_exclude_mm)
      ! virtual_site
      constraint_info%nconst_vsite = 0
      NULLIFY (constraint_info%const_vsite_mol)
      NULLIFY (constraint_info%const_vsite_molname)
      NULLIFY (constraint_info%const_vsite_a)
      NULLIFY (constraint_info%const_vsite_b)
      NULLIFY (constraint_info%const_vsite_c)
      NULLIFY (constraint_info%const_vsite_d)
      NULLIFY (constraint_info%const_vsite_wbc)
      NULLIFY (constraint_info%const_vsite_wdc)
      NULLIFY (constraint_info%vsite_intermolecular)
      NULLIFY (constraint_info%vsite_restraint)
      NULLIFY (constraint_info%vsite_k0)
      NULLIFY (constraint_info%vsite_exclude_qm, constraint_info%vsite_exclude_mm)

   END SUBROUTINE init_constraint

! **************************************************************************************************
!> \brief 1. Just DEALLOCATE all the stuff
!> \param topology ...
!> \par History
!>      none
! **************************************************************************************************
   SUBROUTINE deallocate_topology(topology)
      TYPE(topology_parameters_type), INTENT(INOUT)      :: topology

!-----------------------------------------------------------------------------
! 1. DEALLOCATE things in topology%atom_info
!-----------------------------------------------------------------------------

      IF (ASSOCIATED(topology%atom_info%id_molname)) THEN
         DEALLOCATE (topology%atom_info%id_molname)
      END IF
      IF (ASSOCIATED(topology%atom_info%id_resname)) THEN
         DEALLOCATE (topology%atom_info%id_resname)
      END IF
      IF (ASSOCIATED(topology%atom_info%resid)) THEN
         DEALLOCATE (topology%atom_info%resid)
      END IF
      IF (ASSOCIATED(topology%atom_info%id_atmname)) THEN
         DEALLOCATE (topology%atom_info%id_atmname)
      END IF
      IF (ASSOCIATED(topology%atom_info%id_atom_names)) THEN
         DEALLOCATE (topology%atom_info%id_atom_names)
      END IF
      IF (ASSOCIATED(topology%atom_info%r)) THEN
         DEALLOCATE (topology%atom_info%r)
      END IF
      IF (ASSOCIATED(topology%atom_info%map_mol_typ)) THEN
         DEALLOCATE (topology%atom_info%map_mol_typ)
      END IF
      IF (ASSOCIATED(topology%atom_info%map_mol_num)) THEN
         DEALLOCATE (topology%atom_info%map_mol_num)
      END IF
      IF (ASSOCIATED(topology%atom_info%map_mol_res)) THEN
         DEALLOCATE (topology%atom_info%map_mol_res)
      END IF
      IF (ASSOCIATED(topology%atom_info%atm_charge)) THEN
         DEALLOCATE (topology%atom_info%atm_charge)
      END IF
      IF (ASSOCIATED(topology%atom_info%atm_mass)) THEN
         DEALLOCATE (topology%atom_info%atm_mass)
      END IF
      IF (ASSOCIATED(topology%atom_info%occup)) THEN
         DEALLOCATE (topology%atom_info%occup)
      END IF
      IF (ASSOCIATED(topology%atom_info%beta)) THEN
         DEALLOCATE (topology%atom_info%beta)
      END IF
      IF (ASSOCIATED(topology%atom_info%id_element)) THEN
         DEALLOCATE (topology%atom_info%id_element)
      END IF
      !-----------------------------------------------------------------------------
      ! 2. DEALLOCATE things in topology%conn_info
      !-----------------------------------------------------------------------------
      IF (ASSOCIATED(topology%conn_info%bond_a)) THEN
         DEALLOCATE (topology%conn_info%bond_a)
      END IF
      IF (ASSOCIATED(topology%conn_info%bond_b)) THEN
         DEALLOCATE (topology%conn_info%bond_b)
      END IF
      IF (ASSOCIATED(topology%conn_info%bond_type)) THEN
         DEALLOCATE (topology%conn_info%bond_type)
      END IF
      IF (ASSOCIATED(topology%conn_info%ub_a)) THEN
         DEALLOCATE (topology%conn_info%ub_a)
      END IF
      IF (ASSOCIATED(topology%conn_info%ub_b)) THEN
         DEALLOCATE (topology%conn_info%ub_b)
      END IF
      IF (ASSOCIATED(topology%conn_info%ub_c)) THEN
         DEALLOCATE (topology%conn_info%ub_c)
      END IF
      IF (ASSOCIATED(topology%conn_info%theta_a)) THEN
         DEALLOCATE (topology%conn_info%theta_a)
      END IF
      IF (ASSOCIATED(topology%conn_info%theta_b)) THEN
         DEALLOCATE (topology%conn_info%theta_b)
      END IF
      IF (ASSOCIATED(topology%conn_info%theta_c)) THEN
         DEALLOCATE (topology%conn_info%theta_c)
      END IF
      IF (ASSOCIATED(topology%conn_info%theta_type)) THEN
         DEALLOCATE (topology%conn_info%theta_type)
      END IF
      IF (ASSOCIATED(topology%conn_info%phi_a)) THEN
         DEALLOCATE (topology%conn_info%phi_a)
      END IF
      IF (ASSOCIATED(topology%conn_info%phi_b)) THEN
         DEALLOCATE (topology%conn_info%phi_b)
      END IF
      IF (ASSOCIATED(topology%conn_info%phi_c)) THEN
         DEALLOCATE (topology%conn_info%phi_c)
      END IF
      IF (ASSOCIATED(topology%conn_info%phi_d)) THEN
         DEALLOCATE (topology%conn_info%phi_d)
      END IF
      IF (ASSOCIATED(topology%conn_info%phi_type)) THEN
         DEALLOCATE (topology%conn_info%phi_type)
      END IF
      IF (ASSOCIATED(topology%conn_info%impr_a)) THEN
         DEALLOCATE (topology%conn_info%impr_a)
      END IF
      IF (ASSOCIATED(topology%conn_info%impr_b)) THEN
         DEALLOCATE (topology%conn_info%impr_b)
      END IF
      IF (ASSOCIATED(topology%conn_info%impr_c)) THEN
         DEALLOCATE (topology%conn_info%impr_c)
      END IF
      IF (ASSOCIATED(topology%conn_info%impr_d)) THEN
         DEALLOCATE (topology%conn_info%impr_d)
      END IF
      IF (ASSOCIATED(topology%conn_info%impr_type)) THEN
         DEALLOCATE (topology%conn_info%impr_type)
      END IF
      IF (ASSOCIATED(topology%conn_info%onfo_a)) THEN
         DEALLOCATE (topology%conn_info%onfo_a)
      END IF
      IF (ASSOCIATED(topology%conn_info%onfo_b)) THEN
         DEALLOCATE (topology%conn_info%onfo_b)
      END IF
      IF (ASSOCIATED(topology%conn_info%c_bond_a)) THEN
         DEALLOCATE (topology%conn_info%c_bond_a)
      END IF
      IF (ASSOCIATED(topology%conn_info%c_bond_b)) THEN
         DEALLOCATE (topology%conn_info%c_bond_b)
      END IF
      IF (ASSOCIATED(topology%conn_info%c_bond_type)) THEN
         DEALLOCATE (topology%conn_info%c_bond_type)
      END IF
      !-----------------------------------------------------------------------------
      ! 3. DEALLOCATE things in topology%cons_info
      !-----------------------------------------------------------------------------
      IF (ASSOCIATED(topology%cons_info)) &
         CALL deallocate_constraint(topology%cons_info)
      !-----------------------------------------------------------------------------
      ! 4. DEALLOCATE things in topology
      !-----------------------------------------------------------------------------
      CALL cell_release(topology%cell)
      CALL cell_release(topology%cell_ref)
      CALL cell_release(topology%cell_muc)
      IF (ASSOCIATED(topology%atom_info)) THEN
         DEALLOCATE (topology%atom_info)
      END IF
      IF (ASSOCIATED(topology%conn_info)) THEN
         DEALLOCATE (topology%conn_info)
      END IF
      IF (ASSOCIATED(topology%cons_info)) THEN
         DEALLOCATE (topology%cons_info)
      END IF

   END SUBROUTINE deallocate_topology

! **************************************************************************************************
!> \brief 1. Just DEALLOCATE all the stuff
!> \param constraint_info ...
!> \par History
!>      none
! **************************************************************************************************
   SUBROUTINE deallocate_constraint(constraint_info)
      TYPE(constraint_info_type), INTENT(INOUT)          :: constraint_info

      INTEGER                                            :: i

! Fixed Atoms

      IF (ASSOCIATED(constraint_info%fixed_atoms)) THEN
         DEALLOCATE (constraint_info%fixed_atoms)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_type)) THEN
         DEALLOCATE (constraint_info%fixed_type)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_molnames)) THEN
         DEALLOCATE (constraint_info%fixed_molnames)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_mol_type)) THEN
         DEALLOCATE (constraint_info%fixed_mol_type)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_restraint)) THEN
         DEALLOCATE (constraint_info%fixed_restraint)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_k0)) THEN
         DEALLOCATE (constraint_info%fixed_k0)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_mol_restraint)) THEN
         DEALLOCATE (constraint_info%fixed_mol_restraint)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_mol_k0)) THEN
         DEALLOCATE (constraint_info%fixed_mol_k0)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_exclude_qm)) THEN
         DEALLOCATE (constraint_info%fixed_exclude_qm)
      END IF
      IF (ASSOCIATED(constraint_info%fixed_exclude_mm)) THEN
         DEALLOCATE (constraint_info%fixed_exclude_mm)
      END IF
      ! Collective Constraint
      IF (ASSOCIATED(constraint_info%colvar_set)) THEN
         DO i = 1, SIZE(constraint_info%colvar_set)
            IF (ASSOCIATED(constraint_info%colvar_set(i)%colvar)) THEN
               CALL colvar_release(constraint_info%colvar_set(i)%colvar)
               NULLIFY (constraint_info%colvar_set(i)%colvar)
            END IF
         END DO
         DEALLOCATE (constraint_info%colvar_set)
      END IF
      IF (ASSOCIATED(constraint_info%const_colv_mol)) THEN
         DEALLOCATE (constraint_info%const_colv_mol)
      END IF
      IF (ASSOCIATED(constraint_info%const_colv_molname)) THEN
         DEALLOCATE (constraint_info%const_colv_molname)
      END IF
      IF (ASSOCIATED(constraint_info%const_colv_target)) THEN
         DEALLOCATE (constraint_info%const_colv_target)
      END IF
      IF (ASSOCIATED(constraint_info%const_colv_target_growth)) THEN
         DEALLOCATE (constraint_info%const_colv_target_growth)
      END IF
      IF (ASSOCIATED(constraint_info%colv_intermolecular)) THEN
         DEALLOCATE (constraint_info%colv_intermolecular)
      END IF
      IF (ASSOCIATED(constraint_info%colv_restraint)) THEN
         DEALLOCATE (constraint_info%colv_restraint)
      END IF
      IF (ASSOCIATED(constraint_info%colv_k0)) THEN
         DEALLOCATE (constraint_info%colv_k0)
      END IF
      IF (ASSOCIATED(constraint_info%colv_exclude_qm)) THEN
         DEALLOCATE (constraint_info%colv_exclude_qm)
      END IF
      IF (ASSOCIATED(constraint_info%colv_exclude_mm)) THEN
         DEALLOCATE (constraint_info%colv_exclude_mm)
      END IF
      ! G3x3
      IF (ASSOCIATED(constraint_info%const_g33_mol)) THEN
         DEALLOCATE (constraint_info%const_g33_mol)
      END IF
      IF (ASSOCIATED(constraint_info%const_g33_molname)) THEN
         DEALLOCATE (constraint_info%const_g33_molname)
      END IF
      IF (ASSOCIATED(constraint_info%const_g33_a)) THEN
         DEALLOCATE (constraint_info%const_g33_a)
      END IF
      IF (ASSOCIATED(constraint_info%const_g33_b)) THEN
         DEALLOCATE (constraint_info%const_g33_b)
      END IF
      IF (ASSOCIATED(constraint_info%const_g33_c)) THEN
         DEALLOCATE (constraint_info%const_g33_c)
      END IF
      IF (ASSOCIATED(constraint_info%const_g33_dab)) THEN
         DEALLOCATE (constraint_info%const_g33_dab)
      END IF
      IF (ASSOCIATED(constraint_info%const_g33_dac)) THEN
         DEALLOCATE (constraint_info%const_g33_dac)
      END IF
      IF (ASSOCIATED(constraint_info%const_g33_dbc)) THEN
         DEALLOCATE (constraint_info%const_g33_dbc)
      END IF
      IF (ASSOCIATED(constraint_info%g33_intermolecular)) THEN
         DEALLOCATE (constraint_info%g33_intermolecular)
      END IF
      IF (ASSOCIATED(constraint_info%g33_restraint)) THEN
         DEALLOCATE (constraint_info%g33_restraint)
      END IF
      IF (ASSOCIATED(constraint_info%g33_k0)) THEN
         DEALLOCATE (constraint_info%g33_k0)
      END IF
      IF (ASSOCIATED(constraint_info%g33_exclude_qm)) THEN
         DEALLOCATE (constraint_info%g33_exclude_qm)
      END IF
      IF (ASSOCIATED(constraint_info%g33_exclude_mm)) THEN
         DEALLOCATE (constraint_info%g33_exclude_mm)
      END IF
      ! G4x6
      IF (ASSOCIATED(constraint_info%const_g46_mol)) THEN
         DEALLOCATE (constraint_info%const_g46_mol)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_molname)) THEN
         DEALLOCATE (constraint_info%const_g46_molname)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_a)) THEN
         DEALLOCATE (constraint_info%const_g46_a)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_b)) THEN
         DEALLOCATE (constraint_info%const_g46_b)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_c)) THEN
         DEALLOCATE (constraint_info%const_g46_c)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_d)) THEN
         DEALLOCATE (constraint_info%const_g46_d)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_dab)) THEN
         DEALLOCATE (constraint_info%const_g46_dab)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_dac)) THEN
         DEALLOCATE (constraint_info%const_g46_dac)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_dbc)) THEN
         DEALLOCATE (constraint_info%const_g46_dbc)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_dad)) THEN
         DEALLOCATE (constraint_info%const_g46_dad)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_dbd)) THEN
         DEALLOCATE (constraint_info%const_g46_dbd)
      END IF
      IF (ASSOCIATED(constraint_info%const_g46_dcd)) THEN
         DEALLOCATE (constraint_info%const_g46_dcd)
      END IF
      IF (ASSOCIATED(constraint_info%g46_intermolecular)) THEN
         DEALLOCATE (constraint_info%g46_intermolecular)
      END IF
      IF (ASSOCIATED(constraint_info%g46_restraint)) THEN
         DEALLOCATE (constraint_info%g46_restraint)
      END IF
      IF (ASSOCIATED(constraint_info%g46_k0)) THEN
         DEALLOCATE (constraint_info%g46_k0)
      END IF
      IF (ASSOCIATED(constraint_info%g46_exclude_qm)) THEN
         DEALLOCATE (constraint_info%g46_exclude_qm)
      END IF
      IF (ASSOCIATED(constraint_info%g46_exclude_mm)) THEN
         DEALLOCATE (constraint_info%g46_exclude_mm)
      END IF
      ! virtual_site
      IF (ASSOCIATED(constraint_info%const_vsite_mol)) THEN
         DEALLOCATE (constraint_info%const_vsite_mol)
      END IF
      IF (ASSOCIATED(constraint_info%const_vsite_molname)) THEN
         DEALLOCATE (constraint_info%const_vsite_molname)
      END IF
      IF (ASSOCIATED(constraint_info%const_vsite_a)) THEN
         DEALLOCATE (constraint_info%const_vsite_a)
      END IF
      IF (ASSOCIATED(constraint_info%const_vsite_b)) THEN
         DEALLOCATE (constraint_info%const_vsite_b)
      END IF
      IF (ASSOCIATED(constraint_info%const_vsite_c)) THEN
         DEALLOCATE (constraint_info%const_vsite_c)
      END IF
      IF (ASSOCIATED(constraint_info%const_vsite_d)) THEN
         DEALLOCATE (constraint_info%const_vsite_d)
      END IF
      IF (ASSOCIATED(constraint_info%const_vsite_wbc)) THEN
         DEALLOCATE (constraint_info%const_vsite_wbc)
      END IF
      IF (ASSOCIATED(constraint_info%const_vsite_wdc)) THEN
         DEALLOCATE (constraint_info%const_vsite_wdc)
      END IF
      IF (ASSOCIATED(constraint_info%vsite_intermolecular)) THEN
         DEALLOCATE (constraint_info%vsite_intermolecular)
      END IF
      IF (ASSOCIATED(constraint_info%vsite_restraint)) THEN
         DEALLOCATE (constraint_info%vsite_restraint)
      END IF
      IF (ASSOCIATED(constraint_info%vsite_k0)) THEN
         DEALLOCATE (constraint_info%vsite_k0)
      END IF
      IF (ASSOCIATED(constraint_info%vsite_exclude_qm)) THEN
         DEALLOCATE (constraint_info%vsite_exclude_qm)
      END IF
      IF (ASSOCIATED(constraint_info%vsite_exclude_mm)) THEN
         DEALLOCATE (constraint_info%vsite_exclude_mm)
      END IF
   END SUBROUTINE deallocate_constraint

! **************************************************************************************************
!> \brief Deallocate possibly allocated arrays before reading topology
!> \param topology ...
!> \par History
!>      none
! **************************************************************************************************
   SUBROUTINE pre_read_topology(topology)
      TYPE(topology_parameters_type), INTENT(INOUT)      :: topology

      TYPE(atom_info_type), POINTER                      :: atom_info

      atom_info => topology%atom_info

      IF (ASSOCIATED(atom_info%id_molname)) THEN
         DEALLOCATE (atom_info%id_molname)
      END IF

      IF (ASSOCIATED(atom_info%resid)) THEN
         DEALLOCATE (atom_info%resid)
      END IF

      IF (ASSOCIATED(atom_info%id_resname)) THEN
         DEALLOCATE (atom_info%id_resname)
      END IF

      IF (ASSOCIATED(atom_info%id_atmname)) THEN
         DEALLOCATE (atom_info%id_atmname)
      END IF

      IF (ASSOCIATED(atom_info%atm_charge)) THEN
         DEALLOCATE (atom_info%atm_charge)
      END IF

      IF (ASSOCIATED(atom_info%atm_mass)) THEN
         DEALLOCATE (atom_info%atm_mass)
      END IF

   END SUBROUTINE pre_read_topology

END MODULE topology_types
