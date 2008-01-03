!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!


!!****s* cp2k/topology_types [1.0] *
!!
!!   NAME
!!     topology_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM & JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (30.11.2001) : new entries in setup_parameters_type
!!                        change name from input_file_name to coord_...
!!                        added topology file
!!                        added atom_names
!!
!!   SOURCE
!******************************************************************************

MODULE topology_types
  USE cell_types,                      ONLY: cell_release,&
                                             cell_type
  USE colvar_types,                    ONLY: colvar_p_type,&
                                             colvar_release
  USE f77_blas
  USE input_constants,                 ONLY: do_bondparm_covalent,&
                                             do_conn_generate,&
                                             do_constr_none,&
                                             do_skip_13
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE termination,                     ONLY: stop_memory
#include "cp_common_uses.h"

  IMPLICIT NONE

  TYPE atom_info_type
     CHARACTER (LEN=default_string_length), POINTER :: label_molname(:)
     CHARACTER (LEN=default_string_length), POINTER :: label_resname(:)
     INTEGER,                               POINTER :: label_resid(:)
     CHARACTER (LEN=default_string_length), POINTER :: label_atmname(:)
     CHARACTER (LEN=default_string_length), POINTER :: atom_names(:)
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: r
     INTEGER, POINTER :: map_mol_typ(:)
     INTEGER, POINTER :: map_mol_num(:)
     INTEGER, POINTER :: map_mol_res(:)
     REAL(KIND=dp), POINTER :: atm_charge(:)
     REAL(KIND=dp), POINTER :: atm_mass(:)
     REAL(KIND=dp), POINTER :: occup(:)
     REAL(KIND=dp), POINTER :: beta(:)
     CHARACTER (LEN=default_string_length), POINTER :: element(:)
  END TYPE atom_info_type

  TYPE connectivity_info_type
     INTEGER, POINTER :: bond_a(:),bond_b(:),bond_type(:)
     INTEGER, POINTER :: ub_a(:),ub_b(:),ub_c(:)
     INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:),theta_type(:)
     INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:),phi_type(:)
     INTEGER, POINTER :: impr_a(:),impr_b(:),impr_c(:),impr_d(:),impr_type(:)
     INTEGER, POINTER :: onfo_a(:),onfo_b(:),onfo_type(:)
     INTEGER, POINTER :: c_bond_a(:),c_bond_b(:),c_bond_type(:)
  END TYPE connectivity_info_type

  TYPE constraint_info_type
     ! Bonds involving Hydrogens
     LOGICAL                                       :: hbonds_restraint                        ! Restraints control
     REAL(KIND=dp)                                 :: hbonds_k0                               ! Restraints control
     ! Fixed Atoms
     INTEGER                                       :: nfixed_atoms
     INTEGER, POINTER                              :: fixed_atoms(:)
     LOGICAL, POINTER                              :: fixed_restraint(:)                      ! Restraints control
     REAL(KIND=dp), POINTER                        :: fixed_k0(:)                             ! Restraints control
     ! Freeze QM or MM
     INTEGER                                       :: freeze_qm, freeze_mm
     LOGICAL                                       :: fixed_mm_restraint, fixed_qm_restraint  ! Restraints control
     REAL(KIND=dp)                                 :: fixed_mm_k0, fixed_qm_k0                ! Restraints control
     ! Freeze with molnames
     LOGICAL, POINTER                              :: fixed_mol_restraint(:)                  ! Restraints control
     REAL(KIND=dp), POINTER                        :: fixed_mol_k0(:)                         ! Restraints control
     CHARACTER(LEN=default_string_length), POINTER :: fixed_molnames(:)
     LOGICAL, POINTER, DIMENSION(:)                :: fixed_exclude_qm, fixed_exclude_mm
     ! Collective constraints
     INTEGER                                       :: nconst_colv
     INTEGER, POINTER                              :: const_colv_mol(:)
     CHARACTER(LEN=default_string_length), POINTER :: const_colv_molname(:)
     REAL(KIND=dp), POINTER                        :: const_colv_target(:)
     TYPE(colvar_p_type), POINTER, DIMENSION(:)    :: colvar_set
     LOGICAL, POINTER                              :: colv_intermolecular(:)
     LOGICAL, POINTER                              :: colv_restraint(:)                       ! Restraints control
     REAL(KIND=dp), POINTER                        :: colv_k0(:)                              ! Restraints control
     LOGICAL, POINTER, DIMENSION(:)                :: colv_exclude_qm, colv_exclude_mm
     ! G3x3
     INTEGER                                       :: nconst_g33
     INTEGER, POINTER                              :: const_g33_mol(:)
     CHARACTER(LEN=default_string_length), POINTER :: const_g33_molname(:)
     INTEGER, POINTER                              :: const_g33_a(:)
     INTEGER, POINTER                              :: const_g33_b(:)
     INTEGER, POINTER                              :: const_g33_c(:)
     REAL(KIND=dp), POINTER                        :: const_g33_dab(:)
     REAL(KIND=dp), POINTER                        :: const_g33_dac(:)
     REAL(KIND=dp), POINTER                        :: const_g33_dbc(:)
     LOGICAL, POINTER                              :: g33_intermolecular(:)
     LOGICAL, POINTER                              :: g33_restraint(:)                        ! Restraints control
     REAL(KIND=dp), POINTER                        :: g33_k0(:)                               ! Restraints control
     LOGICAL, POINTER, DIMENSION(:)                :: g33_exclude_qm, g33_exclude_mm
     ! G4x6
     INTEGER                                       :: nconst_g46
     INTEGER, POINTER                              :: const_g46_mol(:)
     CHARACTER(LEN=default_string_length), POINTER :: const_g46_molname(:)
     INTEGER, POINTER                              :: const_g46_a(:)
     INTEGER, POINTER                              :: const_g46_b(:)
     INTEGER, POINTER                              :: const_g46_c(:)
     INTEGER, POINTER                              :: const_g46_d(:)
     REAL(KIND=dp), POINTER                        :: const_g46_dab(:)
     REAL(KIND=dp), POINTER                        :: const_g46_dac(:)
     REAL(KIND=dp), POINTER                        :: const_g46_dbc(:)
     REAL(KIND=dp), POINTER                        :: const_g46_dad(:)
     REAL(KIND=dp), POINTER                        :: const_g46_dbd(:)
     REAL(KIND=dp), POINTER                        :: const_g46_dcd(:)
     LOGICAL, POINTER                              :: g46_intermolecular(:)
     LOGICAL, POINTER                              :: g46_restraint(:)                        ! Restraints control
     REAL(KIND=dp), POINTER                        :: g46_k0(:)                               ! Restraints control
     LOGICAL, POINTER, DIMENSION(:)                :: g46_exclude_qm, g46_exclude_mm
  END TYPE constraint_info_type

  TYPE topology_parameters_type
     TYPE (atom_info_type),POINTER             :: atom_info
     TYPE (connectivity_info_type),POINTER     :: conn_info
     TYPE (constraint_info_type),POINTER       :: cons_info
     TYPE (cell_type), POINTER                 :: cell, cell_ref
     INTEGER                                   :: conn_type
     INTEGER                                   :: coord_type
     INTEGER                                   :: exclude_vdw
     INTEGER                                   :: exclude_ei
     INTEGER                                   :: bondparm_type
     !TRY TO REMOVE THIS FIVE VARIABLE IN THE FUTURE
     INTEGER                                   :: natoms, natom_type
     INTEGER                                   :: nmol, nmol_type, nmol_conn
     !TRY TO REMOVE THIS FIVE VARIABLE IN THE FUTURE
     LOGICAL                                   :: aa_element
     LOGICAL                                   :: molname_generated
     REAL(KIND=dp)                             :: bondparm_factor
     LOGICAL                                   :: create_molecules
     LOGICAL                                   :: reorder_atom
     LOGICAL                                   :: molecules_check
     LOGICAL                                   :: coordinate
     LOGICAL                                   :: use_g96_velocity
     CHARACTER (LEN=default_path_length)       :: coord_file_name
     CHARACTER (LEN=default_path_length)       :: conn_file_name
     LOGICAL                                   :: const_atom
     LOGICAL                                   :: const_hydr
     LOGICAL                                   :: const_colv
     LOGICAL                                   :: const_33
     LOGICAL                                   :: const_46
     LOGICAL                                   :: charge_occup
     LOGICAL                                   :: charge_beta
     LOGICAL                                   :: para_res
  END TYPE topology_parameters_type

  TYPE constr_list_type
     INTEGER, DIMENSION(:), POINTER :: constr
  END TYPE constr_list_type

  PUBLIC :: atom_info_type,&
            connectivity_info_type,&
            constraint_info_type,&
            topology_parameters_type,&
            constr_list_type

  PUBLIC :: init_topology,&
            deallocate_topology

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_types'
  PRIVATE


CONTAINS

  !******************************************************************************
  !!****** topology_types/init_topology [1.0] *
  !!
  !!   NAME
  !!     init_topology
  !!
  !!   FUNCTION
  !!     1. Just NULLIFY and zero all the stuff
  !!
  !!   AUTHOR
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE init_topology (topology)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology

    INTEGER                                  :: istat

!-----------------------------------------------------------------------------
! 1. Nullify and allocate things in topology
!-----------------------------------------------------------------------------

    ALLOCATE(topology%atom_info,STAT=istat)
    IF(istat/=0) CALL stop_memory ('init_topology','atom_info')
    ALLOCATE(topology%conn_info,STAT=istat)
    IF(istat/=0) CALL stop_memory ('init_topology','conn_info')
    ALLOCATE(topology%cons_info,STAT=istat)
    IF(istat/=0) CALL stop_memory ('init_topology','cons_info')
    !-----------------------------------------------------------------------------
    ! 2. Initialize and Nullify things in topology
    !-----------------------------------------------------------------------------
    NULLIFY(topology%cell, topology%cell_ref)
    topology%natoms=0
    topology%natom_type=0
    topology%nmol=0
    topology%nmol_type=0
    topology%nmol_conn=0
    topology%bondparm_type=do_bondparm_covalent
    topology%reorder_atom=.FALSE.
    topology%create_molecules=.FALSE.
    topology%molecules_check=.FALSE.
    topology%coordinate=.FALSE.
    topology%use_g96_velocity=.FALSE.
    topology%coord_type=-1
    topology%coord_file_name=''
    topology%conn_type=do_conn_generate
    topology%conn_file_name='OFF'
    topology%const_atom=.FALSE.
    topology%const_hydr=.FALSE.
    topology%const_colv=.FALSE.
    topology%const_33=.FALSE.
    topology%const_46=.FALSE.
    topology%charge_occup=.FALSE.
    topology%charge_beta=.FALSE.
    topology%para_res=.FALSE.
    topology%molname_generated=.FALSE.
    topology%aa_element=.FALSE.
    topology%exclude_vdw=do_skip_13
    topology%exclude_ei=do_skip_13
    !-----------------------------------------------------------------------------
    ! 3. Initialize and Nullify things in topology%atom_info
    !-----------------------------------------------------------------------------
    NULLIFY(topology%atom_info%label_molname)
    NULLIFY(topology%atom_info%label_resname)
    NULLIFY(topology%atom_info%label_resid)
    NULLIFY(topology%atom_info%label_atmname)
    NULLIFY(topology%atom_info%atom_names)
    NULLIFY(topology%atom_info%r)
    NULLIFY(topology%atom_info%map_mol_typ)
    NULLIFY(topology%atom_info%map_mol_num)
    NULLIFY(topology%atom_info%map_mol_res)
    NULLIFY(topology%atom_info%atm_charge)
    NULLIFY(topology%atom_info%atm_mass)
    NULLIFY(topology%atom_info%occup)
    NULLIFY(topology%atom_info%beta)
    NULLIFY(topology%atom_info%element)
    !-----------------------------------------------------------------------------
    ! 4. Initialize and Nullify things in topology%conn_info
    !-----------------------------------------------------------------------------
    NULLIFY(topology%conn_info%bond_a)
    NULLIFY(topology%conn_info%bond_b)
    NULLIFY(topology%conn_info%bond_type)
    NULLIFY(topology%conn_info%ub_a)
    NULLIFY(topology%conn_info%ub_b)
    NULLIFY(topology%conn_info%ub_c)
    NULLIFY(topology%conn_info%theta_a)
    NULLIFY(topology%conn_info%theta_b)
    NULLIFY(topology%conn_info%theta_c)
    NULLIFY(topology%conn_info%theta_type)
    NULLIFY(topology%conn_info%phi_a)
    NULLIFY(topology%conn_info%phi_b)
    NULLIFY(topology%conn_info%phi_c)
    NULLIFY(topology%conn_info%phi_d)
    NULLIFY(topology%conn_info%phi_type)
    NULLIFY(topology%conn_info%impr_a)
    NULLIFY(topology%conn_info%impr_b)
    NULLIFY(topology%conn_info%impr_c)
    NULLIFY(topology%conn_info%impr_d)
    NULLIFY(topology%conn_info%impr_type)
    NULLIFY(topology%conn_info%onfo_a)
    NULLIFY(topology%conn_info%onfo_b)
    NULLIFY(topology%conn_info%onfo_type)
    NULLIFY(topology%conn_info%c_bond_a)
    NULLIFY(topology%conn_info%c_bond_b)
    NULLIFY(topology%conn_info%c_bond_type)
    !-----------------------------------------------------------------------------
    ! 5. Initialize and Nullify things in topology%cons_info
    !-----------------------------------------------------------------------------
    CALL init_constraint(topology%cons_info)
  END SUBROUTINE init_topology


  SUBROUTINE init_constraint(constraint_info)
    TYPE(constraint_info_type), POINTER      :: constraint_info

! Bonds involving Hydrogens

    constraint_info%hbonds_restraint = .FALSE.
    ! Fixed Atoms
    constraint_info%nfixed_atoms = 0
    constraint_info%freeze_mm = do_constr_none
    constraint_info%freeze_qm = do_constr_none
    NULLIFY(constraint_info%fixed_atoms)
    NULLIFY(constraint_info%fixed_molnames)
    NULLIFY(constraint_info%fixed_restraint)
    NULLIFY(constraint_info%fixed_k0)
    NULLIFY(constraint_info%fixed_mol_restraint)
    NULLIFY(constraint_info%fixed_mol_k0)
    NULLIFY(constraint_info%fixed_exclude_qm, constraint_info%fixed_exclude_mm)
    ! Collective Constraints
    constraint_info%nconst_colv = 0
    NULLIFY(constraint_info%colvar_set)
    NULLIFY(constraint_info%const_colv_mol)
    NULLIFY(constraint_info%const_colv_molname)
    NULLIFY(constraint_info%const_colv_target)
    NULLIFY(constraint_info%colv_intermolecular)
    NULLIFY(constraint_info%colv_restraint)
    NULLIFY(constraint_info%colv_k0)
    NULLIFY(constraint_info%colv_exclude_qm, constraint_info%colv_exclude_mm)
    ! G3x3
    constraint_info%nconst_g33 = 0
    NULLIFY(constraint_info%const_g33_mol)
    NULLIFY(constraint_info%const_g33_molname)
    NULLIFY(constraint_info%const_g33_a)
    NULLIFY(constraint_info%const_g33_b)
    NULLIFY(constraint_info%const_g33_c)
    NULLIFY(constraint_info%const_g33_dab)
    NULLIFY(constraint_info%const_g33_dac)
    NULLIFY(constraint_info%const_g33_dbc)
    NULLIFY(constraint_info%g33_intermolecular)
    NULLIFY(constraint_info%g33_restraint)
    NULLIFY(constraint_info%g33_k0)
    NULLIFY(constraint_info%g33_exclude_qm, constraint_info%g33_exclude_mm)
    ! G4x6
    constraint_info%nconst_g46 = 0
    NULLIFY(constraint_info%const_g46_mol)
    NULLIFY(constraint_info%const_g46_molname)
    NULLIFY(constraint_info%const_g46_a)
    NULLIFY(constraint_info%const_g46_b)
    NULLIFY(constraint_info%const_g46_c)
    NULLIFY(constraint_info%const_g46_d)
    NULLIFY(constraint_info%const_g46_dab)
    NULLIFY(constraint_info%const_g46_dac)
    NULLIFY(constraint_info%const_g46_dbc)
    NULLIFY(constraint_info%const_g46_dad)
    NULLIFY(constraint_info%const_g46_dbd)
    NULLIFY(constraint_info%const_g46_dcd)
    NULLIFY(constraint_info%g46_intermolecular)
    NULLIFY(constraint_info%g46_restraint)
    NULLIFY(constraint_info%g46_k0)
    NULLIFY(constraint_info%g46_exclude_qm, constraint_info%g46_exclude_mm)

  END SUBROUTINE init_constraint

  !******************************************************************************
  !!****** topology_types/deallocate_topology [1.0] *
  !!
  !!   NAME
  !!     deallocate_topology
  !!
  !!   FUNCTION
  !!     1. Just DEALLOCATE all the stuff
  !!
  !!   AUTHOR
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !******************************************************************************
  SUBROUTINE deallocate_topology (topology, error)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: istat

!-----------------------------------------------------------------------------
! 1. DEALLOCATE things in topology%atom_info
!-----------------------------------------------------------------------------

    IF(ASSOCIATED(topology%atom_info%label_molname)) THEN
       DEALLOCATE(topology%atom_info%label_molname,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','label_molname')
    END IF
    IF(ASSOCIATED(topology%atom_info%label_resname)) THEN
       DEALLOCATE(topology%atom_info%label_resname,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','label_resname')
    END IF
    IF(ASSOCIATED(topology%atom_info%label_resid)) THEN
       DEALLOCATE(topology%atom_info%label_resid,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','label_resid')
    END IF
    IF(ASSOCIATED(topology%atom_info%label_atmname)) THEN
       DEALLOCATE(topology%atom_info%label_atmname,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','label_atmname')
    END IF
    IF(ASSOCIATED(topology%atom_info%atom_names)) THEN
       DEALLOCATE(topology%atom_info%atom_names,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','atom_names')
    END IF
    IF(ASSOCIATED(topology%atom_info%r)) THEN
       DEALLOCATE(topology%atom_info%r,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','r')
    END IF
    IF(ASSOCIATED(topology%atom_info%map_mol_typ)) THEN
       DEALLOCATE(topology%atom_info%map_mol_typ,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','map_mol_typ')
    END IF
    IF(ASSOCIATED(topology%atom_info%map_mol_num)) THEN
       DEALLOCATE(topology%atom_info%map_mol_num,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','map_mol_num')
    END IF
    IF(ASSOCIATED(topology%atom_info%map_mol_res)) THEN
       DEALLOCATE(topology%atom_info%map_mol_res,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','map_mol_res')
    END IF
    IF(ASSOCIATED(topology%atom_info%atm_charge)) THEN
       DEALLOCATE(topology%atom_info%atm_charge,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','atm_charge')
    END IF
    IF(ASSOCIATED(topology%atom_info%atm_mass)) THEN
       DEALLOCATE(topology%atom_info%atm_mass,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','atm_mass')
    END IF
    IF(ASSOCIATED(topology%atom_info%occup)) THEN
       DEALLOCATE(topology%atom_info%occup,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','occup')
    END IF
    IF(ASSOCIATED(topology%atom_info%beta)) THEN
       DEALLOCATE(topology%atom_info%beta,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','beta')
    END IF
    IF(ASSOCIATED(topology%atom_info%element)) THEN
       DEALLOCATE(topology%atom_info%element,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','element')
    END IF
    !-----------------------------------------------------------------------------
    ! 2. DEALLOCATE things in topology%conn_info
    !-----------------------------------------------------------------------------
    IF(ASSOCIATED(topology%conn_info%bond_a)) THEN
       DEALLOCATE(topology%conn_info%bond_a,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','bond_a')
    END IF
    IF(ASSOCIATED(topology%conn_info%bond_b)) THEN
       DEALLOCATE(topology%conn_info%bond_b,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','bond_b')
    END IF
    IF(ASSOCIATED(topology%conn_info%bond_type)) THEN
       DEALLOCATE(topology%conn_info%bond_type,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','bond_type')
    END IF
    IF(ASSOCIATED(topology%conn_info%ub_a)) THEN
       DEALLOCATE(topology%conn_info%ub_a,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','ub_a')
    END IF
    IF(ASSOCIATED(topology%conn_info%ub_b)) THEN
       DEALLOCATE(topology%conn_info%ub_b,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','ub_b')
    END IF
    IF(ASSOCIATED(topology%conn_info%ub_c)) THEN
       DEALLOCATE(topology%conn_info%ub_c,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','ub_c')
    END IF
    IF(ASSOCIATED(topology%conn_info%theta_a)) THEN
       DEALLOCATE(topology%conn_info%theta_a,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_a')
    END IF
    IF(ASSOCIATED(topology%conn_info%theta_b)) THEN
       DEALLOCATE(topology%conn_info%theta_b,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_b')
    END IF
    IF(ASSOCIATED(topology%conn_info%theta_c)) THEN
       DEALLOCATE(topology%conn_info%theta_c,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_c')
    END IF
    IF(ASSOCIATED(topology%conn_info%theta_type)) THEN
       DEALLOCATE(topology%conn_info%theta_type,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_type')
    END IF
    IF(ASSOCIATED(topology%conn_info%phi_a)) THEN
       DEALLOCATE(topology%conn_info%phi_a,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_a')
    END IF
    IF(ASSOCIATED(topology%conn_info%phi_b)) THEN
       DEALLOCATE(topology%conn_info%phi_b,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_b')
    END IF
    IF(ASSOCIATED(topology%conn_info%phi_c)) THEN
       DEALLOCATE(topology%conn_info%phi_c,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_c')
    END IF
    IF(ASSOCIATED(topology%conn_info%phi_d)) THEN
       DEALLOCATE(topology%conn_info%phi_d,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_d')
    END IF
    IF(ASSOCIATED(topology%conn_info%phi_type)) THEN
       DEALLOCATE(topology%conn_info%phi_type,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_type')
    END IF
    IF(ASSOCIATED(topology%conn_info%impr_a)) THEN
       DEALLOCATE(topology%conn_info%impr_a,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','impr_a')
    END IF
    IF(ASSOCIATED(topology%conn_info%impr_b)) THEN
       DEALLOCATE(topology%conn_info%impr_b,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','impr_b')
    END IF
    IF(ASSOCIATED(topology%conn_info%impr_c)) THEN
       DEALLOCATE(topology%conn_info%impr_c,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','impr_c')
    END IF
    IF(ASSOCIATED(topology%conn_info%impr_d)) THEN
       DEALLOCATE(topology%conn_info%impr_d,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','impr_d')
    END IF
    IF(ASSOCIATED(topology%conn_info%impr_type)) THEN
       DEALLOCATE(topology%conn_info%impr_type,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','impr_type')
    END IF
    IF(ASSOCIATED(topology%conn_info%onfo_a)) THEN
       DEALLOCATE(topology%conn_info%onfo_a,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','onfo_a')
    END IF
    IF(ASSOCIATED(topology%conn_info%onfo_b)) THEN
       DEALLOCATE(topology%conn_info%onfo_b,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','onfo_b')
    END IF
    IF(ASSOCIATED(topology%conn_info%onfo_type)) THEN
       DEALLOCATE(topology%conn_info%onfo_type,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','onfo_type')
    END IF
    IF(ASSOCIATED(topology%conn_info%c_bond_a)) THEN
       DEALLOCATE(topology%conn_info%c_bond_a,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','c_bond_a')
    END IF
    IF(ASSOCIATED(topology%conn_info%c_bond_b)) THEN
       DEALLOCATE(topology%conn_info%c_bond_b,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','c_bond_b')
    END IF
    IF(ASSOCIATED(topology%conn_info%c_bond_type)) THEN
       DEALLOCATE(topology%conn_info%c_bond_type,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','c_bond_type')
    END IF
    !-----------------------------------------------------------------------------
    ! 3. DEALLOCATE things in topology%cons_info
    !-----------------------------------------------------------------------------
    IF (ASSOCIATED(topology%cons_info)) &
         CALL deallocate_constraint(topology%cons_info, error)
    !-----------------------------------------------------------------------------
    ! 4. DEALLOCATE things in topology
    !-----------------------------------------------------------------------------
    CALL cell_release(topology%cell,error=error)
    CALL cell_release(topology%cell_ref,error=error)
    IF(ASSOCIATED(topology%atom_info)) THEN
       DEALLOCATE(topology%atom_info,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','atom_info')
    END IF
    IF(ASSOCIATED(topology%conn_info)) THEN
       DEALLOCATE(topology%conn_info,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','conn_info')
    END IF
    IF(ASSOCIATED(topology%cons_info)) THEN
       DEALLOCATE(topology%cons_info,STAT=istat)
       IF (istat /= 0) CALL stop_memory('deallocate_topology','cons_info')
    END IF

  END SUBROUTINE deallocate_topology



  SUBROUTINE deallocate_constraint(constraint_info, error)
    TYPE(constraint_info_type), POINTER      :: constraint_info
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_constraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (.NOT.failure) THEN
       ! Fixed Atoms
       IF(ASSOCIATED(constraint_info%fixed_atoms)) THEN
          DEALLOCATE(constraint_info%fixed_atoms, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%fixed_molnames)) THEN
          DEALLOCATE(constraint_info%fixed_molnames, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%fixed_restraint)) THEN
          DEALLOCATE(constraint_info%fixed_restraint, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%fixed_k0)) THEN
          DEALLOCATE(constraint_info%fixed_k0, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%fixed_mol_restraint)) THEN
          DEALLOCATE(constraint_info%fixed_mol_restraint, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%fixed_mol_k0)) THEN
          DEALLOCATE(constraint_info%fixed_mol_k0, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%fixed_exclude_qm)) THEN
          DEALLOCATE(constraint_info%fixed_exclude_qm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%fixed_exclude_mm)) THEN
          DEALLOCATE(constraint_info%fixed_exclude_mm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       ! Collective Constraint
       IF (ASSOCIATED(constraint_info%colvar_set)) THEN
          DO i = 1, SIZE(constraint_info%colvar_set)
             IF (ASSOCIATED(constraint_info%colvar_set(i)%colvar)) THEN
                CALL colvar_release(constraint_info%colvar_set(i)%colvar,error)
                NULLIFY(constraint_info%colvar_set(i)%colvar)
             END IF
          END DO
          DEALLOCATE(constraint_info%colvar_set, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_colv_mol)) THEN
          DEALLOCATE(constraint_info%const_colv_mol,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_colv_molname)) THEN
          DEALLOCATE(constraint_info%const_colv_molname,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_colv_target)) THEN
          DEALLOCATE(constraint_info%const_colv_target,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%colv_intermolecular)) THEN
          DEALLOCATE(constraint_info%colv_intermolecular,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%colv_restraint)) THEN
          DEALLOCATE(constraint_info%colv_restraint,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%colv_k0)) THEN
          DEALLOCATE(constraint_info%colv_k0,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%colv_exclude_qm)) THEN
          DEALLOCATE(constraint_info%colv_exclude_qm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%colv_exclude_mm)) THEN
          DEALLOCATE(constraint_info%colv_exclude_mm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       ! G3x3
       IF(ASSOCIATED(constraint_info%const_g33_mol)) THEN
          DEALLOCATE(constraint_info%const_g33_mol,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g33_molname)) THEN
          DEALLOCATE(constraint_info%const_g33_molname,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g33_a)) THEN
          DEALLOCATE(constraint_info%const_g33_a,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g33_b)) THEN
          DEALLOCATE(constraint_info%const_g33_b,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g33_c)) THEN
          DEALLOCATE(constraint_info%const_g33_c,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g33_dab)) THEN
          DEALLOCATE(constraint_info%const_g33_dab,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g33_dac)) THEN
          DEALLOCATE(constraint_info%const_g33_dac,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g33_dbc)) THEN
          DEALLOCATE(constraint_info%const_g33_dbc,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g33_intermolecular)) THEN
          DEALLOCATE(constraint_info%g33_intermolecular,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g33_restraint)) THEN
          DEALLOCATE(constraint_info%g33_restraint,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g33_k0)) THEN
          DEALLOCATE(constraint_info%g33_k0,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g33_exclude_qm)) THEN
          DEALLOCATE(constraint_info%g33_exclude_qm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g33_exclude_mm)) THEN
          DEALLOCATE(constraint_info%g33_exclude_mm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       ! G4x6
       IF(ASSOCIATED(constraint_info%const_g46_mol)) THEN
          DEALLOCATE(constraint_info%const_g46_mol,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_molname)) THEN
          DEALLOCATE(constraint_info%const_g46_molname,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_a)) THEN
          DEALLOCATE(constraint_info%const_g46_a,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_b)) THEN
          DEALLOCATE(constraint_info%const_g46_b,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_c)) THEN
          DEALLOCATE(constraint_info%const_g46_c,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_d)) THEN
          DEALLOCATE(constraint_info%const_g46_d,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_dab)) THEN
          DEALLOCATE(constraint_info%const_g46_dab,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_dac)) THEN
          DEALLOCATE(constraint_info%const_g46_dac,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_dbc)) THEN
          DEALLOCATE(constraint_info%const_g46_dbc,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_dad)) THEN
          DEALLOCATE(constraint_info%const_g46_dad,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_dbd)) THEN
          DEALLOCATE(constraint_info%const_g46_dbd,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%const_g46_dcd)) THEN
          DEALLOCATE(constraint_info%const_g46_dcd,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g46_intermolecular)) THEN
          DEALLOCATE(constraint_info%g46_intermolecular,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g46_restraint)) THEN
          DEALLOCATE(constraint_info%g46_restraint,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g46_k0)) THEN
          DEALLOCATE(constraint_info%g46_k0,STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g46_exclude_qm)) THEN
          DEALLOCATE(constraint_info%g46_exclude_qm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF(ASSOCIATED(constraint_info%g46_exclude_mm)) THEN
          DEALLOCATE(constraint_info%g46_exclude_mm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE deallocate_constraint

END MODULE topology_types
!******************************************************************************
