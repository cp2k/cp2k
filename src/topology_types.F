!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000,2005  CP2K developers group                            !
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
  USE input_constants,                 ONLY: do_bondparm_covalent,&
                                             do_conn_generate,&
                                             do_nhcopt_global,&
                                             do_skip_ei_14,&
                                             do_skip_vdw_14
  USE kinds,                           ONLY: default_string_length,&
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
    INTEGER, POINTER :: bond_a(:),bond_b(:)
    INTEGER, POINTER :: ub_a(:),ub_b(:),ub_c(:)
    INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)
    INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)
    INTEGER, POINTER :: impr_a(:),impr_b(:),impr_c(:),impr_d(:)
    INTEGER, POINTER :: onfo_a(:),onfo_b(:)
    INTEGER, POINTER :: c_bond_a(:),c_bond_b(:)
  END TYPE connectivity_info_type

  TYPE constraint_info_type
     ! Fixed Atoms
     INTEGER          :: nfixed_atoms
     INTEGER, POINTER :: fixed_atoms(:)
     LOGICAL          :: freeze_qm, freeze_mm
     CHARACTER(LEN=default_string_length), POINTER :: fixed_molnames(:)
     ! Distance
     INTEGER          :: nconst_dist
     INTEGER, POINTER :: const_dist_mol(:)
     INTEGER, POINTER :: const_dist_a(:)
     INTEGER, POINTER :: const_dist_b(:)
     REAL(KIND=dp), POINTER :: const_dist_dab(:)
     ! Angle
     INTEGER          :: nconst_angle
     INTEGER, POINTER :: const_angle_mol(:)
     INTEGER, POINTER :: const_angle_a(:)
     INTEGER, POINTER :: const_angle_b(:)
     INTEGER, POINTER :: const_angle_c(:)
     REAL(KIND=dp), POINTER :: const_angle_abc(:)
     ! Torsion
     INTEGER          :: nconst_tors
     INTEGER, POINTER :: const_tors_mol(:)
     INTEGER, POINTER :: const_tors_a(:)
     INTEGER, POINTER :: const_tors_b(:)
     INTEGER, POINTER :: const_tors_c(:)
     INTEGER, POINTER :: const_tors_d(:)
     REAL(KIND=dp), POINTER :: const_tors_abcd(:)
     ! G3x3
     INTEGER          :: nconst_g33
     INTEGER, POINTER :: const_g33_mol(:)
     INTEGER, POINTER :: const_g33_a(:)
     INTEGER, POINTER :: const_g33_b(:)
     INTEGER, POINTER :: const_g33_c(:)
     REAL(KIND=dp), POINTER :: const_g33_dab(:)
     REAL(KIND=dp), POINTER :: const_g33_dac(:)
     REAL(KIND=dp), POINTER :: const_g33_dbc(:)
     ! G4x6
     INTEGER          :: nconst_g46
     INTEGER, POINTER :: const_g46_mol(:)
     INTEGER, POINTER :: const_g46_a(:)
     INTEGER, POINTER :: const_g46_b(:)
     INTEGER, POINTER :: const_g46_c(:)
     INTEGER, POINTER :: const_g46_d(:)
     REAL(KIND=dp), POINTER :: const_g46_dab(:)
     REAL(KIND=dp), POINTER :: const_g46_dac(:)
     REAL(KIND=dp), POINTER :: const_g46_dbc(:)
     REAL(KIND=dp), POINTER :: const_g46_dad(:)
     REAL(KIND=dp), POINTER :: const_g46_dbd(:)
     REAL(KIND=dp), POINTER :: const_g46_dcd(:)
     ! Collectives
     INTEGER          :: ncollective
  END TYPE constraint_info_type

  TYPE topology_parameters_type
    TYPE (atom_info_type),POINTER                :: atom_info
    TYPE (connectivity_info_type),POINTER        :: conn_info
    TYPE (constraint_info_type),POINTER          :: cons_info
    TYPE ( cell_type ), POINTER                  :: cell, cell_ref
    INTEGER                                   :: conn_type
    INTEGER                                   :: coord_type, nhcopt
    INTEGER                                   :: exclude_vdw
    INTEGER                                   :: exclude_ei
    INTEGER                                   :: bondparm_type
    !TRY TO REMOVE THIS FOUR VARIABLE IN THE FUTURE
    INTEGER                                   :: natoms, natom_type
    INTEGER                                   :: nmol, nmol_type, nmol_conn
    !TRY TO REMOVE THIS FOUR VARIABLE IN THE FUTURE
    REAL(KIND=dp)                             :: bondparm_factor
    LOGICAL                                   :: reorder_atom
    LOGICAL                                   :: molecules_check
    LOGICAL                                   :: coordinate
    LOGICAL                                   :: connectivity
    LOGICAL                                   :: dump_topology
    CHARACTER (LEN=default_string_length)     :: coord_file_name
    CHARACTER (LEN=default_string_length)     :: conn_file_name
    LOGICAL                                   :: const_atom
    LOGICAL                                   :: const_hydr
    LOGICAL                                   :: const_dist
    LOGICAL                                   :: const_33
    LOGICAL                                   :: const_46
    LOGICAL                                   :: charge_occup
    LOGICAL                                   :: charge_beta
    INTEGER                                   :: exclude_check
    LOGICAL                                   :: para_res
  END TYPE topology_parameters_type


  PUBLIC :: atom_info_type,&
            connectivity_info_type,&
            constraint_info_type,&
            topology_parameters_type

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

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
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
  topology%bondparm_factor=1.1
  topology%reorder_atom=.FALSE.
  topology%molecules_check=.FALSE.
  topology%coordinate=.FALSE.
  topology%connectivity=.FALSE.
  topology%dump_topology=.FALSE.
  topology%nhcopt=do_nhcopt_global
  topology%coord_type=-1
  topology%coord_file_name=''
  topology%conn_type=do_conn_generate
  topology%conn_file_name='OFF'
  topology%const_atom=.FALSE.
  topology%const_hydr=.FALSE.
  topology%const_dist=.FALSE.
  topology%const_33=.FALSE.
  topology%const_46=.FALSE.
  topology%charge_occup=.FALSE.
  topology%charge_beta=.FALSE.
  topology%exclude_vdw=do_skip_vdw_14
  topology%exclude_ei=do_skip_ei_14
  topology%exclude_check=3
  topology%para_res=.FALSE.


  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------
  ! 4. Initialize and Nullify things in topology%conn_info
  !-----------------------------------------------------------------------------
  NULLIFY(topology%conn_info%bond_a)
  NULLIFY(topology%conn_info%bond_b)
  NULLIFY(topology%conn_info%ub_a)
  NULLIFY(topology%conn_info%ub_b)
  NULLIFY(topology%conn_info%ub_c)
  NULLIFY(topology%conn_info%theta_a)
  NULLIFY(topology%conn_info%theta_b)
  NULLIFY(topology%conn_info%theta_c)
  NULLIFY(topology%conn_info%phi_a)
  NULLIFY(topology%conn_info%phi_b)
  NULLIFY(topology%conn_info%phi_c)
  NULLIFY(topology%conn_info%phi_d)
  NULLIFY(topology%conn_info%impr_a)
  NULLIFY(topology%conn_info%impr_b)
  NULLIFY(topology%conn_info%impr_c)
  NULLIFY(topology%conn_info%impr_d)
  NULLIFY(topology%conn_info%onfo_a)
  NULLIFY(topology%conn_info%onfo_b)
  NULLIFY(topology%conn_info%c_bond_a)
  NULLIFY(topology%conn_info%c_bond_b)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Initialize and Nullify things in topology%cons_info
  !-----------------------------------------------------------------------------
  CALL init_constraint(topology%cons_info)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE init_topology


SUBROUTINE init_constraint(constraint_info)
    TYPE(constraint_info_type), POINTER      :: constraint_info

! Fixed Atoms

  constraint_info%nfixed_atoms = 0
  constraint_info%freeze_mm = .FALSE.
  constraint_info%freeze_qm = .FALSE.
  NULLIFY(constraint_info%fixed_atoms)
  NULLIFY(constraint_info%fixed_molnames)
  ! Distance
  constraint_info%nconst_dist = 0
  NULLIFY(constraint_info%const_dist_mol)
  NULLIFY(constraint_info%const_dist_a)
  NULLIFY(constraint_info%const_dist_b)
  NULLIFY(constraint_info%const_dist_dab)
  ! Angle
  constraint_info%nconst_angle=0
  NULLIFY(constraint_info%const_angle_mol)
  NULLIFY(constraint_info%const_angle_a)
  NULLIFY(constraint_info%const_angle_b)
  NULLIFY(constraint_info%const_angle_c)
  NULLIFY(constraint_info%const_angle_abc)
  ! Torsion
  constraint_info%nconst_tors = 0
  NULLIFY(constraint_info%const_tors_mol)
  NULLIFY(constraint_info%const_tors_a)
  NULLIFY(constraint_info%const_tors_b)
  NULLIFY(constraint_info%const_tors_c)
  NULLIFY(constraint_info%const_tors_d)
  NULLIFY(constraint_info%const_tors_abcd)
  ! G3x3
  constraint_info%nconst_g33 = 0
  NULLIFY(constraint_info%const_g33_mol)
  NULLIFY(constraint_info%const_g33_a)
  NULLIFY(constraint_info%const_g33_b)
  NULLIFY(constraint_info%const_g33_c)
  NULLIFY(constraint_info%const_g33_dab)
  NULLIFY(constraint_info%const_g33_dac)
  NULLIFY(constraint_info%const_g33_dbc)
  ! G4x6
  constraint_info%nconst_g46 = 0
  NULLIFY(constraint_info%const_g46_mol)
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: istat

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
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
  IF(ASSOCIATED(topology%conn_info%onfo_a)) THEN
    DEALLOCATE(topology%conn_info%onfo_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','onfo_a')
  END IF
  IF(ASSOCIATED(topology%conn_info%onfo_b)) THEN
    DEALLOCATE(topology%conn_info%onfo_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','onfo_b')
  END IF
  IF(ASSOCIATED(topology%conn_info%c_bond_a)) THEN
    DEALLOCATE(topology%conn_info%c_bond_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','c_bond_a')
  END IF
  IF(ASSOCIATED(topology%conn_info%c_bond_b)) THEN
    DEALLOCATE(topology%conn_info%c_bond_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','c_bond_b')
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. DEALLOCATE things in topology%cons_info
  !-----------------------------------------------------------------------------
  IF (ASSOCIATED(topology%cons_info)) &
       CALL deallocate_constraint(topology%cons_info, error)
  
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. DEALLOCATE things in topology
  !-----------------------------------------------------------------------------
  CALL cell_release(topology%cell)
  CALL cell_release(topology%cell_ref)
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



!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE deallocate_topology



SUBROUTINE deallocate_constraint(constraint_info, error)
    TYPE(constraint_info_type), POINTER      :: constraint_info
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_constraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
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
     ! Distance
     IF(ASSOCIATED(constraint_info%const_dist_mol)) THEN
        DEALLOCATE(constraint_info%const_dist_mol,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_dist_a)) THEN
        DEALLOCATE(constraint_info%const_dist_a,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_dist_b)) THEN
        DEALLOCATE(constraint_info%const_dist_b,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_dist_dab)) THEN
        DEALLOCATE(constraint_info%const_dist_dab,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     ! Angle
     IF(ASSOCIATED(constraint_info%const_angle_mol)) THEN
        DEALLOCATE(constraint_info%const_angle_mol,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_angle_a)) THEN
        DEALLOCATE(constraint_info%const_angle_a,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_angle_b)) THEN
        DEALLOCATE(constraint_info%const_angle_b,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_angle_c)) THEN
        DEALLOCATE(constraint_info%const_angle_c,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_angle_abc)) THEN
        DEALLOCATE(constraint_info%const_angle_abc,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF     
     ! Torsion
     IF(ASSOCIATED(constraint_info%const_tors_mol)) THEN
        DEALLOCATE(constraint_info%const_tors_mol,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_tors_a)) THEN
        DEALLOCATE(constraint_info%const_tors_a,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_tors_b)) THEN
        DEALLOCATE(constraint_info%const_tors_b,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF
     IF(ASSOCIATED(constraint_info%const_tors_c)) THEN
        DEALLOCATE(constraint_info%const_tors_c,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF 
     IF(ASSOCIATED(constraint_info%const_tors_d)) THEN
        DEALLOCATE(constraint_info%const_tors_d,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF     
     IF(ASSOCIATED(constraint_info%const_tors_abcd)) THEN
        DEALLOCATE(constraint_info%const_tors_abcd,STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF     
     ! G3x3
     IF(ASSOCIATED(constraint_info%const_g33_mol)) THEN
        DEALLOCATE(constraint_info%const_g33_mol,STAT=stat)
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
     ! G4x6
     IF(ASSOCIATED(constraint_info%const_g46_mol)) THEN
        DEALLOCATE(constraint_info%const_g46_mol,STAT=stat)
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
  END IF
END SUBROUTINE deallocate_constraint

END MODULE topology_types

!!*****
!******************************************************************************
