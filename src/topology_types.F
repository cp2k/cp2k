!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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

  USE kinds,                           ONLY: dp,&
                                             default_string_length
  USE pair_potential_types,            ONLY: pair_potential_type
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  TYPE atom_info_type
    CHARACTER (LEN=default_string_length), POINTER :: label_molname(:)
    CHARACTER (LEN=default_string_length), POINTER :: label_resname(:)
    CHARACTER (LEN=default_string_length), POINTER :: label_atmname(:)
    CHARACTER (LEN=default_string_length), POINTER :: atom_names(:)
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: r
    INTEGER, POINTER :: map_mol_typ(:)
    INTEGER, POINTER :: map_mol_num(:)
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
  END TYPE connectivity_info_type

  TYPE constraint_info_type
    INTEGER          :: nconst_dist
    INTEGER, POINTER :: const_dist_mol(:)
    INTEGER, POINTER :: const_dist_a(:)
    INTEGER, POINTER :: const_dist_b(:)
    REAL(KIND=dp), POINTER :: const_dist_dab(:)
    INTEGER          :: nconst_g33
    INTEGER, POINTER :: const_g33_mol(:)
    INTEGER, POINTER :: const_g33_a(:)
    INTEGER, POINTER :: const_g33_b(:)
    INTEGER, POINTER :: const_g33_c(:)
    REAL(KIND=dp), POINTER :: const_g33_dab(:)
    REAL(KIND=dp), POINTER :: const_g33_dac(:)
    REAL(KIND=dp), POINTER :: const_g33_dbc(:)
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
  END TYPE constraint_info_type

  TYPE topology_parameters_type
    TYPE (atom_info_type),POINTER                :: atom_info
    TYPE (connectivity_info_type),POINTER        :: conn_info
    TYPE (constraint_info_type),POINTER          :: cons_info
    TYPE ( cell_type ), POINTER                  :: cell, cell_ref
    !TRY TO REMOVE THIS FOUR VARIABLE IN THE FUTURE
    INTEGER                                   :: natoms, natom_type
    INTEGER                                   :: nmol, nmol_type
    !TRY TO REMOVE THIS FOUR VARIABLE IN THE FUTURE
    CHARACTER (LEN=default_string_length)     :: bondparm_type
    REAL(KIND=dp)                                :: bondparm_factor
    LOGICAL                                   :: reorder_atom
    LOGICAL                                   :: coordinate
    LOGICAL                                   :: connectivity
    LOGICAL                                   :: dump_topology
    CHARACTER (LEN=default_string_length)     :: nhcopt
    CHARACTER (LEN=default_string_length)     :: coord_type
    CHARACTER (LEN=default_string_length)     :: coord_file_name
    CHARACTER (LEN=default_string_length)     :: conn_type
    CHARACTER (LEN=default_string_length)     :: conn_file_name
    LOGICAL                                   :: constraint
    LOGICAL                                   :: const_atom
    LOGICAL                                   :: const_hydr
    LOGICAL                                   :: const_dist
    LOGICAL                                   :: const_33
    LOGICAL                                   :: const_46
    LOGICAL                                   :: charge_occup
    LOGICAL                                   :: charge_beta
    CHARACTER (LEN=default_string_length)     :: exclude_vdw
    CHARACTER (LEN=default_string_length)     :: exclude_ei
  END TYPE topology_parameters_type


  PUBLIC :: atom_info_type,&
            connectivity_info_type,&
            constraint_info_type,&
            topology_parameters_type

  PUBLIC :: init_topology,&
            deallocate_topology


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
  topology%bondparm_type='COVALENT'
  topology%bondparm_factor=1.1
  topology%reorder_atom=.FALSE.
  topology%coordinate=.FALSE.
  topology%connectivity=.FALSE.
  topology%dump_topology=.FALSE.
  topology%nhcopt='GLOBAL'
  topology%coord_type=''
  topology%coord_file_name=''
  topology%conn_type='GENERATE'
  topology%conn_file_name='OFF'
  topology%constraint=.FALSE.
  topology%const_atom=.FALSE.
  topology%const_hydr=.FALSE.
  topology%const_dist=.FALSE.
  topology%const_33=.FALSE.
  topology%const_46=.FALSE.
  topology%charge_occup=.FALSE.
  topology%charge_beta=.FALSE.
  topology%exclude_vdw="1-4"
  topology%exclude_ei="1-4"


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. Initialize and Nullify things in topology%atom_info
  !-----------------------------------------------------------------------------
  NULLIFY(topology%atom_info%label_molname)
  NULLIFY(topology%atom_info%label_resname)
  NULLIFY(topology%atom_info%label_atmname)
  NULLIFY(topology%atom_info%atom_names)
  NULLIFY(topology%atom_info%r)
  NULLIFY(topology%atom_info%map_mol_typ)
  NULLIFY(topology%atom_info%map_mol_num)
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


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. Initialize and Nullify things in topology%cons_info
  !-----------------------------------------------------------------------------
  topology%cons_info%nconst_dist = 0
  NULLIFY(topology%cons_info%const_dist_mol)
  NULLIFY(topology%cons_info%const_dist_a)
  NULLIFY(topology%cons_info%const_dist_b)
  NULLIFY(topology%cons_info%const_dist_dab)
  topology%cons_info%nconst_g33 = 0
  NULLIFY(topology%cons_info%const_g33_mol)
  NULLIFY(topology%cons_info%const_g33_a)
  NULLIFY(topology%cons_info%const_g33_b)
  NULLIFY(topology%cons_info%const_g33_c)
  NULLIFY(topology%cons_info%const_g33_dab)
  NULLIFY(topology%cons_info%const_g33_dac)
  NULLIFY(topology%cons_info%const_g33_dbc)
  topology%cons_info%nconst_g46 = 0
  NULLIFY(topology%cons_info%const_g46_mol)
  NULLIFY(topology%cons_info%const_g46_a)
  NULLIFY(topology%cons_info%const_g46_b)
  NULLIFY(topology%cons_info%const_g46_c)
  NULLIFY(topology%cons_info%const_g46_d)
  NULLIFY(topology%cons_info%const_g46_dab)
  NULLIFY(topology%cons_info%const_g46_dac)
  NULLIFY(topology%cons_info%const_g46_dbc)
  NULLIFY(topology%cons_info%const_g46_dad)
  NULLIFY(topology%cons_info%const_g46_dbd)
  NULLIFY(topology%cons_info%const_g46_dcd)



!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE init_topology



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
SUBROUTINE deallocate_topology (topology)
    TYPE(topology_parameters_type), &
      INTENT(INOUT)                          :: topology

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


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 3. DEALLOCATE things in topology%cons_info
  !-----------------------------------------------------------------------------
  IF(ASSOCIATED(topology%cons_info%const_dist_mol)) THEN
    DEALLOCATE(topology%cons_info%const_dist_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_dist_mol')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_dist_a)) THEN
    DEALLOCATE(topology%cons_info%const_dist_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_dist_a')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_dist_b)) THEN
    DEALLOCATE(topology%cons_info%const_dist_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_dist_b')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_dist_dab)) THEN
    DEALLOCATE(topology%cons_info%const_dist_dab,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_dist_dab')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g33_mol)) THEN
    DEALLOCATE(topology%cons_info%const_g33_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_mol')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g33_a)) THEN
    DEALLOCATE(topology%cons_info%const_g33_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_a')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g33_b)) THEN
    DEALLOCATE(topology%cons_info%const_g33_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_b')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g33_c)) THEN
    DEALLOCATE(topology%cons_info%const_g33_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_c')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g33_dab)) THEN
    DEALLOCATE(topology%cons_info%const_g33_dab,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_dab')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g33_dac)) THEN
    DEALLOCATE(topology%cons_info%const_g33_dac,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_dac')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g33_dbc)) THEN
    DEALLOCATE(topology%cons_info%const_g33_dbc,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_dbc')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_mol)) THEN
    DEALLOCATE(topology%cons_info%const_g46_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_mol')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_a)) THEN
    DEALLOCATE(topology%cons_info%const_g46_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_a')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_b)) THEN
    DEALLOCATE(topology%cons_info%const_g46_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_b')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_c)) THEN
    DEALLOCATE(topology%cons_info%const_g46_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_c')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_d)) THEN
    DEALLOCATE(topology%cons_info%const_g46_d,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_d')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_dab)) THEN
    DEALLOCATE(topology%cons_info%const_g46_dab,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_dab')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_dac)) THEN
    DEALLOCATE(topology%cons_info%const_g46_dac,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_dac')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_dbc)) THEN
    DEALLOCATE(topology%cons_info%const_g46_dbc,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_dbc')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_dad)) THEN
    DEALLOCATE(topology%cons_info%const_g46_dad,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_dad')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_dbd)) THEN
    DEALLOCATE(topology%cons_info%const_g46_dbd,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_dbd')
  END IF
  IF(ASSOCIATED(topology%cons_info%const_g46_dcd)) THEN
    DEALLOCATE(topology%cons_info%const_g46_dcd,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g46_dcd')
  END IF


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



END MODULE topology_types

!!*****
!******************************************************************************
