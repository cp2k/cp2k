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

  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  TYPE topology_parameters_type
    TYPE ( cell_type ), POINTER               :: cell
    INTEGER                                   :: natoms, natom_type
    INTEGER                                   :: nmol, nmol_type
    LOGICAL                                   :: reorder_atom
    LOGICAL                                   :: coordinate
    LOGICAL                                   :: connectivity
    LOGICAL                                   :: dump_topology
    CHARACTER (LEN=default_string_length)     :: nhcopt
    CHARACTER (LEN=default_string_length)     :: bondparm_type
    REAL (dbl)                                :: bondparm_factor
    CHARACTER (LEN=default_string_length)     :: coord_type
    CHARACTER (LEN=default_string_length)     :: coord_file_name
    CHARACTER (LEN=default_string_length)     :: conn_type
    CHARACTER (LEN=default_string_length)     :: conn_file_name
    CHARACTER (LEN=default_string_length), POINTER :: label_molname(:)
    CHARACTER (LEN=default_string_length), POINTER :: label_resname(:)
    CHARACTER (LEN=default_string_length), POINTER :: label_atmname(:)
    CHARACTER (LEN=default_string_length), POINTER :: atom_names(:)
    REAL(dbl), DIMENSION(:,:), POINTER :: r
    INTEGER, POINTER :: map_mol_typ(:)
    INTEGER, POINTER :: map_mol_num(:)
    INTEGER, POINTER :: bond_a(:),bond_b(:)
    INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)
    INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)
    REAL (dbl), POINTER :: atm_charge(:)
    REAL (dbl), POINTER :: atm_mass(:)
    REAL (dbl), POINTER :: occup(:)
    REAL (dbl), POINTER :: beta(:)

    LOGICAL                                   :: constraint
    LOGICAL                                   :: const_atom
    LOGICAL                                   :: const_hydr
    LOGICAL                                   :: const_dist
    INTEGER          :: nconst_dist
    INTEGER, POINTER :: const_dist_mol(:)
    INTEGER, POINTER :: const_dist_a(:)
    INTEGER, POINTER :: const_dist_b(:)
    REAL (dbl), POINTER :: const_dist_dab(:)
    LOGICAL                                   :: const_33
    INTEGER          :: nconst_g33
    INTEGER, POINTER :: const_g33_mol(:)
    INTEGER, POINTER :: const_g33_a(:)
    INTEGER, POINTER :: const_g33_b(:)
    INTEGER, POINTER :: const_g33_c(:)
    REAL (dbl), POINTER :: const_g33_dab(:)
    REAL (dbl), POINTER :: const_g33_dac(:)
    REAL (dbl), POINTER :: const_g33_dbc(:)

  END TYPE topology_parameters_type

  PUBLIC :: topology_parameters_type,init_topology,deallocate_topology

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

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  NULLIFY(topology%cell)
  topology%natoms=0
  topology%natom_type=0
  topology%nmol=0
  topology%nmol_type=0
  topology%dump_topology=.FALSE.
  topology%reorder_atom=.FALSE.
  topology%coordinate=.FALSE.
  topology%connectivity=.FALSE.
  topology%nhcopt='GLOBAL'
  topology%bondparm_type='COVALENT'
  topology%bondparm_factor=1.1
  topology%coord_type=''
  topology%coord_file_name=''
  topology%conn_type='OFF'
  topology%conn_file_name='OFF'
  NULLIFY(topology%label_molname)
  NULLIFY(topology%label_resname)
  NULLIFY(topology%label_atmname)
  NULLIFY(topology%atom_names)
  NULLIFY(topology%r)
  NULLIFY(topology%map_mol_typ)
  NULLIFY(topology%map_mol_num)
  NULLIFY(topology%bond_a)
  NULLIFY(topology%bond_b)
  NULLIFY(topology%theta_a)
  NULLIFY(topology%theta_b)
  NULLIFY(topology%theta_c)
  NULLIFY(topology%phi_a)
  NULLIFY(topology%phi_b)
  NULLIFY(topology%phi_c)
  NULLIFY(topology%phi_d)
  NULLIFY(topology%atm_charge)
  NULLIFY(topology%atm_mass)
  NULLIFY(topology%occup)
  NULLIFY(topology%beta)

  topology%constraint=.FALSE.
  topology%const_atom=.FALSE.
  topology%const_hydr=.FALSE.
  topology%const_dist=.FALSE.
  topology%nconst_dist = 0
  NULLIFY(topology%const_dist_mol)
  NULLIFY(topology%const_dist_a)
  NULLIFY(topology%const_dist_b)
  NULLIFY(topology%const_dist_dab)
  topology%const_33=.FALSE.
  topology%nconst_g33 = 0
  NULLIFY(topology%const_g33_mol)
  NULLIFY(topology%const_g33_a)
  NULLIFY(topology%const_g33_b)
  NULLIFY(topology%const_g33_c)
  NULLIFY(topology%const_g33_dab)
  NULLIFY(topology%const_g33_dac)
  NULLIFY(topology%const_g33_dbc)
    
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

  IF(ASSOCIATED(topology%cell)) THEN
    DEALLOCATE(topology%cell,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','cell_type')
  END IF

  IF(ASSOCIATED(topology%label_molname)) THEN
    DEALLOCATE(topology%label_molname,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','label_molname')
  END IF

  IF(ASSOCIATED(topology%label_resname)) THEN
    DEALLOCATE(topology%label_resname,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','label_resname')
  END IF

  IF(ASSOCIATED(topology%label_atmname)) THEN
    DEALLOCATE(topology%label_atmname,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','label_atmname')
  END IF

  IF(ASSOCIATED(topology%atom_names)) THEN
    DEALLOCATE(topology%atom_names,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','atom_names')
  END IF

  IF(ASSOCIATED(topology%r)) THEN
    DEALLOCATE(topology%r,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','r')
  END IF

  IF(ASSOCIATED(topology%map_mol_typ)) THEN
    DEALLOCATE(topology%map_mol_typ,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','map_mol_typ')
  END IF

  IF(ASSOCIATED(topology%map_mol_num)) THEN
    DEALLOCATE(topology%map_mol_num,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','map_mol_num')
  END IF

  IF(ASSOCIATED(topology%bond_a)) THEN
    DEALLOCATE(topology%bond_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','bond_a')
  END IF

  IF(ASSOCIATED(topology%bond_b)) THEN
    DEALLOCATE(topology%bond_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','bond_b')
  END IF

  IF(ASSOCIATED(topology%theta_a)) THEN
    DEALLOCATE(topology%theta_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_a')
  END IF

  IF(ASSOCIATED(topology%theta_b)) THEN
    DEALLOCATE(topology%theta_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_b')
  END IF

  IF(ASSOCIATED(topology%theta_c)) THEN
    DEALLOCATE(topology%theta_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_c')
  END IF

  IF(ASSOCIATED(topology%phi_a)) THEN
    DEALLOCATE(topology%phi_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_a')
  END IF

  IF(ASSOCIATED(topology%phi_b)) THEN
    DEALLOCATE(topology%phi_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_b')
  END IF

  IF(ASSOCIATED(topology%phi_c)) THEN
    DEALLOCATE(topology%phi_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_c')
  END IF

  IF(ASSOCIATED(topology%phi_d)) THEN
    DEALLOCATE(topology%phi_d,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_d')
  END IF

  IF(ASSOCIATED(topology%atm_charge)) THEN
    DEALLOCATE(topology%atm_charge,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','atm_charge')
  END IF

  IF(ASSOCIATED(topology%atm_mass)) THEN
    DEALLOCATE(topology%atm_mass,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','atm_mass')
  END IF

  IF(ASSOCIATED(topology%occup)) THEN
    DEALLOCATE(topology%occup,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','occup')
  END IF

  IF(ASSOCIATED(topology%beta)) THEN
    DEALLOCATE(topology%beta,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','beta')
  END IF

  IF(ASSOCIATED(topology%const_dist_a)) THEN
    DEALLOCATE(topology%const_dist_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_dist_a')
  END IF

  IF(ASSOCIATED(topology%const_dist_b)) THEN
    DEALLOCATE(topology%const_dist_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_dist_b')
  END IF

  IF(ASSOCIATED(topology%const_dist_mol)) THEN
    DEALLOCATE(topology%const_dist_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_dist_mol')
  END IF

  IF(ASSOCIATED(topology%const_g33_mol)) THEN
    DEALLOCATE(topology%const_g33_mol,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_mol')
  END IF

  IF(ASSOCIATED(topology%const_g33_a)) THEN
    DEALLOCATE(topology%const_g33_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_a')
  END IF

  IF(ASSOCIATED(topology%const_g33_b)) THEN
    DEALLOCATE(topology%const_g33_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_b')
  END IF

  IF(ASSOCIATED(topology%const_g33_c)) THEN
    DEALLOCATE(topology%const_g33_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_c')
  END IF

  IF(ASSOCIATED(topology%const_g33_dab)) THEN
    DEALLOCATE(topology%const_g33_dab,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_dab')
  END IF

  IF(ASSOCIATED(topology%const_g33_dac)) THEN
    DEALLOCATE(topology%const_g33_dac,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_dac')
  END IF

  IF(ASSOCIATED(topology%const_g33_dbc)) THEN
    DEALLOCATE(topology%const_g33_dbc,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','const_g33_dbc')
  END IF


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE deallocate_topology



END MODULE topology_types

!!*****
!******************************************************************************
