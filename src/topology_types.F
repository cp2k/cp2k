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
    INTEGER                                   :: nmol, natoms, natom_type
    LOGICAL                                   :: coordinate
    LOGICAL                                   :: connectivity
    CHARACTER (LEN=default_string_length)     :: coord_type
    CHARACTER (LEN=default_string_length)     :: coord_file_name
    CHARACTER (LEN=default_string_length)     :: conn_type
    CHARACTER (LEN=default_string_length)     :: conn_file_name
    CHARACTER (LEN=default_string_length),&
               DIMENSION ( : ), POINTER       :: atom_names
    CHARACTER (LEN=20), POINTER :: label_resname(:)
    CHARACTER (LEN=20), POINTER :: label_atmname(:)
    REAL(dbl), DIMENSION(:,:), POINTER :: r
    INTEGER, POINTER :: map_mol_typ(:)
    INTEGER, POINTER :: map_mol_num(:)
    INTEGER, POINTER :: bond_a(:),bond_b(:)
    INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)
    INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)
    REAL (dbl), POINTER :: atm_charge(:)
    REAL (dbl), POINTER :: atm_mass(:)
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
  IMPLICIT NONE
  ! Arguments
  TYPE ( topology_parameters_type ),INTENT (INOUT) :: topology
  ! Locals

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  NULLIFY(topology%cell)
  topology%nmol=0
  topology%natoms=0
  topology%natom_type=0
  topology%coordinate=.FALSE.
  topology%connectivity=.FALSE.
  topology%coord_type=''
  topology%coord_file_name=''
  topology%conn_type=''
  topology%conn_file_name=''
  NULLIFY(topology%atom_names)
  NULLIFY(topology%label_resname)
  NULLIFY(topology%label_atmname)
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
  IMPLICIT NONE
  ! Arguments
  TYPE ( topology_parameters_type ), INTENT (INOUT) :: topology
  ! Locals
  INTEGER :: istat

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

  DEALLOCATE(topology%cell,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','cell_type')

  IF(ASSOCIATED(topology%atom_names)) THEN
    DEALLOCATE(topology%atom_names,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','atom_names')
  END IF

  DEALLOCATE(topology%label_resname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','label_resname')

  DEALLOCATE(topology%label_atmname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','label_atmname')

  DEALLOCATE(topology%r,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','r')

  DEALLOCATE(topology%map_mol_typ,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','map_mol_typ')

  DEALLOCATE(topology%map_mol_num,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','map_mol_num')

  DEALLOCATE(topology%bond_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','bond_a')

  DEALLOCATE(topology%bond_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','bond_b')

  DEALLOCATE(topology%theta_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_a')

  DEALLOCATE(topology%theta_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_b')

  DEALLOCATE(topology%theta_c,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','theta_c')

  DEALLOCATE(topology%phi_a,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_a')

  DEALLOCATE(topology%phi_b,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_b')

  DEALLOCATE(topology%phi_c,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_c')

  DEALLOCATE(topology%phi_d,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','phi_d')

  DEALLOCATE(topology%atm_charge,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','atm_charge')

  DEALLOCATE(topology%atm_mass,STAT=istat)
  IF (istat /= 0) CALL stop_memory('deallocate_topology','atm_mass')

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE deallocate_topology



END MODULE topology_types

!!*****
!******************************************************************************
