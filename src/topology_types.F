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
    REAL(dbl), DIMENSION(:,:), POINTER :: r
    INTEGER, POINTER :: map_mol_typ(:)
    INTEGER, POINTER :: map_mol_num(:)
    REAL (dbl), POINTER :: atm_charge(:)
    REAL (dbl), POINTER :: atm_mass(:)
    REAL (dbl), POINTER :: occup(:)
    REAL (dbl), POINTER :: beta(:)
    CHARACTER (LEN=default_string_length), POINTER :: element(:)
  END TYPE atom_info_type

  TYPE connectivity_info_type
    INTEGER, POINTER :: bond_a(:),bond_b(:)
    INTEGER, POINTER :: theta_a(:),theta_b(:),theta_c(:)
    INTEGER, POINTER :: phi_a(:),phi_b(:),phi_c(:),phi_d(:)
  END TYPE connectivity_info_type

  TYPE constraint_info_type
    INTEGER          :: nconst_dist
    INTEGER, POINTER :: const_dist_mol(:)
    INTEGER, POINTER :: const_dist_a(:)
    INTEGER, POINTER :: const_dist_b(:)
    REAL (dbl), POINTER :: const_dist_dab(:)
    INTEGER          :: nconst_g33
    INTEGER, POINTER :: const_g33_mol(:)
    INTEGER, POINTER :: const_g33_a(:)
    INTEGER, POINTER :: const_g33_b(:)
    INTEGER, POINTER :: const_g33_c(:)
    REAL (dbl), POINTER :: const_g33_dab(:)
    REAL (dbl), POINTER :: const_g33_dac(:)
    REAL (dbl), POINTER :: const_g33_dbc(:)
  END TYPE constraint_info_type

  TYPE force_field_input_info_type
    CHARACTER (LEN=default_string_length), POINTER :: ff_charge_atm(:)
    REAL (dbl), POINTER                            :: ff_charge(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_bond_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_bond_b(:)
    REAL (dbl), POINTER                            :: ff_bond_k(:)
    REAL (dbl), POINTER                            :: ff_bond_r0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_bend_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_bend_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_bend_c(:)
    REAL (dbl), POINTER                            :: ff_bend_k(:)
    REAL (dbl), POINTER                            :: ff_bend_theta0(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_torsion_a(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_torsion_b(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_torsion_c(:)
    CHARACTER (LEN=default_string_length), POINTER :: ff_torsion_d(:)
    REAL (dbl), POINTER                            :: ff_torsion_k(:)
    REAL (dbl), POINTER                            :: ff_torsion_phi0(:)
    TYPE (pair_potential_type), POINTER            :: ff_nonbonded(:)
  END TYPE force_field_input_info_type

  TYPE force_field_charmm_info_type
    INTEGER                                   :: nothing_here_yet
  END TYPE force_field_charmm_info_type

  TYPE force_field_amber_info_type
    INTEGER                                   :: nothing_here_yet
  END TYPE force_field_amber_info_type

  TYPE topology_parameters_type
    TYPE (atom_info_type),POINTER                :: atom_info
    TYPE (connectivity_info_type),POINTER        :: conn_info
    TYPE (constraint_info_type),POINTER          :: cons_info
    TYPE (force_field_input_info_type),POINTER   :: ffin_info
    TYPE (force_field_charmm_info_type),POINTER  :: ffch_info
    TYPE (force_field_charmm_info_type),POINTER  :: ffam_info
    TYPE ( cell_type ), POINTER                  :: cell
    !TRY TO REMOVE THIS FOUR VARIABLE IN THE FUTURE
    INTEGER                                   :: natoms, natom_type
    INTEGER                                   :: nmol, nmol_type
    !TRY TO REMOVE THIS FOUR VARIABLE IN THE FUTURE
    CHARACTER (LEN=default_string_length)     :: bondparm_type
    REAL (dbl)                                :: bondparm_factor
    LOGICAL                                   :: reorder_atom
    LOGICAL                                   :: coordinate
    LOGICAL                                   :: connectivity
    LOGICAL                                   :: dump_topology
    LOGICAL                                   :: ff
    CHARACTER (LEN=default_string_length)     :: ff_file_name
    CHARACTER (LEN=default_string_length)     :: ff_type
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
  END TYPE topology_parameters_type


  PUBLIC :: atom_info_type,&
            connectivity_info_type,&
            constraint_info_type,&
            force_field_input_info_type,&
            force_field_charmm_info_type,&
            force_field_amber_info_type,&
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
  IF(istat/=0) CALL stop_memory &
    ('init_topology','atom_info')
  ALLOCATE(topology%conn_info,STAT=istat)
  IF(istat/=0) CALL stop_memory &
    ('init_topology','conn_info')
  ALLOCATE(topology%cons_info,STAT=istat)
  IF(istat/=0) CALL stop_memory &
    ('init_topology','cons_info')
  ALLOCATE(topology%ffin_info,STAT=istat)
  IF(istat/=0) CALL stop_memory &
    ('init_topology','ffin_info')
  ALLOCATE(topology%ffch_info,STAT=istat)
  IF(istat/=0) CALL stop_memory &
    ('init_topology','ffch_info')
  ALLOCATE(topology%ffam_info,STAT=istat)
  IF(istat/=0) CALL stop_memory &
    ('init_topology','ffam_info')


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. Initialize and Nullify things in topology
  !-----------------------------------------------------------------------------
  NULLIFY(topology%cell)
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
  topology%ff=.FALSE.
  topology%ff_file_name=''
  topology%ff_type='OFF'
  topology%nhcopt='GLOBAL'
  topology%coord_type=''
  topology%coord_file_name=''
  topology%conn_type='OFF'
  topology%conn_file_name='OFF'
  topology%constraint=.FALSE.
  topology%const_atom=.FALSE.
  topology%const_hydr=.FALSE.
  topology%const_dist=.FALSE.
  topology%const_33=.FALSE.


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
  NULLIFY(topology%conn_info%theta_a)
  NULLIFY(topology%conn_info%theta_b)
  NULLIFY(topology%conn_info%theta_c)
  NULLIFY(topology%conn_info%phi_a)
  NULLIFY(topology%conn_info%phi_b)
  NULLIFY(topology%conn_info%phi_c)
  NULLIFY(topology%conn_info%phi_d)


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


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. Initialize and Nullify things in topology%ffin_info
  !-----------------------------------------------------------------------------
  NULLIFY(topology%ffin_info%ff_charge_atm)
  NULLIFY(topology%ffin_info%ff_charge)
  NULLIFY(topology%ffin_info%ff_bond_a)
  NULLIFY(topology%ffin_info%ff_bond_b)
  NULLIFY(topology%ffin_info%ff_bond_k)
  NULLIFY(topology%ffin_info%ff_bond_r0)
  NULLIFY(topology%ffin_info%ff_bend_a)
  NULLIFY(topology%ffin_info%ff_bend_b)
  NULLIFY(topology%ffin_info%ff_bend_c)
  NULLIFY(topology%ffin_info%ff_bend_k)
  NULLIFY(topology%ffin_info%ff_bend_theta0)
  NULLIFY(topology%ffin_info%ff_torsion_a)
  NULLIFY(topology%ffin_info%ff_torsion_b)
  NULLIFY(topology%ffin_info%ff_torsion_c)
  NULLIFY(topology%ffin_info%ff_torsion_d)
  NULLIFY(topology%ffin_info%ff_torsion_k)
  NULLIFY(topology%ffin_info%ff_torsion_phi0)
  NULLIFY(topology%ffin_info%ff_nonbonded)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. Initialize and Nullify things in topology%ffch_info
  !-----------------------------------------------------------------------------
    


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 8. Initialize and Nullify things in topology%ffam_info
  !-----------------------------------------------------------------------------



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


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. DEALLOCATE things in topology%ffin_info
  !-----------------------------------------------------------------------------
  IF(ASSOCIATED(topology%ffin_info%ff_charge_atm)) THEN
    DEALLOCATE(topology%ffin_info%ff_charge_atm,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_charge_atm')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_charge)) THEN
    DEALLOCATE(topology%ffin_info%ff_charge,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_charge')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bond_a)) THEN
    DEALLOCATE(topology%ffin_info%ff_bond_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bond_a')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bond_b)) THEN
    DEALLOCATE(topology%ffin_info%ff_bond_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bond_b')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bond_k)) THEN
    DEALLOCATE(topology%ffin_info%ff_bond_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bond_k')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bond_r0)) THEN
    DEALLOCATE(topology%ffin_info%ff_bond_r0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bond_r0')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bend_a)) THEN
    DEALLOCATE(topology%ffin_info%ff_bend_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bend_a')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bend_b)) THEN
    DEALLOCATE(topology%ffin_info%ff_bend_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bend_b')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bend_c)) THEN
    DEALLOCATE(topology%ffin_info%ff_bend_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bend_c')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bend_k)) THEN
    DEALLOCATE(topology%ffin_info%ff_bend_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bend_k')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_bend_theta0)) THEN
    DEALLOCATE(topology%ffin_info%ff_bend_theta0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_bend_theta0')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_torsion_a)) THEN
    DEALLOCATE(topology%ffin_info%ff_torsion_a,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_torsion_a')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_torsion_b)) THEN
    DEALLOCATE(topology%ffin_info%ff_torsion_b,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_torsion_b')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_torsion_c)) THEN
    DEALLOCATE(topology%ffin_info%ff_torsion_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_torsion_c')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_torsion_d)) THEN
    DEALLOCATE(topology%ffin_info%ff_torsion_d,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_torsion_d')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_torsion_k)) THEN
    DEALLOCATE(topology%ffin_info%ff_torsion_k,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_torsion_k')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_torsion_phi0)) THEN
    DEALLOCATE(topology%ffin_info%ff_torsion_phi0,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_torsion_phi0')
  END IF
  IF(ASSOCIATED(topology%ffin_info%ff_nonbonded)) THEN
    DEALLOCATE(topology%ffin_info%ff_nonbonded,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ff_nonbonded')
  END IF


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 5. DEALLOCATE things in topology%ffch_info
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 6. DEALLOCATE things in topology%ffam_info
  !-----------------------------------------------------------------------------



  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 7. DEALLOCATE things in topology
  !-----------------------------------------------------------------------------
  CALL cell_release(topology%cell)
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
  IF(ASSOCIATED(topology%ffin_info)) THEN
    DEALLOCATE(topology%ffin_info,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ffin_info')
  END IF
  IF(ASSOCIATED(topology%ffch_info)) THEN
    DEALLOCATE(topology%ffch_info,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ffch_info')
  END IF
  IF(ASSOCIATED(topology%ffam_info)) THEN
    DEALLOCATE(topology%ffam_info,STAT=istat)
    IF (istat /= 0) CALL stop_memory('deallocate_topology','ffam_info')
  END IF



!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE deallocate_topology



END MODULE topology_types

!!*****
!******************************************************************************
