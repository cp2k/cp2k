!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/topology_pdb [1.0] *
!!
!!   NAME
!!     topology_pdb
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     IKUO   08.16.03
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
MODULE topology_pdb

  USE atomic_kind_types,        ONLY : allocate_atomic_kind_set,&
                                       atomic_kind_type,&
                                       get_atomic_kind,&
                                       init_atomic_kind_set,&
                                       set_atomic_kind
  USE convert_units,            ONLY : convert_to_cp2k_units
  USE external_potential_types, ONLY : epc_potential_type,&
                                       set_potential
  USE global_types,             ONLY : global_environment_type
  USE input_types,              ONLY : setup_parameters_type
  USE kinds,                    ONLY : dbl
  USE memory_utilities,         ONLY : reallocate
  USE molecule_types,           ONLY : molecule_type
  USE extended_system_dynamics, ONLY : lnhc_parameters_type, &
                                       npt_info_type
  USE particle_types,           ONLY : particle_type, &
                                       allocate_particle_set
  USE parser,                   ONLY : parser_init, &
                                       parser_end, &
                                       read_line, &
                                       test_next, &
                                       p_error, &
                                       get_next, &
                                       search_label
  USE periodic_table,           ONLY : ptable
  USE physcon,                  ONLY : massunit
  USE simulation_cell,          ONLY : cell_type
  USE string_utilities,         ONLY : xstring,&
                                       integer_to_string, &
                                       string_to_integer, &
                                       str_search, &
                                       uppercase
  USE termination,              ONLY : stop_memory, &
                                       stop_program
  USE topology_util,            ONLY : azero,&
                                       topology_coord_pack

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_coordinate_pdb

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** topology_pdb/read_coordinate_pdb [1.0] *
!!
!!   NAME
!!     read_coordinate_pdb
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     IKUO   08.16.03
!!
!!   MODIFICATION HISTORY
!!
!!   OPTIONS
!!
!!   INPUTS
!!*** *************************************************************************

SUBROUTINE read_coordinate_pdb (r,label_resname,label_atmname,setup,globenv)

  IMPLICIT NONE

! Arguments
  REAL(dbl), DIMENSION(:,:), POINTER :: r
  CHARACTER (LEN=20), POINTER :: label_resname(:),label_atmname(:)
  TYPE ( setup_parameters_type ) :: setup
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: ierror, ilen, iw,  ia, ie, i
  INTEGER, PARAMETER :: max_name_length = 20
  INTEGER, PARAMETER :: nblock_atom = 1000,&
                        nblock_kind = 10

  CHARACTER(LEN=max_name_length)  :: current_label,label,string
  CHARACTER(LEN=8)                :: unit_name
  CHARACTER(LEN=3)                :: test_result
  INTEGER                         :: iatom,ikind,istat,natom,nkind,z,n,iat,&
                                     natom_types,nmol_type,k,j
  LOGICAL                         :: new_kind,coords,&
                                       init

  INTEGER, POINTER :: num_atmid(:),num_resid(:)
  REAL (dbl), POINTER :: pos1(:),pos2(:),pos3(:)
  CHARACTER (LEN=20), POINTER :: tmp_resname(:),tmp_atmname(:)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  iw = globenv % scr
  natom = 0
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Change to REALLOCATE later
  NULLIFY(num_atmid,num_resid,pos1,pos2,pos3)
  NULLIFY(tmp_resname,tmp_atmname)
  ALLOCATE(num_atmid(100000),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','num_atmid')
  ALLOCATE(num_resid(100000),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','num_resid')
  ALLOCATE(tmp_resname(100000),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','tmp_resname')
  ALLOCATE(tmp_atmname(100000),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','tmp_atmname')
  ALLOCATE(pos1(100000),pos2(100000),pos3(100000),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','pos')
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Parse the file right now and put all info into temporary storage
  WRITE(iw,*) "READING IN PDB FILE ",setup%coord_file_name
  natom = 0
  label = 'REMARK PDB'
  CALL parser_init(setup%coord_file_name,globenv )
  CALL search_label ( label, ierror, ignore_case=.TRUE. )
  IF ( ierror /= 0 ) THEN
       CALL stop_program ( 'read_coordinate_pdb', 'REMARK PDB'  )
  ELSE
    CALL read_line
    DO WHILE (test_next() /= 'X')
      ilen = 10
      CALL get_next(string,ilen)
      CALL uppercase(string)
      SELECT CASE (string)
      CASE ('ATOM')
        natom = natom + 1
        CALL get_next(num_atmid(natom))
        CALL get_next(tmp_atmname(natom),ilen)
        CALL get_next(tmp_resname(natom),ilen)
        CALL get_next(num_resid(natom))
        CALL get_next(pos1(natom))
        CALL get_next(pos2(natom))
        CALL get_next(pos3(natom))
      CASE ('HEATOM')
        CALL stop_program ('read_coordinate_pdb','HEATOM not implemented yet')
      CASE ('END')
        EXIT
      CASE DEFAULT
      END SELECT
      CALL read_line
    END DO
  END IF
  CALL parser_end
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  !Put the info into array that got passed in
  ALLOCATE(r(3,natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','r(3,natom)')
  ALLOCATE(label_resname(natom),label_atmname(natom),STAT=istat)
  DO iatom=1,natom
    r(1,iatom) = pos1(iatom)
    r(2,iatom) = pos2(iatom)
    r(3,iatom) = pos3(iatom)
    CALL convert_to_cp2k_units("ANGSTROM",length=r(1,iatom))
    CALL convert_to_cp2k_units("ANGSTROM",length=r(2,iatom))
    CALL convert_to_cp2k_units("ANGSTROM",length=r(3,iatom))
    label_resname(iatom)=tmp_resname(iatom)
    label_atmname(iatom)=tmp_atmname(iatom)
  END DO
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  !Deallocate all the temporary arrays
  DEALLOCATE(num_atmid,num_resid,STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','num_atmid,num_resid')
  DEALLOCATE(pos1,pos2,pos3,STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','pos1,pos2,pos3')
  DEALLOCATE(tmp_resname,tmp_atmname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','tmp_[res,atm]name')
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE read_coordinate_pdb

END MODULE topology_pdb
