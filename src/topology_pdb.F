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

  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE memory_utilities,                ONLY: reallocate
  USE parser,                          ONLY: get_next,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             test_next
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: topology_parameters_type

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "topology_pdb"

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

SUBROUTINE read_coordinate_pdb (topology,globenv)


    TYPE(topology_parameters_type)           :: topology
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(len=*), PARAMETER :: routineN = 'read_coordinate_pdb', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 1000

    CHARACTER(LEN=default_string_length)     :: label, string
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: tmp_atmname(:), tmp_resname(:)
    INTEGER                                  :: handle, iatom, ilen, istat, &
                                                iw, natom
    INTEGER, POINTER                         :: num_atmid(:), num_resid(:)
    REAL(dbl), POINTER                       :: pos1(:), pos2(:), pos3(:)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  iw = globenv % scr

  CALL write_checkpoint_information("entering "//routineN,globenv)
  CALL timeset(routineN,'I','',handle)
  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Entering read_coordinate_pdb"
  END IF
  natom = 0
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Change to REALLOCATE later
  NULLIFY(num_atmid,num_resid,pos1,pos2,pos3)
  NULLIFY(tmp_resname,tmp_atmname)
  ALLOCATE(num_atmid(nblock_atom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','num_atmid')
  ALLOCATE(num_resid(nblock_atom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','num_resid')
  ALLOCATE(tmp_resname(nblock_atom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','tmp_resname')
  ALLOCATE(tmp_atmname(nblock_atom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','tmp_atmname')
  ALLOCATE(pos1(nblock_atom),pos2(nblock_atom),pos3(nblock_atom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','pos')

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!Parse the file right now and put all info into temporary storage
  !WRITE(iw,*) "READING IN PDB FILE ",topology%coord_file_name
  natom = 0
  label = 'REMARK PDB'
  CALL parser_init(topology%coord_file_name,globenv )
  CALL search_label ( label, istat, ignore_case=.TRUE. )
  IF ( istat /= 0 ) THEN
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
        IF(natom > SIZE(pos1)) THEN
          CALL reallocate(num_atmid,1,natom+nblock_atom)
          CALL reallocate(num_resid,1,natom+nblock_atom)
          CALL reallocate(tmp_resname,1,natom+nblock_atom)
          CALL reallocate(tmp_atmname,1,natom+nblock_atom)
          CALL reallocate(pos1,1,natom+nblock_atom)
          CALL reallocate(pos2,1,natom+nblock_atom)
          CALL reallocate(pos3,1,natom+nblock_atom)
        END IF
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
  ALLOCATE(topology%r(3,natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','r(3,natom)')
  ALLOCATE(topology%label_resname(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','label_resname')
  ALLOCATE(topology%label_atmname(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','label_atmname')
  ALLOCATE(topology%atm_mass(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','label_atmname')
  ALLOCATE(topology%atm_charge(natom),STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','label_atmname')
  DO iatom=1,natom
    topology%r(1,iatom) = pos1(iatom)
    topology%r(2,iatom) = pos2(iatom)
    topology%r(3,iatom) = pos3(iatom)
    CALL convert_to_cp2k_units("ANGSTROM",length=topology%r(1,iatom))
    CALL convert_to_cp2k_units("ANGSTROM",length=topology%r(2,iatom))
    CALL convert_to_cp2k_units("ANGSTROM",length=topology%r(3,iatom))
    topology%label_resname(iatom)=tmp_resname(iatom)
    topology%label_atmname(iatom)=tmp_atmname(iatom)
    topology%atm_mass(iatom)  =0.0_dbl
    topology%atm_charge(iatom)=0.0_dbl
  END DO

  topology%natoms = natom

!------------------------------------------------------------------------------
  !Deallocate all the temporary arrays
  DEALLOCATE(num_atmid,num_resid,STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','num_atmid,num_resid')
  DEALLOCATE(pos1,pos2,pos3,STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','pos1,pos2,pos3')
  DEALLOCATE(tmp_resname,tmp_atmname,STAT=istat)
  IF (istat /= 0) CALL stop_memory('read_coordinate_pdb','tmp_[res,atm]name')


  IF((globenv%ionode).AND.(globenv%print%level>0))THEN
    WRITE(iw,*) "  Exiting read_coordinate_pdb"
  END IF

  CALL timestop(0.0_dbl,handle)
  CALL write_checkpoint_information("exiting "//routineN,globenv)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END SUBROUTINE read_coordinate_pdb

END MODULE topology_pdb
