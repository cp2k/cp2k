!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/atoms_input [1.0] *
!!
!!   NAME
!!     atoms_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!   cjm, Feb-20-2001 : added all the extended variables to
!!   system_type
!!   gt 23-09-2002 : major changes. Pointer part is allocated/deallocated
!!                   and initialized here. Atomic coordinates can now be 
!!                   read also from &COORD section in the input file.
!!                   If &COORD is not found, .dat file is read.  
!!                   If & coord is found and .NOT. 'INIT', parsing of the .dat
!!                   is performed to get the proper coords/vel/eta variables
!!    CJM 31-7-03  : Major rewrite.  No more atype
!!
!!   SOURCE
!******************************************************************************

MODULE atoms_input

  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             scaled_to_real
  USE checkpoint_handler,              ONLY: Pwarning
  USE cp_linked_list_val,              ONLY: cp_sll_val_type, cp_sll_val_next
  USE input_val_types,                 ONLY: val_type, val_get
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_list_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
  USE input_constants,                 ONLY: do_conn_off,&
                                             do_conn_generate
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_atoms_input
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atoms_input'
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** atoms_input/read_atoms_input [1.0] *
!!
!!   NAME
!!     read_atoms_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    CJM
!!
!!   MODIFICATION HISTORY
!!
!!  OPTIONS
!!
!!  INPUTS
!!*** *************************************************************************

  SUBROUTINE read_atoms_input ( topology, root_section, overwrite, subsys_section, error )


    TYPE(topology_parameters_type)           :: topology
    TYPE(section_vals_type), POINTER         :: root_section
    LOGICAL, INTENT(IN), OPTIONAL            :: overwrite
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_atoms_input', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 10000

    CHARACTER(len=2*default_string_length)   :: line_att
    CHARACTER(len=default_string_length)     :: my_default_index
    INTEGER                                  :: atm_equal_mol, end_c, iatom, &
                                                ILEN, j, natom, start_c, wrd,&
                                                output_unit, handle
    LOGICAL                                  :: explicit, failure, &
                                                my_overwrite, &
                                                scaled_coordinates, ionode, is_ok
    REAL(KIND=dp)                            :: unit_of_length
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: coord_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: list 
    TYPE(val_type), POINTER                  :: val


    failure=.FALSE.
    my_overwrite = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
      output_unit = -1
    END IF
    IF (PRESENT(overwrite)) my_overwrite=overwrite
    NULLIFY(coord_section)
    coord_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
    CALL section_vals_get(coord_section, explicit=explicit, error=error)
    IF (.NOT.explicit) RETURN

    call timeset(routineN,handle)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. get cell and topolgoy%atom_info
    !-----------------------------------------------------------------------------
    atom_info => topology%atom_info
    cell => topology%cell
    CALL get_cell(cell=cell,&
         scaled_coordinates=scaled_coordinates,&
         unit_of_length=unit_of_length)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Read in the coordinates from &COORD section in the input file
    !-----------------------------------------------------------------------------
    CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
         n_rep_val=natom,error=error)
    topology%natoms = natom
    IF (my_overwrite) THEN
       CPPostcondition(SIZE(atom_info%r,2)==natom, cp_failure_level, routineP, error, failure)
       CALL Pwarning("Overwriting coordinates. Active coordinates read from &COORD section."//&
            " Active coordinates READ from &COORD section",&
            root_section, error)
       CALL section_vals_list_get(coord_section,"_DEFAULT_KEYWORD_",list=list,error=error)
       DO iatom=1,natom
          ! we use only the first default_string_length characters of each line
          is_ok=cp_sll_val_next(list,val,error=error)
          CALL val_get(val,c_val=line_att,error=error)
          !CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
          !     i_rep_val=iatom,c_val=line_att,error=error)
          ! *** Read name and atomic coordinates ***
          start_c=1
          DO wrd=1,4
             DO j=start_c,LEN(line_att)
                IF (line_att(j:j)/=' ') THEN
                   start_c=j
                   EXIT
                END IF
             END DO
             end_c=LEN(line_att)+1
             DO j=start_c,LEN(line_att)
                IF (line_att(j:j)==' ') THEN
                   end_c=j
                   EXIT
                END IF
             END DO
             CALL cp_assert(wrd==4.OR.end_c<LEN(line_att)+1,cp_failure_level,&
                  cp_assertion_failed,routineP,&
                  "incorrectly formatted line in coord section'"//line_att//"'",&
                  error,failure)
             IF (wrd==1) THEN
                atom_info%label_atmname(iatom)=line_att(start_c:end_c-1)
             ELSE
                READ (line_att(start_c:end_c-1),*) atom_info%r(wrd-1,iatom)
             END IF
             start_c=end_c
          END DO
       END DO
    ELSE
       CALL reallocate(atom_info%label_molname,1,natom)
       CALL reallocate(atom_info%label_resname,1,natom)
       CALL reallocate(atom_info%label_resid,1,natom)
       CALL reallocate(atom_info%label_atmname,1,natom)
       CALL reallocate(atom_info%element,1,natom)
       CALL reallocate(atom_info%r,1,3,1,natom)
       CALL reallocate(atom_info%atm_mass,1,natom)
       CALL reallocate(atom_info%atm_charge,1,natom)

       CALL section_vals_list_get(coord_section,"_DEFAULT_KEYWORD_",list=list,error=error)
       DO iatom=1,natom
          ! we use only the first default_string_length characters of each line
          is_ok=cp_sll_val_next(list,val,error=error)
          CALL val_get(val,c_val=line_att,error=error)
          !CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
          !     i_rep_val=iatom,c_val=line_att,error=error)
          ILEN = default_string_length
          atom_info%label_molname(iatom) = ""
          atom_info%label_resname(iatom) = ""
          atom_info%label_resid(iatom) = 1
          atom_info%label_atmname(iatom) = ""
          atom_info%element(iatom) = ""
          ! *** Read name and atomic coordinates ***
          start_c=1
          DO wrd=1,4
             DO j=start_c,LEN(line_att)
                IF (line_att(j:j)/=' ') THEN
                   start_c=j
                   EXIT
                END IF
             END DO
             end_c=LEN(line_att)+1
             DO j=start_c,LEN(line_att)
                IF (line_att(j:j)==' ') THEN
                   end_c=j
                   EXIT
                END IF
             END DO
             CALL cp_assert(wrd==4.OR.end_c<LEN(line_att)+1,cp_failure_level,&
                  cp_assertion_failed,routineP,&
                  "incorrectly formatted line in coord section'"//line_att//"'",&
                  error,failure)
             IF (wrd==1) THEN
                atom_info%label_atmname(iatom)=line_att(start_c:end_c-1)
             ELSE
                READ (line_att(start_c:end_c-1),*) atom_info%r(wrd-1,iatom)
             END IF
             start_c=end_c
          END DO
          ! For default, set atom name to residue name to molecule name
          WRITE(my_default_index,'(I0)')iatom
          atom_info%label_molname(iatom) = TRIM(atom_info%label_atmname(iatom))//TRIM(my_default_index)
          atom_info%label_resname(iatom) = TRIM(atom_info%label_molname(iatom))
          atom_info%element(iatom)       = TRIM(atom_info%label_atmname(iatom))
          topology%molname_generated = .TRUE.
          IF (start_c<=LEN_TRIM(line_att)) THEN
             READ(line_att(start_c:LEN_TRIM(line_att)),*)atom_info%label_molname(iatom)
             atom_info%label_resname(iatom) = atom_info%label_molname(iatom)
             topology%molname_generated = .FALSE.
          END IF
          atom_info%atm_mass(iatom)   = 0.0_dp 
          atom_info%atm_charge(iatom) = -HUGE(0.0_dp)
       END DO
    END IF
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Convert coordinates into internal cp2k coordinates
    !-----------------------------------------------------------------------------
    DO iatom=1,natom
       IF (scaled_coordinates) THEN
          CALL scaled_to_real(atom_info%r(:,iatom),atom_info%r(:,iatom),cell)
       ELSE
          atom_info%r(:,iatom) = atom_info%r(:,iatom)*unit_of_length
       END IF
    END DO
    CALL timestop(handle)
  END SUBROUTINE read_atoms_input

END MODULE atoms_input
