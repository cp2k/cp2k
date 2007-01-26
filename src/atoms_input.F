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

  USE atomic_kind_types,               ONLY: atomic_kind_type, &
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell, pbc, &
                                             scaled_to_real
  USE checkpoint_handler,              ONLY: Pwarning
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_type
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_atoms_input, read_shell_coord_input
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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_atoms_input', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 10000

    CHARACTER(len=2*default_string_length)   :: line_att
    CHARACTER(len=default_string_length)     :: my_default_index
    INTEGER                                  :: end_c, handle, iatom, ILEN, &
                                                j, natom, output_unit, &
                                                start_c, wrd
    LOGICAL                                  :: explicit, failure, ionode, &
                                                is_ok, my_overwrite, &
                                                scaled_coordinates
    REAL(KIND=dp)                            :: unit_of_length
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(section_vals_type), POINTER         :: coord_section
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

    CALL timeset(routineN,handle)
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

!******************************************************************************
!!****** atoms_input/read_shell_coord_input [1.0] *
!!
!!   NAME
!!     read_shell_coord_input
!!
!!   FUNCTION
!!
!!   AUTHOR
!!    MI
!!
!!   MODIFICATION HISTORY
!!
!!  OPTIONS
!!
!!  INPUTS
!!*** *************************************************************************

  SUBROUTINE read_shell_coord_input ( particle_set, shell_particle_set, cell, root_section, subsys_section, error )

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: shell_particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_shell_coord_input', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nblock_atom = 10000

    CHARACTER(len=2*default_string_length)   :: line_att
    CHARACTER(len=default_string_length)     :: my_default_index
    CHARACTER(len=default_string_length)     :: name_kind
    CHARACTER(len=default_string_length), &
      DIMENSION(:), ALLOCATABLE              :: at_name
    INTEGER                                  :: end_c, handle, ishell, ILEN, &
                                                istat, j, nshell, output_unit, &
                                                start_c, sh_index, wrd
    INTEGER, DIMENSION(:), ALLOCATABLE       :: at_index
    LOGICAL                                  :: explicit, failure, ionode, &
                                                is_ok, is_shell, &
                                                scaled_coordinates
    REAL(KIND=dp)                            :: unit_of_length
    REAL(KIND=dp)                            :: dab, rab(3)
    REAL(KIND=dp), DIMENSION(:,:), &
      ALLOCATABLE                            :: r
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(section_vals_type), POINTER         :: shell_coord_section
    TYPE(shell_type), POINTER                :: shell
    TYPE(val_type), POINTER                  :: val

    failure=.FALSE.
    NULLIFY(logger, atomic_kind, list, shell_coord_section, shell, val)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
      output_unit = -1
    END IF

    NULLIFY(shell_coord_section)
    shell_coord_section => section_vals_get_subs_vals(subsys_section,"SHELL_COORD",error=error)
    CALL section_vals_get(shell_coord_section, explicit=explicit, error=error)
    IF (.NOT.explicit) RETURN

    CALL timeset(routineN,handle)
    CALL section_vals_val_get(shell_coord_section,"_DEFAULT_KEYWORD_",&
         n_rep_val=nshell,error=error)

    CPPostcondition(SIZE(shell_particle_set,1)==nshell, cp_failure_level, routineP, error, failure)
    IF(.NOT. failure) THEN
      ALLOCATE (r(3,nshell),at_name(nshell),at_index(nshell), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      CALL Pwarning("Overwriting shells coordinates. "//&
           " Active coordinates READ from &SHELL_COORD section",&
           root_section, error)
      CALL section_vals_list_get(shell_coord_section,"_DEFAULT_KEYWORD_",list=list,error=error)
      DO ishell=1,nshell
         ! we use only the first default_string_length characters of each line
         is_ok=cp_sll_val_next(list,val,error=error)
         CALL val_get(val,c_val=line_att,error=error)
         !CALL section_vals_val_get(shell_coord_section,"_DEFAULT_KEYWORD_",&
         !     i_rep_val=iatom,c_val=line_att,error=error)
         ! *** Read name and atomic coordinates ***
         start_c=1
         DO wrd=1,5
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
            CALL cp_assert(wrd==5.OR.end_c<LEN(line_att)+1,cp_failure_level,&
                 cp_assertion_failed,routineP,&
                 "incorrectly formatted line in coord section'"//line_att//"'",&
                 error,failure)
            IF (wrd==1) THEN
               at_name(ishell)=line_att(start_c:end_c-1)
               CALL uppercase(at_name(ishell))
            ELSE IF (wrd==5) THEN
               READ (line_att(start_c:end_c-1),*) at_index(ishell)
            ELSE
               READ (line_att(start_c:end_c-1),*) r(wrd-1,ishell)
            END IF
            start_c=end_c
         END DO
      END DO
      !-----------------------------------------------------------------------------
      !-----------------------------------------------------------------------------
      ! 3. Check corrispondence and convert coordinates into internal cp2k coordinates
      !-----------------------------------------------------------------------------
      DO ishell=1,nshell
         atomic_kind => particle_set(at_index(ishell))%atomic_kind
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              name=name_kind,shell_active=is_shell,shell=shell)
         CALL uppercase(name_kind)
         IF((TRIM(at_name(ishell))==TRIM(name_kind)).AND. is_shell) THEN
           sh_index = particle_set(at_index(ishell))%shell_index
           IF (scaled_coordinates) THEN
             CALL scaled_to_real(r(:,ishell),shell_particle_set(sh_index)%r(:),cell)
           ELSE
             shell_particle_set(sh_index)%r(:) = r(:,ishell)*unit_of_length
           END IF

           rab = pbc(shell_particle_set(sh_index)%r,particle_set(at_index(ishell))%r,cell)
           dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
           IF(shell%max_dist>0.0_dp .AND. shell%max_dist<dab ) THEN
             IF(ionode) THEN
               WRITE(output_unit,*) "WARNING : shell and core for atom ", at_index(ishell), " seem to be too distant. "
             END IF
           END IF
 
         ELSE
           CALL cp_assert(.FALSE.,cp_failure_level,&
                 cp_assertion_failed,routineP,&
                 "shell coordinate assigned to the wrong atom. check the shell indexes in the input",&
                 error,failure)
         END IF
      END DO
      DEALLOCATE(r,at_index,at_name,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    
    END IF

    CALL timestop(handle)

  END SUBROUTINE read_shell_coord_input

END MODULE atoms_input
