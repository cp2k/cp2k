!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>    cjm, Feb-20-2001 : added all the extended variables to
!>    system_type
!>    gt 23-09-2002 : major changes. Pointer part is allocated/deallocated
!>                    and initialized here. Atomic coordinates can now be
!>                    read also from &COORD section in the input file.
!>                    If &COORD is not found, .dat file is read.
!>                    If & coord is found and .NOT. 'INIT', parsing of the .dat
!>                    is performed to get the proper coords/vel/eta variables
!>     CJM 31-7-03  : Major rewrite.  No more atype
! *****************************************************************************
MODULE atoms_input
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc,&
                                             scaled_to_real
  USE cp_linked_list_val,              ONLY: cp_sll_val_next,&
                                             cp_sll_val_type
  USE cp_parser_methods,               ONLY: read_float_object
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_list_get,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: val_get,&
                                             val_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE particle_types,                  ONLY: particle_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_atoms_input, read_shell_coord_input
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atoms_input'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param topology ...
!> \param overwrite ...
!> \param subsys_section ...
!> \param save_mem ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE read_atoms_input ( topology, overwrite, subsys_section, save_mem, error )

    TYPE(topology_parameters_type)           :: topology
    LOGICAL, INTENT(IN), OPTIONAL            :: overwrite
    TYPE(section_vals_type), POINTER         :: subsys_section
    LOGICAL, INTENT(IN), OPTIONAL            :: save_mem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_atoms_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2*default_string_length)   :: line_att
    CHARACTER(len=default_string_length)     :: error_message, &
                                                my_default_index, strtmp, &
                                                unit_str
    INTEGER                                  :: default_id, end_c, handle, &
                                                iatom, j, natom, output_unit, &
                                                start_c, wrd
    LOGICAL                                  :: explicit, failure, is_ok, &
                                                my_overwrite, my_save_mem, &
                                                scaled_coordinates
    REAL(KIND=dp)                            :: r0(3), unit_conv
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(section_vals_type), POINTER         :: coord_section
    TYPE(val_type), POINTER                  :: val

    failure=.FALSE.
    my_overwrite = .FALSE.
    my_save_mem = .FALSE.
    error_message = ""
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    IF (PRESENT(overwrite)) my_overwrite=overwrite
    IF (PRESENT(save_mem)) my_save_mem=save_mem
    NULLIFY(coord_section)
    coord_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
    CALL section_vals_get(coord_section, explicit=explicit, error=error)
    IF (.NOT.explicit) RETURN

    CALL timeset(routineN,handle)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. get cell and topology%atom_info
    !-----------------------------------------------------------------------------
    atom_info => topology%atom_info
    cell => topology%cell_muc
    CALL section_vals_val_get(coord_section,"UNIT",c_val=unit_str,error=error)
    CALL section_vals_val_get(coord_section,"SCALED",l_val=scaled_coordinates,error=error)
    unit_conv = cp_unit_to_cp2k(1.0_dp,TRIM(unit_str),error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Read in the coordinates from &COORD section in the input file
    !-----------------------------------------------------------------------------
    CALL section_vals_val_get(coord_section,"_DEFAULT_KEYWORD_",&
         n_rep_val=natom,error=error)
    topology%natoms = natom
    IF (my_overwrite) THEN
       CPPostcondition(SIZE(atom_info%r,2)==natom, cp_failure_level, routineP, error, failure)
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
            "Overwriting coordinates. Active coordinates read from &COORD section."//&
            " Active coordinates READ from &COORD section "//&
CPSourceFileRef,&
          only_ionode=.TRUE.)
       CALL section_vals_list_get(coord_section,"_DEFAULT_KEYWORD_",list=list,error=error)
       DO iatom=1,natom
          is_ok=cp_sll_val_next(list,val,error=error)
          CALL val_get(val,c_val=line_att,error=error)
          ! Read name and atomic coordinates
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
             CALL cp_assert(LEN_TRIM(line_att(start_c:end_c-1))/=0,cp_failure_level,&
                  cp_assertion_failed,routineP,&
                  "incorrectly formatted line in coord section'"//line_att//"'",&
                  error,failure)
             IF (wrd==1) THEN
                atom_info%id_atmname(iatom)=str2id(s2s(line_att(start_c:end_c-1)))
             ELSE
                READ (line_att(start_c:end_c-1),*) atom_info%r(wrd-1,iatom)
             END IF
             start_c=end_c
          END DO
       END DO
    ELSE
       ! Element is assigned on the basis of the atm_name
       topology%aa_element = .TRUE.

       CALL reallocate(atom_info%id_molname,1,natom)
       CALL reallocate(atom_info%id_resname,1,natom)
       CALL reallocate(atom_info%resid,1,natom)
       CALL reallocate(atom_info%id_atmname,1,natom)
       CALL reallocate(atom_info%id_element,1,natom)
       CALL reallocate(atom_info%r,1,3,1,natom)
       CALL reallocate(atom_info%atm_mass,1,natom)
       CALL reallocate(atom_info%atm_charge,1,natom)

       CALL section_vals_list_get(coord_section,"_DEFAULT_KEYWORD_",list=list,error=error)
       DO iatom=1,natom
          ! we use only the first default_string_length characters of each line
          is_ok=cp_sll_val_next(list,val,error=error)
          CALL val_get(val,c_val=line_att,error=error)
          default_id = str2id(s2s(""))
          atom_info%id_molname(iatom) = default_id
          atom_info%id_resname(iatom) = default_id
          atom_info%resid(iatom) = 1
          atom_info%id_atmname(iatom) = default_id
          atom_info%id_element(iatom) = default_id
          topology%molname_generated = .TRUE.
          ! Read name and atomic coordinates
          start_c=1
          DO wrd=1,6
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
             CALL cp_assert((LEN_TRIM(line_att(start_c:end_c-1)) /= 0),cp_failure_level,&
                            cp_assertion_failed,routineP,&
                            "Incorrectly formatted input line for atom "//&
                            TRIM(ADJUSTL(cp_to_string(iatom)))//&
                            " found in COORD section. Input line: <"//&
                            TRIM(line_att)//"> ",error,failure)
             SELECT CASE (wrd)
             CASE (1)
                atom_info%id_atmname(iatom) = str2id(s2s(line_att(start_c:end_c-1)))
             CASE (2:4)
                CALL read_float_object(line_att(start_c:end_c-1),&
                                       atom_info%r(wrd-1,iatom),error_message)
                CALL cp_assert((LEN_TRIM(error_message) == 0),cp_failure_level,&
                               cp_assertion_failed,routineP,&
                               "Incorrectly formatted input line for atom "//&
                               TRIM(ADJUSTL(cp_to_string(iatom)))//&
                               " found in COORD section. "//TRIM(error_message)//&
                               " Input line: <"//TRIM(line_att)//"> ",error,failure)
             CASE (5)
                READ (line_att(start_c:end_c-1),*) strtmp
                atom_info%id_molname(iatom) = str2id(strtmp)
                atom_info%id_resname(iatom) = atom_info%id_molname(iatom)
                topology%molname_generated = .FALSE.
             CASE (6)
                READ (line_att(start_c:end_c-1),*) strtmp
                atom_info%id_resname(iatom) = str2id(strtmp)
             END SELECT
             start_c = end_c
             IF (start_c > LEN_TRIM(line_att)) EXIT
          END DO
          IF (topology%molname_generated) THEN
             ! Use defaults, if no molname was specified
             WRITE (my_default_index,'(I0)') iatom
             atom_info%id_molname(iatom) = str2id(s2s(TRIM(id2str(atom_info%id_atmname(iatom)))//TRIM(my_default_index)))
             atom_info%id_resname(iatom) = atom_info%id_molname(iatom)
          END IF
          atom_info%id_element(iatom) = atom_info%id_atmname(iatom)
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
          r0 = atom_info%r(:,iatom)
          CALL scaled_to_real(atom_info%r(:,iatom),r0,cell)
       ELSE
          atom_info%r(:,iatom) = atom_info%r(:,iatom)*unit_conv
       END IF
    END DO
    IF(my_save_mem) CALL section_vals_remove_values(coord_section, error)

    CALL timestop(handle)
  END SUBROUTINE read_atoms_input

! *****************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param shell_particle_set ...
!> \param cell ...
!> \param subsys_section ...
!> \param core_particle_set ...
!> \param save_mem ...
!> \param error ...
!> \author MI
! *****************************************************************************
  SUBROUTINE read_shell_coord_input ( particle_set, shell_particle_set, cell,&
             subsys_section, core_particle_set, save_mem, error )

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, &
                                                shell_particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: core_particle_set
    LOGICAL, INTENT(IN), OPTIONAL            :: save_mem
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_shell_coord_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2*default_string_length)   :: line_att
    CHARACTER(len=default_string_length)     :: name_kind, unit_str
    CHARACTER(len=default_string_length), &
      ALLOCATABLE, DIMENSION(:)              :: at_name, at_name_c
    INTEGER                                  :: end_c, handle, ishell, istat, &
                                                j, nshell, output_unit, &
                                                sh_index, start_c, wrd
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: at_index, at_index_c
    LOGICAL :: core_scaled_coordinates, explicit, failure, is_ok, is_shell, &
      my_save_mem, shell_scaled_coordinates
    REAL(KIND=dp)                            :: dab, mass_com, rab(3), &
                                                unit_conv_core, &
                                                unit_conv_shell
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r, rc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_val_type), POINTER           :: list
    TYPE(section_vals_type), POINTER         :: core_coord_section, &
                                                shell_coord_section
    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(val_type), POINTER                  :: val

    failure=.FALSE.
    my_save_mem = .FALSE.
    NULLIFY(logger, atomic_kind, list, shell_coord_section, shell, val)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    IF(PRESENT(save_mem)) my_save_mem=save_mem
    NULLIFY(shell_coord_section, core_coord_section)
    shell_coord_section => section_vals_get_subs_vals(subsys_section,"SHELL_COORD",error=error)
    CALL section_vals_get(shell_coord_section, explicit=explicit, error=error)
    IF (.NOT.explicit) RETURN

    CALL timeset(routineN,handle)
    CPPrecondition(ASSOCIATED(particle_set),cp_failure_level,routineP,error,failure)
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Read in the coordinates from &SHELL_COORD section in the input file
    !-----------------------------------------------------------------------------
    CALL section_vals_val_get(shell_coord_section,"UNIT",c_val=unit_str,error=error)
    CALL section_vals_val_get(shell_coord_section,"SCALED",l_val=shell_scaled_coordinates,error=error)
    unit_conv_shell = cp_unit_to_cp2k(1.0_dp,TRIM(unit_str),error=error)
    CALL section_vals_val_get(shell_coord_section,"_DEFAULT_KEYWORD_",&
         n_rep_val=nshell,error=error)

    IF (ASSOCIATED(shell_particle_set)) THEN
      CPPrecondition((SIZE(shell_particle_set,1) == nshell),cp_failure_level,routineP,error,failure)
      ALLOCATE (r(3,nshell),at_name(nshell),at_index(nshell), STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
           "Overwriting shell coordinates. "//&
           "Active coordinates READ from &SHELL_COORD section. "//&
CPSourceFileRef,&
           only_ionode=.TRUE.)
      CALL section_vals_list_get(shell_coord_section,"_DEFAULT_KEYWORD_",list=list,error=error)
      DO ishell=1,nshell
         ! we use only the first default_string_length characters of each line
         is_ok=cp_sll_val_next(list,val,error=error)
         CALL val_get(val,c_val=line_att,error=error)
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

      IF (PRESENT(core_particle_set)) THEN
        CPPrecondition(ASSOCIATED(core_particle_set),cp_failure_level,routineP,error,failure)
        core_coord_section => section_vals_get_subs_vals(subsys_section,"CORE_COORD",error=error)
        CALL section_vals_get(core_coord_section, explicit=explicit, error=error)
        IF (explicit) THEN
          CALL section_vals_val_get(core_coord_section,"UNIT",c_val=unit_str,error=error)
          CALL section_vals_val_get(core_coord_section,"SCALED",l_val=core_scaled_coordinates,error=error)
          unit_conv_core = cp_unit_to_cp2k(1.0_dp,TRIM(unit_str),error=error)
          CALL section_vals_val_get(core_coord_section,"_DEFAULT_KEYWORD_",&
               n_rep_val=nshell,error=error)

          CPPostcondition((SIZE(core_particle_set,1) == nshell),cp_failure_level,routineP,error,failure)
          ALLOCATE (rc(3,nshell),at_name_c(nshell),at_index_c(nshell), STAT=istat)
          CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Overwriting cores coordinates. "//&
               " Active coordinates READ from &CORE_COORD section. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          CALL section_vals_list_get(core_coord_section,"_DEFAULT_KEYWORD_",list=list,error=error)
          DO ishell=1,nshell
             ! we use only the first default_string_length characters of each line
             is_ok=cp_sll_val_next(list,val,error=error)
             CALL val_get(val,c_val=line_att,error=error)
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
                   at_name_c(ishell)=line_att(start_c:end_c-1)
                   CALL uppercase(at_name_c(ishell))
                ELSE IF (wrd==5) THEN
                   READ (line_att(start_c:end_c-1),*) at_index_c(ishell)
                ELSE
                   READ (line_att(start_c:end_c-1),*) rc(wrd-1,ishell)
                END IF
                start_c=end_c
             END DO
          END DO
          IF(my_save_mem) CALL section_vals_remove_values(core_coord_section, error)
        END IF  ! explicit
      END IF  ! core_particle_set

      !-----------------------------------------------------------------------------
      ! 3. Check corrispondence and convert coordinates into internal cp2k coordinates
      !-----------------------------------------------------------------------------
      DO ishell=1,nshell
         atomic_kind => particle_set(at_index(ishell))%atomic_kind
         CALL get_atomic_kind(atomic_kind=atomic_kind,&
              name=name_kind,shell_active=is_shell,mass=mass_com,shell=shell)
         CALL uppercase(name_kind)
         IF((TRIM(at_name(ishell))==TRIM(name_kind)).AND. is_shell) THEN
           sh_index = particle_set(at_index(ishell))%shell_index
           IF (shell_scaled_coordinates) THEN
             CALL scaled_to_real(r(:,ishell),shell_particle_set(sh_index)%r(:),cell)
           ELSE
             shell_particle_set(sh_index)%r(:) = r(:,ishell)*unit_conv_shell
           END IF
           shell_particle_set(sh_index)%atom_index=at_index(ishell)

           IF(PRESENT(core_particle_set) .AND. .NOT. explicit) THEN
             core_particle_set(sh_index)%r(1) = (mass_com*particle_set(at_index(ishell))%r(1)-&
                   shell%mass_shell*shell_particle_set(sh_index)%r(1))/shell%mass_core
             core_particle_set(sh_index)%r(2) = (mass_com*particle_set(at_index(ishell))%r(2)-&
                   shell%mass_shell*shell_particle_set(sh_index)%r(2))/shell%mass_core
             core_particle_set(sh_index)%r(3) = (mass_com*particle_set(at_index(ishell))%r(3)-&
                   shell%mass_shell*shell_particle_set(sh_index)%r(3))/shell%mass_core
             core_particle_set(sh_index)%atom_index=at_index(ishell)
             rab = pbc(shell_particle_set(sh_index)%r,core_particle_set(sh_index)%r,cell)
           ELSE IF(explicit) THEN
             IF (core_scaled_coordinates) THEN
               CALL scaled_to_real(rc(:,ishell),core_particle_set(sh_index)%r(:),cell)
             ELSE
               core_particle_set(sh_index)%r(:) = rc(:,ishell)*unit_conv_core
             END IF
             core_particle_set(sh_index)%atom_index=at_index_c(ishell)
             rab = pbc(shell_particle_set(sh_index)%r,core_particle_set(sh_index)%r,cell)
             CPPostcondition(TRIM(at_name(ishell))==TRIM(at_name_c(ishell)),cp_failure_level,routineP,error,failure)
             CPPostcondition(at_index(ishell)==at_index_c(ishell),cp_failure_level,routineP,error,failure)
           ELSE
             rab = pbc(shell_particle_set(sh_index)%r,particle_set(at_index(ishell))%r,cell)
           END IF

           dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))
           IF(shell%max_dist>0.0_dp .AND. shell%max_dist<dab ) THEN
             IF(output_unit>0) THEN
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
      DEALLOCATE(rc,at_index_c,at_name_c,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF

    IF (my_save_mem) CALL section_vals_remove_values(shell_coord_section,error)

    CALL timestop(handle)

  END SUBROUTINE read_shell_coord_input

END MODULE atoms_input
