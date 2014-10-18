!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE topology_cp2k

  USE cell_types,                      ONLY: cell_type,&
                                             scaled_to_real
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_object,&
                                             parser_test_next_token
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_get_rval,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'topology_cp2k'

  PRIVATE
  PUBLIC :: read_coordinate_cp2k

CONTAINS

! *****************************************************************************
!> \brief   Read the CP2K &COORD section from an external file, i.e. read
!>          atomic coordinates and molecule/residue information in CP2K format.
!> \param topology ...
!> \param para_env ...
!> \param subsys_section ...
!> \param error ...
!> \date    17.01.2011 (Creation, MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_coordinate_cp2k(topology,para_env,subsys_section,error)

    TYPE(topology_parameters_type)           :: topology
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_coordinate_cp2k', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: string
    INTEGER                                  :: handle, i, iw, natom, &
                                                newsize, number_of_atoms
    LOGICAL                                  :: eof, explicit, failure, &
                                                scaled_coordinates
    REAL(KIND=dp)                            :: pfactor, unit_conv
    REAL(KIND=dp), DIMENSION(3)              :: r
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(section_vals_type), POINTER         :: coord_section

    CALL timeset(routineN,handle)

    NULLIFY (coord_section)
    NULLIFY (logger)
    NULLIFY (parser)
    failure = .FALSE.
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/XYZ_INFO",&
                              extension=".subsysLog",error=error)

    ! Check if there is a &COORD section
    coord_section => section_vals_get_subs_vals(subsys_section,"COORD",error=error)
    CALL section_vals_get(coord_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(coord_section,"UNIT",c_val=string,error=error)
       CALL section_vals_val_get(coord_section,"SCALED",l_val=scaled_coordinates,error=error)
    ELSE
       ! The default is Cartesian coordinates in Angstrom
       scaled_coordinates = .FALSE.
       string = "angstrom"
    END IF
    unit_conv = cp_unit_to_cp2k(1.0_dp,TRIM(string),error=error)

    atom_info => topology%atom_info
    cell => topology%cell_muc

    IF (iw > 0) THEN
       WRITE (UNIT=iw,FMT="(T2,A)")&
        "BEGIN of COORD section data read from file "//TRIM(topology%coord_file_name)
    END IF

    pfactor = section_get_rval(subsys_section,"TOPOLOGY%MEMORY_PROGRESSION_FACTOR",error)
    number_of_atoms = section_get_ival(subsys_section,"TOPOLOGY%NUMBER_OF_ATOMS",error)
    IF (number_of_atoms < 1) THEN
       newsize = 1000
    ELSE
       newsize = number_of_atoms
    END IF

    CALL reallocate(atom_info%id_molname,1,newsize)
    CALL reallocate(atom_info%id_resname,1,newsize)
    CALL reallocate(atom_info%resid,1,newsize)
    CALL reallocate(atom_info%id_atmname,1,newsize)
    CALL reallocate(atom_info%r,1,3,1,newsize)
    CALL reallocate(atom_info%atm_mass,1,newsize)
    CALL reallocate(atom_info%atm_charge,1,newsize)
    CALL reallocate(atom_info%occup,1,newsize)
    CALL reallocate(atom_info%beta,1,newsize)
    CALL reallocate(atom_info%id_element,1,newsize)

    topology%molname_generated = .FALSE.
    ! Element is assigned on the basis of the atm_name
    topology%aa_element = .TRUE.

    CALL parser_create(parser,topology%coord_file_name,para_env=para_env,error=error)

    natom = 0
    DO
       CALL parser_get_object(parser,object=string,newline=.TRUE.,at_end=eof,error=error)
       IF (eof) EXIT
       natom = natom + 1
       IF (natom > SIZE(atom_info%id_atmname)) THEN
          newsize = pfactor*natom
          CALL reallocate(atom_info%id_molname,1,newsize)
          CALL reallocate(atom_info%id_resname,1,newsize)
          CALL reallocate(atom_info%resid,1,newsize)
          CALL reallocate(atom_info%id_atmname,1,newsize)
          CALL reallocate(atom_info%r,1,3,1,newsize)
          CALL reallocate(atom_info%atm_mass,1,newsize)
          CALL reallocate(atom_info%atm_charge,1,newsize)
          CALL reallocate(atom_info%occup,1,newsize)
          CALL reallocate(atom_info%beta,1,newsize)
          CALL reallocate(atom_info%id_element,1,newsize)
       END IF
       atom_info%id_atmname(natom) = str2id(s2s(string))
       DO i=1,3
         CALL parser_get_object(parser,object=r(i),error=error)
       END DO
       IF (scaled_coordinates) THEN
          CALL scaled_to_real(atom_info%r(1:3,natom),r,cell)
       ELSE
          atom_info%r(1:3,natom) = r(1:3)*unit_conv
       END IF
       IF (parser_test_next_token(parser,error=error) /= "EOL") THEN
          CALL parser_get_object(parser,object=string,error=error)
          atom_info%id_molname(natom) = str2id(s2s(string))
          IF (parser_test_next_token(parser,error=error) /= "EOL") THEN
             CALL parser_get_object(parser,object=string,error=error)
             atom_info%id_resname(natom) = str2id(s2s(string))
          ELSE
             atom_info%id_resname(natom) = atom_info%id_molname(natom)
          END IF
       ELSE
          string = ""
          WRITE (UNIT=string,FMT="(I0)") natom
          atom_info%id_molname(natom) = str2id(s2s(TRIM(id2str(atom_info%id_atmname(natom)))//TRIM(string)))
          atom_info%id_resname(natom) = atom_info%id_molname(natom)
          topology%molname_generated = .TRUE.
       END IF
       atom_info%resid(natom) = 1
       atom_info%id_element(natom) = atom_info%id_atmname(natom)
       atom_info%atm_mass(natom)   =  HUGE(0.0_dp)
       atom_info%atm_charge(natom) = -HUGE(0.0_dp)
       IF (iw > 0) THEN
          WRITE (UNIT=iw,FMT="(T2,A,3F20.8,2(2X,A))")&
            TRIM(id2str(atom_info%id_atmname(natom))),atom_info%r(1:3,natom),&
            ADJUSTL(TRIM(id2str(atom_info%id_molname(natom)))),&
            ADJUSTL(TRIM(id2str(atom_info%id_resname(natom))))
       END IF
       IF (natom == number_of_atoms) EXIT
    END DO

    CALL parser_release(parser,error=error)

    topology%natoms = natom
    CALL reallocate(atom_info%id_molname,1,natom)
    CALL reallocate(atom_info%id_resname,1,natom)
    CALL reallocate(atom_info%resid,1,natom)
    CALL reallocate(atom_info%id_atmname,1,natom)
    CALL reallocate(atom_info%r,1,3,1,natom)
    CALL reallocate(atom_info%atm_mass,1,natom)
    CALL reallocate(atom_info%atm_charge,1,natom)
    CALL reallocate(atom_info%occup,1,natom)
    CALL reallocate(atom_info%beta,1,natom)
    CALL reallocate(atom_info%id_element,1,natom)

    CALL section_vals_val_set(subsys_section,"TOPOLOGY%NUMBER_OF_ATOMS",i_val=natom,error=error)

    IF (iw > 0) THEN
       WRITE (UNIT=iw,FMT="(T2,A)")&
        "END of COORD section data read from file "//TRIM(topology%coord_file_name)
    END IF

    CALL cp_print_key_finished_output(iw,logger,subsys_section,&
                                      "PRINT%TOPOLOGY_INFO/XYZ_INFO",&
                                      error=error)

    CALL timestop(handle)

  END SUBROUTINE read_coordinate_cp2k

END MODULE topology_cp2k
