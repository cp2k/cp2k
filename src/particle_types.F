!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/particle_types [1.0] *
!!
!!   NAME
!!     particle_types
!!
!!   FUNCTION
!!     Define the data structure for the particle information.
!!
!!   AUTHOR
!!     CJM, MK
!!
!!   MODIFICATION HISTORY
!!     - Atomic kind added in particle_type (MK,08.01.2002)
!!     - Functionality for particle_type added (MK,14.01.2002)
!!     - Allow for general coordinate input (MK,13.09.2003)
!!     - Molecule concept introduced (MK,26.09.2003)
!!
!!   SOURCE
!******************************************************************************

MODULE particle_types

  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             get_atomic_kind,&
                                             read_atomic_kind_set,&
                                             set_atomic_kind
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl,&
                                             default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE memory_utilities,                ONLY: reallocate
  USE molecule_kind_types,             ONLY: allocate_molecule_kind_set,&
                                             atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             read_molecule_kind_set,&
                                             set_molecule_kind
  USE molecule_types_new,              ONLY: allocate_molecule_set,&
                                             molecule_type,&
                                             set_molecule_set
  USE physcon,                         ONLY: angstrom,&
                                             massunit
  USE qs_parser,                       ONLY: finish_parser,&
                                             integer_object,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             scaled_to_real
  USE string_utilities,                ONLY: ascii_to_string,&
                                             compress,&
                                             string_to_ascii,&
                                             uppercase
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (in this module) ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "particle_types"

! *** Data types ***

  TYPE exclusion_list_type
    INTEGER, DIMENSION(:), POINTER :: list
  END TYPE exclusion_list_type

  TYPE local_particle_type
    TYPE(exclusion_list_type), DIMENSION(:), POINTER :: exclusion
  END TYPE local_particle_type

  TYPE particle_type
     TYPE(atomic_kind_type), POINTER       :: atomic_kind ! atomic kind information
     TYPE(local_particle_type), POINTER    :: lpi ! local particle information
     REAL(KIND = dp), DIMENSION(3)                :: f,& ! force
                                              r,& ! position
                                              v   ! velocity
!MK this components should become obsolete and they are only present
!MK for compatability reasons
     INTEGER                               :: abase,first_cgf,first_cgf_aux,&
                                              first_sgf,first_sgf_aux,grid_id,&
                                              iatom,kind,last_cgf,last_sgf
     INTEGER, DIMENSION(:), POINTER        :: coef_list
  END TYPE particle_type

! *** Public data types ***

  PUBLIC :: exclusion_list_type,&
            local_particle_type,&
            particle_type

! *** Public subroutines ***

  PUBLIC :: allocate_particle_set,&
            deallocate_particle_set,&
            get_particle_set,&
            read_particle_set,&
            write_qs_particle_coordinates,&
            write_particle_distances, &
            write_particle_coordinates

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_particle_set(particle_set,nparticle)

!   Purpose: Allocate a particle set.

!   History: - Creation (14.01.2002,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: nparticle

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_particle_set"

    INTEGER                                  :: iparticle, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) CALL deallocate_particle_set(particle_set)

    ALLOCATE (particle_set(nparticle),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "particle_set",nparticle*int_size)

    DO iparticle=1,nparticle
      NULLIFY (particle_set(iparticle)%atomic_kind)
      NULLIFY (particle_set(iparticle)%coef_list)
      NULLIFY (particle_set(iparticle)%lpi)
      particle_set(iparticle)%abase = 0
      particle_set(iparticle)%grid_id = 0
      particle_set(iparticle)%iatom = iparticle
      particle_set(iparticle)%f(:) = 0.0_dbl
      particle_set(iparticle)%r(:) = 0.0_dbl
      particle_set(iparticle)%v(:) = 0.0_dbl
    END DO

  END SUBROUTINE allocate_particle_set

! *****************************************************************************

  SUBROUTINE deallocate_particle_set(particle_set)

!   Purpose: Deallocate a particle set.

!   History: - Creation (14.01.2002,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_particle_set"

    INTEGER                                  :: iparticle, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) THEN

      DO iparticle=1,SIZE(particle_set)

        IF (ASSOCIATED(particle_set(iparticle)%coef_list)) THEN
          DEALLOCATE (particle_set(iparticle)%coef_list,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                           "particle_set(iparticle)%coef_list")
        END IF

        IF (ASSOCIATED(particle_set(iparticle)%lpi)) THEN
          IF (ASSOCIATED(particle_set(iparticle)%lpi%exclusion)) THEN
            CALL stop_program(routine_name,module_name,__LINE__,&
                              "exclusion list dealloc NYI")
!MK         CALL deallocate_exclusion_list(particle_set(iparticle)%lpi%exclusion)
          END IF
          DEALLOCATE (particle_set(iparticle)%lpi,STAT=istat)
          IF (istat /= 0) THEN
            CALL stop_memory(routine_name,module_name,__LINE__,&
                             "particle_set(iparticle)%lpi")
          END IF
        END IF

      END DO

      DEALLOCATE (particle_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "particle_set")

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer particle_set is not associated and "//&
                        "cannot be deallocated")

    END IF

  END SUBROUTINE deallocate_particle_set

! *****************************************************************************

  SUBROUTINE get_particle_set(particle_set,first_cgf,first_sgf,last_cgf,&
                              last_sgf,ncgf,nsgf,r)

!   Purpose: Get the components of a particle set.

!   History: - Creation (14.01.2002,MK)
!            - particle type cleaned (13.10.2003,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: first_cgf, first_sgf, &
                                                last_cgf, last_sgf, ncgf, nsgf
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT), OPTIONAL                  :: r

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_particle_set"

    INTEGER                                  :: icgf, iparticle, isgf, nc, &
                                                nparticle, ns

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) THEN

      nparticle = SIZE(particle_set)

      IF (PRESENT(first_cgf)) THEN
        IF (SIZE(first_cgf) < nparticle) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <first_cgf>")
        END IF
      END IF

      IF (PRESENT(first_sgf)) THEN
        IF (SIZE(first_sgf) < nparticle) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <first_sgf>")
        END IF
      END IF

      IF (PRESENT(last_cgf)) THEN
        IF (SIZE(last_cgf) < nparticle) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <last_cgf>")
        END IF
      END IF

      IF (PRESENT(last_sgf)) THEN
        IF (SIZE(last_sgf) < nparticle) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <last_sgf>")
        END IF
      END IF

      IF (PRESENT(ncgf)) THEN
        IF (SIZE(ncgf) < nparticle) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <ncgf>")
        END IF
      END IF

      IF (PRESENT(nsgf)) THEN
        IF (SIZE(nsgf) < nparticle) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "An array of insufficient size was supplied "//&
                            "for argument <nsgf>")
        END IF
      END IF

      IF (PRESENT(first_cgf).OR.&
          PRESENT(last_cgf).OR.&
          PRESENT(first_sgf).OR.&
          PRESENT(last_sgf).OR.&
          PRESENT(ncgf).OR.&
          PRESENT(nsgf)) THEN
        icgf = 0
        isgf = 0
        DO iparticle=1,nparticle
          CALL get_atomic_kind(atomic_kind=particle_set(iparticle)%atomic_kind,&
                               ncgf=nc,&
                               nsgf=ns)
          IF (PRESENT(ncgf)) ncgf(iparticle) = nc
          IF (PRESENT(first_cgf)) first_cgf(iparticle) = icgf + 1
          icgf = icgf + nc
          IF (PRESENT(last_cgf)) last_cgf(iparticle) = icgf
          IF (PRESENT(nsgf)) nsgf(iparticle) = ns
          IF (PRESENT(first_sgf)) first_sgf(iparticle) = isgf + 1
          isgf = isgf + ns
          IF (PRESENT(last_sgf)) last_sgf(iparticle) = isgf
        END DO
        IF (PRESENT(first_cgf)) THEN
          IF (SIZE(first_cgf) > nparticle) first_cgf(nparticle+1) = icgf + 1
        END IF
        IF (PRESENT(first_sgf)) THEN
          IF (SIZE(first_sgf) > nparticle) first_sgf(nparticle+1) = isgf + 1
        END IF
      END IF

      IF (PRESENT(r)) THEN
        DO iparticle=1,nparticle
          r(:,iparticle) = particle_set(iparticle)%r(:)
        END DO
      END IF

    ELSE

      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer particle_set is not associated")

    END IF

  END SUBROUTINE get_particle_set

! *****************************************************************************

  SUBROUTINE read_particle_set(particle_set,atomic_kind_set,molecule_set,&
                               molecule_kind_set,cell,globenv)

!   Purpose: Read the input section with the particle (atomic) coordinates.

!   History: - Creation (19.05.2000,MK)
!            - Variable input format including user defined formats
!              (02.09.2003,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_particle_set"
    CHARACTER(LEN=4), PARAMETER              :: required_descriptors = "AXYZ"
    INTEGER, PARAMETER                       :: nblock1 = 1, nblock2 = 1   

    CHARACTER(LEN=1)                         :: descriptor
    CHARACTER(LEN=200)                       :: message
    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=80)                        :: default_format, gro_format, &
                                                input_format, pdb_format, &
                                                xyz_format
    CHARACTER(LEN=default_string_length)     :: label, molecule_name, &
                                                particle_name, string
    INTEGER :: atom_a, handle, i, iatom, idummy, ikind, imolecule, &
      imolecule_kind, input_format_length, iparticle, iparticle_kind, istat, &
      natom, ndigits, nmolecule, nmolecule_kind, nparticle, nparticle_kind, &
      string_length
    INTEGER, DIMENSION(:), POINTER :: first_atom_of_molecule, list, &
      molecule_kind_of, nmolecule_of_kind, nparticle_of_kind, particle_kind_of
    INTEGER, DIMENSION(:, :), POINTER        :: molecule_kind_name, &
                                                particle_kind_name
    LOGICAL                                  :: new_kind, newline, pdb, &
                                                scaled_coordinates, &
                                                section_found
    REAL(KIND=dp)                            :: fdummy, unit_of_length
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: r
    TYPE(atom_type), DIMENSION(:), POINTER   :: atom_list
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind

! blocking for number of particles
! blocking for number of kinds
!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I",'',handle)

    NULLIFY (atom_list)
    NULLIFY (first_atom_of_molecule)
    NULLIFY (list)
    NULLIFY (molecule_kind_name)
    NULLIFY (molecule_kind_of)
    NULLIFY (nmolecule_of_kind)
    NULLIFY (nparticle_of_kind)
    NULLIFY (particle_kind_name)
    NULLIFY (particle_kind_of)
    NULLIFY (r)

    CALL reallocate(first_atom_of_molecule,1,nblock1)
    CALL reallocate(molecule_kind_name,1,default_string_length,1,nblock2)
    CALL reallocate(molecule_kind_of,1,nblock1)
    CALL reallocate(nparticle_of_kind,1,nblock2)
    CALL reallocate(particle_kind_name,1,default_string_length,1,nblock2)
    CALL reallocate(particle_kind_of,1,nblock1)
    CALL reallocate(r,1,3,1,nblock1)

    pdb = .FALSE.

!   *** Define the standard input formats using the following descriptors:

!   *** A)tomic kind
!   *** F)loating point number (not considered)
!   *** I)nteger number (not considered)
!   *** M)olecule kind
!   *** S)tring (not considered)
!   *** X)coordinate
!   *** Y)coordinate
!   *** Z)coordinate

    xyz_format = "A,X,Y,Z"
    pdb_format = "S6,I5,S1,A4,S1,M4,S1,I4,S4,X8,Y8,Z8"
    gro_format = "I5,M5,A5,I5,X,Y,Z"

    default_format = xyz_format

    string_length = LEN_TRIM(globenv%coord_file_name)

    IF (globenv%coord_file_name(MAX(1,string_length-3):string_length) == ".pdb") pdb = .TRUE.

    IF (pdb) THEN

      input_format = pdb_format

      CALL start_parser(file_name=globenv%coord_file_name,&
                        section_char=" ",&
                        globenv=globenv)

    ELSE

      CALL start_parser(file_name=globenv%coord_file_name,&
                        globenv=globenv,&
                        start_section_label="COORD",&
                        section_found=section_found)

      IF (.NOT.section_found) THEN
        CALL stop_parser(routine_name,module_name,__LINE__,"SECTION_NOT_FOUND")
      END IF

      test_result = test_object()

      IF (test_result == "STR") THEN
        CALL read_object(input_format,lower_to_upper=.TRUE.)
        CALL compress(input_format,full=.TRUE.)
        SELECT CASE (TRIM(input_format))
        CASE ("GROMOS87","GRO")
          input_format = gro_format
        CASE ("PDB")
          input_format = pdb_format
          pdb = .TRUE.
        CASE ("XYZ")
          input_format = xyz_format
        END SELECT
      ELSE
        input_format = default_format
      END IF

    END IF

    input_format_length = LEN_TRIM(input_format)
!MKprint*,""""//TRIM(input_format)//"""",input_format_length
!   *** Check for valid input format ***

    DO i=1,LEN_TRIM(required_descriptors)
      descriptor = required_descriptors(i:i)
      IF (INDEX(input_format,descriptor) == 0) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The required descriptor <"//descriptor//&
                          "> is missing in the coordinate format <"//&
                          TRIM(input_format)//">",globenv)
      ELSE IF (INDEX(input_format,descriptor,BACK=.FALSE.) /=&
               INDEX(input_format,descriptor,BACK=.TRUE.)) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "The descriptor <"//descriptor//"> appears more "//&
                          "than once in the coordinate format <"//&
                          TRIM(input_format)//">",globenv)
      END IF
    END DO

    IF (INDEX(input_format,"M",BACK=.FALSE.) /=&
        INDEX(input_format,"M",BACK=.TRUE.)) THEN
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The descriptor <M> appears more than once in the "//&
                        "coordinate format <"//TRIM(input_format)//">",globenv)
    END IF

!   *** Read the particle coordinate and kind information ***

    nmolecule_kind = 0
    nparticle = 0
    nparticle_kind = 0

    line_loop: DO

      test_result = test_object(newline=.TRUE.)

      IF (test_result == "EOS") THEN

        EXIT

      ELSE IF (test_result == "EOF") THEN

        CALL stop_parser(routine_name,module_name,__LINE__,"EOF")

      ELSE

        nparticle = nparticle + 1

        IF (nparticle > SIZE(particle_kind_of)) THEN
          CALL reallocate(particle_kind_of,1,nparticle+nblock1)
          CALL reallocate(r,1,3,1,nparticle+nblock1)
        END IF

!       *** Scan the input format ***

        particle_name = ""
        molecule_name = ""

        newline = .TRUE.

        i = 1

        format_loop: DO

          IF (i > input_format_length) EXIT

          descriptor = input_format(i:i)

          ndigits = INDEX(input_format(i+1:input_format_length)//",",",") - 1

          IF (ndigits > 0) THEN
            IF (integer_object(input_format(i+1:i+ndigits))) THEN
              READ (UNIT=input_format(i+1:i+ndigits),FMT=*) string_length
            ELSE
              CALL stop_parser(routine_name,module_name,__LINE__,&
                               "An integer number was expected for the input "//&
                               "format length, found <"//&
                               input_format(i+1:i+ndigits)//">")
            END IF
          ELSE
            string_length = 0
          END IF
!MKprint*,descriptor," -> string_length=",string_length,"  ndigits=",ndigits
          SELECT CASE (descriptor)
          CASE ("A")
            CALL read_object(object=particle_name,&
                             lower_to_upper=.TRUE.,&
                             newline=newline,&
                             string_length=string_length)
!MKprint*,"A="""//TRIM(particle_name)//"""",string_length
          CASE ("F")
            CALL read_object(object=fdummy,&
                             newline=newline,&
                             string_length=string_length)
!MKprint*,"F=",fdummy,string_length
          CASE ("I")
            CALL read_object(object=idummy,&
                             newline=newline,&
                             string_length=string_length)
!MKprint*,"I=",idummy,string_length
          CASE ("M")
            CALL read_object(object=molecule_name,&
                             lower_to_upper=.TRUE.,&
                             newline=newline,&
                             string_length=string_length)
!MKprint*,"M="""//TRIM(molecule_name)//"""",string_length
          CASE ("S")
            CALL read_object(object=string,&
                             newline=newline,&
                             string_length=string_length)
            IF (pdb.AND.(i == 1)) THEN
              CALL uppercase(string)
              IF ((string /= "ATOM").OR.(string == "HETATM")) THEN
                nparticle = nparticle - 1
                CYCLE line_loop
              END IF
            END IF
!MKprint*,"S="""//TRIM(string)//"""",string_length
          CASE ("X")
            CALL read_object(object=r(1,nparticle),&
                             newline=newline,&
                             string_length=string_length)
!MKprint*,"X=",r(1,nparticle),string_length
          CASE ("Y")
            CALL read_object(object=r(2,nparticle),&
                             newline=newline,&
                             string_length=string_length)
!MKprint*,"Y=",r(2,nparticle),string_length
          CASE ("Z")
            CALL read_object(object=r(3,nparticle),&
                             newline=newline,&
                             string_length=string_length)
!MKprint*,"Z=",r(3,nparticle),string_length
          END SELECT

          i = i + ndigits + 2

          newline = .FALSE.

        END DO format_loop

!       *** Check, if we found a new particle kind ***

        IF (LEN_TRIM(particle_name) > 0) THEN

          label = particle_name
          CALL uppercase(label)

          new_kind = .TRUE.

          DO ikind=1,nparticle_kind
            CALL ascii_to_string(particle_kind_name(:,ikind),string)
            CALL uppercase(string)
            IF (string == label) THEN
              new_kind = .FALSE.
              EXIT
            END IF
          END DO

          IF (new_kind) THEN
            nparticle_kind = nparticle_kind + 1
            IF (nparticle_kind > SIZE(nparticle_of_kind,1)) THEN
              CALL reallocate(particle_kind_name,1,default_string_length,&
                                                 1,nparticle_kind+nblock2)
              CALL reallocate(nparticle_of_kind,1,nparticle_kind+nblock2)
            END IF
            ikind = nparticle_kind
            CALL string_to_ascii(particle_name,particle_kind_name(:,ikind))
          END IF

          particle_kind_of(nparticle) = ikind
          nparticle_of_kind(ikind) = nparticle_of_kind(ikind) + 1

        END IF

!       *** Check, if we found a new molecule kind ***

        IF (LEN_TRIM(molecule_name) > 0) THEN

          label = molecule_name
          CALL uppercase(label)

          new_kind = .TRUE.

          DO ikind=1,nmolecule_kind
            CALL ascii_to_string(molecule_kind_name(:,ikind),string)
            CALL uppercase(string)
            IF (string == label) THEN
              new_kind = .FALSE.
              EXIT
            END IF
          END DO

          IF (new_kind) THEN
            nmolecule_kind = nmolecule_kind + 1
            IF (nmolecule_kind > SIZE(molecule_kind_name,2)) THEN
              CALL reallocate(molecule_kind_name,1,default_string_length,&
                                                  1,nmolecule_kind+nblock2)
            END IF
            ikind = nmolecule_kind
            CALL string_to_ascii(molecule_name,molecule_kind_name(:,ikind))
          END IF

          IF (nparticle > SIZE(molecule_kind_of)) THEN
            CALL reallocate(molecule_kind_of,1,nparticle+nblock1)
          END IF

          molecule_kind_of(nparticle) = ikind

        END IF

      END IF

    END DO line_loop

!   *** Reallocate the work arrays with the final sizes ***

    IF (nparticle > 0) THEN
      IF (nmolecule_kind > 0) THEN
        CALL reallocate(molecule_kind_name,1,default_string_length,1,nmolecule_kind)
      END IF
      CALL reallocate(molecule_kind_of,1,nparticle)
      CALL reallocate(nparticle_of_kind,1,nparticle_kind)
      CALL reallocate(particle_kind_name,1,default_string_length,1,nparticle_kind)
      CALL reallocate(particle_kind_of,1,nparticle)
      CALL reallocate(r,1,3,1,nparticle)
    ELSE
      CALL stop_parser(routine_name,module_name,__LINE__,&
                       "No particle coordinates found")
    END IF

    CALL finish_parser()

!   *** Allocate the data structure for the atomic kind information ***

    CALL allocate_atomic_kind_set(atomic_kind_set,nparticle_kind)

    DO ikind=1,nparticle_kind

      atomic_kind => atomic_kind_set(ikind)

      CALL reallocate(list,1,nparticle_of_kind(ikind))

      i = 0
      DO iparticle=1,nparticle
        IF (particle_kind_of(iparticle) == ikind) THEN
          i = i + 1
          list(i) = iparticle
        END IF
      END DO

      CALL ascii_to_string(particle_kind_name(:,ikind),string)

      CALL set_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=list,&
                           kind_number=ikind,&
                           name=string)

    END DO

!   *** Read the atomic kind information ***

    CALL read_atomic_kind_set(atomic_kind_set,globenv)

!   *** Allocate the data structures for the molecule kind and the molecule ***
!   *** information. Each atomic kind defaults to a molecule kind, if no    ***
!   *** molecule kinds are supplied                                         ***

    IF (nmolecule_kind > 0) THEN

      CALL allocate_molecule_kind_set(molecule_kind_set,nmolecule_kind)

      DO ikind=1,nmolecule_kind

        molecule_kind => molecule_kind_set(ikind)

        CALL ascii_to_string(molecule_kind_name(:,ikind),string)

        CALL set_molecule_kind(molecule_kind=molecule_kind,&
                               kind_number=ikind,&
                               name=string)

      END DO

      CALL read_molecule_kind_set(molecule_kind_set,atomic_kind_set,globenv)

    ELSE

!     *** Each atomic kind becomes a molecule kind, ***
!     *** if no molecule kinds are specified        ***

      nmolecule_kind = nparticle_kind
      molecule_kind_of(:) = particle_kind_of(:)

      CALL allocate_molecule_kind_set(molecule_kind_set,nmolecule_kind)

      ALLOCATE (atom_list(1),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom_list",int_size)

      DO ikind=1,nmolecule_kind

        atomic_kind => atomic_kind_set(ikind)

        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             name=string)

        atom_list(1)%atomic_kind => atomic_kind
        atom_list(1)%name = string

        molecule_kind => molecule_kind_set(ikind)

        CALL set_molecule_kind(molecule_kind=molecule_kind,&
                               atom_list=atom_list,&
                               kind_number=ikind,&
                               name=string)

      END DO

      NULLIFY (atom_list(1)%atomic_kind)
      DEALLOCATE (atom_list,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom_list",int_size)

    END IF

!   *** Allocate the data structure for the particle information ***

    CALL allocate_particle_set(particle_set,nparticle)

!   *** Initialize the data structure for the particle information ***

    CALL get_cell(cell=cell,&
                  scaled_coordinates=scaled_coordinates,&
                  unit_of_length=unit_of_length)

    DO iparticle=1,nparticle
      iparticle_kind = particle_kind_of(iparticle)
      particle_set(iparticle)%atomic_kind => atomic_kind_set(iparticle_kind)
      particle_set(iparticle)%iatom = iparticle
      particle_set(iparticle)%f(:) = 0.0_dp
      IF (scaled_coordinates) THEN
        particle_set(iparticle)%r(:) =&
          scaled_to_real(r(:,iparticle),cell)*unit_of_length
      ELSE
        particle_set(iparticle)%r(:) = r(:,iparticle)*unit_of_length
      END IF
      particle_set(iparticle)%v(:) = 0.0_dp
    END DO

!   *** Count the molecules (or atoms) and prepare the molecule setup ***

    CALL reallocate(nmolecule_of_kind,1,nmolecule_kind)

    nmolecule = 0
    iparticle = 1

    DO

      nmolecule = nmolecule + 1

      IF (nmolecule > SIZE(first_atom_of_molecule)) THEN
        CALL reallocate(first_atom_of_molecule,1,nmolecule+nblock1)
      END IF

      first_atom_of_molecule(nmolecule) = iparticle

      imolecule_kind = molecule_kind_of(iparticle)

      molecule_kind => molecule_kind_set(imolecule_kind)

      CALL get_molecule_kind(molecule_kind=molecule_kind,&
                             atom_list=atom_list,&
                             name=molecule_name,&
                             natom=natom)

      nmolecule_of_kind(imolecule_kind) = nmolecule_of_kind(imolecule_kind) + 1

      IF (natom < 1) THEN
        CALL stop_program(routine_name,module_name,__LINE__,&
                          "Invalid number of atoms for molecule kind found.",&
                          globenv)
      END IF

!     *** Consistency check ***

      DO iatom=1,natom ! loop over all atoms of the current molecule
        atom_a = iparticle + iatom - 1
        IF (atom_a > nparticle) THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
                            "Incomplete coordinate set. Atoms are missing "//&
                            "for the definition of the last molecule.",globenv)
        END IF
        CALL get_atomic_kind(atomic_kind=particle_set(atom_a)%atomic_kind,&
                             name=particle_name)
        CALL get_atomic_kind(atomic_kind=atom_list(iatom)%atomic_kind,&
                             name=string)
        IF (particle_name /= string) THEN
          WRITE (UNIT=message,FMT="(A,3(I6,A))")&
            "The atomic kind <"//TRIM(particle_name)//"> of atom ",atom_a,&
            " in the coordinate set does not coincide with the "//&
            "atomic kind <"//TRIM(string)//"> of atom ",iatom,&
            " of molecule ",nmolecule," which belongs to the molecule "//&
            "kind <"//TRIM(molecule_name)//">."
          CALL compress(message)
          CALL stop_program(routine_name,module_name,__LINE__,message,globenv)
        END IF
      END DO

      iparticle = iparticle + natom

      IF (iparticle > nparticle) EXIT

    END DO

    CALL reallocate(first_atom_of_molecule,1,nmolecule)

!   *** Allocate and initialize the molecule set ***

    CALL allocate_molecule_set(molecule_set,nmolecule)

    CALL set_molecule_set(molecule_set=molecule_set,&
                          first_atom=first_atom_of_molecule)

!   *** Build the molecule index list for each molecule kind ***

    DO imolecule_kind=1,nmolecule_kind
      molecule_kind => molecule_kind_set(imolecule_kind)
      CALL reallocate(list,1,nmolecule_of_kind(imolecule_kind))
      i = 0
      DO imolecule=1,nmolecule
        iparticle = first_atom_of_molecule(imolecule)
        IF (molecule_kind_of(iparticle) == imolecule_kind) THEN
          i = i + 1
          list(i) = imolecule
        END IF
      END DO
      CALL set_molecule_kind(molecule_kind=molecule_kind,&
                             molecule_list=list)
    END DO

!   *** Release work storage ***

    IF (ASSOCIATED(first_atom_of_molecule)) THEN
      DEALLOCATE (first_atom_of_molecule,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "first_atom_of_molecule")
    END IF

    IF (ASSOCIATED(list)) THEN
      DEALLOCATE (list,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "list")
    END IF

    IF (ASSOCIATED(molecule_kind_name)) THEN
      DEALLOCATE (molecule_kind_name,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "molecule_kind_name")
    END IF

    IF (ASSOCIATED(molecule_kind_of)) THEN
      DEALLOCATE (molecule_kind_of,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "molecule_kind_of")
    END IF

    IF (ASSOCIATED(nmolecule_of_kind)) THEN
      DEALLOCATE (nmolecule_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "nmolecule_of_kind")
    END IF

    IF (ASSOCIATED(nparticle_of_kind)) THEN
      DEALLOCATE (nparticle_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "nparticle_of_kind")
    END IF

    IF (ASSOCIATED(particle_kind_name)) THEN
      DEALLOCATE (particle_kind_name,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "particle_kind_name")
    END IF

    IF (ASSOCIATED(particle_kind_of)) THEN
      DEALLOCATE (particle_kind_of,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "particle_kind_of")
    END IF

    IF (ASSOCIATED(r)) THEN
      DEALLOCATE (r,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "r")
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE read_particle_set

! *****************************************************************************
!  should be able to write a few formats e.g. xmol, and some binary format (dcd)
!  some format can be used for x,v,f
!  should work for fist, qs, kg
!  still needs some rethinking of the units
!  FORMAT   CONTENT                                    UNITS x,v,f
!  XMOL     POS, VEL, FORCE, POS_VEL, POS_VEL_FORCE    Angstrom, a.u., a.u. !hmmmm
! *****************************************************************************

  SUBROUTINE write_particle_coordinates(particle_set,cell,iunit,format,&
                                        content,title)

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    INTEGER                                  :: iunit
    CHARACTER(LEN=*)                         :: format, content, title

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: iatom, natom
    REAL(KIND=dp), DIMENSION(3)              :: r

    natom = SIZE(particle_set)

    SELECT CASE (TRIM(format))
    CASE ("XMOL")
      WRITE(iunit,'(I8)') natom
      WRITE(iunit,'(A)')  TRIM(title)
      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             element_symbol=element_symbol)
       SELECT CASE (TRIM(content))
         CASE ("POS")
          r(:) = particle_set(iatom)%r(:)*angstrom
          WRITE (iunit,"(T2,A2,1X,3F20.10)") element_symbol,r(1:3)
        CASE ("VEL")
          WRITE (iunit,"(T2,A2,1X,3F20.10)") element_symbol,&
                                            particle_set(iatom)%v(1:3)
        CASE ("FORCE")
          WRITE (iunit,"(T2,A2,1X,3F20.10)") element_symbol,&
                                            particle_set(iatom)%f(1:3)
        CASE ("POS_VEL")
          r(:) = particle_set(iatom)%r(:)*angstrom
          WRITE (iunit,"(T2,A2,1X,6F20.10)") element_symbol,r(1:3),&
                                            particle_set(iatom)%v(1:3)
        CASE ("POS_VEL_FORCE")
          r(:) = particle_set(iatom)%r(:)*angstrom
          WRITE (iunit,"(T2,A2,1X,9F20.10)") element_symbol,r(1:3),&
                                            particle_set(iatom)%v(1:3),&
                                            particle_set(iatom)%f(1:3)
        END SELECT
      END DO
    CASE ("ATOMIC")
      DO iatom=1,natom
        SELECT CASE (TRIM(content))
        CASE ("POS")
         r(:) = particle_set(iatom)%r(:)
         WRITE (iunit,"(3F20.10)") r(1:3)
        CASE ("VEL")
          WRITE (iunit,"(3F20.10)")particle_set(iatom)%v(1:3)
        CASE ("FORCE")
          WRITE (iunit,"(3F20.10)") particle_set(iatom)%f(1:3)
        CASE ("POS_VEL")
          WRITE (iunit,"(6F20.10)")particle_set(iatom)%r(1:3),&
                                   particle_set(iatom)%v(1:3)
        CASE ("POS_VEL_FORCE")
          r(:) = particle_set(iatom)%r(:)
          WRITE (iunit,"(9F20.10)")r(1:3),&
                                   particle_set(iatom)%v(1:3),&
                                   particle_set(iatom)%f(1:3)
        END SELECT
      END DO
    CASE DEFAULT
    ! others todo
    END SELECT

  END SUBROUTINE write_particle_coordinates

! *****************************************************************************

  SUBROUTINE write_qs_particle_coordinates(particle_set,cell,globenv)

!   Purpose: Write the atomic coordinates to the output unit.

!   History: - Creation (05.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER                                  :: iatom, ikind, natom, &
                                                output_unit, z, zeff
    REAL(KIND=dp)                            :: mass, unit_of_length

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%atomic_coordinates) THEN

      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)

      output_unit = globenv%scr

!     *** Write headlines ***

      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "ATOMIC COORDINATES IN "//TRIM(unit_of_length_name)
      WRITE (UNIT=output_unit,&
           FMT="(/,T3,A,'       ',2(A1,11X),A1,'       ',A6,'    ',A4,/)")&
        "Atom  Kind  Element","X","Y","Z","Z(eff)","Mass"

      natom = SIZE(particle_set)
      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             kind_number=ikind,&
                             element_symbol=element_symbol,&
                             mass=mass,&
                             z=z,&
                             zeff=zeff)
        WRITE (output_unit,&
             "(T2,I5,' ',I4,'   ',A2,'  ',I3,3F12.6,'    ',I3,'   ',F8.4)")&
          iatom,ikind,element_symbol,z,&
          particle_set(iatom)%r(1:3)/unit_of_length,zeff,mass/massunit
      END DO

    END IF

  END SUBROUTINE write_qs_particle_coordinates

! *****************************************************************************

  SUBROUTINE write_particle_distances(particle_set,cell,globenv)

!   Purpose: Write the matrix of the particle distances to the output unit.

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_particle_distances"

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=8)                         :: unit_of_length_name
    INTEGER                                  :: from, iatom, icell, icol, &
                                                ikind, istat, jatom, jcell, &
                                                kcell, natom, output_unit, to
    INTEGER, DIMENSION(3)                    :: periodic
    REAL(KIND=dp)                            :: dab, dab_image, unit_of_length
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: distance_matrix
    REAL(KIND=dp), DIMENSION(3)              :: r, ra, rab, rab_image, rb, &
                                                rb_image

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%interatomic_distances) THEN

      CALL get_cell(cell=cell,&
                    periodic=periodic,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)

      output_unit = globenv%scr

      natom = SIZE(particle_set)

      ALLOCATE (distance_matrix(natom,natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "distance_matrix",dp_size*natom*natom)

      distance_matrix(:,:) = 0.0_dp

      DO iatom=1,natom

        ra(:) = pbc(particle_set(iatom)%r(:),cell)

        DO jatom=iatom+1,natom

          rb(:) = pbc(particle_set(jatom)%r(:),cell)

          rab(:) = rb(:) - ra(:)

          dab = SQRT(rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3))

!         *** Search for the image of atom B which is closest to atome A ***

          DO icell=-periodic(1),periodic(1)
            DO jcell=-periodic(2),periodic(2)
              DO kcell=-periodic(3),periodic(3)
                r = REAL((/icell,jcell,kcell/),dp)
                rb_image(:) = rb(:) + scaled_to_real(r,cell)
                rab_image(:) = rb_image(:) - ra(:)
                dab_image = SQRT(rab_image(1)*rab_image(1) +&
                                 rab_image(2)*rab_image(2) +&
                                 rab_image(3)*rab_image(3))
                IF (dab_image < dab) dab = dab_image
              END DO
            END DO
          END DO

          distance_matrix(iatom,jatom) = dab/unit_of_length
          distance_matrix(jatom,iatom) = distance_matrix(iatom,jatom)

        END DO

      END DO

!     *** Print the distance matrix ***

      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "INTERATOMIC DISTANCES IN "//TRIM(unit_of_length_name)

      DO jatom=1,natom,5
        from = jatom
        to = MIN(from+4,natom)
        WRITE (UNIT=output_unit,FMT="(/,T2,'           ',5('    ',I5,'    '))")&
          (icol,icol=from,to)
        DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               kind_number=ikind,&
                               element_symbol=element_symbol)
          WRITE (UNIT=output_unit,FMT="(T2,I5,'  ',A2,'  ',5(' ',F12.6))")&
            iatom,element_symbol,&
            (distance_matrix(iatom,icol),icol=from,to)
        END DO
      END DO

    END IF

  END SUBROUTINE write_particle_distances

! *****************************************************************************

END MODULE particle_types
