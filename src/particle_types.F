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
!!     Atomic kind added in particle_type (MK, 08.01.2002)
!!     Functionality for particle_type added (MK, 14.01.2002)
!!
!!   SOURCE
!******************************************************************************

MODULE particle_types
  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_type,&
                                             get_atomic_kind,&
                                             set_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp => dp,&
                                             wp_size => dp_size
  USE machine,                         ONLY: m_flush
  USE memory_utilities,                ONLY: reallocate
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: massunit
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             search_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE string_utilities,                ONLY: integer_to_string,&
                                             string_to_integer,&
                                             uppercase
  USE physcon,                         ONLY: angstrom
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: max_name_length = 60

  TYPE particle_type
     TYPE(atomic_kind_type), POINTER       :: atomic_kind
     INTEGER                               :: abase,first_cgf,first_cgf_aux,&
                                              first_sgf,first_sgf_aux,grid_id,&
                                              iatom,kind,last_cgf,last_sgf
     REAL(wp), DIMENSION(3)                :: f,r,v
     INTEGER, DIMENSION(:), POINTER        :: coef_list
  END TYPE particle_type

! *** Public data types ***

  PUBLIC :: particle_type, particle_1d_a_type

! *** Public subroutines ***

  PUBLIC :: allocate_particle_set,&
            deallocate_particle_set,&
            get_particle_set,&
            init_particle_set,&
            read_particle_set,&
            write_qs_particle_coordinates,&
            write_particle_distances, &
            write_particle_coordinates

!!***
! *****************************************************************************

!!****s* particle_types/particle_1d_a_type [1.0] *
!!
!!   NAME
!!     particle_1d_a_type
!!
!!   FUNCTION
!!     represent a pointed to a 1d array of particles (to have arrays of
!!     1d arrays)
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - particle_set: the 1d array of particles
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!   SOURCE
  !***************************************************************************
  TYPE particle_1d_a_type
     TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
  END TYPE particle_1d_a_type
!!***
  !****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE allocate_particle_set(particle_set,natom)

!   Purpose: Allocate a particle set.

!   History: - Creation (14.01.2002, Matthias Krack)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE allocate_particle_set (MODULE particle_types)"

!   *** Local variables ***

    INTEGER :: iatom,istat,natom

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) CALL deallocate_particle_set(particle_set)

    ALLOCATE (particle_set(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"particle_set(natom)",0)

    DO iatom=1,natom
      NULLIFY (particle_set(iatom)%atomic_kind)
      NULLIFY (particle_set(iatom)%coef_list)
    END DO

  END SUBROUTINE allocate_particle_set

! *****************************************************************************

  SUBROUTINE deallocate_particle_set(particle_set)

!   Purpose: Deallocate a particle set.

!   History: - Creation (14.01.2002, Matthias Krack)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE deallocate_particle_set (MODULE particle_types)"

!   *** Local variables ***

    INTEGER :: iatom,istat,natom

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) THEN
      natom = SIZE(particle_set)
      DO iatom=1,natom
        IF (ASSOCIATED(particle_set(iatom)%coef_list)) THEN
          DEALLOCATE (particle_set(iatom)%coef_list,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,&
                                           "particle_set(iatom)%coef_list")
        END IF
      END DO
      DEALLOCATE (particle_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"particle_set")
    ELSE
      CALL stop_program(routine,&
                        "The pointer particle_set is not associated and "//&
                        "cannot be deallocated")
    END IF

  END SUBROUTINE deallocate_particle_set

! *****************************************************************************

  SUBROUTINE get_particle_set(particle_set,r)

!   Purpose: Get the components of a particle set.

!   History: - Creation (14.01.2002, Matthias Krack)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
    REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: r

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE get_particle_set (MODULE particle_types)"

!   *** Local variables ***

    INTEGER :: iatom

!   ---------------------------------------------------------------------------

    IF (PRESENT(r)) THEN
      DO iatom=1,SIZE(particle_set)
        r(:,iatom) = particle_set(iatom)%r(:)
      END DO
    END IF

  END SUBROUTINE get_particle_set

! *****************************************************************************

  SUBROUTINE init_particle_set(particle_set)

!   Purpose: Initialise a particle set.

!   History: - Creation (14.01.2002, Matthias Krack)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_particle_set (MODULE particle_types)"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER :: orb_basis_set

    REAL(wp) :: mass
    INTEGER  :: iatom,icgf,isgf,istat,kind_number,natom,ncgf,nsgf

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(particle_set)) THEN
      natom = SIZE(particle_set)
      icgf = 0
      isgf = 0
      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             orb_basis_set=orb_basis_set,&
                             kind_number=kind_number,&
                             mass=mass)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,ncgf=ncgf,nsgf=nsgf)
        particle_set(iatom)%first_cgf = icgf + 1
        icgf = icgf + ncgf
        particle_set(iatom)%last_cgf = icgf
        particle_set(iatom)%first_sgf = isgf + 1
        isgf = isgf + nsgf
        particle_set(iatom)%last_sgf = isgf
      END DO
    ELSE
      CALL stop_program(routine,"The pointer particle_set is not associated")
    END IF

  END SUBROUTINE init_particle_set

! *****************************************************************************

  SUBROUTINE read_particle_set(particle_set,atomic_kind_set,cell,globenv)

!   Purpose: Read the input section with the atomic coordinates.
!            the format is xmol like : KIND X Y Z

!   History: - Creation (19.05.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_particle_set (MODULE particle_types)"

    INTEGER, PARAMETER :: nblock_atom = 1000,&
                          nblock_kind = 10

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER :: atomic_kind
    CHARACTER(LEN=max_name_length)  :: current_label,label,string
    CHARACTER(LEN=200)              :: message
    CHARACTER(LEN=8)                :: unit_name
    CHARACTER(LEN=3)                :: test_result
    REAL(wp)                        :: unit_of_length
    INTEGER                         :: iatom,ikind,istat,natom,nkind,z
    LOGICAL                         :: found,new_kind,scaled_coordinates,&
                                       section_found

    INTEGER, DIMENSION(:), POINTER :: kind_of,natom_of_kind

    REAL(wp), DIMENSION(:,:), POINTER :: r
    INTEGER, DIMENSION(:,:), POINTER  :: atom_list,kind_label

!   ---------------------------------------------------------------------------

    NULLIFY (atom_list,kind_label,kind_of,natom_of_kind,r)

    nkind = 0
    natom = 0

    CALL reallocate(atom_list,1,nblock_kind,1,nblock_atom)
    CALL reallocate(kind_label,1,max_name_length,1,nblock_atom)
    CALL reallocate(kind_of,1,nblock_atom)
    CALL reallocate(natom_of_kind,1,nblock_kind)
    CALL reallocate(r,1,3,1,nblock_atom)

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="COORD",&
                      section_found=section_found)

    IF (.NOT.section_found) CALL stop_program(routine,&
                                              "No COORD section found.")

    test_result = test_object(newline=.TRUE.)

    DO WHILE (test_result /= "EOS")

      SELECT CASE (test_result)
      CASE ("INT")
        CALL read_object(z,newline=.TRUE.)
        current_label = ptable(z)%symbol
      CASE ("STR")
        CALL read_object(current_label,newline=.TRUE.)
      CASE DEFAULT
        CALL stop_parser(routine,&
                         "Floating point number found for the atomic "//&
                         "label or atomic number")
      END SELECT

      natom = natom + 1

      IF (natom > SIZE(atom_list,2)) THEN
        CALL reallocate(atom_list,1,nkind,1,natom+nblock_atom)
        CALL reallocate(kind_of,1,natom+nblock_atom)
        CALL reallocate(r,1,3,1,natom+nblock_atom)
      END IF

!     *** Read the atomic coordinates ***

      CALL read_object(r(1,natom))
      CALL read_object(r(2,natom))
      CALL read_object(r(3,natom))

!     *** Check, if we found a new atomic kind ***

      label = current_label
      CALL uppercase(label)

      new_kind = .TRUE.

      DO ikind=1,nkind
        CALL integer_to_string(kind_label(:,ikind),string)
        CALL uppercase(string)
        IF (string == label) THEN
          new_kind = .FALSE.
          EXIT
        END IF
      END DO

      IF (new_kind) THEN
        nkind = nkind + 1
        IF (nkind > SIZE(atom_list,1)) THEN
          CALL reallocate(atom_list,1,nkind+nblock_kind,1,natom)
          CALL reallocate(kind_label,1,max_name_length,1,nkind+nblock_kind)
          CALL reallocate(natom_of_kind,1,nkind+nblock_kind)
        END IF
        ikind = nkind
        CALL string_to_integer(current_label,kind_label(:,ikind))
      END IF

      kind_of(natom) = ikind
      natom_of_kind(ikind) = natom_of_kind(ikind) + 1
      atom_list(ikind,natom_of_kind(ikind)) = natom

!     *** Read the next atomic coordinate set ***

      test_result = test_object(newline=.TRUE.)

    END DO

    CALL finish_parser()

    CALL reallocate(atom_list,1,nkind,1,natom)
    CALL reallocate(kind_label,1,max_name_length,1,nkind)
    CALL reallocate(kind_of,1,natom)
    CALL reallocate(natom_of_kind,1,nkind)
    CALL reallocate(r,1,3,1,natom)

!   *** Allocate the data structure for the atomic kind information ***

    CALL allocate_atomic_kind_set(atomic_kind_set,nkind)

!   *** Initialize the data structure for atomic kind information ***

    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL integer_to_string(kind_label(:,ikind),string)
      CALL set_atomic_kind(atomic_kind=atomic_kind,&
                           kind_number=ikind,&
                           name=string,&
                           atom_list=atom_list(ikind,1:natom_of_kind(ikind)))
    END DO

!   *** Allocate the data structure for the atomic information ***

    IF (natom > 0) THEN
      CALL allocate_particle_set(particle_set,natom)
    ELSE
      CALL stop_program(routine,"No atomic coordinates found")
    END IF

!   *** Initialize the data structure for the atomic information ***

    CALL get_cell(cell=cell,&
                  scaled_coordinates=scaled_coordinates,&
                  unit_of_length=unit_of_length)

    DO iatom=1,natom
      ikind = kind_of(iatom)
      particle_set(iatom)%atomic_kind => atomic_kind_set(ikind)
      particle_set(iatom)%iatom = iatom
      particle_set(iatom)%first_cgf = 0
      particle_set(iatom)%first_sgf = 0
      particle_set(iatom)%f(:) = 0.0_wp
      IF (scaled_coordinates) THEN
        particle_set(iatom)%r(:) =&
          scaled_to_real(r(:,iatom),cell)*unit_of_length
      ELSE
        particle_set(iatom)%r(:) = r(:,iatom)*unit_of_length
      END IF
      particle_set(iatom)%v(:) = 0.0_wp
    END DO

!   *** Release work storage ***

    DEALLOCATE (atom_list,kind_label,kind_of,natom_of_kind,r,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine,&
                       "atom_list,kind_label,kind_of,natom_of_kind,r")
    END IF

  END SUBROUTINE read_particle_set

! *****************************************************************************
!  should be able to write a few formats e.g. xmol, and some binary format (dcd)
!  some format can be used for x,v,f
!  should work for fist, qs, kg
!  still needs some rethinking of the units
!  FORMAT   CONTENT                                    UNITS x,v,f
!  XMOL     POS, VEL, FORCE, POS_VEL, POS_VEL_FORCE          Angstrom, a.u., a.u. !hmmmm
! *****************************************************************************
  SUBROUTINE write_particle_coordinates(particle_set,iunit,format,content,title)
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set
    INTEGER                                    :: iunit
    LOGICAL                                    :: ionode
    CHARACTER(LEN=*)                                :: format
    CHARACTER(LEN=*)                                :: content
    CHARACTER(LEN=*)                                :: title

  INTEGER :: natom,iatom
  CHARACTER(LEN=2) :: element_symbol
 
  natom = SIZE(particle_set)

  SELECT CASE (TRIM(format))
  CASE ("XMOL")
       write(iunit,'(I8)') natom
       write(iunit,'(A)')  TRIM(title)
       DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               element_symbol=element_symbol)
          SELECT CASE (TRIM(content))
          CASE ("POS")
               WRITE (iunit,"(T2,A2,X,3F20.10)") element_symbol,particle_set(iatom)%r(1:3)*angstrom
          CASE ("VEL")
               WRITE (iunit,"(T2,A2,X,3F20.10)") element_symbol,particle_set(iatom)%v(1:3)
          CASE ("FORCE")
               WRITE (iunit,"(T2,A2,X,3F20.10)") element_symbol,particle_set(iatom)%f(1:3)
          CASE ("POS_VEL")
               WRITE (iunit,"(T2,A2,X,6F20.10)") element_symbol,particle_set(iatom)%r(1:3)*angstrom, &
                                                                particle_set(iatom)%v(1:3)
          CASE ("POS_VEL_FORCE")
               WRITE (iunit,"(T2,A2,X,9F20.10)") element_symbol,particle_set(iatom)%r(1:3)*angstrom, &
                                                                particle_set(iatom)%v(1:3), &
                                                                particle_set(iatom)%f(1:3)
          END SELECT
       ENDDO
  CASE DEFAULT
    ! others todo
  END SELECT

  END SUBROUTINE write_particle_coordinates

! *****************************************************************************
  SUBROUTINE write_qs_particle_coordinates(particle_set,cell,globenv)

!   Purpose: Write the atomic coordinates to the output unit.

!   History: - Creation (05.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN)  :: globenv
    TYPE(cell_type), POINTER                   :: cell
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set

!   *** Local variables ***

    CHARACTER(LEN=8) :: unit_of_length_name
    CHARACTER(LEN=2) :: element_symbol
    REAL(wp)         :: mass,unit_of_length
    INTEGER          :: iatom,ikind,natom,output_unit,z,zeff

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

      CALL get_cell(cell=cell,&
                    unit_of_length=unit_of_length,&
                    unit_of_length_name=unit_of_length_name)

      output_unit = globenv%scr

!     *** Write headlines ***

      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "ATOMIC COORDINATES IN "//TRIM(unit_of_length_name)
      WRITE (UNIT=output_unit,FMT="(/,T3,A,7X,2(A1,11X),A1,7X,A6,4X,A4,/)")&
        "Atom  Kind  Element","X","Y","Z","Z(eff)","Mass"

      natom = SIZE(particle_set)
      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             kind_number=ikind,&
                             element_symbol=element_symbol,&
                             mass=mass,&
                             z=z,&
                             zeff=zeff)
        WRITE (output_unit,"(T2,I5,1X,I4,3X,A2,2X,I3,3F12.6,4X,I3,3X,F8.4)")&
          iatom,ikind,element_symbol,z,&
          particle_set(iatom)%r(1:3)/unit_of_length,zeff,mass/massunit
      END DO

  END SUBROUTINE write_qs_particle_coordinates

! *****************************************************************************

  SUBROUTINE write_particle_distances(particle_set,cell,globenv)

!   Purpose: Write the matrix of the particle distances to the output unit.

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(cell_type), POINTER                   :: cell
    TYPE(global_environment_type), INTENT(IN)  :: globenv
    TYPE(particle_type), DIMENSION(:), POINTER :: particle_set

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE write_particle_distances (MODULE particle_types)"

!   *** Local variables ***

    CHARACTER(LEN=2) :: element_symbol
    CHARACTER(LEN=8) :: unit_of_length_name
    REAL(wp)         :: dab,dab_image,unit_of_length
    INTEGER          :: from,iatom,icell,icol,ikind,istat,jatom,jcell,kcell,&
                        natom,output_unit,to

    REAL(wp), DIMENSION(3) :: r,ra,rab,rab_image,rb,rb_image
    INTEGER, DIMENSION(3)  :: periodic

    REAL(wp), DIMENSION(:,:), ALLOCATABLE :: distance_matrix

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
      IF (istat /= 0) CALL stop_memory(routine,"distance_matrix(natom,natom)",&
                                       wp_size*natom*natom)

      distance_matrix(:,:) = 0.0_wp

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
                r = REAL((/icell,jcell,kcell/),wp)
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
        WRITE (UNIT=output_unit,FMT="(/,T2,11X,5(4X,I5,4X))")&
          (icol,icol=from,to)
        DO iatom=1,natom
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               kind_number=ikind,&
                               element_symbol=element_symbol)
          WRITE (UNIT=output_unit,FMT="(T2,I5,2X,A2,2X,5(1X,F12.6))")&
            iatom,element_symbol,&
            (distance_matrix(iatom,icol),icol=from,to)
        END DO
      END DO

    END IF

  END SUBROUTINE write_particle_distances

! *****************************************************************************

END MODULE particle_types
