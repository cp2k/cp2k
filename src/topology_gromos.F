!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/topology_gromos [1.0] *
!!
!!   NAME
!!     topology_gromos
!!
!!   FUNCTION
!!     Read and translate a topologuy file generated for gromos
!!     into information needed by cp2k
!!
!!   AUTHOR
!!     JGH 09.12.2001
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE topology_gromos

!  USE global_types,                    ONLY: global_environment_type
!  USE kinds,                           ONLY: dbl
!  USE message_passing,                 ONLY: mp_bcast
!  USE molecule_types,                  ONLY: molecule_type
!  USE parser,                          ONLY: get_next,&
!                                             parser_end,&
!                                             parser_init,&
!                                             read_line,&
!                                             search_label,&
!                                             test_next
!  USE string_utilities,                ONLY: integer_to_string,&
!                                             str_search
!  USE termination,                     ONLY: stop_memory,&
!                                             stop_program
!  USE util,                            ONLY: get_unit
!
!  TYPE solute_atom
!   INTEGER :: atnm                    ! atom number
!   INTEGER :: mres                    ! residue number
!   CHARACTER ( LEN=10 ) :: panm       ! atom name of solute atom
!   INTEGER :: iac                     ! integer (van der Waals) atom type code
!   REAL ( dbl ) :: mass               ! mass of solute atom
!   REAL ( dbl ) :: cg                 ! charge of solute atom
!   INTEGER :: cgc                     ! charge group code (0 or 1)
!   INTEGER :: ine                     ! number of excluded atoms
!   INTEGER, DIMENSION (:), POINTER :: excl
!   INTEGER :: ine14                   ! number of 1-4 interactions
!   INTEGER, DIMENSION (:), POINTER :: int14
!  END TYPE solute_atom
!
!  TYPE solvent_atom
!    INTEGER :: i                  ! solvent atom sequence number
!    INTEGER :: iacs               ! (van der Waals) atom type code
!    CHARACTER (LEN=10) :: anms    ! atom name of solvent atom
!    REAL ( dbl ) :: mass          ! mass of solvent atom
!    REAL ( dbl ) :: cgs           ! charge of solvent atom
!  END TYPE solvent_atom
!
!  TYPE ljtype
!    INTEGER :: iac, jac           ! (van der Waals) atom type code
!    REAL ( dbl ) :: c12           ! r**(-12) term in nonbonded interactions
!    REAL ( dbl ) :: c6            ! r**(-6) term in nonbonded interactions
!    REAL ( dbl ) :: cs12          ! r**(-12) term in 1-4 nonbonded interactions
!    REAL ( dbl ) :: cs6           ! r**(-6) term in 1-4 nonbonded interactions
!  END TYPE ljtype
!
!  TYPE gromos_topology_type
!     CHARACTER (LEN=80), DIMENSION (5) :: title   !title cards
!     REAL ( dbl ) :: frepsi, hbar                 !phys. constants
!     CHARACTER (LEN=10) :: version                !version number
!     INTEGER :: nratt                             !number of atom type names
!     CHARACTER ( LEN=10 ), DIMENSION (:), POINTER :: atomtypename
!     INTEGER :: nraa2                             !number of residues
!     CHARACTER ( LEN=10 ), DIMENSION (:), POINTER :: residue_name
!     INTEGER :: nrp                               !number of solute atoms
!     INTEGER, DIMENSION (:,:), POINTER :: resat   !first and last atom of residue
!     TYPE ( solute_atom ), DIMENSION (:), POINTER :: solute
!     INTEGER :: nbty                              !number of bond types
!     REAL ( dbl ), DIMENSION (:), POINTER :: cb,& !force constant
!                                             b0   !bond length
!     INTEGER :: nbonh                             !nr bonds involving H
!     INTEGER, DIMENSION (:,:), POINTER :: IJBH    !H bonds
!     INTEGER :: nbon                              !nr bonds 
!     INTEGER, DIMENSION (:,:), POINTER :: IJBO    !bonds
!     INTEGER :: ntty                              !number bond angle types
!     REAL ( dbl ), DIMENSION (:), POINTER :: ct,& !force constant
!                                             t0   !bond angle
!     INTEGER :: ntheh                             !bond angles involving H
!     INTEGER, DIMENSION (:,:), POINTER :: IJKTH   !bonds with H
!     INTEGER :: nthe                              !bond angles
!     INTEGER, DIMENSION (:,:), POINTER :: IJKT    !bonds
!     INTEGER :: nqty                              !improper dihedral types
!     REAL ( dbl ), DIMENSION (:), POINTER :: cq,& !force constant
!                                             q0   !improper dihedral angle
!     INTEGER :: nqhih                             !improper DH involving H
!     INTEGER, DIMENSION (:,:), POINTER :: IJKQH   !improper DH with H
!     INTEGER :: nqhi                              !improper DH
!     INTEGER, DIMENSION (:,:), POINTER :: IJKQ    !DH without H
!     INTEGER :: npty                              !dihedrals types
!     REAL ( dbl ), DIMENSION (:), POINTER :: cp,& !force constant
!                                             pd   !dihedral angles
!     INTEGER, DIMENSION (:), POINTER :: np        !multiplicity
!     INTEGER :: nphih                             !DH involving H
!     INTEGER, DIMENSION (:,:), POINTER :: IJKPH   !DH with H
!     INTEGER :: nphi                              !DH
!     INTEGER, DIMENSION (:,:), POINTER :: IJKP    !DH without H
!     INTEGER :: nratt2                            !LJ interaction types
!     TYPE (ljtype), DIMENSION (:), POINTER :: ljparam
!     CHARACTER ( LEN=10 ) :: solvent_name         !name of solvent molecule
!     INTEGER :: nram                              !atoms per solvent molecule
!     TYPE ( solvent_atom ), DIMENSION (:), POINTER :: solvent
!     INTEGER :: ncons                             !solvent constraints
!     INTEGER, DIMENSION (:,:), POINTER :: IJCON   !constraints
!     REAL ( dbl ), DIMENSION (:), POINTER :: cval !constraint value
!  END TYPE gromos_topology_type
!
!  TYPE topology_type
!    LOGICAL :: active
!    INTEGER :: type
!    TYPE ( gromos_topology_type ), POINTER :: gromos
!  END TYPE topology_type
!
!  INTEGER, PARAMETER :: GROMOS_TOPOLOGY = 1001
!
!  TYPE ( topology_type ) :: topo_file
!
!  PRIVATE
!  PUBLIC :: read_topology_gromos, get_molecules_topology, &
!     set_molecule_type_topology, get_molecule_structure_topology
!
!!!*****
!!******************************************************************************
!
!CONTAINS
!
!!******************************************************************************
!!!****** topology/read_topology_file [1.0] *
!!!
!!!   NAME
!!!     read_topology_file
!!!
!!!   FUNCTION
!!!     Read topology file and store the full information in local
!!!     data types
!!!
!!!   AUTHOR
!!!     JGH 09.12.2001
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE read_topology_gromos ( filename, globenv )
!
!
!    CHARACTER(LEN=*), INTENT(IN)             :: filename
!    TYPE(global_environment_type), &
!      INTENT(IN)                             :: globenv
!
!    CHARACTER(LEN=80)                        :: line
!    INTEGER                                  :: ierr, iunit
!    LOGICAL                                  :: found
!
!!------------------------------------------------------------------------------
!! mark topology type as inactive
!
!  topo_file % active = .FALSE.
!
!! Is there a topology file
!  IF ( globenv % ionode ) THEN
!    INQUIRE ( FILE = filename, EXIST = found )
!    IF (found) THEN
!      iunit = get_unit()
!      OPEN ( unit=iunit, file=filename, form="FORMATTED" )
!      READ ( iunit, '(A80)' ) line
!      IF ( line(1:22) == "# GROMOS TOPOLOGY FILE" ) THEN
!        topo_file % active = .TRUE.
!        topo_file % type = GROMOS_TOPOLOGY
!      ELSE
!        CALL stop_program ( " read_topology_file ", &
!                            "Unknown type of topology file" )
!      END IF
!      CLOSE ( iunit )
!    END IF
!  END IF
!! send the outcome of the inquiry to all cpus
!  CALL mp_bcast(topo_file % active, globenv%source, globenv%group)
!  CALL mp_bcast(topo_file % type, globenv%source, globenv%group)
!
!  IF ( topo_file % active ) THEN
!    SELECT CASE (topo_file % type)
!    CASE DEFAULT
!      CALL stop_program ( " read_topology_file ", &
!                          "Unknown type of topology file" )
!    CASE (GROMOS_TOPOLOGY)
!      IF ( globenv % ionode .AND. globenv % print_level > 0 ) THEN
!        WRITE ( globenv % scr , '(A,T61,A20)') & 
!            " Parsing GROMOS Topology File ", TRIM(filename)
!      END IF
!      ALLOCATE ( topo_file % gromos, STAT=ierr )
!      IF ( ierr /= 0 ) CALL stop_program ( "read_topology_file",&
!         "Allocation of gromos type failed" )
!      CALL read_gromos_topology ( filename, topo_file % gromos, globenv )
!    END SELECT
!  END IF
!
!END SUBROUTINE read_topology_gromos
!
!!!*****
!!******************************************************************************
!!!****** topology/read_topology_file [1.0] *
!!!
!!!   NAME
!!!     read_gromos_topology
!!!
!!!   FUNCTION
!!!     Read GROMOS topology file
!!!
!!!   AUTHOR
!!!     JGH 20.01.2002
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE read_gromos_topology ( filename, gromos, globenv )
!
!
!    CHARACTER(LEN=*), INTENT(IN)             :: filename
!    TYPE(gromos_topology_type)               :: gromos
!    TYPE(global_environment_type), &
!      INTENT(IN)                             :: globenv
!
!    CHARACTER(LEN=10)                        :: cnum
!    CHARACTER(LEN=132)                       :: line
!    CHARACTER(LEN=30)                        :: label
!    INTEGER                                  :: i, ierror, ilen, io, j
!
!!------------------------------------------------------------------------------
!
!  CALL parser_init ( filename, globenv, comment="#", section_end="END" )
!
!! title section
!  label="TITLE"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  gromos%title = " "
!  DO i = 1, 5
!    CALL read_line ( line )
!    IF ( line(1:3)=="END" ) EXIT
!    gromos%title(i) = line
!  END DO
!! Phys. constants
!  label="TOPPHYSCON"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%frepsi )
!  CALL read_line
!  CALL get_next ( gromos%hbar )
!! Version
!  label="TOPVERSION"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  ilen=0
!  call get_next(gromos%version,ilen)
!! Atom type names
!  label="ATOMTYPENAME"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nratt )
!  ALLOCATE ( gromos%atomtypename(gromos%nratt), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology",&
!                                 "atomtypename",gromos%nratt)
!  i = 0
!  CALL read_line
!  DO WHILE ( test_next() /= 'X' )
!     i=i+1
!     IF ( i > gromos%nratt ) &
!        CALL stop_program ( "read_gromos_topology", "Too many atoms" )
!     ilen=0
!     CALL get_next ( gromos%atomtypename(i), ilen )
!     CALL read_line
!  END DO
!! Residue names
!  label="RESNAME"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nraa2 )
!  ALLOCATE ( gromos%residue_name(gromos%nraa2), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology",&
!                                 "residue_name",gromos%nraa2)
!  i = 0
!  CALL read_line
!  DO WHILE ( test_next() /= 'X' )
!     i=i+1
!     IF ( i > gromos%nraa2 ) &
!        CALL stop_program ( "read_gromos_topology", "Too many residues" )
!     ilen=0
!     CALL get_next ( gromos%residue_name(i), ilen )
!     CALL integer_to_string(i,cnum)
!     gromos%residue_name(i)=TRIM(gromos%residue_name(i))//TRIM(cnum)
!     CALL read_line
!  END DO
!! Solute atoms
!  label="SOLUTEATOM"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nrp )
!  ALLOCATE ( gromos%solute(gromos%nrp), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology",&
!                                 "solute",gromos%nrp)
!  i = 0
!  CALL read_line
!  DO WHILE ( test_next() /= 'X' )
!     i=i+1
!     IF ( i > gromos%nrp ) &
!        CALL stop_program ( "read_gromos_topology", "Too many solute atoms" )
!     ilen=0
!     CALL get_next ( gromos%solute(i)%atnm )
!     CALL get_next ( gromos%solute(i)%mres )
!     CALL get_next ( gromos%solute(i)%panm, ilen )
!     CALL get_next ( gromos%solute(i)%iac )
!     CALL get_next ( gromos%solute(i)%mass )
!     CALL get_next ( gromos%solute(i)%cg )
!     CALL get_next ( gromos%solute(i)%cgc )
!     CALL get_next ( gromos%solute(i)%ine )
!     ilen = gromos%solute(i)%ine
!     ALLOCATE ( gromos%solute(i)%excl(ilen), STAT=ierror )
!     IF(ierror/=0) CALL stop_memory("read_gromos_topology",&
!                                 "solute%excl",ilen)
!     DO j=1,ilen
!       CALL get_next ( gromos%solute(i)%excl(j) )
!     END DO
!     CALL read_line
!     CALL get_next ( gromos%solute(i)%ine14 )
!     ilen = gromos%solute(i)%ine14
!     ALLOCATE ( gromos%solute(i)%int14(ilen), STAT=ierror )
!     IF(ierror/=0) CALL stop_memory("read_gromos_topology",&
!                                 "solute%int14",ilen)
!     DO j=1,ilen
!       CALL get_next ( gromos%solute(i)%int14(j) )
!     END DO
!     CALL read_line
!  END DO
!! Bond types
!  label="BONDTYPE"
!  CALL find_keyword ( label, filename, ierror )
!!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nbty )
!  ALLOCATE ( gromos%cb(gromos%nbty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","cb",gromos%nbty)
!  ALLOCATE ( gromos%b0(gromos%nbty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","b0",gromos%nbty)
!  DO i = 1, gromos%nbty
!    CALL read_line
!    CALL get_next ( gromos%cb(i) )
!    CALL get_next ( gromos%b0(i) )
!  END DO
!! Bonds with H atoms
!  label="BONDH"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nbonh )
!  ALLOCATE ( gromos%ijbh(gromos%nbonh,3), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ijbh",3*gromos%nbonh)
!  DO i = 1, gromos%nbonh
!    CALL read_line
!    CALL get_next ( gromos%ijbh(i,1) )
!    CALL get_next ( gromos%ijbh(i,2) )
!    CALL get_next ( gromos%ijbh(i,3) )
!  END DO
!! Bonds without H atoms 
!  label="BOND "
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nbon )
!  ALLOCATE ( gromos%ijbo(gromos%nbon,3), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ijbh",3*gromos%nbonh)
!  DO i = 1, gromos%nbon
!    CALL read_line
!    CALL get_next ( gromos%ijbo(i,1) )
!    CALL get_next ( gromos%ijbo(i,2) )
!    CALL get_next ( gromos%ijbo(i,3) )
!  END DO
!! Bond angle types
!  label="BONDANGLETYPE"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%ntty )
!  ALLOCATE ( gromos%ct(gromos%ntty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ct",gromos%ntty)
!  ALLOCATE ( gromos%t0(gromos%ntty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","t0",gromos%ntty)
!  DO i = 1, gromos%ntty
!    CALL read_line
!    CALL get_next ( gromos%ct(i) )
!    CALL get_next ( gromos%t0(i) )
!  END DO
!! Bond angels with H
!  label="BONDANGLEH"
!  CALL search_label ( label, ierror, ignore_case=.TRUE. )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%ntheh )
!  ALLOCATE ( gromos%ijkth(gromos%ntheh,4), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ijkth",4*gromos%ntheh)
!  DO i = 1, gromos%ntheh
!    CALL read_line
!    CALL get_next ( gromos%ijkth(i,1) )
!    CALL get_next ( gromos%ijkth(i,2) )
!    CALL get_next ( gromos%ijkth(i,3) )
!    CALL get_next ( gromos%ijkth(i,4) )
!  END DO
!! Bond angels without H
!  label="BONDANGLE"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nthe )
!  ALLOCATE ( gromos%ijkt(gromos%nthe,4), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ijkt",4*gromos%nthe)
!  DO i = 1, gromos%nthe
!    CALL read_line
!    CALL get_next ( gromos%ijkt(i,1) )
!    CALL get_next ( gromos%ijkt(i,2) )
!    CALL get_next ( gromos%ijkt(i,3) )
!    CALL get_next ( gromos%ijkt(i,4) )
!  END DO
!! improper dihedral types
!  label="IMPDIHEDRALTYPE"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nqty )
!  ALLOCATE ( gromos%cq(gromos%nqty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","cq",gromos%nqty)
!  ALLOCATE ( gromos%q0(gromos%nqty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","q0",gromos%nqty)
!  DO i = 1, gromos%nqty
!    CALL read_line
!    CALL get_next ( gromos%cq(i) )
!    CALL get_next ( gromos%q0(i) )
!  END DO
!! improper dihedrals with H
!  label="IMPDIHEDRALH"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nqhih )
!  ALLOCATE ( gromos%ijkqh(gromos%nqhih,5), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ijkqh",5*gromos%nqhih)
!  DO i = 1, gromos%nqhih
!    CALL read_line
!    CALL get_next ( gromos%ijkqh(i,1) )
!    CALL get_next ( gromos%ijkqh(i,2) )
!    CALL get_next ( gromos%ijkqh(i,3) )
!    CALL get_next ( gromos%ijkqh(i,4) )
!    CALL get_next ( gromos%ijkqh(i,5) )
!  END DO
!! improper dihedrals without H
!  label="IMPDIHEDRAL"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nqhi )
!  ALLOCATE ( gromos%ijkq(gromos%nqhi,5), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ijkq",5*gromos%nqhi)
!  DO i = 1, gromos%nqhi
!    CALL read_line
!    CALL get_next ( gromos%ijkq(i,1) )
!    CALL get_next ( gromos%ijkq(i,2) )
!    CALL get_next ( gromos%ijkq(i,3) )
!    CALL get_next ( gromos%ijkq(i,4) )
!    CALL get_next ( gromos%ijkq(i,5) )
!  END DO
!! dihedral types
!  label="DIHEDRALTYPE"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%npty )
!  ALLOCATE ( gromos%cp(gromos%npty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","cp",gromos%npty)
!  ALLOCATE ( gromos%pd(gromos%npty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","pd",gromos%npty)
!  ALLOCATE ( gromos%np(gromos%npty), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","np",gromos%npty)
!  DO i = 1, gromos%npty
!    CALL read_line
!    CALL get_next ( gromos%cp(i) )
!    CALL get_next ( gromos%pd(i) )
!    CALL get_next ( gromos%np(i) )
!  END DO
!! dihedrals with H
!  label="DIHEDRALH"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nphih )
!  ALLOCATE ( gromos%ijkph(gromos%nphih,5), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ijkph",5*gromos%nphih)
!  DO i = 1, gromos%nphih
!    CALL read_line
!    CALL get_next ( gromos%ijkph(i,1) )
!    CALL get_next ( gromos%ijkph(i,2) )
!    CALL get_next ( gromos%ijkph(i,3) )
!    CALL get_next ( gromos%ijkph(i,4) )
!    CALL get_next ( gromos%ijkph(i,5) )
!  END DO
!! dihedrals without H
!  label="DIHEDRAL"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nphi )
!  ALLOCATE ( gromos%ijkp(gromos%nphi,5), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","ijkp",5*gromos%nphi)
!  DO i = 1, gromos%nphi
!    CALL read_line
!    CALL get_next ( gromos%ijkp(i,1) )
!    CALL get_next ( gromos%ijkp(i,2) )
!    CALL get_next ( gromos%ijkp(i,3) )
!    CALL get_next ( gromos%ijkp(i,4) )
!    CALL get_next ( gromos%ijkp(i,5) )
!  END DO
!! Lennard-Jones Parameters
!  label="LJPARAMETERS"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nratt2 )
!  ALLOCATE ( gromos%ljparam(gromos%nratt2), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology",&
!                                 "ljparam",gromos%nratt2)
!  DO i = 1, gromos%nratt2
!    CALL read_line
!    CALL get_next ( gromos%ljparam(i)%iac )
!    CALL get_next ( gromos%ljparam(i)%jac )
!    CALL get_next ( gromos%ljparam(i)%c12 )
!    CALL get_next ( gromos%ljparam(i)%c6  )
!    CALL get_next ( gromos%ljparam(i)%cs12)
!    CALL get_next ( gromos%ljparam(i)%cs6 )
!  END DO
!! Solvent atoms
!  label="SOLVENTATOM"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%nram )
!  ALLOCATE ( gromos%solvent(gromos%nram), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology",&
!                                 "solvent",gromos%nram)
!  DO i = 1, gromos%nram
!    CALL read_line
!    CALL get_next ( gromos%solvent(i)%i )
!    ilen=0
!    CALL get_next ( gromos%solvent(i)%anms, ilen )
!    CALL get_next ( gromos%solvent(i)%iacs )
!    CALL get_next ( gromos%solvent(i)%mass )
!    CALL get_next ( gromos%solvent(i)%cgs )
!  END DO
!! Solvent constraints
!  label="SOLVENTCONSTR"
!  CALL find_keyword ( label, filename, ierror )
!  IF(ierror/=0) CALL stop_program ( "read_gromos_topology", label)
!  CALL read_line
!  CALL get_next ( gromos%ncons )
!  ALLOCATE ( gromos%ijcon(gromos%ncons,2), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology",&
!                                 "ijcons",2*gromos%ncons)
!  ALLOCATE ( gromos%cval(gromos%ncons), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","cval",gromos%ncons)
!  DO i = 1, gromos%ncons
!    CALL read_line
!    CALL get_next ( gromos%ijcon(i,1) )
!    CALL get_next ( gromos%ijcon(i,2) )
!    CALL get_next ( gromos%cval(i) )
!  END DO
!
!  CALL parser_end ()
!
!! count number of atoms per residue
!  ALLOCATE ( gromos%resat(2,gromos%nraa2), STAT=ierror )
!  IF(ierror/=0) CALL stop_memory("read_gromos_topology","resat",2*gromos%nraa2)
!  gromos%resat(:,:) = 0
!  DO i = 1, gromos%nrp
!     j = gromos%solute(i)%mres
!     IF ( gromos%resat(1,j) == 0 ) gromos%resat(1,j) = i
!     gromos%resat(2,j) = i
!  END DO
!
!! write the data to the output file
!  IF ( globenv % ionode .AND. globenv % print_level > 9 ) THEN
!    io = globenv % scr
!    DO i=1,5
!      IF(LEN_TRIM(gromos%title(i)) > 0 ) &
!         WRITE ( io, '(T2,A)' ) ADJUSTL(gromos%title(i)(1:79))
!    END DO
!    WRITE ( io, '(A,T61,E20.10)' ) " Permittivity of vacuum EPS0", gromos%frepsi
!    WRITE ( io, '(A,T61,E20.10)' ) " Planck's constant", gromos%hbar
!    WRITE ( io, '(A,T71,A10)' ) " Version of Topology file",&
!                                ADJUSTR(gromos%version)
!    WRITE ( io, '(A,T71,I10)' ) " Number of Atom Types",gromos%nratt
!    DO i=1,gromos%nratt
!       WRITE ( io, '(T71,A10)' ) ADJUSTR(gromos%atomtypename(i))
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of Residue Types",gromos%nraa2
!    DO i=1,gromos%nraa2
!       WRITE ( io, '(T71,A10)' ) ADJUSTR(gromos%residue_name(i))
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of Solute Atoms",gromos%nrp
!    WRITE ( io, '(A,A)' ) " ATNM  MRES  PANM  IAC      MASS          CG   CGC ",&
!                          "        INE/INE14 "
!    DO i = 1, gromos%nrp
!      WRITE ( io, '(i5,i6,A6,i5,F10.5,f12.6,i6,i6,10(i4))') &
!         gromos%solute(i)%atnm,gromos%solute(i)%mres,&
!         ADJUSTR(gromos%solute(i)%panm),gromos%solute(i)%iac,&
!         gromos%solute(i)%mass,gromos%solute(i)%cg,&
!         gromos%solute(i)%cgc,gromos%solute(i)%ine,&
!         (gromos%solute(i)%excl(j),j=1,gromos%solute(i)%ine)
!      WRITE ( io, '(T51,i6,10(i4))') gromos%solute(i)%ine14, &
!         (gromos%solute(i)%int14(j),j=1,gromos%solute(i)%ine14)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of Bond types",gromos%nbty
!    WRITE ( io, '(T47,A,T70,A)') "Force constant","Bond length"
!    DO i = 1, gromos%nbty
!      WRITE ( io, '(T41,E20.10,T71,F10.7)') gromos%cb(i),gromos%b0(i)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) &
!      " Number of bonds involving H atoms in solute", gromos%nbonh
!    WRITE ( io, '(A)' ) "    IBH    JBH    ICBH "
!    DO i = 1, gromos%nbonh
!       WRITE ( io, '(i7,i7,i8)') (gromos%ijbh(i,j),j=1,3)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) &
!      " Number of bonds not involving H atoms in solute", gromos%nbon
!    WRITE ( io, '(A)' ) "     IB     JB     ICB "
!    DO i = 1, gromos%nbon
!       WRITE ( io, '(i7,i7,i8)') (gromos%ijbo(i,j),j=1,3)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of Bond Angle Types",gromos%ntty
!    WRITE ( io, '(T47,A,T76,A)') "Force constant","Angle"
!    DO i = 1, gromos%ntty
!      WRITE ( io, '(T41,E20.10,T66,F15.7)') gromos%ct(i),gromos%t0(i)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) &
!      " Number of bond angles involving H atoms in solute", gromos%ntheh
!    WRITE ( io, '(A)' ) "    ITH    JTH    KTH    ICTH "
!    DO i = 1, gromos%ntheh
!       WRITE ( io, '(3(i7),i8)') (gromos%ijkth(i,j),j=1,4)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) &
!      " Number of bond angles not involving H atoms in solute", gromos%nthe
!    WRITE ( io, '(A)' ) "     IT     JT     KT     ICT "
!    DO i = 1, gromos%nthe
!       WRITE ( io, '(3(i7),i8)') (gromos%ijkt(i,j),j=1,4)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of Improper Dihedral Types",gromos%nqty
!    WRITE ( io, '(T47,A,T64,A)') "Force constant","Improper dihedral"
!    DO i = 1, gromos%nqty
!      WRITE ( io, '(T41,E20.10,T71,F10.7)') gromos%cq(i),gromos%q0(i)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) &
!      " Number of improper dihedrals involving H atoms in solute", gromos%nqhih
!    WRITE ( io, '(A)' ) "    IQH    JQH    KQH    LQH    ICQH "
!    DO i = 1, gromos%nqhih
!       WRITE ( io, '(4(i7),i8)') (gromos%ijkqh(i,j),j=1,5)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) &
!      " Number of improper dihedrals not involving H atoms in solute", gromos%nqhi
!    WRITE ( io, '(A)' ) "     IQ     JQ     KQ     LQ     ICQ "
!    DO i = 1, gromos%nqhi
!       WRITE ( io, '(4(i7),i8)') (gromos%ijkq(i,j),j=1,5)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of Dihedral Types",gromos%npty
!    WRITE ( io, '(T21,A,T39,A,T69,A)') &
!      "Force constant","COS phase shift","multiplicity"
!    DO i = 1, gromos%npty
!      WRITE ( io, '(T15,E20.10,T40,F14.7,T71,I10)') &
!        gromos%cp(i),gromos%pd(i),gromos%np(i)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) &
!      " Number of dihedrals involving H atoms in solute", gromos%nphih
!    WRITE ( io, '(A)' ) "    IPH    JPH    KPH    LPH    ICPH "
!    DO i = 1, gromos%nphih
!       WRITE ( io, '(4(i7),i8)') (gromos%ijkph(i,j),j=1,5)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) &
!      " Number of improper dihedrals not involving H atoms in solute", gromos%nphi
!    WRITE ( io, '(A)' ) "     IP     JP     KP     LP     ICP "
!    DO i = 1, gromos%nphi
!       WRITE ( io, '(4(i7),i8)') (gromos%ijkp(i,j),j=1,5)
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of Lennard-Jones interaction types", &
!       gromos%nratt2
!    WRITE ( io, '(T2,A,T7,A,T25,A,T42,A,T61,A,T78,A)' ) &
!      " iac"," jac","C12","C6","CS12","CS6"
!    DO i = 1, gromos%nratt2
!      WRITE ( io, '(I5,I5,T15,E13.7,T31,E13.7,T52,E13.7,T68,E13.7)' ) &
!        gromos%ljparam(i)%iac,gromos%ljparam(i)%jac,gromos%ljparam(i)%c12,&
!        gromos%ljparam(i)%c6,gromos%ljparam(i)%cs12,gromos%ljparam(i)%cs6
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of solvent atoms", gromos%nram
!    WRITE ( io, '(T19,A,T37,A,T43,A,T63,A,T75,A)' ) &
!      " Number","Name","vdW type","Mass","Charge"
!    DO i = 1, gromos%nram
!      WRITE ( io, '(T20,I6,T30,A10,T45,I6,T55,F12.6,T71,F10.5)' ) &
!        gromos%solvent(i)%i,ADJUSTR(gromos%solvent(i)%anms),&
!        gromos%solvent(i)%iacs,gromos%solvent(i)%mass,gromos%solvent(i)%cgs
!    END DO
!    WRITE ( io, '(A,T71,I10)' ) " Number of solvent constraints", gromos%ncons
!    WRITE ( io, '(T40,A,T50,A,T76,A)' ) " icons"," jcons"," cons"
!    DO i = 1, gromos%ncons
!      WRITE ( io, '(T40,I6,T50,I6,T66,F15.6)' ) &
!        gromos%ijcon(i,1),gromos%ijcon(i,2),gromos%cval(i)
!    END DO
!
!  END IF
!
!END SUBROUTINE read_gromos_topology
!
!SUBROUTINE find_keyword ( label, filename, ierr )
!    CHARACTER(LEN=*)                         :: label, filename
!    INTEGER                                  :: ierr
!
!    CHARACTER(len=132)                       :: line
!    INTEGER                                  :: ios, llen, lu
!
!  INQUIRE ( FILE=filename, NUMBER=lu )
!  REWIND ( lu )
!
!  ierr = 0
!  llen = LEN_TRIM ( label )
!  DO
!     READ ( lu, iostat = ios, fmt = '( a )' ) line
!     IF ( ios > 0 ) THEN
!        CALL stop_program ( 'find_keyword','read error')
!     ELSE IF ( ios < 0 ) THEN
!        ierr = 2
!        EXIT
!     ELSE
!        IF ( line(1:llen)==label(1:llen) .AND. &
!             line(llen+1:llen+1)==" " ) EXIT
!     END IF
!  END DO
!
!END SUBROUTINE find_keyword
!
!!!*****
!!******************************************************************************
!!!****** topology/get_molecules_topology [1.0] *
!!!
!!!   NAME
!!!     get_molecules_topology
!!!
!!!   FUNCTION
!!!     Counts number of molecule definitions on topology file
!!!
!!!   AUTHOR
!!!     JGH 22.01.2002
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   SOURCE
!!******************************************************************************
!
!FUNCTION get_molecules_topology ( ) RESULT ( nmol )
!
!
!    INTEGER                                  :: nmol
!
!!------------------------------------------------------------------------------
!
!  IF ( topo_file % active ) THEN
!    SELECT CASE (topo_file % type)
!    CASE DEFAULT
!      CALL stop_program ( " read_topology_file ", &
!                          "Unknown type of topology file" )
!    CASE (GROMOS_TOPOLOGY)
!      nmol = topo_file % gromos % nraa2                      ! residues
!      IF ( topo_file % gromos % nram > 0 ) nmol = nmol + 1   ! solvent
!    END SELECT
!  ELSE
!    nmol = 0
!  ENDIF
!
!END FUNCTION get_molecules_topology 
!
!!!*****
!!******************************************************************************
!!!****** topology/set_molecule_type_topology [1.0] *
!!!
!!!   NAME
!!!     set_molecule_type_topology
!!!
!!!   FUNCTION
!!!     Transfers the setup section information to the mol_setup structure
!!!
!!!   AUTHOR
!!!     JGH 22.01.2002
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE set_molecule_type_topology ( mol_setup, default_thermostat )
!
!
!    TYPE(molecule_type), DIMENSION(:), &
!      INTENT(INOUT)                          :: mol_setup
!    CHARACTER(LEN=*), INTENT(IN)             :: default_thermostat
!
!    CHARACTER(LEN=10)                        :: molname
!    INTEGER                                  :: imol, inserted, istat, j, &
!                                                nlen, nmol
!    LOGICAL, ALLOCATABLE, DIMENSION(:)       :: found
!
!!------------------------------------------------------------------------------
!
!  IF ( topo_file % active ) THEN
!    nmol = SIZE ( mol_setup )
!    ALLOCATE ( found(nmol), STAT=istat )
!    IF ( istat /= 0 ) call stop_memory ( "set_mol_type_topology", &
!       "found", nmol )
!    found = .FALSE.
!    SELECT CASE (topo_file % type)
!    CASE DEFAULT
!      CALL stop_program ( " read_topology_file ", &
!                          "Unknown type of topology file" )
!    CASE (GROMOS_TOPOLOGY)
!      DO imol = 1, topo_file % gromos % nraa2
!        molname = topo_file % gromos % residue_name(imol)
!        inserted = 0
!        DO j = 1, nmol
!          IF ( found ( j ) ) CYCLE
!          nlen = LEN_TRIM ( mol_setup(j) % molpar % name )
!          IF (  mol_setup(j) % molpar % name(1:1) == "#" ) THEN
!            mol_setup(j) % molpar % name = molname
!            mol_setup(j) % num_mol = 1
!            mol_setup(j) % nhcopt = default_thermostat
!            inserted = j
!            EXIT
!          ELSE IF (TRIM(molname)==TRIM(mol_setup(j) % molpar % name)) THEN
!            inserted = j
!            EXIT
!          END IF 
!        END DO
!        IF ( inserted == 0 ) THEN
!          CALL stop_program ( "set_mol_type_topology", &
!            "Could not find molecule  "//molname )
!        ELSE
!          found ( j ) = .TRUE.
!        END IF
!      END DO
!      IF ( topo_file % gromos % nram > 0 ) THEN
!        ! look for the solvent molecule
!        inserted = 0
!        ! run in reverse order to put solvent molecules at the end
!        DO j = nmol, 1, -1
!          IF ( found ( j ) ) CYCLE
!          IF (  mol_setup(j) % molpar % name(1:1) == "#" ) THEN
!            mol_setup(j) % molpar % name = "SOLVENT"
!            mol_setup(j) % num_mol = 1
!            mol_setup(j) % nhcopt = default_thermostat
!            inserted = j
!            EXIT
!          ELSE 
!            ! we give the solvent the name of the first molecule not
!            ! yet used
!            inserted = j
!            topo_file % gromos % solvent_name = mol_setup(j) % molpar % name
!            EXIT
!          END IF 
!        END DO
!        IF ( inserted == 0 ) THEN
!          CALL stop_program ( "set_mol_type_topology", &
!            "Could not find solvent molecule" )
!        ELSE
!          found ( j ) = .TRUE.
!        END IF
!      END IF
!      IF ( .NOT. ALL ( found ) ) THEN
!        CALL stop_program ( "set_mol_type_topology", &
!            "Could not find all molecules" )
!      END IF
!    END SELECT
!    DEALLOCATE ( found, STAT=istat )
!    IF ( istat /= 0 ) call stop_memory ( "set_mol_type_topology", "found" )
!  ENDIF
!
!END SUBROUTINE set_molecule_type_topology
!
!!!*****
!!******************************************************************************
!!!****** topology/get_molecule_structure_topology [1.0] *
!!!
!!!   NAME
!!!     get_molecule_structure_topology
!!!
!!!   FUNCTION
!!!     Transfers the molecule section information to the mol_setup structure
!!!
!!!   AUTHOR
!!!     JGH 25.01.2002
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE get_molecule_structure_topology ( molname, mol_setup )
!
!
!    CHARACTER(LEN=*), DIMENSION(:), &
!      INTENT(IN)                             :: molname
!    TYPE(molecule_type), DIMENSION(:), &
!      INTENT(INOUT)                          :: mol_setup
!
!    CHARACTER(LEN=10)                        :: string
!    INTEGER                                  :: i, imol, ios, iraa, j, na, &
!                                                nmol_type
!    LOGICAL                                  :: solvent
!
!!------------------------------------------------------------------------------
!
!  IF ( topo_file % active ) THEN
!
!    nmol_type = SIZE ( mol_setup )
!
!    SELECT CASE (topo_file % type)
!    CASE DEFAULT
!      CALL stop_program ( " get_molecule_structure_topology ", &
!                          "Unknown type of topology file" )
!    CASE (GROMOS_TOPOLOGY)
!
!      DO iraa = 1, topo_file % gromos % nraa2 + 1
!
!        IF ( iraa <= topo_file % gromos % nraa2 ) THEN
!          string = topo_file % gromos % residue_name ( iraa )
!          solvent = .FALSE.
!        ELSE
!          string = topo_file % gromos % solvent_name
!          solvent = .TRUE.
!        END IF
!        imol = str_search(molname,nmol_type,string)
!        IF (imol<=0) CALL stop_program ( "get_molecule_structure_topology",&
!                      "unknown molecule type "//string )
!
!        IF ( solvent ) THEN
!           na = topo_file % gromos % nram
!        ELSE
!           na = topo_file % gromos % resat(2,iraa) - &
!                topo_file % gromos % resat(1,iraa) + 1
!        END IF
!        mol_setup(imol) %molpar%natom = na
!        ALLOCATE (mol_setup(imol)%molpar%aname(na),STAT=ios)
!        IF ( ios /= 0 ) CALL stop_memory &
!             ( 'read_molecule_section', 'aname', na )
!        ALLOCATE (mol_setup(imol)%molpar%alabel(na),STAT=ios)
!        IF ( ios /= 0 ) CALL stop_memory &
!             ( 'read_molecule_section', 'alabel', na )
!        ALLOCATE (mol_setup(imol)%molpar%aweight(na),STAT=ios)
!        IF ( ios /= 0 ) CALL stop_memory &
!             ( 'read_molecule_section', 'aweight', na )
!        ALLOCATE (mol_setup(imol)%molpar%acharge(na),STAT=ios)
!        IF ( ios /= 0 ) CALL stop_memory &
!             ( 'read_molecule_section', 'acharge', na )
!
!        mol_setup(imol) %molpar%natom_type = 0
!        mol_setup(imol) %molpar%nbonds = 0
!        mol_setup(imol) %molpar%nbends = 0
!        mol_setup(imol) %molpar%ntorsions = 0
!        mol_setup(imol) %molpar%ndcon = 0
!        mol_setup(imol) %molpar%n3x3con = 0
!
!        ! atoms
!        IF ( solvent ) THEN
!           DO i = 1, na
!             mol_setup(imol)%molpar%aname(i) = &
!                               topo_file % gromos % solvent(i) % anms
!             CALL integer_to_string ( topo_file % gromos % solvent(i) % iacs, &
!                                      mol_setup(imol)%molpar%alabel(i) )
!             mol_setup(imol)%molpar%aweight(i) = &
!                               topo_file % gromos % solvent(i) % mass
!             mol_setup(imol)%molpar%acharge(i) = &
!                               topo_file % gromos % solvent(i) % cgs
!           END DO
!!
!! bonds bends torsions in solvent molecule ??????
!!
!           mol_setup(imol) %molpar%ndcon = topo_file % gromos % ncons
!           na = mol_setup(imol) %molpar%ndcon
!           IF (na>=0) THEN
!              ALLOCATE (mol_setup(imol)%molpar%ndc(3,na),STAT=ios)
!              IF ( ios /= 0 ) CALL stop_memory &
!                 ( 'get_molecule_structure_topology', 'ndc', 3 * na )
!              DO i = 1, na
!                mol_setup(imol) %molpar%ndc(1,i) = topo_file % gromos % ijcon (i,1)
!                mol_setup(imol) %molpar%ndc(2,i) = topo_file % gromos % ijcon (i,2)
!                mol_setup(imol) %molpar%ndc(3,i) = i
!              END DO
!           END IF
!
!        ELSE
!
!           DO i = 1, na
!             j = topo_file % gromos % resat(1,iraa) + i - 1
!             mol_setup(imol)%molpar%aname(i) = &
!                               topo_file % gromos % solute(j) % panm
!             CALL integer_to_string ( topo_file % gromos % solute(j) % iac, &
!                                      mol_setup(imol)%molpar%alabel(i) )
!             mol_setup(imol)%molpar%aweight(i) = &
!                               topo_file % gromos % solute(j) % mass
!             mol_setup(imol)%molpar%acharge(i) = &
!                               topo_file % gromos % solute(j) % cg
!           END DO
!           !bonds
!           ! count number of bonds within molecule
!!          na = mol_setup(imol) %molpar%nbonds
!!          IF (na>=0) THEN
!!             ALLOCATE (mol_setup(imol)%molpar%bonds(3,na),STAT=ios)
!!             IF (ios /= 0 ) CALL stop_memory &
!!                ( 'get_molecule_structure_topology', 'bonds', 3 * na )
!!             DO i = 1, na
!!                mol_setup(imol) %molpar%bonds(1,i) 
!!                mol_setup(imol) %molpar%bonds(2,i) 
!!                mol_setup(imol) %molpar%bonds(3,i) 
!!             END DO
!!          END IF
!           !bends
!           !count number of bends within molecule
!!          na = mol_setup(imol) %molpar%nbends
!!          IF (na>=0) THEN
!!             ALLOCATE (mol_setup(imol)%molpar%bends(4,na),STAT=ios)
!!             IF (ios /= 0 ) CALL stop_memory &
!!                ( 'get_molecule_structure_topology', 'bends', 4 * na )
!!             DO i = 1, na
!!                mol_setup(imol) %molpar%bends(1,i) 
!!                mol_setup(imol) %molpar%bends(2,i) 
!!                mol_setup(imol) %molpar%bends(3,i) 
!!                mol_setup(imol) %molpar%bends(4,i) 
!!             END DO
!!          END IF
!           !torsion
!           !count number of torsions
!!          na = mol_setup(imol) %molpar%ntorsions
!!          IF (na>=0) THEN
!!             ALLOCATE ( mol_setup(imol)%molpar%torsions(5,na), STAT=ios )
!!             IF ( ios /= 0 ) CALL stop_memory &
!!                  ( 'get_molecule_structure_topology', 'torsions', 5 * na )
!!             DO i = 1, na
!!                mol_setup(imol) %molpar%torsions(1,i) 
!!                mol_setup(imol) %molpar%torsions(2,i) 
!!                mol_setup(imol) %molpar%torsions(3,i) 
!!                mol_setup(imol) %molpar%torsions(4,i) 
!!             END DO
!!          END IF
!           !distance constraints
!!          na = mol_setup(imol) %molpar%ndcon
!!          IF (na>=0) THEN
!!             ALLOCATE (mol_setup(imol)%molpar%ndc(3,na),STAT=ios)
!!             IF ( ios /= 0 ) CALL stop_memory &
!!                ( 'get_molecule_structure_topology', 'ndc', 3 * na )
!!             DO i = 1, na
!!               mol_setup(imol) %molpar%ndc(1,i) =
!!               mol_setup(imol) %molpar%ndc(2,i) =
!!               mol_setup(imol) %molpar%ndc(3,i) =
!!             END DO
!!          END IF
!
!        END IF
!        
!      END DO
!
!    END SELECT
!
!  ENDIF
!
!END SUBROUTINE get_molecule_structure_topology
!
!!*****
!******************************************************************************

END MODULE topology_gromos

!******************************************************************************
