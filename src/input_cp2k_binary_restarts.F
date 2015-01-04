!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to read the binary restart file of CP2K
!> \author Matthias Krack (MK)
!> \par History
!>      - Creation (17.02.2011,MK)
!> \version 1.0
! *****************************************************************************
MODULE input_cp2k_binary_restarts

  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr,&
                                             debug_print_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE extended_system_types,           ONLY: lnhc_parameters_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE string_table,                    ONLY: id2str,&
                                             s2s,&
                                             str2id
  USE termination,                     ONLY: print_message,&
                                             stop_memory,&
                                             stop_program,&
                                             stop_read
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology_types,                  ONLY: atom_info_type,&
                                             topology_parameters_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_binary_restarts'

  PUBLIC :: read_binary_coordinates,&
            read_binary_cs_coordinates,&
            read_binary_thermostats_nose,&
            read_binary_velocities

CONTAINS

! *****************************************************************************
!> \brief   Read the input section &COORD from an external file written in
!>          binary format.
!> \param topology ...
!> \param root_section ...
!> \param para_env ...
!> \param subsys_section ...
!> \param binary_file_read ...
!> \param error ...
!> \par History
!>      - Creation (10.02.2011,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_binary_coordinates(topology,root_section,para_env,&
                                     subsys_section,binary_file_read,error)

    TYPE(topology_parameters_type)           :: topology
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    LOGICAL, INTENT(OUT)                     :: binary_file_read
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_binary_coordinates', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: binary_restart_file_name
    CHARACTER(LEN=default_string_length)     :: string
    INTEGER :: handle, iatom, ikind, input_unit, istat, iw, natom, natomkind, &
      ncore, nmolecule, nmoleculekind, nshell
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ibuf, id_name
    LOGICAL                                  :: failure
    TYPE(atom_info_type), POINTER            :: atom_info
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    NULLIFY (logger)
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(root_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subsys_section),cp_failure_level,routineP,error,failure)
    logger => cp_error_get_logger(error)

    binary_file_read = .FALSE.

    CALL section_vals_val_get(root_section,"EXT_RESTART%BINARY_RESTART_FILE_NAME",&
                              c_val=binary_restart_file_name,error=error)

    IF (TRIM(ADJUSTL(binary_restart_file_name)) == "") THEN
       CALL timestop(handle)
       RETURN
    END IF

    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/XYZ_INFO",&
                              extension=".subsysLog",error=error)

    natomkind = 0
    natom = 0
    ncore = 0
    nshell = 0
    nmoleculekind = 0
    nmolecule = 0

    ! Open binary restart file and read number atomic kinds, atoms, etc.
    IF (para_env%ionode) THEN
       CALL open_file(file_name=binary_restart_file_name,&
                      file_status="OLD",&
                      file_form="UNFORMATTED",&
                      file_action="READWRITE",&
                      file_position="REWIND",&
                      unit_number=input_unit,&
                      debug=iw)
       READ (UNIT=input_unit,IOSTAT=istat)&
         natomkind,natom,ncore,nshell,nmoleculekind,nmolecule
       IF (istat /= 0) THEN
          CALL stop_read(routineN,moduleN,__LINE__,&
                         "natomkind,natom,ncore,nshell,nmoleculekind,nmolecule "//&
                         "(IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                         input_unit,para_env)
       END IF
       IF (iw > 0) THEN
          WRITE (UNIT=iw,FMT="(T2,A,T71,I10)")&
            "Number of atomic kinds:",natomkind,&
            "Number of atoms:",natom,&
            "Number of cores (only core-shell model):",ncore,&
            "Number of shells (only core-shell model):",nshell,&
            "Number of molecule kinds:",nmoleculekind,&
            "Number of molecules",nmolecule
       END IF
    END IF

    CALL mp_bcast(natomkind,para_env%source,para_env%group)
    CALL mp_bcast(natom,para_env%source,para_env%group)
    CALL mp_bcast(ncore,para_env%source,para_env%group)
    CALL mp_bcast(nshell,para_env%source,para_env%group)
    CALL mp_bcast(nmoleculekind,para_env%source,para_env%group)
    CALL mp_bcast(nmolecule,para_env%source,para_env%group)

    ALLOCATE (id_name(natomkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "id_name",int_size*natomkind,para_env)
    ! Read atomic kind names
    DO ikind=1,natomkind
       IF (para_env%ionode) THEN
          READ (UNIT=input_unit,IOSTAT=istat) string
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "string (IOSTAT = "//&
                                         TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                         input_unit,para_env)
       END IF
       CALL mp_bcast(string,para_env%source,para_env%group)
       id_name(ikind) = str2id(string)
    END DO

    ! Allocate and initialise atom_info array
    atom_info => topology%atom_info
    ALLOCATE (atom_info%id_molname(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%id_molname",&
                                     int_size*natom,para_env)
    atom_info%id_molname(:) = 0
    ALLOCATE (atom_info%id_resname(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%id_resname",&
                                     int_size*natom,para_env)
    atom_info%id_resname(:) = 0
    ALLOCATE (atom_info%resid(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%resid",&
                                     int_size*natom,para_env)
    atom_info%resid = 1
    ALLOCATE (atom_info%id_atmname(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%id_atmname",&
                                     int_size*natom,para_env)
    atom_info%id_atmname = 0
    ALLOCATE (atom_info%r(3,natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%r",&
                                     int_size*3*natom,para_env)
    atom_info%r(:,:) = 0.0_dp
    ALLOCATE (atom_info%atm_mass(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%atm_mass",&
                                     dp_size*natom,para_env)
    atom_info%atm_mass(:) = HUGE(0.0_dp)
    ALLOCATE (atom_info%atm_charge(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%atm_charge",&
                                     dp_size*natom,para_env)
    atom_info%atm_charge(:) = -HUGE(0.0_dp)
    ALLOCATE (atom_info%occup(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%occup",&
                                     dp_size*natom,para_env)
    atom_info%occup(:) = 0.0_dp
    ALLOCATE (atom_info%beta(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%id_beta",&
                                     dp_size*natom,para_env)
    atom_info%beta(:) = 0.0_dp
    ALLOCATE (atom_info%id_element(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_info%id_element",&
                                     int_size*natom,para_env)
    atom_info%id_element(:) = 0
    ALLOCATE (ibuf(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "ibuf",int_size*natom,para_env)

    ! Read atomic kind number of each atom
    IF (para_env%ionode) THEN
       READ (UNIT=input_unit,IOSTAT=istat) ibuf(1:natom)
       IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                      "ibuf (IOSTAT = "//&
                                      TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                      input_unit,para_env)
    END IF
    CALL mp_bcast(ibuf,para_env%source,para_env%group)
    DO iatom=1,natom
       ikind = ibuf(iatom)
       atom_info%id_atmname(iatom) = id_name(ikind)
       atom_info%id_element(iatom) = id_name(ikind)
    END DO
    DEALLOCATE (id_name,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "id_name",para_env)

    ! Read atomic coordinates
    IF (para_env%ionode) THEN
       READ (UNIT=input_unit,IOSTAT=istat) atom_info%r(1:3,1:natom)
       IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                      "atom_info%r(1:3,1:natom) (IOSTAT = "//&
                                      TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                      input_unit,para_env)
    END IF
    CALL mp_bcast(atom_info%r,para_env%source,para_env%group)

    ! Read molecule information if available
    IF (nmolecule > 0) THEN
       ALLOCATE (id_name(nmoleculekind),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "id_name",nmoleculekind*int_size,&
                                        para_env)
       ! Read molecule kind names
       DO ikind=1,nmoleculekind
          IF (para_env%ionode) THEN
             READ (UNIT=input_unit,IOSTAT=istat) string
             IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                            "string (IOSTAT = "//&
                                            TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                            input_unit,para_env)
          END IF
          CALL mp_bcast(string,para_env%source,para_env%group)
          id_name(ikind) = str2id(string)
       END DO
       ! Read molecule kind numbers
       IF (para_env%ionode) THEN
          READ (UNIT=input_unit,IOSTAT=istat) ibuf(1:natom)
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "ibuf(1:natom) (IOSTAT = "//&
                                         TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                         input_unit,para_env)
       END IF
       CALL mp_bcast(ibuf,para_env%source,para_env%group)
       DO iatom=1,natom
          ikind = ibuf(iatom)
          atom_info%id_molname(iatom) = id_name(ikind)
       END DO
       DEALLOCATE (id_name,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"id_name",&
                                        para_env)
       ! Read molecule index which is used also as residue id
       IF (para_env%ionode) THEN
          READ (UNIT=input_unit,IOSTAT=istat) atom_info%resid(1:natom)
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "atom_info%resid(1:natom) (IOSTAT = "//&
                                         TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                         input_unit,para_env)
       END IF
       CALL mp_bcast(atom_info%resid,para_env%source,para_env%group)
       DO iatom=1,natom
          atom_info%id_resname(iatom) = str2id(s2s(cp_to_string(atom_info%resid(iatom))))
       END DO
    END IF
    DEALLOCATE (ibuf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ibuf",para_env)

    !MK to be checked ...
    topology%aa_element = .TRUE.
    topology%molname_generated = .FALSE.
    topology%natoms = natom

    IF (iw > 0) THEN
       WRITE (UNIT=iw,FMT="(T2,A)")&
         "BEGIN of COORD section data [Angstrom] read in binary format from file "//&
         TRIM(binary_restart_file_name)
       DO iatom=1,natom
         WRITE (UNIT=iw,FMT="(T2,A2,3(1X,ES25.16),2(1X,A))")&
           TRIM(ADJUSTL(id2str(atom_info%id_atmname(iatom)))),&
           atom_info%r(1:3,iatom)*angstrom,&
           TRIM(ADJUSTL(id2str(atom_info%id_molname(iatom)))),&
           TRIM(ADJUSTL(id2str(atom_info%id_resname(iatom))))
       END DO
       WRITE (UNIT=iw,FMT="(T2,A)")&
         "END of COORD section data [Angstrom] read from binary restart file "//&
         TRIM(binary_restart_file_name)
    END IF

    IF (para_env%ionode) CALL close_file(unit_number=input_unit,&
                                         keep_preconnection=.TRUE.)

    binary_file_read = .TRUE.

    CALL timestop(handle)

  END SUBROUTINE read_binary_coordinates

! *****************************************************************************
!> \brief   Read the input section &CORE_COORD or &SHELL_COORD from an external
!>          file written in binary format.
!> \param prefix ...
!> \param particle_set ...
!> \param root_section ...
!> \param subsys_section ...
!> \param binary_file_read ...
!> \param error ...
!> \par History
!>      - Creation (17.02.2011,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_binary_cs_coordinates(prefix,particle_set,root_section,&
                                        subsys_section,binary_file_read,error)

    CHARACTER(LEN=*), INTENT(IN)             :: prefix
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: root_section, subsys_section
    LOGICAL, INTENT(OUT)                     :: binary_file_read
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_binary_cs_coordinates', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: binary_restart_file_name, &
                                                message
    CHARACTER(LEN=default_string_length)     :: section_label, section_name
    INTEGER                                  :: handle, input_unit, &
                                                iparticle, istat, iw, nbuf, &
                                                nparticle
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ibuf
    LOGICAL                                  :: exit_routine, failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL timeset(routineN,handle)

    NULLIFY (logger)
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(root_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subsys_section),cp_failure_level,routineP,error,failure)
    logger => cp_error_get_logger(error)
    para_env => logger%para_env

    binary_file_read = .FALSE.

    IF (ASSOCIATED(particle_set)) THEN
       exit_routine = .FALSE.
       nparticle = SIZE(particle_set)
    ELSE
       exit_routine = .TRUE.
       nparticle = 0
    END IF

    CALL section_vals_val_get(root_section,"EXT_RESTART%BINARY_RESTART_FILE_NAME",&
                              c_val=binary_restart_file_name,error=error)

    IF (TRIM(ADJUSTL(binary_restart_file_name)) == "") THEN
       CALL timestop(handle)
       RETURN
    END IF

    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/XYZ_INFO",&
                              extension=".subsysLog",error=error)

    section_name = prefix//" COORDINATES"

    ! Open binary restart file at last position
    IF (para_env%ionode) THEN
       CALL open_file(file_name=TRIM(binary_restart_file_name),&
                      file_status="OLD",&
                      file_form="UNFORMATTED",&
                      file_action="READWRITE",&
                      file_position="ASIS",&
                      unit_number=input_unit,&
                      debug=iw)
       READ (UNIT=input_unit,IOSTAT=istat) section_label,nbuf
       IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                      "section_label, nbuf -> "//TRIM(section_label)//", "//&
                                      TRIM(ADJUSTL(cp_to_string(nbuf)))//&
                                      " (IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//"). "//&
                                      "Section "//TRIM(ADJUSTL(section_name))//" was expected.",&
                                      input_unit,para_env)
       IF (TRIM(section_label) == TRIM(section_name)) THEN
          IF (nbuf /= nparticle) THEN
             IF (iw > 0) THEN
                message = "INFO: The requested number of "//TRIM(section_name)//" ("//&
                          TRIM(ADJUSTL(cp_to_string(nparticle)))//") does not match the "//&
                          "number ("//TRIM(ADJUSTL(cp_to_string(nbuf)))//") available from the "//&
                          "binary restart file <"//TRIM(binary_restart_file_name)//&
                          ">. The restart file information is ignored."
                CALL print_message(message,iw,1,1,1)
             END IF
             ! Ignore this section
             IF (nbuf > 0) THEN
                ! Perform dummy read
                ALLOCATE (rbuf(3,nbuf),STAT=istat)
                IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                 "rbuf",dp_size*3*nbuf,&
                                                 para_env)
                READ (UNIT=input_unit,IOSTAT=istat) rbuf(1:3,1:nbuf)
                IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                               "rbuf(1:3,1:nbuf) -> "//prefix//&
                                               " coordinates (IOSTAT = "//&
                                               TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                               input_unit,para_env)
                DEALLOCATE (rbuf,STAT=istat)
                IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)
                ALLOCATE (ibuf(nbuf),STAT=istat)
                IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                 "ibuf",int_size*nbuf,&
                                                 para_env)
                READ (UNIT=input_unit,IOSTAT=istat) ibuf(1:nbuf)
                IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                               "ibuf(1:nparticle) -> atomic indices of the "//&
                                               TRIM(section_name)//" (IOSTAT = "//&
                                               TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                               input_unit,para_env)
                DEALLOCATE (ibuf,STAT=istat)
                IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ibuf",para_env)
             END IF
             exit_routine = .TRUE.
          ELSE
             IF (iw > 0) THEN
                WRITE (UNIT=iw,FMT="(T2,A,T71,I10)")&
                 "Number of "//prefix//" particles:",nparticle
             END IF
             IF (nparticle == 0) exit_routine = .TRUE.
          END IF
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Section label <"//TRIM(section_label)//"> read from the "//&
                            "binary restart file <"//TRIM(binary_restart_file_name)//&
                            "> does not match the requested section name <"//&
                            TRIM(section_name)//">.")
       END IF
    END IF

    CALL mp_bcast(exit_routine,para_env%source,para_env%group)
    IF (exit_routine) THEN
       IF (para_env%ionode) CALL close_file(unit_number=input_unit,&
                                            keep_preconnection=.TRUE.)
       CALL timestop(handle)
       RETURN
    END IF

    CPPrecondition((nparticle > 0),cp_failure_level,routineP,error,failure)

    ALLOCATE (rbuf(3,nparticle),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "rbuf",dp_size*3*nparticle,&
                                     para_env)

    IF (para_env%ionode) THEN
       READ (UNIT=input_unit,IOSTAT=istat) rbuf(1:3,1:nparticle)
       IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                      "rbuf(1:3,1:nparticle) -> "//prefix//&
                                      " coordinates (IOSTAT = "//&
                                      TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                      input_unit,para_env)
    END IF
    CALL mp_bcast(rbuf,para_env%source,para_env%group)

    DO iparticle=1,nparticle
       particle_set(iparticle)%r(1:3) = rbuf(1:3,iparticle)
    END DO

    DEALLOCATE (rbuf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)

    ALLOCATE (ibuf(nparticle),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "ibuf",int_size*nparticle,&
                                     para_env)

    IF (para_env%ionode) THEN
       READ (UNIT=input_unit,IOSTAT=istat) ibuf(1:nparticle)
       IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                      "ibuf(1:nparticle) -> atomic indices of the "//&
                                      TRIM(section_name)//" (IOSTAT = "//&
                                      TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                      input_unit,para_env)
    END IF

    CALL mp_bcast(ibuf,para_env%source,para_env%group)

    DO iparticle=1,nparticle
       particle_set(iparticle)%atom_index = ibuf(iparticle)
    END DO

    DEALLOCATE (ibuf,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"ibuf",para_env)

    IF (iw > 0) THEN
       WRITE (UNIT=iw,FMT="(T2,A)")&
        "BEGIN of "//TRIM(ADJUSTL(section_name))//&
        " section data [Angstrom] read in binary format from file "//&
        TRIM(binary_restart_file_name)
       DO iparticle=1,nparticle
         WRITE (UNIT=iw,FMT="(T2,A2,3(1X,ES25.16),1X,I0)")&
           TRIM(ADJUSTL(particle_set(iparticle)%atomic_kind%name)),&
           particle_set(iparticle)%r(1:3)*angstrom,&
           particle_set(iparticle)%atom_index
       END DO
       WRITE (UNIT=iw,FMT="(T2,A)")&
        "END of "//TRIM(ADJUSTL(section_name))//&
        " section data [Angstrom] read from binary restart file "//&
        TRIM(binary_restart_file_name)
    END IF

    IF (para_env%ionode) CALL close_file(unit_number=input_unit,&
                                         keep_preconnection=.TRUE.)

    binary_file_read = .TRUE.

    CALL timestop(handle)

  END SUBROUTINE read_binary_cs_coordinates

! *****************************************************************************
!> \brief   Read the input section &VELOCITY, &CORE_VELOCITY, or
!>          &SHELL_VELOCITY from an external file written in binary format.
!> \param prefix ...
!> \param particle_set ...
!> \param root_section ...
!> \param para_env ...
!> \param subsys_section ...
!> \param binary_file_read ...
!> \param error ...
!> \par History
!>      - Creation (17.02.2011,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_binary_velocities(prefix,particle_set,root_section,para_env,&
                                    subsys_section,binary_file_read,error)

    CHARACTER(LEN=*), INTENT(IN)             :: prefix
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: subsys_section
    LOGICAL, INTENT(OUT)                     :: binary_file_read
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_binary_velocities', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: binary_restart_file_name, &
                                                message
    CHARACTER(LEN=default_string_length)     :: section_label, section_name
    INTEGER                                  :: handle, i, input_unit, &
                                                iparticle, istat, iw, nbuf, &
                                                nparticle
    LOGICAL                                  :: failure, have_velocities
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rbuf
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    NULLIFY (logger)
    failure = .FALSE.
    CPPrecondition(ASSOCIATED(root_section),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(subsys_section),cp_failure_level,routineP,error,failure)
    logger => cp_error_get_logger(error)

    binary_file_read = .FALSE.

    CALL section_vals_val_get(root_section,"EXT_RESTART%BINARY_RESTART_FILE_NAME",&
                              c_val=binary_restart_file_name,error=error)

    IF (TRIM(ADJUSTL(binary_restart_file_name)) == "") THEN
       CALL timestop(handle)
       RETURN
    END IF

    iw = cp_print_key_unit_nr(logger,subsys_section,"PRINT%TOPOLOGY_INFO/XYZ_INFO",&
                              extension=".subsysLog",error=error)

    IF (LEN_TRIM(prefix) == 0) THEN
       section_name = "VELOCITIES"
    ELSE
       section_name = prefix//" VELOCITIES"
    END IF

    have_velocities = .FALSE.

    IF (ASSOCIATED(particle_set)) THEN
       nparticle = SIZE(particle_set)
    ELSE
       nparticle = 0
    END IF

    ! Open binary restart file at last position and check if there are
    ! velocities available
    IF (para_env%ionode) THEN
       CALL open_file(file_name=binary_restart_file_name,&
                      file_status="OLD",&
                      file_form="UNFORMATTED",&
                      file_action="READWRITE",&
                      file_position="ASIS",&
                      unit_number=input_unit,&
                      debug=iw)
       DO
          READ (UNIT=input_unit,IOSTAT=istat) section_label,nbuf
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "section_label, nbuf -> "//TRIM(section_label)//", "//&
                                         TRIM(ADJUSTL(cp_to_string(nbuf)))//&
                                         " (IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//"). "//&
                                         "Section "//TRIM(ADJUSTL(section_name))//" was expected.",&
                                         input_unit,para_env)
          IF (INDEX(section_label,"THERMOSTAT") > 0) THEN
             IF (nbuf > 0) THEN
                ! Ignore thermostat information
                ALLOCATE (rbuf(nbuf,1),STAT=istat)
                IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                                 "rbuf",dp_size*nbuf,para_env)
                ! Perform dummy read
                DO i=1,4
                   READ (UNIT=input_unit,IOSTAT=istat) rbuf(1:nbuf,1)
                   IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                                  "rbuf(1:nbuf,1) -> "//&
                                                  TRIM(ADJUSTL(section_label))//&
                                                  " (IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                                  input_unit,para_env)
                END DO
                DEALLOCATE (rbuf,STAT=istat)
                IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)
                IF (iw > 0) THEN
                   message = "INFO: Ignoring section <"//TRIM(ADJUSTL(section_label))//&
                             "> from binary restart file <"//TRIM(binary_restart_file_name)//">."
                   CALL print_message(message,iw,1,1,1)
                END IF
             END IF
             CYCLE
          ELSE IF (INDEX(section_label,"VELOCIT") == 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                               "Section label <"//TRIM(section_label)//"> read from the "//&
                               "binary restart file <"//TRIM(binary_restart_file_name)//&
                               "> does not match the requested section name <"//&
                               TRIM(section_name)//">.")
          ELSE
             IF (nbuf > 0) have_velocities = .TRUE.
             EXIT
          END IF
       END DO
    END IF

    CALL mp_bcast(nbuf,para_env%source,para_env%group)
    CALL mp_bcast(have_velocities,para_env%source,para_env%group)

    IF (have_velocities) THEN

       ALLOCATE (rbuf(3,nbuf),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "rbuf",dp_size*3*nbuf,&
                                        para_env)

       IF (para_env%ionode) THEN
          READ (UNIT=input_unit,IOSTAT=istat) rbuf(1:3,1:nbuf)
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "rbuf(1:3,1:nbuf) -> "//&
                                         TRIM(ADJUSTL(section_name))//&
                                         " (IOSTAT = "//TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                         input_unit,para_env)
       END IF

       IF (nbuf == nparticle) THEN
          CALL mp_bcast(rbuf,para_env%source,para_env%group)
          DO iparticle=1,nparticle
             particle_set(iparticle)%v(1:3) = rbuf(1:3,iparticle)
          END DO
       ELSE
          IF (iw > 0) THEN
             message = "INFO: The requested number of "//TRIM(ADJUSTL(section_name))//&
                       " ("//TRIM(ADJUSTL(cp_to_string(nparticle)))//") does not match the "//&
                       "number ("//TRIM(ADJUSTL(cp_to_string(nbuf)))//") available from the "//&
                       "binary restart file <"//TRIM(binary_restart_file_name)//&
                       ">. The restart file information is ignored."
             CALL print_message(message,iw,1,1,1)
          END IF
       END IF

       DEALLOCATE (rbuf,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)

    END IF

    IF (nbuf == nparticle) THEN
       IF (iw > 0) THEN
          WRITE (UNIT=iw,FMT="(T2,A)")&
           "BEGIN of "//TRIM(ADJUSTL(section_name))//&
           " section data [a.u.] read in binary format from file "//&
           TRIM(binary_restart_file_name)
          IF (have_velocities) THEN
             DO iparticle=1,nparticle
                WRITE (UNIT=iw,FMT="(T2,A2,3(1X,ES25.16))")&
                 TRIM(ADJUSTL(particle_set(iparticle)%atomic_kind%name)),&
                 particle_set(iparticle)%v(1:3)
             END DO
          ELSE
             WRITE (UNIT=iw,FMT="(A)")&
              "# No "//TRIM(ADJUSTL(section_name))//" available"
          END IF
          WRITE (UNIT=iw,FMT="(T2,A)")&
           "END of "//TRIM(ADJUSTL(section_name))//&
           " section data [a.u.] read from binary restart file "//&
           TRIM(binary_restart_file_name)
       END IF
       binary_file_read = .TRUE.
    END IF

    IF (para_env%ionode) CALL close_file(unit_number=input_unit,&
                                         keep_preconnection=.TRUE.)

    CALL timestop(handle)

  END SUBROUTINE read_binary_velocities

! *****************************************************************************
!> \brief   Read the input section &THERMOSTAT for Nose thermostats from an
!>          external file written in binary format.
!> \param prefix ...
!> \param nhc ...
!> \param binary_restart_file_name ...
!> \param restart ...
!> \param para_env ...
!> \param error ...
!> \par History
!>      - Creation (28.02.2011,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE read_binary_thermostats_nose(prefix,nhc,binary_restart_file_name,&
                                          restart,para_env,error)

    CHARACTER(LEN=*), INTENT(IN)             :: prefix
    TYPE(lnhc_parameters_type), POINTER      :: nhc
    CHARACTER(LEN=*), INTENT(IN)             :: binary_restart_file_name
    LOGICAL, INTENT(OUT)                     :: restart
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_binary_thermostats_nose', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: section_label, section_name
    INTEGER                                  :: handle, i, idx, input_unit, &
                                                istat, j, nhc_size
    LOGICAL                                  :: debug, failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rbuf
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(nhc),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(para_env),cp_failure_level,routineP,error,failure)

    ! Set to .TRUE. for debug mode, i.e. all data read are written to stdout
    NULLIFY (logger)
    logger => cp_error_get_logger(error)
    IF (logger%iter_info%print_level >= debug_print_level) THEN
       debug = .TRUE.
    ELSE
       debug = .FALSE.
    END IF

    restart = .FALSE.

    section_name = prefix//" THERMOSTATS"

    ! Open binary restart file at last position
    IF (para_env%ionode) THEN
       CALL open_file(file_name=binary_restart_file_name,&
                      file_status="OLD",&
                      file_form="UNFORMATTED",&
                      file_action="READWRITE",&
                      file_position="ASIS",&
                      unit_number=input_unit)
       READ (UNIT=input_unit,IOSTAT=istat) section_label,nhc_size
       IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                      "nhc_size (IOSTAT = "//&
                                      TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                      input_unit,para_env)
       IF (INDEX(section_label,"THERMOSTAT") == 0) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Section label <"//TRIM(section_label)//"> read from the "//&
                            "binary restart file <"//TRIM(binary_restart_file_name)//&
                            "> does not match the requested section name <"//&
                            TRIM(section_name)//">.")
       END IF
       IF (debug) THEN
          WRITE (UNIT=*,FMT="(T2,A,/,T2,A,I0)")&
            "BEGIN of "//TRIM(ADJUSTL(section_label))//&
            " section data read in binary format from file "//&
            TRIM(binary_restart_file_name),&
            "# nhc_size = ",nhc_size
       END IF
    END IF

    CALL mp_bcast(nhc_size,para_env%source,para_env%group)

    IF (nhc_size > 0) THEN

       ALLOCATE (rbuf(nhc_size),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "rbuf",dp_size*nhc_size,para_env)
       rbuf(:) = 0.0_dp

       ! Read NHC section &COORD
       IF (para_env%ionode) THEN
          READ (UNIT=input_unit,IOSTAT=istat) rbuf(1:nhc_size)
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "eta -> rbuf (IOSTAT = "//&
                                         TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                         input_unit,para_env)
          IF (debug) THEN
             WRITE (UNIT=*,FMT="(T2,A,/,(4(1X,ES25.16)))")&
               "&COORD",rbuf(1:nhc_size)
          END IF
       END IF
       CALL mp_bcast(rbuf,para_env%source,para_env%group)
       DO i=1,SIZE(nhc%nvt,2)
          idx = (nhc%map_info%index(i) - 1)*nhc%nhc_len
          DO j=1,SIZE(nhc%nvt,1)
             idx = idx + 1
             nhc%nvt(j,i)%eta = rbuf(idx)
          END DO
       END DO

       ! Read NHC section &VELOCITY
       IF (para_env%ionode) THEN
          READ (UNIT=input_unit,IOSTAT=istat) rbuf(1:nhc_size)
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "veta -> rbuf (IOSTAT = "//&
                                         TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                         input_unit,para_env)
          IF (debug) THEN
             WRITE (UNIT=*,FMT="(T2,A,/,(4(1X,ES25.16)))")&
               "&VELOCITY",rbuf(1:nhc_size)
          END IF
       END IF
       CALL mp_bcast(rbuf,para_env%source,para_env%group)
       DO i=1,SIZE(nhc%nvt,2)
          idx = (nhc%map_info%index(i) - 1)*nhc%nhc_len
          DO j=1,SIZE(nhc%nvt,1)
             idx = idx + 1
             nhc%nvt(j,i)%v = rbuf(idx)
          END DO
       END DO

       ! Read NHC section &MASS
       IF (para_env%ionode) THEN
          READ (UNIT=input_unit,IOSTAT=istat) rbuf(1:nhc_size)
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "mnhc -> rbuf (IOSTAT = "//&
                                         TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                         input_unit,para_env)
          IF (debug) THEN
             WRITE (UNIT=*,FMT="(T2,A,/,(4(1X,ES25.16)))")&
               "&MASS:",rbuf(1:nhc_size)
          END IF
       END IF
       CALL mp_bcast(rbuf,para_env%source,para_env%group)
       DO i=1,SIZE(nhc%nvt,2)
          idx = (nhc%map_info%index(i) - 1)*nhc%nhc_len
          DO j=1,SIZE(nhc%nvt,1)
             idx = idx + 1
             nhc%nvt(j,i)%mass = rbuf(idx)
          END DO
       END DO

       ! Read NHC section &FORCE
       IF (para_env%ionode) THEN
          READ (UNIT=input_unit,IOSTAT=istat) rbuf(1:nhc_size)
          IF (istat /= 0) CALL stop_read(routineN,moduleN,__LINE__,&
                                         "fnhc -> rbuf (IOSTAT = "//&
                                         TRIM(ADJUSTL(cp_to_string(istat)))//")",&
                                         input_unit,para_env)
          IF (debug) THEN
             WRITE (UNIT=*,FMT="(T2,A,/,(4(1X,ES25.16)))")&
               "&FORCE",rbuf(1:nhc_size)
          END IF
       END IF
       CALL mp_bcast(rbuf,para_env%source,para_env%group)
       DO i=1,SIZE(nhc%nvt,2)
          idx = (nhc%map_info%index(i) - 1)*nhc%nhc_len
          DO j=1,SIZE(nhc%nvt,1)
             idx = idx + 1
             nhc%nvt(j,i)%f = rbuf(idx)
          END DO
       END DO

       DEALLOCATE (rbuf,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"rbuf",para_env)

       restart = .TRUE.

    END IF

    IF (para_env%ionode) THEN
       IF (debug) THEN
          WRITE (UNIT=*,FMT="(T2,A)")&
            "END of"//TRIM(ADJUSTL(section_label))//&
            " section data read in binary format from file "//&
            TRIM(binary_restart_file_name)
       END IF
       CALL close_file(unit_number=input_unit,&
                       keep_preconnection=.TRUE.)
    END IF

    CALL timestop(handle)

  END SUBROUTINE read_binary_thermostats_nose

END MODULE input_cp2k_binary_restarts
