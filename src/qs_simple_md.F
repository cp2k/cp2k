!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_simple_md [1.0] *
!!
!!   NAME
!!     qs_simple_md
!!
!!   FUNCTION
!!     Perform a simple velocity Verlet molecular dynamics (MD) run using
!!     QUICKSTEP (NVE ensemble)
!!
!!   AUTHOR
!!     Matthias Krack (07.09.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_simple_md

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp,&
                                             wp_size => dp_size
  USE md_qs_energies,                  ONLY: md_energy_qs,&
                                             set_qs_energies
  USE message_passing,                 ONLY: mp_bcast
  USE particle_types,                  ONLY: particle_type,&
                                             write_qs_particle_coordinates
  USE physcon,                         ONLY: kelvin,&
                                             femtoseconds
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_forces
  USE qs_parser,                       ONLY: finish_parser,&
                                             open_file,&
                                             read_object,&
                                             search_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE qs_wf_history_methods,           ONLY: wfi_create
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_linear_p_method_nr,&
                                             wfi_release
  USE simulation_cell,                 ONLY: cell_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: gasdev

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_simple_md"

! *** Data types ***

  TYPE simple_md_env_type
    REAL(wp) :: temperature,&
                time_step
    INTEGER  :: first_step,&
                max_step
  END TYPE simple_md_env_type

  PUBLIC :: simple_md

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE simple_md(qs_env,globenv)

!   Purpose: Driver routine for a simple velocity Verlet molecular dynamics
!            run using QUICKSTEP.

!   History: - Creation (05.09.2003,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), POINTER        :: qs_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "simple_md"

!   *** Local variables ***

    TYPE(atomic_kind_type), POINTER   :: atomic_kind
    TYPE(cell_type), POINTER          :: cell
    TYPE(qs_energy_type), POINTER     :: energy
    TYPE(qs_wf_history_type), POINTER :: wf_history

    TYPE(simple_md_env_type) :: simple_md_env
    REAL(wp)                 :: dm,dt,dt2,e_com,e_kinetic,e_thermal,e0,mass,&
                                temperature,scale,total_mass,v2
    INTEGER                  :: atom,dof,handle,i,iatom,idum,ikind,istat,&
                                istep,natom,nkind,output_unit
    LOGICAL                  :: ionode

    REAL(wp), DIMENSION(3) :: v,v_com,v_com_new

    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    INTEGER, DIMENSION(:), POINTER                :: atom_list

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    idum = globenv%idum
    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL read_simple_md_section(simple_md_env,globenv)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    particle_set=particle_set)

    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

!   *** Initialize wavefunction history ***

    CALL wfi_create(wf_history,interpolation_method_nr=wfi_linear_p_method_nr)
    CALL set_qs_env (qs_env=qs_env,wf_history=wf_history)

    dof = 3*natom
    dt = simple_md_env%time_step
    dt2 = 0.5_wp*dt*dt
    e0 = 0.0_wp

!   *** Initialize the velocities ***

    e_kinetic = 0.0_wp
    total_mass = 0.0_wp

    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list,&
                           mass=mass)
      v2 = 0.0_wp
      natom = SIZE(atom_list)
      DO iatom=1,natom
        atom = atom_list(iatom)
        particle_set(atom)%v(1) = gasdev(idum)
        particle_set(atom)%v(2) = gasdev(idum)
        particle_set(atom)%v(3) = gasdev(idum)
        v2 = v2 + particle_set(atom)%v(1)*particle_set(atom)%v(1) +&
                  particle_set(atom)%v(2)*particle_set(atom)%v(2) +&
                  particle_set(atom)%v(3)*particle_set(atom)%v(3)
      END DO
      e_kinetic = e_kinetic + 0.5_wp*mass*v2
      total_mass = total_mass + REAL(natom,wp)*mass
    END DO

    temperature = 2.0_wp*e_kinetic/REAL(dof,wp)

    scale = SQRT(simple_md_env%temperature/temperature)

!   *** Scale the velocities and calculate center of mass motion (COM) ***

    v_com(:) = 0.0_wp

    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list,&
                           mass=mass)
      natom = SIZE(atom_list)
      DO iatom=1,natom
        atom = atom_list(iatom)
        particle_set(atom)%v(:) = scale*particle_set(atom)%v(:)
        v_com(:) = v_com(:) + mass*particle_set(atom)%v(:)
      END DO
    END DO

    v_com(:) = v_com(:)/total_mass

!   *** Remove the center of mass motion and calculate the initial temperature ***

    e_kinetic = 0.0_wp
    v_com_new(:) = 0.0_wp

    DO ikind=1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom_list,&
                           mass=mass)
      v2 = 0.0_wp
      natom = SIZE(atom_list)
      DO iatom=1,natom
        atom = atom_list(iatom)
        particle_set(atom)%v(:) = particle_set(atom)%v(:) - v_com(:)
        v2 = v2 + particle_set(atom)%v(1)*particle_set(atom)%v(1) +&
                  particle_set(atom)%v(2)*particle_set(atom)%v(2) +&
                  particle_set(atom)%v(3)*particle_set(atom)%v(3)
        v_com_new(:) = v_com_new(:) + mass*particle_set(atom)%v(:)
      END DO
      e_kinetic = e_kinetic + 0.5_wp*mass*v2
    END DO

    v_com_new(:) = v_com_new(:)/total_mass
    temperature = 2.0_wp*e_kinetic/REAL(dof,wp)

    IF (ionode) THEN
      WRITE (UNIT=output_unit,FMT="(/,T2,A,F10.2,A,/,/,(T2,A,F20.12))")&
        "Initial temperature = ",temperature*kelvin," K",&
        "Center of mass motion: x =",v_com_new(1),&
        "                       y =",v_com_new(2),&
        "                       z =",v_com_new(3)
    END IF

!   *** Main MD loop ***

    DO istep=simple_md_env%first_step,simple_md_env%max_step

      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,T2,A,I6)")&
          "MD timestep number: ",istep
      END IF

!     *** Velocity Verlet (part 1) ***

      IF (istep > simple_md_env%first_step) THEN
        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               mass=mass)
          dm = 0.5_wp*dt/mass
          natom = SIZE(atom_list)
          DO iatom=1,natom
            atom = atom_list(iatom)
            particle_set(atom)%v(:) = particle_set(atom)%v(:) +&
                                      dm*particle_set(atom)%f(:)
            particle_set(atom)%r(:) = particle_set(atom)%r(:) +&
                                      dt*particle_set(atom)%v(:)
          END DO
        END DO
        IF (ionode) THEN
          WRITE (UNIT=output_unit,FMT="(/,T2,A,F10.3,A,F10.3,A)")&
            "MD time step: ",dt," a.u. = ",dt*femtoseconds," fs"
        END IF
!MK     CALL write_qs_particle_coordinates(particle_set,cell,globenv)
      END IF

!     *** Update energies and forces ***

      CALL qs_forces(qs_env,globenv)

!     *** Velocity Verlet (part 2) ***

      IF (istep > simple_md_env%first_step) THEN
        DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,&
                               atom_list=atom_list,&
                               mass=mass)
          dm = 0.5_wp*dt/mass
          natom = SIZE(atom_list)
          DO iatom=1,natom
            atom = atom_list(iatom)
            particle_set(atom)%v(:) = particle_set(atom)%v(:) +&
                                      dm*particle_set(atom)%f(:)
          END DO
        END DO
      END IF

!     *** Calculate center of mass motion (COM) of the current MD step ***

      e_kinetic = 0.0_wp
      v_com(:) = 0.0_wp

      DO ikind=1,nkind
        atomic_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             atom_list=atom_list,&
                             mass=mass)
        v2 = 0.0_wp
        natom = SIZE(atom_list)
        DO iatom=1,natom
          atom = atom_list(iatom)
          v_com(:) = v_com(:) + mass*particle_set(atom)%v(:)
          v2 = v2 + particle_set(atom)%v(1)*particle_set(atom)%v(1) +&
                    particle_set(atom)%v(2)*particle_set(atom)%v(2) +&
                    particle_set(atom)%v(3)*particle_set(atom)%v(3)
        END DO
        e_kinetic = e_kinetic + 0.5_wp*mass*v2
      END DO

      e_com = 0.5_wp*(v_com(1)*v_com(1) +&
                      v_com(2)*v_com(2) +&
                      v_com(3)*v_com(3))/total_mass
      e_thermal = e_kinetic - e_com

      temperature = 2.0_wp*e_thermal/REAL(dof,wp)

      CALL get_qs_env(qs_env=qs_env,energy=energy)

      IF (istep == simple_md_env%first_step) e0 = e_kinetic + energy%total

      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,(T2,A,F20.12))")&
          "Center of mass motion: x =",v_com(1),&
          "                       y =",v_com(2),&
          "                       z =",v_com(3)
        WRITE (UNIT=output_unit,FMT="(/,T2,A,F10.2,A,/,/,(T2,A,F20.12,A))")&
          "Temperature     = ",temperature*kelvin," K",&
          "E(thermal)      = ",e_thermal," a.u.",&
          "E(kinetic)      = ",e_kinetic," a.u.",&
          "E(potential)    = ",energy%total," a.u.",&
          "E(total)        = ",e_kinetic + energy%total," a.u.",&
          "E(total) - E(0) = ",e_kinetic + energy%total - e0," a.u."
      END IF

    END DO

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE simple_md

! *****************************************************************************

  SUBROUTINE read_simple_md_section(simple_md_env,globenv)

!   Purpose: Read the Driver routine for a simple velocity Verlet molecular dynamics
!            run using QUICKSTEP.

!   History: - Creation (05.09.2003,MK)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(simple_md_env_type), INTENT(OUT)     :: simple_md_env

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_simple_md"
    CHARACTER(LEN=*), PARAMETER :: start_section_label = "SIMPLE_MD"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword
    CHARACTER(LEN=3)  :: test_result
    LOGICAL           :: section_found

!   ---------------------------------------------------------------------------

    simple_md_env%first_step = 0
    simple_md_env%max_step = 100
    simple_md_env%temperature = 300.0_wp/kelvin ! 300 K -> a.u.
    simple_md_env%time_step = 1.0_wp/femtoseconds ! 1 fs -> a.u.

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label=start_section_label,&
                      section_found=section_found)

    IF (section_found) THEN
      DO
        test_result = test_object()
        IF (test_result /= "EOL") THEN
          CALL read_object(keyword)
          CALL stop_parser(routine_name,module_name,__LINE__,"INVALID_KEYWORD")
        END IF
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine_name,module_name,__LINE__,"EOF")
        ELSE
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("FIRST_STEP")
            CALL read_object(simple_md_env%first_step,lower_to_upper=.TRUE.)
          CASE ("MAX_STEP")
            CALL read_object(simple_md_env%max_step,lower_to_upper=.TRUE.)
          CASE ("TEMPERATURE","T")
            CALL read_object(simple_md_env%temperature,lower_to_upper=.TRUE.)
            simple_md_env%temperature = simple_md_env%temperature/kelvin
          CASE ("TIME_STEP")
            CALL read_object(simple_md_env%time_step,lower_to_upper=.TRUE.)
            simple_md_env%time_step = simple_md_env%time_step/femtoseconds
          CASE DEFAULT
            CALL stop_parser(routine_name,module_name,__LINE__,"INVALID_KEYWORD")
          END SELECT
        END IF
      END DO
    END IF

    CALL finish_parser()

  END SUBROUTINE read_simple_md_section

! *****************************************************************************

END MODULE qs_simple_md
