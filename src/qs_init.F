!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_init [1.0] *
!!
!!   NAME
!!     qs_init
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MK (18.05.2000)
!!
!!   MODIFICATION HISTORY
!!     - Merged with the Quickstep MODULE method_specification (17.01.2002, MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_init

  USE kinds, ONLY: int_size,&
                   wp => dp

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE dft_types,                       ONLY: dft_control_type,&
                                             read_dft_control,&
                                             read_qs_control,&
                                             write_dft_control,&
                                             write_qs_control
  USE gamma,                           ONLY: init_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE machine,                         ONLY: m_flush
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE particle_types,                  ONLY: particle_type,&
                                             init_particle_set,&
                                             read_particle_set,&
                                             write_particle_coordinates,&
                                             write_particle_distances
  USE qs_environment_types,            ONLY: init_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE simulation_cell,                 ONLY: cell_type,&
                                             read_cell,&
                                             write_cell
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE pao_types,                       ONLY: pao_env_type
  USE pao_read_section,                ONLY: pao_env_init_from_file
  use cp_para_types,                   ONLY: cp_para_env_type
  use cp_para_env,                     only: cp_para_env_retain, &
       cp_para_env_release, cp_para_env_from_globenv
  use scf_control_types,               only: scf_control_type,&
       scf_c_create, scf_c_read_parameters, scf_c_write_parameters,&
       scf_c_release

  IMPLICIT NONE

  PRIVATE

! *** Public subroutines ***

  PUBLIC :: read_qs_env

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_qs_env(qs_env,globenv)

!   Purpose: Read the input and the database files for the setup of the
!            QUICKSTEP environment.

!   History: - Creation (22.05.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)    :: qs_env
    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_qs_env (MODULE qs_init)"

!   *** Local variables ***

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(dft_control_type), POINTER               :: dft_control
    TYPE(pao_env_type), POINTER                   :: pao_env
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    TYPE(qs_energy_type), POINTER                 :: energy
    TYPE(qs_force_type), DIMENSION(:), POINTER    :: force
    type(scf_control_type), pointer               :: scf_control
    TYPE(cp_para_env_type), pointer               :: para_env

    INTEGER :: handle,ikind,istat,lmax,maxl,maxlppl,maxlppnl,nkind,output_unit
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: natom_of_kind

!   ---------------------------------------------------------------------------

    CALL timeset("read_qs_env","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Initialise the QUICKSTEP environment ***

    NULLIFY (atomic_kind_set)
    NULLIFY (cell)
    NULLIFY (dft_control)
    NULLIFY (pao_env)
    NULLIFY (particle_set)
    NULLIFY (energy)
    NULLIFY (force)
    nullify (scf_control)

    para_env=>cp_para_env_from_globenv(globenv) 
    CALL init_qs_env(qs_env, para_env=para_env)

!   *** Print the Quickstep program banner (copyright and version number) ***

    CALL write_qs_program_banner(globenv)

!   *** Read the input section with the DFT control parameters ***

    CALL read_dft_control(dft_control,globenv)

!   *** Print the DFT control parameters ***

    CALL write_dft_control(dft_control,globenv)

!   *** Read the input section with the QUICKSTEP control parameters ***

    CALL read_qs_control(dft_control%qs_control,globenv)

!   *** Print the QUICKSTEP control parameters ***

    CALL write_qs_control(dft_control%qs_control,globenv)

!   *** Read the input section with the cell parameters ***

    CALL read_cell(cell,globenv)

!   *** Print the cell parameters ***

    CALL write_cell(cell,globenv)

!   *** Read the input section with the atomic coordinates ***

    CALL read_particle_set(particle_set,atomic_kind_set,cell,globenv)

!   *** Read the atomic kind information ***

    CALL read_atomic_kind_set(atomic_kind_set,globenv)

!   *** Print the unnormalized basis set information (input data) ***

    CALL write_gto_basis_sets(atomic_kind_set,globenv)

!   *** Initialize the spherical harmonics and ***
!   *** the orbital transformation matrices    ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxl=maxl,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl)

    lmax = MAX(maxl,maxlppl,maxlppnl)

    CALL init_spherical_harmonics(lmax,globenv)

!   *** Initialize the pretabulation for the calculation of the   ***
!   *** incomplete Gamma function F_n(t) after McMurchie-Davidson ***

    lmax = 3*maxl + 1

    CALL init_md_ftable(lmax)

!   *** Initialise the atomic kind set ***

    CALL init_atomic_kind_set(atomic_kind_set)

!   *** Print the atomic kind set ***

    CALL write_atomic_kind_set(atomic_kind_set,globenv)

!   *** Initialize the atomic data sets completely ***

    CALL init_particle_set(particle_set)

!   *** Print the total number of kinds, atoms, basis functions etc. ***

    CALL write_total_numbers(atomic_kind_set,particle_set,globenv)

!   *** Print the atomic coordinates ***

    CALL write_particle_coordinates(particle_set,cell,globenv)

!   *** Print the interatomic distances ***

    CALL write_particle_distances(particle_set,cell,globenv)

!   *** Initialize the atomic interaction radii ***

    CALL init_interaction_radii(dft_control%qs_control,cell,atomic_kind_set,&
                                globenv)

!   *** SCF parameters ***
    
    call scf_c_create(scf_control)
    CALL scf_c_read_parameters(scf_control,globenv)
    CALL scf_c_write_parameters(scf_control,globenv)

!   *** Allocate the data structure for Quickstep energies ***

    CALL allocate_qs_energy(energy)

!   *** Allocate the force data structure ***

    IF (dft_control%forces) THEN
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"natom_of_kind",nkind*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom_of_kind=natom_of_kind) 
      CALL allocate_qs_force(force,natom_of_kind)
      DEALLOCATE (natom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"natom_of_kind")
    END IF

!   *** Get the current QUICKSTEP environment ***

    CALL set_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    pao_env=pao_env,&
                    particle_set=particle_set,&
                    energy=energy,&
                    force=force,&
                    scf_control=scf_control)

    call scf_c_release(scf_control)

!   *** Read the pao section ***

    if (dft_control%qs_control%pao) then
       allocate(pao_env) ! check status ?
       call set_qs_env(qs_env=qs_env,pao_env=pao_env)
       call pao_env_init_from_file(qs_env,globenv)
       !call cp_write_description(pao_env, unit_nr=globenv%scr)
    end if

    IF (ionode) CALL m_flush(output_unit)

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE read_qs_env

! *****************************************************************************

  SUBROUTINE write_qs_program_banner(globenv)

!   Purpose: Write the Quickstep program banner with the copyright and the
!            program version information to the output unit.

!   History: - Creation (15.06.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(IN) :: globenv

!   *** Local variables ***

    INTEGER :: output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%program_banner) THEN

      output_unit = globenv%scr

      WRITE (UNIT=output_unit,FMT="(T2,A79)")&
        "*******************************************************************************",&
        "*******************************************************************************",&
        "**                                                                           **",&
        "**     #####                         ##              ##                      **",&
        "**    ##   ##            ##          ##              ##                      **",&
        "**   ##     ##                       ##            ######                    **",&
        "**   ##     ##  ##   ##  ##   #####  ##  ##   ####   ##    #####    #####    **",&
        "**   ##     ##  ##   ##  ##  ##      ## ##   ##      ##   ##   ##  ##   ##   **",&
        "**   ##  ## ##  ##   ##  ##  ##      ####     ###    ##   ######   ######    **",&
        "**    ##  ###   ##   ##  ##  ##      ## ##      ##   ##   ##       ##        **",&
        "**     #######   #####   ##   #####  ##  ##  ####    ##    #####   ##        **",&
        "**           ##                                                    ##        **",&
        "**                                                                           **",&
        "**                                                ... make the atoms dance   **",&
        "**                                                                           **",&
        "**                         Version 3.1 (June 2002)                           **",&
        "**                                                                           **",&
        "**            Copyright (C) by CP2K Developers Group (2000 - 2002)           **",&
        "**                                                                           **",&
        "*******************************************************************************",&
        "*******************************************************************************"

    END IF

  END SUBROUTINE write_qs_program_banner

! *****************************************************************************

  SUBROUTINE write_total_numbers(atomic_kind_set,particle_set,globenv)

!   Purpose: Write the total number of kinds, atoms, etc. to the logical unit
!            number lunit.

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(global_environment_type), INTENT(IN)     :: globenv
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set

!   *** Local variables ***

    INTEGER :: iatom,maxl,maxlppl,maxlppnl,natom,ncgf,nkind,npgf,nset,nsgf,&
               nshell,output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%total_numbers) THEN

      output_unit = globenv%scr

      natom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               maxl=maxl,&
                               maxlppl=maxlppl,&
                               maxlppnl=maxlppnl,&
                               ncgf=ncgf,&
                               npgf=npgf,&
                               nset=nset,&
                               nsgf=nsgf,&
                               nshell=nshell)

      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "TOTAL NUMBERS AND MAXIMUM NUMBERS"

      WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(T3,A,I6))")&
        "Total number of",&
        "- Atomic kinds:                  ",nkind,&
        "- Atoms:                         ",natom,&
        "- Shell sets:                    ",nset,&
        "- Shells:                        ",nshell,&
        "- Primitive Cartesian functions: ",npgf,&
        "- Cartesian basis functions:     ",ncgf,&
        "- Spherical basis functions:     ",nsgf

      IF ((maxlppl > -1).OR.(maxlppnl > -1)) THEN
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(T3,A,I6))")&
          "Maximum angular momentum quantum number of the",&
          "- Orbital basis functions:                   ",maxl,&
          "- Local part of the GTH pseudopotential:     ",maxlppl,&
          "- Non-local part of the GTH pseudopotential: ",maxlppnl
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(T3,A,I6))")&
          "Maximum angular momentum quantum number of the",&
          "- Orbital basis functions: ",maxl
      END IF

    END IF

  END SUBROUTINE write_total_numbers

! *****************************************************************************

END MODULE qs_init
