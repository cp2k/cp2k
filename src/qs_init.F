!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
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
!!     - Merged with the Quickstep MODULE method_specification (17.01.2002,MK)
!!     - USE statements cleaned, write_checkpoint_information added
!!       (25.09.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_init
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_blacs_env,                    ONLY: cp_blacs_env_from_globenv
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_f_matrix_struct,              ONLY: cp_f_matrix_struct_type,&
                                             cp_fmstruct_create,&
                                             cp_fmstruct_release,&
                                             cp_fmstruct_retain
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dft_types,                       ONLY: dft_control_type,&
                                             read_dft_control,&
                                             read_qs_control,&
                                             write_dft_control,&
                                             write_qs_control
  USE gamma,                           ONLY: init_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: int_size,&
                                             wp => dp
  USE machine,                         ONLY: m_flush
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE pao_read_section,                ONLY: pao_env_init_from_file
  USE pao_types,                       ONLY: pao_env_type
  USE particle_types,                  ONLY: init_particle_set,&
                                             particle_type,&
                                             read_particle_set,&
                                             write_particle_coordinates,&
                                             write_particle_distances
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_rebuild_aomo_pool,&
                                             qs_env_rebuild_pw_pools
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             init_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE scf_control_types,               ONLY: scf_c_create,&
                                             scf_c_read_parameters,&
                                             scf_c_release,&
                                             scf_c_write_parameters,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             read_cell,&
                                             write_cell
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_init"

! *** Public subroutines ***

  PUBLIC :: read_qs_env

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_qs_env(qs_env,globenv,error)

!   Purpose: Read the input and the database files for the setup of the
!            QUICKSTEP environment.

!   History: - Creation (22.05.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(qs_environment_type), INTENT(OUT)       :: qs_env
    TYPE(global_environment_type), INTENT(IN)    :: globenv
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_qs_env"
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

!   *** Local variables ***

    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(cell_type), POINTER                      :: cell
    TYPE(dft_control_type), POINTER               :: dft_control
    TYPE(pao_env_type), POINTER                   :: pao_env
    TYPE(particle_type), DIMENSION(:), POINTER    :: particle_set
    TYPE(qs_energy_type), POINTER                 :: energy
    TYPE(qs_force_type), DIMENSION(:), POINTER    :: force
    TYPE(scf_control_type), POINTER               :: scf_control
    TYPE(cp_para_env_type), POINTER               :: para_env
    TYPE(cp_blacs_env_type), POINTER              :: blacs_env
    TYPE(cp_f_matrix_struct_type), POINTER        :: fmstruct

    INTEGER :: handle,ikind,istat,lmax,maxl,maxlppl,maxlppnl,nkind,&
         output_unit, n_mo, n_ao
    LOGICAL :: ionode

    INTEGER, DIMENSION(:), ALLOCATABLE :: natom_of_kind

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Initialise the Quickstep environment ***

    NULLIFY (atomic_kind_set)
    NULLIFY (cell)
    NULLIFY (dft_control)
    NULLIFY (pao_env)
    NULLIFY (particle_set)
    NULLIFY (energy)
    NULLIFY (force)
    NULLIFY (scf_control,fmstruct, para_env, blacs_env)

    para_env => cp_para_env_from_globenv(globenv)
    blacs_env => cp_blacs_env_from_globenv(globenv)

    CALL init_qs_env(qs_env,para_env=para_env)

!   *** Print the Quickstep program banner (copyright and version number) ***

    CALL write_qs_program_banner(globenv)

!   *** Read the input section with the DFT control parameters ***

    CALL read_dft_control(dft_control,globenv)

!   *** Print the DFT control parameters ***

    CALL write_dft_control(dft_control,globenv)

!   *** Read the input section with the Quickstep control parameters ***

    CALL read_qs_control(dft_control%qs_control,globenv)

!   *** Print the Quickstep control parameters ***

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

    CALL init_atomic_kind_set(atomic_kind_set,globenv)

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

    CALL scf_c_create(scf_control)
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

!   *** Get the current Quickstep environment ***

    CALL set_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    pao_env=pao_env,&
                    particle_set=particle_set,&
                    energy=energy,&
                    force=force,&
                    scf_control=scf_control)

    CALL scf_c_release(scf_control)

!   *** Initialize pools ***
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, &
           dft_control=dft_control, scf_control=scf_control)
      CALL get_atomic_kind_set(atomic_kind_set, nsgf=n_ao, nelectron=n_mo)
      n_mo = n_mo-dft_control%charge

     IF (globenv%print%mo_eigenvalues.OR.&
         (scf_control%level_shift /= 0.0_wp).OR.&
         (scf_control%smear /= 0.0_wp)) THEN
       n_mo = n_ao
     ELSE
       n_mo = n_mo/2
     END IF

     CALL cp_fmstruct_create(fmstruct, nrow_global=n_ao,&
         ncol_global=n_mo, para_env=para_env,&
         context=blacs_env,&
         nrow_block=scf_control%nrow_block,&
         ncol_block=scf_control%ncol_block)
     CALL qs_env_rebuild_aomo_pool(qs_env,fmstruct,error=error)
     CALL cp_fmstruct_release(fmstruct,error=error)

     CALL qs_env_rebuild_pw_pools(qs_env,error=error)

!   *** Read the pao section ***

    IF (dft_control%qs_control%pao) THEN
      ALLOCATE (pao_env) ! check status ?
      CALL set_qs_env(qs_env=qs_env,pao_env=pao_env)
      CALL pao_env_init_from_file(qs_env,globenv)
!     CALL cp_write_description(pao_env, unit_nr=globenv%scr)
    END IF

    IF (ionode) CALL m_flush(output_unit)

    CALL timestop(0.0_wp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

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

      WRITE (UNIT=output_unit,FMT="(/,(T2,A79))")&
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
