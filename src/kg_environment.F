!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_environment [1.0] *
!!
!!   NAME
!!     kg_environment
!!
!!   FUNCTION
!!    initialize kg environment             
!!   AUTHOR
!!     gt SEPT-23-2002
!!   HISTORY 
!!     GT 11-13-2003 :            
!!     moved the main driver to cp2k 
!!     performs only the initialization of the environment
!!*****
!******************************************************************************

MODULE kg_environment
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: allocate_atomic_kind_set,&
                                             atomic_kind_p_type,&
                                             atomic_kind_type,&
                                             check_atomic_kind_set,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_type
  USE cp_control_utils,                ONLY: read_dft_control,&
                                             read_qs_control,&
                                             write_dft_control,&
                                             write_qs_control
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_retain
  USE cp_subsystem_types,              ONLY: cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_create,&
                                             dyn_coeff_set_release,&
                                             dyn_coeff_set_type
  USE gamma,                           ONLY: init_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE kg_energy_types,                 ONLY: allocate_kg_energy,&
                                             kg_energy_type
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type,&
                                             set_kg_env
  USE kg_environment_methods,          ONLY: kg_env_setup
  USE kg_force_types,                  ONLY: allocate_kg_force,&
                                             kg_force_type
  USE kg_rho_methods,                  ONLY: kg_rho_create
  USE kg_rho_types,                    ONLY: kg_rho_release,&
                                             kg_rho_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: allocate_molecule_kind_set,&
                                             atom_type,&
                                             get_molecule_kind,&
                                             molecule_kind_type,&
                                             set_molecule_kind,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: allocate_molecule_set,&
                                             molecule_type,&
                                             set_molecule,&
                                             set_molecule_set
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: allocate_particle_set,&
                                             particle_type,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates,&
                                             write_structure_data
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_type
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
  USE qs_environment_types,            ONLY: qs_environment_p_type,&
                                             qs_environment_type,&
                                             qs_env_create, &
                                             qs_env_release, &
                                             get_qs_env,&
                                             set_qs_env
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type
  USE qs_wf_history_methods,           ONLY: wfi_create
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_release
  USE scf_control_types,               ONLY: scf_c_create,&
                                             scf_c_read_parameters,&
                                             scf_c_release,&
                                             scf_c_write_parameters,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type,&
                                             read_cell,&
                                             write_cell
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kg_init

  CHARACTER(LEN=*), PARAMETER :: moduleN = "kg_environment"
  
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!****** kg_environment/kg_init [1.0] *
!!
!!   NAME
!!     kg_init
!!
!!   FUNCTION
!!   reads the input and database file for KG
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!****
!!*****************************************************************************

  SUBROUTINE kg_init (kg_env,globenv,error)

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "kg_init"

    INTEGER                                  :: handle, istat, maxl, maxlgto, &
                                                maxlppl, maxlppnl, nkind, &
                                                output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    LOGICAL                                  :: ionode
    LOGICAL                                  :: use_ref_cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(kg_energy_type), POINTER            :: energy
    TYPE(kg_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(scf_control_type), POINTER          :: scf_control

!------------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I",'',handle)

    ionode = globenv%ionode
    output_unit = globenv%scr

    ! nullifying  pointers
    NULLIFY ( atomic_kind_set )
    NULLIFY ( cell )
    NULLIFY ( dft_control )
    NULLIFY ( dyn_coeff_set )
    NULLIFY ( energy )
    NULLIFY ( force )
    NULLIFY ( local_molecules )
    NULLIFY ( local_particles )
    NULLIFY ( molecule_set )
    NULLIFY ( particle_set )
    NULLIFY ( scf_control )

    ! Read/write the input section with the dft and qs control parameters
    CALL read_dft_control ( dft_control, globenv, error )
    CALL read_qs_control ( dft_control % qs_control, globenv )
    CALL write_dft_control ( dft_control, globenv )
    CALL write_qs_control ( dft_control % qs_control, globenv )

    ! Read/write the input section with the cell parameters
    CALL read_cell ( cell, globenv, cell_ref=cell_ref, &
                     use_ref_cell=use_ref_cell )
    CALL write_cell ( cell, globenv )

    ! topology
    CALL topology_control (atomic_kind_set,particle_set,&
                           molecule_kind_set,molecule_set,&
                           globenv)

    ! Read/write the atomic kind set
    CALL read_atomic_kind_set(atomic_kind_set,globenv)
    CALL check_atomic_kind_set(atomic_kind_set,dft_control,globenv)
    CALL write_gto_basis_sets(atomic_kind_set,globenv)

    ! Initialize the spherical harmonics
    ! the orbital transformation matrices
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl)

    maxl = MAX(maxlgto,maxlppl,maxlppnl)
    IF (dft_control%forces) maxl = maxl + 1

    CALL init_orbital_pointers(maxl)
    CALL init_spherical_harmonics(maxl,globenv)

    ! Initialize the pretabulation for the calculation of the
    ! incomplete Gamma function F_n(t) after McMurchie-Davidson
    maxl = 3*maxlgto + 1
    CALL init_md_ftable(maxl)

    ! Initialise/print the atomic/molecular kind set
    CALL init_atomic_kind_set(atomic_kind_set,globenv)
    CALL write_atomic_kind_set(atomic_kind_set,globenv)
    CALL write_molecule_kind_set(molecule_kind_set,globenv)
    CALL write_total_numbers(atomic_kind_set,particle_set,globenv)

    ! Print the atomic coordinates 
    IF (globenv%print%atomic_coordinates) THEN
      CALL write_qs_particle_coordinates(particle_set,cell,globenv)
    END IF

    ! Print the interatomic distances ***
    CALL write_particle_distances(particle_set,cell,globenv)
    CALL write_structure_data(particle_set,cell,globenv)

    ! Initialize the atomic interaction radii 
    CALL init_interaction_radii ( dft_control%qs_control, cell, &
                                  atomic_kind_set, globenv)

    ! Distribute molecules and atoms using the new data structures
    CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules,&
                                 error=error)

    ! If polarization is also present, allocate and initialize coefs pointer  
    IF ( dft_control % qs_control % polarization ) THEN
       CALL dyn_coeff_set_create(dyn_coeff_set=dyn_coeff_set,&
                                 atomic_kind_set=atomic_kind_set,&
                                 distribution=local_particles,&
                                 error=error)
       CALL set_kg_env (kg_env=kg_env, dyn_coeff_set=dyn_coeff_set) 
       CALL dyn_coeff_set_release(dyn_coeff_set, error=error)
    END IF

    ! SCF parameters
    IF ( dft_control % qs_control % method /= "KG_NOPOL" ) THEN
      CALL scf_c_create(scf_control)
      CALL scf_c_read_parameters(scf_control,globenv)
      CALL scf_c_write_parameters(scf_control,globenv)
      CALL set_kg_env(kg_env=kg_env,scf_control=scf_control)
      CALL scf_c_release(scf_control,error=error)  
    END IF

    ! Allocate the data structure for energies ***
    CALL allocate_kg_energy(energy)

    ! Allocate the force data structure
    IF (dft_control%forces) THEN
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                       "natom_of_kind",nkind*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom_of_kind=natom_of_kind)
      CALL allocate_kg_force(force,natom_of_kind)
      DEALLOCATE (natom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,moduleN,__LINE__,&
                                       "natom_of_kind")
    END IF

    ! set the kg_env
    CALL set_kg_env ( kg_env=kg_env,&
                      atomic_kind_set=atomic_kind_set,&
                      cell=cell, &
                      cell_ref=cell_ref, &
                      use_ref_cell=use_ref_cell, &
                      dft_control=dft_control, &
                      energy=energy,&
                      force=force,&
                      local_molecules=local_molecules,&
                      local_particles=local_particles,&
                      molecule_kind_set=molecule_kind_set,&
                      molecule_set=molecule_set,&
                      particle_set=particle_set)  
    
    CALL cell_release(cell,error=error)
    CALL distribution_1d_release(local_particles,error=error)
    CALL distribution_1d_release(local_molecules,error=error)
    CALL dft_control_release(dft_control, error=error)

    !Sets up pw, grids, densities              
    CALL kg_env_setup(kg_env,error)

    !Sets up qs_env for subsystems             
    CALL get_kg_env ( kg_env=kg_env,&
                      dft_control=dft_control)
    IF ( dft_control % qs_control % method == "KG_GPW" ) THEN
      CALL kg_subsys_setup(kg_env,error)
    END IF

    IF (ionode) CALL m_flush(output_unit)

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE kg_init

!***************************************************************************

  SUBROUTINE write_total_numbers(atomic_kind_set,particle_set,globenv)

!   Purpose: Write the total number of kinds, atoms, etc. to the logical unit
!            number lunit.

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    INTEGER                                  :: maxlgto, maxlppl, maxlppnl, &
                                                natom, ncgf, nkind, npgf, &
                                                nset, nsgf, nshell, &
                                                output_unit

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%ionode) RETURN

    IF (globenv%print%total_numbers) THEN

      output_unit = globenv%scr

      natom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               maxlgto=maxlgto,&
                               maxlppl=maxlppl,&
                               maxlppnl=maxlppnl,&
                               ncgf=ncgf,&
                               npgf=npgf,&
                               nset=nset,&
                               nsgf=nsgf,&
                               nshell=nshell)

      WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
        "TOTAL NUMBERS AND MAXIMUM NUMBERS"

      WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T71,I10))")&
        "Total number of",&
        "- Atomic kinds:                  ",nkind,&
        "- Atoms:                         ",natom,&
        "- Shell sets:                    ",nset,&
        "- Shells:                        ",nshell,&
        "- Primitive Cartesian functions: ",npgf,&
        "- Cartesian basis functions:     ",ncgf,&
        "- Spherical basis functions:     ",nsgf

      IF ( maxlppnl > -1 ) THEN
        WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T75,I6))")&
          "Maximum angular momentum of the",&
          "- Orbital basis functions:                   ",maxlgto,&
          "- Local part of the GTH pseudopotential:     ",maxlppl,&
          "- Non-local part of the GTH pseudopotential: ",maxlppnl
      ELSE IF ( maxlppl > -1 ) THEN
        WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T75,I6))")&
          "Maximum angular momentum of the",&
          "- Orbital basis functions:                   ",maxlgto,&
          "- Local part of the GTH pseudopotential:     ",maxlppl
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T3,A,T75,I6))")&
          "Maximum angular momentum of the orbital basis functions: ",maxlgto
      END IF

    END IF

  END SUBROUTINE write_total_numbers

!******************************************************************************

  SUBROUTINE kg_subsys_setup(kg_env,error)

!   We generate a qs_environment for the system

    TYPE(kg_environment_type), INTENT(inout) :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_subsys_setup', &
      routineP = moduleN//':'//routineN

    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL                                  :: use_ref_cell
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force

    LOGICAL                                  :: failure
    INTEGER                                  :: nkind, istat
    INTEGER, DIMENSION(:), ALLOCATABLE       :: natom_of_kind

!   ---------------------------------------------------------------------------

    failure=.FALSE.

    NULLIFY(qs_env,subsys,scf_control,cell,cell_ref,atomic_kind_set,&
            para_env,dft_control,wf_history,qs_energy,qs_force)

    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    cell_ref=cell_ref,&
                    use_ref_cell=use_ref_cell,&
                    para_env=para_env,&
                    scf_control=scf_control,&
                    subsys=subsys,&
                    error=error)

    CALL qs_env_create(qs_env, para_env=para_env, error=error)

    CALL wfi_create(wf_history,&
                    interpolation_method_nr = &
                      dft_control%qs_control%wf_interpolation_method_nr,&
                    extrapolation_order = &
                      dft_control%qs_control%wf_extrapolation_order,&
                    error=error)

    CALL allocate_qs_energy ( qs_energy )

    IF ( dft_control%forces ) THEN
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "natom_of_kind",nkind*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom_of_kind=natom_of_kind)
      CALL allocate_qs_force ( qs_force,natom_of_kind )
      DEALLOCATE (natom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "natom_of_kind")
    END IF

    CALL set_qs_env(qs_env=qs_env,&
                    cell=cell                 ,&
                    cell_ref=cell_ref         ,&
                    use_ref_cell=use_ref_cell ,&
                    dft_control=dft_control   ,&
                    energy=qs_energy          ,&
                    force=qs_force            ,&
                    scf_control=scf_control   ,&
                    wf_history=wf_history     ,&
                    subsys=subsys             ,&
                    error=error                 )

    CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
    CALL set_qs_env(qs_env, qs_charges=qs_charges,error=error)
    CALL qs_charges_release(qs_charges,error=error)

    CALL qs_env_rebuild_pw_env(qs_env, error=error)

    CALL set_kg_env(kg_env=kg_env,&
                    sub_qs_env=qs_env)
    CALL qs_env_release(qs_env,error=error)

  END SUBROUTINE kg_subsys_setup

!******************************************************************************

END MODULE kg_environment

!******************************************************************************

