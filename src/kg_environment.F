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
  USE cp_output_handling,              ONLY: cp_print_key_should_output,&
                                             cp_p_file,&
                                             cp_print_key_unit_nr,&
                                             cp_print_key_finished_output,&
                                             cp_iter_string,&
                                             cp_iterate,&
                                             cp_add_iter_level,&
                                             cp_rm_iter_level
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
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
                                             read_mgrid_section,&
                                             read_qs_section,&
                                             write_dft_control,&
                                             write_qs_control
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_type
  USE cp_subsystem_methods,            ONLY: cp_subsys_read_colvar
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_create,&
                                             dyn_coeff_set_release,&
                                             dyn_coeff_set_type
  USE gamma,                           ONLY: init_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get
  USE kg_energy_types,                 ONLY: allocate_kg_energy,&
                                             kg_energy_type
  USE kg_environment_methods,          ONLY: kg_env_setup
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type,&
                                             set_kg_env
  USE kg_force_types,                  ONLY: allocate_kg_force,&
                                             kg_force_type
  USE kg_gpw_fm_mol_methods,           ONLY: build_local_fm_mol
  USE kg_gpw_fm_mol_types,             ONLY: kg_fm_mol_set_create,&
                                             kg_fm_mol_set_release,&
                                             kg_fm_mol_set_type,&
                                             kg_fm_p_type
  USE kg_gpw_pw_env_methods,           ONLY: build_molbox_env
  USE kg_gpw_pw_env_types,             ONLY: kg_molbox_env_type,&
                                             kg_sub_pw_env_create,&
                                             kg_sub_pw_env_release,&
                                             kg_sub_pw_env_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_max,mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             num_ao_el_per_molecule,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: molecule_type
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates,&
                                             write_structure_data
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_type
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_env_release,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE qs_wf_history_methods,           ONLY: wfi_create
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,wfi_release
  USE scf_control_types,               ONLY: scf_c_create,&
                                             scf_c_read_parameters,&
                                             scf_c_release,&
                                             scf_c_write_parameters,&
                                             scf_control_type
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_type,&
                                             read_cell,&
                                             write_cell,&
                                             get_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
  USE virial_types,                    ONLY: virial_create
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
#include "cp_common_uses.h"

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
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "kg_init"

    INTEGER                                  :: handle, istat, maxl, maxlgto, &
                                                maxlppl, maxlppnl, nkind, &
                                                output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    LOGICAL                                  :: ionode, kg_gpw, use_ref_cell, was_present
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
    TYPE(section_vals_type), POINTER         :: cell_section,&
                                                force_env_section,&
                                                poisson_section
    REAL(kind=dp), DIMENSION(3)              :: abc
    TYPE(cp_logger_type), pointer            :: logger

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I",'',handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    was_present = .FALSE.
    
    ! nullifying  pointers
    NULLIFY ( atomic_kind_set )
    NULLIFY ( cell_section, force_env_section, poisson_section)
    NULLIFY ( cell)
    NULLIFY ( dft_control )
    NULLIFY ( dyn_coeff_set )
    NULLIFY ( energy )
    NULLIFY ( force )
    NULLIFY ( local_molecules )
    NULLIFY ( local_particles )
    NULLIFY ( molecule_set )
    NULLIFY ( particle_set )
    NULLIFY ( scf_control )
    NULLIFY ( cell_ref )
    NULLIFY ( cell_section, force_env_section, poisson_section )
    NULLIFY ( logger )

    logger => cp_error_get_logger(error)
    force_env_section => section_vals_get_subs_vals(globenv%input_file,&
         "FORCE_EVAL",error=error)
    cell_section => section_vals_get_subs_vals(force_env_section,&
         "SUBSYS%CELL",error=error) 
    poisson_section => section_vals_get_subs_vals(force_env_section,&
         "DFT%POISSON_DFT",error=error)


    CALL set_kg_env(kg_env,input=force_env_section,&
         error=error)

    ! Read/write the input section with the dft and qs control parameters
    CALL read_dft_control ( dft_control, globenv, error )

!   *** Read the input section with the Quickstep control parameters ***
    CALL read_qs_section(dft_control%qs_control,globenv,error=error)
    CALL read_mgrid_section(dft_control%qs_control,globenv,error=error)

    CALL write_dft_control ( dft_control, globenv )
    CALL write_qs_control ( dft_control % qs_control, globenv )

!   *** Read the input section with the cell parameters
    CALL read_cell ( cell, cell_ref, globenv,  use_ref_cell=use_ref_cell,&
         cell_section=cell_section, error=error)
    CALL get_cell( cell, abc=abc)    

!   *** Print the cell parameters ***
    CALL write_cell ( cell, globenv )

    ! topology
    CALL topology_control (atomic_kind_set,particle_set,&
                           molecule_kind_set,molecule_set,&
                           globenv, error=error)

    ! Read/write the atomic kind set
    CALL read_atomic_kind_set(atomic_kind_set,globenv)
    CALL check_atomic_kind_set(atomic_kind_set,dft_control,globenv)
    CALL write_gto_basis_sets(atomic_kind_set,globenv)

    CALL num_ao_el_per_molecule(molecule_kind_set)


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
    CALL init_atomic_kind_set(atomic_kind_set,globenv,&
                              method=dft_control%qs_control%method)
    CALL write_atomic_kind_set(atomic_kind_set,globenv)
    CALL write_molecule_kind_set(molecule_kind_set,globenv)
    CALL write_total_numbers(atomic_kind_set,particle_set,globenv)

    ! Print the atomic coordinates 
    CALL write_qs_particle_coordinates(particle_set,cell,globenv,error)

    ! Print the interatomic distances ***
    CALL write_particle_distances(particle_set,cell,globenv,error)
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

    ! Flag for the KG_GPW method: if true the sub-qs_enviroment is created
    kg_gpw = ( dft_control % qs_control % method == "KG_GPW" )

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
    CALL cell_release(cell_ref,error=error)
    CALL distribution_1d_release(local_particles,error=error)
    CALL distribution_1d_release(local_molecules,error=error)
    CALL dft_control_release(dft_control, error=error)
 
    !Sets up pw, grids, densities              
    CALL kg_env_setup(kg_env,kg_gpw,error)

    IF ( kg_gpw ) THEN
      !Sets up qs_env for subsystems             
      CALL kg_subsys_setup(kg_env,globenv,error)
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
      pointer                             :: globenv

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

  SUBROUTINE kg_subsys_setup(kg_env,globenv,error)

!   We generate a qs_environment for the system

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(global_environment_type), &
      pointer                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_subsys_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imolecule_kind, istat, &
                                                nao_global, nao_max, nelectron_global,&
                                                nkind, nmolecule_kind, nmo_max, nspins
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    REAL(dp)                                 :: maxocc_global(2)
    LOGICAL                                  :: failure, nmo_eq_nao, use_ref_cell
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(kg_fm_p_type), POINTER              :: kg_fm
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(kg_sub_pw_env_type), POINTER        :: kg_sub_pw_env
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(scf_control_type), POINTER          :: scf_control
    INTEGER, DIMENSION(2)                    :: added_mos
    TYPE(section_vals_type), POINTER         :: poisson_section,&
                                                force_env_section,&
                                                colvar_section
    TYPE(cp_logger_type), pointer            :: logger
    
    failure=.FALSE.

    logger => cp_error_get_logger(error)
    NULLIFY(qs_env,subsys,scf_control,cell,cell_ref,atomic_kind_set,&
            para_env,dft_control,distribution_2d,wf_history,qs_energy,qs_force,&
            poisson_section)

    NULLIFY(distribution_1d,kg_fm,molecule_kind_set,molecule_set,particle_set)

    poisson_section => section_vals_get_subs_vals(globenv%input_file,&
         "FORCE_EVAL%DFT%POISSON_DFT",error=error)
    force_env_section => section_vals_get_subs_vals(globenv%input_file,&
         "FORCE_EVAL",error=error)

    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    cell_ref=cell_ref,&
                    distribution_2d=distribution_2d,&
                    use_ref_cell=use_ref_cell,&
                    para_env=para_env,&
                    scf_control=scf_control,&
                    subsys=subsys,&
                    error=error)

    ! read colvars and add to subsystem
    colvar_section => section_vals_get_subs_vals(globenv%input_file,&
                              "FORCE_EVAL%SUBSYS%COLVAR",error=error)
    CALL cp_subsys_read_colvar( subsys, colvar_section, error=error)


    CALL qs_env_create(qs_env, para_env=para_env, error=error)

    CALL virial_create ( qs_env % virial, virial_in_use = .FALSE., numer_virial = .FALSE. )

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


!   No idea whether OT can work with KG
    CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.FALSE.)

    CALL set_qs_env(qs_env=qs_env,&
                    cell=cell                 ,&
                    input=force_env_section   ,&
                    cell_ref=cell_ref         ,&
                    use_ref_cell=use_ref_cell ,&
                    dft_control=dft_control   ,&
                    energy=qs_energy          ,&
                    force=qs_force            ,&
                    scf_control=scf_control   ,&
                    subsys=subsys             ,&
                    wf_history=wf_history     ,&
                    error=error                 )

    CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
    CALL set_qs_env(qs_env, qs_charges=qs_charges,error=error)
    CALL qs_charges_release(qs_charges,error=error)

    CALL wfi_release(wf_history,error=error)

    CALL qs_env_rebuild_pw_env(qs_env, error=error)

!MI
!   Here add the allocations and initialization of the separated blocks of the fm
!   They should be arguments of the kg_env
    CALL get_kg_env(kg_env=kg_env,&
                    kg_fm_set=kg_fm,&
                    local_molecules=distribution_1d,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    particle_set=particle_set)

    nmolecule_kind = SIZE(molecule_kind_set,1)

    nspins = dft_control%nspins
 
!   *** Problems with the distribution of the electrons ***

!   I suppose, the charge of each molecule should be given from input
!    to assign the correct number of electrons. Is this possible?
!   It cannot be the global charge of the system, because each molecule
!   must have a fixed number of electron from the beginning, 
!   and charge transfers are not admitted
    IF(dft_control%charge /= 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
         "Automatic collocation of the charge on molecules not implemented")
    END IF
!   the multiplicity is given by the even or odd number of electrons
!   for each molecule separatedly. I don't whether it is possible to do differently
    IF (dft_control%multiplicity /= 0) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
         "The multiplicity is assigned molecule by molecule")
    END IF
    
!   *** Some options require that all MOs are computed ... ***
    nmo_eq_nao = .FALSE.
    IF (BTEST(cp_print_key_should_output(logger%iter_info,globenv%input_file,&
            "FORCE_EVAL%DFT%PRINT%MO",error=error),cp_p_file).OR.&
        (scf_control%level_shift /= 0.0_dp).OR.&
        (scf_control%smear /= 0.0_dp).OR.&
        (scf_control%eps_jacobi /= 0.0_dp)) THEN
        nmo_eq_nao = .TRUE.
    END IF

!   *** In case of Added_mos /= 0 the same number of states is added to each molecule
    added_mos(1:2) = scf_control%added_mos(1:2)

    CALL kg_fm_mol_set_create(kg_fm,nmolecule_kind,error=error)
    CALL kg_sub_pw_env_create(kg_sub_pw_env,nmolecule_kind,error=error)
    nelectron_global = 0
    nao_global = 0
    nao_max = 0
    nmo_max = 0
    maxocc_global(1:2) = 0.0_dp
    DO imolecule_kind = 1,nmolecule_kind
      molecule_kind => molecule_kind_set(imolecule_kind)
      CALL build_local_fm_mol(kg_fm%kg_fm_mol_set(imolecule_kind),&
                              molecule_kind,&
                              imolecule_kind,distribution_1d,&
                              molecule_set,particle_set,&
                              nspins,nmo_eq_nao,added_mos,&
                              scf_control%use_cholesky,&
                              nelectron_global,nao_global,&
                              nao_max, nmo_max, &
                              maxocc_global=maxocc_global,&
                              error=error)
      ! Initialize the local cell and the pw_env for each molecule kind
      fm_mol_set => kg_fm%kg_fm_mol_set(imolecule_kind) 
      CALL build_molbox_env(kg_sub_pw_env%molbox_env_set(imolecule_kind),&
                            fm_mol_set,qs_env,poisson_section=poisson_section,&
                            error=error)
    END DO

    CALL set_qs_env(qs_env=qs_env,kg_sub_pw_env=kg_sub_pw_env)
    CALL kg_sub_pw_env_release(kg_sub_pw_env,error=error)

    CALL mp_sum(nelectron_global, para_env%group)
    CALL mp_sum(nao_global, para_env%group)
    CALL mp_max(nao_max, para_env%group)
    CALL mp_max(nmo_max, para_env%group)
    kg_fm%nelectron_global = nelectron_global
    kg_fm%nao_global = nao_global
    kg_fm%nao_max = nao_max
    kg_fm%nmo_max = nmo_max
    kg_fm%maxocc_global = maxocc_global
    kg_fm%nmolecule_global = SIZE(molecule_set,1)
    CALL set_kg_env(kg_env=kg_env,&
                    kg_fm_set=kg_fm,&
                    sub_qs_env=qs_env)
    CALL qs_env_release(qs_env,error=error)
    CALL kg_fm_mol_set_release(kg_fm,error=error)

  END SUBROUTINE kg_subsys_setup

!******************************************************************************

END MODULE kg_environment

!******************************************************************************

