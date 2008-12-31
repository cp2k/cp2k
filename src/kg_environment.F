!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief initialize kg environment
!> \par History
!>      GT 11-13-2003 :
!>      moved the main driver to cp2k
!>      performs only the initialization of the environment
!> \author gt SEPT-23-2002
! *****************************************************************************
MODULE kg_environment
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE cell_types,                      ONLY: cell_release,&
                                             cell_type,&
                                             get_cell,&
                                             read_cell,&
                                             write_cell
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_type
  USE cp_control_utils,                ONLY: read_dft_control,&
                                             read_mgrid_section,&
                                             read_qs_section,&
                                             write_dft_control,&
                                             write_qs_control
  USE cp_ddapc_types,                  ONLY: cp_ddapc_ewald_create
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_methods,               ONLY: cp_subsys_read_colvar
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE cp_symmetry,                     ONLY: write_symmetry
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_create,&
                                             dyn_coeff_set_initialize,&
                                             dyn_coeff_set_release,&
                                             dyn_coeff_set_type
  USE f77_blas
  USE gamma,                           ONLY: init_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: kg_header
  USE input_constants,                 ONLY: xc_vdw_fun_pairpot
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kg_energy_types,                 ONLY: allocate_kg_energy,&
                                             kg_energy_type
  USE kg_environment_methods,          ONLY: kg_env_setup
  USE kg_environment_types,            ONLY: get_kg_env,&
                                             kg_environment_type,&
                                             set_kg_env
  USE kg_force_types,                  ONLY: allocate_kg_force,&
                                             kg_force_type
  USE kg_gpw_fm_mol_methods,           ONLY: build_local_fm_mol
  USE kg_gpw_fm_mol_types,             ONLY: get_kg_fm_mol_set,&
                                             kg_fm_mol_set_create,&
                                             kg_fm_mol_set_release,&
                                             kg_fm_mol_set_type,&
                                             kg_fm_p_type
  USE kg_gpw_pw_env_methods,           ONLY: build_molbox_env
  USE kg_gpw_pw_env_types,             ONLY: kg_sub_pw_env_create,&
                                             kg_sub_pw_env_release,&
                                             kg_sub_pw_env_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_max,&
                                             mp_sum,&
                                             mp_sync
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
  USE qs_dispersion_types,             ONLY: qs_dispersion_type
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
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_release
  USE rel_control_types,               ONLY: rel_c_create,&
                                             rel_c_read_parameters,&
                                             rel_c_release,&
                                             rel_control_type
  USE scf_control_types,               ONLY: scf_c_create,&
                                             scf_c_read_parameters,&
                                             scf_c_release,&
                                             scf_c_write_parameters,&
                                             scf_control_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
  USE virial_types,                    ONLY: virial_create
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kg_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_environment'

!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief reads the input and database file for KG
! *****************************************************************************
  SUBROUTINE kg_init (kg_env,root_section,para_env,globenv,force_env_section,&
       subsys_section, use_motion_section, error)

    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, iw, maxl, &
                                                maxlgto, maxlppl, maxlppnl, &
                                                nkind, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    LOGICAL                                  :: ionode, kg_gpw, use_ref_cell, &
                                                was_present
    REAL(kind=dp), DIMENSION(3)              :: abc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
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
    TYPE(section_vals_type), POINTER         :: cell_section, colvar_section, &
                                                dft_section, kind_section, &
                                                qs_section, work_section

    CALL timeset(routineN,handle)
    ionode = para_env%ionode
    logger => cp_error_get_logger(error)
    IF(ionode) output_unit = cp_logger_get_default_unit_nr(logger)
    was_present = .FALSE.

    ! nullifying  pointers
    NULLIFY ( atomic_kind_set )
    NULLIFY ( cell_section )
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
    NULLIFY ( cell_section, kind_section )
    NULLIFY ( logger )

    logger => cp_error_get_logger(error)
    IF (.NOT.ASSOCIATED(subsys_section)) THEN
       subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    END IF
    cell_section => section_vals_get_subs_vals(subsys_section,"CELL",error=error)
    dft_section =>  section_vals_get_subs_vals(force_env_section,"DFT",error=error)
    qs_section =>  section_vals_get_subs_vals(dft_section,"QS",error=error)

    iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PROGRAM_BANNER",&
         extension=".Log",error=error)
    CALL kg_header(iw, iw>0)
    CALL cp_print_key_finished_output(iw,logger,dft_section,&
         "PRINT%PROGRAM_BANNER",error=error)

    CALL set_kg_env(kg_env,input=force_env_section,&
         error=error)

    ! Read/write the input section with the dft and qs control parameters
    CALL read_dft_control ( dft_control, dft_section,error )

!   *** Read the input section with the Quickstep control parameters ***
    CALL read_qs_section(dft_control%qs_control,qs_section,error=error)
    dft_control%qs_control%becke_restraint=.FALSE.
    CALL read_mgrid_section(dft_control%qs_control,dft_section,error=error)

    CALL write_dft_control ( dft_control, dft_section, error)
    CALL write_qs_control ( dft_control % qs_control, dft_section, error)

    ! read colvars and add to subsys
    colvar_section => section_vals_get_subs_vals(subsys_section,"COLVAR",error=error)
    CALL cp_subsys_read_colvar( kg_env%subsys, colvar_section, error=error)

!   *** Read the input section with the cell parameters
    CALL read_cell ( cell, cell_ref, use_ref_cell=use_ref_cell,&
         cell_section=cell_section, para_env=para_env, error=error)
    CALL get_cell( cell, abc=abc)

!   *** Print the cell parameters ***
    CALL write_cell ( cell, subsys_section, error=error)

    ! topology
    CALL topology_control (atomic_kind_set,particle_set, molecule_kind_set,molecule_set,&
         kg_env%subsys%colvar_p, kg_env%subsys%gci, root_section,para_env, &
         force_env_section=force_env_section, subsys_section=subsys_section, &
         use_motion_section=use_motion_section, error=error)

    ! Read/write the atomic kind set
    kind_section => section_vals_get_subs_vals(subsys_section,"KIND",error=error)
    CALL read_atomic_kind_set(atomic_kind_set,kind_section,para_env,force_env_section,error)
    CALL write_gto_basis_sets(atomic_kind_set,subsys_section,error)

    CALL num_ao_el_per_molecule(molecule_kind_set)

    ! Initialize the spherical harmonics
    ! the orbital transformation matrices
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl)

    maxl = MAX(maxlgto,maxlppl,maxlppnl) + 1

    CALL init_orbital_pointers(maxl)
    CALL init_spherical_harmonics(maxl,root_section,error)

    ! Initialize the pretabulation for the calculation of the
    ! incomplete Gamma function F_n(t) after McMurchie-Davidson
    maxl = 3*maxlgto + 1
    CALL init_md_ftable(maxl)

    ! Initialise/print the atomic/molecular kind set
    CALL init_atomic_kind_set(atomic_kind_set,para_env,force_env_section=force_env_section,&
                              method=dft_control%qs_control%method,error=error)
    CALL write_atomic_kind_set(atomic_kind_set,subsys_section,error)
    CALL write_molecule_kind_set(molecule_kind_set,subsys_section,error)
    CALL write_total_numbers(atomic_kind_set,particle_set,kg_env%input,error)

    ! Print the atomic coordinates
    CALL write_qs_particle_coordinates(particle_set,cell,subsys_section,label="KIM-GORDON",error=error)

    ! Print the interatomic distances ***
    CALL write_particle_distances(particle_set,cell,subsys_section,error)

    ! Print the requested structure data
    CALL write_structure_data(particle_set,cell,subsys_section,error)

    ! Print symmetry information
    CALL write_symmetry(particle_set,cell,subsys_section,error)

    ! Initialize the atomic interaction radii
    CALL init_interaction_radii ( dft_control%qs_control, cell, &
                                  atomic_kind_set, subsys_section, error)

    ! Distribute molecules and atoms using the new data structures
    CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules,&
                                 force_env_section=kg_env%input,&
                                 error=error)

    ! If polarization is also present, allocate and initialize coefs pointer
    IF ( dft_control % qs_control % polarization ) THEN
       CALL dyn_coeff_set_create(dyn_coeff_set=dyn_coeff_set,&
                                 atomic_kind_set=atomic_kind_set,&
                                 distribution=local_particles,&
                                 error=error)
       work_section => section_vals_get_subs_vals(root_section,"MOTION%MD%KG_COEFF",error=error)
       CALL dyn_coeff_set_initialize(dyn_coeff_set,root_section,para_env,globenv,work_section,skipvel=.TRUE.,error=error)
       CALL set_kg_env (kg_env=kg_env, dyn_coeff_set=dyn_coeff_set,error=error)
       CALL dyn_coeff_set_release(dyn_coeff_set, error=error)
    END IF

    ! Flag for the KG_GPW method: if true the sub-qs_enviroment is created
    kg_gpw = ( dft_control % qs_control % method == "KG_GPW" )

    ! SCF parameters
    IF ( dft_control % qs_control % method /= "KG_NOPOL" ) THEN
      CALL scf_c_create(scf_control,error=error)
      CALL scf_c_read_parameters(scf_control,dft_section,error=error)
      CALL scf_c_write_parameters(scf_control,dft_section,error=error)
      CALL set_kg_env(kg_env=kg_env,scf_control=scf_control,error=error)
      CALL scf_c_release(scf_control,error=error)
    END IF

    ! Allocate the data structure for energies ***
    CALL allocate_kg_energy(energy)

    ! Allocate the force data structure
    nkind = SIZE(atomic_kind_set)
    ALLOCATE (natom_of_kind(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "natom_of_kind",nkind*int_size)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                           natom_of_kind=natom_of_kind)
    CALL allocate_kg_force(force,natom_of_kind)
    DEALLOCATE (natom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "natom_of_kind")

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
                      particle_set=particle_set,error=error)

    CALL cell_release(cell,error=error)
    CALL cell_release(cell_ref,error=error)
    CALL distribution_1d_release(local_particles,error=error)
    CALL distribution_1d_release(local_molecules,error=error)

    !Sets up pw, grids, densities
    CALL kg_env_setup(kg_env,kg_gpw,error)

    IF ( kg_gpw ) THEN
      !Sets up qs_env for subsys
      CALL kg_subsys_setup(kg_env,force_env_section,subsys_section,error)
    ELSE
      CALL dft_control_release(dft_control, error=error)
    END IF

    IF (ionode) CALL m_flush(output_unit)

    CALL timestop(handle)

  END SUBROUTINE kg_init

! *****************************************************************************
!> \brief   Write the total number of kinds, atoms, etc. to the logical unit
!>          number lunit.
!> \author  Matthias Krack
!> \date    06.10.2000
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_total_numbers(atomic_kind_set,particle_set,force_env_section,error)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: maxlgto, maxlppl, maxlppnl, &
                                                natom, ncgf, nkind, npgf, &
                                                nset, nsgf, nshell, &
                                                output_unit
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_print_key_unit_nr(logger,force_env_section,"PRINT%TOTAL_NUMBERS",&
         extension=".Log",error=error)

    IF (output_unit>0) THEN
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
        WRITE (UNIT=output_unit,FMT="(/,T3,A,T75,I6)")&
          "Maximum angular momentum of the orbital basis functions: ",maxlgto
      END IF

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "PRINT%TOTAL_NUMBERS",error=error)

  END SUBROUTINE write_total_numbers

! *****************************************************************************
!> \brief Distribute the fm_mol_set and allocate the mos for each molecule
!>       Distribute the electrons among the molecules
!>       Generate a qs_environment for the system
!> \note
!>       Charged molecule kinds need to be specified here, reading the info from input somehow
!> \author MI
! *****************************************************************************
  SUBROUTINE kg_subsys_setup(kg_env,force_env_section,subsys_section,error)
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_subsys_setup', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: molch_name
    INTEGER :: ich, iglobal, imol, imolecule_kind, istat, molch_kind, &
      molch_val, nao_global, nao_max, nelectron_global, nkind, nmo_max, &
      nmol_local, nmolcharge, nmolecule_kind, nrep, nspins, sumcharge
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    INTEGER, DIMENSION(2)                    :: added_mos
    INTEGER, DIMENSION(:), POINTER           :: charge_x_mol, mult_x_mol
    LOGICAL                                  :: explicit, failure, &
                                                nmo_eq_nao, use_ref_cell
    REAL(dp)                                 :: maxocc_global(2)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: distribution_1d
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(kg_fm_mol_set_type), POINTER        :: fm_mol_set
    TYPE(kg_fm_p_type), POINTER              :: kg_fm
    TYPE(kg_sub_pw_env_type), POINTER        :: kg_sub_pw_env
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, &
                                                mol_charge_section, &
                                                poisson_section, vdw_section, &
                                                xc_fun_section

    failure=.FALSE.

    logger => cp_error_get_logger(error)
    NULLIFY(qs_env,subsys,scf_control,cell,cell_ref,atomic_kind_set,&
            para_env,dft_control,distribution_2d,wf_history,qs_energy,qs_force,&
            poisson_section)
    NULLIFY(mol_charge_section, charge_x_mol)

    NULLIFY(distribution_1d,kg_fm,molecule_kind_set,molecule_set,particle_set)

    poisson_section => section_vals_get_subs_vals(force_env_section,&
         "DFT%POISSON",error=error)

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

    CALL qs_env_create(qs_env, para_env=para_env, error=error)

    CALL virial_create ( qs_env % virial, error=error)

    CALL wfi_create(wf_history,&
                    interpolation_method_nr = &
                    dft_control%qs_control%wf_interpolation_method_nr,&
                    extrapolation_order = &
                    dft_control%qs_control%wf_extrapolation_order,&
                    error=error)

    CALL allocate_qs_energy ( qs_energy )

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

    !   Create relativistic control section

    NULLIFY(dft_section)
    dft_section =>  section_vals_get_subs_vals(kg_env%input,"DFT",error=error)
    CALL rel_c_create(rel_control,error=error)
    CALL rel_c_read_parameters(rel_control,dft_section,error=error)
    CALL set_qs_env(qs_env,rel_control=rel_control,error=error)
    CALL rel_c_release(rel_control,error=error)

!   No idea whether OT can work with KG
    CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.FALSE.,error=error)

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

    !   *** Setup the grids for the G-space Interpolation if any
    CALL cp_ddapc_ewald_create(qs_env%cp_ddapc_ewald, .FALSE., .FALSE.,&
         cell, force_env_section, subsys_section, para_env, error)

    CALL dft_control_release(dft_control, error=error)
    CALL wfi_release(wf_history,error=error)

    CALL get_kg_env(kg_env=kg_env, dft_control=dft_control, error=error)
    CALL qs_charges_create(qs_charges,nspins=dft_control%nspins,error=error)
    CALL set_qs_env(qs_env, qs_charges=qs_charges,error=error)
    CALL qs_charges_release(qs_charges,error=error)

    ! setup the dispersion part (stop not implemented)
    ALLOCATE(dispersion_env,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(vdw_section,xc_fun_section)
    vdw_section =>  section_vals_get_subs_vals(dft_section,"XC%vdw_potential",error=error)
    xc_fun_section => section_vals_get_subs_vals(dft_section,"XC%XC_FUNCTIONAL",error=error)
    CALL section_vals_val_get(vdw_section, "POTENTIAL_TYPE", i_val=dispersion_env%type, error=error)
    IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN
      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END IF
    CALL set_qs_env(qs_env, dispersion_env=dispersion_env, error=error)

    CALL qs_env_rebuild_pw_env(qs_env, error=error)

!   Here add the allocations and initialization of the separated blocks of the fm
!   They should be arguments of the kg_env
    CALL get_kg_env(kg_env=kg_env,&
                    kg_fm_set=kg_fm,&
                    local_molecules=distribution_1d,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    particle_set=particle_set,error=error)

    nmolecule_kind = SIZE(molecule_kind_set,1)

    nspins = dft_control%nspins

!   *** Some options require that all MOs are computed ... ***
    CALL get_kg_env(kg_env=kg_env, scf_control=scf_control, error=error)
    nmo_eq_nao = .FALSE.
    IF (BTEST(cp_print_key_should_output(logger%iter_info,force_env_section,&
            "DFT%PRINT%MO",error=error),cp_p_file).OR.&
        (scf_control%level_shift /= 0.0_dp).OR.&
        (scf_control%smear%window_size /= 0.0_dp).OR.&
        (scf_control%diagonalization%eps_jacobi /= 0.0_dp)) THEN
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

    kg_fm%nmolecule_global = SIZE(molecule_set,1)
    ALLOCATE(kg_fm%imol_pe_pos(kg_fm%nmolecule_global),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    kg_fm%imol_pe_pos = 0
    ALLOCATE(kg_fm%imol_local_name(kg_fm%nmolecule_global),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    kg_fm%imol_local_name = 0

!   read from input the molecular charges if given 
    sumcharge = 0
    ALLOCATE(charge_x_mol(nmolecule_kind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    charge_x_mol = 0
    ALLOCATE(mult_x_mol(nmolecule_kind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    mult_x_mol = 0

    mol_charge_section => section_vals_get_subs_vals(force_env_section,&
         "DFT%KG%MOL_CHARGE",error=error)
    CALL section_vals_get(mol_charge_section,explicit=explicit,n_repetition=nmolcharge,error=error)
    IF(explicit) THEN
      DO ich =  1, nmolcharge
        molch_kind = 0
        molch_name = "UNDEF"
        CALL section_vals_val_get(mol_charge_section,"MOLECULE",i_rep_section=ich,&
             n_rep_val=nrep,error=error)
        IF (nrep/=0) THEN
            CALL section_vals_val_get(mol_charge_section,"MOLECULE",i_rep_section=ich,&
                 i_val=molch_kind,error=error)
        END IF

        CALL section_vals_val_get(mol_charge_section,"MOLNAME",i_rep_section=ich,&
                n_rep_val=nrep,error=error)
        IF (nrep/=0) THEN
            CALL section_vals_val_get(mol_charge_section,"MOLNAME",i_rep_section=ich,&
                   c_val=molch_name,error=error)
        END IF
        IF  (((molch_kind==0).AND.(molch_name=="UNDEF")).OR.&
                ((molch_kind/=0).AND.(molch_name/="UNDEF"))) THEN
           CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF

        CALL section_vals_val_get(mol_charge_section,"CHARGE",i_rep_section=ich,&
                i_val=molch_val,error=error)

        IF(molch_kind/=0) THEN
          imolecule_kind = molch_kind
          charge_x_mol(imolecule_kind) = molch_val
        ELSE
          DO imolecule_kind = 1, nmolecule_kind
             molecule_kind => molecule_kind_set(imolecule_kind)
             IF (molecule_kind%name==molch_name) THEN
               charge_x_mol(imolecule_kind) = molch_val
             END IF
          END DO 
        END IF
      END DO 

      CALL num_ao_el_per_molecule(molecule_kind_set,charge_x_mol=charge_x_mol)

    END IF

    ! Build the full matrix for each molecule
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
                              mol_charge=charge_x_mol(imolecule_kind),&
                              mol_multiplicity=mult_x_mol(imolecule_kind),&
                              error=error)
      ! Initialize the local cell and the pw_env for each molecule kind
      fm_mol_set => kg_fm%kg_fm_mol_set(imolecule_kind)

      CALL get_kg_fm_mol_set(kg_fm_mol_set=fm_mol_set,&
          nmolecule_local=nmol_local)
      IF(nmol_local>0) THEN 

        CALL build_molbox_env(kg_sub_pw_env%molbox_env_set(imolecule_kind),&
                              fm_mol_set,qs_env,poisson_section=poisson_section,&
                              error=error)
      END IF
      CALL mp_sync(para_env%group)
   
      DO imol = 1,distribution_1d%n_el(imolecule_kind)
        iglobal = distribution_1d%list(imolecule_kind)%array(imol)
        kg_fm%imol_local_name(iglobal) = imol
        kg_fm%imol_pe_pos(iglobal) = para_env%mepos
      END DO

!      sumcharge = sumcharge + molecule_kind%nmolecule*charge_x_mol(imolecule_kind) 
      sumcharge = sumcharge + nmol_local*charge_x_mol(imolecule_kind) 

    END DO

!   check that the total charge is correct
    CALL mp_sum(sumcharge,para_env%group)
    IF(dft_control%charge /= sumcharge) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
         " The charge distributed over the molecules is not equal to the total charge")
    END IF

!   deallocate temporary arrays
    DEALLOCATE(charge_x_mol, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(mult_x_mol, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL set_qs_env(qs_env=qs_env,kg_sub_pw_env=kg_sub_pw_env,error=error)
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
    CALL mp_sum(kg_fm%imol_local_name, para_env%group)
    CALL mp_sum(kg_fm%imol_pe_pos, para_env%group)

    CALL set_kg_env(kg_env=kg_env,&
                    kg_fm_set=kg_fm,&
                    sub_qs_env=qs_env,error=error)
    CALL qs_env_release(qs_env,error=error)
    CALL kg_fm_mol_set_release(kg_fm,error=error)

  END SUBROUTINE kg_subsys_setup

END MODULE kg_environment

