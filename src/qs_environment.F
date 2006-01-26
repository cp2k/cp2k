!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2004  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_environment [1.0] *
!!
!!   NAME
!!     qs_environment
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
!!     - Added more LSD structure (01.2003,Joost VandeVondele)
!!     - New molecule data types introduced (Sep. 2003,MK)
!!     - Cleaning; getting rid of pnode (02.10.2003,MK)
!!     - Sub-system setup added (08.10.2003,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_environment
  USE cp_output_handling,              ONLY: cp_print_key_should_output,&
                                             cp_p_file,&
                                             cp_print_key_unit_nr,&
                                             cp_print_key_finished_output,&
                                             cp_iter_string,&
                                             cp_iterate,&
                                             cp_add_iter_level,&
                                             cp_rm_iter_level
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             check_atomic_kind_set,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             init_gapw_basis_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_type,&
                                             gapw_control_type,&
                                             qs_control_type,&
                                             semi_empirical_control_type
  USE cp_control_utils,                ONLY: read_dft_control,&
                                             read_mgrid_section,&
                                             read_qs_section,&
                                             read_tddfpt_control,&
                                             write_dft_control,&
                                             write_qs_control
  USE cp_fm_struct,                    ONLY: optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_type
  USE cp_subsystem_methods,            ONLY: cp_subsys_read_colvar
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE gamma,                           ONLY: init_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE harris_energy_types,             ONLY: harris_energy_type, &
                                             harris_energy_create
  USE harris_env_types,                ONLY: harris_env_type, &
                                             harris_env_create, &
                                             harris_env_set
  USE harris_force_types,              ONLY: harris_force_type, &
                                             harris_force_create
  USE hartree_local_methods,           ONLY: init_coulomb_local
  USE input_constants,                 ONLY: do_method_gapw,&
                                             do_method_gapw_xc,&
                                             state_loc_all,&
                                             state_loc_none,&
                                             xas_tp_hh, xas_tp_fh
  USE header,                          only: qs_header
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_val_get, &
                                             section_vals_get,&
                                             section_vals_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             num_ao_el_per_molecule,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: molecule_type
  USE nddo_methods,                    ONLY: init_rho0_nddo
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_setup
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE qs_loc_control,                  ONLY: localized_wfn_control_type,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             mo_set_p_type
  USE qs_rho0_ggrid,                   ONLY: rho0_s_grid_create
  USE qs_rho0_methods,                 ONLY: init_rho0
  USE qs_rho0_types,                   ONLY: rho0_mpole_type
  USE qs_rho_atom_methods,             ONLY: init_rho_atom
  USE qs_wf_history_methods,           ONLY: wfi_create
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_release
  USE scf_control_types,               ONLY: scf_c_create,&
                                             scf_c_read_parameters,&
                                             scf_c_release,&
                                             scf_c_write_parameters,&
                                             scf_control_type
  USE cell_types,                 ONLY: cell_release,&
                                             cell_type,&
                                             read_cell,&
                                             write_cell,&
                                             get_cell,&
                                             cell_create,&
                                             cell_clone
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
  USE xas_control,                     ONLY: write_xas_control
  USE cp_ddapc_types,                  ONLY: cp_ddapc_ewald_create
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_environment"

! *** Public subroutines ***

  PUBLIC :: qs_init, qs_init_subsys

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_init(qs_env,globenv,subsystem,cell,cell_ref,qmmm,error)

!   Purpose: Read the input and the database files for the setup of the
!            QUICKSTEP environment.

!   History: - Creation (22.05.2000,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER         :: qs_env
    TYPE(global_environment_type), POINTER     :: globenv
    TYPE(cp_subsystem_type), POINTER, OPTIONAL :: subsystem
    TYPE(cell_type), POINTER, OPTIONAL         :: cell, cell_ref
    LOGICAL, INTENT(IN), OPTIONAL              :: qmmm
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                              :: error

    CHARACTER(len=*), PARAMETER :: routineN = "qs_init", &
      routineP = module_name//":"//routineN

    LOGICAL                                  :: use_ref_cell, &
                                                harris_flag
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: my_cell, my_cell_ref
    TYPE(cp_subsystem_type), POINTER         :: my_subsystem
    TYPE(mol_kind_new_list_type), POINTER    :: mol_kinds
    TYPE(mol_new_list_type), POINTER         :: mols
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: force_env_section,&
                                                cell_section,&
                                                harris_section,&
                                                colvar_section
    REAL(kind=dp), dimension(3)              :: abc
    LOGICAL                                  :: my_qmmm

    NULLIFY(my_subsystem, molecule_kind_set, molecule_set, atomic_kind_set,&
         particle_set,atomic_kinds,particles,mols,mol_kinds,my_cell,my_cell_ref,&
         cell_section, harris_section)
    
    force_env_section => section_vals_get_subs_vals(globenv%input_file,&
         "FORCE_EVAL",error=error)
    cell_section => section_vals_get_subs_vals(force_env_section,&
         "SUBSYS%CELL",error=error)
    harris_section => section_vals_get_subs_vals(force_env_section, &
         "DFT%QS%HARRIS", error=error)
    
    my_qmmm     = .FALSE.
    IF (PRESENT(qmmm)) my_qmmm=qmmm

    harris_flag = .FALSE.
    CALL section_vals_val_get(harris_section, "ACTIVATE", &
                              l_val=harris_flag, error=error)

    CALL set_qs_env(qs_env,input=force_env_section,&
         use_harris=harris_flag,error=error)

    IF (PRESENT(cell)) THEN
       my_cell => cell
       IF (PRESENT(cell_ref)) THEN
          my_cell_ref => cell_ref
          use_ref_cell = .TRUE.
       ELSE
          CALL cell_create(my_cell_ref,error=error)
          CALL cell_clone(my_cell,my_cell_ref,error)
          use_ref_cell = .FALSE.
       END IF
    ELSE
       !   *** Read the input section with the cell parameters ***
       CALL read_cell(my_cell, my_cell_ref, globenv, use_ref_cell=use_ref_cell, &
                      cell_section=cell_section, error=error)
    END IF
    CALL get_cell ( my_cell, abc=abc)

!   *** Setup the grids for the G-space Interpolation if any
    CALL cp_ddapc_ewald_create(qs_env%cp_ddapc_ewald, my_qmmm, my_cell, force_env_section,&
         globenv, error)

!   *** Print the cell parameters ***

    CALL write_cell(my_cell,force_env_section,error)
    IF (PRESENT(subsystem)) THEN 
       my_subsystem => subsystem
       CALL read_atomic_kind_set(my_subsystem%atomic_kinds%els, globenv,error)
    ELSE
       CALL cp_subsys_create(my_subsystem,para_env=globenv%para_env,error=error)

       colvar_section => section_vals_get_subs_vals(globenv%input_file,&
                                 "FORCE_EVAL%SUBSYS%COLVAR",error=error)
       CALL cp_subsys_read_colvar( my_subsystem, colvar_section, error=error)


       !   *** Read the particle coordinates and allocate the atomic kind, ***
       !   *** the molecule kind, and the molecule data structures         ***
       CALL topology_control (atomic_kind_set,particle_set,&
                              molecule_kind_set,molecule_set,&
                              globenv,error=error)
  
       CALL read_atomic_kind_set(atomic_kind_set,globenv,error)
       
       CALL num_ao_el_per_molecule(molecule_kind_set)
       
       CALL particle_list_create(particles,els_ptr=particle_set,error=error)
       CALL atomic_kind_list_create(atomic_kinds,els_ptr=atomic_kind_set,error=error)
       CALL mol_new_list_create(mols,els_ptr=molecule_set,error=error)
       CALL mol_kind_new_list_create(mol_kinds,els_ptr=molecule_kind_set,error=error)
       CALL cp_subsys_set(my_subsystem,particles=particles,atomic_kinds=atomic_kinds,&
            molecules_new=mols,molecule_kinds_new=mol_kinds,error=error)
       CALL particle_list_release(particles,error=error)
       CALL atomic_kind_list_release(atomic_kinds,error=error)
       CALL mol_new_list_release(mols,error=error)
       CALL mol_kind_new_list_release(mol_kinds,error=error)
    END IF

    CALL qs_init_subsys(qs_env,my_subsystem,my_cell,my_cell_ref,use_ref_cell,&
         globenv,my_qmmm,harris=harris_flag,error=error)

    IF (.NOT.PRESENT(cell))      CALL cell_release(my_cell,error=error)
    IF (.NOT.PRESENT(cell_ref))  CALL cell_release(my_cell_ref,error=error)
    IF (.NOT.PRESENT(subsystem)) CALL cp_subsys_release(my_subsystem,error=error)

  END SUBROUTINE qs_init

! *****************************************************************************

  SUBROUTINE qs_init_subsys(qs_env,subsystem,cell,cell_ref,use_ref_cell,&
       globenv,qmmm,harris,error)

!   Purpose: Read the input and the database files for the setup of the
!            QUICKSTEP environment.

!   History: - Creation (22.05.2000,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL, INTENT(in)                      :: use_ref_cell
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(in)                      :: qmmm
    LOGICAL, INTENT(IN), OPTIONAL            :: harris
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "qs_init_subsys"

    INTEGER :: handle, ispin, istat, lmax_sphere, maxl, maxlgto, maxlppl, &
      maxlppnl, multiplicity, n_ao, nelectron, nkind, output_unit, natom, iw
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    INTEGER, DIMENSION(2)                    :: n_mo, nelectron_spin, nmoloc
    LOGICAL                                  :: all_potential_present, ionode,&
                                                xas_calculation, was_present
    REAL(dp)                                 :: maxocc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(gapw_control_type), POINTER         :: gapw_control
    TYPE(qs_control_type), POINTER           :: qs_control
    TYPE(localized_wfn_control_type), &
      POINTER                                :: localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(rho0_mpole_type), POINTER           :: nddo_mpole, rho0_mpole
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: xas_section, dft_section
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control

    LOGICAL :: harris_flag
    TYPE(cp_logger_type), POINTER            :: logger
    
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL write_checkpoint_information("entering "//routine_name,globenv,error)

    CALL timeset(routine_name,"I",'',handle)

    ionode = globenv%ionode
    IF(ionode) output_unit = globenv%scr
    was_present = .FALSE.

    IF (PRESENT(harris)) THEN
      harris_flag = harris
    ELSE
      harris_flag = .FALSE.
    END IF

!   *** Initialise the Quickstep environment ***

    NULLIFY (mos)
    NULLIFY (dft_control)
    NULLIFY (energy)
    NULLIFY (force)
    NULLIFY (harris_energy)
    NULLIFY (harris_env)
    NULLIFY (harris_force)
    NULLIFY (local_molecules)
    NULLIFY (local_particles)
    NULLIFY (scf_control)
    NULLIFY (localized_wfn_control)
    NULLIFY (xas_section, dft_section)
    dft_section =>  section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    particle_set => subsystem%particles%els
    atomic_kind_set => subsystem%atomic_kinds%els
    molecule_set => subsystem%molecules_new%els
    molecule_kind_set => subsystem%molecule_kinds_new%els

!   *** Print the Quickstep program banner (copyright and version number) ***
    iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PROGRAM_BANNER",&
         extension=".Log",error=error)    
    CALL qs_header(iw, iw>0)
    CALL cp_print_key_finished_output(iw,logger,dft_section,&
         "PRINT%PROGRAM_BANNER",error=error)    

!   *** Read the input section with the DFT control parameters ***

    CALL read_dft_control(dft_control,globenv)
    ! ONLY WORKS WITH NEW INPUT
    CALL read_tddfpt_control(dft_control%tddfpt_control,dft_control%nspins,globenv)

!   *** Read the input section with the Quickstep control parameters ***

    CALL read_qs_section(dft_control%qs_control,globenv,error=error)
    CALL read_mgrid_section(dft_control%qs_control,globenv,error=error)

!   *** Print the DFT control parameters ***

    CALL write_dft_control(dft_control,dft_section,error)

!   *** Print the Quickstep control parameters ***

    CALL write_qs_control(dft_control%qs_control,dft_section,error)

!   *** Print XES/XAS   control parameters *** ***
    if (dft_control % do_xas_calculation) then
      CALL write_xas_control(dft_control%xas_control,dft_section,error=error)
    end if

!   *** Check basis and fill in missing parts ***

    CALL check_atomic_kind_set(atomic_kind_set,dft_control,globenv,&
         force_env_section=qs_env%input,error=error)

!   *** Check that no all-electron potential is present if GPW or GAPW_XC

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                            all_potential_present=all_potential_present)
    IF(all_potential_present .AND. .NOT.(dft_control%qs_control%method == "GAPW"&
          .OR. dft_control%qs_control%semi_empirical))&
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "all-el calculations with GPW and GAPW_XC are not implemented ")      

!   *** Initialize the spherical harmonics and ***
!   *** the orbital transformation matrices    ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl)

    lmax_sphere = dft_control%qs_control%gapw_control%lmax_sphere
    IF(lmax_sphere .LT.0) THEN
       lmax_sphere = 2*maxlgto
       dft_control%qs_control%gapw_control%lmax_sphere= lmax_sphere
    END IF
    maxl = MAX(2*maxlgto,maxlppl,maxlppnl,lmax_sphere)

    IF (dft_control%forces) maxl = maxl + 1

    CALL init_orbital_pointers(maxl)
    CALL init_spherical_harmonics(maxl,globenv,error)

!   *** Initialise the atomic kind set ***

    CALL init_atomic_kind_set(atomic_kind_set,globenv)

!   *** Print the unnormalized basis set information (input data) ***

    CALL write_gto_basis_sets(atomic_kind_set,qs_env%input,error)

!   *** Initialise GAPW soft basis and projectors

    IF(dft_control%qs_control%method == "GAPW" .OR. &
       dft_control%qs_control%method == "GAPW_XC")  THEN
      qs_control => dft_control%qs_control
      gapw_control => dft_control%qs_control%gapw_control
      CALL init_gapw_basis_set(atomic_kind_set,qs_control,qs_env%input,error)
    ENDIF

!   *** Initialize the pretabulation for the calculation of the   ***
!   *** incomplete Gamma function F_n(t) after McMurchie-Davidson ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto)
    maxl = 3*maxlgto + 1
    CALL init_md_ftable(maxl)

!   *** Print the atomic kind set ***

    CALL write_atomic_kind_set(atomic_kind_set,qs_env%input,error)

!   *** Print the molecule kind set ***

    CALL write_molecule_kind_set(molecule_kind_set,qs_env%input,error)

!   *** Print the total number of kinds, atoms, basis functions etc. ***

    CALL write_total_numbers(atomic_kind_set,particle_set,qs_env%input,error)

!   *** Print the atomic coordinates ***

    CALL write_qs_particle_coordinates(particle_set,cell,qs_env%input,error)

!   *** Print the interatomic distances ***

    CALL write_particle_distances(particle_set,cell,qs_env%input,error)

!   *** Initialize the atomic interaction radii ***
    CALL init_interaction_radii(dft_control%qs_control,cell,atomic_kind_set,&
                                qs_env%input,error)

!   *** Distribute molecules and atoms using the new data structures ***

    CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules,&
                                 force_env_section=qs_env%input,&
                                 error=error)

!   *** SCF parameters ***

    CALL scf_c_create(scf_control)
    CALL scf_c_read_parameters(scf_control,globenv,qs_env%input,error=error)
    CALL scf_c_write_parameters(scf_control,globenv,qs_env%input,error=error)

!   *** Use the just read block sizes to set the optimal values in cp_fm_struct ***

    optimal_blacs_row_block_size=scf_control%nrow_block
    optimal_blacs_col_block_size=scf_control%ncol_block

!   *** Allocate the data structure for Quickstep energies ***

    CALL allocate_qs_energy(energy)

!   *** Allocate the data structure for the Harris energies ***

    IF (harris_flag) THEN
      CALL harris_energy_create(harris_energy=harris_energy,error=error)
    END IF

!   *** Activate the interpolation ***

    CALL wfi_create(wf_history,&
                    interpolation_method_nr=&
                      dft_control%qs_control%wf_interpolation_method_nr,&
                    extrapolation_order = dft_control%qs_control%wf_extrapolation_order,&
                    error=error)

!   *** Set the actual Harris environment ***
    
    IF (harris_flag) THEN
      CALL harris_env_create(harris_env=harris_env, natom=natom, &
                             nspins=dft_control%nspins, error=error)
      CALL harris_env_set(harris_env=harris_env, harris_energy=harris_energy, &
                          harris_force=harris_force, error=error)
    END IF

!   *** Set the current Quickstep environment ***

    CALL set_qs_env(qs_env,subsys=subsystem,error=error)
    CALL set_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    cell_ref=cell_ref,&
                    use_ref_cell=use_ref_cell,&
                    qmmm=qmmm,&
                    dft_control=dft_control,&
                    energy=energy,&
                    force=force,&
                    scf_control=scf_control,&
                    wf_history=wf_history)

    IF (harris_flag) THEN
      CALL set_qs_env(qs_env=qs_env, harris_env=harris_env, error=error)
    END IF

    CALL cp_subsys_set(subsystem,local_molecules_new=local_molecules,&
         local_particles=local_particles,error=error)

    CALL distribution_1d_release(local_particles,error=error)
    CALL distribution_1d_release(local_molecules,error=error)
    CALL scf_c_release(scf_control,error=error)
    CALL wfi_release(wf_history,error=error)
    CALL dft_control_release(dft_control, error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    scf_control=scf_control)

    ! decide what conditions need mo_derivs
    ! right now, this only appears to be OT
    IF (scf_control%use_ot) THEN
      CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.TRUE.)
    ELSE
      CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.FALSE.)
    ENDIF

!   Initialize the GAPW local densities and potentials
    IF(dft_control%qs_control%method_id == do_method_gapw .OR. &
      dft_control%qs_control%method_id == do_method_gapw_xc) THEN
!     *** Allocate and initialize the set of atomic densities ***
      CALL init_rho_atom(qs_env,gapw_control)
      IF(dft_control%qs_control%method_id /= do_method_gapw_xc) THEN
!       *** Allocate and initialize the compensation density rho0 ***
        CALL init_rho0(qs_env,gapw_control,cell,error=error)
!       *** Allocate and Initialize the local coulomb terms (1,2 and 3 centers) ***
        CALL init_coulomb_local(qs_env,gapw_control)
      END IF
    ELSE IF(dft_control%qs_control%semi_empirical) THEN
      se_control => dft_control%qs_control%se_control
      SELECT CASE ( dft_control%qs_control%method )
        CASE DEFAULT
        CASE ("AM1","PM3","MNDO")
          CALL init_rho0_nddo(qs_env, se_control, error)
      END SELECT
    ENDIF

!   *** Allocate the MO data types ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             nsgf=n_ao,&
                             nelectron=nelectron)

    nelectron = nelectron - dft_control%charge ! the total number of electrons

    IF (dft_control%multiplicity == 0) THEN
      IF (MODULO(nelectron,2) == 0) THEN
        dft_control%multiplicity = 1
      ELSE
        dft_control%multiplicity = 2
      END IF
    END IF

    multiplicity = dft_control%multiplicity

    IF ((dft_control%nspins < 1).OR.(dft_control%nspins > 2)) THEN
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "nspins should be 1 or 2 for the time being...")
    END IF

    IF ((MODULO(nelectron,2) /= 0).AND.(dft_control%nspins == 1)) THEN
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "Use the LSD option for an odd number of electrons.")
    END IF

    ! The transition potential method to calculate XAS needs LSD
    xas_section =>  section_vals_get_subs_vals(globenv%input_file, &
                                                 "FORCE_EVAL%DFT%XAS")
    CALL section_vals_val_get(xas_section,"_SECTION_PARAMETERS_",&
         l_val=xas_calculation,error=error)

    IF(xas_calculation .AND. (dft_control%nspins == 1) ) THEN
      CALL stop_program(routine_name,module_name,__LINE__,&
                   "Use the LSD option for XAS with transition potential.")
    END IF
!    IF(xas_calculation .AND.  (scf_control%smear /= 0.0_dp)) THEN
!      CALL stop_program(routine_name,module_name,__LINE__,&
!          "transition potential XAS and smearing not implemented.")
!    END IF

!   assigning the number of states per spin initial version, not yet very
!   general. Should work for an even number of electrons and a single
!   additional electron this set of options that requires full matrices,
!   however, makes things a bit ugly right now.... we try to make a
!   distinction between the number of electrons per spin and the number of
!   MOs per spin this should allow the use of fractional occupations later
!   on

   IF (dft_control%nspins == 1) THEN

     maxocc = 2.0_dp
     nelectron_spin(1) = nelectron
     nelectron_spin(2) = 0
     n_mo(1) = nelectron/2
     n_mo(2) = 0

   ELSE

     maxocc=1.0_dp

!    *** the simplist spin distribution is written here. Special cases will
!    *** need additional user input

     IF (MODULO(nelectron + multiplicity - 1,2) /= 0) THEN
       CALL stop_program(routine_name,module_name,__LINE__,&
                         "LSD: try to use a different multiplicity.")
     END IF

     nelectron_spin(1) = (nelectron + multiplicity - 1)/2
     nelectron_spin(2) = (nelectron - multiplicity + 1)/2

     IF (nelectron_spin(2) < 0) THEN
       CALL stop_program(routine_name,module_name,__LINE__,&
                         "LSD: too few electrons for this multiplicity.")
     END IF

     n_mo(1) = nelectron_spin(1)
     n_mo(2) = nelectron_spin(2)

   END IF

    scf_control%added_mos(1) = MIN(scf_control%added_mos(1),n_ao-n_mo(1))

    n_mo(1) = n_mo(1) + scf_control%added_mos(1)
    IF (n_mo(2) /= 0) n_mo(2) = n_mo(2) + scf_control%added_mos(2)

    IF(scf_control%smear /= 0.0_dp .AND. scf_control%added_mos(1) == 0) THEN
      CALL stop_program(routine_name,module_name,__LINE__,&
                         "More states are needed to apply smeraing.")
    END IF

!   *** Some options require that all MOs are computed ... ***

    IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,&
            "PRINT%MO/CARTESIAN",error=error),cp_p_file).OR.&
        (scf_control%level_shift /= 0.0_dp).OR.(scf_control%eps_jacobi /= 0.0_dp)) THEN
       n_mo(1) = n_ao
       IF (n_mo(2) /= 0) n_mo(2) = n_ao
    END IF

    ALLOCATE (mos(dft_control%nspins),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "mos",dft_control%nspins*int_size)
    END IF

    ! in principle the restricted calculation could be performed 
    ! using just one set of MOs and special casing most of the code
    ! right now we'll just take care of what is effectively an additional constraint
    ! at as few places as possible, just duplicating the beta orbitals
    IF (dft_control%restricted) THEN
       ! it is really not yet tested till the end ! Joost
       IF (globenv%ionode) WRITE(globenv%scr,*) ""
       IF (globenv%ionode) WRITE(globenv%scr,*) " **************************************"
       IF (globenv%ionode) WRITE(globenv%scr,*) " restricted calculation cutting corners"
       IF (globenv%ionode) WRITE(globenv%scr,*) " experimental feature, check code      "
       IF (globenv%ionode) WRITE(globenv%scr,*) " **************************************"
    ENDIF

    DO ispin=1,dft_control%nspins
      NULLIFY (mos(ispin)%mo_set)
      CALL allocate_mo_set(mo_set=mos(ispin)%mo_set,&
                           nao=n_ao,&
                           nmo=n_mo(ispin),&
                           nelectron=nelectron_spin(ispin),&
                           maxocc=maxocc,&
                           error=error)
    END DO

    CALL set_qs_env(qs_env,mos=mos,error=error)

!   *** Some preparetion for the localization of the MOS if required
    IF (dft_control%qs_control%localized_wfn) THEN
      NULLIFY (localized_wfn_control)
      localized_wfn_control =>  dft_control%qs_control%localized_wfn_control
      nmoloc = 0
      IF(localized_wfn_control%set_of_states == state_loc_all) THEN
        DO ispin = 1,dft_control%nspins
          nmoloc(ispin) = NINT(nelectron_spin(ispin)/maxocc)
        END DO
      ELSEIF(localized_wfn_control%set_of_states == state_loc_none) THEN
      ELSE
        DO ispin = 1,dft_control%nspins
          nmoloc(ispin) = MIN(localized_wfn_control%nloc_states(1) ,n_mo(ispin))
        END DO
      END IF
      NULLIFY (localized_wfn_control)
      CALL set_loc_wfn_lists(dft_control%qs_control%localized_wfn_control,&
           nmoloc,n_mo,dft_control%nspins,error=error)
      CALL set_loc_centers(dft_control%qs_control%localized_wfn_control,nmoloc,&
           dft_control%nspins,globenv,error=error)
    END IF 


!   *** Sets up pw_env, qs_charges, mpools... ***

    CALL qs_env_setup(qs_env,error)

!   Allocate and Initialie rho0 soft on the global grid
    IF(dft_control%qs_control%method == "GAPW") THEN 
      CALL get_qs_env(qs_env=qs_env,rho0_mpole=rho0_mpole)
      CALL rho0_s_grid_create(qs_env, rho0_mpole, error=error)
    END IF
!   Allocate grid 
    IF(dft_control%qs_control%se_control%ewald) THEN 
      CALL get_qs_env(qs_env=qs_env,nddo_mpole=nddo_mpole)
      CALL rho0_s_grid_create(qs_env, nddo_mpole, error=error)
    END IF

    IF (ionode) CALL m_flush(output_unit)

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv,error)

  END SUBROUTINE qs_init_subsys

! *****************************************************************************

  SUBROUTINE write_total_numbers(atomic_kind_set,particle_set,force_env_section,error)

!   Purpose: Write the total number of kinds, atoms, etc. to the logical unit
!            number lunit.

!   History: - Creation (06.10.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    
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

       IF ((maxlppl > -1).OR.(maxlppnl > -1)) THEN
          WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T75,I6))")&
               "Maximum angular momentum of the",&
               "- Orbital basis functions:                   ",maxlgto,&
               "- Local part of the GTH pseudopotential:     ",maxlppl,&
               "- Non-local part of the GTH pseudopotential: ",maxlppnl
       ELSE
          WRITE (UNIT=output_unit,FMT="(/,T3,A,T75,I6))")&
               "Maximum angular momentum of the orbital basis functions: ",maxlgto
       END IF

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "PRINT%TOTAL_NUMBERS",error=error)

  END SUBROUTINE write_total_numbers

! *****************************************************************************

END MODULE qs_environment
