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
                                             attach_atomic_kind_set,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             set_atomic_kind,&
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
  USE cp_para_env,                     ONLY: cp_para_env_create
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d,&
                                             distribute_molecules_2d
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
  USE kg_force_types,                  ONLY: allocate_kg_force,&
                                             kg_force_type
  USE kg_rho_methods,                  ONLY: kg_rho_create
  USE kg_rho_types,                    ONLY: kg_rho_release,&
                                             kg_rho_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_comm_split
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
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_environment_p_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type
  USE qs_interactions,                 ONLY: init_interaction_radii
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
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
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

    ! Read the input section with the dft control parameters *
    CALL read_dft_control ( dft_control, globenv, error )
    CALL write_dft_control ( dft_control, globenv )

    ! Read the input section with the KG control parameters 
    CALL read_qs_control ( dft_control % qs_control, globenv )
    CALL write_qs_control ( dft_control % qs_control, globenv )

    ! Read/write the input section with the cell parameters
    CALL read_cell ( cell, globenv )
    CALL write_cell ( cell, globenv )

    ! topology
    CALL topology_control (atomic_kind_set,particle_set,&
                           molecule_kind_set,molecule_set,&
                           globenv)

    ! Read/write the atomic kind set
    CALL read_atomic_kind_set(atomic_kind_set,globenv)
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

  SUBROUTINE kg_env_setup(kg_env,error)

    TYPE(kg_environment_type), INTENT(inout) :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_env_setup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(kg_rho_type), POINTER               :: rho
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure=.FALSE.

    NULLIFY(rho)

    CALL get_kg_env(kg_env=kg_env,&
                    atomic_kind_set=atomic_kind_set,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    particle_set=particle_set,&
                    rho=rho)
    CALL kg_rho_create(kg_rho=rho,kg_env=kg_env,error=error)
    CALL set_kg_env(kg_env,rho=rho)
    CALL kg_rho_release(kg_rho=rho,error=error)
     ! create 2d distribution

    CALL distribute_molecules_2d(particle_kind_set=atomic_kind_set,&
                                  particle_set=particle_set,&
                                  molecule_kind_set=molecule_kind_set,&
                                  molecule_set=molecule_set,&
                                  distribution_2d=distribution_2d,&
                                  blacs_env=kg_env%blacs_env,&
                                  error=error)

    CALL set_kg_env(kg_env,distribution_2d=distribution_2d,error=error)
    CALL distribution_2d_release(distribution_2d, error=error)

  END SUBROUTINE kg_env_setup

!******************************************************************************

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

      WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(T15,A,T75,I6))")&
        "Total number of",&
        "- Atomic kinds:                  ",nkind,&
        "- Atoms:                         ",natom,&
        "- Shell sets:                    ",nset,&
        "- Shells:                        ",nshell,&
        "- Primitive Cartesian functions: ",npgf,&
        "- Cartesian basis functions:     ",ncgf,&
        "- Spherical basis functions:     ",nsgf

      IF ( maxlppnl > -1 ) THEN
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(T15,A,T75,I6))")&
          "Maximum angular momentum quantum number of the",&
          "- Orbital basis functions:                   ",maxlgto,&
          "- Local part of the GTH pseudopotential:     ",maxlppl,&
          "- Non-local part of the GTH pseudopotential: ",maxlppnl
      ELSE IF ( maxlppl > -1 ) THEN
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(T15,A,T75,I6))")&
          "Maximum angular momentum quantum number of the",&
          "- Orbital basis functions:                   ",maxlgto,&
          "- Local part of the GTH pseudopotential:     ",maxlppl
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(T15,A,T75,I6))")&
          "Maximum angular momentum quantum number of the",&
          "- Orbital basis functions: ",maxlgto
      END IF

    END IF

  END SUBROUTINE write_total_numbers

!******************************************************************************

  SUBROUTINE kg_subsys_setup(kg_env,error)

!   We generate a qs_environment for each subsystem
!   We make the assignment molecule==subsystem

    TYPE(kg_environment_type), INTENT(inout) :: kg_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kg_subsys_setup', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, i, iatom, ierror, ikind, imol, iset, j, jmol, &
      last_atom, n, natom, ngroups, nkind, nmol, nmol_tot, nset, sub_comm
    INTEGER, DIMENSION(:), POINTER           :: atom_list, &
                                                group_distribution, &
                                                molecule_list, natom_of_kind
    LOGICAL                                  :: failure
    TYPE(atom_type), DIMENSION(:), POINTER   :: mol_atom_list
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_p_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set_list
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: sub_atomic_kind
    TYPE(cell_type), POINTER                 :: cell, sub_cell
    TYPE(cp_para_env_type), POINTER          :: para_env, sub_para_env
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(mol_kind_new_list_type), POINTER    :: mol_kinds
    TYPE(mol_new_list_type), POINTER         :: mols
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set, &
                                                sub_molecule_kind_set
    TYPE(molecule_kind_type), POINTER        :: molecule_kind, &
                                                sub_molecule_kind
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set, sub_molecule_set
    TYPE(molecule_type), POINTER             :: sub_molecule
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set, sub_particle_set
    TYPE(qs_energy_type), POINTER            :: sub_energy
    TYPE(qs_environment_p_type), &
      DIMENSION(:), POINTER                  :: qs_env
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: sub_force

!   ---------------------------------------------------------------------------

    failure=.FALSE.

    CALL get_kg_env(kg_env=kg_env,&
                    dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    para_env=para_env,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    particle_set=particle_set   ,&
                    local_molecules=local_molecules)

    ! Each processor works alone on subsystems
    ! This assumes #subsystems >> #cpus
    ! Loadbalancing is given by the local_molecules distribution
    ! A more general code would take all possible cases into account
    ! and make use of better estimates for work needed for different 
    ! subsystems
    ALLOCATE ( group_distribution(0:para_env%num_pe-1), STAT=ierror )
    IF(ierror/=0) CALL stop_memory('kg_subsys_setup',&
       'group_distribution',para_env%num_pe)
    ! generate the new group (just the local proc for the moment)
    CALL mp_comm_split(para_env%group,sub_comm,ngroups,group_distribution, &
                       n_subgroups=para_env%num_pe)
    CALL cp_para_env_create(sub_para_env,sub_comm,error=error)
    DEALLOCATE ( group_distribution, STAT=ierror )
    IF(ierror/=0) CALL stop_memory('kg_subsys_setup','group_distribution')

    ! allocate and generate a qs_env for each subsystem

    nkind = SIZE (local_molecules%n_el)
    nmol_tot = SUM ( local_molecules%n_el(1:nkind) )
    ALLOCATE ( qs_env(nmol_tot), STAT=ierror )
    IF(ierror/=0) CALL stop_memory('kg_subsys_setup','qs_env',nmol_tot)
    DO imol=1,nmol_tot
      CALL qs_env_create(qs_env(imol)%qs_env, para_env=sub_para_env, error=error)
    END DO

    ALLOCATE ( atomic_kind_set_list(nmol_tot), STAT=ierror )
    IF(ierror/=0) CALL stop_memory('kg_subsys_setup',&
                                   'atomic_kind_set_list',nmol_tot)
    nset = SIZE (atomic_kind_set)
    DO imol=1,nmol_tot
      CALL allocate_atomic_kind_set(atomic_kind_set_list(imol)%atomic_kind_set,nset)
      CALL attach_atomic_kind_set(atomic_kind_set_list(imol)%atomic_kind_set,&
                                  atomic_kind_set)
    END DO

    ALLOCATE ( molecule_list(1), STAT=ierror )
    IF(ierror/=0) CALL stop_memory('kg_subsys_setup','molecule_list',1)
    molecule_list = 1

    jmol = 0
    DO ikind=1,nkind
      nmol = local_molecules%n_el(ikind)
      DO imol=1,nmol
        jmol = jmol + 1
        i = local_molecules%list(ikind)%array(imol)
        molecule_kind => molecule_set(i)%molecule_kind
        CALL get_molecule_kind(molecule_kind=molecule_kind,&
                               atom_list=mol_atom_list,&
                               natom=natom)
        first_atom=molecule_set(i)%first_atom
        last_atom=first_atom+natom-1

        NULLIFY(sub_molecule_kind_set)
        CALL allocate_molecule_kind_set(sub_molecule_kind_set,1)
        sub_molecule_kind_set(1) = molecule_kind
        sub_molecule_kind => sub_molecule_kind_set(1)
        CALL set_molecule_kind(molecule_kind=sub_molecule_kind,&
                               molecule_list=molecule_list)
        NULLIFY(sub_molecule_set)
        CALL allocate_molecule_set(sub_molecule_set,1)
        CALL set_molecule_set(sub_molecule_set,molecule_list)
        sub_molecule => sub_molecule_set(1)
        CALL set_molecule(sub_molecule,sub_molecule_kind,&
                 molecule_set(i)%lmi)

        NULLIFY(sub_particle_set)
        CALL allocate_particle_set(sub_particle_set,natom)
        sub_particle_set(:) = particle_set(first_atom:last_atom)

        ALLOCATE ( atom_list(natom), STAT=ierror )
        IF(ierror/=0) CALL stop_memory('kg_subsys_setup','atom_list',natom)
        ALLOCATE ( natom_of_kind(nset), STAT=ierror )
        IF(ierror/=0) CALL stop_memory('kg_subsys_setup','natom_of_kind',nset)
        natom_of_kind = 0

        DO iset = 1, nset
           atom_list = 0
           n = 0
           DO iatom = 1,natom
             CALL get_atomic_kind(mol_atom_list(iatom)%atomic_kind,kind_number=j)
             IF ( j == iset ) THEN
               n = n + 1
               atom_list(n) = iatom
             END IF
           END DO
           natom_of_kind(iset) = n
           DO iatom = 1,n
             j = atom_list(iatom)
             sub_particle_set(j)%atomic_kind => &
                              atomic_kind_set_list(jmol)%atomic_kind_set(iset)
           END DO
           sub_atomic_kind => atomic_kind_set_list(jmol)%atomic_kind_set(iset)
           CALL set_atomic_kind(sub_atomic_kind,atom_list=atom_list(1:n))
        END DO

        NULLIFY(sub_energy)
        CALL allocate_qs_energy(sub_energy)
        NULLIFY(sub_force)
        CALL allocate_qs_force(sub_force,natom_of_kind)

        CALL distribute_molecules_2d(particle_kind_set=&
                                     atomic_kind_set_list(jmol)%atomic_kind_set,&
                                     particle_set=sub_particle_set,&
                                     molecule_kind_set=sub_molecule_kind_set,&
                                     molecule_set=sub_molecule_set,&
                                     distribution_2d=distribution_2d,&
                                     blacs_env=qs_env(jmol)%qs_env%blacs_env,&
                                     error=error)

        !initialize the sub_cell
        sub_cell => cell
        CALL cell_retain(cell,error)

        !FM initializes the subsys
        NULLIFY(subsystem,particles,atomic_kinds,mols,mol_kinds)
        CALL cp_subsys_create(subsystem,para_env=para_env,error=error)
        CALL particle_list_create(particles,els_ptr=sub_particle_set,error=error)
        CALL atomic_kind_list_create(atomic_kinds,&
             els_ptr=atomic_kind_set_list(jmol)%atomic_kind_set,&
             owns_els=.TRUE.,& !FM correct? should that be false?
             error=error)
        CALL mol_new_list_create(mols,els_ptr=sub_molecule_set,error=error)
        CALL mol_kind_new_list_create(mol_kinds,els_ptr=sub_molecule_kind_set,error=error)
        CALL cp_subsys_set(subsystem,particles=particles,&
             atomic_kinds=atomic_kinds,&
             molecules_new=mols,molecule_kinds_new=mol_kinds,error=error)
        CALL particle_list_release(particles,error=error)
        CALL atomic_kind_list_release(atomic_kinds,error=error)
        CALL mol_new_list_release(mols,error=error)
        CALL mol_kind_new_list_release(mol_kinds,error=error)
       !FM you might want to have a list of subsys, and just one qs_env per proc.
       !FM where you swap the subsys, in that case qs_environment:qs_init_subsys
       !FM would be a starting point, but you should further split up
       !FM the initializazion that don't need to be done a second time.

        CALL set_qs_env(qs_env=qs_env(jmol)%qs_env,&
                        cell=sub_cell             ,&
                        dft_control=dft_control   ,&
!d                      c=                        ,&
                        energy=sub_energy         ,&
                        force=sub_force           ,&
!d                      h=                        ,&
!d                      k=                        ,&
!d                      s=                        ,&
!d                      w=                        ,&
!d                      pw_env=                   ,&
!d                      mpools=                   ,&
!d                      rho=                      ,&
!d                      kinetic=                  ,&
!d                      scf_control=              ,&
!d                      qs_charges=               ,&
!d                      rho_core=                 ,&
!d                      ks_env=                   ,&
!d                      wf_history=               ,&
!d                      scf_env=                  ,&
!d                      local_particles=          ,&
!d                      local_molecules=          ,&
                        distribution_2d=distribution_2d,&
                        subsys=subsystem          ,&
                        error=error                 )

        CALL cp_subsys_release(subsystem,error=error)
        CALL distribution_2d_release(distribution_2d, error=error)

        DEALLOCATE ( atom_list, STAT=ierror )
        IF(ierror/=0) CALL stop_memory('kg_subsys_setup','atom_list')
        DEALLOCATE ( natom_of_kind, STAT=ierror )
        IF(ierror/=0) CALL stop_memory('kg_subsys_setup','natom_of_kind')

      END DO
    END DO

    CALL set_kg_env(kg_env=kg_env,&
                    sub_qs_env=qs_env)

    DEALLOCATE ( atomic_kind_set_list, STAT=ierror )
    IF(ierror/=0) CALL stop_memory('kg_subsys_setup','atomic_kind_set_list')
    DEALLOCATE ( molecule_list, STAT=ierror )
    IF(ierror/=0) CALL stop_memory('kg_subsys_setup','molecule_list')

  END SUBROUTINE kg_subsys_setup

!******************************************************************************

END MODULE kg_environment

!******************************************************************************

