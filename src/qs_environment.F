!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set,&
                                             init_atomic_kind_set,&
                                             read_atomic_kind_set,&
                                             write_atomic_kind_set,&
                                             write_gto_basis_sets
  USE checkpoint_handler,              ONLY: write_checkpoint_information
  USE cntl_input,                      ONLY: read_cntl_section
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_fm_struct,                    ONLY: optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE dft_types,                       ONLY: dft_control_release,&
                                             dft_control_type,&
                                             read_dft_control,&
                                             read_qs_control,&
                                             write_dft_control,&
                                             write_qs_control
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE gamma,                           ONLY: init_md_ftable
  USE global_types,                    ONLY: global_environment_type
  USE input_types,                     ONLY: setup_parameters_type
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE machine,                         ONLY: m_flush
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: molecule_type
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE pao_read_section,                ONLY: pao_env_init_from_file
  USE pao_types,                       ONLY: pao_env_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             read_particle_set,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env,&
                                             qs_env_setup
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             init_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type
  USE qs_interactions,                 ONLY: init_interaction_radii
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_wf_history_methods,           ONLY: wfi_create
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_release
  USE scf_control_types,               ONLY: scf_c_create,&
                                             scf_c_read_parameters,&
                                             scf_c_release,&
                                             scf_c_write_parameters,&
                                             scf_control_type
  USE setup_input,                     ONLY: read_setup_section
  USE simulation_cell,                 ONLY: cell_release,&
                                             cell_type,&
                                             read_cell,&
                                             write_cell
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_environment"

! *** Public subroutines ***

  PUBLIC :: qs_init

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE qs_init(qs_env,globenv,error,new_h_matrix)

!   Purpose: Read the input and the database files for the setup of the
!            QUICKSTEP environment.

!   History: - Creation (22.05.2000,MK)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER           :: qs_env
    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    REAL(dp),DIMENSION(1:3,1:3),OPTIONAL,INTENT(IN) :: new_h_matrix

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "qs_init"

!   *** Local variables ***

    TYPE(cell_type), POINTER            :: cell
    TYPE(dft_control_type), POINTER     :: dft_control
    TYPE(distribution_1d_type), POINTER :: local_molecules,&
                                           local_particles
    TYPE(qs_energy_type), POINTER       :: energy
    TYPE(qs_wf_history_type), POINTER   :: wf_history
    TYPE(pao_env_type), POINTER         :: pao_env
    TYPE(scf_control_type), POINTER     :: scf_control

    TYPE(setup_parameters_type) :: setup

    REAL(dp) :: maxocc
    INTEGER  :: handle,ikind,ispin,istat,maxl,maxlgto,maxlppl,maxlppnl,&
                multiplicity,n_ao,natom,nelectron,nkind,output_unit
    LOGICAL  :: ionode

    TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
    TYPE(mo_set_p_type), DIMENSION(:), POINTER      :: c
    TYPE(molecule_kind_type), DIMENSION(:), POINTER :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), POINTER      :: molecule_set
    TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force

    INTEGER, DIMENSION(2) :: n_mo,nelectron_spin

    INTEGER, DIMENSION(:), ALLOCATABLE :: natom_of_kind

!   ---------------------------------------------------------------------------

    CALL write_checkpoint_information("entering "//routine_name,globenv)

    CALL timeset(routine_name,"I",'',handle)


    ionode = globenv%ionode
    output_unit = globenv%scr

!   *** Initialise the Quickstep environment ***

    NULLIFY (atomic_kind_set)
    NULLIFY (c)
    NULLIFY (cell)
    NULLIFY (dft_control)
    NULLIFY (energy)
    NULLIFY (force)
    NULLIFY (local_molecules)
    NULLIFY (local_particles)
    NULLIFY (molecule_set)
    NULLIFY (pao_env)
    NULLIFY (particle_set)
    NULLIFY (scf_control)

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

    IF (PRESENT(new_h_matrix)) THEN
        CALL read_cell(cell,globenv,new_h_matrix)
    ELSE
        CALL read_cell(cell,globenv)
    ENDIF

!   *** Print the cell parameters ***

    CALL write_cell(cell,globenv)

    CALL read_cntl_section(setup,globenv)

    CALL read_setup_section(setup,globenv)

!   *** Read the particle coordinates and allocate the atomic kind, ***
!   *** the molecule kind, and the molecule data structures         ***
! topology_control
    CALL topology_control (atomic_kind_set,particle_set,&
                         molecule_kind_set,molecule_set,&
                         globenv)
  
    CALL read_atomic_kind_set(atomic_kind_set,globenv)

!   *** Print the unnormalized basis set information (input data) ***

    CALL write_gto_basis_sets(atomic_kind_set,globenv)

!   *** Initialize the spherical harmonics and ***
!   *** the orbital transformation matrices    ***

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             maxlppl=maxlppl,&
                             maxlppnl=maxlppnl)

    maxl = MAX(maxlgto,maxlppl,maxlppnl)

    IF (dft_control%forces) maxl = maxl + 1

    CALL init_orbital_pointers(maxl)
    CALL init_spherical_harmonics(maxl,globenv)

!   *** Initialize the pretabulation for the calculation of the   ***
!   *** incomplete Gamma function F_n(t) after McMurchie-Davidson ***

    maxl = 3*maxlgto + 1

    CALL init_md_ftable(maxl)

!   *** Initialise the atomic kind set ***

    CALL init_atomic_kind_set(atomic_kind_set,globenv)

!   *** Print the atomic kind set ***

    CALL write_atomic_kind_set(atomic_kind_set,globenv)

!   *** Print the molecule kind set ***

    CALL write_molecule_kind_set(molecule_kind_set,globenv)

!   *** Print the total number of kinds, atoms, basis functions etc. ***

    CALL write_total_numbers(atomic_kind_set,particle_set,globenv)

!   *** Print the atomic coordinates ***

    IF (globenv%print%atomic_coordinates) THEN
      CALL write_qs_particle_coordinates(particle_set,cell,globenv)
    END IF

!   *** Print the interatomic distances ***

    CALL write_particle_distances(particle_set,cell,globenv)

!   *** Initialize the atomic interaction radii ***

    CALL init_interaction_radii(dft_control%qs_control,cell,atomic_kind_set,&
                                globenv)

!   *** Distribute molecules and atoms using the new data structures ***

    CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set,&
                                 particle_set=particle_set,&
                                 local_particles=local_particles,&
                                 molecule_kind_set=molecule_kind_set,&
                                 molecule_set=molecule_set,&
                                 local_molecules=local_molecules,&
                                 error=error)

!   *** SCF parameters ***

    CALL scf_c_create(scf_control)
    CALL scf_c_read_parameters(scf_control,globenv)
    CALL scf_c_write_parameters(scf_control,globenv)

!   *** Use the just read block sizes to set the optimal values in cp_fm_struct ***

    optimal_blacs_row_block_size=scf_control%nrow_block
    optimal_blacs_col_block_size=scf_control%ncol_block

!   *** Allocate the data structure for Quickstep energies ***

    CALL allocate_qs_energy(energy)

!   *** Allocate the force data structure ***

    IF (dft_control%forces) THEN
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "natom_of_kind",nkind*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom_of_kind=natom_of_kind)
      CALL allocate_qs_force(force,natom_of_kind)
      DEALLOCATE (natom_of_kind,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "natom_of_kind")
    END IF

!   *** Activate the interpolation ***

    CALL wfi_create(wf_history,&
                    interpolation_method_nr=&
                      dft_control%qs_control%wf_interpolation_method_nr,&
                    extrapolation_order = dft_control%qs_control%wf_extrapolation_order,&
                    error=error)

!   *** Set the current Quickstep environment ***

    CALL set_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    dft_control=dft_control,&
                    energy=energy,&
                    force=force,&
                    local_molecules=local_molecules,&
                    local_particles=local_particles,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    pao_env=pao_env,&
                    particle_set=particle_set,&
                    scf_control=scf_control,&
                    wf_history=wf_history)

    CALL cell_release(cell,error=error)
    CALL distribution_1d_release(local_particles,error=error)
    CALL distribution_1d_release(local_molecules,error=error)
    CALL scf_c_release(scf_control,error=error)
    CALL wfi_release(wf_history,error=error)
    CALL dft_control_release(dft_control, error=error)

!   *** Allocate the MO data types ***

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    scf_control=scf_control)

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

     IF (nelectron_spin(2) <= 0) THEN
       CALL stop_program(routine_name,module_name,__LINE__,&
                         "LSD: too few electrons for this multiplicity.")
     END IF

     n_mo(1) = nelectron_spin(1)
     n_mo(2) = nelectron_spin(2)

    END IF

!   *** Some options require that all MOs are computed ... ***

    IF (globenv%print%mo_eigenvalues.OR.&
        (scf_control%level_shift /= 0.0_dp).OR.&
        (scf_control%smear /= 0.0_dp).OR.&
        (scf_control%eps_jacobi /= 0.0_dp)) THEN
       n_mo(1) = n_ao
       IF (n_mo(2) /= 0) n_mo(2) = n_ao
    END IF

    ALLOCATE (c(dft_control%nspins),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "c",dft_control%nspins*int_size)
    END IF

    DO ispin=1,dft_control%nspins
      NULLIFY (c(ispin)%mo_set)
      CALL allocate_mo_set(mo_set=c(ispin)%mo_set,&
                           nao=n_ao,&
                           nmo=n_mo(ispin),&
                           nelectron=nelectron_spin(ispin),&
                           maxocc=maxocc,&
                           error=error)
    END DO

    CALL set_qs_env(qs_env,c=c,error=error)

!   *** Sets up pw_env, qs_charges, mpools... ***

    CALL qs_env_setup(qs_env,error)

!   *** Read the pao section ***

    IF (dft_control%qs_control%pao) THEN
      ALLOCATE (pao_env) ! check status ?
      CALL set_qs_env(qs_env=qs_env,pao_env=pao_env)
      CALL pao_env_init_from_file(qs_env,globenv)
!     CALL cp_write_description(pao_env, unit_nr=globenv%scr)
    END IF

    IF (ionode) CALL m_flush(output_unit)

    CALL timestop(0.0_dp,handle)

    CALL write_checkpoint_information("leaving "//routine_name,globenv)

  END SUBROUTINE qs_init

! *****************************************************************************

  SUBROUTINE write_qs_program_banner(globenv)

!   Purpose: Write the Quickstep program banner with the copyright and the
!            program version information to the output unit.
!   JVdV it is also a great place to provide references to QS related papers

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
        "**            Copyright (C) by CP2K Developers Group (2000 - 2003)           **",&
        "**                                                                           **",&
        "*******************************************************************************",&
        "** G. Lippert, J. Hutter, and M. Parrinello,                                 **",&
        "** A hybrid Gaussian and plane wave density functional scheme,               **",&
        "** Molecular Physics, Vol. 92, Issue 3, 477-487 (1997)                       **",&
        ! GAPW references
        ! "*******************************************************************************",&
        ! "** G. Lippert, J. Hutter, and M. Parrinello                                  **",&
        ! "** The Gaussian and augmented-plane-wave density functional method for       **",&
        ! "** ab initio molecular dynamics simulations                                  **",&
        ! "** Theoretical Chemistry Accounts, Vol. 103, Issue 2, 124-140 (1999)         **",&
        ! "*******************************************************************************",&
        ! "** M. Krack and M. Parrinello                                                **",&
        ! "** All-electron ab-initio molecular dynamics                                 **",&
        ! "** Physical Chemistry Chemical Physics, Vol. 2, Issue 10, 2105-2112 (2000)   **",&
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

    INTEGER :: iatom,maxlgto,maxlppl,maxlppnl,natom,ncgf,nkind,npgf,nset,&
               nsgf,nshell,output_unit

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
          "- Orbital basis functions:                   ",maxlgto,&
          "- Local part of the GTH pseudopotential:     ",maxlppl,&
          "- Non-local part of the GTH pseudopotential: ",maxlppnl
      ELSE
        WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,(T3,A,I6))")&
          "Maximum angular momentum quantum number of the",&
          "- Orbital basis functions: ",maxlgto
      END IF

    END IF

  END SUBROUTINE write_total_numbers

! *****************************************************************************

!*************************************************************
END MODULE qs_environment
