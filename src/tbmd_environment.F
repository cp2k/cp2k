!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_environment [1.0] *
!!
!!   NAME
!!     tbmd_environment
!!
!!   FUNCTION
!!     Methods and functions on the tight-binding environment
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_environment
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                             atomic_kind_list_release,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set,&
                                             write_atomic_kind_set
  USE cell_types,                      ONLY: cell_release,&
                                             cell_type,&
                                             get_cell,&
                                             read_cell,&
                                             write_cell
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create
  USE cp_fm_struct,                    ONLY: optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE cp_subsystem_methods,            ONLY: cp_subsys_read_colvar
  USE cp_subsystem_types,              ONLY: cp_subsys_create,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: use_gsp,&
                                             use_kwon,&
                                             use_lenosky,&
                                             use_sawada,&
                                             use_udp
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_create,&
                                             mol_kind_new_list_release,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_create,&
                                             mol_new_list_release,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates,&
                                             write_structure_data
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             mo_set_type
  USE tbmd_energy_types,               ONLY: tbmd_energy_create,&
                                             tbmd_energy_type
  USE tbmd_environment_methods,        ONLY: tbmd_env_setup
  USE tbmd_environment_types,          ONLY: tbmd_env_set,&
                                             tbmd_environment_type
  USE tbmd_force_types,                ONLY: tbmd_force_create,&
                                             tbmd_force_type
  USE tbmd_model_types,                ONLY: tbmd_model_create,&
                                             tbmd_model_type
  USE tbmd_neighbor_list_env_types,    ONLY: tbmd_neighbor_list_env_create,&
                                             tbmd_neighbor_list_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tbmd_environment'

! *** Public subroutines ***

  PUBLIC :: tbmd_init

!!***
! *****************************************************************************

CONTAINS

!!****f* tbmd_environment/tbmd_init [1.0] *
!!
!!   NAME
!!     tbmd_init
!!
!!   FUNCTION
!!     Initialize the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - globenv: The global environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_init(tbmd_env, globenv, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, use_ref_cell
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: cell_section, colvar_section, &
                                                force_env_section, &
                                                tbmd_section

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      ! nullifying  pointers
      NULLIFY(force_env_section, cell_section, colvar_section, tbmd_section, &
              cell, cell_ref, subsystem, molecule_kind_set, molecule_set, &
              particle_set, atomic_kind_set, particles, molecules, &
              molecule_kinds, atomic_kinds)

      force_env_section => section_vals_get_subs_vals(globenv%input_file, &
                                                      "FORCE_EVAL", error=error)
      cell_section => section_vals_get_subs_vals(force_env_section,&
                                                 "SUBSYS%CELL", error=error)
      colvar_section => section_vals_get_subs_vals(force_env_section, &
                                                   "SUBSYS%COLVAR", error=error)
      tbmd_section => section_vals_get_subs_vals(force_env_section, &
                                                 "TBMD", error=error)

      CALL tbmd_env_set(tbmd_env=tbmd_env, tbmd_input=tbmd_section, &
                        force_env_input=force_env_section, para_env=globenv%para_env, &
                        error=error)

      CALL cp_blacs_env_create(blacs_env=tbmd_env%blacs_env, para_env=tbmd_env%para_env, &
                               error=error)

      CALL read_cell(cell=cell, cell_ref=cell_ref, &
                     use_ref_cell=use_ref_cell, cell_section=cell_section, &
                     error=error)
      CALL get_cell(cell=cell, abc=abc)
      CALL write_cell(cell=cell, force_env_section=force_env_section,error=error)

      CALL cp_subsys_create(subsystem, para_env=globenv%para_env, error=error)
      CALL cp_subsys_read_colvar(subsystem, colvar_section, globenv, error=error)

      CALL topology_control(atomic_kind_set, particle_set, &
                            molecule_kind_set, molecule_set, &
                            subsystem%colvar_p,globenv, error=error)

      !CALL read_atomic_kind_set(atomic_kind_set, globenv)

      CALL particle_list_create(particles, els_ptr=particle_set, error=error)
      CALL atomic_kind_list_create(atomic_kinds, els_ptr=atomic_kind_set, error=error)
      CALL mol_new_list_create(molecules, els_ptr=molecule_set, error=error)
      CALL mol_kind_new_list_create(molecule_kinds, els_ptr=molecule_kind_set, error=error)
      CALL cp_subsys_set(subsystem, particles=particles, atomic_kinds=atomic_kinds, &
                         molecules_new=molecules, molecule_kinds_new=molecule_kinds, &
                         error=error)
      CALL particle_list_release(particles, error=error)
      CALL atomic_kind_list_release(atomic_kinds, error=error)
      CALL mol_new_list_release(molecules, error=error)
      CALL mol_kind_new_list_release(molecule_kinds, error=error)

      CALL tbmd_init_subsys(tbmd_env=tbmd_env, subsystem=subsystem, cell=cell, &
                            cell_ref=cell_ref, use_ref_cell=use_ref_cell, &
                            globenv=globenv, error=error)

      CALL cell_release(cell, error=error)
      CALL cell_release(cell_ref, error=error)
      CALL cp_subsys_release(subsystem, error=error)
    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_init
!***************************************************************************

!!****f* tbmd_environment/tbmd_init_subsys [1.0] *
!!
!!   NAME
!!     tbmd_init_subsys
!!
!!   FUNCTION
!!     Initialize the tbmd environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment of matter
!!     - subsystem: the subsystem
!!     - cell: Pointer to the actual simulation cell
!!     - cell_ref: Pointer to the reference cell, used e.g. in NPT simulations
!!     - use_ref_cell: Logical which indicates if cell_ref is in use
!!     - globenv: The global environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_init_subsys(tbmd_env, subsystem, cell, cell_ref, &
                              use_ref_cell, globenv, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL, INTENT(in)                      :: use_ref_cell
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_init_subsys', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: handle, ikind, n_ao, natom, &
                                                nelectron, stat, temp_zeff
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: temp_natom_of_kind
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_ptr
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(mo_set_type), POINTER               :: mos
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    TYPE(tbmd_model_type), POINTER           :: tbmd_model
    TYPE(tbmd_neighbor_list_env_type), &
      POINTER                                :: tbmd_neighbor_list_env

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(atomic_kind_set, molecule_kind_set, molecule_set, particle_set, &
            tbmd_energy, tbmd_force, tbmd_model, local_molecules, &
            local_particles, tbmd_neighbor_list_env, mos, atomic_kind_ptr)

    failure = .FALSE.

    nelectron = 0
    n_ao      = 0
    temp_zeff = 0

    particle_set => subsystem%particles%els
    atomic_kind_set => subsystem%atomic_kinds%els
    molecule_set => subsystem%molecules_new%els
    molecule_kind_set => subsystem%molecule_kinds_new%els

!   *** Print the atomic kind set ***
    CALL write_atomic_kind_set(atomic_kind_set, tbmd_env%force_env_input, &
                               error=error)

!   *** Print the molecule kind set ***
    CALL write_molecule_kind_set(molecule_kind_set, tbmd_env%force_env_input, &
                                 error=error)

!   *** Print the atomic coordinates 
    CALL write_qs_particle_coordinates(particle_set, cell=cell, &
                                       force_env_section=tbmd_env%force_env_input, &
                                       error=error)
    CALL write_particle_distances(particle_set, cell=cell, &
                                  force_env_section=tbmd_env%force_env_input, &
                                  error=error)
    CALL write_structure_data(particle_set, cell=cell, &
                              force_env_section=tbmd_env%force_env_input, &
                              error=error)

!   *** Distribute molecules and atoms using the new data structures ***
    CALL distribute_molecules_1d(particle_kind_set=atomic_kind_set, &
                                 particle_set=particle_set, &
                                 local_particles=local_particles, &
                                 molecule_kind_set=molecule_kind_set, &
                                 molecule_set=molecule_set, &
                                 local_molecules=local_molecules, &
                                 force_env_section=tbmd_env%force_env_input, &
                                 error=error)

    natom = SIZE(particle_set)
    CALL tbmd_energy_create(tbmd_energy=tbmd_energy, error=error)
    CALL tbmd_force_create(tbmd_force=tbmd_force, natom=natom, error=error)
    CALL tbmd_model_create(tbmd_model=tbmd_model, error=error)
    CALL tbmd_neighbor_list_env_create(tbmd_neighbor_list_env=tbmd_neighbor_list_env, &
                                       error=error)

    ALLOCATE(temp_natom_of_kind(SIZE(atomic_kind_set)), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineN, error, failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom_of_kind=temp_natom_of_kind)
    DO ikind = 1,SIZE(atomic_kind_set)
      atomic_kind_ptr => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind_ptr, element_symbol=element_symbol)
      SELECT CASE(element_symbol)
        CASE("SI", "Si", "C")
          temp_zeff = 4
        CASE default
          temp_zeff = 0
      END SELECT
      nelectron = nelectron + temp_natom_of_kind(ikind)*temp_zeff
    END DO

    DEALLOCATE(temp_natom_of_kind, STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineN, error, failure)

    CALL allocate_mo_set(mo_set=mos, nao=4*natom, nmo=nelectron/2, nelectron=nelectron, &
                         maxocc=2.0_dp, error=error)

    CALL tbmd_env_set(tbmd_env=tbmd_env, tbmd_energy=tbmd_energy, &
                      tbmd_force=tbmd_force, tbmd_model=tbmd_model, &
                      tbmd_neighbor_list_env=tbmd_neighbor_list_env, &
                      subsystem=subsystem, local_molecules=local_molecules, &
                      local_particles=local_particles, cell=cell, &
                      cell_ref=cell_ref, use_ref_cell=use_ref_cell, &
                      para_env=globenv%para_env, mos=mos, error=error)

    CALL tbmd_init_tb_model(tbmd_env=tbmd_env, globenv=globenv, error=error)
    CALL tbmd_init_tb_neighbor_list(tbmd_env=tbmd_env, globenv=globenv, &
                                    error=error)

    CALL distribution_1d_release(local_particles, error=error)
    CALL distribution_1d_release(local_molecules, error=error)

    optimal_blacs_row_block_size=32
    optimal_blacs_col_block_size=32

    CALL tbmd_env_setup(tbmd_env=tbmd_env, error=error)

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_init_subsys
!***************************************************************************

!!****f* tbmd_environment/tbmd_init_tb_neighbor_list [1.0] *
!!
!!   NAME
!!     tbmd_init_tb_neighbor_list
!!
!!   FUNCTION
!!     Initialize the tight-binding neighbor-list environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment of matter
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_init_tb_neighbor_list(tbmd_env, globenv, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_init_tb_neighbor_list', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nkinds, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: r_cutoff
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                tbmd_neighbor_list_section, &
                                                tbmd_section
    TYPE(tbmd_neighbor_list_env_type), &
      POINTER                                :: tbmd_neighbor_list_env

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(tbmd_neighbor_list_env, force_env_section, tbmd_section, &
            tbmd_neighbor_list_section, atomic_kind_set)

    tbmd_neighbor_list_env => tbmd_env%tbmd_neighbor_list_env

    force_env_section => section_vals_get_subs_vals(globenv%input_file, &
                                                    "FORCE_EVAL", error=error)
    tbmd_section => section_vals_get_subs_vals(force_env_section, &
                                               "TBMD", error=error)
    tbmd_neighbor_list_section => section_vals_get_subs_vals(tbmd_section, &
                                         "TBMD_NEIGHBOR_LIST", error=error)

    CALL section_vals_val_get(section_vals=tbmd_neighbor_list_section, &
                              keyword_name="VERLET_SKIN", &
                              r_val=tbmd_neighbor_list_env%verlet_skin, &
                              error=error)
    CALL section_vals_val_get(section_vals=tbmd_neighbor_list_section, &
                              keyword_name="R_CUTOFF", &
                              r_val=r_cutoff, error=error)

    atomic_kind_set => tbmd_env%subsystem%atomic_kinds%els
    nkinds = SIZE(atomic_kind_set)
    ALLOCATE(tbmd_neighbor_list_env%r_cutoff(nkinds,nkinds), STAT=stat)
    CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
    tbmd_neighbor_list_env%r_cutoff(:,:) = r_cutoff

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_init_tb_neighbor_list
!***************************************************************************

!!****f* tbmd_environment/tbmd_init_tb_model [1.0] *
!!
!!   NAME
!!     tbmd_init_tb_model
!!
!!   FUNCTION
!!     Initialize the tight-binding model environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment of matter
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_init_tb_model(tbmd_env, globenv, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_init_tb_model', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: tb_atomic_kind_name
    INTEGER                                  :: handle, i
    REAL(KIND=dp)                            :: a_latt, Am, eV
    REAL(KIND=dp), DIMENSION(:), POINTER     :: vector_ptr
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_ptr
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                tbmd_model_section, &
                                                tbmd_section
    TYPE(tbmd_model_type), POINTER           :: tbmd_model

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(tbmd_model, atomic_kind_set, atomic_kind_ptr, force_env_section, &
            tbmd_section, tbmd_model_section, vector_ptr)

    tbmd_model => tbmd_env%tbmd_model

    force_env_section => section_vals_get_subs_vals(globenv%input_file, &
                                                    "FORCE_EVAL", error=error)
    tbmd_section => section_vals_get_subs_vals(force_env_section, &
                                               "TBMD", error=error)
    tbmd_model_section => section_vals_get_subs_vals(tbmd_section, &
                                                     "TB_MODEL", error=error)

    atomic_kind_set => tbmd_env%subsystem%atomic_kinds%els

    eV = 1./27.2116  !1/evolt
    Am = 1./0.529177 !1/angstrom

    ! loop over all kinds
    DO i=1,SIZE(atomic_kind_set)
      atomic_kind_ptr => tbmd_env%subsystem%atomic_kinds%els(i)
      CALL get_atomic_kind(atomic_kind=atomic_kind_ptr, &
                           name=tb_atomic_kind_name)
      SELECT CASE(tb_atomic_kind_name)
        CASE("SI", "Si")
          a_latt = 5.43*Am
          CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                    keyword_name="TB-Model", &
                                    i_val=tbmd_model%tb_model_name, &
                                    error=error)
          SELECT CASE(tbmd_model%tb_model_name)
            CASE(use_gsp)
              CALL section_vals_val_get(section_vals=tbmd_section, &
                                        keyword_name="Noisy_TB", &
                                        l_val=tbmd_model%noisy_tb, &
                                        error=error)

              tbmd_model%h_ssS         = -1.82_dp*eV
              tbmd_model%h_spS         =  1.96_dp*eV
              tbmd_model%h_ppS         =  3.06_dp*eV
              tbmd_model%h_ppP         = -0.87_dp*eV
              tbmd_model%n             =  2.0_dp
              tbmd_model%r_0           =  (SQRT(3.0_dp)/4.0_dp)*a_latt
              tbmd_model%Phi           =  3.4581_dp*eV
              tbmd_model%m             =  4.54_dp
              tbmd_model%d_c           =  3.67_dp*Am
              tbmd_model%m_c           =  6.48_dp
              tbmd_model%dEpsilon      =  8.295_dp*eV
              tbmd_model%E_s           =  0.0_dp
              tbmd_model%E_0           =  0.0_dp
              tbmd_model%r_calpha(:)   =  3.67_dp*Am
              tbmd_model%n_calpha(:)   =  6.48_dp
              tbmd_model%C14(2:4)      =  0.0_dp
              tbmd_model%C14(1)        =  1.0_dp
            CASE(use_sawada)
              CALL cp_unimplemented_error(fromWhere=routineP, &
                   message="The Sawada tight-binding" //&
                   "parametrization isn't implemented yet.", &
                   error=error, error_level=cp_failure_level)

              tbmd_model%tb_model_name = use_sawada
              CALL section_vals_val_get(section_vals=tbmd_section, &
                                        keyword_name="Noisy_TB", &
                                        l_val=tbmd_model%noisy_tb, &
                                        error=error)

            CASE(use_kwon)
              CALL cp_unimplemented_error(fromWhere=routineP, &
                   message="The Kwon tight-binding" //&
                   "parametrization isn't implemented yet.", &
                   error=error, error_level=cp_failure_level)

              CALL section_vals_val_get(section_vals=tbmd_section, &
                                        keyword_name="Noisy_TB", &
                                        l_val=tbmd_model%noisy_tb, &
                                        error=error)

            CASE(use_lenosky)
              CALL cp_unimplemented_error(fromWhere=routineP, &
                   message="The Lenosky tight-binding" //&
                   "parametrization isn't implemented yet.", &
                   error=error, error_level=cp_failure_level)

              CALL section_vals_val_get(section_vals=tbmd_section, &
                                        keyword_name="Noisy_TB", &
                                        l_val=tbmd_model%noisy_tb, &
                                        error=error)

            CASE(use_udp)
              CALL section_vals_val_get(section_vals=tbmd_section, &
                                        keyword_name="Noisy_TB", &
                                        l_val=tbmd_model%noisy_tb, &
                                        error=error)

              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="h_ssS", &
                                        r_val=tbmd_model%h_ssS, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="h_spS", &
                                        r_val=tbmd_model%h_spS, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="h_ppS", &
                                        r_val=tbmd_model%h_ppS, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="h_ppP", &
                                        r_val=tbmd_model%h_ppP, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="n", &
                                        r_val=tbmd_model%n, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="r_0", &
                                        r_val=tbmd_model%r_0, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="Phi(i)", &
                                        r_val=tbmd_model%Phi, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="m", &
                                        r_val=tbmd_model%m, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="d_c", &
                                        r_val=tbmd_model%d_c, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="m_c", &
                                        r_val=tbmd_model%m_c, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="dEpsilon", &
                                        r_val=tbmd_model%dEpsilon, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="E_s", &
                                        r_val=tbmd_model%E_s, &
                                        error=error)
              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="E_0", &
                                        r_val=tbmd_model%E_0, &
                                        error=error)

              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="r_calpha", &
                                        r_vals=vector_ptr, &
                                        error=error)
              tbmd_model%r_calpha(1:4) = vector_ptr(1:4)

              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="n_calpha", &
                                        r_vals=vector_ptr, &
                                        error=error)
              tbmd_model%n_calpha(1:4) = vector_ptr(1:4)

              CALL section_vals_val_get(section_vals=tbmd_model_section, &
                                        keyword_name="C14", &
                                        r_vals=vector_ptr, &
                                        error=error)
              tbmd_model%C14(1:4) = vector_ptr(1:4)
            CASE DEFAULT
              CALL cp_unimplemented_error(fromWhere=routineP, &
                   message="Other tight-binding parametrizations"//&
                   "for Silicon than GSP, Sawada, Kwon & Lenosky"//&
                   "aren't implemented yet.", &
                    error=error, error_level=cp_failure_level)
          END SELECT
        CASE DEFAULT
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Tight-binding for other elements" //&
               "than Si isn't implemented yet.", &
               error=error, error_level=cp_failure_level)
      END SELECT
    END DO

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_init_tb_model
!***************************************************************************

END MODULE tbmd_environment
