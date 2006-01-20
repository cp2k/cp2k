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
  USE atomic_kind_list_types,                 ONLY: atomic_kind_list_create,&
                                                    atomic_kind_list_release,&
                                                    atomic_kind_list_type
  USE atomic_kind_types,                      ONLY: atomic_kind_type, &
                                                    read_atomic_kind_set, &
                                                    get_atomic_kind
  USE cp_subsystem_methods,                   ONLY: cp_subsys_read_colvar
  USE cp_subsystem_types,                     ONLY: cp_subsystem_type, &
                                                    cp_subsys_create, &
                                                    cp_subsys_release, &
                                                    cp_subsys_set
  USE global_types,                           ONLY: global_environment_type
  USE input_constants,                        ONLY: use_gsp, &
                                                    use_sawada, &
                                                    use_kwon, &
                                                    use_lenosky, &
                                                    use_udp
  USE input_section_types,                    ONLY: section_vals_get_subs_vals,&
                                                    section_vals_type
  USE kinds,                                  ONLY: dp, &
                                                    default_string_length
  USE mol_kind_new_list_types,                ONLY: mol_kind_new_list_create,&
                                                    mol_kind_new_list_release,&
                                                    mol_kind_new_list_type
  USE mol_new_list_types,                     ONLY: mol_new_list_create,&
                                                    mol_new_list_release,&
                                                    mol_new_list_type
  USE molecule_kind_types,                    ONLY: molecule_kind_type
  USE molecule_types_new,                     ONLY: molecule_type
  USE particle_list_types,                    ONLY: particle_list_create, &
                                                    particle_list_type, &
                                                    particle_list_release
  USE particle_types,                         ONLY: particle_type
  USE physcon,                                ONLY: evolt, &
                                                    angstrom
  USE cell_types,                        ONLY: cell_type, &
                                                    read_cell, &
                                                    get_cell, &
                                                    cell_release,&
                                                    write_cell
  USE tbmd_energy_types,                      ONLY: tbmd_energy_type, &
                                                    tbmd_energy_create
  USE tbmd_environment_types,                 ONLY: tbmd_environment_type, &
                                                    tbmd_env_set
  USE tbmd_force_types,                       ONLY: tbmd_force_type, &
                                                    tbmd_force_create
  USE tbmd_model_types,                       ONLY: tbmd_model_type, &
                                                    tbmd_model_create
  USE topology,                               ONLY: topology_control

#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: moduleN = "tbmd_environment"

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

    LOGICAL                                  :: failure, &
                                                use_ref_cell
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                cell_section, &
                                                colvar_section, &
                                                tbmd_section
    TYPE(cell_type), POINTER                 :: cell, &
                                                cell_ref
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(mol_new_list_type), POINTER         :: molecules
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
 
!   ------------------------------------------------------------------------

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      ! nullifying  pointers
      NULLIFY(force_env_section, cell_section, colvar_section, tbmd_section, &
              cell, cell_ref, subsystem, molecule_kind_set, molecule_set, &
              particle_set, atomic_kind_set, particles, molecules, molecule_kinds)

      force_env_section => section_vals_get_subs_vals(globenv%input_file, &
                                                      "FORCE_EVAL", error=error)
      cell_section => section_vals_get_subs_vals(force_env_section,&
                                                 "SUBSYS%CELL", error=error)
      colvar_section => section_vals_get_subs_vals(force_env_section, &
                                                   "SUBSYS%COLVAR", error=error)
      tbmd_section => section_vals_get_subs_vals(force_env_section, &
                                                 "TBMD", error=error)

      CALL tbmd_env_set(tbmd_env=tbmd_env, tbmd_input=tbmd_section, error=error)

      CALL read_cell(cell=cell, cell_ref=cell_ref, globenv=globenv, &
                     use_ref_cell=use_ref_cell, cell_section=cell_section, &
                     error=error)
      CALL get_cell(cell=cell, abc=abc)
      CALL write_cell(cell=cell, globenv=globenv,error=error)

      CALL cp_subsys_create(subsystem, para_env=globenv%para_env, error=error)
      CALL cp_subsys_read_colvar(subsystem, colvar_section, error=error)

      CALL topology_control(atomic_kind_set, particle_set, &
                            molecule_kind_set, molecule_set, &
                            globenv, error=error)

      CALL read_atomic_kind_set(atomic_kind_set, globenv)

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
!!     - tbmd_env: The tbmd environment to retain
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

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: natom, &
                                                kinds, &
                                                stat, &
                                                i
    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    TYPE(tbmd_model_type), POINTER           :: tbmd_model

!   ------------------------------------------------------------------------

    NULLIFY(atomic_kind_set, molecule_kind_set, molecule_set, particle_set, &
            tbmd_energy, tbmd_force, tbmd_model)

    particle_set => subsystem%particles%els
    atomic_kind_set => subsystem%atomic_kinds%els
    molecule_set => subsystem%molecules_new%els
    molecule_kind_set => subsystem%molecule_kinds_new%els

    natom = SIZE(particle_set)
    CALL tbmd_energy_create(tbmd_energy=tbmd_energy, error=error)
    CALL tbmd_force_create(tbmd_force=tbmd_force, natom=natom, error=error)
    CALL tbmd_model_create(tbmd_model=tbmd_model, error=error)

    CALL tbmd_init_tb_model(tbmd_env=tbmd_env, error=error)

    CALL tbmd_env_set(tbmd_env=tbmd_env, tbmd_energy=tbmd_energy, &
                      tbmd_force=tbmd_force, tbmd_model=tbmd_model, &
                      subsystem=subsystem, error=error)

  END SUBROUTINE tbmd_init_subsys
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
!!     - tbmd_env: The tbmd environment to retain
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
  SUBROUTINE tbmd_init_tb_model(tbmd_env, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_init_tb_model', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: tbmd_input
    TYPE(tbmd_model_type), POINTER           :: tbmd_model
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_ptr
    INTEGER                                  :: i
    CHARACTER(LEN=default_string_length)     :: tb_model_name
    REAL(KIND=dp)                            :: eV, &
                                                Am, &
                                                a_latt

!   ------------------------------------------------------------------------

    NULLIFY(tbmd_input, tbmd_model, atomic_kind_set, atomic_kind_ptr)

    tbmd_input => tbmd_env%tbmd_input
    tbmd_model => tbmd_env%tbmd_model

    atomic_kind_set => tbmd_env%subsystem%atomic_kinds%els
    
    eV = 1/evolt
    Am = 1/angstrom
    a_latt = 5.43*Am

    ! loop over all kinds
    DO i=1,SIZE(atomic_kind_set)
      atomic_kind_ptr => tbmd_env%subsystem%atomic_kinds%els(i)
      CALL get_atomic_kind(atomic_kind=atomic_kind_ptr, &
                           name=tb_model_name)
      SELECT CASE(tb_model_name)
        CASE("SI", "Si")
          SELECT CASE(tbmd_model%tb_model_name)
            CASE(use_gsp)
              tbmd_model%tb_model_name = use_gsp
              tbmd_model%h_ssS         = -1.82_dp*eV
              tbmd_model%h_spS         =  1.96_dp*eV
              tbmd_model%h_ppS         =  3.06_dp*eV
              tbmd_model%h_ppP         = -0.87_dp*eV
              tbmd_model%n             =  2.0_dp
              tbmd_model%r_0           =  SQRT(3.0_dp)/4.0_dp*a_latt
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
            CASE(use_kwon)
              CALL cp_unimplemented_error(fromWhere=routineP, &
                   message="The Kwon tight-binding" //&
                   "parametrization isn't implemented yet.", &
                   error=error, error_level=cp_failure_level)
            CASE(use_lenosky)
              CALL cp_unimplemented_error(fromWhere=routineP, &
                   message="The Lenosky tight-binding" //&
                   "parametrization isn't implemented yet.", &
                   error=error, error_level=cp_failure_level)
            CASE(use_udp)
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

  END SUBROUTINE tbmd_init_tb_model
!***************************************************************************

END MODULE tbmd_environment
