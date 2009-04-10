!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - Merged with the Quickstep MODULE method_specification (17.01.2002,MK)
!>      - USE statements cleaned, added
!>        (25.09.2002,MK)
!>      - Added more LSD structure (01.2003,Joost VandeVondele)
!>      - New molecule data types introduced (Sep. 2003,MK)
!>      - Cleaning; getting rid of pnode (02.10.2003,MK)
!>      - Sub-system setup added (08.10.2003,MK)
!> \author MK (18.05.2000)
! *****************************************************************************
MODULE qs_environment
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
  USE bibliography,                    ONLY: Iannuzzi2007,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             get_cell,&
                                             read_cell,&
                                             write_cell
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_type,&
                                             dftb_control_type,&
                                             gapw_control_type,&
                                             qs_control_type,&
                                             semi_empirical_control_type
  USE cp_control_utils,                ONLY: read_becke_section,&
                                             read_ddapc_section,&
                                             read_dft_control,&
                                             read_mgrid_section,&
                                             read_qs_section,&
                                             read_tddfpt_control,&
                                             write_dft_control,&
                                             write_qs_control
  USE cp_ddapc_types,                  ONLY: cp_ddapc_ewald_create
  USE cp_fm_struct,                    ONLY: optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_methods,               ONLY: cp_subsys_read_colvar
  USE cp_subsys_types,                 ONLY: cp_subsys_create,&
                                             cp_subsys_release,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE cp_symmetry,                     ONLY: write_symmetry
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE et_coupling_types,               ONLY: et_coupling_create
  USE ewald_environment_types,         ONLY: ewald_env_create,&
                                             ewald_env_get,&
                                             ewald_env_release,&
                                             ewald_env_set,&
                                             ewald_environment_type,&
                                             read_ewald_section
  USE ewald_pw_methods,                ONLY: ewald_pw_grid_change
  USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                             ewald_pw_release,&
                                             ewald_pw_type
  USE f77_blas
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_create,&
                                             fist_nonbond_env_type
  USE gamma,                           ONLY: init_md_ftable
  USE harris_energy_types,             ONLY: harris_energy_create,&
                                             harris_energy_type
  USE harris_env_types,                ONLY: harris_env_create,&
                                             harris_env_set,&
                                             harris_env_type
  USE harris_force_types,              ONLY: harris_force_type
  USE hartree_local_methods,           ONLY: init_coulomb_local
  USE header,                          ONLY: dftb_header,&
                                             qs_header,&
                                             se_header
  USE hfx_types,                       ONLY: hfx_create
  USE input_constants,                 ONLY: &
       do_et_becke, do_et_ddapc, do_method_am1, do_method_dftb, &
       do_method_gapw, do_method_gapw_xc, do_method_gpw, do_method_mndo, &
       do_method_mndod, do_method_pdg, do_method_pm3, do_method_pm6, &
       do_method_pnnl, do_method_rm1, do_multipole_none, general_roks, &
       restart_guess, use_aux_fit_basis_set, xc_vdw_fun_pairpot
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
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
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates,&
                                             write_structure_data
  USE qs_dftb_parameters,              ONLY: qs_dftb_param_init
  USE qs_dftb_types,                   ONLY: qs_dftb_pairpot_type
  USE qs_dispersion_pairpot,           ONLY: qs_dispersion_pairpot_init,&
                                             qs_scaling_init,&
                                             write_dispersion
  USE qs_dispersion_types,             ONLY: qs_dispersion_type
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_setup
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_interactions,                 ONLY: init_interaction_radii,&
                                             init_se_nlradius
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             mo_set_p_type
  USE qs_rho0_ggrid,                   ONLY: rho0_s_grid_create
  USE qs_rho0_methods,                 ONLY: init_rho0
  USE qs_rho0_types,                   ONLY: rho0_mpole_type
  USE qs_rho_atom_methods,             ONLY: init_rho_atom
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
  USE semi_empirical_expns3_methods,   ONLY: semi_empirical_expns3_setup
  USE semi_empirical_int_arrays,       ONLY: init_se_intd_array
  USE semi_empirical_mpole_methods,    ONLY: nddo_mpole_setup
  USE semi_empirical_mpole_types,      ONLY: nddo_mpole_type
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_create,&
                                             semi_empirical_si_type
  USE semi_empirical_types,            ONLY: se_taper_create,&
                                             se_taper_type
  USE semi_empirical_utils,            ONLY: se_cutoff_compatible
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE topology,                        ONLY: topology_control
  USE xas_control,                     ONLY: write_xas_control
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_environment'

  ! *** Public subroutines ***
  PUBLIC :: qs_init

CONTAINS

! *****************************************************************************
!> \brief Read the input and the database files for the setup of the
!>      QUICKSTEP environment.
!> \author Creation (22.05.2000,MK)
! *****************************************************************************
  SUBROUTINE qs_init(qs_env,root_section,subsys,cell,cell_ref,qmmm,&
       qmmm_periodic,force_env_section,subsys_section,use_motion_section,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm, qmmm_periodic
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: method_id, natom
    LOGICAL                                  :: do_et, do_hfx, harris_flag, &
                                                my_qmmm, my_qmmm_periodic, &
                                                use_ref_cell
    REAL(kind=dp), DIMENSION(3)              :: abc
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: my_cell, my_cell_ref
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: my_subsys
    TYPE(mol_kind_new_list_type), POINTER    :: mol_kinds
    TYPE(mol_new_list_type), POINTER         :: mols
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER :: cell_section, colvar_section, &
      et_coupling_section, harris_section, hfx_section, kind_section

    NULLIFY(my_subsys, molecule_kind_set, molecule_set, atomic_kind_set,&
         particle_set,atomic_kinds,particles,mols,mol_kinds,my_cell,my_cell_ref,&
         cell_section, harris_section, kind_section)
    para_env=>qs_env%para_env
    IF (.NOT.ASSOCIATED(subsys_section)) THEN
       subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    END IF
    cell_section => section_vals_get_subs_vals(subsys_section,"CELL",error=error)
    harris_section => section_vals_get_subs_vals(force_env_section, "DFT%QS%HARRIS", error=error)
    my_qmmm     = .FALSE.
    IF (PRESENT(qmmm)) my_qmmm=qmmm
    my_qmmm_periodic = .FALSE.
    IF (PRESENT(qmmm_periodic)) my_qmmm_periodic=qmmm_periodic

    harris_flag = .FALSE.
    CALL section_vals_val_get(harris_section, "ACTIVATE", l_val=harris_flag, error=error)
    CALL set_qs_env(qs_env,input=force_env_section,use_harris=harris_flag,error=error)

    ! Possibly initialize arrays for SE
    CALL section_vals_val_get(force_env_section,"DFT%QS%METHOD",i_val=method_id,error=error)
    SELECT CASE (method_id)
    CASE ( do_method_rm1, do_method_am1, do_method_mndo, do_method_pdg,&
           do_method_pm3, do_method_pm6, do_method_mndod, do_method_pnnl )
       CALL init_se_intd_array(error)
    CASE DEFAULT
       ! Do nothing
    END SELECT
    
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
       CALL read_cell(my_cell, my_cell_ref, use_ref_cell=use_ref_cell, &
            cell_section=cell_section, para_env=para_env, error=error)
    END IF
    CALL get_cell ( my_cell, abc=abc)

    !   *** Setup the grids for the G-space Interpolation if any
    CALL cp_ddapc_ewald_create(qs_env%cp_ddapc_ewald, my_qmmm, my_qmmm_periodic,&
         my_cell, force_env_section, subsys_section, para_env, error)

    !   *** Print the cell parameters ***
    CALL write_cell(my_cell,subsys_section,cell_ref=my_cell_ref,error=error)

    !   Kind section
    kind_section => section_vals_get_subs_vals(subsys_section,"KIND",error=error)
    IF (PRESENT(subsys)) THEN
       my_subsys => subsys
       CALL read_atomic_kind_set(my_subsys%atomic_kinds%els,kind_section,para_env,&
            force_env_section,error)
    ELSE
       CALL cp_subsys_create(my_subsys,para_env=para_env,error=error)

       colvar_section => section_vals_get_subs_vals(subsys_section,"COLVAR",error=error)
       CALL cp_subsys_read_colvar( my_subsys, colvar_section, error=error)

       !   *** Read the particle coordinates and allocate the atomic kind, ***
       !   *** the molecule kind, and the molecule data structures         ***
       CALL topology_control (atomic_kind_set,particle_set, molecule_kind_set, molecule_set,&
            my_subsys%colvar_p,my_subsys%gci, root_section, para_env, &
            force_env_section=force_env_section, subsys_section=subsys_section,&
            use_motion_section=use_motion_section, error=error)

       CALL read_atomic_kind_set(atomic_kind_set,kind_section,para_env,&
            force_env_section,error)

       CALL num_ao_el_per_molecule(molecule_kind_set)

       CALL particle_list_create(particles,els_ptr=particle_set,error=error)
       CALL atomic_kind_list_create(atomic_kinds,els_ptr=atomic_kind_set,error=error)
       CALL mol_new_list_create(mols,els_ptr=molecule_set,error=error)
       CALL mol_kind_new_list_create(mol_kinds,els_ptr=molecule_kind_set,error=error)
       CALL cp_subsys_set(my_subsys,particles=particles,atomic_kinds=atomic_kinds,&
            molecules_new=mols,molecule_kinds_new=mol_kinds,error=error)
       CALL particle_list_release(particles,error=error)
       CALL atomic_kind_list_release(atomic_kinds,error=error)
       CALL mol_new_list_release(mols,error=error)
       CALL mol_kind_new_list_release(mol_kinds,error=error)
    END IF

    CALL qs_init_subsys(qs_env,my_subsys,my_cell,my_cell_ref,use_ref_cell,&
         root_section,subsys_section,my_qmmm,harris=harris_flag,error=error)

    IF (.NOT.PRESENT(cell))      CALL cell_release(my_cell,error=error)
    IF (.NOT.PRESENT(cell_ref))  CALL cell_release(my_cell_ref,error=error)
    IF (.NOT.PRESENT(subsys)) CALL cp_subsys_release(my_subsys,error=error)

    do_hfx =.FALSE.
    hfx_section => section_vals_get_subs_vals(qs_env%input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_section,explicit=do_hfx,error=error)
    IF (do_hfx) THEN
      natom=SIZE(particle_set)
      CALL hfx_create(qs_env%x_data, para_env, hfx_section,natom,atomic_kind_set, &
                      qs_env%dft_control, qs_env%cell, error=error)
    END IF
    et_coupling_section => section_vals_get_subs_vals(qs_env%input,&
                           "PROPERTIES%ET_COUPLING",error=error)
    CALL section_vals_get(et_coupling_section,explicit=do_et,error=error)
    IF (do_et) CALL et_coupling_create(qs_env%et_coupling,error=error)
    
  END SUBROUTINE qs_init

! *****************************************************************************
!> \brief Initialize the qs environment (subsys)
!> \author Creation (22.05.2000,MK)
! *****************************************************************************
  SUBROUTINE qs_init_subsys(qs_env,subsys,cell,cell_ref,use_ref_cell,&
       root_section,subsys_section,qmmm,harris,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL, INTENT(in)                      :: use_ref_cell
    TYPE(section_vals_type), POINTER         :: root_section, subsys_section
    LOGICAL, INTENT(in)                      :: qmmm
    LOGICAL, INTENT(IN), OPTIONAL            :: harris
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_init_subsys', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ispin, istat, iw, lmax_sphere, maxl, maxlgto, maxlppl, &
      maxlppnl, method_id, multiplicity, my_ival, n_ao, n_ao_aux_fit, natom, &
      nelectron, output_unit
    INTEGER, DIMENSION(2)                    :: n_mo, nelectron_spin
    LOGICAL                                  :: all_potential_present, &
                                                explicit, harris_flag, &
                                                ionode, was_present
    REAL(dp)                                 :: ewald_rcut, maxocc, &
                                                verlet_skin
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(dftb_control_type), POINTER         :: dftb_control
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: se_nonbond_env
    TYPE(gapw_control_type), POINTER         :: gapw_control
    TYPE(harris_energy_type), POINTER        :: harris_energy
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(nddo_mpole_type), POINTER           :: se_nddo_mpole
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_control_type), POINTER           :: qs_control
    TYPE(qs_dftb_pairpot_type), &
      DIMENSION(:, :), POINTER               :: dftb_potential
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(section_vals_type), POINTER :: dft_section, et_becke_section, &
      et_coupling_section, et_ddapc_section, ewald_section, poisson_section, &
      print_section, qs_section, se_section, vdw_section, xc_fun_section
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_si_type), POINTER    :: se_store_int_env

    CALL timeset(routineN,handle)
    NULLIFY(logger)
    para_env=>qs_env%para_env
    logger => cp_error_get_logger(error)

    ionode = para_env%ionode
    output_unit = -1
    IF(ionode) output_unit = cp_logger_get_default_unit_nr(logger)
    was_present = .FALSE.

    IF (PRESENT(harris)) THEN
       harris_flag = harris
    ELSE
       harris_flag = .FALSE.
    END IF

    ! Initialise the Quickstep environment
    NULLIFY (mos, se_taper, mos_aux_fit)
    NULLIFY (dft_control)
    NULLIFY (energy)
    NULLIFY (force)
    NULLIFY (harris_energy)
    NULLIFY (harris_env)
    NULLIFY (harris_force)
    NULLIFY (local_molecules)
    NULLIFY (local_particles)
    NULLIFY (scf_control)
    NULLIFY (dft_section)
    NULLIFY (et_coupling_section)
    dft_section =>  section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    qs_section =>  section_vals_get_subs_vals(dft_section,"QS",error=error)
    et_coupling_section =>  section_vals_get_subs_vals(qs_env%input,"PROPERTIES%ET_COUPLING",error=error)

    particle_set => subsys%particles%els
    atomic_kind_set => subsys%atomic_kinds%els
    molecule_set => subsys%molecules_new%els
    molecule_kind_set => subsys%molecule_kinds_new%els

    !   *** Read the input section with the DFT control parameters ***
            CALL read_dft_control(dft_control,dft_section,error=error)

    IF (dft_control % do_tddfpt_calculation) THEN
       CALL read_tddfpt_control(dft_control%tddfpt_control,dft_control%nspins, &
            dft_section,error)
    END IF

    !   *** Print the Quickstep program banner (copyright and version number) ***
    iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PROGRAM_BANNER",extension=".Log",error=error)
    CALL section_vals_val_get(qs_section,"METHOD",i_val=method_id,error=error)
    SELECT CASE (method_id)
    CASE DEFAULT
       CALL qs_header(iw, iw>0)
    CASE ( do_method_rm1, do_method_am1, do_method_mndo, do_method_pdg,&
           do_method_pm3, do_method_pm6, do_method_mndod, do_method_pnnl )
       CALL se_header(iw, iw>0)
    CASE ( do_method_dftb )
       CALL dftb_header(iw, iw>0)
    END SELECT

    CALL cp_print_key_finished_output(iw,logger,dft_section,&
         "PRINT%PROGRAM_BANNER",error=error)

    !   *** Read the input section with the Quickstep control parameters ***
    CALL read_qs_section(dft_control%qs_control,qs_section,error=error)

    !   *******  check if any kind of electron transfer calculation has to be performed
    CALL  section_vals_val_get(et_coupling_section,"TYPE_OF_CONSTRAINT",i_val=my_ival,error=error)
    dft_control%qs_control%et_coupling_calc=.FALSE.
    IF (my_ival==do_et_ddapc)THEN
       et_ddapc_section =>  section_vals_get_subs_vals(et_coupling_section,"DDAPC_RESTRAINT_A",error=error)
       dft_control%qs_control%et_coupling_calc=.TRUE.
       dft_control%qs_control%ddapc_restraint=.TRUE.
       CALL  read_ddapc_section(dft_control%qs_control,ddapc_restraint_section=et_ddapc_section,error=error)
    ENDIF

    IF (my_ival==do_et_becke)THEN
       dft_control%qs_control%becke_restraint=.TRUE.
       dft_control%qs_control%et_coupling_calc=.TRUE.       
       et_becke_section =>  section_vals_get_subs_vals(et_coupling_section,"BECKE_RESTRAINT_A",error=error)      
       CALL  read_becke_section(dft_control%qs_control,et_becke_section,error)
    END IF
    CALL read_mgrid_section(dft_control%qs_control,dft_section,error=error)

    !   Create relativistic control section
    CALL rel_c_create(rel_control,error=error)
    CALL rel_c_read_parameters(rel_control,dft_section,error=error)
    CALL set_qs_env(qs_env,rel_control=rel_control,error=error)
    CALL rel_c_release(rel_control,error=error)
  
    !   *** Read DFTB parameter files ***
    IF ( dft_control%qs_control%method == "DFTB" ) THEN
       NULLIFY (ewald_env,ewald_pw,dftb_potential)
       dftb_control => dft_control%qs_control%dftb_control
       CALL qs_dftb_param_init (atomic_kind_set,dftb_control,dftb_potential,&
            subsys_section=subsys_section,para_env=para_env,error=error)
       CALL set_qs_env(qs_env,dftb_potential=dftb_potential,error=error)
       ! check for Ewald
       IF ( dftb_control%do_ewald ) THEN
          CALL ewald_env_create(ewald_env,para_env,error=error)
          poisson_section => section_vals_get_subs_vals(dft_section,"POISSON",error=error)
          CALL ewald_env_set(ewald_env,poisson_section=poisson_section,error=error)
          ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD",error=error)
          print_section => section_vals_get_subs_vals(qs_env%input,&
               "PRINT%GRID_INFORMATION",error=error)
          CALL read_ewald_section(ewald_env,ewald_section,error=error)
          CALL ewald_pw_create(ewald_pw,ewald_env,cell,cell_ref,&
               print_section=print_section,error=error)
          CALL set_qs_env(qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
          CALL ewald_env_release(ewald_env,error=error)
          CALL ewald_pw_release(ewald_pw,error=error)
       END IF
    END IF

    ! DFT+U
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             dft_plus_u_atom_present=dft_control%dft_plus_u)

    !   *** Check basis and fill in missing parts ***
    CALL check_atomic_kind_set(atomic_kind_set,dft_control,para_env,&
         subsys_section=subsys_section,error=error)

    !   *** Check that no all-electron potential is present if GPW or GAPW_XC
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         all_potential_present=all_potential_present)
    IF ( (dft_control%qs_control%method == "GPW") .OR. &
         (dft_control%qs_control%method == "GAPW_XC") ) THEN
       IF( all_potential_present ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "all-el calculations with GPW and GAPW_XC are not implemented ")
       END IF
    END IF

    !   *** Initialize the spherical harmonics and ***
    !   *** the orbital transformation matrices    ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxlgto=maxlgto,maxlppl=maxlppl,maxlppnl=maxlppnl)

    lmax_sphere = dft_control%qs_control%gapw_control%lmax_sphere
    IF(lmax_sphere .LT.0) THEN
       lmax_sphere = 2*maxlgto
       dft_control%qs_control%gapw_control%lmax_sphere= lmax_sphere
    END IF
    maxl = MAX(2*maxlgto,maxlppl,maxlppnl,lmax_sphere) + 1

    CALL init_orbital_pointers(maxl)
    CALL init_spherical_harmonics(maxl,root_section,error)
    !   *** Initialise the atomic kind set ***
    CALL init_atomic_kind_set(atomic_kind_set,para_env,qs_env%input,error=error)

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

    !   *** Initialize the atomic interaction radii ***
    CALL init_interaction_radii(dft_control%qs_control,cell,atomic_kind_set,&
         subsys_section,error)

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
    CALL scf_c_create(scf_control,error=error)
    CALL scf_c_read_parameters(scf_control,dft_section,error=error)

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
       natom = SIZE(particle_set)
       CALL harris_env_create(harris_env=harris_env, natom=natom, &
            nspins=dft_control%nspins, error=error)
       CALL harris_env_set(harris_env=harris_env, harris_energy=harris_energy, &
            harris_force=harris_force, error=error)
    END IF

    !   *** Set the current Quickstep environment ***
    CALL set_qs_env(qs_env,subsys=subsys,error=error)
    CALL set_qs_env(qs_env=qs_env,&
                    cell=cell,&
                    cell_ref=cell_ref,&
                    use_ref_cell=use_ref_cell,&
                    qmmm=qmmm,&
                    dft_control=dft_control,&
                    energy=energy,&
                    force=force,&
                    scf_control=scf_control,&
                    wf_history=wf_history,&
                    error=error)

    IF (harris_flag) THEN
       CALL set_qs_env(qs_env=qs_env, harris_env=harris_env, error=error)
    END IF

    CALL cp_subsys_set(subsys,local_molecules_new=local_molecules,&
         local_particles=local_particles,error=error)

    CALL distribution_1d_release(local_particles,error=error)
    CALL distribution_1d_release(local_molecules,error=error)
    CALL scf_c_release(scf_control,error=error)
    CALL wfi_release(wf_history,error=error)
    CALL dft_control_release(dft_control, error=error)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,&
                    scf_control=scf_control,error=error)

    ! decide what conditions need mo_derivs
    ! right now, this only appears to be OT
    IF (scf_control%use_ot) THEN
       CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.TRUE.,error=error)
    ELSE
       CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.FALSE.,error=error)
    ENDIF

    ! XXXXXXX this is backwards XXXXXXXX
    dft_control%smear = scf_control%smear%do_smear

    !   Initialize the GAPW local densities and potentials
    IF  (dft_control%qs_control%method_id == do_method_gapw .OR. &
         dft_control%qs_control%method_id == do_method_gapw_xc) THEN
       !     *** Allocate and initialize the set of atomic densities ***
       CALL init_rho_atom(qs_env,gapw_control,error=error)
       IF(dft_control%qs_control%method_id /= do_method_gapw_xc) THEN
          !       *** Allocate and initialize the compensation density rho0 ***
          CALL init_rho0(qs_env,gapw_control,cell,error=error)
          !       *** Allocate and Initialize the local coulomb terms (1,2 and 3 centers) ***
          CALL init_coulomb_local(qs_env,gapw_control,error=error)
       END IF
    ELSE IF(dft_control%qs_control%semi_empirical) THEN
       NULLIFY(se_store_int_env, se_nddo_mpole, se_nonbond_env)
       natom = SIZE(particle_set)
       se_section => section_vals_get_subs_vals(qs_section,"SE",error=error)
       se_control => dft_control%qs_control%se_control

       ! Make the cutoff radii choice a bit smarter
       CALL se_cutoff_compatible(se_control, se_section, cell, output_unit, error)

       SELECT CASE ( dft_control%qs_control%method_id)
       CASE DEFAULT
       CASE (do_method_rm1,do_method_am1,do_method_mndo,do_method_pm3,&
             do_method_pm6,do_method_mndod,do_method_pnnl)
          ! Neighbor lists have to be MAX(interaction range, orbital range)
          ! set new kind radius
          CALL init_se_nlradius(se_control,cell,atomic_kind_set,subsys_section,error)
       END SELECT
       ! Initialize to zero the max multipole to treat in the EWALD scheme..
       se_control%max_multipole = do_multipole_none
       ! check for Ewald
       IF (se_control%do_ewald .OR. se_control%do_ewald_gks) THEN
          CALL ewald_env_create(ewald_env,para_env,error=error)
          poisson_section => section_vals_get_subs_vals(dft_section,"POISSON",error=error)
          CALL ewald_env_set(ewald_env,poisson_section=poisson_section,error=error)
          ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD",error=error)
          print_section => section_vals_get_subs_vals(qs_env%input,&
               "PRINT%GRID_INFORMATION",error=error)
          CALL read_ewald_section(ewald_env,ewald_section,error=error)
          ! Create ewald grids
          CALL ewald_pw_create(ewald_pw,ewald_env,cell,cell_ref,&
               print_section=print_section,error=error)
          ! Initialize ewald grids
          CALL ewald_pw_grid_change ( ewald_pw, ewald_env, cell, error )
          ! Setup the nonbond environment (real space part of Ewald) and the maximum
          ! level of multipoles to be treated in the periodic SE scheme
          CALL ewald_env_get(ewald_env, rcut=ewald_rcut, max_multipole=se_control%max_multipole,&
               error=error)
          CALL section_vals_val_get(se_section,"NEIGHBOR_LISTS%VERLET_SKIN",&
               r_val=verlet_skin,error=error)
          CALL fist_nonbond_env_create(se_nonbond_env, atomic_kind_set, nonbonded14=.FALSE.,&
               nonbonded=.TRUE., verlet_skin=verlet_skin, global_rcut=ewald_rcut, error=error)
          ! Create and Setup NDDO multipole environment
          CALL nddo_mpole_setup(se_nddo_mpole, natom, error)
          CALL set_qs_env(qs_env,ewald_env=ewald_env, ewald_pw=ewald_pw,&
               se_nonbond_env=se_nonbond_env, se_nddo_mpole=se_nddo_mpole,&
               error=error)
          CALL ewald_env_release(ewald_env,error=error)
          CALL ewald_pw_release(ewald_pw,error=error)
          ! Handle the residual integral part 1/R^3
          CALL semi_empirical_expns3_setup(atomic_kind_set, se_control,&
               dft_control%qs_control%method_id,error)
       END IF
       ! Taper function
       CALL se_taper_create(se_taper, se_control%integral_screening, &
            se_control%taper_cou, se_control%range_cou, &
            se_control%taper_exc, se_control%range_exc, &
            se_control%taper_scr, se_control%range_scr, error)
       CALL set_qs_env(qs_env, se_taper=se_taper, error=error)
       ! Store integral environment
       CALL semi_empirical_si_create(se_store_int_env, se_section, error)
       CALL set_qs_env(qs_env, se_store_int_env=se_store_int_env, error=error)
    ENDIF

    !   Initialize possible dispersion parameters
    IF (dft_control%qs_control%method_id == do_method_gpw .OR. &
        dft_control%qs_control%method_id == do_method_gapw .OR. &
        dft_control%qs_control%method_id == do_method_gapw_xc) THEN
      ALLOCATE(dispersion_env,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"dispersion_env",1)
      NULLIFY(vdw_section,xc_fun_section)
      vdw_section =>  section_vals_get_subs_vals(dft_section,"XC%vdw_potential",error=error)
      xc_fun_section => section_vals_get_subs_vals(dft_section,"XC%XC_FUNCTIONAL",error=error)
      CALL section_vals_val_get(vdw_section, "POTENTIAL_TYPE", i_val=dispersion_env%type, error=error)
      IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN
        CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%TYPE", c_val=dispersion_env%pp_type, error=error)
        CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%EXP_PRE", r_val=dispersion_env%exp_pre, error=error)
        CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%SCALING" ,explicit=explicit,error=error)
        IF ( .NOT. explicit ) THEN
          CALL qs_scaling_init(dispersion_env%scaling,dft_section,error)
        ELSE
          CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%SCALING", r_val=dispersion_env%scaling, error=error)
        END IF
        CALL section_vals_val_get(vdw_section, "PAIR_POTENTIAL%R_CUTOFF", r_val=dispersion_env%rc_disp, error=error)
        CALL qs_dispersion_pairpot_init(atomic_kind_set,vdw_section,para_env,error)
      END IF
      CALL set_qs_env(qs_env, dispersion_env=dispersion_env, error=error)
    END IF

    !   *** Allocate the MO data types ***
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, nsgf=n_ao, nelectron=nelectron)

    ! the total number of electrons
    nelectron = nelectron - dft_control%charge

    IF (dft_control%multiplicity == 0) THEN
       IF (MODULO(nelectron,2) == 0) THEN
          dft_control%multiplicity = 1
       ELSE
          dft_control%multiplicity = 2
       END IF
    END IF

    multiplicity = dft_control%multiplicity

    IF ((dft_control%nspins < 1).OR.(dft_control%nspins > 2)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            "nspins should be 1 or 2 for the time being...")
    END IF

    IF ((MODULO(nelectron,2) /= 0).AND.(dft_control%nspins == 1)) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
            "Use the LSD option for an odd number of electrons.")
    END IF

    ! The transition potential method to calculate XAS needs LSD
    IF(dft_control%do_xas_calculation ) THEN
       IF(dft_control%nspins == 1) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Use the LSD option for XAS with transition potential.")
       END IF
       IF(dft_control%xas_control%xas_restart .AND. &
            scf_control%density_guess/= restart_guess ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Restart of XAS orbitals possible only if the GS orbitals"//&
               " are also read from restart file.")
       END IF
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
          CALL stop_program(routineN,moduleN,__LINE__,&
               "LSD: try to use a different multiplicity.")
       END IF

       nelectron_spin(1) = (nelectron + multiplicity - 1)/2
       nelectron_spin(2) = (nelectron - multiplicity + 1)/2

       IF (nelectron_spin(2) < 0) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "LSD: too few electrons for this multiplicity.")
       END IF

       n_mo(1) = nelectron_spin(1)
       n_mo(2) = nelectron_spin(2)

    END IF

    scf_control%added_mos(1) = MIN(scf_control%added_mos(1),n_ao-n_mo(1))

    n_mo(1) = n_mo(1) + scf_control%added_mos(1)
    IF (n_mo(2) /= 0 .AND. scf_control%added_mos(2).GT.0) THEN
        scf_control%added_mos(2) = MIN(scf_control%added_mos(2),n_ao-n_mo(2))
        n_mo(2) = n_mo(2) + scf_control%added_mos(2)
        ! we require n_mo(2)<=n_mo(1)
        n_mo(2) = MIN(n_mo(1),n_mo(2))
    ENDIF

    IF (scf_control%smear%do_smear.AND.(scf_control%added_mos(1) == 0)) THEN
      CALL stop_program(routineN,moduleN,__LINE__,&
                        "Added MOs are required for smearing.")
    END IF

    !   *** Some options require that all MOs are computed ... ***
    IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,&
         "PRINT%MO/CARTESIAN",error=error),&
         cp_p_file).OR.&
         (scf_control%level_shift /= 0.0_dp).OR.&
         (scf_control%diagonalization%eps_jacobi /= 0.0_dp).OR.&
         (dft_control%roks.AND.(.NOT.scf_control%use_ot))) THEN
       n_mo(:) = n_ao
    END IF

    ! Compatibility checks for ROKS
    IF (dft_control%roks.AND.(.NOT.scf_control%use_ot)) THEN
       CALL cp_assert((scf_control%roks_scheme /= general_roks),cp_warning_level,&
            cp_assertion_failed,routineP,"General ROKS scheme is not yet tested!"//&
CPSourceFileRef,&
            only_ionode=.TRUE.)

       IF (scf_control%smear%do_smear) THEN
         CALL stop_program(routineN,moduleN,__LINE__,&
                           "The options ROKS and SMEAR are not compatible. "//&
                           "Try UKS instead of ROKS")
       END IF
    END IF

    ALLOCATE (mos(dft_control%nspins),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
            "mos",dft_control%nspins*int_size)
    END IF

    ! in principle the restricted calculation could be performed
    ! using just one set of MOs and special casing most of the code
    ! right now we'll just take care of what is effectively an additional constraint
    ! at as few places as possible, just duplicating the beta orbitals
    IF (dft_control%restricted .AND. logger%para_env%mepos==logger%para_env%source) THEN
       ! it is really not yet tested till the end ! Joost
       iw = cp_logger_get_default_unit_nr(logger,local=.FALSE.)
       WRITE(iw,*) ""
       WRITE(iw,*) " **************************************"
       WRITE(iw,*) " restricted calculation cutting corners"
       WRITE(iw,*) " experimental feature, check code      "
       WRITE(iw,*) " **************************************"
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
    ! If we do wfn fitting, set mo_set_aux_fit
    IF( dft_control%do_wfn_fitting ) THEN
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             nsgf=n_ao_aux_fit,&
                             nelectron=nelectron,&
                             basis_set_id=use_aux_fit_basis_set)
      ALLOCATE (mos_aux_fit(dft_control%nspins),STAT=istat)
      IF (istat /= 0) THEN
         CALL stop_memory(routineN,moduleN,__LINE__,&
              "mos_aux_fit",dft_control%nspins*int_size)
      END IF

      DO ispin=1,dft_control%nspins
         NULLIFY (mos_aux_fit(ispin)%mo_set)
         CALL allocate_mo_set(mo_set=mos_aux_fit(ispin)%mo_set,&
                              nao=n_ao_aux_fit,&
                              nmo=n_mo(ispin),&
                              nelectron=nelectron_spin(ispin),&
                              maxocc=maxocc,&
                              error=error)
      END DO
      CALL set_qs_env(qs_env,mos_aux_fit=mos_aux_fit,error=error)
    END IF


    ! Print the DFT control parameters
    CALL write_dft_control(dft_control,dft_section,error)

    ! Print the DFT control parameters
    IF (dft_control%qs_control%method_id == do_method_gpw .OR. &
        dft_control%qs_control%method_id == do_method_gapw .OR. &
        dft_control%qs_control%method_id == do_method_gapw_xc) THEN
      CALL write_dispersion(qs_env,error)
    END IF

    ! Print the Quickstep control parameters
    CALL write_qs_control(dft_control%qs_control,dft_section,error)

    ! Print XES/XAS control parameters
    IF (dft_control%do_xas_calculation) THEN
       CALL cite_reference(Iannuzzi2007)
       CALL write_xas_control(dft_control%xas_control,dft_section,error=error)
    END IF

    ! Print the unnormalized basis set information (input data)
    CALL write_gto_basis_sets(atomic_kind_set,subsys_section,error)

    ! Print the atomic kind set
    CALL write_atomic_kind_set(atomic_kind_set,subsys_section,error)

    ! Print the molecule kind set
    CALL write_molecule_kind_set(molecule_kind_set,subsys_section,error)

    ! Print the total number of kinds, atoms, basis functions etc.
    CALL write_total_numbers(atomic_kind_set,particle_set,qs_env%input,error)

    ! Print the atomic coordinates
    CALL write_qs_particle_coordinates(particle_set,cell,subsys_section,label="QUICKSTEP",error=error)

    ! Print the interatomic distances
    CALL write_particle_distances(particle_set,cell,subsys_section,error)

    ! Print the requested structure data
    CALL write_structure_data(particle_set,cell,subsys_section,error)

    ! Print symmetry information
    CALL write_symmetry(particle_set,cell,subsys_section,error)

    ! Print the SCF parameters
    CALL scf_c_write_parameters(scf_control,dft_section,error=error)

    ! Sets up pw_env, qs_charges, mpools ...
    CALL qs_env_setup(qs_env,error)

    ! Allocate and Initialie rho0 soft on the global grid
    IF(dft_control%qs_control%method == "GAPW") THEN
       CALL get_qs_env(qs_env=qs_env,rho0_mpole=rho0_mpole,error=error)
       CALL rho0_s_grid_create(qs_env, rho0_mpole, error=error)
    END IF

    IF (output_unit>0) CALL m_flush(output_unit)
    CALL timestop(handle)

  END SUBROUTINE qs_init_subsys

! *****************************************************************************
!> \brief Write the total number of kinds, atoms, etc. to the logical unit
!>      number lunit.
!> \author Creation (06.10.2000)
! *****************************************************************************
  SUBROUTINE write_total_numbers(atomic_kind_set,particle_set,force_env_section,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_total_numbers', &
      routineP = moduleN//':'//routineN

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

       IF ( nset+npgf+ncgf > 0 ) THEN
         WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T71,I10))")&
            "Total number of",&
            "- Atomic kinds:                  ",nkind,&
            "- Atoms:                         ",natom,&
            "- Shell sets:                    ",nset,&
            "- Shells:                        ",nshell,&
            "- Primitive Cartesian functions: ",npgf,&
            "- Cartesian basis functions:     ",ncgf,&
            "- Spherical basis functions:     ",nsgf
       ELSE IF ( nshell+nsgf > 0 ) THEN
         WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T71,I10))")&
            "Total number of",&
            "- Atomic kinds:                  ",nkind,&
            "- Atoms:                         ",natom,&
            "- Shells:                        ",nshell,&
            "- Spherical basis functions:     ",nsgf
       ELSE
         WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T71,I10))")&
            "Total number of",&
            "- Atomic kinds:                  ",nkind,&
            "- Atoms:                         ",natom
       END IF

       IF ((maxlppl > -1).OR.(maxlppnl > -1)) THEN
          WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T75,I6))")&
               "Maximum angular momentum of the",&
               "- Orbital basis functions:                   ",maxlgto,&
               "- Local part of the GTH pseudopotential:     ",maxlppl,&
               "- Non-local part of the GTH pseudopotential: ",maxlppnl
       ELSE
          WRITE (UNIT=output_unit,FMT="(/,T3,A,T75,I6)")&
               "Maximum angular momentum of the orbital basis functions: ",maxlgto
       END IF

    END IF
    CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
         "PRINT%TOTAL_NUMBERS",error=error)

  END SUBROUTINE write_total_numbers

END MODULE qs_environment
