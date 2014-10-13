!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE atom_kind_orbitals,              ONLY: calculate_atomic_relkin
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE basis_set_types,                 ONLY: copy_gto_basis_set,&
                                             gto_basis_set_type
  USE bibliography,                    ONLY: Iannuzzi2006,&
                                             Iannuzzi2007,&
                                             cite_reference
  USE cell_methods,                    ONLY: cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_retain,&
                                             cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_type,&
                                             dftb_control_type,&
                                             gapw_control_type,&
                                             qs_control_type,&
                                             scptb_control_type,&
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
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE cp_symmetry,                     ONLY: write_symmetry
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_type
  USE distribution_methods,            ONLY: distribute_molecules_1d
  USE dm_ls_scf,                       ONLY: ls_scf_create
  USE et_coupling_types,               ONLY: et_coupling_create
  USE ewald_environment_types,         ONLY: ewald_env_create,&
                                             ewald_env_get,&
                                             ewald_env_release,&
                                             ewald_env_set,&
                                             ewald_environment_type,&
                                             read_ewald_section
  USE ewald_pw_methods,                ONLY: ewald_pw_grid_update
  USE ewald_pw_types,                  ONLY: ewald_pw_create,&
                                             ewald_pw_release,&
                                             ewald_pw_type
  USE f77_blas
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_create,&
                                             fist_nonbond_env_type
  USE ga_environment_types,            ONLY: ga_environment_type,&
                                             init_ga_env
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
                                             scptb_header,&
                                             se_header
  USE hfx_ri_methods,                  ONLY: hfx_ri_env_create
  USE hfx_types,                       ONLY: hfx_create
  USE input_constants,                 ONLY: &
       dispersion_d3, do_et_becke, do_et_ddapc, do_method_am1, &
       do_method_dftb, do_method_gapw, do_method_gapw_xc, do_method_gpw, &
       do_method_lrigpw, do_method_mndo, do_method_mndod, do_method_ofgpw, &
       do_method_pdg, do_method_pm3, do_method_pm6, do_method_pnnl, &
       do_method_rm1, do_method_scptb, do_multipole_none, general_roks, &
       rel_none, rel_trans_atom, restart_guess, use_aux_fit_basis_set, &
       use_lri_basis_set, vdw_pairpot_dftd3, xc_vdw_fun_none, &
       xc_vdw_fun_nonloc, xc_vdw_fun_pairpot
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kg_environment,                  ONLY: kg_env_create
  USE kinds,                           ONLY: dp
  USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                             kpoint_initialize,&
                                             kpoint_initialize_mos
  USE kpoint_types,                    ONLY: kpoint_type,&
                                             read_kpoint_section,&
                                             write_kpoint_info
  USE lri_environment_types,           ONLY: lri_basis_init,&
                                             lri_env_create,&
                                             lri_environment_type
  USE machine,                         ONLY: m_flush
  USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                             write_molecule_kind_set
  USE molecule_types_new,              ONLY: molecule_type
  USE mp2_setup,                       ONLY: read_mp2_section
  USE mp2_types,                       ONLY: mp2_env_create
  USE orbital_pointers,                ONLY: init_orbital_pointers
  USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
  USE particle_types,                  ONLY: particle_type,&
                                             write_particle_distances,&
                                             write_qs_particle_coordinates,&
                                             write_structure_data
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE qs_dftb_parameters,              ONLY: qs_dftb_param_init
  USE qs_dftb_types,                   ONLY: qs_dftb_pairpot_type
  USE qs_dispersion_nonloc,            ONLY: qs_dispersion_nonloc_init
  USE qs_dispersion_pairpot,           ONLY: qs_dispersion_pairpot_init
  USE qs_dispersion_types,             ONLY: qs_dispersion_type
  USE qs_dispersion_utils,             ONLY: qs_dispersion_env_set,&
                                             qs_write_dispersion
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_methods,          ONLY: qs_env_setup
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_interactions,                 ONLY: init_interaction_radii,&
                                             init_se_nlradius,&
                                             write_core_charge_radii,&
                                             write_geminal_radii,&
                                             write_paw_radii,&
                                             write_pgf_orb_radii,&
                                             write_ppl_radii,&
                                             write_ppnl_radii
  USE qs_kind_types,                   ONLY: &
       check_qs_kind_set, get_qs_kind, get_qs_kind_set, init_gapw_basis_set, &
       init_qs_kind_set, qs_kind_type, set_qs_kind, write_gto_basis_sets, &
       write_qs_kind_set
  USE qs_ks_types,                     ONLY: qs_ks_env_create,&
                                             qs_ks_env_type,&
                                             qs_ks_release,&
                                             set_ks_env
  USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                             mo_set_p_type
  USE qs_rho0_ggrid,                   ONLY: rho0_s_grid_create
  USE qs_rho0_methods,                 ONLY: init_rho0
  USE qs_rho0_types,                   ONLY: rho0_mpole_type
  USE qs_rho_atom_methods,             ONLY: init_rho_atom
  USE qs_subsys_methods,               ONLY: qs_subsys_create
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_release,&
                                             qs_subsys_set,&
                                             qs_subsys_type
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
  USE scptb_utils,                     ONLY: scptb_parameter_init
  USE se_ga_tools,                     ONLY: se_ga_pair_list_init
  USE semi_empirical_expns3_methods,   ONLY: semi_empirical_expns3_setup
  USE semi_empirical_int_arrays,       ONLY: init_se_intd_array
  USE semi_empirical_mpole_methods,    ONLY: nddo_mpole_setup
  USE semi_empirical_mpole_types,      ONLY: nddo_mpole_type
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_create,&
                                             semi_empirical_si_type
  USE semi_empirical_types,            ONLY: se_taper_create,&
                                             se_taper_type
  USE semi_empirical_utils,            ONLY: se_cutoff_compatible
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xas_control,                     ONLY: write_xas_control
#include "./common/cp_common_uses.f90"

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
!> \param qs_env ...
!> \param root_section ...
!> \param cp_subsys ...
!> \param cell ...
!> \param cell_ref ...
!> \param qmmm ...
!> \param qmmm_periodic ...
!> \param qmmm_env_qm ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \param error ...
!> \author Creation (22.05.2000,MK)
! *****************************************************************************
  SUBROUTINE qs_init(qs_env,root_section,cp_subsys,cell,cell_ref,qmmm,&
       qmmm_periodic,qmmm_env_qm,force_env_section,subsys_section,&
       use_motion_section,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: cp_subsys
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm, qmmm_periodic
    TYPE(qmmm_env_qm_type), OPTIONAL, &
      POINTER                                :: qmmm_env_qm
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, istat, method_id, &
                                                natom, nkind
    LOGICAL :: do_et, do_exx, do_hfx, do_hfx_ri, do_kp, failure, harris_flag, &
      mp2_present, my_qmmm, my_qmmm_multipole, my_qmmm_periodic, use_ref_cell
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rtmat
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: my_cell, my_cell_ref
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(ga_environment_type), POINTER       :: ga_env
    TYPE(gto_basis_set_type), POINTER        :: lri_basis_set, orb_basis_set
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(lri_environment_type), POINTER      :: lri_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(section_vals_type), POINTER :: dft_section, et_coupling_section, &
      harris_section, hfx_ri_section, hfx_section, kpoint_section, mp2_section

    failure = .FALSE.
    NULLIFY(my_cell, my_cell_ref, atomic_kind_set, particle_set, &
         qs_kind_set, harris_section, kpoint_section, dft_section, ga_env, &
         subsys, ks_env)


    para_env=>qs_env%para_env
    IF (.NOT.ASSOCIATED(subsys_section)) THEN
       subsys_section => section_vals_get_subs_vals(force_env_section,"SUBSYS",error=error)
    END IF

    harris_section => section_vals_get_subs_vals(force_env_section, "DFT%QS%HARRIS", error=error)
    my_qmmm     = .FALSE.
    IF (PRESENT(qmmm)) my_qmmm=qmmm
    my_qmmm_periodic = .FALSE.
    IF (PRESENT(qmmm_periodic)) my_qmmm_periodic=qmmm_periodic
    my_qmmm_multipole = .FALSE.
    IF (PRESENT(qmmm_env_qm)) my_qmmm_multipole=qmmm_env_qm%multipole
    IF(PRESENT(qmmm_env_qm)) THEN
     qs_env%qmmm_env_qm => qmmm_env_qm
    END IF

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

    CALL qs_subsys_create(subsys, para_env, root_section=root_section,&
                          force_env_section=force_env_section,&
                          subsys_section=subsys_section,&
                          use_motion_section=use_motion_section,&
                          cp_subsys=cp_subsys, cell=cell, cell_ref=cell_ref, &
                          error=error)

    CALL qs_ks_env_create(ks_env, error=error)
    CALL set_ks_env(ks_env, subsys=subsys, error=error)
    CALL set_qs_env(qs_env, ks_env=ks_env, error=error)

    CALL qs_subsys_get(subsys,&
                       cell=my_cell,&
                       cell_ref=my_cell_ref,&
                       use_ref_cell=use_ref_cell,&
                       atomic_kind_set=atomic_kind_set,&
                       qs_kind_set=qs_kind_set,&
                       particle_set=particle_set,&
                       error=error)

    !   *** Setup the grids for the G-space Interpolation if any
    CALL cp_ddapc_ewald_create(qs_env%cp_ddapc_ewald, my_qmmm, my_qmmm_periodic,&
         my_qmmm_multipole, my_cell, force_env_section, subsys_section, para_env, error)

    CALL qs_init_subsys(qs_env,subsys,my_cell,my_cell_ref,use_ref_cell,&
         root_section,subsys_section,my_qmmm,harris=harris_flag,error=error)

    ! kpoints
    CALL get_qs_env(qs_env=qs_env,kpoints=kpoints,error=error)
    kpoint_section => section_vals_get_subs_vals(qs_env%input,"DFT%KPOINTS",error=error)
    CALL read_kpoint_section(kpoints,kpoint_section,error)
    CALL kpoint_initialize(kpoints, particle_set, my_cell, error)
    dft_section => section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
    CALL write_kpoint_info(kpoints,dft_section,error)
    kpoints%para_env => para_env
    CALL cp_para_env_retain(para_env, error)
    CALL get_qs_env(qs_env=qs_env,blacs_env=blacs_env,error=error)
    kpoints%blacs_env_all => blacs_env
    CALL cp_blacs_env_retain(blacs_env, error)
    CALL get_qs_env(qs_env=qs_env,do_kpoints=do_kp,error=error)
    IF(do_kp) THEN
       CALL kpoint_env_initialize(kpoints, error)
       CALL kpoint_initialize_mos(kpoints, qs_env%mos, error)
    END IF

    do_hfx =.FALSE.
    hfx_section => section_vals_get_subs_vals(qs_env%input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_section,explicit=do_hfx,error=error)
    IF (do_hfx) THEN
      ! Retrieve particle_set and atomic_kind_set (needed for both kinds of initialization)
      natom=SIZE(particle_set)
      CALL hfx_create(qs_env%x_data, para_env, hfx_section, natom, atomic_kind_set,&
                      qs_kind_set, qs_env%dft_control, my_cell, error=error)
      hfx_ri_section => section_vals_get_subs_vals(hfx_section,"HFX_RI",error=error)
      CALL section_vals_get(hfx_ri_section,explicit=do_hfx_ri,error=error)
      IF (do_hfx_ri) THEN
        CALL hfx_ri_env_create(qs_env,error)
      END IF
    END IF

    mp2_section => section_vals_get_subs_vals(qs_env%input,"DFT%XC%WF_CORRELATION",error=error)
    CALL section_vals_get(mp2_section,explicit=mp2_present,error=error)
    IF (mp2_present) THEN
       CALL mp2_env_create(qs_env%mp2_env,error)
       CALL read_mp2_section(qs_env%input,qs_env%mp2_env,error)
       ! create the EXX section if necessary
       do_exx =.FALSE.
       hfx_section => section_vals_get_subs_vals(qs_env%input,"DFT%XC%WF_CORRELATION%RI_RPA%HF",error=error)
       CALL section_vals_get(hfx_section,explicit=do_exx,error=error)
       IF (do_exx) THEN
         ! Retrieve particle_set and atomic_kind_set (needed for both kinds of initialization)
         natom=SIZE(particle_set)
         CALL hfx_create(qs_env%mp2_env%ri_rpa%x_data, para_env, hfx_section, natom, atomic_kind_set,&
                         qs_kind_set, qs_env%dft_control, my_cell, error=error)
       END IF
    END IF

    IF (qs_env%dft_control%qs_control%do_kg) THEN
      CALL cite_reference(Iannuzzi2006)
      CALL kg_env_create(qs_env%kg_env, qs_env%input, error)
    END IF

    et_coupling_section => section_vals_get_subs_vals(qs_env%input,&
                           "PROPERTIES%ET_COUPLING",error=error)
    CALL section_vals_get(et_coupling_section,explicit=do_et,error=error)
    IF (do_et) CALL et_coupling_create(qs_env%et_coupling,error=error)

    ! lri env
    IF (method_id == do_method_lrigpw.OR.qs_env%dft_control%qs_control%lri_optbas) THEN
       NULLIFY(lri_env)
       CALL lri_env_create(lri_env,error)
       ! initialize the basic basis sets (orb and ri)
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,error=error)
       nkind = SIZE(atomic_kind_set)
       ALLOCATE(lri_env%orb_basis(nkind),lri_env%ri_basis(nkind),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DO ikind=1,nkind
          NULLIFY(orb_basis_set,lri_basis_set)
          CALL get_qs_kind(qs_kind_set(ikind),&
                           orb_basis_set=orb_basis_set,&
                           lri_basis_set=lri_basis_set)
          NULLIFY(lri_env%orb_basis(ikind)%gto_basis_set)
          NULLIFY(lri_env%ri_basis(ikind)%gto_basis_set)
          IF (ASSOCIATED(orb_basis_set)) THEN
             CALL copy_gto_basis_set(orb_basis_set,lri_env%orb_basis(ikind)%gto_basis_set,error)
             CALL copy_gto_basis_set(lri_basis_set,lri_env%ri_basis(ikind)%gto_basis_set,error)
          END IF
       END DO
       CALL lri_basis_init(lri_env,atomic_kind_set,error)
       !
       ! check for debug
       CALL section_vals_val_get(force_env_section,"DFT%QS%DEBUG_LRI_INTEGRALS",&
                                 l_val=lri_env%debug,error=error)
       !
       CALL set_qs_env (qs_env,lri_env=lri_env,error=error)
    END IF

    ! GA option
    CALL init_ga_env ( ga_env, error )
    CALL set_qs_env ( qs_env, ga_env=ga_env, error=error )
    CALL se_ga_pair_list_init ( qs_env, error )

    IF (qs_env%dft_control%qs_control%do_ls_scf) THEN
       CALL ls_scf_create(qs_env,error)
    ENDIF

    ! see if we have atomic relativistic corrections
    CALL get_qs_env(qs_env,rel_control=rel_control,error=error)
    IF (rel_control%rel_method /= rel_none) THEN
       IF (rel_control%rel_transformation == rel_trans_atom) THEN
          nkind = SIZE(atomic_kind_set)
          DO ikind=1,nkind
             NULLIFY(rtmat)
             CALL calculate_atomic_relkin(atomic_kind_set(ikind), qs_kind_set(ikind),rel_control,rtmat,error)
             IF(ASSOCIATED(rtmat)) CALL set_qs_kind(qs_kind_set(ikind), reltmat=rtmat)
          END DO
       END IF
    END IF

    CALL qs_subsys_release(subsys, error=error)
    CALL qs_ks_release(ks_env, error=error)
  END SUBROUTINE qs_init

! *****************************************************************************
!> \brief Initialize the qs environment (subsys)
!> \param qs_env ...
!> \param subsys ...
!> \param cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param root_section ...
!> \param subsys_section ...
!> \param qmmm ...
!> \param harris ...
!> \param error ...
!> \author Creation (22.05.2000,MK)
! *****************************************************************************
  SUBROUTINE qs_init_subsys(qs_env,subsys,cell,cell_ref,use_ref_cell,&
       root_section,subsys_section,qmmm,harris,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(cell_type), POINTER                 :: cell, cell_ref
    LOGICAL, INTENT(in)                      :: use_ref_cell
    TYPE(section_vals_type), POINTER         :: root_section, subsys_section
    LOGICAL, INTENT(in)                      :: qmmm
    LOGICAL, INTENT(IN), OPTIONAL            :: harris
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_init_subsys', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ispin, istat, iw, lmax_sphere, maxl, maxlgto, &
      maxlgto_lri, maxlppl, maxlppnl, method_id, multiplicity, my_ival, n_ao, &
      n_ao_aux_fit, n_mo_add, natom, nelectron, output_unit
    INTEGER, DIMENSION(2)                    :: n_mo, nelectron_spin
    LOGICAL                                  :: all_potential_present, &
                                                failure, harris_flag, &
                                                has_unit_metric, was_present
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
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_wf_history_type), POINTER        :: wf_history
    TYPE(rel_control_type), POINTER          :: rel_control
    TYPE(rho0_mpole_type), POINTER           :: rho0_mpole
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(scptb_control_type), POINTER        :: scptb_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(section_vals_type), POINTER :: dft_section, et_becke_section, &
      et_coupling_section, et_ddapc_section, ewald_section, nl_section, &
      poisson_section, pp_section, print_section, qs_section, se_section, &
      xc_section
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(semi_empirical_si_type), POINTER    :: se_store_int_env

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(logger)
    para_env=>qs_env%para_env
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)
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

    CALL qs_subsys_get(subsys, particle_set=particle_set,&
                       qs_kind_set=qs_kind_set,&
                       atomic_kind_set=atomic_kind_set,&
                       molecule_set=molecule_set,&
                       molecule_kind_set=molecule_kind_set,&
                       error=error)

    !   *** Read the input section with the DFT control parameters ***
    CALL read_dft_control(dft_control,dft_section,error=error)

    IF (dft_control % do_tddfpt_calculation) THEN
       CALL read_tddfpt_control(dft_control%tddfpt_control, &
            dft_section,error)
    END IF

    !   *** Print the Quickstep program banner (copyright and version number) ***
    iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%PROGRAM_BANNER",extension=".Log",error=error)
    CALL section_vals_val_get(qs_section,"METHOD",i_val=method_id,error=error)
    SELECT CASE (method_id)
    CASE DEFAULT
       CALL qs_header(iw)
    CASE ( do_method_rm1, do_method_am1, do_method_mndo, do_method_pdg,&
           do_method_pm3, do_method_pm6, do_method_mndod, do_method_pnnl )
       CALL se_header(iw)
    CASE ( do_method_dftb )
       CALL dftb_header(iw)
    CASE ( do_method_scptb )
       CALL scptb_header(iw)
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
    CALL read_mgrid_section(dft_control%qs_control,dft_section,para_env=para_env,error=error)

    !   Create relativistic control section
    CALL rel_c_create(rel_control,error=error)
    CALL rel_c_read_parameters(rel_control,dft_section,error=error)
    CALL set_qs_env(qs_env,rel_control=rel_control,error=error)
    CALL rel_c_release(rel_control,error=error)

    !   *** Read DFTB parameter files ***
    IF ( dft_control%qs_control%method == "DFTB" ) THEN
       NULLIFY (ewald_env,ewald_pw,dftb_potential)
       dftb_control => dft_control%qs_control%dftb_control
       CALL qs_dftb_param_init (atomic_kind_set,qs_kind_set,dftb_control,dftb_potential,&
            subsys_section=subsys_section,para_env=para_env,error=error)
       CALL set_qs_env(qs_env,dftb_potential=dftb_potential,error=error)
       ! check for Ewald
       IF ( dftb_control%do_ewald ) THEN
          CALL ewald_env_create(ewald_env,para_env,error=error)
          poisson_section => section_vals_get_subs_vals(dft_section,"POISSON",error=error)
          CALL ewald_env_set(ewald_env,poisson_section=poisson_section,error=error)
          ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD",error=error)
          print_section => section_vals_get_subs_vals(qs_env%input,"PRINT%GRID_INFORMATION",error=error)
          CALL read_ewald_section(ewald_env,ewald_section,error=error)
          CALL ewald_pw_create(ewald_pw,ewald_env,cell,cell_ref,print_section=print_section,error=error)
          CALL set_qs_env(qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
          CALL ewald_env_release(ewald_env,error=error)
          CALL ewald_pw_release(ewald_pw,error=error)
       END IF
    ELSE IF ( dft_control%qs_control%method == "SCPTB" ) THEN
       scptb_control => dft_control%qs_control%scptb_control
       print_section => section_vals_get_subs_vals(subsys_section,"PRINT",error=error)
       CALL scptb_parameter_init(atomic_kind_set,qs_kind_set,scptb_control,print_section,para_env,error)
       ! check for Ewald
       IF ( scptb_control%do_ewald ) THEN
          CALL ewald_env_create(ewald_env,para_env,error=error)
          poisson_section => section_vals_get_subs_vals(dft_section,"POISSON",error=error)
          CALL ewald_env_set(ewald_env,poisson_section=poisson_section,error=error)
          ewald_section => section_vals_get_subs_vals(poisson_section,"EWALD",error=error)
          print_section => section_vals_get_subs_vals(qs_env%input,"PRINT%GRID_INFORMATION",error=error)
          CALL read_ewald_section(ewald_env,ewald_section,error=error)
          CALL ewald_pw_create(ewald_pw,ewald_env,cell,cell_ref,print_section=print_section,error=error)
          CALL set_qs_env(qs_env,ewald_env=ewald_env,ewald_pw=ewald_pw,error=error)
          CALL ewald_env_release(ewald_env,error=error)
          CALL ewald_pw_release(ewald_pw,error=error)
       END IF
    END IF

    ! DFT+U
    CALL get_qs_kind_set(qs_kind_set, dft_plus_u_atom_present=dft_control%dft_plus_u)

    !   *** Check basis and fill in missing parts ***
    CALL check_qs_kind_set(qs_kind_set,dft_control,para_env,&
         subsys_section=subsys_section,error=error)

    !   *** Check that no all-electron potential is present if GPW or GAPW_XC
    CALL get_qs_kind_set(qs_kind_set, all_potential_present=all_potential_present)
    IF ( (dft_control%qs_control%method == "GPW") .OR. &
         (dft_control%qs_control%method == "GAPW_XC") .OR. &
         (dft_control%qs_control%method == "OFGPW") ) THEN
       IF( all_potential_present ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "all-el calculations with GPW, GAPW_XC, and OFGPW are not implemented ",&
               para_env)
       END IF
    END IF

    !   *** Initialize the spherical harmonics and ***
    !   *** the orbital transformation matrices    ***
    CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto, maxlppl=maxlppl, maxlppnl=maxlppnl)

    lmax_sphere = dft_control%qs_control%gapw_control%lmax_sphere
    IF(lmax_sphere .LT.0) THEN
       lmax_sphere = 2*maxlgto
       dft_control%qs_control%gapw_control%lmax_sphere= lmax_sphere
    END IF
    IF(dft_control%qs_control%method == "LRIGPW".OR.dft_control%qs_control%lri_optbas) THEN
       CALL get_qs_kind_set(qs_kind_set,maxlgto=maxlgto_lri,basis_set_id=use_lri_basis_set)
       !take maxlgto from lri basis if larger (usually) 
       maxlgto = MAX(maxlgto,maxlgto_lri)
    END IF
    maxl = MAX(2*maxlgto,maxlppl,maxlppnl,lmax_sphere) + 1

    CALL init_orbital_pointers(maxl)

    output_unit = cp_print_key_unit_nr(logger,root_section,"GLOBAL%PRINT/SPHERICAL_HARMONICS",&
           extension=".Log",error=error)
    CALL init_spherical_harmonics(maxl,output_unit,error)
    CALL cp_print_key_finished_output(output_unit,logger,root_section,&
            "GLOBAL%PRINT/SPHERICAL_HARMONICS",error=error)

    !   *** Initialise the qs_kind_set ***
    CALL init_qs_kind_set(qs_kind_set,error)

    !   *** Initialise GAPW soft basis and projectors
    IF(dft_control%qs_control%method == "GAPW" .OR. &
         dft_control%qs_control%method == "GAPW_XC")  THEN
       qs_control => dft_control%qs_control
       gapw_control => dft_control%qs_control%gapw_control
       CALL init_gapw_basis_set(qs_kind_set,qs_control,qs_env%input,error)
    ENDIF

    !   *** Initialize the pretabulation for the calculation of the   ***
    !   *** incomplete Gamma function F_n(t) after McMurchie-Davidson ***
    CALL get_qs_kind_set(qs_kind_set, maxlgto=maxlgto)
    maxl = MAX(3*maxlgto + 1,0)
    CALL init_md_ftable(maxl)

    !   *** Initialize the atomic interaction radii ***
    CALL init_interaction_radii(dft_control%qs_control,atomic_kind_set,qs_kind_set,error)

    CALL write_pgf_orb_radii("orb",atomic_kind_set,qs_kind_set,subsys_section,error)
    CALL write_geminal_radii(atomic_kind_set,qs_kind_set,subsys_section,error)
    CALL write_pgf_orb_radii("aux",atomic_kind_set,qs_kind_set,subsys_section,error)
    CALL write_pgf_orb_radii("lri",atomic_kind_set,qs_kind_set,subsys_section,error)
    CALL write_core_charge_radii(atomic_kind_set,qs_kind_set,subsys_section,error)
    CALL write_ppl_radii(atomic_kind_set,qs_kind_set,subsys_section,error)
    CALL write_ppnl_radii(atomic_kind_set,qs_kind_set,subsys_section,error)
    CALL write_paw_radii(atomic_kind_set,qs_kind_set,subsys_section,error)


    !   *** Distribute molecules and atoms using the new data structures ***
    CALL distribute_molecules_1d(atomic_kind_set=atomic_kind_set,&
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

    !   *** Allocate the data structure for Quickstep energies ***
    CALL allocate_qs_energy(energy)

    !   *** Allocate the data structure for the Harris energies ***
    IF (harris_flag) THEN
       CALL harris_energy_create(harris_energy=harris_energy,error=error)
    END IF

    ! check for orthogoanl basis
    has_unit_metric = .FALSE.
    IF (dft_control%qs_control%semi_empirical) THEN
      IF (dft_control%qs_control%se_control%orthogonal_basis) has_unit_metric = .TRUE.
    END IF
    IF (dft_control%qs_control%dftb) THEN
      IF (dft_control%qs_control%dftb_control%orthogonal_basis) has_unit_metric = .TRUE.
    END IF
    CALL set_qs_env(qs_env,has_unit_metric=has_unit_metric,error=error)

    !   *** Activate the interpolation ***
    CALL wfi_create(wf_history,&
                    interpolation_method_nr=&
                    dft_control%qs_control%wf_interpolation_method_nr,&
                    extrapolation_order = dft_control%qs_control%wf_extrapolation_order,&
                    has_unit_metric = has_unit_metric, &
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
    CALL set_qs_env(qs_env=qs_env,&
                    qmmm=qmmm,&
                    dft_control=dft_control,&
                    scf_control=scf_control,&
                    wf_history=wf_history,&
                    error=error)

    CALL qs_subsys_set(subsys,&
                       cell_ref=cell_ref,&
                       use_ref_cell=use_ref_cell,&
                       energy=energy,&
                       force=force,&
                       error=error)

    IF (harris_flag) THEN
       CALL set_qs_env(qs_env=qs_env, harris_env=harris_env, error=error)
    END IF

    CALL qs_subsys_set(subsys,local_molecules_new=local_molecules,&
         local_particles=local_particles,cell=cell,error=error)

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
    IF (qs_env%dft_control%qs_control%do_ls_scf .OR. &
        qs_env%dft_control%qs_control%do_almo_scf) THEN
      CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.FALSE.,error=error)
    ELSE
      IF (scf_control%use_ot) THEN
         CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.TRUE.,error=error)
      ELSE
         CALL set_qs_env(qs_env=qs_env,requires_mo_derivs=.FALSE.,error=error)
      ENDIF
    ENDIF

    ! XXXXXXX this is backwards XXXXXXXX
    dft_control%smear = scf_control%smear%do_smear

    !   Initialize the GAPW local densities and potentials
    IF  (dft_control%qs_control%method_id == do_method_gapw .OR. &
         dft_control%qs_control%method_id == do_method_gapw_xc) THEN
       !     *** Allocate and initialize the set of atomic densities ***
       CALL init_rho_atom(qs_env,gapw_control,error=error)
       IF(dft_control%qs_control%method_id /= do_method_gapw_xc) THEN
          CALL get_qs_env(qs_env=qs_env,natom=natom,error=error)
          !       *** Allocate and initialize the compensation density rho0 ***
          CALL init_rho0(qs_env,gapw_control,error=error)
          !       *** Allocate and Initialize the local coulomb term ***
          CALL init_coulomb_local(qs_env%hartree_local,natom,error=error)
       END IF
    ELSE IF (dft_control%qs_control%method_id == do_method_lrigpw) THEN
       ! allocate local ri environment
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
          CALL init_se_nlradius(se_control,atomic_kind_set,qs_kind_set,subsys_section,error)
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
          CALL ewald_pw_grid_update(ewald_pw, ewald_env, cell%hmat, error)
          ! Setup the nonbond environment (real space part of Ewald)
          CALL ewald_env_get(ewald_env, rcut=ewald_rcut, error=error)
          ! Setup the maximum level of multipoles to be treated in the periodic SE scheme
          IF (se_control%do_ewald) THEN
             CALL ewald_env_get(ewald_env, max_multipole=se_control%max_multipole, error=error)
          ENDIF
          CALL section_vals_val_get(se_section,"NEIGHBOR_LISTS%VERLET_SKIN",&
               r_val=verlet_skin,error=error)
          CALL fist_nonbond_env_create(se_nonbond_env, atomic_kind_set, &
               do_nonbonded=.TRUE., verlet_skin=verlet_skin, ewald_rcut=ewald_rcut, &
               ei_scale14=0.0_dp, vdw_scale14=0.0_dp, shift_cutoff=.FALSE., &
               error=error)
          ! Create and Setup NDDO multipole environment
          CALL nddo_mpole_setup(se_nddo_mpole, natom, error)
          CALL set_qs_env(qs_env,ewald_env=ewald_env, ewald_pw=ewald_pw,&
               se_nonbond_env=se_nonbond_env, se_nddo_mpole=se_nddo_mpole,&
               error=error)
          CALL ewald_env_release(ewald_env,error=error)
          CALL ewald_pw_release(ewald_pw,error=error)
          ! Handle the residual integral part 1/R^3
          CALL semi_empirical_expns3_setup(qs_kind_set, se_control,&
               dft_control%qs_control%method_id,error)
       END IF
       ! Taper function
       CALL se_taper_create(se_taper, se_control%integral_screening, se_control%do_ewald,&
            se_control%taper_cou, se_control%range_cou, &
            se_control%taper_exc, se_control%range_exc, &
            se_control%taper_scr, se_control%range_scr, &
            se_control%taper_lrc, se_control%range_lrc, error)
       CALL set_qs_env(qs_env, se_taper=se_taper, error=error)
       ! Store integral environment
       CALL semi_empirical_si_create(se_store_int_env, se_section, error=error)
       CALL set_qs_env(qs_env, se_store_int_env=se_store_int_env, error=error)
    ENDIF

    !   Initialize possible dispersion parameters
    IF (dft_control%qs_control%method_id == do_method_gpw .OR. &
        dft_control%qs_control%method_id == do_method_gapw .OR. &
        dft_control%qs_control%method_id == do_method_gapw_xc .OR. &
        dft_control%qs_control%method_id == do_method_lrigpw .OR. &
        dft_control%qs_control%method_id == do_method_ofgpw) THEN
      ALLOCATE(dispersion_env,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      NULLIFY(xc_section)
      xc_section => section_vals_get_subs_vals(dft_section,"XC",error=error)
      CALL qs_dispersion_env_set(dispersion_env,xc_section,error)
      IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN
        NULLIFY(pp_section)
        pp_section =>  section_vals_get_subs_vals(xc_section,"VDW_POTENTIAL%PAIR_POTENTIAL",error=error)
        CALL qs_dispersion_pairpot_init(atomic_kind_set,qs_kind_set,dispersion_env,pp_section,para_env,error)
      ELSE IF ( dispersion_env%type == xc_vdw_fun_nonloc ) THEN
        NULLIFY(nl_section)
        nl_section =>  section_vals_get_subs_vals(xc_section,"VDW_POTENTIAL%NON_LOCAL",error=error)
        CALL qs_dispersion_nonloc_init(dispersion_env,nl_section,para_env,error)
      END IF
      CALL set_qs_env(qs_env, dispersion_env=dispersion_env, error=error)
    ELSE IF (dft_control%qs_control%method_id == do_method_scptb) THEN
      ALLOCATE(dispersion_env,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ! set general defaults
      dispersion_env%doabc=.FALSE.
      dispersion_env%c9cnst=.FALSE.
      dispersion_env%lrc=.FALSE.
      dispersion_env%verbose=.FALSE.
      NULLIFY(dispersion_env%c6ab,dispersion_env%maxci,dispersion_env%r0ab,dispersion_env%rcov,&
              dispersion_env%r2r4,dispersion_env%cn,dispersion_env%cnkind,dispersion_env%cnlist)
      NULLIFY(dispersion_env%q_mesh,dispersion_env%kernel,dispersion_env%d2phi_dk2,&
              dispersion_env%d2y_dx2,dispersion_env%dftd_section)
      NULLIFY(dispersion_env%sab_vdw,dispersion_env%sab_cn)
      IF (scptb_control%dispersion) THEN
         dispersion_env%type = xc_vdw_fun_pairpot
         dispersion_env%pp_type = vdw_pairpot_dftd3
         dispersion_env%eps_cn = scptb_control%epscn
         dispersion_env%s6  = scptb_control%sd3(1)
         dispersion_env%sr6 = scptb_control%sd3(2)
         dispersion_env%s8  = scptb_control%sd3(3)
         dispersion_env%rc_disp = scptb_control%rcdisp
         dispersion_env%exp_pre = 0._dp
         dispersion_env%scaling = 0._dp
         dispersion_env%parameter_file_name = scptb_control%dispersion_parameter_file
         CALL qs_dispersion_pairpot_init(atomic_kind_set,qs_kind_set,dispersion_env,para_env=para_env,error=error)
      ELSE
         dispersion_env%type = xc_vdw_fun_none
      END IF
      CALL set_qs_env(qs_env, dispersion_env=dispersion_env, error=error)
    ELSE IF (dft_control%qs_control%method_id == do_method_dftb) THEN
      ALLOCATE(dispersion_env,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ! set general defaults
      dispersion_env%doabc=.FALSE.
      dispersion_env%c9cnst=.FALSE.
      dispersion_env%lrc=.FALSE.
      dispersion_env%verbose=.FALSE.
      NULLIFY(dispersion_env%c6ab,dispersion_env%maxci,dispersion_env%r0ab,dispersion_env%rcov,&
              dispersion_env%r2r4,dispersion_env%cn,dispersion_env%cnkind,dispersion_env%cnlist)
      NULLIFY(dispersion_env%q_mesh,dispersion_env%kernel,dispersion_env%d2phi_dk2,&
              dispersion_env%d2y_dx2,dispersion_env%dftd_section)
      NULLIFY(dispersion_env%sab_vdw,dispersion_env%sab_cn)
      IF (dftb_control%dispersion .AND. dftb_control%dispersion_type == dispersion_d3) THEN
         dispersion_env%type = xc_vdw_fun_pairpot
         dispersion_env%pp_type = vdw_pairpot_dftd3
         dispersion_env%eps_cn = dftb_control%epscn
         dispersion_env%s6  = dftb_control%sd3(1)
         dispersion_env%sr6 = dftb_control%sd3(2)
         dispersion_env%s8  = dftb_control%sd3(3)
         dispersion_env%rc_disp = dftb_control%rcdisp
         dispersion_env%exp_pre = 0._dp
         dispersion_env%scaling = 0._dp
         dispersion_env%parameter_file_name = dftb_control%dispersion_parameter_file
         CALL qs_dispersion_pairpot_init(atomic_kind_set,qs_kind_set,dispersion_env,para_env=para_env,error=error)
      ELSE
         dispersion_env%type = xc_vdw_fun_none
      END IF
      CALL set_qs_env(qs_env, dispersion_env=dispersion_env, error=error)
    ELSE IF(dft_control%qs_control%semi_empirical) THEN
      ALLOCATE(dispersion_env,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      ! set general defaults
      dispersion_env%doabc=.FALSE.
      dispersion_env%c9cnst=.FALSE.
      dispersion_env%lrc=.FALSE.
      dispersion_env%verbose=.FALSE.
      NULLIFY(dispersion_env%c6ab,dispersion_env%maxci,dispersion_env%r0ab,dispersion_env%rcov,&
              dispersion_env%r2r4,dispersion_env%cn,dispersion_env%cnkind,dispersion_env%cnlist)
      NULLIFY(dispersion_env%q_mesh,dispersion_env%kernel,dispersion_env%d2phi_dk2,&
              dispersion_env%d2y_dx2,dispersion_env%dftd_section)
      NULLIFY(dispersion_env%sab_vdw,dispersion_env%sab_cn)
      IF (se_control%dispersion) THEN
         dispersion_env%type = xc_vdw_fun_pairpot
         dispersion_env%pp_type = vdw_pairpot_dftd3
         dispersion_env%eps_cn = se_control%epscn
         dispersion_env%s6  = se_control%sd3(1)
         dispersion_env%sr6 = se_control%sd3(2)
         dispersion_env%s8  = se_control%sd3(3)
         dispersion_env%rc_disp = se_control%rcdisp
         dispersion_env%exp_pre = 0._dp
         dispersion_env%scaling = 0._dp
         dispersion_env%parameter_file_name = se_control%dispersion_parameter_file
         CALL qs_dispersion_pairpot_init(atomic_kind_set,qs_kind_set,dispersion_env,para_env=para_env,error=error)
      ELSE
         dispersion_env%type = xc_vdw_fun_none
      END IF
      CALL set_qs_env(qs_env, dispersion_env=dispersion_env, error=error)
    END IF

    !   *** Allocate the MO data types ***
    CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao, nelectron=nelectron)

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
            "nspins should be 1 or 2 for the time being ...",para_env)
    END IF

    IF ((MODULO(nelectron,2) /= 0).AND.(dft_control%nspins == 1)) THEN
       IF ( .NOT. dft_control%qs_control%ofgpw .AND.  .NOT. dft_control%smear) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Use the LSD option for an odd number of electrons",para_env)
       END IF
    END IF

    ! The transition potential method to calculate XAS needs LSD
    IF (dft_control%do_xas_calculation) THEN
       IF (dft_control%nspins == 1) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Use the LSD option for XAS with transition potential",para_env)
       END IF
       IF (dft_control%xas_control%xas_restart.AND.&
           scf_control%density_guess/= restart_guess) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Restart of XAS orbitals possible only if the GS orbitals"//&
               " are also read from restart file.",para_env)
       END IF
    END IF

    !   assigning the number of states per spin initial version, not yet very
    !   general. Should work for an even number of electrons and a single
    !   additional electron this set of options that requires full matrices,
    !   however, makes things a bit ugly right now.... we try to make a
    !   distinction between the number of electrons per spin and the number of
    !   MOs per spin this should allow the use of fractional occupations later
    !   on
    IF ( dft_control%qs_control%ofgpw ) THEN

       IF (dft_control%nspins == 1) THEN
          maxocc = nelectron
          nelectron_spin(1) = nelectron
          nelectron_spin(2) = 0
          n_mo(1) = 1
          n_mo(2) = 0
       ELSE
          IF (MODULO(nelectron + multiplicity - 1,2) /= 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "LSD: try to use a different multiplicity",para_env)
          END IF
          nelectron_spin(1) = (nelectron + multiplicity - 1)/2
          nelectron_spin(2) = (nelectron - multiplicity + 1)/2
          IF (nelectron_spin(2) < 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "LSD: too few electrons for this multiplicity",para_env)
          END IF
          maxocc=MAXVAL(nelectron_spin)
          n_mo(1) = MIN(nelectron_spin(1),1)
          n_mo(2) = MIN(nelectron_spin(2),1)
       END IF

    ELSE

       IF (dft_control%nspins == 1) THEN
          maxocc = 2.0_dp
          nelectron_spin(1) = nelectron
          nelectron_spin(2) = 0
          IF (MODULO(nelectron,2) ==0) THEN
            n_mo(1) = nelectron/2
          ELSE
            n_mo(1) = INT(nelectron/2._dp) + 1
          END IF
          n_mo(2) = 0
       ELSE
          maxocc=1.0_dp

          ! The simplist spin distribution is written here. Special cases will
          ! need additional user input
          IF (MODULO(nelectron + multiplicity - 1,2) /= 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "LSD: try to use a different multiplicity",para_env)
          END IF

          nelectron_spin(1) = (nelectron + multiplicity - 1)/2
          nelectron_spin(2) = (nelectron - multiplicity + 1)/2

          IF (nelectron_spin(2) < 0) THEN
             CALL stop_program(routineN,moduleN,__LINE__,&
                  "LSD: too few electrons for this multiplicity",para_env)
          END IF

          n_mo(1) = nelectron_spin(1)
          n_mo(2) = nelectron_spin(2)

       END IF

    END IF

    ! store the number of electrons once an for all
    CALL set_qs_env(qs_env,nelectron_total=nelectron,nelectron_spin=nelectron_spin,error=error)

    ! Check and set number of added (unoccupied) MOs
    CALL cp_assert((scf_control%added_mos(1) <= n_ao - n_mo(1)),cp_warning_level,&
                   cp_assertion_failed,routineP,&
                   "More added MOs requested than available. "//&
                   "The full set of unoccupied MOs will be used.",&
                   only_ionode=.TRUE.)
    scf_control%added_mos(1) = MIN(scf_control%added_mos(1),n_ao - n_mo(1))
    n_mo(1) = n_mo(1) + scf_control%added_mos(1)

    IF (dft_control%nspins == 2) THEN
      IF (scf_control%added_mos(2) > 0) THEN
        n_mo_add = scf_control%added_mos(2)
      ELSE
        n_mo_add = scf_control%added_mos(1)
      END IF
      CALL cp_assert((n_mo_add <= n_ao - n_mo(2)),cp_warning_level,&
                     cp_assertion_failed,routineP,&
                     "More added MOs requested for beta spin than available.",&
                     only_ionode=.TRUE.)
      scf_control%added_mos(2) = MIN(n_mo_add,n_ao - n_mo(2))
      n_mo(2) = n_mo(2) + scf_control%added_mos(2)
      CALL cp_assert((n_mo(2) <= n_mo(1)),cp_warning_level,&
                     cp_assertion_failed,routineP,&
                     "More beta than alpha MOs requested. "//&
                     "The number of beta MOs will be reduced to the number alpha MOs.",&
                     only_ionode=.TRUE.)
      n_mo(2) = MIN(n_mo(1),n_mo(2))
    END IF

    ! Compatibility checks for smearing

    IF (scf_control%smear%do_smear) THEN
      IF (scf_control%added_mos(1) == 0) THEN
        CALL stop_program(routineN,moduleN,__LINE__,&
                          "Extra MOs (ADDED_MOS) are required for smearing",para_env)
      END IF
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
            cp_assertion_failed,routineP,&
            "General ROKS scheme is not yet tested!",&
            only_ionode=.TRUE.)

       IF (scf_control%smear%do_smear) THEN
         CALL stop_program(routineN,moduleN,__LINE__,&
                           "The options ROKS and SMEAR are not compatible. "//&
                           "Try UKS instead of ROKS",para_env)
       END IF
    END IF


    ! in principle the restricted calculation could be performed
    ! using just one set of MOs and special casing most of the code
    ! right now we'll just take care of what is effectively an additional constraint
    ! at as few places as possible, just duplicating the beta orbitals
    IF (dft_control%restricted .AND. (output_unit>0)) THEN
       ! it is really not yet tested till the end ! Joost
       WRITE(output_unit,*) ""
       WRITE(output_unit,*) " **************************************"
       WRITE(output_unit,*) " restricted calculation cutting corners"
       WRITE(output_unit,*) " experimental feature, check code      "
       WRITE(output_unit,*) " **************************************"
    ENDIF

    ! no point in allocating these things here ?
    IF (qs_env%dft_control%qs_control%do_ls_scf .OR. &
        qs_env%dft_control%qs_control%do_almo_scf) THEN
      NULLIFY(mos)
    ELSE
      ALLOCATE (mos(dft_control%nspins),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DO ispin=1,dft_control%nspins
         NULLIFY (mos(ispin)%mo_set)
         CALL allocate_mo_set(mo_set=mos(ispin)%mo_set,&
                              nao=n_ao,&
                              nmo=n_mo(ispin),&
                              nelectron=nelectron_spin(ispin),&
                              n_el_f=REAL(nelectron_spin(ispin),dp),&
                              maxocc=maxocc,&
                              flexible_electron_count=dft_control%relax_multiplicity,&
                              error=error)
      END DO
    END IF

    CALL set_qs_env(qs_env,mos=mos,error=error)

    ! If we use auxiliary density matrix methods , set mo_set_aux_fit
    IF( dft_control%do_admm ) THEN
      CALL get_qs_kind_set(qs_kind_set, nelectron=nelectron, nsgf=n_ao_aux_fit, &
                           basis_set_id=use_aux_fit_basis_set)
      ALLOCATE (mos_aux_fit(dft_control%nspins),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      DO ispin=1,dft_control%nspins
         NULLIFY (mos_aux_fit(ispin)%mo_set)
         CALL allocate_mo_set(mo_set=mos_aux_fit(ispin)%mo_set,&
                              nao=n_ao_aux_fit,&
                              nmo=n_mo(ispin),&
                              nelectron=nelectron_spin(ispin),&
                              n_el_f=REAL(nelectron_spin(ispin),dp),&
                              maxocc=maxocc,&
                              flexible_electron_count=dft_control%relax_multiplicity,&
                              error=error)
      END DO
      CALL set_qs_env(qs_env,mos_aux_fit=mos_aux_fit,error=error)
    END IF


    ! Print the DFT control parameters
    CALL write_dft_control(dft_control,dft_section,error)

    ! Print the vdW control parameters
    IF (dft_control%qs_control%method_id == do_method_gpw .OR. &
        dft_control%qs_control%method_id == do_method_gapw .OR. &
        dft_control%qs_control%method_id == do_method_gapw_xc .OR. &
        dft_control%qs_control%method_id == do_method_lrigpw .OR. &
        dft_control%qs_control%method_id == do_method_scptb .OR. &
        dft_control%qs_control%method_id == do_method_dftb .OR. &
        dft_control%qs_control%method_id == do_method_ofgpw) THEN
      CALL get_qs_env(qs_env,dispersion_env=dispersion_env,error=error)
      CALL qs_write_dispersion(qs_env,dispersion_env,error=error)
    END IF

    ! Print the Quickstep control parameters
    CALL write_qs_control(dft_control%qs_control,dft_section,error)

    ! Print XES/XAS control parameters
    IF (dft_control%do_xas_calculation) THEN
       CALL cite_reference(Iannuzzi2007)
       CALL write_xas_control(dft_control%xas_control,dft_section,error=error)
    END IF

    ! Print the unnormalized basis set information (input data)
    CALL write_gto_basis_sets(qs_kind_set,subsys_section,error=error)

    ! Print the atomic kind set
    CALL write_qs_kind_set(qs_kind_set,subsys_section,error)

    ! Print the molecule kind set
    CALL write_molecule_kind_set(molecule_kind_set,subsys_section,error)

    ! Print the total number of kinds, atoms, basis functions etc.
    CALL write_total_numbers(atomic_kind_set, qs_kind_set,particle_set,qs_env%input,error)

    ! Print the atomic coordinates
    CALL write_qs_particle_coordinates(particle_set, qs_kind_set, subsys_section,label="QUICKSTEP",error=error)

    ! Print the interatomic distances
    CALL write_particle_distances(particle_set,cell,subsys_section,error)

    ! Print the requested structure data
    CALL write_structure_data(particle_set,cell,subsys_section,error)

    ! Print symmetry information
    CALL write_symmetry(particle_set,cell,subsys_section,error)

    ! Print the SCF parameters
    IF ((.NOT. qs_env%dft_control%qs_control%do_ls_scf) .AND. &
        (.NOT. qs_env%dft_control%qs_control%do_almo_scf) ) THEN
       CALL scf_c_write_parameters(scf_control,dft_section,error=error)
    ENDIF

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
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param force_env_section ...
!> \param error ...
!> \author Creation (06.10.2000)
! *****************************************************************************
  SUBROUTINE write_total_numbers(atomic_kind_set,qs_kind_set,particle_set,force_env_section,error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
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
       nkind = SIZE(qs_kind_set)

       CALL get_qs_kind_set(qs_kind_set,&
                            maxlgto=maxlgto,&
                            ncgf=ncgf,&
                            npgf=npgf,&
                            nset=nset,&
                            nsgf=nsgf,&
                            nshell=nshell,&
                            maxlppl=maxlppl,&
                            maxlppnl=maxlppnl)

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

       IF ((maxlppl > -1).AND.(maxlppnl > -1)) THEN
          WRITE (UNIT=output_unit,FMT="(/,T3,A,(T30,A,T75,I6))")&
               "Maximum angular momentum of the",&
               "- Orbital basis functions:                   ",maxlgto,&
               "- Local part of the GTH pseudopotential:     ",maxlppl,&
               "- Non-local part of the GTH pseudopotential: ",maxlppnl
       ELSEIF (maxlppl > -1) THEN
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

END MODULE qs_environment
