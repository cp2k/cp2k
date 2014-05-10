!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - mo_set_p_type added to qs_env (23.04.02,MK)
!>      - qs_force_type added to qs_env (05.06.02,MK)
!> \author MK (23.01.2002)
! *****************************************************************************
MODULE qs_environment_types
  USE admm_types,                      ONLY: admm_env_release,&
                                             admm_type
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release
  USE cp_control_types,                ONLY: dft_control_release,&
                                             dft_control_retain,&
                                             dft_control_type,&
                                             tddfpt_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_distribution_release,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_release_p
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_deallocate_matrix_set
  USE cp_ddapc_types,                  ONLY: cp_ddapc_ewald_release,&
                                             cp_ddapc_ewald_type,&
                                             cp_ddapc_release,&
                                             cp_ddapc_retain,&
                                             cp_ddapc_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_type
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_create,&
                                             cp_result_release,&
                                             cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_create,&
                                             cp_subsys_get,&
                                             cp_subsys_release,&
                                             cp_subsys_retain,&
                                             cp_subsys_type
  USE dbcsr_api,                       ONLY: dbcsr_distribution_obj
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE dm_ls_scf_types,                 ONLY: ls_scf_env_type,&
                                             ls_scf_release
  USE ep_qs_types,                     ONLY: ep_qs_release,&
                                             ep_qs_retain,&
                                             ep_qs_type
  USE et_coupling_types,               ONLY: et_coupling_release,&
                                             et_coupling_type
  USE ewald_environment_types,         ONLY: ewald_env_release,&
                                             ewald_env_retain,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_release,&
                                             ewald_pw_retain,&
                                             ewald_pw_type
  USE f77_blas
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_release,&
                                             fist_nonbond_env_type
  USE ga_environment_types,            ONLY: ga_env_release,&
                                             ga_env_retain,&
                                             ga_environment_type
  USE global_types,                    ONLY: global_environment_type
  USE harris_env_types,                ONLY: harris_env_release,&
                                             harris_env_type
  USE hartree_local_types,             ONLY: ecoul_1center_type,&
                                             get_hartree_local,&
                                             hartree_local_create,&
                                             hartree_local_release,&
                                             hartree_local_type,&
                                             set_hartree_local
  USE hfx_types,                       ONLY: hfx_release,&
                                             hfx_type
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type
  USE kg_environment_types,            ONLY: kg_env_release,&
                                             kg_environment_type
  USE kinds,                           ONLY: dp
  USE lri_environment_types,           ONLY: lri_density_release,&
                                             lri_density_type,&
                                             lri_env_release,&
                                             lri_environment_type
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE mp2_types,                       ONLY: mp2_env_release,&
                                             mp2_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_release,&
                                             pw_env_retain,&
                                             pw_env_type
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE qs_charges_types,                ONLY: qs_charges_release,&
                                             qs_charges_retain,&
                                             qs_charges_type
  USE qs_dftb_types,                   ONLY: qs_dftb_pairpot_release,&
                                             qs_dftb_pairpot_type
  USE qs_dispersion_types,             ONLY: qs_dispersion_release,&
                                             qs_dispersion_type
  USE qs_energy_types,                 ONLY: deallocate_qs_energy,&
                                             qs_energy_type
  USE qs_force_types,                  ONLY: deallocate_qs_force,&
                                             qs_force_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type,&
                                             qs_ks_qmmm_release,&
                                             qs_ks_qmmm_retain
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             qs_ks_release,&
                                             qs_ks_retain
  USE qs_linres_types,                 ONLY: linres_control_release,&
                                             linres_control_retain,&
                                             linres_control_type
  USE qs_local_rho_types,              ONLY: get_local_rho,&
                                             local_rho_set_create,&
                                             local_rho_set_release,&
                                             local_rho_type,&
                                             rhoz_type,&
                                             set_local_rho
  USE qs_matrix_pools,                 ONLY: mpools_get,&
                                             mpools_release,&
                                             mpools_retain,&
                                             qs_matrix_pools_type
  USE qs_mo_methods,                   ONLY: make_basis_cholesky,&
                                             make_basis_lowdin,&
                                             make_basis_simple,&
                                             make_basis_sm
  USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                             neighbor_list_set_p_type
  USE qs_oce_types,                    ONLY: deallocate_oce_set,&
                                             oce_matrix_type
  USE qs_period_efield_types,          ONLY: efield_berry_release,&
                                             efield_berry_type
  USE qs_rho0_types,                   ONLY: rho0_atom_type,&
                                             rho0_mpole_type
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_retain,&
                                             qs_rho_type
  USE qs_scf_types,                    ONLY: ot_method_nr,&
                                             qs_scf_env_type,&
                                             scf_env_release,&
                                             scf_env_retain
  USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                             wfi_release,&
                                             wfi_retain
  USE rel_control_types,               ONLY: rel_c_release,&
                                             rel_c_retain,&
                                             rel_control_type
  USE ri_environment_types,            ONLY: ri_env_release,&
                                             ri_environment_type
  USE rt_propagation_types,            ONLY: rt_prop_release,&
                                             rt_prop_type
  USE scf_control_types,               ONLY: scf_c_release,&
                                             scf_c_retain,&
                                             scf_control_type
  USE scp_environment_types,           ONLY: scp_env_release,&
                                             scp_env_retain,&
                                             scp_environment_type
  USE semi_empirical_mpole_types,      ONLY: nddo_mpole_release,&
                                             nddo_mpole_type
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_release,&
                                             semi_empirical_si_type
  USE semi_empirical_types,            ONLY: se_taper_release,&
                                             se_taper_type
  USE task_list_types,                 ONLY: deallocate_task_list,&
                                             task_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_release,&
                                             virial_retain,&
                                             virial_type
  USE wannier_states_types,            ONLY: wannier_centres_type
  USE xas_env_types,                   ONLY: xas_env_release,&
                                             xas_env_retain,&
                                             xas_environment_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_environment_types'
  INTEGER, PRIVATE, SAVE :: last_qs_env_id_nr=0

! *** Public data types ***

  PUBLIC :: qs_environment_p_type,&
            qs_environment_type

! *** Public subroutines ***

  PUBLIC :: get_qs_env,&
            init_qs_env,&
            qs_env_create,&
            qs_env_release,&
            qs_env_reorthogonalize_vectors,&
            qs_env_retain,&
            set_qs_env

! *****************************************************************************
!> \brief contains all the info needed by quickstep to do dft scf
!> \param id_nr identification number, guarateed to be unique
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param cell the simulation cell (box)
!> \param image_matrix matrix needed for the image charge approach within QM/MM
!> \param image_coeff coefficients of the image charges within QM/MM
!> \param ipiv auxiliary matrix needed for solving linear system of equations
!>        within the QM/MM image charge approach
!> \param calc_image_preconditioner LOGICAL which indicates if the preconditioner
!>        matrix is calculated when determining the coefficients of the image
!>        charges iteratively by CG
!> \param dft_control object that contains the values of various parameters
!>        relevant to a dft calculation
!> \param mos the molecular orbitals, in LSD array with the different spins
!> \param mos_aux_fit molecular orbitals for admm
!> \param sab_orb ,sac_ae, sac_ppl, sap_ppnl:
!>        neighbor lists for the calculation of the core Hamiltonian matrix
!> \param sap_oce: neighbor lists for the calculation of the expansion
!>        coefficients of the local atomic densities rho1_hard and rho1_soft
!> \param sab_se: neighbor lists for the calculation of the 2 centers
!>        hartree term in semi-empirical methods
!> \param sac_tbe: neighbor lists for the calculation of the tight binding
!>        Ewald terms (DFTB)
!> \param sab_core neighbor lists for the calculation of the core interactions
!> \param sab_all neighbor lists for the calculation of the  matrix element of
!>        non-symmetric operators
!> \param sab_vdw: neighbor lists for the calculation of dispersion interactions
!> \param sab_scp: neighbor lists for the calculation of the dispersion interactions
!>        when self-consistent polarization (SCP) is included
!> \param sab_almo: neighbor lists to create ALMO delocalization template
!> \param energy place where the energies calculated are stored
!> \param matrix_h core hamiltonian (matrix_h(1)), if needed by the calculation also
!>        its derivatives wrt. x,y, and z (matrix_h(2:4))
!> \param matrix_ks Khon-Sham matrix (filtred by the structure of S).
!>        with LSD an array with the differents orbitals
!> \param matrix_s the overlap matrix (matrix_s(1)), if needed by the calculation also
!>        its derivatives wrt. x,y, and z (matrix_s(2:4))
!> \param kinetic kinetic part of h
!> \param gamma_matrix gamma (Coulomb) matrix in DFTB, if needed by the calculation also
!>        its derivatives wrt. x,y, and z
!> \param rho the (old) density, in various representations (ao+grid)
!> \param rho_xc the (old) soft density, used for the GAPW_XC method only
!> \param scf_control various parameters relevant to the scf cycle
!> \param rel_control various parameters relevant to the relativistic calclation
!> \param pw_env various things to create and get info on the grids used
!>        by qs
!> \param ewald_env information to perform ewald summation (used by DFTB)
!> \param ewald_pw information on pw used to perform ewald summation (used by DFTB)
!> \param para_env parallel environement that should be  used to do the qs
!>        related calculations
!> \param blacs_env parallel environement needed by scalapack (full matrixes)
!> \param rho_core a grid with the charges of the cores of the atoms in the
!>        reciprocal (g) space
!> \param vppl a realspace grid with the local pseudopotential
!> \param vee a realspace grid with the external electrostatic potential
!>
!> \param rhoext a realspace grid with the external density 
!> \param external_vxc a realspace grid with the external vxc 
!> \param mask a realspace mask function to exclude pseuo region
!>
!> \param qs_charges an object where to put the charges that are present on
!>        the various grids
!> \param ks_env environement for the calculation of the Kohn Sham matrix
!> \param wf_history place where the snapshots of the previous states are
!>        stored (used during the md to extrapolate the new starting point)
!> \param scf_env place where temporary information used in the scf is stored
!> \param mpools all the pools of matrixes
!> \param distribution_ 2d: distribution of the atom pairs between the processors
!> \param subsys the particles, molecules,... of this environment
!> \param local_rho_set contains the atomic, compensations and core densities
!>                       and the local parts of the xc terms
!> \param hartree_local contains the 1, 2 and 3 centers coulomb terms
!> \param requires_mo_derivs logical, true if dE/dC is required (e.g. OT)
!> \param has_unit_metric logical, true if the S matrix is considered unity for the SCF
!> \param mo_derivs the actual derivatives of the total energy wrt to MO coeffs (divided by 2*f_i)
!> \param neighbor_list_id the current version of the neighbor_list
!> \param xas_env temporary information for xas calculation
!> \param dftb_potential pair potentials for use with DFTB
!> \param scp_env SCP environment for use with QS
!> \param dispersion_env environment for use with QS dispersion
!>
!>      compatibility get (things that you should get from the subsys):
!> \param atomic_kind_set array with infos about the species (atomic_kinds)
!>        present in the system
!> \param particle_set info on the atoms you simulate, pos,...
!> \param local_particles which particles ar local to this processor
!>      new:
!> \param local_molecules which molecules are local to this processor
!> \param molecule_kind_set description of the molecule kinds
!> \param molecule_set all the molecule description
!> \param use_harris LOGICAL which indicates, if the Harris energy functional
!>                    should be used.
!> \param harris_env The harris environment, which consist of the harris energy
!>                    and force type and the rho structure to pass the density of
!>                    the harris energy functional to the force calculation
!> \param ep_qs_env environment to perform an ep calculation
!> \param rtp all data needed for real time propagation
!> \param x contains data used in Hartree-Fock-Exchange calculations
!> \param task_list the list of tasks used in collocate and integrate
!> \param task_list_soft the list of tasks used in collocate and integrate in case of soft basis functions
!> \param mo_loc_history if a history of localized wfn is kept, they are stored here.
!> \par History
!>      11.2002 added doc and attribute description [fawzi]
!>      08.2004 renamed some of the very short names (s,c,k,h) for easier grepping
!> \author Matthias Krack & fawzi
! *****************************************************************************
  TYPE qs_environment_type
    INTEGER :: id_nr, ref_count
    LOGICAL :: use_ref_cell
    LOGICAL :: qmmm, qmmm_periodic
    LOGICAL :: requires_mo_derivs
    LOGICAL :: requires_matrix_vxc
    LOGICAL :: has_unit_metric
    LOGICAL :: use_harris
    LOGICAL :: run_rtp
    LOGICAL :: linres_run
    LOGICAL :: calc_image_preconditioner
    INTEGER :: nelectron_total
    INTEGER :: nelectron_spin(2)
    REAL(KIND=dp)                                         :: sim_time
    REAL(KIND=dp) :: start_time, target_time
    REAL(KIND=dp),DIMENSION(:,:),POINTER                  :: image_matrix
    REAL(KIND=dp),DIMENSION(:),POINTER                    :: image_coeff
    INTEGER,DIMENSION(:),POINTER                          :: ipiv
    INTEGER                                               :: sim_step
    TYPE(ga_environment_type), POINTER                    :: ga_env
    TYPE(ls_scf_env_type), POINTER                        :: ls_scf_env
    TYPE(cell_type), POINTER                              :: cell, super_cell
    TYPE(cell_type), POINTER                              :: cell_ref
    TYPE(dft_control_type), POINTER                       :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), POINTER            :: mos, mos_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER             :: mo_derivs_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER             :: mo_loc_history
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs
    INTEGER                                               :: neighbor_list_id
    TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb,sac_ae,&
                                                             sac_ppl,sap_ppnl,&
                                                             sap_oce,sab_se,sab_vdw,&
                                                             sab_tbe,sab_core,sab_all,sab_scp,&
                                                             sab_aux_fit, sab_aux_fit_vs_orb,&
                                                             sab_aux_fit_asymm, sab_lrc,&
                                                             sab_almo
    TYPE(qs_energy_type), POINTER                         :: energy
    TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
    TYPE(harris_env_type), POINTER                        :: harris_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_h,&
                                                             matrix_ks,&
                                                             matrix_vxc,&
                                                             matrix_s,&
                                                             kinetic,&
                                                             matrix_w,&
                                                             matrix_w_mp2,&
                                                             matrix_p_mp2,&
                                                             gamma_matrix,&
                                                             matrix_s_aux_fit,&
                                                             matrix_s_aux_fit_vs_orb,&
                                                             matrix_ks_aux_fit
    TYPE(qs_rho_type), POINTER                            :: rho,&
                                                             rho_buffer,&
                                                             rho_xc,&
                                                             rho_aux_fit
    TYPE(scf_control_type), POINTER                       :: scf_control
    TYPE(rel_control_type), POINTER                       :: rel_control
    TYPE(scp_environment_type), POINTER                   :: scp_env
    TYPE(pw_env_type), POINTER                            :: pw_env
    TYPE(cp_para_env_type), POINTER                       :: para_env
    TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
    TYPE(pw_p_type), POINTER                              :: rho_core
    TYPE(pw_p_type), POINTER                              :: vppl
    TYPE(pw_p_type), POINTER                              :: rho_nlcc
    TYPE(pw_p_type), POINTER                              :: rho_nlcc_g
    TYPE(pw_p_type), POINTER                              :: vee
! ZMP adding variables
    TYPE(qs_rho_type), POINTER                            :: rho_external
    TYPE(pw_p_type), POINTER                              :: external_vxc
    TYPE(pw_p_type), POINTER                              :: mask
    TYPE(task_list_type), POINTER                         :: task_list_zmp
    TYPE(qs_charges_type), POINTER                        :: qs_charges
    TYPE(qs_ks_env_type), POINTER                         :: ks_env
    TYPE(qs_ks_qmmm_env_type), POINTER                    :: ks_qmmm_env
    TYPE(qmmm_env_qm_type),POINTER                        :: qmmm_env_qm
    TYPE(qs_wf_history_type), POINTER                     :: wf_history
    TYPE(qs_scf_env_type), POINTER                        :: scf_env
    TYPE(qs_matrix_pools_type), POINTER                   :: mpools
    TYPE(qs_matrix_pools_type), POINTER                   :: mpools_aux_fit
    TYPE(distribution_2d_type), POINTER                   :: distribution_2d
    TYPE(dbcsr_distribution_obj), POINTER                 :: dbcsr_dist
    TYPE(cp_subsys_type), POINTER                         :: subsys
    TYPE(oce_matrix_type), POINTER                        :: oce
    TYPE(local_rho_type), POINTER                         :: local_rho_set
    TYPE(hartree_local_type),  POINTER                    :: hartree_local
    TYPE(section_vals_type), POINTER                      :: input
    TYPE(linres_control_type), POINTER                    :: linres_control
    TYPE(xas_environment_type), POINTER                   :: xas_env
    TYPE(cp_ddapc_type), POINTER                          :: cp_ddapc_env
    TYPE(cp_ddapc_ewald_type), POINTER                    :: cp_ddapc_ewald
    TYPE(virial_type), POINTER                            :: virial
    REAL(KIND = dp), DIMENSION(:,:), POINTER              :: outer_scf_history
    INTEGER                                               :: outer_scf_ihistory
    TYPE(ep_qs_type), POINTER                             :: ep_qs_env
    TYPE(hfx_type),DIMENSION(:,:), POINTER                :: x_data
    TYPE(et_coupling_type),POINTER                        :: et_coupling
    TYPE(qs_dftb_pairpot_type), DIMENSION(:,:), POINTER   :: dftb_potential
    TYPE(cp_result_type) ,POINTER                         :: results
    TYPE(task_list_type), POINTER                         :: task_list
    TYPE(task_list_type), POINTER                         :: task_list_aux_fit
    TYPE(task_list_type), POINTER                         :: task_list_soft
    TYPE(admm_type), POINTER                              :: admm_env
    ! LRI
    TYPE(lri_environment_type), POINTER                   :: lri_env
    TYPE(lri_density_type), POINTER                       :: lri_density
    ! HFX RI
    TYPE(ri_environment_type), POINTER                    :: hfx_ri_env
    ! Empirical dispersion
    TYPE(qs_dispersion_type), POINTER                     :: dispersion_env
    ! Semi-empirical and DFTB types
    TYPE(ewald_environment_type),POINTER                  :: ewald_env
    TYPE(ewald_pw_type),POINTER                           :: ewald_pw
    ! Semi-empirical types
    TYPE(se_taper_type), POINTER                          :: se_taper
    TYPE(semi_empirical_si_type),POINTER                  :: se_store_int_env
    TYPE(nddo_mpole_type), POINTER                        :: se_nddo_mpole
    TYPE(fist_nonbond_env_type), POINTER                  :: se_nonbond_env
    TYPE(rt_prop_type),POINTER                            :: rtp
    TYPE(efield_berry_type),POINTER                       :: efield
    ! a history for the broyden ot
    REAL(KIND = dp)                                       :: broyden_adaptive_sigma
    TYPE(mp2_type), POINTER                               :: mp2_env
    TYPE(kg_environment_type), POINTER                    :: kg_env
    TYPE(wannier_centres_type), POINTER,DIMENSION(:)      :: WannierCentres=>NULL()
    ! atomic properties
    TYPE(atprop_type),POINTER                             :: atprop
  END TYPE qs_environment_type

! *****************************************************************************
!> \brief to build arrays of pointers
!> \param qs_env the pointer to the qs_env
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE qs_environment_p_type
     TYPE(qs_environment_type), POINTER :: qs_env
  END TYPE qs_environment_p_type

CONTAINS

! *****************************************************************************
!> \brief   Get the QUICKSTEP environment.
!> \author  MK
!> \date    23.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_qs_env(qs_env,atomic_kind_set,cell,super_cell,cell_ref,use_ref_cell,&
       dft_control,mos,mos_aux_fit,sab_orb,sab_aux_fit,sab_aux_fit_asymm, sab_aux_fit_vs_orb,&
       sab_all,qmmm,qmmm_periodic,sac_ae,sac_ppl,sap_ppnl,sab_vdw,sab_scp,sap_oce,sab_lrc,&
       sab_se,sab_tbe,sab_core,sab_almo,particle_set,energy,force,harris_env,&
       matrix_h,matrix_ks,matrix_vxc,run_rtp,rtp,&
       matrix_ks_aux_fit,matrix_s,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb,matrix_w,matrix_w_mp2,matrix_p_mp2,&
       gamma_matrix, rho, rho_aux_fit,rho_buffer, rho_xc, pw_env, ewald_env, ewald_pw, &
       mpools,mpools_aux_fit,input, para_env, blacs_env, scf_control,rel_control,kinetic,qs_charges,&
       vppl, rho_core, rho_nlcc, rho_nlcc_g, ks_env, ks_qmmm_env, wf_history,scf_env,use_harris, id_nr, local_particles,&
       local_molecules, distribution_2d, dbcsr_dist, molecule_kind_set,molecule_set,subsys,oce,rho_atom_set,&
       task_list,task_list_zmp,task_list_aux_fit,task_list_soft,rho0_atom_set,rho0_mpole,rhoz_set,ecoul_1c,&
       rho0_s_rs,rho0_s_gs,has_unit_metric,requires_mo_derivs,mo_derivs,&
       mo_derivs_aux_fit, mo_loc_history, tddfpt_control,nkind,natom,nelectron_total,nelectron_spin,efield,&
       neighbor_list_id,linres_control,xas_env,virial,cp_ddapc_env,cp_ddapc_ewald,&
       outer_scf_history,outer_scf_ihistory,ep_qs_env,x_data,et_coupling,dftb_potential,results,&
       scp_env,se_taper,se_store_int_env,se_nddo_mpole,se_nonbond_env, admm_env, &
       lri_env,lri_density,hfx_ri_env,dispersion_env, vee,rho_external,external_vxc,mask,&
       mp2_env,kg_env,WannierCentres,ga_env,atprop,ls_scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, super_cell, cell_ref
    LOGICAL, OPTIONAL                        :: use_ref_cell
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mos, mos_aux_fit
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: sab_orb, sab_aux_fit, &
                                                sab_aux_fit_asymm, &
                                                sab_aux_fit_vs_orb, sab_all
    LOGICAL, OPTIONAL                        :: qmmm, qmmm_periodic
    TYPE(neighbor_list_set_p_type), DIMENSION(:), OPTIONAL, POINTER :: &
      sac_ae, sac_ppl, sap_ppnl, sab_vdw, sab_scp, sap_oce, sab_lrc, sab_se, &
      sab_tbe, sab_core, sab_almo
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(harris_env_type), OPTIONAL, POINTER :: harris_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: matrix_h, matrix_ks, &
                                                matrix_vxc
    LOGICAL, OPTIONAL                        :: run_rtp
    TYPE(rt_prop_type), OPTIONAL, POINTER    :: rtp
    TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, POINTER :: &
      matrix_ks_aux_fit, matrix_s, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, &
      matrix_w, matrix_w_mp2, matrix_p_mp2, gamma_matrix
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho, rho_aux_fit, rho_buffer, &
                                                rho_xc
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(ewald_pw_type), OPTIONAL, POINTER   :: ewald_pw
    TYPE(qs_matrix_pools_type), OPTIONAL, &
      POINTER                                :: mpools, mpools_aux_fit
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                :: scf_control
    TYPE(rel_control_type), OPTIONAL, &
      POINTER                                :: rel_control
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: kinetic
    TYPE(qs_charges_type), OPTIONAL, POINTER :: qs_charges
    TYPE(pw_p_type), OPTIONAL, POINTER       :: vppl, rho_core, rho_nlcc, &
                                                rho_nlcc_g
    TYPE(qs_ks_env_type), OPTIONAL, POINTER  :: ks_env
    TYPE(qs_ks_qmmm_env_type), OPTIONAL, &
      POINTER                                :: ks_qmmm_env
    TYPE(qs_wf_history_type), OPTIONAL, &
      POINTER                                :: wf_history
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    LOGICAL, OPTIONAL                        :: use_harris
    INTEGER, INTENT(out), OPTIONAL           :: id_nr
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles, &
                                                local_molecules
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(dbcsr_distribution_obj), OPTIONAL, &
      POINTER                                :: dbcsr_dist
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(oce_matrix_type), OPTIONAL, POINTER :: oce
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_atom_set
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list, task_list_zmp, &
                                                task_list_aux_fit, &
                                                task_list_soft
    TYPE(rho0_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho0_atom_set
    TYPE(rho0_mpole_type), OPTIONAL, POINTER :: rho0_mpole
    TYPE(rhoz_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rhoz_set
    TYPE(ecoul_1center_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ecoul_1c
    TYPE(pw_p_type), OPTIONAL, POINTER       :: rho0_s_rs, rho0_s_gs
    LOGICAL, OPTIONAL                        :: has_unit_metric, &
                                                requires_mo_derivs
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mo_derivs
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mo_derivs_aux_fit, &
                                                mo_loc_history
    TYPE(tddfpt_control_type), OPTIONAL, &
      POINTER                                :: tddfpt_control
    INTEGER, OPTIONAL                        :: nkind, natom, nelectron_total
    INTEGER, DIMENSION(2), OPTIONAL          :: nelectron_spin
    TYPE(efield_berry_type), OPTIONAL, &
      POINTER                                :: efield
    INTEGER, OPTIONAL                        :: neighbor_list_id
    TYPE(linres_control_type), OPTIONAL, &
      POINTER                                :: linres_control
    TYPE(xas_environment_type), OPTIONAL, &
      POINTER                                :: xas_env
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_ddapc_type), OPTIONAL, POINTER   :: cp_ddapc_env
    TYPE(cp_ddapc_ewald_type), OPTIONAL, &
      POINTER                                :: cp_ddapc_ewald
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: outer_scf_history
    INTEGER, INTENT(out), OPTIONAL           :: outer_scf_ihistory
    TYPE(ep_qs_type), OPTIONAL, POINTER      :: ep_qs_env
    TYPE(hfx_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: x_data
    TYPE(et_coupling_type), OPTIONAL, &
      POINTER                                :: et_coupling
    TYPE(qs_dftb_pairpot_type), &
      DIMENSION(:, :), OPTIONAL, POINTER     :: dftb_potential
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(scp_environment_type), OPTIONAL, &
      POINTER                                :: scp_env
    TYPE(se_taper_type), OPTIONAL, POINTER   :: se_taper
    TYPE(semi_empirical_si_type), OPTIONAL, &
      POINTER                                :: se_store_int_env
    TYPE(nddo_mpole_type), OPTIONAL, POINTER :: se_nddo_mpole
    TYPE(fist_nonbond_env_type), OPTIONAL, &
      POINTER                                :: se_nonbond_env
    TYPE(admm_type), OPTIONAL, POINTER       :: admm_env
    TYPE(lri_environment_type), OPTIONAL, &
      POINTER                                :: lri_env
    TYPE(lri_density_type), OPTIONAL, &
      POINTER                                :: lri_density
    TYPE(ri_environment_type), OPTIONAL, &
      POINTER                                :: hfx_ri_env
    TYPE(qs_dispersion_type), OPTIONAL, &
      POINTER                                :: dispersion_env
    TYPE(pw_p_type), OPTIONAL, POINTER       :: vee
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho_external
    TYPE(pw_p_type), OPTIONAL, POINTER       :: external_vxc, mask
    TYPE(mp2_type), OPTIONAL, POINTER        :: mp2_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(wannier_centres_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: WannierCentres
    TYPE(ga_environment_type), OPTIONAL, &
      POINTER                                :: ga_env
    TYPE(atprop_type), OPTIONAL, POINTER     :: atprop
    TYPE(ls_scf_env_type), OPTIONAL, POINTER :: ls_scf_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_qs_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(mol_kind_new_list_type), POINTER    :: molecule_kinds_new
    TYPE(mol_new_list_type), POINTER         :: molecules_new
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(rho0_mpole_type), POINTER           :: rho0_m

    failure=.FALSE.
    NULLIFY(atomic_kinds, particles, molecules_new,&
         molecule_kinds_new)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env%subsys),cp_failure_level,routineP,error,failure)
    IF (PRESENT(nelectron_total)) nelectron_total=qs_env%nelectron_total
    IF (PRESENT(nelectron_spin))  nelectron_spin=qs_env%nelectron_spin
    IF (PRESENT(outer_scf_history)) outer_scf_history=>qs_env%outer_scf_history
    IF (PRESENT(outer_scf_ihistory)) outer_scf_ihistory=qs_env%outer_scf_ihistory
    IF (PRESENT(ga_env)) ga_env => qs_env%ga_env
    IF (PRESENT(cell)) cell => qs_env%cell
    IF (PRESENT(mp2_env)) mp2_env => qs_env%mp2_env
    IF (PRESENT(kg_env)) kg_env => qs_env%kg_env   
    IF (PRESENT(super_cell)) super_cell => qs_env%super_cell
    IF (PRESENT(cell_ref)) cell_ref => qs_env%cell_ref
    IF (PRESENT(use_ref_cell)) use_ref_cell = qs_env%use_ref_cell
    IF (PRESENT(qmmm)) qmmm = qs_env%qmmm
    IF (PRESENT(qmmm_periodic)) qmmm_periodic = qs_env%qmmm_periodic
    IF (PRESENT(dft_control)) dft_control => qs_env%dft_control
    IF (PRESENT(mos)) mos => qs_env%mos
    IF (PRESENT(mos_aux_fit)) mos_aux_fit => qs_env%mos_aux_fit
    IF (PRESENT(neighbor_list_id)) neighbor_list_id=qs_env%neighbor_list_id
    IF (PRESENT(sab_orb)) sab_orb => qs_env%sab_orb
    IF (PRESENT(sab_aux_fit)) sab_aux_fit => qs_env%sab_aux_fit
    IF (PRESENT(sab_aux_fit_asymm)) sab_aux_fit_asymm => qs_env%sab_aux_fit_asymm
    IF (PRESENT(sab_aux_fit_vs_orb)) sab_aux_fit_vs_orb => qs_env%sab_aux_fit_vs_orb
    IF (PRESENT(sab_all)) sab_all => qs_env%sab_all
    IF (PRESENT(sab_scp)) sab_scp => qs_env%sab_scp
    IF (PRESENT(sab_vdw)) sab_vdw => qs_env%sab_vdw
    IF (PRESENT(sac_ae)) sac_ae => qs_env%sac_ae
    IF (PRESENT(sac_ppl)) sac_ppl => qs_env%sac_ppl
    IF (PRESENT(sap_ppnl)) sap_ppnl => qs_env%sap_ppnl
    IF (PRESENT(sap_oce)) sap_oce =>  qs_env%sap_oce
    IF (PRESENT(sab_se)) sab_se =>  qs_env%sab_se
    IF (PRESENT(sab_lrc)) sab_lrc =>  qs_env%sab_lrc
    IF (PRESENT(sab_tbe)) sab_tbe =>  qs_env%sab_tbe
    IF (PRESENT(sab_core)) sab_core =>  qs_env%sab_core
    IF (PRESENT(sab_almo)) sab_almo =>  qs_env%sab_almo
    IF (PRESENT(energy)) energy => qs_env%energy
    IF (PRESENT(force)) force => qs_env%force
    IF (PRESENT(use_harris)) use_harris = qs_env%use_harris
    IF (PRESENT(harris_env)) harris_env => qs_env%harris_env
    IF (PRESENT(matrix_h)) matrix_h => qs_env%matrix_h
    IF (PRESENT(matrix_ks)) matrix_ks => qs_env%matrix_ks
    IF (PRESENT(matrix_vxc)) matrix_vxc => qs_env%matrix_vxc
    IF (PRESENT(matrix_ks_aux_fit)) matrix_ks_aux_fit => qs_env%matrix_ks_aux_fit
    IF (PRESENT(kinetic)) kinetic => qs_env%kinetic
    IF (PRESENT(matrix_s)) matrix_s => qs_env%matrix_s
    IF (PRESENT(matrix_s_aux_fit)) matrix_s_aux_fit => qs_env%matrix_s_aux_fit
    IF (PRESENT(matrix_s_aux_fit_vs_orb)) matrix_s_aux_fit_vs_orb => qs_env%matrix_s_aux_fit_vs_orb
    IF (PRESENT(matrix_w)) matrix_w => qs_env%matrix_w
    IF (PRESENT(matrix_w_mp2)) matrix_w_mp2 => qs_env%matrix_w_mp2
    IF (PRESENT(matrix_p_mp2)) matrix_p_mp2 => qs_env%matrix_p_mp2
    IF (PRESENT(gamma_matrix)) gamma_matrix => qs_env%gamma_matrix
    IF (PRESENT(rho)) rho => qs_env%rho
    IF (PRESENT(rho_aux_fit)) rho_aux_fit => qs_env%rho_aux_fit
    IF (PRESENT(rho_buffer)) rho_buffer => qs_env%rho_buffer
    IF (PRESENT(rho_xc)) rho_xc => qs_env%rho_xc
    IF (PRESENT(pw_env)) pw_env => qs_env%pw_env
    IF (PRESENT(ewald_env)) ewald_env => qs_env%ewald_env
    IF (PRESENT(ewald_pw)) ewald_pw => qs_env%ewald_pw
    IF (PRESENT(mpools)) mpools => qs_env%mpools
    IF (PRESENT(mpools_aux_fit)) mpools_aux_fit => qs_env%mpools_aux_fit
    IF (PRESENT(blacs_env)) blacs_env => qs_env%blacs_env
    IF (PRESENT(scf_control)) scf_control => qs_env%scf_control
    IF (PRESENT(rel_control)) rel_control => qs_env%rel_control
    IF (PRESENT(rho_core)) rho_core => qs_env%rho_core
    IF (PRESENT(rho_nlcc)) rho_nlcc => qs_env%rho_nlcc
    IF (PRESENT(rho_nlcc_g)) rho_nlcc_g => qs_env%rho_nlcc_g
    IF (PRESENT(vppl)) vppl => qs_env%vppl
    IF (PRESENT(vee)) vee => qs_env%vee
! ZMP pointing vectors
    IF (PRESENT(rho_external)) rho_external => qs_env%rho_external
    IF (PRESENT(external_vxc)) external_vxc => qs_env%external_vxc
    IF (PRESENT(mask)) mask => qs_env%mask

    IF (PRESENT(qs_charges)) qs_charges => qs_env%qs_charges
    IF (PRESENT(ks_env)) ks_env => qs_env%ks_env
    IF (PRESENT(ks_qmmm_env)) ks_qmmm_env => qs_env%ks_qmmm_env
    IF (PRESENT(wf_history)) wf_history => qs_env%wf_history
    IF (PRESENT(scf_env)) scf_env => qs_env%scf_env
    IF (PRESENT(id_nr)) id_nr=qs_env%id_nr
    IF (PRESENT(subsys)) subsys => qs_env%subsys
    IF (PRESENT(para_env)) para_env => qs_env%para_env
    IF (PRESENT(distribution_2d)) distribution_2d => qs_env%distribution_2d
    IF (PRESENT(dbcsr_dist)) dbcsr_dist => qs_env%dbcsr_dist
    IF (PRESENT(oce)) oce =>  qs_env%oce
    IF (PRESENT(requires_mo_derivs)) requires_mo_derivs = qs_env%requires_mo_derivs
    IF (PRESENT(has_unit_metric)) has_unit_metric = qs_env%has_unit_metric
    IF (PRESENT(mo_derivs)) mo_derivs => qs_env%mo_derivs
    IF (PRESENT(mo_derivs_aux_fit)) mo_derivs_aux_fit => qs_env%mo_derivs_aux_fit
    IF (PRESENT(mo_loc_history)) mo_loc_history => qs_env%mo_loc_history
    IF (PRESENT(tddfpt_control)) tddfpt_control => qs_env%dft_control%tddfpt_control
    IF (PRESENT(linres_control)) linres_control => qs_env%linres_control
    IF (PRESENT(virial)) virial => qs_env%virial
    IF (PRESENT(se_taper)) se_taper => qs_env%se_taper
    IF (PRESENT(se_store_int_env)) se_store_int_env => qs_env%se_store_int_env
    IF (PRESENT(se_nddo_mpole)) se_nddo_mpole => qs_env%se_nddo_mpole
    IF (PRESENT(se_nonbond_env)) se_nonbond_env => qs_env%se_nonbond_env
    IF (PRESENT(admm_env)) admm_env => qs_env%admm_env
    IF (PRESENT(lri_env)) lri_env => qs_env%lri_env
    IF (PRESENT(lri_density)) lri_density => qs_env%lri_density
    IF (PRESENT(hfx_ri_env)) hfx_ri_env => qs_env%hfx_ri_env
    IF (PRESENT(dispersion_env)) dispersion_env => qs_env%dispersion_env
    IF (PRESENT(run_rtp)) run_rtp=qs_env%run_rtp
    IF (PRESENT(rtp)) rtp=>qs_env%rtp
    IF (PRESENT(ls_scf_env)) ls_scf_env=>qs_env%ls_scf_env

    IF (PRESENT(rho_atom_set)) &
       CALL get_local_rho(qs_env%local_rho_set, rho_atom_set=rho_atom_set)
    IF (PRESENT(rho0_atom_set)) &
       CALL get_local_rho(qs_env%local_rho_set, rho0_atom_set=rho0_atom_set)
    IF (PRESENT(rho0_mpole)) &
         CALL get_local_rho(qs_env%local_rho_set,  rho0_mpole=rho0_mpole)
    IF (PRESENT(rhoz_set)) &
         CALL get_local_rho(qs_env%local_rho_set, rhoz_set=rhoz_set)
    IF (PRESENT(ecoul_1c)) &
         CALL get_hartree_local(qs_env%hartree_local,ecoul_1c=ecoul_1c)
    IF (PRESENT(rho0_s_rs)) THEN
       CALL get_local_rho(qs_env%local_rho_set,  rho0_mpole=rho0_m)
       IF(ASSOCIATED(rho0_m)) THEN
         rho0_s_rs => rho0_m%rho0_s_rs
       END IF
    END IF
    IF (PRESENT(rho0_s_gs)) THEN
       CALL get_local_rho(qs_env%local_rho_set,  rho0_mpole=rho0_m)
       IF(ASSOCIATED(rho0_m)) THEN
         rho0_s_gs => rho0_m%rho0_s_gs
       END IF
    END IF

    IF(PRESENT(xas_env)) xas_env => qs_env%xas_env

    CALL cp_subsys_get(qs_env%subsys,&
                      atomic_kinds=atomic_kinds,&
                      local_molecules_new=local_molecules,&
                      local_particles=local_particles,&
                      particles=particles,&
                      molecule_kinds_new=molecule_kinds_new,&
                      molecules_new=molecules_new,&
                      error=error)

    IF (PRESENT(natom)) natom=SIZE(particles%els)
    IF (PRESENT(nkind)) nkind=SIZE(atomic_kinds%els)
    IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
    IF (PRESENT(particle_set)) particle_set => particles%els
    IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds_new%els
    IF (PRESENT(molecule_set)) molecule_set => molecules_new%els
    IF (PRESENT(input)) input => qs_env%input
    IF (PRESENT(cp_ddapc_env)) cp_ddapc_env => qs_env%cp_ddapc_env
    IF (PRESENT(cp_ddapc_ewald)) cp_ddapc_ewald => qs_env%cp_ddapc_ewald
    IF (PRESENT(ep_qs_env)) ep_qs_env => qs_env%ep_qs_env
    IF (PRESENT(x_data)) x_data => qs_env%x_data
    IF (PRESENT(et_coupling)) et_coupling => qs_env%et_coupling
    IF (PRESENT(dftb_potential)) dftb_potential => qs_env%dftb_potential
    IF (PRESENT(scp_env)) scp_env => qs_env%scp_env
    IF (PRESENT(efield)) efield => qs_env%efield

    IF (PRESENT(atprop)) atprop => qs_env%atprop

    IF (PRESENT(results)) results => qs_env%results
    IF (PRESENT(task_list)) task_list => qs_env%task_list
    IF (PRESENT(task_list_aux_fit)) task_list_aux_fit => qs_env%task_list_aux_fit
    IF (PRESENT(task_list_soft)) task_list_soft => qs_env%task_list_soft

    IF (PRESENT(WannierCentres)) WannierCentres => qs_env%WannierCentres

  END SUBROUTINE get_qs_env

! *****************************************************************************
!> \brief  Initialise the QUICKSTEP environment.
!> \author  MK
!> \date    25.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_qs_env(qs_env,para_env,globenv,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    NULLIFY (qs_env%blacs_env)
    NULLIFY (qs_env%ls_scf_env)
    NULLIFY (qs_env%image_matrix)
    NULLIFY (qs_env%ipiv)
    NULLIFY (qs_env%image_coeff)
    NULLIFY (qs_env%ga_env)
    NULLIFY (qs_env%cell)
    NULLIFY (qs_env%super_cell)
    NULLIFY (qs_env%cell_ref)
    NULLIFY (qs_env%dft_control)
    NULLIFY (qs_env%mos)
    NULLIFY (qs_env%mos_aux_fit)
    NULLIFY (qs_env%sab_orb)
    NULLIFY (qs_env%sab_aux_fit)
    NULLIFY (qs_env%sab_aux_fit_asymm)
    NULLIFY (qs_env%sab_aux_fit_vs_orb)
    NULLIFY (qs_env%sab_scp)
    NULLIFY (qs_env%sab_vdw)
    NULLIFY (qs_env%sab_all)
    NULLIFY (qs_env%sac_ae)
    NULLIFY (qs_env%sac_ppl)
    NULLIFY (qs_env%sap_ppnl)
    NULLIFY (qs_env%sap_oce)
    NULLIFY (qs_env%sab_se)
    NULLIFY (qs_env%sab_lrc)
    NULLIFY (qs_env%sab_tbe)
    NULLIFY (qs_env%sab_core)
    NULLIFY (qs_env%sab_almo)
    NULLIFY (qs_env%energy)
    NULLIFY (qs_env%force)
    NULLIFY (qs_env%harris_env)
    NULLIFY (qs_env%matrix_h)
    NULLIFY (qs_env%matrix_ks)
    NULLIFY (qs_env%matrix_vxc)
    NULLIFY (qs_env%matrix_ks_aux_fit)
    NULLIFY (qs_env%kinetic)
    NULLIFY (qs_env%matrix_s)
    NULLIFY (qs_env%matrix_s_aux_fit)
    NULLIFY (qs_env%matrix_s_aux_fit_vs_orb)
    NULLIFY (qs_env%matrix_w)
    NULLIFY (qs_env%matrix_w_mp2)
    NULLIFY (qs_env%matrix_p_mp2)
    NULLIFY (qs_env%gamma_matrix)
    NULLIFY (qs_env%mpools)
    NULLIFY (qs_env%mpools_aux_fit)
    NULLIFY (qs_env%rho)
    NULLIFY (qs_env%rho_aux_fit)
    NULLIFY (qs_env%rho_buffer)
    NULLIFY (qs_env%rho_xc)
    NULLIFY (qs_env%pw_env)
    NULLIFY (qs_env%ewald_env)
    NULLIFY (qs_env%ewald_pw)
    NULLIFY (qs_env%scf_control)
    NULLIFY (qs_env%rel_control)
    NULLIFY (qs_env%qs_charges)
    NULLIFY (qs_env%rho_nlcc)
    NULLIFY (qs_env%rho_nlcc_g)
    NULLIFY (qs_env%rho_core)
    NULLIFY (qs_env%vppl)
    NULLIFY (qs_env%vee)
! ZMP initializing arrays
    NULLIFY (qs_env%rho_external)
    NULLIFY (qs_env%external_vxc)
    NULLIFY (qs_env%mask)

    NULLIFY (qs_env%ks_env)
    NULLIFY (qs_env%ks_qmmm_env)
    NULLIFY (qs_env%wf_history)
    NULLIFY (qs_env%scf_env)
    NULLIFY (qs_env%distribution_2d)
    NULLIFY (qs_env%dbcsr_dist)
    NULLIFY (qs_env%subsys)
    NULLIFY (qs_env%blacs_env)
    NULLIFY (qs_env%oce)
    NULLIFY (qs_env%local_rho_set)
    NULLIFY (qs_env%hartree_local)
    NULLIFY (qs_env%input)
    NULLIFY (qs_env%linres_control)
    NULLIFY (qs_env%xas_env)
    NULLIFY (qs_env%virial)
    NULLIFY (qs_env%cp_ddapc_env)
    NULLIFY (qs_env%cp_ddapc_ewald)
    NULLIFY (qs_env%outer_scf_history)
    NULLIFY (qs_env%ep_qs_env)
    NULLIFY (qs_env%x_data)
    NULLIFY (qs_env%et_coupling)
    NULLIFY (qs_env%dftb_potential)
    NULLIFY (qs_env%results)
    NULLIFY (qs_env%task_list)
    NULLIFY (qs_env%task_list_aux_fit)
    NULLIFY (qs_env%task_list_soft)
! ZMP zmp task list
    NULLIFY (qs_env%task_list_zmp)

    NULLIFY (qs_env%scp_env)
    NULLIFY (qs_env%se_taper)
    NULLIFY (qs_env%se_store_int_env)
    NULLIFY (qs_env%se_nddo_mpole)
    NULLIFY (qs_env%se_nonbond_env)
    NULLIFY (qs_env%admm_env)
    NULLIFY (qs_env%efield)
    NULLIFY (qs_env%lri_env)
    NULLIFY (qs_env%lri_density)
    NULLIFY (qs_env%hfx_ri_env)
    NULLIFY (qs_env%dispersion_env)
    NULLIFY (qs_env%rtp)
    NULLIFY (qs_env%mp2_env)
    NULLIFY (qs_env%kg_env)
    NULLIFY (qs_env%WannierCentres)
    NULLIFY (qs_env%atprop)

    qs_env%outer_scf_ihistory=0
    qs_env%broyden_adaptive_sigma=-1.0_dp
    qs_env%neighbor_list_id=-1
    qs_env%nelectron_total=-1
    qs_env%nelectron_spin=-1

    CALL local_rho_set_create(qs_env%local_rho_set, error=error)
    CALL hartree_local_create(qs_env%hartree_local, error=error)
    CALL cp_subsys_create(qs_env%subsys, para_env=para_env, error=error)
    CALL cp_para_env_retain(para_env,error=error)
    qs_env%para_env => para_env
    CALL cp_blacs_env_create(qs_env%blacs_env,qs_env%para_env,&
                             globenv%blacs_grid_layout,&
                             globenv%blacs_repeatable,error=error)
    qs_env%ref_count=1
    last_qs_env_id_nr=last_qs_env_id_nr+1
    qs_env%id_nr=last_qs_env_id_nr
    qs_env%run_rtp=.FALSE.
    qs_env%linres_run=.FALSE.
    qs_env%use_ref_cell=.FALSE.
    qs_env%qmmm=.FALSE.
    qs_env%qmmm_periodic=.FALSE.
    qs_env%requires_mo_derivs=.FALSE.
    qs_env%requires_matrix_vxc=.FALSE.
    qs_env%has_unit_metric=.FALSE.
    qs_env%use_harris=.FALSE.
    qs_env%calc_image_preconditioner=.TRUE.
    qs_env%target_time = globenv%cp2k_target_time
    qs_env%start_time = globenv%cp2k_start_time

    qs_env%sim_time=0._dp
    qs_env%sim_step=0
    ! Zero all variables containing results
    CALL cp_result_create(qs_env%results, error)
    NULLIFY(qs_env%mo_derivs)
    NULLIFY(qs_env%mo_derivs_aux_fit)
    NULLIFY(qs_env%mo_loc_history)

  END SUBROUTINE init_qs_env

! *****************************************************************************
!> \brief   Set the QUICKSTEP environment.
!> \author  MK
!> \date    23.01.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_qs_env(qs_env,cell,super_cell,cell_ref,use_ref_cell,dft_control,&
       mos,mos_aux_fit,qmmm,qmmm_periodic,sab_orb,sab_aux_fit,sab_aux_fit_asymm, &
       sab_aux_fit_vs_orb,sab_all,sac_ae,sac_ppl,sap_ppnl,sap_oce,sab_lrc,&
       sab_se,sab_tbe,sab_core,sab_scp,sab_vdw,sab_almo,energy,force,harris_env,&
       matrix_h,matrix_ks,matrix_vxc,&
       matrix_ks_aux_fit, matrix_s,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb,matrix_w,matrix_w_mp2,matrix_p_mp2,&
       gamma_matrix,pw_env,ewald_env,ewald_pw,mpools,mpools_aux_fit,rho,vee,&
       rho_external,external_vxc,mask,&
       rho_aux_fit,rho_buffer,rho_xc,kinetic,scf_control,rel_control,qs_charges,vppl,rho_core,rho_nlcc, rho_nlcc_g, ks_env,&
       ks_qmmm_env,wf_history,scf_env,use_harris,task_list,task_list_aux_fit,task_list_soft,&
       task_list_zmp,&
       distribution_2d,dbcsr_dist,input, subsys, oce,rho_atom_set,rho0_atom_set,rho0_mpole,run_rtp,rtp,&
       rhoz_set,rhoz_tot,ecoul_1c,has_unit_metric,requires_mo_derivs,mo_derivs,&
       mo_derivs_aux_fit, mo_loc_history, nelectron_total, nelectron_spin, efield,&
       neighbor_list_id,linres_control,xas_env,virial,cp_ddapc_env,cp_ddapc_ewald,&
       outer_scf_history,outer_scf_ihistory,ep_qs_env,x_data,et_coupling,dftb_potential,&
       scp_env,se_taper,se_store_int_env,se_nddo_mpole,se_nonbond_env,admm_env,ls_scf_env,&
       lri_env,lri_density,hfx_ri_env,dispersion_env,mp2_env,kg_env,WannierCentres,ga_env,atprop,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, super_cell, cell_ref
    LOGICAL, OPTIONAL                        :: use_ref_cell
    TYPE(dft_control_type), OPTIONAL, &
      POINTER                                :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mos, mos_aux_fit
    LOGICAL, OPTIONAL                        :: qmmm, qmmm_periodic
    TYPE(neighbor_list_set_p_type), DIMENSION(:), OPTIONAL, POINTER :: &
      sab_orb, sab_aux_fit, sab_aux_fit_asymm, sab_aux_fit_vs_orb, sab_all, &
      sac_ae, sac_ppl, sap_ppnl, sap_oce, sab_lrc, sab_se, sab_tbe, sab_core, &
      sab_scp, sab_vdw, sab_almo
    TYPE(qs_energy_type), OPTIONAL, POINTER  :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: force
    TYPE(harris_env_type), OPTIONAL, POINTER :: harris_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), OPTIONAL, POINTER :: matrix_h, &
      matrix_ks, matrix_vxc, matrix_ks_aux_fit, matrix_s, matrix_s_aux_fit, &
      matrix_s_aux_fit_vs_orb, matrix_w, matrix_w_mp2, matrix_p_mp2, &
      gamma_matrix
    TYPE(pw_env_type), OPTIONAL, POINTER     :: pw_env
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(ewald_pw_type), OPTIONAL, POINTER   :: ewald_pw
    TYPE(qs_matrix_pools_type), OPTIONAL, &
      POINTER                                :: mpools, mpools_aux_fit
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho
    TYPE(pw_p_type), OPTIONAL, POINTER       :: vee
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho_external
    TYPE(pw_p_type), OPTIONAL, POINTER       :: external_vxc, mask
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho_aux_fit, rho_buffer, &
                                                rho_xc
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: kinetic
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                :: scf_control
    TYPE(rel_control_type), OPTIONAL, &
      POINTER                                :: rel_control
    TYPE(qs_charges_type), OPTIONAL, POINTER :: qs_charges
    TYPE(pw_p_type), OPTIONAL, POINTER       :: vppl, rho_core, rho_nlcc, &
                                                rho_nlcc_g
    TYPE(qs_ks_env_type), OPTIONAL, POINTER  :: ks_env
    TYPE(qs_ks_qmmm_env_type), OPTIONAL, &
      POINTER                                :: ks_qmmm_env
    TYPE(qs_wf_history_type), OPTIONAL, &
      POINTER                                :: wf_history
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    LOGICAL, OPTIONAL                        :: use_harris
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list, task_list_aux_fit, &
                                                task_list_soft, task_list_zmp
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    TYPE(dbcsr_distribution_obj), OPTIONAL, &
      POINTER                                :: dbcsr_dist
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(oce_matrix_type), OPTIONAL, POINTER :: oce
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_atom_set
    TYPE(rho0_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho0_atom_set
    TYPE(rho0_mpole_type), OPTIONAL, POINTER :: rho0_mpole
    LOGICAL, OPTIONAL                        :: run_rtp
    TYPE(rt_prop_type), OPTIONAL, POINTER    :: rtp
    TYPE(rhoz_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rhoz_set
    REAL(dp), OPTIONAL                       :: rhoz_tot
    TYPE(ecoul_1center_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: ecoul_1c
    LOGICAL, OPTIONAL                        :: has_unit_metric, &
                                                requires_mo_derivs
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mo_derivs
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mo_derivs_aux_fit, &
                                                mo_loc_history
    INTEGER, OPTIONAL                        :: nelectron_total
    INTEGER, DIMENSION(2), OPTIONAL          :: nelectron_spin
    TYPE(efield_berry_type), OPTIONAL, &
      POINTER                                :: efield
    INTEGER, OPTIONAL                        :: neighbor_list_id
    TYPE(linres_control_type), OPTIONAL, &
      POINTER                                :: linres_control
    TYPE(xas_environment_type), OPTIONAL, &
      POINTER                                :: xas_env
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_ddapc_type), OPTIONAL, POINTER   :: cp_ddapc_env
    TYPE(cp_ddapc_ewald_type), OPTIONAL, &
      POINTER                                :: cp_ddapc_ewald
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: outer_scf_history
    INTEGER, INTENT(IN), OPTIONAL            :: outer_scf_ihistory
    TYPE(ep_qs_type), OPTIONAL, POINTER      :: ep_qs_env
    TYPE(hfx_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: x_data
    TYPE(et_coupling_type), OPTIONAL, &
      POINTER                                :: et_coupling
    TYPE(qs_dftb_pairpot_type), &
      DIMENSION(:, :), OPTIONAL, POINTER     :: dftb_potential
    TYPE(scp_environment_type), OPTIONAL, &
      POINTER                                :: scp_env
    TYPE(se_taper_type), OPTIONAL, POINTER   :: se_taper
    TYPE(semi_empirical_si_type), OPTIONAL, &
      POINTER                                :: se_store_int_env
    TYPE(nddo_mpole_type), OPTIONAL, POINTER :: se_nddo_mpole
    TYPE(fist_nonbond_env_type), OPTIONAL, &
      POINTER                                :: se_nonbond_env
    TYPE(admm_type), OPTIONAL, POINTER       :: admm_env
    TYPE(ls_scf_env_type), OPTIONAL, POINTER :: ls_scf_env
    TYPE(lri_environment_type), OPTIONAL, &
      POINTER                                :: lri_env
    TYPE(lri_density_type), OPTIONAL, &
      POINTER                                :: lri_density
    TYPE(ri_environment_type), OPTIONAL, &
      POINTER                                :: hfx_ri_env
    TYPE(qs_dispersion_type), OPTIONAL, &
      POINTER                                :: dispersion_env
    TYPE(mp2_type), OPTIONAL, POINTER        :: mp2_env
    TYPE(kg_environment_type), OPTIONAL, &
      POINTER                                :: kg_env
    TYPE(wannier_centres_type), &
      DIMENSION(:), OPTIONAL, POINTER        :: WannierCentres
    TYPE(ga_environment_type), OPTIONAL, &
      POINTER                                :: ga_env
    TYPE(atprop_type), OPTIONAL, POINTER     :: atprop
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_qs_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
!    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (PRESENT(ga_env)) THEN
       CALL ga_env_retain(ga_env, error)
       CALL ga_env_release(qs_env%ga_env,error)
       qs_env%ga_env => ga_env
    END IF
    IF (PRESENT(cell)) THEN
       CALL cell_retain(cell, error=error)
       CALL cell_release(qs_env%cell,error=error)
       qs_env%cell => cell
    END IF
    IF (PRESENT(nelectron_total)) qs_env%nelectron_total=nelectron_total
    IF (PRESENT(mp2_env)) qs_env%mp2_env=>mp2_env
    IF (PRESENT(kg_env)) qs_env%kg_env => kg_env
    IF (PRESENT(nelectron_spin))  qs_env%nelectron_spin=nelectron_spin
    IF (PRESENT(super_cell)) THEN
       CALL cell_retain(super_cell, error=error)
       CALL cell_release(qs_env%super_cell,error=error)
       qs_env%super_cell => super_cell
    END IF
    IF (PRESENT(cell_ref)) THEN
       CALL cell_retain(cell_ref, error=error)
       CALL cell_release(qs_env%cell_ref,error=error)
       qs_env%cell_ref => cell_ref
    END IF
    IF (PRESENT(dft_control)) THEN
       CALL dft_control_retain(dft_control, error=error)
       CALL dft_control_release(qs_env%dft_control, error=error)
       qs_env%dft_control => dft_control
    END IF
    IF ( PRESENT ( use_ref_cell ) ) qs_env%use_ref_cell = use_ref_cell
    IF ( PRESENT ( qmmm ) ) qs_env%qmmm = qmmm
    IF ( PRESENT ( qmmm_periodic ) ) qs_env%qmmm_periodic = qmmm_periodic
    IF (PRESENT(mos)) qs_env%mos => mos
    IF (PRESENT(mos_aux_fit)) qs_env%mos_aux_fit => mos_aux_fit
    IF (PRESENT(neighbor_list_id)) qs_env%neighbor_list_id=neighbor_list_id
    IF (PRESENT(sab_orb)) qs_env%sab_orb => sab_orb
    IF (PRESENT(sab_aux_fit)) qs_env%sab_aux_fit => sab_aux_fit
    IF (PRESENT(sab_aux_fit_asymm)) qs_env%sab_aux_fit_asymm => sab_aux_fit_asymm
    IF (PRESENT(sab_aux_fit_vs_orb)) qs_env%sab_aux_fit_vs_orb => sab_aux_fit_vs_orb
    IF (PRESENT(sab_scp)) qs_env%sab_scp => sab_scp
    IF (PRESENT(sab_vdw)) qs_env%sab_vdw => sab_vdw
    IF (PRESENT(sab_all)) qs_env%sab_all => sab_all
    IF (PRESENT(sac_ae)) qs_env%sac_ae => sac_ae
    IF (PRESENT(sac_ppl)) qs_env%sac_ppl => sac_ppl
    IF (PRESENT(sap_ppnl)) qs_env%sap_ppnl => sap_ppnl
    IF (PRESENT(ls_scf_env)) qs_env%ls_scf_env => ls_scf_env
    IF (PRESENT(sap_oce)) qs_env%sap_oce => sap_oce
    IF (PRESENT(sab_se)) qs_env%sab_se => sab_se
    IF (PRESENT(sab_lrc)) qs_env%sab_lrc => sab_lrc
    IF (PRESENT(sab_tbe)) qs_env%sab_tbe => sab_tbe
    IF (PRESENT(sab_core)) qs_env%sab_core => sab_core
    IF (PRESENT(sab_almo)) qs_env%sab_almo => sab_almo
    IF (PRESENT(energy)) qs_env%energy => energy
    ! if intels checking (-C) complains here, you have rediscovered a bug in the intel
    ! compiler (present in at least 10.0.025). A testcase has been submitted to intel.
    IF (PRESENT(force)) qs_env%force => force
    IF (PRESENT(use_harris)) qs_env%use_harris = use_harris
    IF (PRESENT(harris_env)) qs_env%harris_env => harris_env
    IF (PRESENT(matrix_h)) qs_env%matrix_h => matrix_h
    IF (PRESENT(matrix_ks)) qs_env%matrix_ks => matrix_ks
    IF (PRESENT(matrix_vxc)) qs_env%matrix_vxc => matrix_vxc
    IF (PRESENT(matrix_ks_aux_fit)) qs_env%matrix_ks_aux_fit => matrix_ks_aux_fit
    IF (PRESENT(matrix_s)) qs_env%matrix_s => matrix_s
    IF (PRESENT(matrix_s_aux_fit)) qs_env%matrix_s_aux_fit => matrix_s_aux_fit
    IF (PRESENT(matrix_s_aux_fit_vs_orb)) qs_env%matrix_s_aux_fit_vs_orb => matrix_s_aux_fit_vs_orb
    IF (PRESENT(matrix_w)) qs_env%matrix_w => matrix_w
    IF (PRESENT(matrix_w_mp2)) qs_env%matrix_w_mp2 => matrix_w_mp2
    IF (PRESENT(matrix_p_mp2)) qs_env%matrix_p_mp2 => matrix_p_mp2
    IF (PRESENT(gamma_matrix)) qs_env%gamma_matrix => gamma_matrix
    IF (PRESENT(oce)) qs_env%oce => oce
    IF (PRESENT(outer_scf_history)) qs_env%outer_scf_history=>outer_scf_history
    IF (PRESENT(outer_scf_ihistory)) qs_env%outer_scf_ihistory=outer_scf_ihistory
    IF (PRESENT(requires_mo_derivs)) qs_env%requires_mo_derivs=requires_mo_derivs
    IF (PRESENT(has_unit_metric)) qs_env%has_unit_metric=has_unit_metric
    IF (PRESENT(mo_derivs)) qs_env%mo_derivs=>mo_derivs
    IF (PRESENT(mo_derivs_aux_fit)) qs_env%mo_derivs_aux_fit=>mo_derivs_aux_fit
    IF (PRESENT(mo_loc_history)) qs_env%mo_loc_history=>mo_loc_history
    IF (PRESENT(task_list)) qs_env%task_list=>task_list
    IF (PRESENT(task_list_aux_fit)) qs_env%task_list_aux_fit=>task_list_aux_fit
    IF (PRESENT(task_list_soft)) qs_env%task_list_soft=>task_list_soft
! ZMP pointer
    IF (PRESENT(task_list_zmp)) qs_env%task_list_zmp=>task_list_zmp
    IF (PRESENT(run_rtp)) qs_env%run_rtp=run_rtp
    IF (PRESENT(rtp)) qs_env%rtp=>rtp
    IF (PRESENT(efield)) qs_env%efield=>efield
    IF (PRESENT(atprop)) qs_env%atprop=>atprop
    IF (PRESENT(rho)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho)) CALL qs_rho_retain(rho,error=error)
      CALL qs_rho_release(qs_env%rho,error=error)
      qs_env%rho => rho
    END IF
    IF (PRESENT(rho_aux_fit)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho_aux_fit)) CALL qs_rho_retain(rho_aux_fit,error=error)
      CALL qs_rho_release(qs_env%rho_aux_fit,error=error)
      qs_env%rho_aux_fit => rho_aux_fit
    END IF
    IF (PRESENT(rho_buffer)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho_buffer)) CALL qs_rho_retain(rho_buffer,error=error)
      CALL qs_rho_release(qs_env%rho_buffer,error=error)
      qs_env%rho_buffer => rho_buffer
    END IF
    IF (PRESENT(rho_xc)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho_xc)) CALL qs_rho_retain(rho_xc,error=error)
      CALL qs_rho_release(qs_env%rho_xc,error=error)
      qs_env%rho_xc => rho_xc
    END IF

    IF (PRESENT(pw_env)) THEN ! accept also null pointers?
      CALL pw_env_retain(pw_env,error=error)
      CALL pw_env_release(qs_env%pw_env,error=error)
      qs_env%pw_env => pw_env
    END IF
    IF (PRESENT(ewald_env)) THEN ! accept also null pointers?
      CALL ewald_env_retain(ewald_env,error=error)
      CALL ewald_env_release(qs_env%ewald_env,error=error)
      qs_env%ewald_env => ewald_env
    END IF
    IF (PRESENT(ewald_pw)) THEN ! accept also null pointers?
      CALL ewald_pw_retain(ewald_pw,error=error)
      CALL ewald_pw_release(qs_env%ewald_pw,error=error)
      qs_env%ewald_pw => ewald_pw
    END IF
    IF (PRESENT(scf_control)) THEN ! accept also null pointers?
       CALL scf_c_retain(scf_control,error=error)
       CALL scf_c_release(qs_env%scf_control,error=error)
       qs_env%scf_control => scf_control
    END IF
    IF (PRESENT(rel_control)) THEN ! accept also null pointers?
       CALL rel_c_retain(rel_control,error=error)
       CALL rel_c_release(qs_env%rel_control,error=error)
       qs_env%rel_control => rel_control
    END IF
    IF (PRESENT(linres_control)) THEN ! accept also null pointers?
       CALL linres_control_retain(linres_control,error=error)
       CALL linres_control_release(qs_env%linres_control,error=error)
       qs_env%linres_control => linres_control
    END IF
    IF (PRESENT(kinetic)) qs_env%kinetic => kinetic
    IF (PRESENT(rho_core)) qs_env%rho_core => rho_core
    IF (PRESENT(rho_nlcc)) qs_env%rho_nlcc => rho_nlcc
    IF (PRESENT(rho_nlcc_g)) qs_env%rho_nlcc_g => rho_nlcc_g
    IF (PRESENT(vppl)) qs_env%vppl => vppl
    IF (PRESENT(vee)) qs_env%vee => vee
! ZMP associating variables
    IF (PRESENT(rho_external)) THEN ! accepts also null pointers !
       IF (ASSOCIATED(rho_external)) CALL qs_rho_retain(rho_external,error=error)
       CALL qs_rho_release(qs_env%rho_external,error=error)
       qs_env%rho_external => rho_external
    END IF
    IF (PRESENT(external_vxc)) qs_env%external_vxc => external_vxc
    IF (PRESENT(mask)) qs_env%mask => mask

    IF (PRESENT(qs_charges)) THEN
       CALL qs_charges_retain(qs_charges,error=error)
       CALL qs_charges_release(qs_env%qs_charges,error=error)
       qs_env%qs_charges => qs_charges
    END IF
    IF (PRESENT(ks_qmmm_env)) THEN
       CALL qs_ks_qmmm_retain(ks_qmmm_env, error=error)
       CALL qs_ks_qmmm_release(qs_env%ks_qmmm_env, error=error)
       qs_env%ks_qmmm_env => ks_qmmm_env
    END IF
    IF (PRESENT(ks_env)) THEN ! accept also null pointers?
       CALL qs_ks_retain(ks_env, error=error)
       CALL qs_ks_release(qs_env%ks_env, error=error)
       qs_env%ks_env => ks_env
    END IF
    IF (PRESENT(wf_history)) THEN ! accept also null pointers ?
       CALL wfi_retain(wf_history,error=error)
       CALL wfi_release(qs_env%wf_history,error=error)
       qs_env%wf_history => wf_history
    END IF
    IF (PRESENT(scf_env)) THEN ! accept also null pointers ?
       CALL scf_env_retain(scf_env,error=error)
       CALL scf_env_release(qs_env%scf_env, error=error)
       qs_env%scf_env => scf_env
    END IF
    IF (PRESENT(xas_env)) THEN ! accept also null pointers?
       CALL xas_env_retain(xas_env, error=error)
       CALL xas_env_release(qs_env%xas_env, error=error)
       qs_env%xas_env => xas_env
    END IF
    IF (PRESENT(virial)) THEN
       IF (ASSOCIATED(virial)) THEN
          CALL virial_retain(virial,error=error)
       END IF
       CALL virial_release(qs_env%virial,error=error)
       qs_env%virial => virial
     END IF
    IF (PRESENT(mpools)) THEN
       CALL mpools_retain(mpools,error=error)
       CALL mpools_release(qs_env%mpools, error=error)
       qs_env%mpools => mpools
    END IF
    IF (PRESENT(mpools_aux_fit)) THEN
       CALL mpools_retain(mpools_aux_fit,error=error)
       CALL mpools_release(qs_env%mpools_aux_fit, error=error)
       qs_env%mpools_aux_fit => mpools_aux_fit
    END IF
    IF (PRESENT(rho_atom_set)) THEN
        CALL set_local_rho(qs_env%local_rho_set,rho_atom_set=rho_atom_set,error=error)
     END IF
    IF (PRESENT(rho0_atom_set)) THEN
        CALL set_local_rho(qs_env%local_rho_set,rho0_atom_set=rho0_atom_set,error=error)
    END IF
    IF (PRESENT(rho0_mpole)) THEN
        CALL set_local_rho(qs_env%local_rho_set,rho0_mpole=rho0_mpole,error=error)
    END IF
    IF (PRESENT(rhoz_set)) THEN
        CALL set_local_rho(qs_env%local_rho_set,rhoz_set=rhoz_set,error=error)
    END IF
    IF (PRESENT(rhoz_tot)) qs_env%local_rho_set%rhoz_tot = rhoz_tot
    IF (PRESENT(ecoul_1c)) THEN
        CALL set_hartree_local(qs_env%hartree_local,ecoul_1c=ecoul_1c)
    END IF

    IF (PRESENT(subsys)) THEN
       CALL cp_subsys_retain(subsys,error=error)
       CALL cp_subsys_release(qs_env%subsys,error=error)
       qs_env%subsys => subsys
    END IF
    IF (PRESENT(distribution_2d)) THEN
       CALL distribution_2d_retain(distribution_2d,error=error)
       CALL distribution_2d_release(qs_env%distribution_2d,error=error)
       qs_env%distribution_2d => distribution_2d
    ENDIF
    IF (PRESENT(dbcsr_dist)) THEN
       IF(ASSOCIATED(qs_env%dbcsr_dist)) THEN
          CALL cp_dbcsr_distribution_release(qs_env%dbcsr_dist)
          DEALLOCATE(qs_env%dbcsr_dist)
       ENDIF
       qs_env%dbcsr_dist => dbcsr_dist
    ENDIF
    IF (PRESENT(input)) THEN
       CALL section_vals_retain(input,error=error)
       CALL section_vals_release(qs_env%input,error=error)
       qs_env%input => input
    END IF
    IF (PRESENT(cp_ddapc_env)) THEN
       CALL cp_ddapc_retain(cp_ddapc_env, error=error)
       CALL cp_ddapc_release(qs_env%cp_ddapc_env, error=error)
       qs_env%cp_ddapc_env => cp_ddapc_env
    END IF
    IF (PRESENT(cp_ddapc_ewald)) THEN
       qs_env%cp_ddapc_ewald => cp_ddapc_ewald
    END IF
    IF (PRESENT(ep_qs_env)) THEN
       IF (ASSOCIATED(ep_qs_env)) CALL ep_qs_retain(ep_qs_env,error=error)
       CALL ep_qs_release(qs_env%ep_qs_env,error=error)
       qs_env%ep_qs_env => ep_qs_env
    END IF
    IF (PRESENT(x_data)) qs_env%x_data => x_data
    IF (PRESENT(et_coupling))qs_env%et_coupling => et_coupling
    IF (PRESENT(dftb_potential))qs_env%dftb_potential => dftb_potential
    IF (PRESENT(scp_env)) THEN
      CALL scp_env_retain(scp_env,error=error)
      CALL scp_env_release(qs_env%scp_env,error=error)
      qs_env%scp_env => scp_env
    END IF
    IF (PRESENT(se_taper)) THEN
       CALL se_taper_release(qs_env%se_taper,error=error)
       qs_env%se_taper => se_taper
    END IF
    IF (PRESENT(se_store_int_env)) THEN
       CALL semi_empirical_si_release(qs_env%se_store_int_env,error=error)
       qs_env%se_store_int_env => se_store_int_env
    END IF
    IF (PRESENT(se_nddo_mpole)) THEN
       CALL nddo_mpole_release(qs_env%se_nddo_mpole,error=error)
       qs_env%se_nddo_mpole => se_nddo_mpole
    END IF
    IF (PRESENT(se_nonbond_env)) THEN
       CALL fist_nonbond_env_release(qs_env%se_nonbond_env,error)
       qs_env%se_nonbond_env => se_nonbond_env
    END IF
    IF( PRESENT(admm_env) ) qs_env%admm_env => admm_env
    IF( PRESENT(lri_env) ) qs_env%lri_env => lri_env
    IF( PRESENT(lri_density) ) qs_env%lri_density => lri_density
    IF( PRESENT(hfx_ri_env) ) qs_env%hfx_ri_env => hfx_ri_env
    IF( PRESENT(dispersion_env) ) qs_env%dispersion_env => dispersion_env
    IF( PRESENT(WannierCentres) ) qs_env%WannierCentres => WannierCentres
  END SUBROUTINE set_qs_env

! *****************************************************************************
!> \brief allocates and intitializes a qs_env
!> \param qs_env the object to create
!> \param para_env the parallel environement for the qs_env
!> \param blacs_env the blacs env for the qs_env (for full matrixes,
!>        it must be coerent with the para_env)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_env_create(qs_env,para_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(qs_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL init_qs_env(qs_env,para_env=para_env,&
            globenv=globenv, error=error)
    END IF
  END SUBROUTINE qs_env_create

! *****************************************************************************
!> \brief retains the given qs_env (see doc/ReferenceCounting.html)
!> \param qs_env the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_env_retain(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
       qs_env%ref_count=qs_env%ref_count+1
    END IF
  END SUBROUTINE qs_env_retain

! *****************************************************************************
!> \brief releases the given qs_env (see doc/ReferenceCounting.html)
!> \param qs_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      12.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_env_release(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(qs_env)) THEN
       CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
       qs_env%ref_count=qs_env%ref_count-1
       IF (qs_env%ref_count<1) THEN
          CALL ga_env_release(qs_env%ga_env,error=error)
          CALL cell_release(qs_env%cell,error=error)
          CALL cell_release(qs_env%super_cell,error=error)
          CALL cell_release(qs_env%cell_ref,error=error)
          CALL dft_control_release(qs_env%dft_control, error=error)
          IF (ASSOCIATED(qs_env%mos)) THEN
             DO i=1,SIZE(qs_env%mos)
                CALL deallocate_mo_set(qs_env%mos(i)%mo_set,error=error)
             END DO
             DEALLOCATE(qs_env%mos, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%mos_aux_fit)) THEN
             DO i=1,SIZE(qs_env%mos_aux_fit)
                CALL deallocate_mo_set(qs_env%mos_aux_fit(i)%mo_set,error=error)
             END DO
             DEALLOCATE(qs_env%mos_aux_fit, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF

          IF (ASSOCIATED(qs_env%mo_derivs)) THEN
             DO I=1,SIZE(qs_env%mo_derivs)
                CALL cp_dbcsr_release_p(qs_env%mo_derivs(I)%matrix, error=error)
             ENDDO
             DEALLOCATE(qs_env%mo_derivs, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          ENDIF

          IF (ASSOCIATED(qs_env%mo_derivs_aux_fit)) THEN
             DO I=1,SIZE(qs_env%mo_derivs_aux_fit)
                CALL cp_fm_release(qs_env%mo_derivs_aux_fit(I)%matrix,error=error)
             ENDDO
             DEALLOCATE(qs_env%mo_derivs_aux_fit, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          ENDIF

          IF (ASSOCIATED(qs_env%mo_loc_history)) THEN
             DO I=1,SIZE(qs_env%mo_loc_history)
                CALL cp_fm_release(qs_env%mo_loc_history(I)%matrix,error=error)
             ENDDO
             DEALLOCATE(qs_env%mo_loc_history, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          ENDIF

          IF (ASSOCIATED(qs_env%sab_scp)) THEN
             DO i=1,SIZE(qs_env%sab_scp)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_scp(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_scp,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_vdw)) THEN
             DO i=1,SIZE(qs_env%sab_vdw)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_vdw(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_vdw,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_orb)) THEN
             DO i=1,SIZE(qs_env%sab_orb)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_orb(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_orb,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_aux_fit)) THEN
             DO i=1,SIZE(qs_env%sab_aux_fit)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_aux_fit(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_aux_fit,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_aux_fit_asymm)) THEN
             DO i=1,SIZE(qs_env%sab_aux_fit_asymm)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_aux_fit_asymm(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_aux_fit_asymm,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_aux_fit_vs_orb)) THEN
             DO i=1,SIZE(qs_env%sab_aux_fit_vs_orb)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_aux_fit_vs_orb(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_aux_fit_vs_orb,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF

          IF (ASSOCIATED(qs_env%sab_all)) THEN
             DO i=1,SIZE(qs_env%sab_all)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_all(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_all,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sac_ae)) THEN
             DO i=1,SIZE(qs_env%sac_ae)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sac_ae(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sac_ae,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sac_ppl)) THEN
             DO i=1,SIZE(qs_env%sac_ppl)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sac_ppl(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sac_ppl,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sap_ppnl)) THEN
             DO i=1,SIZE(qs_env%sap_ppnl)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sap_ppnl(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sap_ppnl,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sap_oce)) THEN
             DO i=1,SIZE(qs_env%sap_oce)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sap_oce(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sap_oce,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_se)) THEN
             DO i=1,SIZE(qs_env%sab_se)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_se(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_se,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_lrc)) THEN
             DO i=1,SIZE(qs_env%sab_lrc)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_lrc(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_lrc,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_tbe)) THEN
             DO i=1,SIZE(qs_env%sab_tbe)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_tbe(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_tbe,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_core)) THEN
             DO i=1,SIZE(qs_env%sab_core)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_core(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_core,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%sab_almo)) THEN
             DO i=1,SIZE(qs_env%sab_almo)
                CALL deallocate_neighbor_list_set(&
                     qs_env%sab_almo(i)%neighbor_list_set)
             END DO
             DEALLOCATE(qs_env%sab_almo,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(qs_env%rtp))THEN
             CALL rt_prop_release(qs_env%rtp,error)
             DEALLOCATE (qs_env%rtp)
          END IF
          IF (ASSOCIATED(qs_env%outer_scf_history)) THEN
             DEALLOCATE(qs_env%outer_scf_history,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             qs_env%outer_scf_ihistory=0
          ENDIF
          IF (ASSOCIATED(qs_env%energy)) &
               CALL deallocate_qs_energy(qs_env%energy)
          IF (ASSOCIATED(qs_env%force))&
               CALL deallocate_qs_force(qs_env%force)
          IF (ASSOCIATED(qs_env%harris_env)) THEN
             CALL harris_env_release(qs_env%harris_env, error=error)
          END IF
          IF (ASSOCIATED(qs_env%matrix_h)) THEN
             CALL cp_dbcsr_deallocate_matrix_set(qs_env%matrix_h,error)
          ENDIF
          IF (ASSOCIATED(qs_env%matrix_ks)) THEN
             CALL cp_dbcsr_deallocate_matrix_set(qs_env % matrix_ks,error)
          ENDIF
          IF (ASSOCIATED(qs_env%matrix_vxc)) THEN
             CALL cp_dbcsr_deallocate_matrix_set(qs_env % matrix_vxc,error)
          ENDIF
          IF (ASSOCIATED(qs_env%matrix_ks_aux_fit)) THEN
             CALL cp_dbcsr_deallocate_matrix_set(qs_env % matrix_ks_aux_fit,error)
          ENDIF
          IF (ASSOCIATED(qs_env%matrix_s)) CALL cp_dbcsr_deallocate_matrix_set(qs_env%matrix_s,error)
          IF (ASSOCIATED(qs_env%matrix_s_aux_fit)) CALL cp_dbcsr_deallocate_matrix_set(qs_env%matrix_s_aux_fit,error)
          IF (ASSOCIATED(qs_env%matrix_s_aux_fit_vs_orb)) CALL cp_dbcsr_deallocate_matrix_set(qs_env%matrix_s_aux_fit_vs_orb,error)
          IF (ASSOCIATED(qs_env%matrix_w)) CALL cp_dbcsr_deallocate_matrix_set(qs_env%matrix_w,error)
          IF (ASSOCIATED(qs_env%matrix_w_mp2)) CALL cp_dbcsr_deallocate_matrix_set(qs_env%matrix_w_mp2,error)
          IF (ASSOCIATED(qs_env%matrix_p_mp2)) CALL cp_dbcsr_deallocate_matrix_set(qs_env%matrix_p_mp2,error)
          IF (ASSOCIATED(qs_env%gamma_matrix)) CALL cp_dbcsr_deallocate_matrix_set(qs_env%gamma_matrix,error)
          IF (ASSOCIATED(qs_env%oce)) CALL deallocate_oce_set(qs_env%oce,error)

          IF (ASSOCIATED(qs_env%local_rho_set)) THEN
             CALL local_rho_set_release(qs_env%local_rho_set,error=error)
          END IF
          IF (ASSOCIATED(qs_env%hartree_local)) THEN
             CALL hartree_local_release(qs_env%hartree_local,error=error)
          END IF
          IF (ASSOCIATED(qs_env%kinetic))&
               CALL cp_dbcsr_deallocate_matrix_set(qs_env%kinetic,error)

          IF(ASSOCIATED(qs_env%rho_xc)) THEN
             NULLIFY(qs_env%rho_xc%rho_ao)
             CALL qs_rho_release(qs_env%rho_xc, error=error)
          END IF
          CALL qs_rho_release(qs_env%rho, error=error)
          CALL qs_rho_release(qs_env%rho_aux_fit, error=error)
          CALL qs_rho_release(qs_env%rho_buffer,error=error)
          CALL scf_c_release(qs_env%scf_control, error=error)
          CALL rel_c_release(qs_env%rel_control, error=error)

          IF(ASSOCIATED(qs_env%linres_control)) THEN
             CALL linres_control_release(qs_env%linres_control, error=error)
          END IF

          IF(ASSOCIATED(qs_env%task_list)) THEN
             CALL deallocate_task_list(qs_env%task_list,error)
          ENDIF
          IF(ASSOCIATED(qs_env%task_list_aux_fit)) THEN
             CALL deallocate_task_list(qs_env%task_list_aux_fit,error)
          ENDIF
          IF(ASSOCIATED(qs_env%task_list_soft)) THEN
             CALL deallocate_task_list(qs_env%task_list_soft,error)
          ENDIF
          IF (ASSOCIATED(qs_env%ls_scf_env)) THEN
             CALL ls_scf_release(qs_env%ls_scf_env,error)
          ENDIF
! ZMP 
          IF(ASSOCIATED(qs_env%task_list_zmp)) THEN
             CALL deallocate_task_list(qs_env%task_list_zmp,error)
          ENDIF

          !Only if do_xas_calculation
          IF(ASSOCIATED(qs_env%xas_env)) THEN
             CALL xas_env_release(qs_env%xas_env,error=error)
          END IF
          IF(ASSOCIATED(qs_env%virial)) THEN
             CALL virial_release(qs_env%virial,error=error)
          END IF
          CALL pw_env_release(qs_env%pw_env, error=error)
          CALL ewald_env_release(qs_env%ewald_env, error=error)
          CALL ewald_pw_release(qs_env%ewald_pw, error=error)
          CALL cp_para_env_release(qs_env%para_env,error=error)
          CALL cp_blacs_env_release(qs_env%blacs_env, error=error)
          IF (ASSOCIATED(qs_env%rho_nlcc_g)) THEN
             CALL pw_release(qs_env%rho_nlcc_g%pw,error=error)
             DEALLOCATE(qs_env%rho_nlcc_g)
          ENDIF
          IF (ASSOCIATED(qs_env%rho_nlcc)) THEN
             CALL pw_release(qs_env%rho_nlcc%pw,error=error)
             DEALLOCATE(qs_env%rho_nlcc)
          ENDIF
          IF (ASSOCIATED(qs_env%rho_core)) THEN
             CALL pw_release(qs_env%rho_core%pw,error=error)
             DEALLOCATE(qs_env%rho_core)
          ENDIF
          IF (ASSOCIATED(qs_env%image_matrix)) THEN
             DEALLOCATE(qs_env%image_matrix)
          ENDIF
          IF (ASSOCIATED(qs_env%ipiv)) THEN
             DEALLOCATE(qs_env%ipiv)
          ENDIF
          IF (ASSOCIATED(qs_env%image_coeff)) THEN
             DEALLOCATE(qs_env%image_coeff)
          ENDIF
          IF (ASSOCIATED(qs_env%vppl)) THEN
             CALL pw_release(qs_env%vppl%pw,error=error)
             DEALLOCATE(qs_env%vppl)
          ENDIF
          IF (ASSOCIATED(qs_env%vee)) THEN
             CALL pw_release(qs_env%vee%pw,error=error)
             DEALLOCATE(qs_env%vee)
          ENDIF
! ZMP
          IF(ASSOCIATED(qs_env%rho_external)) THEN
             NULLIFY(qs_env%rho_external%rho_ao)
             CALL qs_rho_release(qs_env%rho_external, error=error)
          END IF
          IF (ASSOCIATED(qs_env%external_vxc)) THEN
             CALL pw_release(qs_env%external_vxc%pw,error=error)
             DEALLOCATE(qs_env%external_vxc)
          ENDIF
          IF (ASSOCIATED(qs_env%mask)) THEN
             CALL pw_release(qs_env%mask%pw,error=error)
             DEALLOCATE(qs_env%mask)
          ENDIF

          CALL qs_charges_release(qs_env%qs_charges, error=error)
          CALL qs_ks_release(qs_env%ks_env, error=error)
          CALL qs_ks_qmmm_release(qs_env%ks_qmmm_env, error=error)
          CALL wfi_release(qs_env%wf_history,error=error)
          CALL scf_env_release(qs_env%scf_env, error=error)
          CALL mpools_release(qs_env%mpools,error=error)
          CALL mpools_release(qs_env%mpools_aux_fit,error=error)
          IF (ASSOCIATED(qs_env%distribution_2d)) THEN
             CALL distribution_2d_release(qs_env%distribution_2d,error=error)
          ENDIF
          IF(ASSOCIATED(qs_env%dbcsr_dist)) THEN
             CALL cp_dbcsr_distribution_release(qs_env%dbcsr_dist)
             DEALLOCATE(qs_env%dbcsr_dist)
          ENDIF
          CALL cp_subsys_release(qs_env%subsys,error=error)
          CALL section_vals_release(qs_env%input,error=error)
          CALL cp_ddapc_release(qs_env%cp_ddapc_env, error=error)
          CALL cp_ddapc_ewald_release(qs_env%cp_ddapc_ewald, error=error)
          CALL ep_qs_release(qs_env%ep_qs_env,error=error)
          CALL efield_berry_release(qs_env%efield,error=error)
          CALL cp_result_release(qs_env%results,error=error)
          IF(ASSOCIATED(qs_env%x_data)) THEN
             CALL hfx_release(qs_env%x_data, error=error)
          END IF
          IF(ASSOCIATED(qs_env%et_coupling)) THEN
             CALL et_coupling_release(qs_env%et_coupling,error)
          END IF
          IF (ASSOCIATED(qs_env%dftb_potential)) THEN
             CALL qs_dftb_pairpot_release(qs_env%dftb_potential,error)
          END IF
          IF (ASSOCIATED(qs_env%scp_env)) THEN
             CALL scp_env_release ( qs_env % scp_env, error )
          END IF
          IF (ASSOCIATED(qs_env%se_taper)) THEN
             CALL se_taper_release(qs_env%se_taper, error)
          END IF
          IF (ASSOCIATED(qs_env%se_store_int_env)) THEN
             CALL semi_empirical_si_release(qs_env%se_store_int_env, error)
          END IF
          IF (ASSOCIATED(qs_env%se_nddo_mpole)) THEN
             CALL nddo_mpole_release(qs_env%se_nddo_mpole, error)
          END IF
          IF (ASSOCIATED(qs_env%se_nonbond_env)) THEN
             CALL fist_nonbond_env_release(qs_env%se_nonbond_env,error)
          END IF
          IF (ASSOCIATED(qs_env%admm_env)) THEN
            CALL admm_env_release(qs_env%admm_env, error)
          END IF
          IF (ASSOCIATED(qs_env%lri_env)) THEN
            CALL lri_env_release(qs_env%lri_env, error)
          END IF
          IF (ASSOCIATED(qs_env%lri_density)) THEN
             CALL lri_density_release(qs_env%lri_density, error)
          END IF
          IF (ASSOCIATED(qs_env%hfx_ri_env)) THEN
            CALL ri_env_release(qs_env%hfx_ri_env, error)
          END IF
          IF (ASSOCIATED(qs_env%mp2_env)) THEN
            CALL mp2_env_release(qs_env%mp2_env, error)
          END IF
          IF (ASSOCIATED(qs_env%kg_env)) THEN
            CALL kg_env_release(qs_env%kg_env, error)
          END IF

          ! dispersion
          CALL qs_dispersion_release(qs_env%dispersion_env,error)

          IF( ASSOCIATED(qs_env%WannierCentres)) THEN
            DO i=1,SIZE(qs_env%WannierCentres)
              DEALLOCATE(qs_env%WannierCentres(i)%WannierHamDiag, stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
              DEALLOCATE(qs_env%WannierCentres(i)%centres, stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
            ENDDO
            DEALLOCATE(qs_env%WannierCentres, stat=stat)
            CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          ENDIF
          ! dereference atprop (we assume it will be deallocated in force_env!
          NULLIFY(qs_env%atprop)
          ! now we are ready to deallocate the full structure
          DEALLOCATE(qs_env, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(qs_env)
  END SUBROUTINE qs_env_release

! *****************************************************************************
!> \brief reorthogonalizes the mos
!> \param qs_env the qs_env in which to orthogonalize
!> \param v the vectors to orthogonalize
!> \param n_col number of column of v to orthogonalize
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      04.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qs_env_reorthogonalize_vectors(qs_env, v_matrix, n_col,&
       error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: v_matrix
    INTEGER, INTENT(in), OPTIONAL            :: n_col
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'qs_env_reorthogonalize_vectors', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, my_n_col
    LOGICAL                                  :: failure, has_unit_metric, &
                                                ortho_contains_cholesky, &
                                                smearing_is_used
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_pool_type), POINTER           :: maxao_maxmo_fm_pool
    TYPE(qs_matrix_pools_type), POINTER      :: mpools
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control

    failure=.FALSE.
    NULLIFY(scf_env, scf_control, maxao_maxmo_fm_pool, matrix_s, mpools)
    CALL timeset(routineN,handle)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(v_matrix),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN

       CALL cp_fm_get_info(v_matrix,ncol_global=my_n_col,error=error)
       IF (PRESENT(n_col)) my_n_col=n_col
       CALL get_qs_env(qs_env,mpools=mpools,&
            scf_env=scf_env,&
            scf_control=scf_control,&
            matrix_s=matrix_s,&
            error=error)
       CALL mpools_get(mpools,maxao_maxmo_fm_pool=maxao_maxmo_fm_pool,&
            error=error)
       IF (ASSOCIATED(scf_env)) THEN
          ortho_contains_cholesky=(scf_env%method /= ot_method_nr).AND.&
               (scf_env%cholesky_method>0 ).AND.&
               ASSOCIATED(scf_env%ortho)
       ELSE
          ortho_contains_cholesky=.FALSE.
       END IF

       CALL get_qs_env(qs_env,has_unit_metric=has_unit_metric,error=error)
       smearing_is_used = .FALSE.
       IF (qs_env%dft_control%smear) THEN
           smearing_is_used = .TRUE.
       END IF


       IF (has_unit_metric) THEN
          CALL make_basis_simple(v_matrix,my_n_col,error=error)
       ELSE IF (smearing_is_used) THEN
          CALL make_basis_lowdin(vmatrix=v_matrix,ncol=my_n_col,&
                 matrix_s=matrix_s(1)%matrix,error=error)
       ELSE IF (ortho_contains_cholesky) THEN
          CALL make_basis_cholesky(vmatrix=v_matrix,ncol=my_n_col,&
               ortho=scf_env%ortho,error=error)
       ELSE
          CALL make_basis_sm(v_matrix,my_n_col, matrix_s(1)%matrix ,error=error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE qs_env_reorthogonalize_vectors

END MODULE qs_environment_types
