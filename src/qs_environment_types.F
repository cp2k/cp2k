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
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atprop_types,                    ONLY: atprop_type
  USE cell_methods,                    ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_release_p,&
                                             dbcsr_distribution_obj
  USE cp_ddapc_types,                  ONLY: cp_ddapc_ewald_release,&
                                             cp_ddapc_ewald_type,&
                                             cp_ddapc_release,&
                                             cp_ddapc_retain,&
                                             cp_ddapc_type
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE distribution_2d_types,           ONLY: distribution_2d_type
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
  USE kpoint_types,                    ONLY: kpoint_type
  USE lri_environment_types,           ONLY: lri_density_release,&
                                             lri_density_type,&
                                             lri_env_release,&
                                             lri_environment_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE mp2_types,                       ONLY: mp2_env_release,&
                                             mp2_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_type
  USE pw_types,                        ONLY: pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE qs_charges_types,                ONLY: qs_charges_release,&
                                             qs_charges_retain,&
                                             qs_charges_type
  USE qs_dftb_types,                   ONLY: qs_dftb_pairpot_release,&
                                             qs_dftb_pairpot_type
  USE qs_dispersion_types,             ONLY: qs_dispersion_release,&
                                             qs_dispersion_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type,&
                                             qs_ks_qmmm_release,&
                                             qs_ks_qmmm_retain
  USE qs_ks_types,                     ONLY: get_ks_env,&
                                             qs_ks_env_type,&
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
  USE qs_matrix_pools,                 ONLY: mpools_release,&
                                             mpools_retain,&
                                             qs_matrix_pools_type
  USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                             mo_set_p_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
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
  USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                             scf_env_release,&
                                             scf_env_retain
  USE qs_subsys_types,                 ONLY: qs_subsys_type
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
  USE virial_types,                    ONLY: virial_type
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

  PUBLIC :: qs_environment_type

! *** Public subroutines ***

  PUBLIC :: get_qs_env,&
            qs_env_create,&
            qs_env_release,&
            qs_env_retain,&
            set_qs_env

! *****************************************************************************
!> \param local_rho_set contains the atomic, compensations and core densities
!>                       and the local parts of the xc terms
!> \param hartree_local contains the 1, 2 and 3 centers coulomb terms
!> \param requires_mo_derivs logical, true if dE/dC is required (e.g. OT)
!> \param has_unit_metric logical, true if the S matrix is considered unity for the SCF
!> \param mo_derivs the actual derivatives of the total energy wrt to MO coeffs (divided by 2*f_i)
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
    TYPE(cell_type), POINTER                              :: super_cell
    TYPE(mo_set_p_type), DIMENSION(:), POINTER            :: mos, mos_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER             :: mo_derivs_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER             :: mo_loc_history
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER          :: mo_derivs
    TYPE(harris_env_type), POINTER                        :: harris_env
    TYPE(scf_control_type), POINTER                       :: scf_control
    TYPE(rel_control_type), POINTER                       :: rel_control
    TYPE(scp_environment_type), POINTER                   :: scp_env
    TYPE(cp_para_env_type), POINTER                       :: para_env
    TYPE(cp_blacs_env_type), POINTER                      :: blacs_env

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
    TYPE(oce_matrix_type), POINTER                        :: oce
    TYPE(local_rho_type), POINTER                         :: local_rho_set
    TYPE(hartree_local_type),  POINTER                    :: hartree_local
    TYPE(section_vals_type), POINTER                      :: input
    TYPE(linres_control_type), POINTER                    :: linres_control
    TYPE(xas_environment_type), POINTER                   :: xas_env
    TYPE(cp_ddapc_type), POINTER                          :: cp_ddapc_env
    TYPE(cp_ddapc_ewald_type), POINTER                    :: cp_ddapc_ewald
    REAL(KIND = dp), DIMENSION(:,:), POINTER              :: outer_scf_history
    INTEGER                                               :: outer_scf_ihistory
    TYPE(ep_qs_type), POINTER                             :: ep_qs_env
    TYPE(hfx_type),DIMENSION(:,:), POINTER                :: x_data
    TYPE(et_coupling_type),POINTER                        :: et_coupling
    TYPE(qs_dftb_pairpot_type), DIMENSION(:,:), POINTER   :: dftb_potential
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
!> \param qs_env ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param super_cell ...
!> \param cell_ref ...
!> \param use_ref_cell ...
!> \param kpoints ...
!> \param dft_control ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param sab_orb ...
!> \param sab_aux_fit ...
!> \param sab_aux_fit_asymm ...
!> \param sab_aux_fit_vs_orb ...
!> \param sab_all ...
!> \param qmmm ...
!> \param qmmm_periodic ...
!> \param sac_ae ...
!> \param sac_ppl ...
!> \param sap_ppnl ...
!> \param sab_vdw ...
!> \param sab_scp ...
!> \param sap_oce ...
!> \param sab_lrc ...
!> \param sab_se ...
!> \param sab_tbe ...
!> \param sab_core ...
!> \param sab_almo ...
!> \param particle_set ...
!> \param energy ...
!> \param force ...
!> \param harris_env ...
!> \param matrix_h ...
!> \param matrix_ks ...
!> \param matrix_vxc ...
!> \param run_rtp ...
!> \param rtp ...
!> \param matrix_ks_aux_fit ...
!> \param matrix_s ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param matrix_w ...
!> \param matrix_w_mp2 ...
!> \param matrix_p_mp2 ...
!> \param gamma_matrix ...
!> \param rho ...
!> \param rho_aux_fit ...
!> \param rho_buffer ...
!> \param rho_xc ...
!> \param pw_env ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param mpools ...
!> \param mpools_aux_fit ...
!> \param input ...
!> \param para_env ...
!> \param blacs_env ...
!> \param scf_control ...
!> \param rel_control ...
!> \param kinetic ...
!> \param qs_charges ...
!> \param vppl ...
!> \param rho_core ...
!> \param rho_nlcc ...
!> \param rho_nlcc_g ...
!> \param ks_env ...
!> \param ks_qmmm_env ...
!> \param wf_history ...
!> \param scf_env ...
!> \param use_harris ...
!> \param id_nr ...
!> \param local_particles ...
!> \param local_molecules ...
!> \param distribution_2d ...
!> \param dbcsr_dist ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param subsys ...
!> \param cp_subsys ...
!> \param oce ...
!> \param rho_atom_set ...
!> \param task_list ...
!> \param task_list_zmp ...
!> \param task_list_aux_fit ...
!> \param task_list_soft ...
!> \param rho0_atom_set ...
!> \param rho0_mpole ...
!> \param rhoz_set ...
!> \param ecoul_1c ...
!> \param rho0_s_rs ...
!> \param rho0_s_gs ...
!> \param do_kpoints ...
!> \param has_unit_metric ...
!> \param requires_mo_derivs ...
!> \param mo_derivs ...
!> \param mo_derivs_aux_fit ...
!> \param mo_loc_history ...
!> \param nkind ...
!> \param natom ...
!> \param nelectron_total ...
!> \param nelectron_spin ...
!> \param efield ...
!> \param neighbor_list_id ...
!> \param linres_control ...
!> \param xas_env ...
!> \param virial ...
!> \param cp_ddapc_env ...
!> \param cp_ddapc_ewald ...
!> \param outer_scf_history ...
!> \param outer_scf_ihistory ...
!> \param ep_qs_env ...
!> \param x_data ...
!> \param et_coupling ...
!> \param dftb_potential ...
!> \param results ...
!> \param scp_env ...
!> \param se_taper ...
!> \param se_store_int_env ...
!> \param se_nddo_mpole ...
!> \param se_nonbond_env ...
!> \param admm_env ...
!> \param lri_env ...
!> \param lri_density ...
!> \param hfx_ri_env ...
!> \param dispersion_env ...
!> \param vee ...
!> \param rho_external ...
!> \param external_vxc ...
!> \param mask ...
!> \param mp2_env ...
!> \param kg_env ...
!> \param WannierCentres ...
!> \param ga_env ...
!> \param atprop ...
!> \param ls_scf_env ...
!> \param v_hartree_rspace ...
!> \param s_mstruct_changed ...
!> \param rho_changed ...
!> \param potential_changed ...
!> \param forces_up_to_date ...
!> \param error ...
!> \date    23.01.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_qs_env(qs_env,atomic_kind_set,qs_kind_set,cell,super_cell,cell_ref,use_ref_cell,kpoints,&
       dft_control,mos,mos_aux_fit,sab_orb,sab_aux_fit,sab_aux_fit_asymm, sab_aux_fit_vs_orb,&
       sab_all,qmmm,qmmm_periodic,sac_ae,sac_ppl,sap_ppnl,sab_vdw,sab_scp,sap_oce,sab_lrc,&
       sab_se,sab_tbe,sab_core,sab_almo,particle_set,energy,force,harris_env,&
       matrix_h,matrix_ks,matrix_vxc,run_rtp,rtp,&
       matrix_ks_aux_fit,matrix_s,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb,matrix_w,matrix_w_mp2,matrix_p_mp2,&
       gamma_matrix, rho, rho_aux_fit,rho_buffer, rho_xc, pw_env, ewald_env, ewald_pw, &
       mpools,mpools_aux_fit,input, para_env, blacs_env, scf_control,rel_control,kinetic,qs_charges,&
       vppl, rho_core, rho_nlcc, rho_nlcc_g, ks_env, ks_qmmm_env, wf_history,scf_env,use_harris, id_nr, local_particles,&
       local_molecules, distribution_2d, dbcsr_dist, molecule_kind_set,molecule_set,subsys,cp_subsys,oce,rho_atom_set,&
       task_list,task_list_zmp,task_list_aux_fit,task_list_soft,rho0_atom_set,rho0_mpole,rhoz_set,ecoul_1c,&
       rho0_s_rs,rho0_s_gs,do_kpoints,has_unit_metric,requires_mo_derivs,mo_derivs,&
       mo_derivs_aux_fit, mo_loc_history,nkind,natom,nelectron_total,nelectron_spin,efield,&
       neighbor_list_id,linres_control,xas_env,virial,cp_ddapc_env,cp_ddapc_ewald,&
       outer_scf_history,outer_scf_ihistory,ep_qs_env,x_data,et_coupling,dftb_potential,results,&
       scp_env,se_taper,se_store_int_env,se_nddo_mpole,se_nonbond_env, admm_env, &
       lri_env,lri_density,hfx_ri_env,dispersion_env, vee,rho_external,external_vxc,mask,&
       mp2_env,kg_env,WannierCentres,ga_env,atprop,ls_scf_env,v_hartree_rspace,&
       s_mstruct_changed,rho_changed,potential_changed,forces_up_to_date,&
       error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(qs_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: qs_kind_set
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, super_cell, cell_ref
    LOGICAL, OPTIONAL                        :: use_ref_cell
    TYPE(kpoint_type), OPTIONAL, POINTER     :: kpoints
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
    TYPE(qs_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: cp_subsys
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
    LOGICAL, OPTIONAL                        :: do_kpoints, has_unit_metric, &
                                                requires_mo_derivs
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mo_derivs
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mo_derivs_aux_fit, &
                                                mo_loc_history
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
    TYPE(pw_type), OPTIONAL, POINTER         :: v_hartree_rspace
    LOGICAL, OPTIONAL                        :: s_mstruct_changed, &
                                                rho_changed, &
                                                potential_changed, &
                                                forces_up_to_date
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_qs_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(rho0_mpole_type), POINTER           :: rho0_m

    failure=.FALSE.
    NULLIFY(rho0_m)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qs_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env%ks_env),cp_failure_level,routineP,error,failure)

    IF (PRESENT(nelectron_total)) nelectron_total=qs_env%nelectron_total
    IF (PRESENT(nelectron_spin))  nelectron_spin=qs_env%nelectron_spin
    IF (PRESENT(outer_scf_history)) outer_scf_history=>qs_env%outer_scf_history
    IF (PRESENT(outer_scf_ihistory)) outer_scf_ihistory=qs_env%outer_scf_ihistory
    IF (PRESENT(ga_env)) ga_env => qs_env%ga_env
    IF (PRESENT(mp2_env)) mp2_env => qs_env%mp2_env
    IF (PRESENT(kg_env)) kg_env => qs_env%kg_env
    IF (PRESENT(super_cell)) super_cell => qs_env%super_cell
    IF (PRESENT(qmmm)) qmmm = qs_env%qmmm
    IF (PRESENT(qmmm_periodic)) qmmm_periodic = qs_env%qmmm_periodic
    IF (PRESENT(mos)) mos => qs_env%mos
    IF (PRESENT(mos_aux_fit)) mos_aux_fit => qs_env%mos_aux_fit
    IF (PRESENT(use_harris)) use_harris = qs_env%use_harris
    IF (PRESENT(harris_env)) harris_env => qs_env%harris_env
    IF (PRESENT(ewald_env)) ewald_env => qs_env%ewald_env
    IF (PRESENT(ewald_pw)) ewald_pw => qs_env%ewald_pw
    IF (PRESENT(mpools)) mpools => qs_env%mpools
    IF (PRESENT(mpools_aux_fit)) mpools_aux_fit => qs_env%mpools_aux_fit
    IF (PRESENT(blacs_env)) blacs_env => qs_env%blacs_env
    IF (PRESENT(scf_control)) scf_control => qs_env%scf_control
    IF (PRESENT(rel_control)) rel_control => qs_env%rel_control
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
    IF (PRESENT(para_env)) para_env => qs_env%para_env
    IF (PRESENT(oce)) oce =>  qs_env%oce
    IF (PRESENT(requires_mo_derivs)) requires_mo_derivs = qs_env%requires_mo_derivs
    IF (PRESENT(has_unit_metric)) has_unit_metric = qs_env%has_unit_metric
    IF (PRESENT(mo_derivs)) mo_derivs => qs_env%mo_derivs
    IF (PRESENT(mo_derivs_aux_fit)) mo_derivs_aux_fit => qs_env%mo_derivs_aux_fit
    IF (PRESENT(mo_loc_history)) mo_loc_history => qs_env%mo_loc_history
    IF (PRESENT(linres_control)) linres_control => qs_env%linres_control
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
    IF (PRESENT(input)) input => qs_env%input
    IF (PRESENT(cp_ddapc_env)) cp_ddapc_env => qs_env%cp_ddapc_env
    IF (PRESENT(cp_ddapc_ewald)) cp_ddapc_ewald => qs_env%cp_ddapc_ewald
    IF (PRESENT(ep_qs_env)) ep_qs_env => qs_env%ep_qs_env
    IF (PRESENT(x_data)) x_data => qs_env%x_data
    IF (PRESENT(et_coupling)) et_coupling => qs_env%et_coupling
    IF (PRESENT(dftb_potential)) dftb_potential => qs_env%dftb_potential
    IF (PRESENT(scp_env)) scp_env => qs_env%scp_env
    IF (PRESENT(efield)) efield => qs_env%efield
    IF (PRESENT(task_list)) task_list => qs_env%task_list
    IF (PRESENT(task_list_aux_fit)) task_list_aux_fit => qs_env%task_list_aux_fit
    IF (PRESENT(task_list_soft)) task_list_soft => qs_env%task_list_soft
    IF (PRESENT(WannierCentres)) WannierCentres => qs_env%WannierCentres

    CALL get_ks_env(qs_env%ks_env,&
                    v_hartree_rspace=v_hartree_rspace,&
                    s_mstruct_changed=s_mstruct_changed,&
                    rho_changed=rho_changed,&
                    potential_changed=potential_changed,&
                    forces_up_to_date=forces_up_to_date,&
                    matrix_h=matrix_h,&
                    matrix_ks=matrix_ks,&
                    matrix_vxc=matrix_vxc,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    kinetic=kinetic,&
                    matrix_s=matrix_s,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    matrix_w=matrix_w,&
                    matrix_w_mp2=matrix_w_mp2,&
                    matrix_p_mp2=matrix_p_mp2,&
                    gamma_matrix=gamma_matrix,&
                    rho=rho,&
                    rho_buffer=rho_buffer,&
                    rho_xc=rho_xc,&
                    rho_aux_fit=rho_aux_fit,&
                    rho_core=rho_core,&
                    rho_nlcc=rho_nlcc,&
                    rho_nlcc_g=rho_nlcc_g,&
                    vppl=vppl,&
                    vee=vee,&
                    neighbor_list_id=neighbor_list_id,&
                    sab_orb=sab_orb,&
                    sab_aux_fit=sab_aux_fit,&
                    sab_aux_fit_asymm=sab_aux_fit_asymm,&
                    sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,&
                    sab_all=sab_all,&
                    sab_scp=sab_scp,&
                    sab_vdw=sab_vdw,&
                    sac_ae=sac_ae,&
                    sac_ppl=sac_ppl,&
                    sap_ppnl=sap_ppnl,&
                    sap_oce=sap_oce,&
                    sab_se=sab_se,&
                    sab_lrc=sab_lrc,&
                    sab_tbe=sab_tbe,&
                    sab_core=sab_core,&
                    sab_almo=sab_almo,&
                    kpoints=kpoints,&
                    do_kpoints=do_kpoints,&
                    local_molecules=local_molecules,&
                    local_particles=local_particles,&
                    atprop=atprop,&
                    virial=virial,&
                    results=results,&
                    cell=cell,&
                    cell_ref=cell_ref,&
                    use_ref_cell=use_ref_cell,&
                    energy=energy,&
                    force=force,&
                    qs_kind_set=qs_kind_set,&
                    subsys=subsys,&
                    cp_subsys=cp_subsys,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    molecule_kind_set=molecule_kind_set,&
                    molecule_set=molecule_set,&
                    natom=natom,&
                    nkind=nkind,&
                    dft_control=dft_control,&
                    dbcsr_dist=dbcsr_dist,&
                    distribution_2d=distribution_2d,&
                    pw_env=pw_env,&
                    error=error)

  END SUBROUTINE get_qs_env

! *****************************************************************************
!> \brief  Initialise the QUICKSTEP environment.
!> \param qs_env ...
!> \param para_env ...
!> \param globenv ...
!> \param error ...
!> \date    25.01.2002
!> \author  MK
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
    NULLIFY (qs_env%super_cell)
    NULLIFY (qs_env%mos)
    NULLIFY (qs_env%mos_aux_fit)
    NULLIFY (qs_env%harris_env)
    NULLIFY (qs_env%mpools)
    NULLIFY (qs_env%mpools_aux_fit)
    NULLIFY (qs_env%ewald_env)
    NULLIFY (qs_env%ewald_pw)
    NULLIFY (qs_env%scf_control)
    NULLIFY (qs_env%rel_control)
    NULLIFY (qs_env%qs_charges)
    ! ZMP initializing arrays
    NULLIFY (qs_env%rho_external)
    NULLIFY (qs_env%external_vxc)
    NULLIFY (qs_env%mask)

    NULLIFY (qs_env%ks_env)
    NULLIFY (qs_env%ks_qmmm_env)
    NULLIFY (qs_env%wf_history)
    NULLIFY (qs_env%scf_env)
    NULLIFY (qs_env%blacs_env)
    NULLIFY (qs_env%oce)
    NULLIFY (qs_env%local_rho_set)
    NULLIFY (qs_env%hartree_local)
    NULLIFY (qs_env%input)
    NULLIFY (qs_env%linres_control)
    NULLIFY (qs_env%xas_env)
    NULLIFY (qs_env%cp_ddapc_env)
    NULLIFY (qs_env%cp_ddapc_ewald)
    NULLIFY (qs_env%outer_scf_history)
    NULLIFY (qs_env%ep_qs_env)
    NULLIFY (qs_env%x_data)
    NULLIFY (qs_env%et_coupling)
    NULLIFY (qs_env%dftb_potential)
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

    qs_env%outer_scf_ihistory=0
    qs_env%broyden_adaptive_sigma=-1.0_dp
    qs_env%nelectron_total=-1
    qs_env%nelectron_spin=-1

    CALL local_rho_set_create(qs_env%local_rho_set, error=error)
    CALL hartree_local_create(qs_env%hartree_local, error=error)
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
    NULLIFY(qs_env%mo_derivs)
    NULLIFY(qs_env%mo_derivs_aux_fit)
    NULLIFY(qs_env%mo_loc_history)

  END SUBROUTINE init_qs_env

! *****************************************************************************
!> \brief   Set the QUICKSTEP environment.
!> \param qs_env ...
!> \param super_cell ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param qmmm ...
!> \param qmmm_periodic ...
!> \param harris_env ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param mpools ...
!> \param mpools_aux_fit ...
!> \param rho_external ...
!> \param external_vxc ...
!> \param mask ...
!> \param kinetic ...
!> \param scf_control ...
!> \param rel_control ...
!> \param qs_charges ...
!> \param ks_env ...
!> \param ks_qmmm_env ...
!> \param wf_history ...
!> \param scf_env ...
!> \param use_harris ...
!> \param task_list ...
!> \param task_list_aux_fit ...
!> \param task_list_soft ...
!> \param task_list_zmp ...
!> \param input ...
!> \param oce ...
!> \param rho_atom_set ...
!> \param rho0_atom_set ...
!> \param rho0_mpole ...
!> \param run_rtp ...
!> \param rtp ...
!> \param rhoz_set ...
!> \param rhoz_tot ...
!> \param ecoul_1c ...
!> \param has_unit_metric ...
!> \param requires_mo_derivs ...
!> \param mo_derivs ...
!> \param mo_derivs_aux_fit ...
!> \param mo_loc_history ...
!> \param nelectron_total ...
!> \param nelectron_spin ...
!> \param efield ...
!> \param linres_control ...
!> \param xas_env ...
!> \param cp_ddapc_env ...
!> \param cp_ddapc_ewald ...
!> \param outer_scf_history ...
!> \param outer_scf_ihistory ...
!> \param ep_qs_env ...
!> \param x_data ...
!> \param et_coupling ...
!> \param dftb_potential ...
!> \param scp_env ...
!> \param se_taper ...
!> \param se_store_int_env ...
!> \param se_nddo_mpole ...
!> \param se_nonbond_env ...
!> \param admm_env ...
!> \param ls_scf_env ...
!> \param lri_env ...
!> \param lri_density ...
!> \param hfx_ri_env ...
!> \param dispersion_env ...
!> \param mp2_env ...
!> \param kg_env ...
!> \param WannierCentres ...
!> \param ga_env ...
!> \param error ...
!> \date    23.01.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_qs_env(qs_env,super_cell,&
       mos,mos_aux_fit,qmmm,qmmm_periodic,harris_env,&
       ewald_env,ewald_pw,mpools,mpools_aux_fit,&
       rho_external,external_vxc,mask,&
       kinetic,scf_control,rel_control,qs_charges,ks_env,&
       ks_qmmm_env,wf_history,scf_env,use_harris,task_list,task_list_aux_fit,task_list_soft,&
       task_list_zmp,&
       input,oce,rho_atom_set,rho0_atom_set,rho0_mpole,run_rtp,rtp,&
       rhoz_set,rhoz_tot,ecoul_1c,has_unit_metric,requires_mo_derivs,mo_derivs,&
       mo_derivs_aux_fit, mo_loc_history, nelectron_total, nelectron_spin, efield,&
       linres_control,xas_env,cp_ddapc_env,cp_ddapc_ewald,&
       outer_scf_history,outer_scf_ihistory,ep_qs_env,x_data,et_coupling,dftb_potential,&
       scp_env,se_taper,se_store_int_env,se_nddo_mpole,se_nonbond_env,admm_env,ls_scf_env,&
       lri_env,lri_density,hfx_ri_env,dispersion_env,mp2_env,kg_env,WannierCentres,ga_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cell_type), OPTIONAL, POINTER       :: super_cell
    TYPE(mo_set_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: mos, mos_aux_fit
    LOGICAL, OPTIONAL                        :: qmmm, qmmm_periodic
    TYPE(harris_env_type), OPTIONAL, POINTER :: harris_env
    TYPE(ewald_environment_type), OPTIONAL, &
      POINTER                                :: ewald_env
    TYPE(ewald_pw_type), OPTIONAL, POINTER   :: ewald_pw
    TYPE(qs_matrix_pools_type), OPTIONAL, &
      POINTER                                :: mpools, mpools_aux_fit
    TYPE(qs_rho_type), OPTIONAL, POINTER     :: rho_external
    TYPE(pw_p_type), OPTIONAL, POINTER       :: external_vxc, mask
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: kinetic
    TYPE(scf_control_type), OPTIONAL, &
      POINTER                                :: scf_control
    TYPE(rel_control_type), OPTIONAL, &
      POINTER                                :: rel_control
    TYPE(qs_charges_type), OPTIONAL, POINTER :: qs_charges
    TYPE(qs_ks_env_type), OPTIONAL, POINTER  :: ks_env
    TYPE(qs_ks_qmmm_env_type), OPTIONAL, &
      POINTER                                :: ks_qmmm_env
    TYPE(qs_wf_history_type), OPTIONAL, &
      POINTER                                :: wf_history
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    LOGICAL, OPTIONAL                        :: use_harris
    TYPE(task_list_type), OPTIONAL, POINTER  :: task_list, task_list_aux_fit, &
                                                task_list_soft, task_list_zmp
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
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
    TYPE(linres_control_type), OPTIONAL, &
      POINTER                                :: linres_control
    TYPE(xas_environment_type), OPTIONAL, &
      POINTER                                :: xas_env
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
    IF (PRESENT(nelectron_total)) qs_env%nelectron_total=nelectron_total
    IF (PRESENT(mp2_env)) qs_env%mp2_env=>mp2_env
    IF (PRESENT(kg_env)) qs_env%kg_env => kg_env
    IF (PRESENT(nelectron_spin))  qs_env%nelectron_spin=nelectron_spin
    IF (PRESENT(super_cell)) THEN
       CALL cell_retain(super_cell, error=error)
       CALL cell_release(qs_env%super_cell,error=error)
       qs_env%super_cell => super_cell
    END IF
    !
    IF (PRESENT(qmmm)) qs_env%qmmm = qmmm
    IF (PRESENT(qmmm_periodic)) qs_env%qmmm_periodic = qmmm_periodic
    IF (PRESENT(mos)) qs_env%mos => mos
    IF (PRESENT(mos_aux_fit)) qs_env%mos_aux_fit => mos_aux_fit
    IF (PRESENT(ls_scf_env)) qs_env%ls_scf_env => ls_scf_env
    ! if intels checking (-C) complains here, you have rediscovered a bug in the intel
    ! compiler (present in at least 10.0.025). A testcase has been submitted to intel.
    IF (PRESENT(use_harris)) qs_env%use_harris = use_harris
    IF (PRESENT(harris_env)) qs_env%harris_env => harris_env
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
!> \param globenv ...
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
          CALL cell_release(qs_env%super_cell,error=error)
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
          IF (ASSOCIATED(qs_env%rtp))THEN
             CALL rt_prop_release(qs_env%rtp,error)
             DEALLOCATE (qs_env%rtp)
          END IF
          IF (ASSOCIATED(qs_env%outer_scf_history)) THEN
             DEALLOCATE(qs_env%outer_scf_history,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             qs_env%outer_scf_ihistory=0
          ENDIF
          IF (ASSOCIATED(qs_env%harris_env)) THEN
             CALL harris_env_release(qs_env%harris_env, error=error)
          END IF
          IF (ASSOCIATED(qs_env%oce)) CALL deallocate_oce_set(qs_env%oce,error)
          IF (ASSOCIATED(qs_env%local_rho_set)) THEN
             CALL local_rho_set_release(qs_env%local_rho_set,error=error)
          END IF
          IF (ASSOCIATED(qs_env%hartree_local)) THEN
             CALL hartree_local_release(qs_env%hartree_local,error=error)
          END IF
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
          CALL ewald_env_release(qs_env%ewald_env, error=error)
          CALL ewald_pw_release(qs_env%ewald_pw, error=error)
          CALL cp_para_env_release(qs_env%para_env,error=error)
          CALL cp_blacs_env_release(qs_env%blacs_env, error=error)
          IF (ASSOCIATED(qs_env%image_matrix)) THEN
             DEALLOCATE(qs_env%image_matrix)
          ENDIF
          IF (ASSOCIATED(qs_env%ipiv)) THEN
             DEALLOCATE(qs_env%ipiv)
          ENDIF
          IF (ASSOCIATED(qs_env%image_coeff)) THEN
             DEALLOCATE(qs_env%image_coeff)
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
          CALL section_vals_release(qs_env%input,error=error)
          CALL cp_ddapc_release(qs_env%cp_ddapc_env, error=error)
          CALL cp_ddapc_ewald_release(qs_env%cp_ddapc_ewald, error=error)
          CALL ep_qs_release(qs_env%ep_qs_env,error=error)
          CALL efield_berry_release(qs_env%efield,error=error)
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
          ! now we are ready to deallocate the full structure
          DEALLOCATE(qs_env, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(qs_env)
  END SUBROUTINE qs_env_release

END MODULE qs_environment_types
