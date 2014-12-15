!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for hfx and admm methods
!>
!>
!> \par History
!>     refactoring 03-2011 [MI]
!> \author MI
! *****************************************************************************
MODULE hfx_admm_utils
  USE admm_dm_types,                   ONLY: admm_dm_create,&
                                             admm_dm_type
  USE admm_methods,                    ONLY: scale_dm
  USE admm_types,                      ONLY: admm_env_create,&
                                             admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_p_type,&
                                             cp_dbcsr_set,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE hfx_derivatives,                 ONLY: derivatives_four_center
  USE hfx_energy_potential,            ONLY: integrate_four_center
  USE hfx_ri_methods,                  ONLY: hfx_ri_energy_potential
  USE input_constants,                 ONLY: &
       do_admm_aux_exch_func_bee, do_admm_aux_exch_func_default, &
       do_admm_aux_exch_func_none, do_admm_aux_exch_func_opt, &
       do_admm_aux_exch_func_pbex, do_hfx_potential_coulomb, &
       do_hfx_potential_short, do_hfx_potential_truncated, &
       xc_funct_no_shortcut
  USE input_section_types,             ONLY: section_vals_duplicate,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_transfer
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_create,&
                                             pw_p_type,&
                                             pw_release
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                             set_ks_env
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE rt_propagation_types,            ONLY: rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
  USE xc_adiabatic_utils,              ONLY: rescale_xc_potential
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: hfx_ks_matrix, hfx_admm_init, create_admm_xc_section

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_admm_utils'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
 SUBROUTINE  hfx_admm_init(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_admm_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n_rep_hf, natoms, &
                                                nspins
    LOGICAL                                  :: do_hfx, failure, &
                                                s_mstruct_changed
    TYPE(admm_dm_type), POINTER              :: admm_dm
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux_fit
    TYPE(section_vals_type), POINTER         :: hfx_sections, input, &
                                                xc_section

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(admm_env, admm_dm, hfx_sections, mos, mos_aux_fit, para_env, &
         particle_set, xc_section, matrix_s_aux_fit, &
         matrix_s_aux_fit_vs_orb, rho, rho_aux_fit, ks_env, dft_control, &
         input)

    CALL get_qs_env(qs_env,&
                    mos_aux_fit=mos_aux_fit,&
                    mos=mos,&
                    admm_env=admm_env,&
                    admm_dm=admm_dm,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    para_env=para_env,&
                    s_mstruct_changed=s_mstruct_changed,&
                    rho=rho,&
                    rho_aux_fit=rho_aux_fit,&
                    ks_env=ks_env,&
                    dft_control=dft_control,&
                    input=input,&
                    error=error)

    nspins=dft_control%nspins

    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
    !! ** ADMM can only be used with HFX
    IF ( .NOT. do_hfx ) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "Wavefunction fitting requested without Hartree-Fock."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END IF

    ! ** Method only implemented for GPW
    IF( dft_control%qs_control%gapw ) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
                               message="ADMM only implemented for GPW", &
                               error=error, error_level=cp_failure_level)
    END IF


    CALL section_vals_get(hfx_sections,n_repetition=n_rep_hf,error=error)
    IF ( n_rep_hf > 1 ) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "ADMM can handle only one HF section."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END IF

    IF (.NOT. ASSOCIATED(admm_env)) THEN
      ! setup admm environment
      CALL get_qs_env(qs_env, input=input, particle_set=particle_set, error=error)
      natoms = SIZE(particle_set,1)
      CALL admm_env_create(admm_env, dft_control%admm_control, mos, mos_aux_fit, &
                           para_env, natoms, error)
      CALL set_qs_env(qs_env,admm_env=admm_env,error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env=qs_env, xc_section=xc_section, &
                                  admm_env=admm_env, error=error)
    ENDIF

    IF(dft_control%do_admm_dm .AND. .NOT. ASSOCIATED(admm_dm)) THEN
        CALL admm_dm_create(admm_dm, dft_control%admm_control, nspins=nspins, natoms=natoms)
        CALL set_ks_env(ks_env, admm_dm=admm_dm, error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE hfx_admm_init

! *****************************************************************************
!> \brief Add the hfx contributions to the Hamiltonian
!>
!> \param ks_env ...
!> \param qs_env ...
!> \param matrix_ks ...
!> \param rho ...
!> \param energy ...
!> \param calculate_forces ...
!> \param just_energy ...
!> \param v_rspace_new ...
!> \param v_tau_rspace ...
!> \param error ...
!> \par History
!>     refactoring 03-2011 [MI]
! *****************************************************************************

  SUBROUTINE hfx_ks_matrix(ks_env,qs_env,matrix_ks,rho,energy,calculate_forces,&
              just_energy, v_rspace_new,v_tau_rspace,error)

    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces, just_energy
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, v_tau_rspace
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, irep, ispin, &
                                                mspin, n_rep_hf, nspins, stat
    LOGICAL :: distribute_fock_matrix, do_adiabatic_rescaling, do_hfx_ri, &
      failure, hfx_treat_lsd_in_core, s_mstruct_changed, use_virial
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: hf_energy
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_h, &
      matrix_ks_aux_fit, matrix_ks_aux_fit_hfx, matrix_ks_aux_fit_im, &
      matrix_ks_im, matrix_ks_orb, rho_ao_orb
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_rho_type), POINTER               :: rho_orb
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(section_vals_type), POINTER         :: adiabatic_rescaling_section, &
                                                hfx_ri_section, hfx_sections, &
                                                input
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN, handle)

    failure = .FALSE.
    NULLIFY(auxbas_pw_pool, dft_control, force, hfx_sections, hfx_ri_section, input, &
            para_env, poisson_env, pw_env, virial,&
            matrix_ks_im, matrix_ks_orb, rho_ao_orb,&
            matrix_h, matrix_ks_aux_fit,matrix_ks_aux_fit_im, matrix_ks_aux_fit_hfx)

    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,&
         input=input,&
         matrix_h=matrix_h,&
         para_env=para_env,&
         pw_env=pw_env,&
         virial=virial,&
         matrix_ks_im=matrix_ks_im,&
         matrix_ks_aux_fit=matrix_ks_aux_fit,&
         matrix_ks_aux_fit_im=matrix_ks_aux_fit_im,&
         matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx,&
         s_mstruct_changed=s_mstruct_changed,&
         error=error)

    nspins=dft_control%nspins

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_sections,n_repetition=n_rep_hf,error=error)
    hfx_ri_section => section_vals_get_subs_vals(hfx_sections,"HFX_RI",error=error)
    CALL section_vals_get(hfx_ri_section,explicit=do_hfx_ri,error=error)
    CALL section_vals_val_get(hfx_sections, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core,&
                                i_rep_section=1,error=error)
    adiabatic_rescaling_section => section_vals_get_subs_vals(input,"DFT%XC%ADIABATIC_RESCALING",error=error)
    CALL section_vals_get(adiabatic_rescaling_section,explicit=do_adiabatic_rescaling,error=error)



    ! *** Initialize the auxiliary ks matrix to zero if required
    IF( dft_control%do_admm ) THEN
      DO ispin = 1,nspins
         CALL cp_dbcsr_set(matrix_ks_aux_fit(ispin)%matrix,0.0_dp,error=error)
         CALL cp_dbcsr_set(matrix_ks_aux_fit_hfx(ispin)%matrix,0.0_dp,error=error)
      END DO
    END IF
    DO ispin = 1,nspins
       CALL cp_dbcsr_set(matrix_ks(ispin)%matrix,0.0_dp,error=error)
    END DO

    CALL section_vals_get(hfx_sections,n_repetition=n_rep_hf,error=error)
    !! set hf exchange energy to zero
    energy%ex = 0.0_dp

    IF ( do_hfx_ri ) THEN
       IF( calculate_forces ) THEN
          !! initalize force array to zero
          CALL get_qs_env(qs_env=qs_env, force=force, error=error)
          DO ikind = 1,SIZE(force)
             force(ikind)%hfx_ri(:,:) = 0.0_dp
          END DO
       END IF
       CALL hfx_ri_energy_potential (qs_env,calculate_forces,error)
    ELSE
       IF( calculate_forces ) THEN
          !! initalize force array to zero
          CALL get_qs_env(qs_env=qs_env, force=force, error=error)
          DO ikind = 1,SIZE(force)
             force(ikind)%fock_4c(:,:) = 0.0_dp
          END DO
       END IF
       ALLOCATE(hf_energy(n_rep_hf), STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO irep = 1,n_rep_hf

          IF( do_adiabatic_rescaling .AND. hfx_treat_lsd_in_core)&
             CALL cp_unimplemented_error(fromWhere=routineP, &
                     message="HFX_TREAT_LSD_IN_CORE not implemented for adiabatically rescaled hybrids",&
                     error=error, error_level=cp_failure_level)
          ! everything is calulated with adiabatic rescaling but the potential is not added in a first step
          distribute_fock_matrix=.NOT.do_adiabatic_rescaling

          mspin=1
          IF ( hfx_treat_lsd_in_core ) mspin=nspins
          
          ! fetch the correct matrices for normal HFX or ADMM
          IF( dft_control%do_admm ) THEN
             CALL get_qs_env(qs_env=qs_env,matrix_ks_aux_fit=matrix_ks_orb,&
                             rho_aux_fit=rho_orb,error=error)
          ELSE
             CALL get_qs_env(qs_env=qs_env,matrix_ks=matrix_ks_orb,rho=rho_orb,error=error)
          END IF
          CALL qs_rho_get(rho_struct=rho_orb,rho_ao=rho_ao_orb,error=error)
          ! Finally the real hfx calulation
          DO ispin = 1,mspin
             CALL integrate_four_center(qs_env, matrix_ks_orb ,energy,rho_ao_orb,hfx_sections,&
                                        para_env, s_mstruct_changed, irep, distribute_fock_matrix,&
                                        ispin=ispin, error=error)
          END DO
       
          IF( calculate_forces .AND. .NOT. do_adiabatic_rescaling ) THEN
             !Scale auxiliary density matrix for ADMMP (see Merlot2014) with gsi(ispin) to scale force
             IF( dft_control%do_admm ) THEN
               CALL scale_dm(qs_env, rho_ao_orb, scale_back=.FALSE., error=error) 
             END IF
             CALL derivatives_four_center(qs_env, rho_ao_orb, hfx_sections, &
                                          para_env, irep, use_virial, error=error)
             !Scale auxiliary density matrix for ADMMP back with 1/gsi(ispin) 
             IF( dft_control%do_admm ) THEN
               CALL scale_dm(qs_env, rho_ao_orb, scale_back=.TRUE., error=error)
             END IF

          END IF

          !! If required, the calculation of the forces will be done later with adiabatic rescaling
          IF( do_adiabatic_rescaling ) hf_energy(irep) = energy%ex

          ! special case RTP/EMD we have a full complex density and HFX has a contrinution from the imaginary part
          IF( qs_env%run_rtp)THEN

             CALL get_qs_env(qs_env=qs_env,rtp=rtp,error=error)
             DO ispin = 1,nspins
                CALL cp_dbcsr_set(matrix_ks_im(ispin)%matrix,0.0_dp,error=error)
             END DO
             IF( dft_control%do_admm ) THEN
                matrix_ks_orb => matrix_ks_aux_fit_im
                DO ispin = 1,nspins
                   CALL cp_dbcsr_set(matrix_ks_aux_fit_im(ispin)%matrix,0.0_dp,error=error)
                END DO
             ELSE
                matrix_ks_orb => matrix_ks_im
             END IF

             CALL qs_rho_get(rho_orb, rho_ao_im=rho_ao_orb, error=error)

             DO ispin = 1,mspin
                CALL integrate_four_center(qs_env, matrix_ks_orb ,energy,rho_ao_orb,hfx_sections,&
                                           para_env, .FALSE., irep, distribute_fock_matrix,&
                                           ispin=ispin, error=error)
             END DO
   
             IF( calculate_forces .AND. .NOT. do_adiabatic_rescaling ) THEN
                CALL derivatives_four_center(qs_env, rho_ao_orb, hfx_sections, &
                                             para_env, irep, use_virial, error=error)
             END IF
   
             !! If required, the calculation of the forces will be done later with adiabatic rescaling
             IF( do_adiabatic_rescaling ) hf_energy(irep) = energy%ex
          END IF

          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                       poisson_env=poisson_env,error=error)
          CALL pw_hfx(qs_env,energy,hfx_sections,poisson_env,auxbas_pw_pool,irep,error)
         
       END DO
    END IF

    ! *** Add Core-Hamiltonian-Matrix ***
    DO ispin=1,nspins
         CALL cp_dbcsr_add(matrix_ks(ispin)%matrix,matrix_h(1)%matrix,&
              1.0_dp,1.0_dp,error=error)
    END DO
    IF (use_virial .AND. calculate_forces) THEN
         virial%pv_virial = virial%pv_virial - virial%pv_fock_4c
         virial%pv_calculate = .FALSE.
    ENDIF

    !! If we perform adiabatic rescaling we are now able to rescale the xc-potential
    IF( do_adiabatic_rescaling ) THEN
       CALL rescale_xc_potential(qs_env,matrix_ks,rho,energy,v_rspace_new,v_tau_rspace,&
            hf_energy, just_energy,calculate_forces,use_virial,error=error)
    END IF   ! do_adiabatic_rescaling

    CALL timestop(handle)

  END SUBROUTINE  hfx_ks_matrix

! *****************************************************************************
!> \brief computes the Hartree-Fock energy brute force in a pw basis
!> \param qs_env ...
!> \param energy ...
!> \param hfx_section ...
!> \param poisson_env ...
!> \param auxbas_pw_pool ...
!> \param irep ...
!> \param error ...
!> \par History
!>      12.2007 created [Joost VandeVondele]
!> \note
!>      only computes the HFX energy, no derivatives as yet
! *****************************************************************************
  SUBROUTINE pw_hfx(qs_env,energy,hfx_section,poisson_env,auxbas_pw_pool,irep,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(section_vals_type), POINTER         :: hfx_section
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    INTEGER                                  :: irep
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(*), PARAMETER :: routineN = 'pw_hfx', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blocksize, handle, iloc, &
                                                iorb, iorb_block, ispin, iw, &
                                                jloc, jorb, jorb_block, norb
    LOGICAL                                  :: do_pw_hfx, failure
    REAL(KIND=dp)                            :: exchange_energy, fraction, &
                                                pair_energy, scaling
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_b
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: pot_g, rho_g, rho_r
    TYPE(pw_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: rho_i, rho_j
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)

    CALL section_vals_val_get(hfx_section, "PW_HFX", l_val=do_pw_hfx, i_rep_section=irep, error=error)

    IF (do_pw_hfx) THEN
       CALL section_vals_val_get(hfx_section, "FRACTION", r_val=fraction, error=error)
       CALL section_vals_val_get(hfx_section, "PW_HFX_BLOCKSIZE", i_val=blocksize, error=error)

       CALL get_qs_env(qs_env,mos=mo_array,pw_env=pw_env, dft_control=dft_control, &
                       cell=cell, particle_set=particle_set, &
                       atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, error=error)

       ! limit the blocksize by the number of orbitals
       CALL get_mo_set(mo_set=mo_array(1)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,ncol_global=norb,error=error)
       blocksize=MAX(1,MIN(blocksize,norb))

       CALL pw_pool_create_pw(auxbas_pw_pool,rho_r%pw,&
                               use_data=REALDATA3D,&
                               in_space=REALSPACE,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,rho_g%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,pot_g%pw,&
                               use_data=COMPLEXDATA1D,&
                               in_space=RECIPROCALSPACE,error=error)

       ALLOCATE(rho_i(blocksize))
       ALLOCATE(rho_j(blocksize))

       DO iorb_block=1,blocksize
          NULLIFY(rho_i(iorb_block)%pw)
          CALL pw_create(rho_i(iorb_block)%pw,rho_r%pw%pw_grid,&
                                  use_data=REALDATA3D,&
                                  in_space=REALSPACE,error=error)
          NULLIFY(rho_j(iorb_block)%pw)
          CALL pw_create(rho_j(iorb_block)%pw,rho_r%pw%pw_grid,&
                                  use_data=REALDATA3D,&
                                  in_space=REALSPACE,error=error)
       ENDDO

       exchange_energy = 0.0_dp

       DO ispin=1,SIZE(mo_array)
          CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff,mo_coeff_b=mo_coeff_b)

          IF(mo_array(ispin)%mo_set%use_mo_coeff_b) THEN!fm->dbcsr
             CALL copy_dbcsr_to_fm(mo_coeff_b,mo_coeff,error=error)!fm->dbcsr
          ENDIF!fm->dbcsr

          CALL cp_fm_get_info(mo_coeff,ncol_global=norb,error=error)

          DO iorb_block=1,norb,blocksize

             DO iorb=iorb_block,MIN(iorb_block+blocksize-1,norb)

                 iloc=iorb-iorb_block+1
                 CALL calculate_wavefunction(mo_coeff,iorb,rho_i(iloc),rho_g, &
                          atomic_kind_set,qs_kind_set, cell,dft_control,particle_set, &
                          pw_env,error=error)

             ENDDO

             DO jorb_block=iorb_block,norb,blocksize

                DO jorb=jorb_block,MIN(jorb_block+blocksize-1,norb)

                   jloc=jorb-jorb_block+1
                   CALL calculate_wavefunction(mo_coeff,jorb,rho_j(jloc),rho_g, &
                            atomic_kind_set,qs_kind_set, cell,dft_control,particle_set, &
                            pw_env,error=error)

                ENDDO

                DO iorb=iorb_block,MIN(iorb_block+blocksize-1,norb)
                   iloc=iorb-iorb_block+1
                   DO jorb=jorb_block,MIN(jorb_block+blocksize-1,norb)
                      jloc=jorb-jorb_block+1
                      IF (jorb<iorb) CYCLE

                      ! compute the pair density
                      rho_r%pw%cr3d = rho_i(iloc)%pw%cr3d * rho_j(jloc)%pw%cr3d

                      ! go the g-space and compute hartree energy
                      CALL pw_transfer(rho_r%pw, rho_g%pw, error=error)
                      CALL pw_poisson_solve(poisson_env,rho_g%pw, pair_energy, pot_g%pw,error=error)

                      ! sum up to the full energy
                      scaling=fraction
                      IF (SIZE(mo_array)==1) scaling=scaling*2.0_dp
                      IF (iorb/=jorb) scaling=scaling*2.0_dp

                      exchange_energy=exchange_energy - scaling * pair_energy

                   ENDDO
                ENDDO

             ENDDO
          ENDDO
       ENDDO

       DO iorb_block=1,blocksize
          CALL pw_release(rho_i(iorb_block)%pw,error=error)
          CALL pw_release(rho_j(iorb_block)%pw,error=error)
       ENDDO

       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_r%pw,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_g%pw,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,pot_g%pw,error=error)

       iw = cp_print_key_unit_nr(logger,hfx_section,"HF_INFO",&
              extension=".scfLog",error=error)

       IF (iw>0) THEN
           WRITE (UNIT=iw,FMT="((T3,A,T61,F20.10))")&
                    "HF_PW_HFX| PW exchange energy:",exchange_energy
           WRITE (UNIT=iw,FMT="((T3,A,T61,F20.10),/)")&
                    "HF_PW_HFX| Gaussian exchange energy:",energy%ex
       ENDIF

       CALL cp_print_key_finished_output(iw,logger,hfx_section,&
                                         "HF_INFO", error=error)

    ENDIF

    CALL timestop(handle)

  END SUBROUTINE pw_hfx


! *****************************************************************************
!> \brief This routine modifies the xc section depending on the potential type
!>        used for the HF exchange and the resulting correction term. Currently
!>        three types of corrections are implemented:
!>
!>        coulomb:     Ex,hf = Ex,hf' + (PBEx-PBEx')
!>        shortrange:  Ex,hf = Ex,hf' + (XWPBEX-XWPBEX')
!>        truncated:   Ex,hf = Ex,hf' + ( (XWPBEX0-PBE_HOLE_TC_LR) -(XWPBEX0-PBE_HOLE_TC_LR)' )
!>
!>        with ' denoting the auxiliary basis set and
!>
!>          PBEx:           PBE exchange functional
!>          XWPBEX:         PBE exchange hole for short-range potential (erfc(omega*r)/r)
!>          XWPBEX0:        PBE exchange hole for standard coulomb potential
!>          PBE_HOLE_TC_LR: PBE exchange hole for longrange truncated coulomb potential
!>
!>
!> \param qs_env the qs environment
!> \param xc_section the original xc_section
!> \param admm_env the ADMM environment
!> \param error ...
!>
!> \param udate_only update gsi scaling  UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      12.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE create_admm_xc_section(qs_env, xc_section, admm_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_admm_xc_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=20)                        :: name_x_func
    INTEGER                                  :: hfx_potential_type, ifun, nfun
    LOGICAL                                  :: failure, funct_found
    REAL(dp)                                 :: cutoff_radius, hfx_fraction, &
                                                omega, scale_x
    TYPE(section_vals_type), POINTER         :: xc_fun, xc_fun_section

    NULLIFY(admm_env%xc_section_aux, admm_env%xc_section_primary)

    failure=.FALSE.


    CALL get_qs_env(qs_env, error=error)

    !! ** Duplicate existing xc-section
        CALL section_vals_duplicate(xc_section,admm_env%xc_section_aux,error=error)
        CALL section_vals_duplicate(xc_section,admm_env%xc_section_primary,error=error)
        !** Now modify the auxiliary basis
        !** First remove all functionals
        xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux,"XC_FUNCTIONAL",error=error)

        !* Overwrite possible shortcut
        CALL section_vals_val_set(xc_fun_section,"_SECTION_PARAMETERS_",&
                                  i_val=xc_funct_no_shortcut,error=error)

    !** Get number of Functionals in the list
    ifun = 0
    nfun = 0
    DO
      ifun = ifun+1
      xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
      IF (.NOT.ASSOCIATED(xc_fun)) EXIT
      nfun = nfun + 1
    END DO

    ifun = 0
    DO ifun = 1,nfun
      xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=1,error=error)
      IF (.NOT.ASSOCIATED(xc_fun)) EXIT
      CALL section_vals_remove_values(xc_fun, error=error)
    END DO

    hfx_potential_type = qs_env%x_data(1,1)%potential_parameter%potential_type
    hfx_fraction = qs_env%x_data(1,1)%general_parameter%fraction

    !in case of no admm exchange corr., no auxiliary exchange functional needed
    IF(admm_env%aux_exch_func == do_admm_aux_exch_func_none) THEN
      hfx_fraction = 0.0_dp
    END IF

    ! default PBE Functional
    IF(admm_env%aux_exch_func== do_admm_aux_exch_func_default .OR. &
       admm_env%aux_exch_func== do_admm_aux_exch_func_none ) THEN 
  
      !! ** Add functionals evaluated with auxiliary basis
      SELECT CASE (hfx_potential_type)
        CASE (do_hfx_potential_coulomb)
          CALL section_vals_val_set(xc_fun_section,"PBE%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                    r_val=-hfx_fraction,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_C",&
                                    r_val=0.0_dp,error=error)
        CASE (do_hfx_potential_short)
          omega =  qs_env%x_data(1,1)%potential_parameter%omega
          CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                    r_val=-hfx_fraction,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                    r_val=0.0_dp,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%OMEGA",&
                                    r_val=omega,error=error)
        CASE (do_hfx_potential_truncated)
          cutoff_radius = qs_env%x_data(1,1)%potential_parameter%cutoff_radius
          CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                    r_val=hfx_fraction,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                    r_val=cutoff_radius,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                    r_val=0.0_dp,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                    r_val=-hfx_fraction,error=error)
        CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

      !** Now modify the functionals for the primary basis
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary,"XC_FUNCTIONAL",error=error)
      !* Overwrite possible shortcut
      CALL section_vals_val_set(xc_fun_section,"_SECTION_PARAMETERS_",&
                              i_val=xc_funct_no_shortcut,error=error)

      SELECT CASE (hfx_potential_type)
      CASE (do_hfx_potential_coulomb)
        ifun = 0
        funct_found = .FALSE.
        DO
          ifun = ifun+1
          xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
          IF (.NOT.ASSOCIATED(xc_fun)) EXIT 
          IF( xc_fun%section%name == "PBE" ) THEN
            funct_found = .TRUE.
          END IF
        END DO
        IF( .NOT. funct_found ) THEN
          CALL section_vals_val_set(xc_fun_section,"PBE%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                    r_val=hfx_fraction,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_C",&
                                    r_val=0.0_dp,error=error)
        ELSE
          CALL section_vals_val_get(xc_fun_section,"PBE%SCALE_X",&
                                    r_val=scale_x,error=error)
          scale_x = scale_x + hfx_fraction
          CALL section_vals_val_set(xc_fun_section,"PBE%SCALE_X",&
                                    r_val=scale_x,error=error)
        END IF
      CASE (do_hfx_potential_short)
        omega =  qs_env%x_data(1,1)%potential_parameter%omega
        ifun = 0
        funct_found = .FALSE.
        DO
          ifun = ifun+1
          xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
          IF (.NOT.ASSOCIATED(xc_fun)) EXIT
          IF( xc_fun%section%name == "XWPBE" ) THEN
            funct_found = .TRUE.
          END IF
        END DO
        IF( .NOT. funct_found ) THEN
          CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                    r_val=hfx_fraction,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                    r_val=0.0_dp,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%OMEGA",&
                                    r_val=omega,error=error)
        ELSE
          CALL section_vals_val_get(xc_fun_section,"XWPBE%SCALE_X",&
                                    r_val=scale_x,error=error)
          scale_x = scale_x + hfx_fraction
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                    r_val=scale_x,error=error)
        END IF
      CASE (do_hfx_potential_truncated)
        cutoff_radius =  qs_env%x_data(1,1)%potential_parameter%cutoff_radius
        ifun = 0
        funct_found = .FALSE.
        DO
          ifun = ifun+1
          xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
          IF (.NOT.ASSOCIATED(xc_fun)) EXIT
          IF( xc_fun%section%name == "PBE_HOLE_T_C_LR" ) THEN
            funct_found = .TRUE.
          END IF
        END DO
        IF( .NOT. funct_found ) THEN
          CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                    r_val=-hfx_fraction,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                    r_val=cutoff_radius,error=error)
  
        ELSE
          CALL section_vals_val_get(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                    r_val=scale_x,error=error)
          scale_x = scale_x - hfx_fraction
          CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%SCALE_X",&
                                    r_val=scale_x,error=error)
          CALL section_vals_val_set(xc_fun_section,"PBE_HOLE_T_C_LR%CUTOFF_RADIUS",&
                                    r_val=cutoff_radius,error=error)
        END IF
        ifun = 0
        funct_found = .FALSE.
        DO
          ifun = ifun+1
          xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
          IF (.NOT.ASSOCIATED(xc_fun)) EXIT
          IF( xc_fun%section%name == "XWPBE" ) THEN
            funct_found = .TRUE.
          END IF
        END DO
        IF( .NOT. funct_found ) THEN
          CALL section_vals_val_set(xc_fun_section,"XWPBE%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                    r_val=hfx_fraction,error=error)
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X",&
                                    r_val=0.0_dp,error=error)

        ELSE
          CALL section_vals_val_get(xc_fun_section,"XWPBE%SCALE_X0",&
                                    r_val=scale_x,error=error)
          scale_x = scale_x + hfx_fraction
          CALL section_vals_val_set(xc_fun_section,"XWPBE%SCALE_X0",&
                                    r_val=scale_x,error=error)
        END IF

       END SELECT

    ! PBEX (always bare form), OPTX and Becke88 functional
    ELSE IF(admm_env%aux_exch_func== do_admm_aux_exch_func_pbex .OR. &
            admm_env%aux_exch_func== do_admm_aux_exch_func_opt .OR. &
            admm_env%aux_exch_func== do_admm_aux_exch_func_bee) THEN 
        IF(admm_env%aux_exch_func== do_admm_aux_exch_func_pbex) THEN
          name_x_func = 'PBE'
        ELSE IF(admm_env%aux_exch_func== do_admm_aux_exch_func_opt) THEN
          name_x_func = 'OPTX'      
        ELSE IF(admm_env%aux_exch_func== do_admm_aux_exch_func_bee) THEN
          name_x_func = 'BECKE88'
        END IF
        !primary basis
          CALL section_vals_val_set(xc_fun_section, TRIM(name_x_func) // "%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,TRIM(name_x_func) //"%SCALE_X",&
                                    r_val=-hfx_fraction,error=error)

          IF(admm_env%aux_exch_func== do_admm_aux_exch_func_pbex) THEN
            CALL section_vals_val_set(xc_fun_section,TRIM(name_x_func) //"%SCALE_C",&
                                    r_val=0.0_dp,error=error)
          END IF
          
      !** Now modify the functionals for the primary basis
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary,"XC_FUNCTIONAL",error=error)
      !* Overwrite possible L",error=error)
      !* Overwrite possible shortcut
      CALL section_vals_val_set(xc_fun_section,"_SECTION_PARAMETERS_",&
                              i_val=xc_funct_no_shortcut,error=error)

        ifun = 0
        funct_found = .FALSE.
        DO
          ifun = ifun+1
          xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
          IF (.NOT.ASSOCIATED(xc_fun)) EXIT
          IF( xc_fun%section%name == TRIM(name_x_func) ) THEN
            funct_found = .TRUE.
          END IF
        END DO
        IF( .NOT. funct_found ) THEN
          CALL section_vals_val_set(xc_fun_section,TRIM(name_x_func)//"%_SECTION_PARAMETERS_",&
                                    l_val=.TRUE.,error=error)
          CALL section_vals_val_set(xc_fun_section,TRIM(name_x_func)//"%SCALE_X",&
                                    r_val=hfx_fraction,error=error)
          IF(admm_env%aux_exch_func== do_admm_aux_exch_func_pbex) THEN
            CALL section_vals_val_set(xc_fun_section,TRIM(name_x_func) //"%SCALE_C",&
                                    r_val=0.0_dp,error=error)
          END IF


        ELSE
          CALL section_vals_val_get(xc_fun_section,TRIM(name_x_func)//"%SCALE_X",&
                                    r_val=scale_x,error=error)
          scale_x = scale_x + hfx_fraction
          CALL section_vals_val_set(xc_fun_section,TRIM(name_x_func)//"%SCALE_X",&
                                    r_val=scale_x,error=error)


        END IF

    END IF 




    IF( 1==0 ) THEN
      WRITE(*,*) "primary"
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_primary,"XC_FUNCTIONAL",error=error)
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT

        scale_x=-1000.0_dp
        IF(xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X",&
                                      r_val=scale_x,error=error)
        END IF
        IF( xc_fun%section%name == "XWPBE" ) THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X0",&
                                    r_val=hfx_fraction,error=error)

           WRITE(*,*) xc_fun%section%name, scale_x, hfx_fraction
        ELSE
          WRITE(*,*) xc_fun%section%name, scale_x
        END IF
      END DO

      WRITE(*,*) "auxiliary"
      xc_fun_section => section_vals_get_subs_vals(admm_env%xc_section_aux,"XC_FUNCTIONAL",error=error)
      ifun = 0
      funct_found = .FALSE.
      DO
        ifun = ifun+1
        xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,error=error)
        IF (.NOT.ASSOCIATED(xc_fun)) EXIT
        scale_x=-1000.0_dp
        IF(xc_fun%section%name /= "LYP" .AND. xc_fun%section%name /= "VWN") THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X",&
                                      r_val=scale_x,error=error)
        END IF
        IF( xc_fun%section%name == "XWPBE" ) THEN
          CALL section_vals_val_get(xc_fun,"SCALE_X0",&
                                    r_val=hfx_fraction,error=error)

           WRITE(*,*) xc_fun%section%name, scale_x, hfx_fraction
        ELSE
          WRITE(*,*) xc_fun%section%name, scale_x
        END IF
      END DO
    END IF

  END SUBROUTINE create_admm_xc_section

END MODULE hfx_admm_utils
