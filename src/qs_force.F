!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Quickstep force driver routine
!> \author MK (12.06.2002)
! *****************************************************************************
MODULE qs_force

  USE admm_methods,                    ONLY: calc_aux_mo_derivs_none,&
                                             calc_mixed_overlap_force
  USE admm_types,                      ONLY: admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_add,&
                                             cp_dbcsr_copy,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_deallocate_matrix_set
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_types,                     ONLY: cp_fm_type,&
                                             cp_fm_write
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dft_plus_u,                      ONLY: plus_u
  USE efield_utils,                    ONLY: calculate_ecore_efield
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE harris_env_types,                ONLY: harris_env_set,&
                                             harris_env_type
  USE harris_force,                    ONLY: harris_force_correction
  USE harris_force_types,              ONLY: harris_force_create,&
                                             harris_force_type
  USE input_constants,                 ONLY: do_admm_purify_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE message_passing,                 ONLY: mp_sum
  USE mulliken,                        ONLY: mulliken_restraint
  USE particle_types,                  ONLY: particle_type
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_dftb_dispersion,              ONLY: calculate_dftb_dispersion
  USE qs_dftb_matrices,                ONLY: build_dftb_matrices
  USE qs_energy,                       ONLY: qs_energies
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_external_potential,           ONLY: external_c_potential,&
                                             external_e_potential
  USE qs_force_types,                  ONLY: allocate_qs_force,&
                                             qs_force_type,&
                                             zero_qs_force
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE rt_propagation_utils,            ONLY: calc_c_mat_force
  USE scf_control_types,               ONLY: scf_control_type
  USE scptb_core_interactions,         ONLY: scptb_core_interaction
  USE scptb_core_matrix,               ONLY: build_scptb_core_matrix
  USE se_core_core,                    ONLY: se_core_core_interaction
  USE se_core_matrix,                  ONLY: build_se_core_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_force'

! *** Public subroutines ***

  PUBLIC :: qs_forces, write_forces

CONTAINS

! *****************************************************************************
!> \brief   Calculate the Quickstep forces.
!> \author  MK
!> \date    29.10.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE qs_forces(qs_env,globenv,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dir, handle, i, iatom, ikind, &
                                                ispin, istat, iw, natom, &
                                                nkind, nspin, output_unit
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: failure, gapw, harris_flag, &
                                                has_unit_metric, use_virial
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s, &
                                                matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb, &
                                                matrix_w, mo_derivs
    TYPE(cp_fm_type), POINTER                :: dp_imag, mo_coeff, &
                                                mo_coeff_aux_fit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY (logger)
    logger => cp_error_get_logger(error)
    para_env=>qs_env%para_env

    NULLIFY (atomic_kind_set)
    NULLIFY (dft_control)
    NULLIFY (force)
    NULLIFY (harris_env)
    NULLIFY (ks_env)
    NULLIFY (mos)
    NULLIFY (matrix_ks)
    NULLIFY (matrix_s)
    NULLIFY (matrix_w)
    NULLIFY (particle_set)
    NULLIFY (rho)
    NULLIFY (scf_control)
    NULLIFY (admm_env)
    NULLIFY (matrix_s_aux_fit_vs_orb)

    harris_flag = qs_env%use_harris
    CALL get_qs_env(qs_env=qs_env,&
         atomic_kind_set=atomic_kind_set,&
         dft_control=dft_control,&
         force=force,harris_env=harris_env,&
         particle_set=particle_set,&
         scf_control=scf_control,&
         virial=virial,&
         has_unit_metric=has_unit_metric,&
         energy=energy,error=error)

    natom = SIZE(particle_set)

    ! zero out the forces
    DO iatom=1,natom
       particle_set(iatom)%f=0.0_dp
    END DO

    gapw = (dft_control%qs_control%method=="GAPW")

    ALLOCATE (atom_of_kind(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom_of_kind",natom*int_size)

    ALLOCATE (kind_of(natom),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "kind_of",natom*int_size)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         atom_of_kind=atom_of_kind,&
         kind_of=kind_of)

    IF (.NOT.ASSOCIATED(force)) THEN
       !   *** Allocate the force data structure ***
       nkind = SIZE(atomic_kind_set)
       ALLOCATE (natom_of_kind(nkind),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            natom_of_kind=natom_of_kind)
       CALL allocate_qs_force(force,natom_of_kind)
       DEALLOCATE (natom_of_kind,STAT=istat)
       CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
       CALL set_qs_env(qs_env,force=force,error=error)
       IF (harris_flag) THEN
          NULLIFY(harris_force)
          CALL harris_force_create(harris_force=harris_force,natom=natom,error=error)
          CALL harris_env_set(harris_env=harris_env, harris_force=harris_force, error=error)
       END IF
    END IF
    CALL zero_qs_force(force)

    CALL get_qs_env(qs_env=qs_env,&
         matrix_s=matrix_s,&
         has_unit_metric=has_unit_metric,&
         error=error)

    CALL qs_energies(qs_env,calc_forces=.TRUE.,error=error)

    CALL get_qs_env(qs_env=qs_env,&
         ks_env=ks_env,&
         matrix_ks=matrix_ks,&
         matrix_s=matrix_s,&
         matrix_w=matrix_w,&
         mo_derivs=mo_derivs,&
         rho=rho,error=error)

    IF (qs_env%run_rtp) THEN
       CALL cp_dbcsr_allocate_matrix_set(matrix_w,dft_control%nspins,error=error)
       DO ispin=1,dft_control%nspins
          ALLOCATE(matrix_w(ispin)%matrix)
          CALL cp_dbcsr_init(matrix_w(ispin)%matrix,error=error)
          CALL cp_dbcsr_copy(matrix_w(ispin)%matrix,matrix_s(1)%matrix,&
               name="W MATRIX",error=error)
          CALL cp_dbcsr_set(matrix_w(ispin)%matrix,0.0_dp,error=error)
       END DO
       CALL set_qs_env(qs_env=qs_env,matrix_w=matrix_w,error=error)

       CALL calc_c_mat_force(qs_env,dp_imag,error)
    END IF

    IF (.NOT.has_unit_metric) THEN
       nspin = SIZE(matrix_w)
       DO ispin=1,nspin
          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX",error=error),cp_p_file)) THEN
             iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX",&
                  extension=".Log",error=error)
             CALL cp_dbcsr_write_sparse_matrix(matrix_w(ispin)%matrix,4,6,qs_env,para_env,output_unit=iw,error=error)
             CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
                  "DFT%PRINT%AO_MATRICES/W_MATRIX", error=error)
          END IF
       END DO
    ENDIF

    ! from an eventual Mulliken restraint
    IF (dft_control%qs_control%mulliken_restraint) THEN
       CALL mulliken_restraint(dft_control%qs_control%mulliken_restraint_control, &
            para_env,matrix_s(1)%matrix, rho%rho_ao,w_matrix=matrix_w,error=error)
    END IF

    ! Add non-Pulay contribution of DFT+U to W matrix, since it has also to be
    ! digested with overlap matrix derivatives
    IF (dft_control%dft_plus_u) THEN
       CALL plus_u(qs_env=qs_env,matrix_w=matrix_w,error=error)
    END IF

    IF (.NOT. harris_flag) THEN ! harris computes the forces differenctly
       ! Compute core forces (also overwrites matrix_w)
       IF (dft_control%qs_control%semi_empirical) THEN
          CALL build_se_core_matrix(qs_env=qs_env,para_env=para_env,&
               calculate_forces=.TRUE.,error=error)
          CALL se_core_core_interaction(qs_env, para_env, calculate_forces=.TRUE.,&
               error=error)
       ELSEIF (dft_control%qs_control%dftb) THEN
          CALL build_dftb_matrices(qs_env=qs_env,para_env=para_env,&
               calculate_forces=.TRUE.,error=error)
          CALL calculate_dftb_dispersion(qs_env=qs_env,para_env=para_env,&
               calculate_forces=.TRUE.,error=error)
       ELSEIF ( dft_control%qs_control%scptb ) THEN
          CALL build_scptb_core_matrix(qs_env=qs_env,calculate_forces=.TRUE.,error=error)
          CALL scptb_core_interaction(qs_env,calculate_forces=.TRUE.,error=error)
       ELSE
          CALL build_core_hamiltonian_matrix(qs_env=qs_env,calculate_forces=.TRUE.,error=error)
          CALL calculate_ecore_self(qs_env,error=error)
          CALL calculate_ecore_overlap(qs_env, para_env, &
               calculate_forces=.TRUE.,&
               error=error)
          CALL calculate_ecore_efield(qs_env,calculate_forces=.TRUE.,error=error)
          CALL external_c_potential(qs_env,calculate_forces=.TRUE.,error=error)
          CALL external_e_potential(qs_env,error)
       END IF

          ! Compute grid-based forces
          CALL qs_ks_update_qs_env(ks_env=ks_env,&
                                   qs_env=qs_env,&
                                   calculate_forces=.TRUE.,error=error)
   
       IF(ASSOCIATED(qs_env%mp2_env)) THEN
         ! with MP2 we have to recalculate the SCF energy with the
         ! correct density
         DO ispin=1, qs_env%dft_control%nspins
           CALL cp_dbcsr_add(qs_env%rho%rho_ao(ispin)%matrix, qs_env%matrix_p_mp2(ispin)%matrix, 1.0_dp,  -1.0_dp, error)
         END DO
         CALL qs_rho_update_rho(qs_env%rho, qs_env=qs_env, error=error)
         CALL qs_ks_did_change(ks_env,rho_changed=.TRUE.,&
                               s_mstruct_changed=qs_env%mp2_env%ri_mp2%free_hfx_buffer,&
                               error=error)
         CALL qs_ks_update_qs_env(ks_env=ks_env,&
                                  qs_env=qs_env,&
                                  just_energy=.TRUE.,error=error)
         energy%total = energy%total + energy%mp2
         ! deallocate mp2_W 
         CALL cp_dbcsr_deallocate_matrix_set(qs_env%matrix_w_mp2,error=error)
       END IF

       ! Add forces resulting from wavefunction fitting

       IF (dft_control%do_admm ) THEN
          CALL get_qs_env(qs_env=qs_env,&
               matrix_s_aux_fit=matrix_s_aux_fit,&
               matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
               mos=mos,&
               admm_env=admm_env,&
               error=error)
          DO ispin=1,nspin             
             mo_set => mos(ispin)%mo_set
             CALL get_mo_set(mo_set=mo_set,mo_coeff=mo_coeff)

             ! if no purification we need to calculate the H matrix for forces
             IF (admm_env%purification_method==do_admm_purify_none) THEN
                CALL get_qs_env(qs_env, mos_aux_fit=mos_aux_fit,&
                     error=error)
                CALL get_mo_set(mo_set=mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
                CALL calc_aux_mo_derivs_none(ispin, qs_env%admm_env, mo_set, &
                     mo_coeff, mo_coeff_aux_fit, &
                     qs_env%matrix_ks_aux_fit, error=error)
                IF(.FALSE.) CALL cp_fm_write(admm_env%H(ispin)%matrix, handle, .TRUE., .FALSE., error)            
             END IF

             use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

             CALL calc_mixed_overlap_force(qs_env, para_env, ispin, admm_env, mo_coeff, &
                                           matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, logger, &
                                           iw, use_virial, error)
          END DO
       END IF

       !  *** replicate forces ***
       DO ikind=1,SIZE(force)
          CALL mp_sum(force(ikind)%overlap,para_env%group)
          CALL mp_sum(force(ikind)%overlap_admm,para_env%group)
          CALL mp_sum(force(ikind)%kinetic,para_env%group)
          CALL mp_sum(force(ikind)%gth_ppl,para_env%group)
          CALL mp_sum(force(ikind)%gth_nlcc,para_env%group)
          CALL mp_sum(force(ikind)%gth_ppnl,para_env%group)
          CALL mp_sum(force(ikind)%all_potential,para_env%group)
          CALL mp_sum(force(ikind)%core_overlap,para_env%group)
          CALL mp_sum(force(ikind)%rho_core,para_env%group)
          CALL mp_sum(force(ikind)%rho_elec,para_env%group)
          CALL mp_sum(force(ikind)%vhxc_atom,para_env%group)
          CALL mp_sum(force(ikind)%g0s_Vh_elec,para_env%group)
          CALL mp_sum(force(ikind)%fock_4c,para_env%group)
          CALL mp_sum(force(ikind)%mp2_non_sep,para_env%group)
          CALL mp_sum(force(ikind)%mp2_sep,para_env%group)
          CALL mp_sum(force(ikind)%hfx_ri,para_env%group)
          CALL mp_sum(force(ikind)%repulsive,para_env%group)
          CALL mp_sum(force(ikind)%dispersion,para_env%group)
          CALL mp_sum(force(ikind)%ehrenfest,para_env%group)

          force(ikind)%total(:,:) = force(ikind)%total(:,:) +&
               force(ikind)%core_overlap(:,:) +&
               force(ikind)%gth_ppl(:,:) +&
               force(ikind)%gth_nlcc(:,:) +&
               force(ikind)%gth_ppnl(:,:) +&
               force(ikind)%all_potential(:,:) +&
               force(ikind)%kinetic(:,:) +&
               force(ikind)%overlap(:,:) +&
               force(ikind)%overlap_admm(:,:) +&
               force(ikind)%rho_core(:,:) +&
               force(ikind)%rho_elec(:,:) +&
               force(ikind)%vhxc_atom(:,:) +&
               force(ikind)%g0s_Vh_elec(:,:) +&
               force(ikind)%fock_4c(:,:) +&
               force(ikind)%mp2_non_sep(:,:) +&
               force(ikind)%mp2_sep(:,:) +&
               force(ikind)%hfx_ri(:,:) +&
               force(ikind)%repulsive(:,:)  +&
               force(ikind)%dispersion(:,:) +&
               force(ikind)%ehrenfest(:,:) +&
               force(ikind)%efield(:,:)+&
               force(ikind)%eev(:,:)

       END DO


       DO iatom=1,natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          ! the force is - dE/dR, what is called force is actually the gradient
          ! Things should have the right name
          ! The minus sign below is a hack
          ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          force(ikind)%other(1:3,i)=-particle_set(iatom)%f(1:3)+force(ikind)%ch_pulay(1:3,i)
          force(ikind)%total(1:3,i)=force(ikind)%total(1:3,i)+force(ikind)%other(1:3,i)
          particle_set(iatom)%f = -force(ikind)%total(1:3,i)
       END DO

       !   *** distribute virial ***
       IF (virial%pv_availability) THEN
          CALL mp_sum(virial%pv_virial,para_env%group)
          !  *** add the volume terms of the virial ***
          IF ((.NOT.virial%pv_numer) .AND. &
              (.NOT.(dft_control%qs_control%dftb .OR. &
                     qs_env%dft_control%qs_control%scptb .OR. &
                     dft_control%qs_control%semi_empirical))) THEN
             DO dir=1, 3
                virial%pv_virial(dir,dir) = virial%pv_virial(dir,dir) - energy%exc  &
                     - 2.0_dp*energy%hartree
                IF( dft_control%do_admm ) THEN
                  virial%pv_virial(dir,dir) = virial%pv_virial(dir,dir) - energy%exc_aux_fit
                END IF
                ! The factor 2 is a hack. It compensates the plus sign in h_stress/pw_poisson_solve.
                ! The sign in pw_poisson_solve is correct for FIST, but not for QS.
                ! There should be a more elegant solution to that...
             END DO
          END IF
       END IF

       output_unit = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%DERIVATIVES",&
            extension=".Log",error=error)
       print_section => section_vals_get_subs_vals(qs_env%input,"DFT%PRINT%DERIVATIVES",error=error)
       IF (dft_control%qs_control%semi_empirical.AND.has_unit_metric) THEN
          CALL write_forces(force,atomic_kind_set,2,output_unit=output_unit,&
               print_section=print_section,error=error)
       ELSE IF (dft_control%qs_control%dftb.AND.(.NOT.has_unit_metric)) THEN
          CALL write_forces(force,atomic_kind_set,4,output_unit=output_unit,&
               print_section=print_section,error=error)
       ELSE IF (gapw) THEN
          CALL write_forces(force,atomic_kind_set,1,output_unit=output_unit,&
               print_section=print_section,error=error)
       ELSE IF (dft_control%qs_control%scptb) THEN
          CALL write_forces(force,atomic_kind_set,5,output_unit=output_unit,&
               print_section=print_section,error=error)
       ELSE
          CALL write_forces(force,atomic_kind_set,0,output_unit=output_unit,&
               print_section=print_section,error=error)
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,qs_env%input,&
            "DFT%PRINT%DERIVATIVES",error=error)
    ELSE
       CALL get_qs_env(qs_env=qs_env, harris_env=harris_env,error=error)
       CALL harris_force_correction(qs_env=qs_env, harris_env=harris_env, &
            globenv=globenv,error=error)
       DO iatom = 1,natom
          i = atom_of_kind(iatom)
          particle_set(iatom)%f = -harris_env%harris_force%f_harris(iatom, 1:3)
       END DO
    END IF

    CALL cp_dbcsr_deallocate_matrix_set(matrix_w,error=error)
    CALL set_qs_env(qs_env=qs_env,matrix_w=matrix_w,error=error)

    DEALLOCATE (atom_of_kind,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"atom_of_kind")

    DEALLOCATE (kind_of,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"kind_of")

    CALL timestop(handle)

  END SUBROUTINE qs_forces

! *****************************************************************************
!> \brief   Write a Quickstep force data structure to output unit
!> \author  MK
!> \date    05.06.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_forces(qs_force,atomic_kind_set,ftype,output_unit,&
                          print_section, error)

    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: qs_force
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    INTEGER, INTENT(IN)                      :: ftype, output_unit
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_forces', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=13)                        :: fmtstr5
    CHARACTER(LEN=15)                        :: fmtstr4
    CHARACTER(LEN=20)                        :: fmtstr3
    CHARACTER(LEN=35)                        :: fmtstr2
    CHARACTER(LEN=48)                        :: fmtstr1
    INTEGER                                  :: i, iatom, ikind, istat, &
                                                my_ftype, natom, ndigits
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of
    REAL(KIND=dp), DIMENSION(3)              :: grand_total

    IF (output_unit>0) THEN

       IF (.NOT.ASSOCIATED(qs_force)) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "The qs_force pointer is not associated "//&
                            "and cannot be printed")
       END IF

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            natom=natom)
       ALLOCATE (atom_of_kind(natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "atom_of_kind",natom*int_size)
       ALLOCATE (kind_of(natom),STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                        "kind_of",natom*int_size)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            atom_of_kind=atom_of_kind,&
            kind_of=kind_of)

       ! Variable precision output of the forces
       CALL section_vals_val_get(print_section,"NDIGITS",&
            i_val=ndigits,error=error)

       fmtstr1 = "(/,/,T2,A,/,/,T3,A,T11,A,T23,A,T40,A1,2(  X,A1))"
       WRITE (UNIT=fmtstr1(41:42),FMT="(I2)") ndigits + 5

       fmtstr2 = "(/,(T2,I5,4X,I4,T18,A,T34,3F  .  ))"
       WRITE (UNIT=fmtstr2(32:33),FMT="(I2)") ndigits
       WRITE (UNIT=fmtstr2(29:30),FMT="(I2)") ndigits + 6

       fmtstr3 = "(/,T3,A,T34,3F  .  )"
       WRITE (UNIT=fmtstr3(18:19),FMT="(I2)") ndigits
       WRITE (UNIT=fmtstr3(15:16),FMT="(I2)") ndigits + 6

       fmtstr4 = "((T34,3F  .  ))"
       WRITE (UNIT=fmtstr4(12:13),FMT="(I2)") ndigits
       WRITE (UNIT=fmtstr4(9:10),FMT="(I2)") ndigits + 6

       fmtstr5 = "(/T2,A//T3,A)"

       WRITE (UNIT=output_unit,FMT=fmtstr1)&
         "FORCES [a.u.]","Atom","Kind","Component","X","Y","Z"

       grand_total(:) = 0.0_dp

       my_ftype = ftype

       SELECT CASE (my_ftype)
       CASE DEFAULT
          DO iatom=1,natom
             ikind = kind_of(iatom)
             i = atom_of_kind(iatom)
             WRITE (UNIT=output_unit,FMT=fmtstr2)&
                  iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
             grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
          END DO
       CASE (0)
          DO iatom=1,natom
             ikind = kind_of(iatom)
             i = atom_of_kind(iatom)
             WRITE (UNIT=output_unit,FMT=fmtstr2)&
                  iatom,ikind,"       overlap",qs_force(ikind)%overlap(1:3,i),&
                  iatom,ikind,"  overlap_admm",qs_force(ikind)%overlap_admm(1:3,i),&
                  iatom,ikind,"       kinetic",qs_force(ikind)%kinetic(1:3,i),&
                  iatom,ikind,"       gth_ppl",qs_force(ikind)%gth_ppl(1:3,i),&
                  iatom,ikind,"      gth_nlcc",qs_force(ikind)%gth_nlcc(1:3,i),&
                  iatom,ikind,"      gth_ppnl",qs_force(ikind)%gth_ppnl(1:3,i),&
                  iatom,ikind,"  core_overlap",qs_force(ikind)%core_overlap(1:3,i),&
                  iatom,ikind,"      rho_core",qs_force(ikind)%rho_core(1:3,i),&
                  iatom,ikind,"      rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
                  iatom,ikind,"      ch_pulay",qs_force(ikind)%ch_pulay(1:3,i),&
                  iatom,ikind,"    dispersion",qs_force(ikind)%dispersion(1:3,i),&
                  iatom,ikind,"         other",qs_force(ikind)%other(1:3,i),&
                  iatom,ikind,"       fock_4c",qs_force(ikind)%fock_4c(1:3,i),&
                  iatom,ikind,"        hfx_ri",qs_force(ikind)%hfx_ri(1:3,i),&
                  iatom,ikind,"     ehrenfest",qs_force(ikind)%ehrenfest(1:3,i),&
                  iatom,ikind,"        efield",qs_force(ikind)%efield(1:3,i),&
                  iatom,ikind,"           eev",qs_force(ikind)%eev(1:3,i),&
                  iatom,ikind,"   mp2_non_sep",qs_force(ikind)%mp2_non_sep(1:3,i),&
                  iatom,ikind,"       mp2_sep",qs_force(ikind)%mp2_sep(1:3,i),&
                  iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
             grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
          END DO
       CASE (1)
          DO iatom=1,natom
             ikind = kind_of(iatom)
             i = atom_of_kind(iatom)
             WRITE (UNIT=output_unit,FMT=fmtstr2)&
                  iatom,ikind,"       overlap",qs_force(ikind)%overlap(1:3,i),&
                  iatom,ikind,"  overlap_admm",qs_force(ikind)%overlap_admm(1:3,i),&
                  iatom,ikind,"       kinetic",qs_force(ikind)%kinetic(1:3,i),&
                  iatom,ikind,"       gth_ppl",qs_force(ikind)%gth_ppl(1:3,i),&
                  iatom,ikind,"      gth_nlcc",qs_force(ikind)%gth_nlcc(1:3,i),&
                  iatom,ikind,"      gth_ppnl",qs_force(ikind)%gth_ppnl(1:3,i),&
                  iatom,ikind," all_potential",qs_force(ikind)%all_potential(1:3,i),&
                  iatom,ikind,"  core_overlap",qs_force(ikind)%core_overlap(1:3,i),&
                  iatom,ikind,"      rho_core",qs_force(ikind)%rho_core(1:3,i),&
                  iatom,ikind,"      rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
                  iatom,ikind,"     vhxc_atom",qs_force(ikind)%vhxc_atom(1:3,i),&
                  iatom,ikind,"   g0s_Vh_elec",qs_force(ikind)%g0s_Vh_elec(1:3,i),&
                  iatom,ikind,"      ch_pulay",qs_force(ikind)%ch_pulay(1:3,i),&
                  iatom,ikind,"    dispersion",qs_force(ikind)%dispersion(1:3,i),&
                  iatom,ikind,"       fock_4c",qs_force(ikind)%fock_4c(1:3,i),&
                  iatom,ikind,"        hfx_ri",qs_force(ikind)%hfx_ri(1:3,i),&
                  iatom,ikind,"        efield",qs_force(ikind)%efield(1:3,i),&
                  iatom,ikind,"           eev",qs_force(ikind)%eev(1:3,i),&
                  iatom,ikind,"   mp2_non_sep",qs_force(ikind)%mp2_non_sep(1:3,i),&
                  iatom,ikind,"       mp2_sep",qs_force(ikind)%mp2_sep(1:3,i),&
                  iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
             grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
          END DO
       CASE (2)
          DO iatom=1,natom
             ikind = kind_of(iatom)
             i = atom_of_kind(iatom)
             WRITE (UNIT=output_unit,FMT=fmtstr2)&
                  iatom,ikind," all_potential",qs_force(ikind)%all_potential(1:3,i),&
                  iatom,ikind,"      rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
                  iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
             grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
          END DO
       CASE (3)
          DO iatom=1,natom
             ikind = kind_of(iatom)
             i = atom_of_kind(iatom)
             WRITE (UNIT=output_unit,FMT=fmtstr2)&
                  iatom,ikind,"        overlap",qs_force(ikind)%overlap(1:3,i),&
                  iatom,ikind,"overlap_admm",qs_force(ikind)%overlap_admm(1:3,i),&
                  iatom,ikind,"        kinetic",qs_force(ikind)%kinetic(1:3,i),&
                  iatom,ikind,"        gth_ppl",qs_force(ikind)%gth_ppl(1:3,i),&
                  iatom,ikind,"       gth_nlcc",qs_force(ikind)%gth_nlcc(1:3,i),&
                  iatom,ikind,"       gth_ppnl",qs_force(ikind)%gth_ppnl(1:3,i),&
                  iatom,ikind,"   core_overlap",qs_force(ikind)%core_overlap(1:3,i),&
                  iatom,ikind,"       rho_core",qs_force(ikind)%rho_core(1:3,i),&
                  iatom,ikind,"       rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
                  iatom,ikind,"       ch_pulay",qs_force(ikind)%ch_pulay(1:3,i),&
                  iatom,ikind,"        fock_4c",qs_force(ikind)%fock_4c(1:3,i),&
                  iatom,ikind,"        hfx_ri",qs_force(ikind)%hfx_ri(1:3,i),&
                  iatom,ikind,"   mp2_non_sep",qs_force(ikind)%mp2_non_sep(1:3,i),&
                  iatom,ikind,"       mp2_sep",qs_force(ikind)%mp2_sep(1:3,i),&
                  iatom,ikind,"          total",qs_force(ikind)%total(1:3,i)
             grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
          END DO
       CASE (4)
          DO iatom=1,natom
             ikind = kind_of(iatom)
             i = atom_of_kind(iatom)
             WRITE (UNIT=output_unit,FMT=fmtstr2)&
                  iatom,ikind,"  all_potential",qs_force(ikind)%all_potential(1:3,i),&
                  iatom,ikind,"        overlap",qs_force(ikind)%overlap(1:3,i),&
                  iatom,ikind,"   overlap_admm",qs_force(ikind)%overlap_admm(1:3,i),&
                  iatom,ikind,"       rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
                  iatom,ikind,"      repulsive",qs_force(ikind)%repulsive(1:3,i),&
                  iatom,ikind,"     dispersion",qs_force(ikind)%dispersion(1:3,i),&
                  iatom,ikind,"     ehrenfest",qs_force(ikind)%ehrenfest(1:3,i),&
                  iatom,ikind,"          total",qs_force(ikind)%total(1:3,i)
             grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
          END DO
       CASE (5)
          DO iatom=1,natom
             ikind = kind_of(iatom)
             i = atom_of_kind(iatom)
             WRITE (UNIT=output_unit,FMT=fmtstr2)&
                  iatom,ikind,"       overlap",qs_force(ikind)%overlap(1:3,i),&
                  iatom,ikind,"       kinetic",qs_force(ikind)%kinetic(1:3,i),&
                  iatom,ikind,"      rho_elec",qs_force(ikind)%rho_elec(1:3,i),&
                  iatom,ikind,"    dispersion",qs_force(ikind)%dispersion(1:3,i),&
                  iatom,ikind," all potential",qs_force(ikind)%all_potential(1:3,i),&
                  iatom,ikind,"         other",qs_force(ikind)%other(1:3,i),&
                  iatom,ikind,"         total",qs_force(ikind)%total(1:3,i)
             grand_total(1:3) = grand_total(1:3) + qs_force(ikind)%total(1:3,i)
          END DO
       END SELECT

       WRITE (UNIT=output_unit,FMT=fmtstr3) "Sum of total",grand_total(1:3)

       DEALLOCATE (atom_of_kind,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"atom_of_kind")
       DEALLOCATE (kind_of,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"kind_of")

    END IF

  END SUBROUTINE write_forces

END MODULE qs_force
