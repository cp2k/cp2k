!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Density Derived atomic point charges from a QM calculation
!>      (see Bloechl, J. Chem. Phys. Vol. 103 pp. 7422-7428)
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE cp_ddapc
  USE bibliography,                    ONLY: Blochl1995,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: ddapc_restraint_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_ddapc_forces,                 ONLY: ewald_ddapc_force,&
                                             reset_ch_pulay,&
                                             restraint_functional_force,&
                                             solvation_ddapc_force
  USE cp_ddapc_util,                   ONLY: get_ddapc,&
                                             modify_hartree_pot,&
                                             restraint_functional_potential
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants,                 ONLY: do_spin_density
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE pw_methods,                      ONLY: pw_integral_ab,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc'

  PUBLIC :: cp_ddapc_apply_CD,&  ! Apply Coupling/Decoupling to Periodic Images
            cp_ddapc_apply_RS,&  ! Apply Restraints/Constraints
            cp_ddapc_apply_RF,&  ! Apply an SCRF (Solvation Scheme)
            qs_ks_ddapc

CONTAINS

! *****************************************************************************
!> \brief Set of methods using DDAPC charges
!> \par History
!>      08.2005 created [tlaino]
!>      08.2008 extended to restraint/constraint DDAPC charges [fschiff]
! *****************************************************************************
  SUBROUTINE qs_ks_ddapc(qs_env,auxbas_pw_pool,rho_tot_gspace,v_hartree_gspace,&
       v_spin_ddapc_rest_r,energy,calculate_forces, ks_env, ks_matrix,just_energy,&
       error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_p_type)                          :: rho_tot_gspace, &
                                                v_hartree_gspace, &
                                                v_spin_ddapc_rest_r
    TYPE(qs_energy_type), POINTER            :: energy
    LOGICAL, INTENT(in)                      :: calculate_forces
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    LOGICAL, INTENT(in)                      :: just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_ks_ddapc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ddapc_size, handle, i, my_id, &
                                                stat
    LOGICAL                                  :: ddapc_restraint_is_spin, &
                                                et_coupling_calc, &
                                                explicit_potential, failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(pw_p_type)                          :: v_spin_ddapc_rest_g

    CALL timeset(routineN,handle)
    CALL cite_reference(Blochl1995)
    ! In case decouple periodic images and/or apply restraints to charges
    failure        = .FALSE.
    logger         => cp_error_get_logger(error)
    ddapc_restraint_is_spin =.FALSE.
    et_coupling_calc        =.FALSE.

    IF(qs_env%dft_control%qs_control%ddapc_restraint)THEN
       ddapc_size=SIZE(qs_env%dft_control%qs_control%ddapc_restraint_control)
       IF(SIZE(energy%ddapc_restraint).NE.ddapc_size)THEN
          DEALLOCATE(energy%ddapc_restraint,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(energy%ddapc_restraint(ddapc_size),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF

       DO i=1,SIZE( qs_env%dft_control%qs_control%ddapc_restraint_control )
          my_id= qs_env%dft_control%qs_control%ddapc_restraint_control(i)%ddapc_restraint_control%density_type
          IF (my_id==do_spin_density.OR.ddapc_restraint_is_spin) ddapc_restraint_is_spin=.TRUE.
       END DO
       et_coupling_calc= qs_env%dft_control%qs_control%et_coupling_calc
    ENDIF

    explicit_potential=ddapc_restraint_is_spin.OR.et_coupling_calc
    qs_env%dft_control%qs_control%ddapc_explicit_potential = explicit_potential
    qs_env%dft_control%qs_control%ddapc_restraint_is_spin = ddapc_restraint_is_spin
    IF (explicit_potential) THEN
       CALL pw_pool_create_pw(auxbas_pw_pool,v_spin_ddapc_rest_g%pw,&
            use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE,error=error)
       CALL pw_zero(v_spin_ddapc_rest_g%pw, error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,v_spin_ddapc_rest_r%pw,&
            use_data=REALDATA3D,in_space=REALSPACE,error=error)
    ENDIF

    IF (calculate_forces) CALL reset_ch_pulay(qs_env, error=error)

    ! Decoupling/Recoupling
    CALL cp_ddapc_apply_CD(qs_env, rho_tot_gspace, energy%hartree, v_hartree_gspace,&
         calculate_forces, Itype_of_density="FULL DENSITY", error=error)
    IF(qs_env%dft_control%qs_control%ddapc_restraint)THEN
       ! Restraints/Constraints
       DO i=1,ddapc_size
          NULLIFY(ddapc_restraint_control)
          ddapc_restraint_control=> qs_env%dft_control%qs_control%ddapc_restraint_control(i)%ddapc_restraint_control

          CALL cp_ddapc_apply_RS(qs_env, energy%hartree, energy%ddapc_restraint(i), v_hartree_gspace,&
               v_spin_ddapc_rest_g, i,ddapc_restraint_control,calculate_forces, error)
       END DO
    END IF
    CALL cp_ddapc_apply_RF(qs_env, rho_tot_gspace, energy%hartree, v_hartree_gspace,&
         calculate_forces, Itype_of_density="FULL DENSITY", error=error)

    ! CJM Copying the real-space Hartree potential to KS_ENV
    IF ((.NOT. just_energy).OR.et_coupling_calc) THEN
       CALL pw_transfer(v_hartree_gspace%pw, ks_env%v_hartree_rspace%pw, error=error)
       CALL pw_scale(ks_env%v_hartree_rspace%pw, ks_env%v_hartree_rspace%pw%pw_grid%dvol, error=error)
       IF (explicit_potential) THEN
          CALL pw_transfer(v_spin_ddapc_rest_g%pw,v_spin_ddapc_rest_r%pw, error=error)
          CALL pw_scale(v_spin_ddapc_rest_r%pw,v_spin_ddapc_rest_r%pw%pw_grid%dvol, error=error)
          IF(et_coupling_calc)THEN
             IF(qs_env%et_coupling%keep_matrix)THEN
                IF(qs_env%et_coupling%first_run)THEN
                   NULLIFY(qs_env%et_coupling%rest_mat(1)%matrix)
                   ALLOCATE(qs_env%et_coupling%rest_mat(1)%matrix)
                   CALL cp_dbcsr_init(qs_env%et_coupling%rest_mat(1)%matrix, error=error)
                   CALL cp_dbcsr_copy(qs_env%et_coupling%rest_mat(1)%matrix, ks_matrix(1)%matrix, &
                        name="ET_RESTRAINT_MATRIX_B", error=error)
                   CALL cp_dbcsr_set(qs_env%et_coupling%rest_mat(1)%matrix, 0.0_dp, error=error)
                   CALL integrate_v_rspace(v_spin_ddapc_rest_r,&
                        h=qs_env%et_coupling%rest_mat(1),&
                        qs_env=qs_env,calculate_forces=.FALSE.,error=error)
                   qs_env%et_coupling%order_p=&
                        qs_env%dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%ddapc_order_p
                   qs_env%et_coupling%e1=qs_env%dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%strength
                   qs_env%et_coupling%keep_matrix=.FALSE.
                ELSE
                   NULLIFY(qs_env%et_coupling%rest_mat(2)%matrix)
                   ALLOCATE(qs_env%et_coupling%rest_mat(2)%matrix)
                   CALL cp_dbcsr_init(qs_env%et_coupling%rest_mat(2)%matrix, error=error)
                   CALL cp_dbcsr_copy(qs_env%et_coupling%rest_mat(2)%matrix, ks_matrix(1)%matrix, &
                        name="ET_RESTRAINT_MATRIX_B", error=error)
                   CALL cp_dbcsr_set(qs_env%et_coupling%rest_mat(2)%matrix, 0.0_dp, error=error)
                   CALL integrate_v_rspace(v_spin_ddapc_rest_r,&
                        h=qs_env%et_coupling%rest_mat(2),&
                        qs_env=qs_env,calculate_forces=.FALSE.,error=error)
                END IF
             END IF
          END IF
       ENDIF
    ENDIF

    IF (explicit_potential) THEN
       CALL pw_pool_give_back_pw(auxbas_pw_pool,v_spin_ddapc_rest_g%pw,error=error)
    ENDIF
    CALL timestop(handle)

  END SUBROUTINE qs_ks_ddapc

! *****************************************************************************
!> \brief Routine to couple/decouple periodic images with the Bloechl scheme
!>
!>      The coupling/decoupling is obtaines evaluating terms E2 and E3 in
!>      J. Chem. Phys. Vol. 103 pp. 7422-7428.. The E2 terms is just a
!>      Ewald summation, and for performance reason I'm writing a specific
!>      driver instead of using and setting-up the environment of the already
!>      available routines
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE cp_ddapc_apply_CD(qs_env, rho_tot_gspace, energy, v_hartree_gspace, &
       calculate_forces, Itype_of_density, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type)                          :: rho_tot_gspace
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    TYPE(pw_p_type)                          :: v_hartree_gspace
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    CHARACTER(LEN=*)                         :: Itype_of_density
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_apply_CD', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, stat
    LOGICAL                                  :: apply_decpl, failure, need_f
    REAL(KINd=dp)                            :: e_decpl, e_recpl
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges, radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dq
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER :: density_fit_section, &
      force_env_section, multipole_section, poisson_section, &
      qmmm_periodic_section

    CALL timeset(routineN,handle)
    failure        = .FALSE.
    need_f         = .FALSE.
    IF (PRESENT(calculate_forces)) need_f   = calculate_forces
    logger         => cp_error_get_logger(error)
    apply_decpl    = qs_env%cp_ddapc_ewald%do_decoupling.OR.qs_env%cp_ddapc_ewald%do_qmmm_periodic
    IF ((.NOT.failure).AND.(apply_decpl)) THEN
       ! Initialize
       NULLIFY(multipole_section,&
               poisson_section,&
               force_env_section,&
               particle_set,&
               qmmm_periodic_section,&
               density_fit_section,&
               charges,&
               radii,&
               dq,&
               cell,&
               super_cell)

       CALL get_qs_env(qs_env=qs_env,&
                       input=force_env_section,&
                       particle_set=particle_set,&
                       cell=cell,&
                       super_cell=super_cell,&
                       error=error)
       CPPostcondition(ASSOCIATED(qs_env%cp_ddapc_ewald),cp_failure_level,routineP,error,failure)
       poisson_section => section_vals_get_subs_vals(force_env_section,"DFT%POISSON", error=error)

       density_fit_section =>section_vals_get_subs_vals(force_env_section,"DFT%DENSITY_FITTING", error=error)

       IF (qs_env%cp_ddapc_ewald%do_decoupling) THEN
          multipole_section => section_vals_get_subs_vals(poisson_section,"MULTIPOLE", error=error)
       END IF
       IF (qs_env%cp_ddapc_ewald%do_qmmm_periodic) THEN
          qmmm_periodic_section   => section_vals_get_subs_vals(force_env_section,"QMMM%PERIODIC", error=error)
          multipole_section       => section_vals_get_subs_vals(qmmm_periodic_section,"MULTIPOLE", error=error)
       END IF
       ! Start the real calculation
       iw=cp_print_key_unit_nr(logger,multipole_section,"PROGRAM_RUN_INFO",&
            extension=".fitChargeLog",error=error)
       ! First we evaluate the charges at the corresponding SCF STEP
       IF (need_f) THEN
          CALL get_ddapc(qs_env,&
                         need_f,&
                         density_fit_section,&
                         qout1=charges,&
                         out_radii=radii,&
                         dq_out=dq,&
                         ext_rho_tot_g=rho_tot_gspace%pw,&
                         Itype_of_density=Itype_of_density,&
                         error=error)
       ELSE
          CALL get_ddapc(qs_env,&
                         need_f,&
                         density_fit_section,&
                         qout1=charges,&
                         out_radii=radii,&
                         ext_rho_tot_g=rho_tot_gspace%pw,&
                         Itype_of_density=Itype_of_density,&
                         error=error)
       END IF
       ! Evaluate the Ewald contribution to the decoupling/coupling E2 and E3
       IF (iw>0) THEN
          e_decpl = 0.5_dp*DOT_PRODUCT(charges,MATMUL(qs_env%cp_ddapc_env%Md,charges))
          WRITE(iw,FMT="(T3,A,T60,F20.10)")"Decoupling Energy: ", e_decpl
       END IF
       IF (qs_env%cp_ddapc_ewald%do_qmmm_periodic.AND.(iw>0)) THEN
          e_recpl = 0.5_dp*DOT_PRODUCT(charges,MATMUL(qs_env%cp_ddapc_env%Mr,charges))
          WRITE(iw,FMT="(T3,A,T60,F20.10)")"Recoupling Energy: ", e_recpl
       END IF
       CALL modify_hartree_pot(v_hartree_gspace,&
                               density_fit_section,&
                               particle_set,&
                               qs_env%cp_ddapc_env%Mt,&
                               qs_env%cp_ddapc_env%AmI,&
                               radii,&
                               charges,&
                               error)
       ! Modify the Hartree potential due to the decoupling/recoupling
       energy = 0.5_dp * pw_integral_ab ( rho_tot_gspace%pw, v_hartree_gspace%pw,error=error)
       IF (need_f) THEN
          CALL ewald_ddapc_force(qs_env, qs_env%cp_ddapc_ewald%coeff_qm,&
               .FALSE., 1.0_dp, multipole_section, cell, particle_set,&
               radii, dq, charges, error=error)
          IF (qs_env%cp_ddapc_ewald%do_qmmm_periodic) THEN
             CALL ewald_ddapc_force(qs_env, qs_env%cp_ddapc_ewald%coeff_mm,&
                  .TRUE., -1.0_dp, multipole_section, super_cell, particle_set, &
                  radii, dq, charges, error=error)
          END IF
       END IF
       ! Clean the allocated arrays
       DEALLOCATE(charges, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (ASSOCIATED(dq)) THEN
          DEALLOCATE(dq, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL cp_print_key_finished_output(iw,logger,multipole_section,&
            "PROGRAM_RUN_INFO",error=error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE cp_ddapc_apply_CD

! *****************************************************************************
!> \brief Routine to apply RESTRAINT/CONSTRAINTS to the density
!>      with the Bloechl scheme
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE cp_ddapc_apply_RS(qs_env, energy, energy_res,v_hartree_gspace,&
       v_spin_ddapc_rest_g,section_id,ddapc_restraint_control,calculate_forces,&
       error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: energy_res
    TYPE(pw_p_type)                          :: v_hartree_gspace, &
                                                v_spin_ddapc_rest_g
    INTEGER                                  :: section_id
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_apply_RS', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, my_id, stat
    LOGICAL                                  :: apply_restrain, failure, &
                                                need_f
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges, radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dq
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: density_fit_section, &
                                                force_env_section, &
                                                restraint_section

    CALL timeset(routineN,handle)
    failure        = .FALSE.
    need_f         = .FALSE.
    IF (PRESENT(calculate_forces)) need_f = calculate_forces
    apply_restrain = qs_env%dft_control%qs_control%ddapc_restraint
    logger         => cp_error_get_logger(error)
    IF ((.NOT.failure).AND.apply_restrain) THEN
       ! Initialize
       NULLIFY(restraint_section, force_env_section, particle_set, charges, &
            radii, dq, cell, density_fit_section, super_cell)

       CALL get_qs_env(qs_env=qs_env,&
                       input=force_env_section,&
                       particle_set=particle_set,&
                       cell=cell,&
                       super_cell=super_cell,&
                       error=error)

       density_fit_section => section_vals_get_subs_vals(force_env_section,"DFT%DENSITY_FITTING",error=error)
       restraint_section => section_vals_get_subs_vals(force_env_section,"DFT%QS%DDAPC_RESTRAINT",&
                            error=error)
       iw=cp_print_key_unit_nr(logger,restraint_section,"PROGRAM_RUN_INFO",&
            extension=".fitChargeLog",error=error)
       ! First we evaluate the charges at the corresponding SCF STEP
       my_id= ddapc_restraint_control%density_type
       IF (need_f) THEN
          CALL get_ddapc(qs_env,&
                         need_f,&
                         density_fit_section,&
                         density_type=my_id,&
                         qout1=charges,&
                         out_radii=radii,&
                         dq_out=dq,&
                         error=error)
       ELSE
          CALL get_ddapc(qs_env,&
                         need_f,&
                         density_fit_section,&
                         density_type=my_id,&
                         qout1=charges,&
                         out_radii=radii,&
                         error=error)
       END IF

       ! Modify the Hartree potential due to the restrain or the v_spin_ddapc_rest_g
       IF ((my_id==do_spin_density).OR.qs_env%dft_control%qs_control%et_coupling_calc) THEN
          CALL restraint_functional_potential(v_spin_ddapc_rest_g, density_fit_section,&
                                           particle_set, qs_env%cp_ddapc_env%AmI, radii, charges, &
                                           ddapc_restraint_control, energy_res, error)
       ELSE
          CALL restraint_functional_potential(v_hartree_gspace, density_fit_section,&
                                           particle_set, qs_env%cp_ddapc_env%AmI, radii, charges, &
                                           ddapc_restraint_control, energy_res,error)
       ENDIF

       IF (need_f) THEN
          CALL restraint_functional_force(qs_env,&
                                          ddapc_restraint_control,&
                                          dq,&
                                          charges,&
                                          SIZE(radii),&
                                          particle_set,&
                                          error)
       END IF
       ! Clean the allocated arrays
       DEALLOCATE(charges, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (ASSOCIATED(dq)) THEN
          DEALLOCATE(dq, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL cp_print_key_finished_output(iw,logger,restraint_section,&
            "PROGRAM_RUN_INFO",error=error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE cp_ddapc_apply_RS

! *****************************************************************************
!> \brief Routine to apply a reaction field during SCF (SCRF) with the Bloechl scheme
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE cp_ddapc_apply_RF(qs_env,rho_tot_gspace,energy,&
       v_hartree_gspace,calculate_forces, Itype_of_density, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type)                          :: rho_tot_gspace
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    TYPE(pw_p_type)                          :: v_hartree_gspace
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    CHARACTER(LEN=*)                         :: Itype_of_density
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_apply_RF', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, stat
    LOGICAL                                  :: apply_solvation, failure, &
                                                need_f
    REAL(KINd=dp)                            :: e_recpl
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges, radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dq
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: density_fit_section, &
                                                force_env_section, &
                                                solvation_section

    CALL timeset(routineN,handle)
    failure         = .FALSE.
    need_f          = .FALSE.
    IF (PRESENT(calculate_forces)) need_f = calculate_forces
    logger          => cp_error_get_logger(error)
    apply_solvation = qs_env%cp_ddapc_ewald%do_solvation
    IF ((.NOT.failure).AND.(apply_solvation)) THEN
       ! Initialize
       NULLIFY(force_env_section, particle_set, charges, &
            radii, dq, cell, super_cell)

       CALL get_qs_env(qs_env=qs_env,&
                       input=force_env_section,&
                       particle_set=particle_set,&
                       cell=cell,&
                       super_cell=super_cell,&
                       error=error)

       solvation_section => section_vals_get_subs_vals(force_env_section,"DFT%SCRF", error=error)
       ! Start the real calculation
       iw=cp_print_key_unit_nr(logger,solvation_section,"PROGRAM_RUN_INFO",&
            extension=".fitChargeLog",error=error)
       density_fit_section =>  section_vals_get_subs_vals(force_env_section,"DFT%DENSITY_FITTING", error=error)
       ! First we evaluate the charges at the corresponding SCF STEP
       IF (need_f) THEN
          CALL get_ddapc(qs_env,&
                         need_f,&
                         density_fit_section,&
                         qout1=charges,&
                         out_radii=radii,&
                         dq_out=dq,&
                         ext_rho_tot_g=rho_tot_gspace%pw,&
                         Itype_of_density=Itype_of_density,&
                         error=error)
       ELSE
          CALL get_ddapc(qs_env,&
                         need_f,&
                         density_fit_section,&
                         qout1=charges,&
                         out_radii=radii,&
                         ext_rho_tot_g=rho_tot_gspace%pw,&
                         Itype_of_density=Itype_of_density,&
                         error=error)
       END IF
       ! Evaluate the Ewald contribution to the decoupling/coupling E2 and E3
       IF (iw>0) THEN
          e_recpl = 0.5_dp*DOT_PRODUCT(charges,MATMUL(qs_env%cp_ddapc_env%Ms,charges))
          WRITE(iw,FMT="(T3,A,T60,F20.10)")"Solvation  Energy: ", e_recpl
       END IF
       CALL modify_hartree_pot(v_hartree_gspace,&
                               density_fit_section,&
                               particle_set,&
                               qs_env%cp_ddapc_env%Ms,&
                               qs_env%cp_ddapc_env%AmI,&
                               radii,&
                               charges,&
                               error)
       ! Modify the Hartree potential due to the reaction field
       energy = 0.5_dp * pw_integral_ab ( rho_tot_gspace%pw, v_hartree_gspace%pw, error=error)
       IF (need_f) THEN
          CALL solvation_ddapc_force(qs_env, cell, solvation_section, particle_set,&
               radii, dq, charges, error)
       END IF
       ! Clean the allocated arrays
       DEALLOCATE(charges, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (ASSOCIATED(dq)) THEN
          DEALLOCATE(dq, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL cp_print_key_finished_output(iw,logger,solvation_section,&
            "PROGRAM_RUN_INFO",error=error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE cp_ddapc_apply_RF

END MODULE cp_ddapc
