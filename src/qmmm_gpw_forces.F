!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines to compute energy and forces in a QM/MM calculation
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_gpw_forces
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_spline_utils,                 ONLY: pw_restrict_s3
  USE cube_utils,                      ONLY: cube_info_type
  USE input_constants,                 ONLY: do_par_atom,&
                                             do_qmmm_coulomb,&
                                             do_qmmm_gauss,&
                                             do_qmmm_none,&
                                             do_qmmm_pcharge,&
                                             do_qmmm_swave,&
                                             spline3_nopbc_interp,&
                                             spline3_pbc_interp
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_irecv,&
                                             mp_isend,&
                                             mp_sum,&
                                             mp_wait
  USE mm_collocate_potential,          ONLY: collocate_gf_rspace_NoPBC,&
                                             integrate_gf_rspace_NoPBC
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_integral_ab,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_gpw_energy,                 ONLY: qmmm_elec_with_gaussian,&
                                             qmmm_elec_with_gaussian_LG,&
                                             qmmm_elec_with_gaussian_LR
  USE qmmm_se_forces,                  ONLY: deriv_se_qmmm_matrix
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_per_pot_p_type,&
                                             qmmm_per_pot_type,&
                                             qmmm_pot_p_type,&
                                             qmmm_pot_type
  USE qmmm_util,                       ONLY: spherical_cutoff_factor
  USE qs_dftb_matrices,                ONLY: deriv_dftb_qmmm_matrix,&
                                             deriv_dftb_qmmm_matrix_pc
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  LOGICAL,       PARAMETER, PRIVATE    :: debug_this_module=.FALSE.
  REAL(KIND=dp), PARAMETER, PRIVATE    :: Dx = 0.01_dp     ! Debug Variables
  REAL(KIND=dp), PARAMETER, PRIVATE    :: MaxErr = 10.0_dp ! Debug Variables
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_gpw_forces'
  PUBLIC :: qmmm_forces

CONTAINS

! *****************************************************************************
!> \brief General driver to Compute the contribution
!>      to the forces due to the QM/MM potential
!> \param qs_env ...
!> \param qmmm_env ...
!> \param mm_particles ...
!> \param calc_force ...
!> \param mm_cell ...
!> \param error ...
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_forces(qs_env,qmmm_env,mm_particles,calc_force,mm_cell,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, image_IndMM, &
                                                Imm, IndMM, ispin, iw, stat
    LOGICAL                                  :: failure, gapw, need_f, &
                                                periodic
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces, Forces_added_charges
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_r
    TYPE(pw_p_type), POINTER                 :: rho0_s_gs, rho_core
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pool
    TYPE(pw_type), POINTER                   :: rho_tot_r, rho_tot_r2
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input_section, &
                                                interp_section, print_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    need_f = .TRUE.
    periodic = qmmm_env%periodic
    IF (PRESENT(calc_force)) need_f = calc_force
    NULLIFY(dft_control, ks_qmmm_env_loc, rho, pw_env, rho_tot_r, energy, Forces,&
            Forces_added_charges,input_section,rho0_s_gs, rho_r)
    CALL get_qs_env(qs_env=qs_env,&
         rho=rho,&
         rho_core=rho_core,&
         pw_env=pw_env,&
         energy=energy,&
         para_env=para_env,&
         input=input_section,&
         rho0_s_gs=rho0_s_gs,&
         dft_control=dft_control,&
         error=error)

    CALL qs_rho_get(rho, rho_r=rho_r, error=error)

    logger => cp_error_get_logger(error)
    ks_qmmm_env_loc => qs_env%ks_qmmm_env
    interp_section => section_vals_get_subs_vals(input_section,"QMMM%INTERPOLATOR",error=error)
    print_section => section_vals_get_subs_vals(input_section,"QMMM%PRINT",error=error)
    iw=cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",&
         extension=".qmmmLog",error=error)
    gapw = dft_control%qs_control%gapw
    ! If forces are required allocate these temporary arrays
    IF (need_f) THEN
       ALLOCATE(Forces(3,qmmm_env%num_mm_atoms),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Forces_added_charges(3,qmmm_env%added_charges%num_mm_atoms),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       Forces(:,:)               = 0.0_dp
       Forces_added_charges(:,:) = 0.0_dp
    END IF
    IF (dft_control%qs_control%semi_empirical) THEN
       ! SEMIEMPIRICAL
       SELECT CASE(qmmm_env%qmmm_coupl_type)
       CASE(do_qmmm_coulomb)
          CALL deriv_se_qmmm_matrix(qs_env,qmmm_env,mm_particles,mm_cell,para_env,&
               need_f,Forces,Forces_added_charges,error)
       CASE(do_qmmm_pcharge)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Point Charge QM/MM electrostatic coupling not yet implemented for SE.", &
               error=error, error_level=cp_failure_level)
       CASE(do_qmmm_gauss,do_qmmm_swave)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="GAUSS or SWAVE QM/MM electrostatic coupling not yet implemented for SE.", &
               error=error, error_level=cp_failure_level)
       CASE(do_qmmm_none)
          IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')&
               "- No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
       CASE DEFAULT
          IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')"Unknown Coupling..."
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    ELSEIF (dft_control%qs_control%dftb) THEN
       ! DFTB
       SELECT CASE(qmmm_env%qmmm_coupl_type)
       CASE(do_qmmm_none)
          IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')&
               "- No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
       CASE(do_qmmm_coulomb)
          CALL deriv_dftb_qmmm_matrix(qs_env,qmmm_env,mm_particles,mm_cell,para_env,&
               need_f,Forces,Forces_added_charges,error)
       CASE(do_qmmm_pcharge)
          CALL deriv_dftb_qmmm_matrix_pc(qs_env,qmmm_env,mm_particles,mm_cell,para_env,&
               need_f,Forces,Forces_added_charges,error)
       CASE(do_qmmm_gauss,do_qmmm_swave)
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="GAUSS or SWAVE QM/MM electrostatic coupling not yet implemented for DFTB.", &
               error=error, error_level=cp_failure_level)
       CASE DEFAULT
          IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')"Unknown Coupling..."
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    ELSE
       ! GPW/GAPW
       CALL pw_env_get(pw_env=pw_env,&
                       pw_pools=pw_pools,&
                       auxbas_pw_pool=auxbas_pool,&
                       error=error)
       CALL pw_pool_create_pw(auxbas_pool,rho_tot_r,&
                              in_space=REALSPACE,&
                              use_data=REALDATA3D,&
                              error=error)
       ! IF GAPW the core charge is replaced by the compensation charge
       IF(gapw) THEN
          IF( dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
            CALL pw_transfer(rho_core%pw,rho_tot_r,error=error)
            energy%qmmm_nu = pw_integral_ab ( rho_tot_r, ks_qmmm_env_loc%v_qmmm_rspace%pw,error=error)
            CALL pw_pool_create_pw(auxbas_pool,rho_tot_r2,&
                              in_space=REALSPACE,&
                              use_data=REALDATA3D,&
                              error=error)
            CALL pw_transfer(rho0_s_gs%pw,rho_tot_r2,error=error)
            CALL pw_axpy(rho_tot_r2,rho_tot_r,error=error)
            CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_r2,accept_non_compatible=.TRUE.,error=error)
          ELSE
            CALL pw_transfer(rho0_s_gs%pw,rho_tot_r,error=error)
          !
          ! QM/MM Nuclear Electrostatic Potential already included through rho0
          !
            energy%qmmm_nu = 0.0_dp
          END IF
       ELSE
          CALL pw_transfer(rho_core%pw,rho_tot_r,error=error)
          !
          ! Computes the QM/MM Nuclear Electrostatic Potential
          !
          energy%qmmm_nu = pw_integral_ab ( rho_tot_r, ks_qmmm_env_loc%v_qmmm_rspace%pw,error=error)
       END IF
       IF (need_f) THEN
          !
          DO ispin=1,SIZE(rho_r)
             CALL pw_axpy(rho_r(ispin)%pw,rho_tot_r,error=error)
          END DO
          IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')"Evaluating forces on MM atoms due to the:"
          ! Electrostatic Interaction type...
          SELECT CASE(qmmm_env%qmmm_coupl_type)
          CASE(do_qmmm_coulomb)
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Coulomb QM/MM electrostatic coupling not implemented for GPW/GAPW.", &
                  error=error, error_level=cp_failure_level)
          CASE(do_qmmm_pcharge)
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Point Charge QM/MM electrostatic coupling not yet implemented for GPW/GAPW.", &
                  error=error, error_level=cp_failure_level)
          CASE(do_qmmm_gauss,do_qmmm_swave)
             IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')&
                  "- QM/MM Coupling computed collocating the Gaussian Potential Functions."
             CALL qmmm_forces_with_gaussian(rho=rho_tot_r,&
                                            qmmm_env=qmmm_env,&
                                            mm_particles=mm_particles,&
                                            aug_pools=qmmm_env%aug_pools,&
                                            auxbas_grid=qmmm_env%gridlevel_info%auxbas_grid,&
                                            coarser_grid=qmmm_env%gridlevel_info%coarser_grid,&
                                            para_env=para_env,&
                                            pw_pools=pw_pools,&
                                            eps_mm_rspace=qmmm_env%eps_mm_rspace,&
                                            cube_info=ks_qmmm_env_loc%cube_info,&
                                            Forces=Forces,&
                                            Forces_added_charges=Forces_added_charges,&
                                            interp_section=interp_section,&
                                            iw=iw,&
                                            mm_cell=mm_cell,&
                                            error=error)
          CASE(do_qmmm_none)
             IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')&
                  "- No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
          CASE DEFAULT
             IF (iw>0) WRITE(iw,'(T2,"QMMM|",1X,A)')"- Unknown Coupling..."
 CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT
       END IF
    END IF
    ! Correct Total Energy adding the contribution of the QM/MM nuclear interaction
    energy%total = energy%total +  energy%qmmm_nu
    ! Proceed if gradients are requested..
    IF (need_f) THEN
       !ikuo Temporary change to alleviate compiler problems on Intel with
       !array dimension of 0
       IF(qmmm_env%num_mm_atoms/=0)               CALL mp_sum(Forces,para_env%group)
       IF(qmmm_env%added_charges%num_mm_atoms/=0) CALL mp_sum(Forces_added_charges,para_env%group)
       ! Debug Forces
       IF (debug_this_module) THEN
          IF ( (dft_control%qs_control%semi_empirical).OR.&
               (dft_control%qs_control%dftb))THEN
             WRITE(iw,*)"NO DEBUG AVAILABLE in module"//TRIM(routineN)
          ELSE
             ! Print Out Forces
             IF (iw>0) THEN
                DO Imm = 1, SIZE(qmmm_env%mm_atom_index)
                   WRITE(iw,*)"ANALYTICAL FORCES:"
                   IndMM = qmmm_env%mm_atom_index(Imm)
                   WRITE(iw,'(I6,3F15.9)') IndMM, Forces(:,Imm)
                END DO
             END IF
             CALL qmmm_debug_forces(rho=rho_tot_r,&
                                    qs_env=qs_env,&
                                    qmmm_env=qmmm_env,&
                                    Analytical_Forces=Forces,&
                                    mm_particles=mm_particles,&
                                    mm_atom_index=qmmm_env%mm_atom_index,&
                                    num_mm_atoms=qmmm_env%num_mm_atoms,&
                                    interp_section=interp_section,&
                                    mm_cell=mm_cell,&
                                    error=error)
          END IF
       ENDIF
    END IF
    ! Give back rho_tot_t to auxbas_pool only for GPW/GAPW
    IF  ((.NOT.dft_control%qs_control%semi_empirical).AND.&
         (.NOT.dft_control%qs_control%dftb)) THEN
       CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_r,&
            accept_non_compatible=.TRUE.,error=error)
    END IF
    IF (iw>0) THEN
       IF(.NOT. gapw) WRITE(iw,'(T2,"QMMM|",1X,A,T66,F15.9)')&
            "QM/MM Nuclear Electrostatic Potential :",energy%qmmm_nu
       WRITE(iw,'(T2,"QMMM|",1X,A,T66,F15.9)')&
            "QMMM Total Energy (QM + QMMM electronic + QMMM nuclear):",energy%total
       WRITE(iw,'(T2,"QMMM|",1X,A)')"MM energy NOT included in the above term!"//&
            " Check for:  FORCE_EVAL ( QMMM )"
       WRITE(iw,'(T2,"QMMM|",1X,A)')"that includes both QM, QMMM and MM energy terms!"
    END IF
    IF (need_f) THEN
       ! Transfer Forces
       DO Imm = 1, qmmm_env%num_mm_atoms
          IndMM = qmmm_env%mm_atom_index(Imm)

          !add image forces to Forces
          IF(qmmm_env%image_charge) THEN
             DO iatom=1,qmmm_env%num_image_mm_atoms
               image_IndMM=qmmm_env%image_charge_pot%image_mm_list(iatom)
               IF(image_IndMM.eq.IndMM) THEN
               Forces(:,Imm)=Forces(:,Imm)&
                             +qmmm_env%image_charge_pot%image_forcesMM(:,iatom)
               ENDIF
             ENDDO
          ENDIF

          ! Hack: In Forces there the gradients indeed...
          ! Minux sign to take care of this misunderstanding...
          mm_particles(IndMM)%f(:) = - Forces(:,Imm) + mm_particles(IndMM)%f(:)
       END DO
       DEALLOCATE(Forces, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
          DO Imm = 1, qmmm_env%added_charges%num_mm_atoms
             IndMM = qmmm_env%added_charges%mm_atom_index(Imm)
             ! Hack: In Forces there the gradients indeed...
             ! Minux sign to take care of this misunderstanding...
             qmmm_env%added_charges%added_particles(IndMM)%f(:) = - Forces_added_charges(:,Imm)
          END DO
       END IF
       DEALLOCATE(Forces_added_charges, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL cp_print_key_finished_output(iw,logger,print_section,"PROGRAM_RUN_INFO",&
         error=error)
    CALL timestop(handle)

  END SUBROUTINE qmmm_forces

! *****************************************************************************
!> \brief Evaluates the contribution to the forces due to the
!>      QM/MM potential computed collocating the Electrostatic
!>      Gaussian Potential.
!> \param rho ...
!> \param qmmm_env ...
!> \param mm_particles ...
!> \param aug_pools ...
!> \param auxbas_grid ...
!> \param coarser_grid ...
!> \param cube_info ...
!> \param para_env ...
!> \param eps_mm_rspace ...
!> \param pw_pools ...
!> \param Forces ...
!> \param Forces_added_charges ...
!> \param interp_section ...
!> \param iw ...
!> \param mm_cell ...
!> \param error ...
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_forces_with_gaussian(rho, qmmm_env, mm_particles,  &
       aug_pools, auxbas_grid, coarser_grid, cube_info, para_env,    &
       eps_mm_rspace, pw_pools, Forces, Forces_added_charges,        &
       interp_section, iw, mm_cell, error)
    TYPE(pw_type), POINTER                   :: rho
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    INTEGER, INTENT(IN)                      :: auxbas_grid, coarser_grid
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces, Forces_added_charges
    TYPE(section_vals_type), POINTER         :: interp_section
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_with_gaussian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, handle, i, igrid, j, &
                                                k, kind_interp, ngrids, &
                                                request, stat
    INTEGER, DIMENSION(3)                    :: glb, gub, lb, ub
    INTEGER, DIMENSION(:), POINTER           :: pos_of_x
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tmp
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids

! Statements

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(grids,tmp)
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qmmm_env%mm_atom_chrg),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qmmm_env%mm_atom_index),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(Forces),cp_failure_level,routineP,error,failure)
    !Statements
    ngrids=SIZE(pw_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw,error=error)
    END DO
    ! Collocate Density on multigrids
    lb = rho%pw_grid%bounds_local(1,:)
    ub = rho%pw_grid%bounds_local(2,:)
    grids(auxbas_grid) % pw % cr3d(lb(1):ub(1),&
                                   lb(2):ub(2),&
                                   lb(3):ub(3) )= rho % cr3d
    ! copy the boundaries
    DO i=lb(1),ub(1)
       grids(auxbas_grid) % pw %cr3d(i,ub(2)+1,ub(3)+1)=rho%cr3d(i,lb(2),lb(3))
    END DO
    DO k=lb(3),ub(3)
       DO i=lb(1),ub(1)
          grids(auxbas_grid) % pw %cr3d(i,ub(2)+1,k)=rho%cr3d(i,lb(2),k)
       END DO
    END DO
    DO j=lb(2),ub(2)
       DO i=lb(1),ub(1)
          grids(auxbas_grid) % pw %cr3d(i,j,ub(3)+1)=rho%cr3d(i,j,lb(3))
       END DO
    END DO
    pos_of_x => grids(auxbas_grid)%pw%pw_grid%para%pos_of_x
    group=grids(auxbas_grid)%pw%pw_grid%para%group
    glb=rho%pw_grid%bounds(1,:)
    gub=rho%pw_grid%bounds(2,:)
    IF (lb(1)==glb(1).AND.ub(1)==gub(1)) THEN
       DO k=lb(3),ub(3)
          DO j=lb(2),ub(2)
             grids(auxbas_grid) % pw %cr3d(ub(1)+1,j,k)=rho%cr3d(lb(1),j,k)
          END DO
          grids(auxbas_grid) % pw %cr3d(ub(1)+1,ub(2)+1,k)=rho%cr3d(lb(1),lb(2),k)
       END DO
       DO j=lb(2),ub(2)
          grids(auxbas_grid) % pw %cr3d(ub(1)+1,j,ub(3)+1)=rho%cr3d(lb(1),j,lb(3))
       END DO
       grids(auxbas_grid) % pw %cr3d(ub(1)+1,ub(2)+1,ub(3)+1)=rho%cr3d(lb(1),lb(2),lb(3))
    ELSE IF (lb(1)==glb(1)) THEN
       ALLOCATE(tmp(rho%pw_grid%bounds_local(1,2):rho%pw_grid%bounds_local(2,2),&
            rho%pw_grid%bounds_local(1,3):rho%pw_grid%bounds_local(2,3)),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       tmp=rho%cr3d(lb(1),:,:)
       CALL mp_isend(msgin=tmp,dest=pos_of_x(rho%pw_grid%bounds(2,1)),comm=group,&
            request=request,tag=112)
       CALL mp_wait(request)
    ELSE IF (ub(1)==gub(1)) THEN
       ALLOCATE(tmp(rho%pw_grid%bounds_local(1,2):rho%pw_grid%bounds_local(2,2),&
            rho%pw_grid%bounds_local(1,3):rho%pw_grid%bounds_local(2,3)),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL mp_irecv(msgout=tmp,source=pos_of_x(rho%pw_grid%bounds(1,1)),&
            comm=group,request=request,tag=112)
       CALL mp_wait(request)

       DO k=lb(3),ub(3)
          DO j=lb(2),ub(2)
             grids(auxbas_grid) % pw %cr3d(ub(1)+1,j,k)=tmp(j,k)
          END DO
          grids(auxbas_grid) % pw %cr3d(ub(1)+1,ub(2)+1,k)=tmp(lb(2),k)
       END DO
       DO j=lb(2),ub(2)
          grids(auxbas_grid) % pw %cr3d(ub(1)+1,j,ub(3)+1)=tmp(j,lb(3))
       END DO
       grids(auxbas_grid) % pw %cr3d(ub(1)+1,ub(2)+1,ub(3)+1)=tmp(lb(2),lb(3))
    END IF
    IF (ASSOCIATED(tmp)) THEN
       DEALLOCATE(tmp, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       NULLIFY(tmp)
    END IF
    ! Further setup of parallelization scheme
    IF (qmmm_env%par_scheme==do_par_atom) THEN
       CALL mp_sum(grids(auxbas_grid)%pw%cr3d,para_env%group)
    END IF
    ! RealSpace Interpolation
    CALL section_vals_val_get(interp_section,"kind", i_val=kind_interp, error=error)
    SELECT CASE(kind_interp)
    CASE(spline3_nopbc_interp, spline3_pbc_interp)
       ! Spline Interpolator
       DO Igrid =  auxbas_grid, SIZE(grids)-1
          CALL pw_restrict_s3(grids(Igrid  )%pw,&
                              grids(Igrid+1)%pw,&
                              aug_pools(Igrid+1)%pool,&
                              param_section=interp_section,&
                              error=error)
       END DO
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL qmmm_force_with_gaussian_low(grids, mm_particles,&
         qmmm_env%mm_atom_chrg, qmmm_env%mm_el_pot_radius, qmmm_env%mm_atom_index,&
         qmmm_env%num_mm_atoms, cube_info, para_env, eps_mm_rspace, auxbas_grid, &
         coarser_grid, qmmm_env%pgfs, qmmm_env%potentials, Forces, aug_pools,    &
         mm_cell, qmmm_env%dOmmOqm, qmmm_env%periodic, qmmm_env%per_potentials,  &
         iw, qmmm_env%par_scheme, qmmm_env%spherical_cutoff, error)

    IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
       CALL qmmm_force_with_gaussian_low(grids, qmmm_env%added_charges%added_particles,&
            qmmm_env%added_charges%mm_atom_chrg, qmmm_env%added_charges%mm_el_pot_radius,&
            qmmm_env%added_charges%mm_atom_index, qmmm_env%added_charges%num_mm_atoms,&
            cube_info, para_env, eps_mm_rspace, auxbas_grid, coarser_grid, qmmm_env%added_charges%pgfs,&
            qmmm_env%added_charges%potentials, Forces_added_charges, aug_pools, mm_cell, &
            qmmm_env%dOmmOqm, qmmm_env%periodic, qmmm_env%per_potentials, iw, qmmm_env%par_scheme,&
            qmmm_env%spherical_cutoff, error)
    END IF

    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)
    CALL timestop(handle)

  END SUBROUTINE qmmm_forces_with_gaussian

! *****************************************************************************
!> \brief Evaluates the contribution to the forces due to the
!>      QM/MM potential computed collocating the Electrostatic
!>      Gaussian Potential. Low Level
!> \param grids ...
!> \param mm_particles ...
!> \param mm_charges ...
!> \param mm_el_pot_radius ...
!> \param mm_atom_index ...
!> \param num_mm_atoms ...
!> \param cube_info ...
!> \param para_env ...
!> \param eps_mm_rspace ...
!> \param auxbas_grid ...
!> \param coarser_grid ...
!> \param pgfs ...
!> \param potentials ...
!> \param Forces ...
!> \param aug_pools ...
!> \param mm_cell ...
!> \param dOmmOqm ...
!> \param periodic ...
!> \param per_potentials ...
!> \param iw ...
!> \param par_scheme ...
!> \param qmmm_spherical_cutoff ...
!> \param error ...
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_force_with_gaussian_low(grids, mm_particles, mm_charges, &
       mm_el_pot_radius, mm_atom_index, num_mm_atoms, cube_info, para_env, &
       eps_mm_rspace, auxbas_grid, coarser_grid, pgfs, potentials, Forces, &
       aug_pools, mm_cell, dOmmOqm, periodic, per_potentials, iw, par_scheme,&
       qmmm_spherical_cutoff, error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    INTEGER, INTENT(IN)                      :: auxbas_grid, coarser_grid
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: Potentials
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    LOGICAL, INTENT(in)                      :: periodic
    TYPE(qmmm_per_pot_p_type), &
      DIMENSION(:), POINTER                  :: per_potentials
    INTEGER, INTENT(IN)                      :: iw, par_scheme
    REAL(KIND=dp), INTENT(IN)                :: qmmm_spherical_cutoff(2)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_force_with_gaussian_low', &
      routineNb = 'qmmm_forces_gaussian_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, IGauss, &
                                                ilevel, Imm, IndMM, IRadTyp, &
                                                LIndMM, myind, n_rep_real(3), &
                                                stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, dvol, height, &
                                                sph_chrg_factor, W
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: xdat, ydat, zdat
    REAL(KIND=dp), DIMENSION(3)              :: force, ra
    TYPE(qmmm_gaussian_type), POINTER        :: pgf
    TYPE(qmmm_per_pot_type), POINTER         :: per_pot
    TYPE(qmmm_pot_type), POINTER             :: pot

    CALL timeset(routineN,handle)
    CALL timeset(routineNb//"_G",handle2)
    NULLIFY(pgf, pot, per_pot)
    IF (par_scheme==do_par_atom) myind=0
    Radius: DO IRadTyp = 1, SIZE(pgfs)
       pgf => pgfs(IRadTyp)%pgf
       pot => potentials(IRadTyp)%pot
       n_rep_real = 0
       IF (periodic) THEN
          per_pot => per_potentials(IRadTyp)%pot
          n_rep_real = per_pot%n_rep_real
       END IF
       Gaussian: DO IGauss = 1, pgf%Number_of_Gaussians
          alpha     = 1.0_dp / pgf%Gk(IGauss)
          alpha     = alpha * alpha
          height    = pgf%Ak(IGauss)
          ilevel    = pgf%grid_level(IGauss)
          dvol      = grids(ilevel)%pw%pw_grid%dvol
          bo        = grids(ilevel) % pw % pw_grid % bounds_local
          ALLOCATE (xdat(2,bo(1,1):bo(2,1)), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (ydat(2,bo(1,2):bo(2,2)), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE (zdat(2,bo(1,3):bo(2,3)), STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          Atoms: DO Imm = 1, SIZE(pot%mm_atom_index)
             IF (par_scheme==do_par_atom) THEN
                myind = myind + 1
                IF (MOD(myind,para_env%num_pe)/=para_env%mepos) CYCLE Atoms
             END IF
             LIndMM    =   pot%mm_atom_index(Imm)
             IndMM     =   mm_atom_index(LIndMM)
             ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm,mm_cell)+dOmmOqm
             W         =   mm_charges(LIndMM) * height
             force     =   0.0_dp
             ! Possible Spherical Cutoff
             IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
                CALL spherical_cutoff_factor(qmmm_spherical_cutoff, ra, sph_chrg_factor, error)
                W = W * sph_chrg_factor
             END IF
             IF (ABS(W)<= EPSILON(0.0_dp)) CYCLE Atoms
             CALL integrate_gf_rspace_NoPBC(zetp=alpha,&
                                            rp=ra,&
                                            scale=-1.0_dp,&
                                            W=W,&
                                            pwgrid=grids(ilevel)%pw,&
                                            cube_info=cube_info(ilevel),&
                                            eps_mm_rspace=eps_mm_rspace,&
                                            xdat=xdat,&
                                            ydat=ydat,&
                                            zdat=zdat,&
                                            bo=bo,&
                                            force=force,&
                                            n_rep_real=n_rep_real,&
                                            mm_cell=mm_cell)
             force = force * dvol
             Forces(:,LIndMM) = Forces(:,LIndMM) + force(:)
             !
             ! Debug Statement
             !
             IF (debug_this_module) THEN
                CALL debug_integrate_gf_rspace_NoPBC(ilevel=ilevel,&
                                                     zetp=alpha,&
                                                     rp=ra,&
                                                     W=W,&
                                                     pwgrid=grids(ilevel)%pw,&
                                                     cube_info=cube_info(ilevel),&
                                                     eps_mm_rspace=eps_mm_rspace,&
                                                     aug_pools=aug_pools,&
                                                     debug_force=force,&
                                                     mm_cell=mm_cell,&
                                                     auxbas_grid=auxbas_grid,&
                                                     n_rep_real=n_rep_real,&
                                                     iw=iw,&
                                                     error=error)
             END IF
          END DO Atoms
          DEALLOCATE (xdat, STAT=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (ydat, STAT=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE (zdat, STAT=stat)
 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END DO Gaussian
    END DO Radius
    CALL timestop(handle2)
    CALL timeset(routineNb//"_R",handle2)
    IF (periodic) THEN
       CALL qmmm_forces_with_gaussian_LG  (pgfs=pgfs,&
                                           cgrid=grids(coarser_grid)%pw,&
                                           num_mm_atoms=num_mm_atoms,&
                                           mm_charges=mm_charges,&
                                           mm_atom_index=mm_atom_index,&
                                           mm_particles=mm_particles,&
                                           mm_el_pot_radius=mm_el_pot_radius,&
                                           para_env=para_env,&
                                           coarser_grid_level=coarser_grid,&
                                           Forces=Forces,&
                                           per_potentials=per_potentials,&
                                           aug_pools=aug_pools,&
                                           mm_cell=mm_cell,&
                                           dOmmOqm=dOmmOqm,&
                                           iw=iw,&
                                           par_scheme=par_scheme,&
                                           qmmm_spherical_cutoff=qmmm_spherical_cutoff,&
                                           error=error)
    ELSE
       CALL qmmm_forces_with_gaussian_LR  (pgfs=pgfs,&
                                           cgrid=grids(coarser_grid)%pw,&
                                           num_mm_atoms=num_mm_atoms,&
                                           mm_charges=mm_charges,&
                                           mm_atom_index=mm_atom_index,&
                                           mm_particles=mm_particles,&
                                           mm_el_pot_radius=mm_el_pot_radius,&
                                           para_env=para_env,&
                                           coarser_grid_level=coarser_grid,&
                                           Forces=Forces,&
                                           potentials=potentials,&
                                           aug_pools=aug_pools,&
                                           mm_cell=mm_cell,&
                                           dOmmOqm=dOmmOqm,&
                                           iw=iw,&
                                           par_scheme=par_scheme,&
                                           qmmm_spherical_cutoff=qmmm_spherical_cutoff,&
                                           error=error)
    END IF
    CALL timestop(handle2)
    CALL timestop(handle)
  END SUBROUTINE qmmm_force_with_gaussian_low

! *****************************************************************************
!> \brief Evaluates the contribution to the forces due to the Long Range
!>      part of the QM/MM potential computed collocating the Electrostatic
!>      Gaussian Potential.
!> \param pgfs ...
!> \param cgrid ...
!> \param num_mm_atoms ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_particles ...
!> \param mm_el_pot_radius ...
!> \param para_env ...
!> \param coarser_grid_level ...
!> \param Forces ...
!> \param per_potentials ...
!> \param aug_pools ...
!> \param mm_cell ...
!> \param dOmmOqm ...
!> \param iw ...
!> \param par_scheme ...
!> \param qmmm_spherical_cutoff ...
!> \param error ...
!> \par History
!>      08.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_forces_with_gaussian_LG  (pgfs,cgrid,num_mm_atoms,mm_charges,mm_atom_index,&
       mm_particles,mm_el_pot_radius,para_env, coarser_grid_level,Forces, per_potentials,&
       aug_pools, mm_cell, dOmmOqm, iw, par_scheme, qmmm_spherical_cutoff, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_type), POINTER                   :: cgrid
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: coarser_grid_level
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces
    TYPE(qmmm_per_pot_p_type), &
      DIMENSION(:), POINTER                  :: per_potentials
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    INTEGER, INTENT(IN)                      :: iw, par_scheme
    REAL(KIND=dp), DIMENSION(2), INTENT(IN)  :: qmmm_spherical_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_with_gaussian_LG', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, ii1, ii2, ii3, ii4, ij1, ij2, ij3, ij4, ik1, ik2, &
      ik3, ik4, Imm, IndMM, IRadTyp, ivec(3), j, k, LIndMM, my_i, my_j, my_k, &
      myind, npts(3), stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: a1, a2, a3, abc_X(4,4), abc_X_Y(4), b1, b2, b3, c1, c2, &
      c3, d1, d2, d3, dr1, dr1c, dr1i, dr2, dr2c, dr2i, dr3, dr3c, dr3i, &
      dvol, e1, e2, e3, f1, f2, f3, fac, ft1, ft2, ft3, g1, g2, g3, h1, h2, &
      h3, p1, p2, p3, q1, q2, q3, qt, r1, r2, r3, rt1, rt2, rt3, rv1, rv2, &
      rv3, s1, s1d, s1o, s2, s2d, s2o, s3, s3d, s3o, s4, s4d, s4o, &
      sph_chrg_factor, t1, t1d, t1o, t2, t2d, t2o, t3, t3d, t3o, t4, t4d, &
      t4o, u1, u2, u3, v1, v1d, v1o, v2, v2d, v2o, v3, v3d, v3o, v4, v4d, &
      v4o, xd1, xd2, xd3, xs1, xs2, xs3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: LForces
    REAL(KIND=dp), DIMENSION(3)              :: ra, val, vec
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid, grid2
    TYPE(pw_type), POINTER                   :: pw
    TYPE(qmmm_per_pot_type), POINTER         :: per_pot

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(grid)
    ALLOCATE(LForces(3,num_mm_atoms), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    LForces = 0.0_dp
    dr1c  = cgrid%pw_grid%dr(1)
    dr2c  = cgrid%pw_grid%dr(2)
    dr3c  = cgrid%pw_grid%dr(3)
    dvol  = cgrid%pw_grid%dvol
    gbo   = cgrid%pw_grid%bounds
    bo    =  cgrid%pw_grid%bounds_local
    grid  => cgrid%cr3d
    IF (par_scheme==do_par_atom) myind = 0
    Radius: DO IRadTyp = 1, SIZE(pgfs)
       per_pot => per_potentials(IRadTyp)%pot
       pw    => per_pot%TabLR
       grid2 => pw%cr3d(:,:,:)
       npts = pw%pw_grid%npts
       dr1  = pw%pw_grid%dr(1)
       dr2  = pw%pw_grid%dr(2)
       dr3  = pw%pw_grid%dr(3)
       dr1i = 1.0_dp / dr1
       dr2i = 1.0_dp / dr2
       dr3i = 1.0_dp / dr3
       Atoms: DO Imm = 1, SIZE(per_pot%mm_atom_index)
          IF (par_scheme==do_par_atom) THEN
             myind = myind + 1
             IF (MOD(myind,para_env%num_pe)/=para_env%mepos) CYCLE
          END IF
          LIndMM    =   per_pot%mm_atom_index(Imm)
          IndMM     =   mm_atom_index(LIndMM)
          ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm,mm_cell)+dOmmOqm
          qt        =   mm_charges(LIndMM)
          ! Possible Spherical Cutoff
          IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
             CALL spherical_cutoff_factor(qmmm_spherical_cutoff, ra, sph_chrg_factor, error)
             qt = qt * sph_chrg_factor
          END IF
          IF (ABS(qt)<= EPSILON(0.0_dp)) CYCLE Atoms
          rt1 = ra(1)
          rt2 = ra(2)
          rt3 = ra(3)
          ft1 = 0.0_dp
          ft2 = 0.0_dp
          ft3 = 0.0_dp
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k = k-gbo(1,3)
             xs3  = REAL(my_k,dp)*dr3c
             my_j = bo(1,2)-gbo(1,2)
             xs2 = REAL(my_j,dp)*dr2c
             rv3 = rt3 - xs3
             vec(3) = rv3
             ivec(3) = FLOOR(vec(3)/pw%pw_grid%dr(3))
                   ik1 = MODULO(ivec(3)-1,npts(3))+1
                   ik2 = MODULO(ivec(3)  ,npts(3))+1
                   ik3 = MODULO(ivec(3)+1,npts(3))+1
                   ik4 = MODULO(ivec(3)+2,npts(3))+1
             xd3  = (vec(3)/dr3)-REAL(ivec(3),kind=dp)
             p1  = 3.0_dp + xd3
             p2  = p1*p1
             p3  = p2*p1
             q1  = 2.0_dp + xd3
             q2  = q1*q1
             q3  = q2*q1
             r1  = 1.0_dp + xd3
             r2  = r1*r1
             r3  = r2*r1
             u1  = xd3
             u2  = u1*u1
             u3  = u2*u1
             v1o =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*p1 + 12.0_dp*p2 - p3)
             v2o = -22.0_dp/3.0_dp + 10.0_dp*q1 - 4.0_dp*q2 + 0.5_dp*q3
             v3o =   2.0_dp/3.0_dp -  2.0_dp*r1 + 2.0_dp*r2 - 0.5_dp*r3
             v4o =   1.0_dp/6.0_dp*u3
             v1d =  -8.0_dp + 4.0_dp*p1 - 0.5_dp*p2
             v2d =  10.0_dp - 8.0_dp*q1 + 1.5_dp*q2
             v3d =  -2.0_dp + 4.0_dp*r1 - 1.5_dp*r2
             v4d =   0.5_dp*u2
             DO j =  bo(1,2), bo(2,2)
                my_i= bo(1,1)-gbo(1,1)
                xs1 = REAL(my_i,dp)*dr1c
                rv2 = rt2 - xs2
                vec(2) = rv2
                ivec(2) = FLOOR(vec(2)/pw%pw_grid%dr(2))
                   ij1 = MODULO(ivec(2)-1,npts(2))+1
                   ij2 = MODULO(ivec(2)  ,npts(2))+1
                   ij3 = MODULO(ivec(2)+1,npts(2))+1
                   ij4 = MODULO(ivec(2)+2,npts(2))+1
                xd2  = (vec(2)/dr2)-REAL(ivec(2),kind=dp)
                e1  = 3.0_dp + xd2
                e2  = e1*e1
                e3  = e2*e1
                f1  = 2.0_dp + xd2
                f2  = f1*f1
                f3  = f2*f1
                g1  = 1.0_dp + xd2
                g2  = g1*g1
                g3  = g2*g1
                h1  = xd2
                h2  = h1*h1
                h3  = h2*h1
                s1o =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*e1 + 12.0_dp*e2 - e3)
                s2o = -22.0_dp/3.0_dp + 10.0_dp*f1 - 4.0_dp*f2 + 0.5_dp*f3
                s3o =   2.0_dp/3.0_dp -  2.0_dp*g1 + 2.0_dp*g2 - 0.5_dp*g3
                s4o =   1.0_dp/6.0_dp*h3
                s1d =  -8.0_dp + 4.0_dp*e1 - 0.5_dp*e2
                s2d =  10.0_dp - 8.0_dp*f1 + 1.5_dp*f2
                s3d =  -2.0_dp + 4.0_dp*g1 - 1.5_dp*g2
                s4d =   0.5_dp*h2
                DO i =  bo(1,1), bo(2,1)
                   rv1 = rt1 - xs1
                   vec(1) = rv1
                   ivec(1) = FLOOR(vec(1)/pw%pw_grid%dr(1))
                   ii1 = MODULO(ivec(1)-1,npts(1))+1
                   ii2 = MODULO(ivec(1)  ,npts(1))+1
                   ii3 = MODULO(ivec(1)+1,npts(1))+1
                   ii4 = MODULO(ivec(1)+2,npts(1))+1
                   xd1  = (vec(1)/dr1)-REAL(ivec(1),kind=dp)
                   a1  = 3.0_dp + xd1
                   a2  = a1*a1
                   a3  = a2*a1
                   b1  = 2.0_dp + xd1
                   b2  = b1*b1
                   b3  = b2*b1
                   c1  = 1.0_dp + xd1
                   c2  = c1*c1
                   c3  = c2*c1
                   d1  = xd1
                   d2  = d1*d1
                   d3  = d2*d1
                   t1o =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*a1 + 12.0_dp*a2 - a3)
                   t2o = -22.0_dp/3.0_dp + 10.0_dp*b1 - 4.0_dp*b2 + 0.5_dp*b3
                   t3o =   2.0_dp/3.0_dp -  2.0_dp*c1 + 2.0_dp*c2 - 0.5_dp*c3
                   t4o =   1.0_dp/6.0_dp*d3
                   t1d =  -8.0_dp + 4.0_dp*a1 - 0.5_dp*a2
                   t2d =  10.0_dp - 8.0_dp*b1 + 1.5_dp*b2
                   t3d =  -2.0_dp + 4.0_dp*c1 - 1.5_dp*c2
                   t4d =   0.5_dp*d2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! # v then t then s
!
! for ii in 1 2 3; do
! if [[ $ii == 1 ]]; then ld=t; fi
! if [[ $ii == 2 ]]; then ld=s; fi
! if [[ $ii == 3 ]]; then ld=v; fi
! #
! for l in t s v; do
! for i in 1 2 3 4; do
! if [[ $ld == $l ]]; then
!  echo "$l$i = $l${i}d*dr${ii}i"
! else
!  echo "$l$i = $l${i}o"
! fi
! done
! done
! #
! for i in 1 2 3 4; do
! for j in 1 2 3 4; do
! echo -n "abc_X($i,$j) = "
! for k in 1 2 3 4; do
! if [ $k == 4 ]; then
! echo "grid2(ii$i,ij$j,ik$k)*v$k"
! else
! echo -n "grid2(ii$i,ij$j,ik$k)*v$k + "
! fi
! done
! done
! done
! echo ""
! for j in 1 2 3 4; do
! echo -n "abc_X_Y($j) = "
! for k in 1 2 3 4; do
! if [ $k == 4 ]; then
! echo "abc_X($k,$j)*t$k"
! else
! echo -n "abc_X($k,$j)*t$k + "
! fi
! done
! done
! echo ""
! echo -n "val($ii) = "
! for k in 1 2 3 4; do
! if [ $k == 4 ]; then
! echo "abc_X_Y($k)*s$k"
! else
! echo -n "abc_X_Y($k)*s$k + "
! fi
! done
! echo ""
! done

                   t1     = t1d*dr1i
                   t2     = t2d*dr1i
                   t3     = t3d*dr1i
                   t4     = t4d*dr1i
                   s1     = s1o
                   s2     = s2o
                   s3     = s3o
                   s4     = s4o
                   v1     = v1o
                   v2     = v2o
                   v3     = v3o
                   v4     = v4o
                   abc_X(1,1) = grid2(ii1,ij1,ik1)*v1 + grid2(ii1,ij1,ik2)*v2 + grid2(ii1,ij1,ik3)*v3 + grid2(ii1,ij1,ik4)*v4
                   abc_X(1,2) = grid2(ii1,ij2,ik1)*v1 + grid2(ii1,ij2,ik2)*v2 + grid2(ii1,ij2,ik3)*v3 + grid2(ii1,ij2,ik4)*v4
                   abc_X(1,3) = grid2(ii1,ij3,ik1)*v1 + grid2(ii1,ij3,ik2)*v2 + grid2(ii1,ij3,ik3)*v3 + grid2(ii1,ij3,ik4)*v4
                   abc_X(1,4) = grid2(ii1,ij4,ik1)*v1 + grid2(ii1,ij4,ik2)*v2 + grid2(ii1,ij4,ik3)*v3 + grid2(ii1,ij4,ik4)*v4
                   abc_X(2,1) = grid2(ii2,ij1,ik1)*v1 + grid2(ii2,ij1,ik2)*v2 + grid2(ii2,ij1,ik3)*v3 + grid2(ii2,ij1,ik4)*v4
                   abc_X(2,2) = grid2(ii2,ij2,ik1)*v1 + grid2(ii2,ij2,ik2)*v2 + grid2(ii2,ij2,ik3)*v3 + grid2(ii2,ij2,ik4)*v4
                   abc_X(2,3) = grid2(ii2,ij3,ik1)*v1 + grid2(ii2,ij3,ik2)*v2 + grid2(ii2,ij3,ik3)*v3 + grid2(ii2,ij3,ik4)*v4
                   abc_X(2,4) = grid2(ii2,ij4,ik1)*v1 + grid2(ii2,ij4,ik2)*v2 + grid2(ii2,ij4,ik3)*v3 + grid2(ii2,ij4,ik4)*v4
                   abc_X(3,1) = grid2(ii3,ij1,ik1)*v1 + grid2(ii3,ij1,ik2)*v2 + grid2(ii3,ij1,ik3)*v3 + grid2(ii3,ij1,ik4)*v4
                   abc_X(3,2) = grid2(ii3,ij2,ik1)*v1 + grid2(ii3,ij2,ik2)*v2 + grid2(ii3,ij2,ik3)*v3 + grid2(ii3,ij2,ik4)*v4
                   abc_X(3,3) = grid2(ii3,ij3,ik1)*v1 + grid2(ii3,ij3,ik2)*v2 + grid2(ii3,ij3,ik3)*v3 + grid2(ii3,ij3,ik4)*v4
                   abc_X(3,4) = grid2(ii3,ij4,ik1)*v1 + grid2(ii3,ij4,ik2)*v2 + grid2(ii3,ij4,ik3)*v3 + grid2(ii3,ij4,ik4)*v4
                   abc_X(4,1) = grid2(ii4,ij1,ik1)*v1 + grid2(ii4,ij1,ik2)*v2 + grid2(ii4,ij1,ik3)*v3 + grid2(ii4,ij1,ik4)*v4
                   abc_X(4,2) = grid2(ii4,ij2,ik1)*v1 + grid2(ii4,ij2,ik2)*v2 + grid2(ii4,ij2,ik3)*v3 + grid2(ii4,ij2,ik4)*v4
                   abc_X(4,3) = grid2(ii4,ij3,ik1)*v1 + grid2(ii4,ij3,ik2)*v2 + grid2(ii4,ij3,ik3)*v3 + grid2(ii4,ij3,ik4)*v4
                   abc_X(4,4) = grid2(ii4,ij4,ik1)*v1 + grid2(ii4,ij4,ik2)*v2 + grid2(ii4,ij4,ik3)*v3 + grid2(ii4,ij4,ik4)*v4

                   abc_X_Y(1) = abc_X(1,1)*t1 + abc_X(2,1)*t2 + abc_X(3,1)*t3 + abc_X(4,1)*t4
                   abc_X_Y(2) = abc_X(1,2)*t1 + abc_X(2,2)*t2 + abc_X(3,2)*t3 + abc_X(4,2)*t4
                   abc_X_Y(3) = abc_X(1,3)*t1 + abc_X(2,3)*t2 + abc_X(3,3)*t3 + abc_X(4,3)*t4
                   abc_X_Y(4) = abc_X(1,4)*t1 + abc_X(2,4)*t2 + abc_X(3,4)*t3 + abc_X(4,4)*t4

                   val(1) = abc_X_Y(1)*s1 + abc_X_Y(2)*s2 + abc_X_Y(3)*s3 + abc_X_Y(4)*s4

                   t1     = t1o
                   t2     = t2o
                   t3     = t3o
                   t4     = t4o
                   s1     = s1d*dr2i
                   s2     = s2d*dr2i
                   s3     = s3d*dr2i
                   s4     = s4d*dr2i
                   !! v1 = v1o
                   !! v2 = v2o
                   !! v3 = v3o
                   !! v4 = v4o
                   !! abc_X(1,1) = grid2(ii1,ij1,ik1)*v1 + grid2(ii1,ij1,ik2)*v2 + grid2(ii1,ij1,ik3)*v3 + grid2(ii1,ij1,ik4)*v4
                   !! abc_X(1,2) = grid2(ii1,ij2,ik1)*v1 + grid2(ii1,ij2,ik2)*v2 + grid2(ii1,ij2,ik3)*v3 + grid2(ii1,ij2,ik4)*v4
                   !! abc_X(1,3) = grid2(ii1,ij3,ik1)*v1 + grid2(ii1,ij3,ik2)*v2 + grid2(ii1,ij3,ik3)*v3 + grid2(ii1,ij3,ik4)*v4
                   !! abc_X(1,4) = grid2(ii1,ij4,ik1)*v1 + grid2(ii1,ij4,ik2)*v2 + grid2(ii1,ij4,ik3)*v3 + grid2(ii1,ij4,ik4)*v4
                   !! abc_X(2,1) = grid2(ii2,ij1,ik1)*v1 + grid2(ii2,ij1,ik2)*v2 + grid2(ii2,ij1,ik3)*v3 + grid2(ii2,ij1,ik4)*v4
                   !! abc_X(2,2) = grid2(ii2,ij2,ik1)*v1 + grid2(ii2,ij2,ik2)*v2 + grid2(ii2,ij2,ik3)*v3 + grid2(ii2,ij2,ik4)*v4
                   !! abc_X(2,3) = grid2(ii2,ij3,ik1)*v1 + grid2(ii2,ij3,ik2)*v2 + grid2(ii2,ij3,ik3)*v3 + grid2(ii2,ij3,ik4)*v4
                   !! abc_X(2,4) = grid2(ii2,ij4,ik1)*v1 + grid2(ii2,ij4,ik2)*v2 + grid2(ii2,ij4,ik3)*v3 + grid2(ii2,ij4,ik4)*v4
                   !! abc_X(3,1) = grid2(ii3,ij1,ik1)*v1 + grid2(ii3,ij1,ik2)*v2 + grid2(ii3,ij1,ik3)*v3 + grid2(ii3,ij1,ik4)*v4
                   !! abc_X(3,2) = grid2(ii3,ij2,ik1)*v1 + grid2(ii3,ij2,ik2)*v2 + grid2(ii3,ij2,ik3)*v3 + grid2(ii3,ij2,ik4)*v4
                   !! abc_X(3,3) = grid2(ii3,ij3,ik1)*v1 + grid2(ii3,ij3,ik2)*v2 + grid2(ii3,ij3,ik3)*v3 + grid2(ii3,ij3,ik4)*v4
                   !! abc_X(3,4) = grid2(ii3,ij4,ik1)*v1 + grid2(ii3,ij4,ik2)*v2 + grid2(ii3,ij4,ik3)*v3 + grid2(ii3,ij4,ik4)*v4
                   !! abc_X(4,1) = grid2(ii4,ij1,ik1)*v1 + grid2(ii4,ij1,ik2)*v2 + grid2(ii4,ij1,ik3)*v3 + grid2(ii4,ij1,ik4)*v4
                   !! abc_X(4,2) = grid2(ii4,ij2,ik1)*v1 + grid2(ii4,ij2,ik2)*v2 + grid2(ii4,ij2,ik3)*v3 + grid2(ii4,ij2,ik4)*v4
                   !! abc_X(4,3) = grid2(ii4,ij3,ik1)*v1 + grid2(ii4,ij3,ik2)*v2 + grid2(ii4,ij3,ik3)*v3 + grid2(ii4,ij3,ik4)*v4
                   !! abc_X(4,4) = grid2(ii4,ij4,ik1)*v1 + grid2(ii4,ij4,ik2)*v2 + grid2(ii4,ij4,ik3)*v3 + grid2(ii4,ij4,ik4)*v4

                   abc_X_Y(1) = abc_X(1,1)*t1 + abc_X(2,1)*t2 + abc_X(3,1)*t3 + abc_X(4,1)*t4
                   abc_X_Y(2) = abc_X(1,2)*t1 + abc_X(2,2)*t2 + abc_X(3,2)*t3 + abc_X(4,2)*t4
                   abc_X_Y(3) = abc_X(1,3)*t1 + abc_X(2,3)*t2 + abc_X(3,3)*t3 + abc_X(4,3)*t4
                   abc_X_Y(4) = abc_X(1,4)*t1 + abc_X(2,4)*t2 + abc_X(3,4)*t3 + abc_X(4,4)*t4

                   val(2) = abc_X_Y(1)*s1 + abc_X_Y(2)*s2 + abc_X_Y(3)*s3 + abc_X_Y(4)*s4

                   t1     = t1o
                   t2     = t2o
                   t3     = t3o
                   t4     = t4o
                   s1     = s1o
                   s2     = s2o
                   s3     = s3o
                   s4     = s4o
                   v1     = v1d*dr3i
                   v2     = v2d*dr3i
                   v3     = v3d*dr3i
                   v4     = v4d*dr3i
                   abc_X(1,1) = grid2(ii1,ij1,ik1)*v1 + grid2(ii1,ij1,ik2)*v2 + grid2(ii1,ij1,ik3)*v3 + grid2(ii1,ij1,ik4)*v4
                   abc_X(1,2) = grid2(ii1,ij2,ik1)*v1 + grid2(ii1,ij2,ik2)*v2 + grid2(ii1,ij2,ik3)*v3 + grid2(ii1,ij2,ik4)*v4
                   abc_X(1,3) = grid2(ii1,ij3,ik1)*v1 + grid2(ii1,ij3,ik2)*v2 + grid2(ii1,ij3,ik3)*v3 + grid2(ii1,ij3,ik4)*v4
                   abc_X(1,4) = grid2(ii1,ij4,ik1)*v1 + grid2(ii1,ij4,ik2)*v2 + grid2(ii1,ij4,ik3)*v3 + grid2(ii1,ij4,ik4)*v4
                   abc_X(2,1) = grid2(ii2,ij1,ik1)*v1 + grid2(ii2,ij1,ik2)*v2 + grid2(ii2,ij1,ik3)*v3 + grid2(ii2,ij1,ik4)*v4
                   abc_X(2,2) = grid2(ii2,ij2,ik1)*v1 + grid2(ii2,ij2,ik2)*v2 + grid2(ii2,ij2,ik3)*v3 + grid2(ii2,ij2,ik4)*v4
                   abc_X(2,3) = grid2(ii2,ij3,ik1)*v1 + grid2(ii2,ij3,ik2)*v2 + grid2(ii2,ij3,ik3)*v3 + grid2(ii2,ij3,ik4)*v4
                   abc_X(2,4) = grid2(ii2,ij4,ik1)*v1 + grid2(ii2,ij4,ik2)*v2 + grid2(ii2,ij4,ik3)*v3 + grid2(ii2,ij4,ik4)*v4
                   abc_X(3,1) = grid2(ii3,ij1,ik1)*v1 + grid2(ii3,ij1,ik2)*v2 + grid2(ii3,ij1,ik3)*v3 + grid2(ii3,ij1,ik4)*v4
                   abc_X(3,2) = grid2(ii3,ij2,ik1)*v1 + grid2(ii3,ij2,ik2)*v2 + grid2(ii3,ij2,ik3)*v3 + grid2(ii3,ij2,ik4)*v4
                   abc_X(3,3) = grid2(ii3,ij3,ik1)*v1 + grid2(ii3,ij3,ik2)*v2 + grid2(ii3,ij3,ik3)*v3 + grid2(ii3,ij3,ik4)*v4
                   abc_X(3,4) = grid2(ii3,ij4,ik1)*v1 + grid2(ii3,ij4,ik2)*v2 + grid2(ii3,ij4,ik3)*v3 + grid2(ii3,ij4,ik4)*v4
                   abc_X(4,1) = grid2(ii4,ij1,ik1)*v1 + grid2(ii4,ij1,ik2)*v2 + grid2(ii4,ij1,ik3)*v3 + grid2(ii4,ij1,ik4)*v4
                   abc_X(4,2) = grid2(ii4,ij2,ik1)*v1 + grid2(ii4,ij2,ik2)*v2 + grid2(ii4,ij2,ik3)*v3 + grid2(ii4,ij2,ik4)*v4
                   abc_X(4,3) = grid2(ii4,ij3,ik1)*v1 + grid2(ii4,ij3,ik2)*v2 + grid2(ii4,ij3,ik3)*v3 + grid2(ii4,ij3,ik4)*v4
                   abc_X(4,4) = grid2(ii4,ij4,ik1)*v1 + grid2(ii4,ij4,ik2)*v2 + grid2(ii4,ij4,ik3)*v3 + grid2(ii4,ij4,ik4)*v4

                   abc_X_Y(1) = abc_X(1,1)*t1 + abc_X(2,1)*t2 + abc_X(3,1)*t3 + abc_X(4,1)*t4
                   abc_X_Y(2) = abc_X(1,2)*t1 + abc_X(2,2)*t2 + abc_X(3,2)*t3 + abc_X(4,2)*t4
                   abc_X_Y(3) = abc_X(1,3)*t1 + abc_X(2,3)*t2 + abc_X(3,3)*t3 + abc_X(4,3)*t4
                   abc_X_Y(4) = abc_X(1,4)*t1 + abc_X(2,4)*t2 + abc_X(3,4)*t3 + abc_X(4,4)*t4

                   val(3) = abc_X_Y(1)*s1 + abc_X_Y(2)*s2 + abc_X_Y(3)*s3 + abc_X_Y(4)*s4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                   fac = grid(i,j,k)
                   ft1 = ft1 + val(1) * fac
                   ft2 = ft2 + val(2) * fac
                   ft3 = ft3 + val(3) * fac
                   xs1 = xs1 + dr1c
                END DO
                xs2 = xs2 + dr2c
             END DO
          END DO LoopOnGrid
          qt = - qt * dvol
          LForces(1,LindMM) = ft1 * qt
          LForces(2,LindMM) = ft2 * qt
          LForces(3,LindMM) = ft3 * qt

          Forces(1,LIndMM) = Forces(1,LIndMM) + LForces(1,LindMM)
          Forces(2,LIndMM) = Forces(2,LIndMM) + LForces(2,LindMM)
          Forces(3,LIndMM) = Forces(3,LIndMM) + LForces(3,LindMM)
       END DO Atoms
    END DO Radius
    !
    ! Debug Statement
    !
    IF (debug_this_module) THEN
       CALL debug_qmmm_forces_with_gauss_LG(pgfs=pgfs,&
                                            aug_pools=aug_pools,&
                                            rho=cgrid,&
                                            num_mm_atoms=num_mm_atoms,&
                                            mm_charges=mm_charges,&
                                            mm_atom_index=mm_atom_index,&
                                            mm_particles=mm_particles,&
                                            mm_el_pot_radius=mm_el_pot_radius,&
                                            coarser_grid_level=coarser_grid_level,&
                                            debug_force=LForces,&
                                            per_potentials=per_potentials,&
                                            para_env=para_env,&
                                            mm_cell=mm_cell,&
                                            dOmmOqm=dOmmOqm,&
                                            iw=iw,&
                                            par_scheme=par_scheme,&
                                            qmmm_spherical_cutoff=qmmm_spherical_cutoff,&
                                            error=error)
    END IF
    DEALLOCATE(LForces, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE qmmm_forces_with_gaussian_LG

! *****************************************************************************
!> \brief Evaluates the contribution to the forces due to the Long Range
!>      part of the QM/MM potential computed collocating the Electrostatic
!>      Gaussian Potential.
!> \param pgfs ...
!> \param cgrid ...
!> \param num_mm_atoms ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_particles ...
!> \param mm_el_pot_radius ...
!> \param para_env ...
!> \param coarser_grid_level ...
!> \param Forces ...
!> \param potentials ...
!> \param aug_pools ...
!> \param mm_cell ...
!> \param dOmmOqm ...
!> \param iw ...
!> \param par_scheme ...
!> \param qmmm_spherical_cutoff ...
!> \param error ...
!> \par History
!>      08.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_forces_with_gaussian_LR  (pgfs,cgrid,num_mm_atoms,mm_charges,mm_atom_index,&
       mm_particles,mm_el_pot_radius,para_env, coarser_grid_level,Forces, potentials,&
       aug_pools, mm_cell, dOmmOqm, iw, par_scheme, qmmm_spherical_cutoff,error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_type), POINTER                   :: cgrid
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: coarser_grid_level
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: Potentials
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    INTEGER, INTENT(IN)                      :: iw, par_scheme
    REAL(KIND=dp), DIMENSION(2), INTENT(IN)  :: qmmm_spherical_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_with_gaussian_LR', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, Imm, IndMM, &
                                                IRadTyp, ix, j, k, LIndMM, &
                                                my_i, my_j, my_k, myind, n1, &
                                                n2, n3, stat
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    LOGICAL                                  :: failure
    REAL(KIND=dp) :: dr1, dr2, dr3, dvol, dx, fac, ft1, ft2, ft3, qt, r, r2, &
      rd1, rd2, rd3, rt1, rt2, rt3, rv1, rv2, rv3, rx, rx2, sph_chrg_factor, &
      Term, xs1, xs2, xs3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: LForces
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pot0_2
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid
    TYPE(qmmm_pot_type), POINTER             :: pot

    CALL timeset(routineN,handle)
    failure = .FALSE.
    ALLOCATE(LForces(3,num_mm_atoms), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    LForces = 0.0_dp
    n1   = cgrid%pw_grid%npts(1)
    n2   = cgrid%pw_grid%npts(2)
    n3   = cgrid%pw_grid%npts(3)
    dr1  = cgrid%pw_grid%dr(1)
    dr2  = cgrid%pw_grid%dr(2)
    dr3  = cgrid%pw_grid%dr(3)
    dvol = cgrid%pw_grid%dvol
    gbo  = cgrid%pw_grid%bounds
    bo    =  cgrid%pw_grid%bounds_local
    grid  => cgrid%cr3d
    IF (par_scheme==do_par_atom) myind = 0
    Radius: DO IRadTyp = 1, SIZE(pgfs)
       pot => potentials(IRadTyp)%pot
       dx     =  Pot%dx
       pot0_2 => Pot%pot0_2
       Atoms: DO Imm = 1, SIZE(pot%mm_atom_index)
          IF (par_scheme==do_par_atom) THEN
             myind = myind + 1
             IF (MOD(myind,para_env%num_pe)/=para_env%mepos) CYCLE
          END IF
          LIndMM    =   pot%mm_atom_index(Imm)
          IndMM     =   mm_atom_index(LIndMM)
          ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm,mm_cell)+dOmmOqm
          qt        =   mm_charges(LIndMM)
          ! Possible Spherical Cutoff
          IF (qmmm_spherical_cutoff(1)>0.0_dp) THEN
             CALL spherical_cutoff_factor(qmmm_spherical_cutoff, ra, sph_chrg_factor, error)
             qt = qt * sph_chrg_factor
          END IF
          IF (ABS(qt)<= EPSILON(0.0_dp)) CYCLE Atoms
          rt1 = ra(1)
          rt2 = ra(2)
          rt3 = ra(3)
          ft1 = 0.0_dp
          ft2 = 0.0_dp
          ft3 = 0.0_dp
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k = k-gbo(1,3)
             xs3  = REAL(my_k,dp)*dr3
             my_j = bo(1,2)-gbo(1,2)
             xs2 = REAL(my_j,dp)*dr2
             rv3 = rt3 - xs3
             DO j =  bo(1,2), bo(2,2)
                my_i= bo(1,1)-gbo(1,1)
                xs1 = REAL(my_i,dp)*dr1
                rv2 = rt2 - xs2
                DO i =  bo(1,1), bo(2,1)
                   rv1 = rt1 - xs1
                   r2  = rv1*rv1 + rv2*rv2 + rv3*rv3
                   r   = SQRT(r2)
                   ix  = FLOOR(r/dx)+1
                   rx  = (r-REAL(ix-1,dp)*dx)/dx
                   rx2 = rx*rx
                   Term = pot0_2(1,ix  )*(-6.0_dp*(rx-rx2))                &
                         +pot0_2(2,ix  )*(1.0_dp-4.0_dp*rx+3.0_dp*rx2)     &
                         +pot0_2(1,ix+1)*(6.0_dp*(rx-rx2))                 &
                         +pot0_2(2,ix+1)*(-2.0_dp*rx+3.0_dp*rx2)
                   fac = grid(i,j,k) * Term
                   IF ( r == 0.0_dp ) THEN
                      rd1 = 1.0_dp
                      rd2 = 1.0_dp
                      rd3 = 1.0_dp
                   ELSE
                      rd1 = rv1 / r
                      rd2 = rv2 / r
                      rd3 = rv3 / r
                   END IF
                   ft1 = ft1 + fac * rd1
                   ft2 = ft2 + fac * rd2
                   ft3 = ft3 + fac * rd3
                   xs1 = xs1 + dr1
                END DO
                xs2 = xs2 + dr2
             END DO
          END DO LoopOnGrid
          qt = - qt * dvol / dx
          LForces(1,LindMM) = ft1 * qt
          LForces(2,LindMM) = ft2 * qt
          LForces(3,LindMM) = ft3 * qt

          Forces(1,LIndMM) = Forces(1,LIndMM) + LForces(1,LindMM)
          Forces(2,LIndMM) = Forces(2,LIndMM) + LForces(2,LindMM)
          Forces(3,LIndMM) = Forces(3,LIndMM) + LForces(3,LindMM)
       END DO Atoms
    END DO Radius
    !
    ! Debug Statement
    !
    IF (debug_this_module) THEN
       CALL debug_qmmm_forces_with_gauss_LR(pgfs=pgfs,&
                                            aug_pools=aug_pools,&
                                            rho=cgrid,&
                                            num_mm_atoms=num_mm_atoms,&
                                            mm_charges=mm_charges,&
                                            mm_atom_index=mm_atom_index,&
                                            mm_particles=mm_particles,&
                                            mm_el_pot_radius=mm_el_pot_radius,&
                                            coarser_grid_level=coarser_grid_level,&
                                            debug_force=LForces,&
                                            potentials=potentials,&
                                            para_env=para_env,&
                                            mm_cell=mm_cell,&
                                            dOmmOqm=dOmmOqm,&
                                            iw=iw,&
                                            par_scheme=par_scheme,&
                                            qmmm_spherical_cutoff=qmmm_spherical_cutoff,&
                                            error=error)
    END IF

    DEALLOCATE(LForces, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
  END SUBROUTINE qmmm_forces_with_gaussian_LR

! *****************************************************************************
!> \brief Evaluates numerically QM/MM forces and compares them with
!>      the analytically computed ones.
!>      It is evaluated only when debug_this_module is set to .TRUE.
!> \param rho ...
!> \param qs_env ...
!> \param qmmm_env ...
!> \param Analytical_Forces ...
!> \param mm_particles ...
!> \param mm_atom_index ...
!> \param num_mm_atoms ...
!> \param interp_section ...
!> \param mm_cell ...
!> \param error ...
!> \par History
!>      08.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_debug_forces(rho,qs_env,qmmm_env,Analytical_Forces,&
                               mm_particles,mm_atom_index,num_mm_atoms,&
                               interp_section,mm_cell,error)
    TYPE(pw_type), POINTER                   :: rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Analytical_Forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(section_vals_type), POINTER         :: interp_section
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_debug_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, I, IndMM, iw, J, K, &
                                                stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: Coord_save
    REAL(KIND=dp), DIMENSION(2)              :: energy
    REAL(KIND=dp), DIMENSION(3)              :: Err
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Num_Forces
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: v_qmmm_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(section_vals_type), POINTER         :: input_section, print_section

    CALL timeset(routineN,handle)
    NULLIFY( Num_Forces )
    CALL get_qs_env(qs_env=qs_env,&
                    pw_env=pw_env,&
                    input=input_section,&
                    para_env=para_env,&
                    error=error)

    print_section => section_vals_get_subs_vals(input_section,"QMMM%PRINT",error=error)
    logger => cp_error_get_logger(error)
    iw=cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",extension=".qmmmLog",error=error)
    CALL pw_env_get(pw_env=pw_env, pw_pools=pw_pools, error=error)
    CALL pw_pool_create_pw(pw_pools(1)%pool, v_qmmm_rspace%pw,&
         use_data=REALDATA3D, in_space=REALSPACE,error=error)
    ALLOCATE(Num_Forces(3,num_mm_atoms),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ks_qmmm_env_loc =>  qs_env%ks_qmmm_env
    IF (iw>0) WRITE(iw,'(/A)')"DEBUG SECTION:"
    Atoms: DO I = 1, num_mm_atoms
       IndMM = mm_atom_index(I)
       Coords: DO J = 1, 3
          Coord_save = mm_particles(IndMM)%r(J)
          energy = 0.0_dp
          Diff: DO K = 1, 2
             mm_particles(IndMM)%r(J) = Coord_save + (-1)**K * Dx
             CALL pw_zero(v_qmmm_rspace%pw,error=error)
             SELECT CASE(qmmm_env%qmmm_coupl_type)
             CASE(do_qmmm_coulomb)
                CALL cp_unimplemented_error(fromWhere=routineP, &
                     message="Coulomb QM/MM electrostatic coupling not implemented for GPW/GAPW.", &
                     error=error, error_level=cp_failure_level)
             CASE(do_qmmm_pcharge)
                CALL cp_unimplemented_error(fromWhere=routineP, &
                     message="Point Charge  QM/MM electrostatic coupling not implemented for GPW/GAPW.", &
                     error=error, error_level=cp_failure_level)
             CASE(do_qmmm_gauss,do_qmmm_swave)
                CALL    qmmm_elec_with_gaussian(qmmm_env=qmmm_env,&
                                                v_qmmm=v_qmmm_rspace,&
                                                mm_particles=mm_particles,&
                                                aug_pools=qmmm_env%aug_pools,&
                                                para_env=para_env,&
                                                eps_mm_rspace=qmmm_env%eps_mm_rspace,&
                                                cube_info=ks_qmmm_env_loc%cube_info,&
                                                pw_pools=pw_pools,&
                                                auxbas_grid=qmmm_env%gridlevel_info%auxbas_grid,&
                                                coarser_grid=qmmm_env%gridlevel_info%coarser_grid,&
                                                interp_section=interp_section,&
                                                mm_cell=mm_cell,&
                                                error=error)
             CASE(do_qmmm_none)
                CYCLE Diff
             CASE DEFAULT
                IF (iw>0) WRITE(iw,'(T3,A)')"Unknown Coupling..."
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             energy(K) =  pw_integral_ab ( rho, v_qmmm_rspace%pw ,error=error)
          END DO Diff
          IF (iw>0) THEN
             WRITE(iw,'(A,I6,A,I3,A,2F15.9)')&
                  "DEBUG :: MM Atom = ",IndMM," Coord = ",J," Energies (+/-) :: ",energy(2), energy(1)
          END IF
          Num_Forces(J,I)  = ( energy(2) - energy(1) ) / (2.0_dp * Dx)
          mm_particles(IndMM)%r(J) = Coord_save
       END DO Coords
    END DO Atoms

    SELECT CASE(qmmm_env%qmmm_coupl_type)
    CASE(do_qmmm_coulomb)
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Coulomb QM/MM electrostatic coupling not implemented for GPW/GAPW.", &
            error=error, error_level=cp_failure_level)
    CASE(do_qmmm_pcharge)
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Point Charge QM/MM electrostatic coupling not implemented for GPW/GAPW.", &
            error=error, error_level=cp_failure_level)
    CASE(do_qmmm_gauss,do_qmmm_swave)
       IF (iw>0) WRITE(iw,'(/A/)')"CHECKING NUMERICAL Vs ANALYTICAL FORCES (Err%):"
       DO I =  1, num_mm_atoms
          IndMM = mm_atom_index(I)
          Err = 0.0_dp
          DO K = 1, 3
             IF (ABS(Num_Forces(K,I)) >= 5.0E-5_dp) THEN
                Err(K) = (Analytical_Forces(K,I)-Num_Forces(K,I))/Num_Forces(K,I)*100.0_dp
             END IF
          END DO
          IF (iw>0)&
               WRITE(iw,100)IndMM,Analytical_Forces(1,I),Num_Forces(1,I),Err(1),&
                                  Analytical_Forces(2,I),Num_Forces(2,I),Err(2),&
                                  Analytical_Forces(3,I),Num_Forces(3,I),Err(3)
          CPPostcondition(ABS(Err(1))<=MaxErr,cp_failure_level,routineP,error,failure)
          CPPostcondition(ABS(Err(2))<=MaxErr,cp_failure_level,routineP,error,failure)
          CPPostcondition(ABS(Err(3))<=MaxErr,cp_failure_level,routineP,error,failure)
       END DO
    CASE(do_qmmm_none)
       IF (iw>0) WRITE(iw,'(T3,A)')"No QM/MM Derivatives to debug. Just Mechanical Coupling!"
    CASE DEFAULT
       IF (iw>0) WRITE(iw,'(T3,A)')"Unknown Coupling..."
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    CALL cp_print_key_finished_output(iw,logger,print_section,"PROGRAM_RUN_INFO",error=error)

    CALL pw_pool_give_back_pw ( pw_pools(1)%pool, v_qmmm_rspace%pw, error=error)
    DEALLOCATE(Num_Forces,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(handle)
100 FORMAT(I5,2F15.9," ( ",F7.2," ) ",2F15.9," ( ",F7.2," ) ",2F15.9," ( ",F7.2," ) ")
  END SUBROUTINE qmmm_debug_forces

! *****************************************************************************
!> \brief Debugs the integrate_gf_rspace_NoPBC.. It may helps ;-P
!> \param ilevel ...
!> \param zetp ...
!> \param rp ...
!> \param W ...
!> \param pwgrid ...
!> \param cube_info ...
!> \param eps_mm_rspace ...
!> \param aug_pools ...
!> \param debug_force ...
!> \param mm_cell ...
!> \param auxbas_grid ...
!> \param n_rep_real ...
!> \param iw ...
!> \param error ...
!> \par History
!>      08.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE debug_integrate_gf_rspace_NoPBC(ilevel, zetp, rp, W, pwgrid, cube_info,&
                                             eps_mm_rspace, aug_pools, debug_force,&
                                             mm_cell,auxbas_grid, n_rep_real, iw, error)
    INTEGER, INTENT(IN)                      :: ilevel
    REAL(KIND=dp), INTENT(IN)                :: zetp
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rp
    REAL(KIND=dp), INTENT(IN)                :: W
    TYPE(pw_type), POINTER                   :: pwgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: debug_force
    TYPE(cell_type), POINTER                 :: mm_cell
    INTEGER, INTENT(IN)                      :: auxbas_grid
    INTEGER, DIMENSION(3), INTENT(IN)        :: n_rep_real
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'debug_integrate_gf_rspace_NoPBC', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, igrid, k, ngrids, &
                                                stat
    INTEGER, DIMENSION(2, 3)                 :: bo2
    INTEGER, SAVE                            :: Icount
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(2)              :: energy
    REAL(KIND=dp), DIMENSION(3)              :: Err, force, myrp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdat, ydat, zdat
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids

    DATA Icount /0/
    ! Statements
    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(grids)
    !Statements
    ngrids = SIZE(aug_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw,error=error)
    END DO
    bo2 = grids(auxbas_grid)%pw%pw_grid%bounds
    ALLOCATE (xdat(bo2(1,1):bo2(2,1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (ydat(bo2(1,2):bo2(2,2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zdat(bo2(1,3):bo2(2,3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    Icount = Icount + 1
    DO i = 1, 3
       DO k = 1, 2
          myrp = rp
          myrp(i) = myrp(i) + (-1.0_dp)**k * Dx
          CALL pw_zero(grids(ilevel)%pw,error=error)
          CALL collocate_gf_rspace_NoPBC(zetp=zetp,&
                                         rp=myrp,&
                                         scale=-1.0_dp,&
                                         W=W,&
                                         pwgrid=grids(ilevel)%pw,&
                                         cube_info=cube_info,&
                                         eps_mm_rspace=eps_mm_rspace,&
                                         xdat=xdat,&
                                         ydat=ydat,&
                                         zdat=zdat,&
                                         bo2=bo2,&
                                         n_rep_real=n_rep_real,&
                                         mm_cell=mm_cell)

          energy(k) = pw_integral_ab(pwgrid, grids(ilevel)%pw,error=error)
       END DO
       force(i) = ( energy(2) - energy(1) ) / (2.0_dp * Dx)
    END DO
    Err = 0.0_dp
    IF (ALL(force /= 0.0_dp)) THEN
       Err(1) = (debug_force(1)-force(1))/force(1)*100.0_dp
       Err(2) = (debug_force(2)-force(2))/force(2)*100.0_dp
       Err(3) = (debug_force(3)-force(3))/force(3)*100.0_dp
    END IF
    IF (iw>0) &
         WRITE(iw,100)Icount, debug_force(1), force(1), Err(1),&
                              debug_force(2), force(2), Err(2),&
                              debug_force(3), force(3), Err(3)
    CPPostcondition(ABS(Err(1))<=MaxErr,cp_failure_level,routineP,error,failure)
    CPPostcondition(ABS(Err(2))<=MaxErr,cp_failure_level,routineP,error,failure)
    CPPostcondition(ABS(Err(3))<=MaxErr,cp_failure_level,routineP,error,failure)

    IF (ASSOCIATED(xdat)) THEN
       DEALLOCATE (xdat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(ydat)) THEN
       DEALLOCATE (ydat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(zdat)) THEN
       DEALLOCATE (zdat, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)
    CALL timestop(handle)
100 FORMAT("Collocation   : ",I5,2F15.9," ( ",F7.2," ) ",2F15.9," ( ",F7.2," ) ",2F15.9," ( ",F7.2," ) ")
  END SUBROUTINE debug_integrate_gf_rspace_NoPBC

! *****************************************************************************
!> \brief Debugs qmmm_forces_with_gaussian_LG.. It may helps too ... ;-]
!> \param pgfs ...
!> \param aug_pools ...
!> \param rho ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_particles ...
!> \param mm_el_pot_radius ...
!> \param num_mm_atoms ...
!> \param coarser_grid_level ...
!> \param per_potentials ...
!> \param debug_force ...
!> \param para_env ...
!> \param mm_cell ...
!> \param dOmmOqm ...
!> \param iw ...
!> \param par_scheme ...
!> \param qmmm_spherical_cutoff ...
!> \param error ...
!> \par History
!>      08.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE debug_qmmm_forces_with_gauss_LG(pgfs,aug_pools, rho, mm_charges, mm_atom_index,&
       mm_particles, mm_el_pot_radius, num_mm_atoms, coarser_grid_level, per_potentials,&
       debug_force, para_env, mm_cell,dOmmOqm,iw,par_scheme,qmmm_spherical_cutoff,error)

    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(pw_type), POINTER                   :: rho
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius
    INTEGER, INTENT(IN)                      :: num_mm_atoms, &
                                                coarser_grid_level
    TYPE(qmmm_per_pot_p_type), &
      DIMENSION(:), POINTER                  :: per_potentials
    REAL(KIND=dp), DIMENSION(:, :)           :: debug_force
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    INTEGER, INTENT(IN)                      :: iw, par_scheme
    REAL(KIND=dp), DIMENSION(2), INTENT(IN)  :: qmmm_spherical_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'debug_qmmm_forces_with_gauss_LG', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, I, igrid, IndMM, J, &
                                                K, ngrids, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: Coord_save
    REAL(KIND=dp), DIMENSION(2)              :: energy
    REAL(KIND=dp), DIMENSION(3)              :: Err
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Num_Forces
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids

    ALLOCATE(Num_Forces(3,num_mm_atoms),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(grids)
    CALL timeset(routineN,handle)
    ngrids = SIZE(aug_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw,error=error)
    END DO
    Atoms: DO I = 1, num_mm_atoms
       IndMM = mm_atom_index(I)
       Coords: DO J = 1, 3
          Coord_save = mm_particles(IndMM)%r(J)
          energy = 0.0_dp
          Diff: DO K = 1, 2
             mm_particles(IndMM)%r(J) = Coord_save + (-1)**K * Dx
             CALL pw_zero(grids(coarser_grid_level)%pw,error=error)

             CALL qmmm_elec_with_gaussian_LG    (pgfs=pgfs,&
                                                 cgrid=grids(coarser_grid_level)%pw,&
                                                 mm_charges=mm_charges,&
                                                 mm_atom_index=mm_atom_index,&
                                                 mm_particles=mm_particles,&
                                                 para_env=para_env,&
                                                 coarser_grid_level=coarser_grid_level,&
                                                 per_potentials=per_potentials,&
                                                 mm_cell=mm_cell,&
                                                 dOmmOqm=dOmmOqm,&
                                                 par_scheme=par_scheme,&
                                                 qmmm_spherical_cutoff=qmmm_spherical_cutoff,&
                                                 error=error)

             energy(K) =  pw_integral_ab ( rho, grids(coarser_grid_level)%pw ,error=error)
          END DO Diff
          IF (iw>0)&
               WRITE(iw,'(A,I6,A,I3,A,2F15.9)')&
               "DEBUG LR:: MM Atom = ",IndMM," Coord = ",J," Energies (+/-) :: ",energy(2), energy(1)
          Num_Forces(J,I)  = ( energy(2) - energy(1) ) / (2.0_dp * Dx)
          mm_particles(IndMM)%r(J) = Coord_save
       END DO Coords
    END DO Atoms

    DO I =  1, num_mm_atoms
       IndMM = mm_atom_index(I)
       Err = 0.0_dp
       IF (ALL(Num_Forces /= 0.0_dp)) THEN
          Err(1) = (debug_force(1,I)-Num_Forces(1,I))/Num_Forces(1,I)*100.0_dp
          Err(2) = (debug_force(2,I)-Num_Forces(2,I))/Num_Forces(2,I)*100.0_dp
          Err(3) = (debug_force(3,I)-Num_Forces(3,I))/Num_Forces(3,I)*100.0_dp
       END IF
       IF (iw>0) &
            WRITE(iw,100)IndMM,debug_force(1,I),Num_Forces(1,I),Err(1),&
                               debug_force(2,I),Num_Forces(2,I),Err(2),&
                               debug_force(3,I),Num_Forces(3,I),Err(3)
       CPPostcondition(ABS(Err(1))<=MaxErr,cp_failure_level,routineP,error,failure)
       CPPostcondition(ABS(Err(2))<=MaxErr,cp_failure_level,routineP,error,failure)
       CPPostcondition(ABS(Err(3))<=MaxErr,cp_failure_level,routineP,error,failure)
    END DO

    DEALLOCATE(Num_Forces,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)
    CALL timestop(handle)
100 FORMAT("MM Atom LR    : ",I5,2F15.9," ( ",F7.2," ) ",2F15.9," ( ",F7.2," ) ",2F15.9," ( ",F7.2," ) ")
  END SUBROUTINE debug_qmmm_forces_with_gauss_LG

! *****************************************************************************
!> \brief Debugs qmmm_forces_with_gaussian_LR.. It may helps too ... ;-]
!> \param pgfs ...
!> \param aug_pools ...
!> \param rho ...
!> \param mm_charges ...
!> \param mm_atom_index ...
!> \param mm_particles ...
!> \param mm_el_pot_radius ...
!> \param num_mm_atoms ...
!> \param coarser_grid_level ...
!> \param potentials ...
!> \param debug_force ...
!> \param para_env ...
!> \param mm_cell ...
!> \param dOmmOqm ...
!> \param iw ...
!> \param par_scheme ...
!> \param qmmm_spherical_cutoff ...
!> \param error ...
!> \par History
!>      08.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE debug_qmmm_forces_with_gauss_LR(pgfs,aug_pools, rho, mm_charges, mm_atom_index,&
       mm_particles, mm_el_pot_radius, num_mm_atoms, coarser_grid_level, potentials,&
       debug_force, para_env, mm_cell,dOmmOqm,iw, par_scheme, qmmm_spherical_cutoff, error)

    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(pw_type), POINTER                   :: rho
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius
    INTEGER, INTENT(IN)                      :: num_mm_atoms, &
                                                coarser_grid_level
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: Potentials
    REAL(KIND=dp), DIMENSION(:, :)           :: debug_force
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cell_type), POINTER                 :: mm_cell
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    INTEGER, INTENT(IN)                      :: iw, par_scheme
    REAL(KIND=dp), DIMENSION(2), INTENT(IN)  :: qmmm_spherical_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'debug_qmmm_forces_with_gauss_LR', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, I, igrid, IndMM, J, &
                                                K, ngrids, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: Coord_save
    REAL(KIND=dp), DIMENSION(2)              :: energy
    REAL(KIND=dp), DIMENSION(3)              :: Err
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Num_Forces
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids

    ALLOCATE(Num_Forces(3,num_mm_atoms),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(grids)
    CALL timeset(routineN,handle)
    ngrids = SIZE(aug_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw,error=error)
    END DO
    Atoms: DO I = 1, num_mm_atoms
       IndMM = mm_atom_index(I)
       Coords: DO J = 1, 3
          Coord_save = mm_particles(IndMM)%r(J)
          energy = 0.0_dp
          Diff: DO K = 1, 2
             mm_particles(IndMM)%r(J) = Coord_save + (-1)**K * Dx
             CALL pw_zero(grids(coarser_grid_level)%pw,error=error)

             CALL qmmm_elec_with_gaussian_LR    (pgfs=pgfs,&
                                                 grid=grids(coarser_grid_level)%pw,&
                                                 mm_charges=mm_charges,&
                                                 mm_atom_index=mm_atom_index,&
                                                 mm_particles=mm_particles,&
                                                 mm_el_pot_radius=mm_el_pot_radius,&
                                                 para_env=para_env,&
                                                 coarser_grid_level=coarser_grid_level,&
                                                 potentials=potentials,&
                                                 mm_cell=mm_cell,&
                                                 dOmmOqm=dOmmOqm,&
                                                 par_scheme=par_scheme,&
                                                 qmmm_spherical_cutoff=qmmm_spherical_cutoff,&
                                                 error=error)

             energy(K) =  pw_integral_ab ( rho, grids(coarser_grid_level)%pw ,error=error)
          END DO Diff
          IF (iw>0)&
               WRITE(iw,'(A,I6,A,I3,A,2F15.9)')&
               "DEBUG LR:: MM Atom = ",IndMM," Coord = ",J," Energies (+/-) :: ",energy(2), energy(1)
          Num_Forces(J,I)  = ( energy(2) - energy(1) ) / (2.0_dp * Dx)
          mm_particles(IndMM)%r(J) = Coord_save
       END DO Coords
    END DO Atoms

    DO I =  1, num_mm_atoms
       IndMM = mm_atom_index(I)
       Err = 0.0_dp
       IF (ALL(Num_Forces(:,I) /= 0.0_dp)) THEN
          Err(1) = (debug_force(1,I)-Num_Forces(1,I))/Num_Forces(1,I)*100.0_dp
          Err(2) = (debug_force(2,I)-Num_Forces(2,I))/Num_Forces(2,I)*100.0_dp
          Err(3) = (debug_force(3,I)-Num_Forces(3,I))/Num_Forces(3,I)*100.0_dp
       END IF
       IF (iw>0)&
            WRITE(iw,100)IndMM,debug_force(1,I),Num_Forces(1,I),Err(1),&
                               debug_force(2,I),Num_Forces(2,I),Err(2),&
                               debug_force(3,I),Num_Forces(3,I),Err(3)
       CPPostcondition(ABS(Err(1))<=MaxErr,cp_failure_level,routineP,error,failure)
       CPPostcondition(ABS(Err(2))<=MaxErr,cp_failure_level,routineP,error,failure)
       CPPostcondition(ABS(Err(3))<=MaxErr,cp_failure_level,routineP,error,failure)
    END DO

    DEALLOCATE(Num_Forces,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)
    CALL timestop(handle)
100 FORMAT("MM Atom LR    : ",I5,2F15.9," ( ",F7.2," ) ",2F15.9," ( ",F7.2," ) ",2F15.9," ( ",F7.2," ) ")
  END SUBROUTINE debug_qmmm_forces_with_gauss_LR

END MODULE qmmm_gpw_forces
