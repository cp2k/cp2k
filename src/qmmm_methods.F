
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002,2003,2004  CP2K developers group                       !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_methods [1.0] *
!!
!!   NAME
!!     qmmm_methods
!!
!!   FUNCTION
!!     Routines to compute energy and forces in a QM/MM calculation
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!
!!   SOURCE
!****************************************************************************

MODULE qmmm_methods
  USE coefficient_types,               ONLY: coeff_type,&
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE global_types,                    ONLY: global_environment_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: cube_info_type
  USE input_constants,                 ONLY: do_qmmm_gauss,&
                                             do_qmmm_grid,&
                                             do_qmmm_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mm_collocate_potential,          ONLY: collocate_gf_rspace_NoPBC,&
                                             integrate_gf_rspace_NoPBC
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws
  USE pw_spline_utils,                 ONLY: pw_inject_s3,&
                                             pw_project_s3
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_inject_NoPBC,&
                                             pw_p_type,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type,&
                                             pw_zero
  USE pws,                             ONLY: pw_integral_ab
  USE qmmm_electrostatic_methods,      ONLY: qmmm_elec_with_gaussian,&
                                             qmmm_elec_with_gaussian_LR,&
                                             qmmm_elec_with_gaussian_low,&
                                             qmmm_electrostatic_on_grid
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_pot_p_type,&
                                             qmmm_pot_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc

  IMPLICIT NONE
  
  PRIVATE
  LOGICAL,        PARAMETER, PRIVATE   :: debug_this_module=.FALSE.
  REAL(KIND=dp), PARAMETER, PRIVATE   :: Dx = 0.01_dp    ! Debug Variables
  REAL(KIND=dp), PARAMETER, PRIVATE   :: MaxErr = 1.0_dp ! Debug Variables
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_methods'
  PUBLIC :: qmmm_calculate_energy,&
            qmmm_modify_hartree_pot,&
            qmmm_forces
  
CONTAINS
!!****f* qmmm_methods/qmmm_calculate_energy [1.0] *
!!
!!   NAME
!!     qmmm_calculate_energy
!!
!!   FUNCTION
!!     Computes the contribution to the total energy of the QM/MM 
!!     electrostatic coupling
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_calculate_energy(qs_env, rho, v_qmmm, qmmm_energy, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(coeff_type), DIMENSION(:), POINTER  :: rho
    TYPE(coeff_type), INTENT(IN)             :: v_qmmm
    REAL(KIND=dp), INTENT(INOUT)             :: qmmm_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_calculate_energy', &
      routineP = moduleN//':'//routineN   

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,"I","",handle)    
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    NULLIFY(dft_control)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    error=error)

    qmmm_energy = 0.0_dp
    DO ispin=1, dft_control%nspins
       qmmm_energy = qmmm_energy + pw_integral_ab ( rho(ispin)%pw, v_qmmm%pw )
    END DO

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_calculate_energy

!!****f* qmmm_methods/qmmm_modify_hartree_pot [1.0] *
!!
!!   NAME
!!     qmmm_modify_hartree_pot
!!
!!   FUNCTION
!!     Modify the hartree potential in order to include the QM/MM correction
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_modify_hartree_pot(v_hartree,v_qmmm, error)
    TYPE(coeff_type), INTENT(INOUT)          :: v_hartree
    TYPE(coeff_type), INTENT(IN)             :: v_qmmm
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_modify_hartree_pot', &
      routineP = moduleN//':'//routineN   

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,"I","",handle)    
    failure=.FALSE.

    v_hartree%pw%cr3d = v_hartree%pw%cr3d + &
                        v_qmmm%pw%pw_grid%dvol*v_qmmm%pw%cr3d

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_modify_hartree_pot

!!****f* qmmm_methods/qmmm_forces [1.0] *
!!
!!   NAME
!!     qmmm_forces
!!
!!   FUNCTION
!!     General driver to Compute the contribution 
!!     to the forces due to the QM/MM potential
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces(qs_env,qmmm_env,mm_particles,calc_force,mm_cell,&
       globenv, error)
    TYPE(global_environment_type), &
      INTENT(in), OPTIONAL, TARGET           :: globenv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces', &
      routineP = moduleN//':'//routineN    

    INTEGER                                  :: handle, Imm, IndMM, ispin, &
                                                stat, iw
    LOGICAL                                  :: failure, need_f
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces, Forces_added_charges
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pool
    TYPE(pw_type), POINTER                   :: rho_tot_r
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input_section, interp_section

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    need_f = .TRUE.
    iw = globenv%scr
    IF (PRESENT(calc_force)) need_f = calc_force
    NULLIFY(ks_qmmm_env_loc, rho, pw_env, rho_tot_r, energy, Forces, Forces_added_charges,&
            input_section)
    CALL get_qs_env(qs_env=qs_env,&
         rho=rho,&
         rho_core=rho_core,&
         pw_env=pw_env,&
         energy=energy,&
         para_env=para_env,&
         input=input_section,&
         error=error)
    ks_qmmm_env_loc => qs_env%ks_qmmm_env
    interp_section => section_vals_get_subs_vals(input_section,"DFT",error=error)
    interp_section => section_vals_get_subs_vals(interp_section,"MGRID",&
         error=error)
    interp_section => section_vals_get_subs_vals(interp_section,"INTERPOLATOR",&
         error=error)
    !
    CALL pw_env_get(pw_env=pw_env,&
         pw_pools=pw_pools,& 
         auxbas_pw_pool=auxbas_pool,&
         error=error)
    CALL pw_pool_create_pw(auxbas_pool,rho_tot_r,&
         in_space=REALSPACE,use_data=REALDATA3D, error=error)
    CALL pw_transfer(rho_core%pw,rho_tot_r)
    !
    ! Computes the QM/MM Nuclear Electrostatic Potential
    !
    energy%qmmm_nu = pw_integral_ab ( rho_tot_r, ks_qmmm_env_loc%v_qmmm_rspace%pw )
    energy%total = energy%total +  energy%qmmm_nu
    WRITE(iw,'(T3,A,F15.9)')" QM/MM Nuclear Electrostatic Potential :",energy%qmmm_nu
    WRITE(iw,'(T3,A,F15.9)')" Total Energy after the QM/MM Calculation :",energy%total
    IF (need_f) THEN
       ALLOCATE(Forces(qmmm_env%num_mm_atoms,3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Forces_added_charges(qmmm_env%added_charges%num_mm_atoms,3),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !
       DO ispin=1,SIZE(rho%rho_r)
          CALL pw_sumup(rho%rho_r(ispin)%pw,rho_tot_r)
       END DO
       WRITE(iw,'(T3,A)')" Evaluating forces on MM atoms due to the:"
       !
       ! Electrostatic Interaction type...
       !
       Forces               = 0.0_dp
       Forces_added_charges = 0.0_dp
       SELECT CASE(qmmm_env%qmmm_coupl_type)
       CASE(do_qmmm_grid)
          WRITE(iw,'(T3,A)')"QM/MM Coupling computed on GRID."
          ! Compute the MM electrostatic potential on GRID
          CALL qmmm_forces_on_grid     ( rho=rho_tot_r,&
               v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
               qmmm_env=qmmm_env,&
               mm_particles=mm_particles,&
               dvol=rho%rho_r(1)%pw%pw_grid%dvol,&
               Forces=Forces,&
               Forces_added_charges=Forces_added_charges,&
               mm_cell=mm_cell,&
               error=error)
       CASE(do_qmmm_gauss)
          WRITE(iw,'(T3,A)')"QM/MM Coupling computed collocating the Gaussian Potential Functions."
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
               mm_cell=mm_cell,&
               error=error)
       CASE(do_qmmm_none)
          WRITE(iw,'(T3,A)')"No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
       CASE DEFAULT
          WRITE(iw,'(T3,A)')"Unknown Coupling..."
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       !
       ! Debug Forces
       !
       IF (debug_this_module) THEN
          !
          ! Print Out Forces
          !
          WRITE(iw,*)"ANALYTICAL FORCES:"
          DO Imm = 1, SIZE(qmmm_env%mm_atom_index)
             IndMM = qmmm_env%mm_atom_index(Imm)
             WRITE(iw,'(I6,3F15.9)') IndMM, Forces(Imm,:)
          END DO
          CALL qmmm_debug_forces(rho=rho_tot_r,&
               qs_env=qs_env,&
               qmmm_env=qmmm_env,&
               Analytical_Forces=Forces,&
               mm_particles=mm_particles,&
               mm_atom_index=qmmm_env%mm_atom_index,&
               num_mm_atoms=qmmm_env%num_mm_atoms,&
               Potentials=qmmm_env%Potentials,&
               interp_section=interp_section,&
               mm_cell=mm_cell,&
               globenv=globenv,&
               error=error)                              
       ENDIF
       !
       ! Transfer Forces 
       !
       DO Imm = 1, qmmm_env%num_mm_atoms
          IndMM = qmmm_env%mm_atom_index(Imm)
          ! Hack: In Forces there the gradients indeed...
          ! Minux sign to take care of this misunderstanding...
          mm_particles(IndMM)%f(:) = - Forces(Imm,:) + mm_particles(IndMM)%f(:)
       END DO
       DEALLOCATE(Forces, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN
          DO Imm = 1, qmmm_env%added_charges%num_mm_atoms
             IndMM = qmmm_env%added_charges%mm_atom_index(Imm)
             ! Hack: In Forces there the gradients indeed...
             ! Minux sign to take care of this misunderstanding...
             qmmm_env%added_charges%added_particles(IndMM)%f(:) = - Forces_added_charges(Imm,:)
          END DO
       END IF
       DEALLOCATE(Forces_added_charges, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    ! Give back rho_tot_t to auxbas_pool
    CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_r,&
         accept_non_compatible=.TRUE.,error=error)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_forces

!!****f* qmmm_methods/qmmm_forces_on_grid [1.0] *
!!
!!   NAME
!!     qmmm_forces_on_grid
!!
!!   FUNCTION
!!     Evaluates the contribution 
!!     to the forces due to the QM/MM potential computed on GRID
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_on_grid(rho,v_qmmm,qmmm_env,mm_particles,&
       dvol, Forces, Forces_added_charges, mm_cell, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(pw_type), POINTER                   :: rho
    TYPE(coeff_type), INTENT(INOUT)          :: v_qmmm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), INTENT(IN)                :: dvol
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces, Forces_added_charges
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_on_grid', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    CALL qmmm_forces_on_grid_low(v_qmmm,&
                                 rho%cr3d,&
                                 mm_particles,&
                                 qmmm_env%mm_atom_chrg,&
                                 qmmm_env%mm_el_pot_radius,&
                                 qmmm_env%mm_atom_index,&
                                 qmmm_env%num_mm_atoms,&
                                 dvol,&
                                 Forces,&
                                 qmmm_env%potentials,&
                                 mm_cell,&
                                 qmmm_env%dOmmOqm,&
                                 error)

    IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN

       CALL qmmm_forces_on_grid_low(v_qmmm,&
                                    rho%cr3d,&
                                    qmmm_env%added_charges%added_particles,&
                                    qmmm_env%added_charges%mm_atom_chrg,&
                                    qmmm_env%added_charges%mm_el_pot_radius,&
                                    qmmm_env%added_charges%mm_atom_index,&
                                    qmmm_env%added_charges%num_mm_atoms,&
                                    dvol,&
                                    Forces_added_charges,&
                                    qmmm_env%added_charges%potentials,&
                                    mm_cell,&
                                    qmmm_env%dOmmOqm,&
                                    error)       

    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_forces_on_grid

!!****f* qmmm_methods/qmmm_forces_on_grid_low [1.0] *
!!
!!   NAME
!!     qmmm_forces_on_grid_low
!!
!!   FUNCTION
!!     Evaluates the contribution to the forces due to the 
!!     QM/MM potential computed on GRID. Low Level Routine.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_on_grid_low(v_qmmm,Rho,mm_particles,mm_charges,&
       mm_el_pot_radius, mm_atom_index,num_mm_atoms,dvol,Forces,Potentials,&
       mm_cell, dOmmOqm, error)
    TYPE(coeff_type), INTENT(INOUT)          :: v_qmmm
    REAL(KIND=dp), DIMENSION(*)              :: Rho
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    REAL(KIND=dp), INTENT(IN)                :: dvol
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Forces
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: Potentials
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_on_grid_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, Imm, Imp, IndMM, &
                                                Ipot, ix, j, jj, k, my_j, &
                                                my_k, n1, n2, n3
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dpot, dr1, dr2, dr3, dx, fac, &
                                                ft1, ft2, ft3, qt, r, rc, &
                                                rt1, rt2, rt3, rv1, rv2, rv3, &
                                                rx, rx2, rx3, rx4, xs1, xs2, &
                                                xs3
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pot0_2
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc
    TYPE(qmmm_pot_type), POINTER             :: Pot

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_charges),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_atom_index),cp_failure_level,routineP,error,failure)

    bo = v_qmmm%pw%pw_grid%bounds_local  
    n1 = v_qmmm%pw%pw_grid%npts(1)   
    n2 = v_qmmm%pw%pw_grid%npts(2)   
    n3 = v_qmmm%pw%pw_grid%npts(3)   
    dr1 = v_qmmm%pw%pw_grid%dr(1)
    dr2 = v_qmmm%pw%pw_grid%dr(2)
    dr3 = v_qmmm%pw%pw_grid%dr(3)
    CALL timeset(routineN,"I","",handle)
    MainLoopPot: DO Ipot = 1, SIZE(Potentials)
       Pot    => Potentials(Ipot)%Pot
       dx     =  Pot%dx 
       pot0_2 => Pot%pot0_2
       LoopMM: DO Imp = 1, SIZE(Pot%mm_atom_index)
          Imm = Pot%mm_atom_index(Imp)
          IndMM = mm_atom_index(Imm) 
          rc = mm_el_pot_radius(Imm) 
          r_pbc = pbc(mm_particles(IndMM)%r-dOmmOqm, mm_cell)+dOmmOqm
          rt1= r_pbc(1) 
          rt2= r_pbc(2) 
          rt3= r_pbc(3)
          ft1 =0.0_dp
          ft2 =0.0_dp
          ft3 =0.0_dp
          qt=mm_charges(Imm)*dvol/dx
          !$omp parallel do default (none) &
          !$omp private(i,j,k,Imm,IndMM,jj,my_k,my_j,&
          !$omp         r, rc, qt,ix,rx,rx2,rx3,rx4, & 
          !$omp         dr1,dr2,dr3, xs1,xs2,xs3, rv1,rv2,rv3, rt1,rt2,rt3,&
          !$omp         ft1,ft2,ft3,dpot,fac)&
          !$omp shared(Rho,bo,n1,n2,n3,dx,pot0_2)
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k=(k+n3/2)
             xs3 = REAL(my_k,dp)*dr3
             my_j=bo(1,2)+n2/2
             jj=my_k*n1*n2+my_j*n1
             xs2 = REAL(my_j,dp)*dr2
             rv3 = rt3 - xs3
             DO j = bo(1,2), bo(2,2)
                xs1 = (bo(1,1)+n1/2)*dr1
                rv2 = rt2 - xs2
                DO i = bo(1,1), bo(2,1)
                   jj=jj+1
                   rv1  = rt1 - xs1
                   r    = SQRT(rv1*rv1 + rv2*rv2 + rv3*rv3)
                   ix   = FLOOR(r/dx)+1
                   rx   = (r-REAL(ix-1,dp)*dx)/dx
                   rx2  = rx*rx
                   rx3  = rx2*rx
                   rx4  = rx3*rx
                   dpot = &
                         pot0_2(1,ix  )*(-30.0_dp*(rx2+rx4)+60.0_dp*rx3)               &
                        +pot0_2(2,ix  )*(1.0_dp-18.0_dp*rx2+32.0_dp*rx3-15.0_dp*rx4)   &
                        +pot0_2(3,ix  )*(rx-4.5_dp*rx2+6.0_dp*rx3-2.5_dp*rx4)        &
                        +pot0_2(1,ix+1)*(30.0_dp*(rx2+rx4)-60.0_dp*rx3)                &
                        +pot0_2(2,ix+1)*(-12.0_dp*rx2+28.0_dp*rx3-15.0_dp*rx4)         &
                        +pot0_2(3,ix+1)*(1.5_dp*rx2-4.0_dp*rx3+2.5_dp*rx4)

                   dpot = dpot / r
                   fac = - Rho(jj)* dpot
                   ft1 = ft1 + fac * rv1
                   ft2 = ft2 + fac * rv2
                   ft3 = ft3 + fac * rv3
                   xs1 = xs1 + dr1
                END DO
                xs2 = xs2 + dr2
             END DO
          END DO LoopOnGrid
          Forces(Imm,1) =  ft1 * qt
          Forces(Imm,2) =  ft2 * qt
          Forces(Imm,3) =  ft3 * qt
       END DO LoopMM
    END DO MainLoopPot
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_forces_on_grid_low
!!****f* qmmm_methods/qmmm_forces_with_gaussian [1.0] *
!!
!!   NAME
!!     qmmm_forces_with_gaussian
!!
!!   FUNCTION
!!     Evaluates the contribution to the forces due to the 
!!     QM/MM potential computed collocating the Electrostatic
!!     Gaussian Potential.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_with_gaussian(rho, qmmm_env, mm_particles,  &
       aug_pools, auxbas_grid, coarser_grid, cube_info, para_env,    &
       eps_mm_rspace, pw_pools, Forces, Forces_added_charges,        &
       interp_section, mm_cell,error)
    IMPLICIT NONE
    ! Arguments
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
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_with_gaussian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, igrid, j, k, ngrids, stat
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids

! Statements

    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
    NULLIFY(grids)
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qmmm_env%mm_atom_chrg),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qmmm_env%mm_atom_index),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(Forces),cp_failure_level,routineP,error,failure)
    !Statements
    ngrids=SIZE(pw_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw)
    END DO

    ! Collocate Density on multigrids

    lb = rho%pw_grid%bounds_local(1,:)
    ub = rho%pw_grid%bounds_local(2,:)
    grids(auxbas_grid) % pw % cr3d(lb(1):ub(1),&
                                   lb(2):ub(2),&
                                   lb(3):ub(3) )= rho % cr3d
    ! copy the boundaries
     DO k=lb(3),ub(3)
        DO j=lb(2),ub(2)
           grids(auxbas_grid) % pw %cr3d(ub(1)+1,j,k)=rho%cr3d(lb(1),j,k)
        END DO
        grids(auxbas_grid) % pw %cr3d(ub(1)+1,ub(2)+1,k)=rho%cr3d(lb(1),lb(2),k)
     END DO
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
        grids(auxbas_grid) % pw %cr3d(ub(1)+1,j,ub(3)+1)=rho%cr3d(lb(1),j,lb(3))
     END DO
     grids(auxbas_grid) % pw %cr3d(ub(1)+1,ub(2)+1,ub(3)+1)=rho%cr3d(lb(1),lb(2),lb(3))

    ! RealSpace Interpolation
    IF (.FALSE.) THEN
       ! Linear Interpolator
       DO Igrid = auxbas_grid, SIZE(grids)-1
          CALL pw_inject_NoPBC(grids(Igrid  ) % pw,&
                               grids(Igrid+1) % pw,&
                               error)
          grids(igrid+1) % pw % cr3d = grids(igrid+1) % pw % cr3d / REAL(8,dp)
       END DO
    ELSE
       ! Spline Interpolator
       
       DO Igrid =  auxbas_grid, SIZE(grids)-1
          CALL pw_project_s3(grids(Igrid  ) % pw,&
                             grids(Igrid+1) % pw,&
                             aug_pools(Igrid+1)%pool,&
                             param_section=interp_section,&
                             error=error)
       END DO
    END IF

    CALL qmmm_forces_with_gaussian_low(grids, mm_particles,&
         qmmm_env%mm_atom_chrg, qmmm_env%mm_el_pot_radius, qmmm_env%mm_atom_index,&
         qmmm_env%num_mm_atoms, cube_info, para_env, eps_mm_rspace, auxbas_grid, &
         coarser_grid, qmmm_env%pgfs, qmmm_env%potentials, Forces, aug_pools,    &
         mm_cell, qmmm_env%dOmmOqm, error)

    IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN

       CALL qmmm_forces_with_gaussian_low(grids, qmmm_env%added_charges%added_particles,&
            qmmm_env%added_charges%mm_atom_chrg, qmmm_env%added_charges%mm_el_pot_radius,&
            qmmm_env%added_charges%mm_atom_index, qmmm_env%added_charges%num_mm_atoms,&
            cube_info, para_env, eps_mm_rspace, auxbas_grid, coarser_grid, qmmm_env%added_charges%pgfs,&
            qmmm_env%added_charges%potentials, Forces_added_charges, aug_pools, mm_cell, &
            qmmm_env%dOmmOqm, error)

    END IF

    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_forces_with_gaussian

!!****f* qmmm_methods/qmmm_forces_with_gaussian_low [1.0] *
!!
!!   NAME
!!     qmmm_forces_with_gaussian_low
!!
!!   FUNCTION
!!     Evaluates the contribution to the forces due to the 
!!     QM/MM potential computed collocating the Electrostatic
!!     Gaussian Potential. Low Level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_with_gaussian_low(grids, mm_particles,&
       mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms,&
       cube_info, para_env, eps_mm_rspace, auxbas_grid, &
       coarser_grid, pgfs, potentials, Forces, aug_pools, mm_cell, &
       dOmmOqm, error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'qmmm_forces_with_gaussian_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, Icg, IGauss, ilevel, &
                                                Imm, IndMM, IRadTyp, ithread, &
                                                LIndMM, stat
    INTEGER, DIMENSION(3)                    :: ndim
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, dvol, height, W
    REAL(KIND=dp), DIMENSION(3)              :: force, ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: wrk, xdat, xder, ydat, yder, &
                                                zdat, zder
    TYPE(qmmm_gaussian_type), POINTER        :: pgf
    TYPE(qmmm_pot_type), POINTER             :: pot

    NULLIFY(pgf, xdat, ydat, zdat, xder, yder, zder, wrk)    
    CALL timeset(routineN,"I","",handle)
    ithread = 0
    Icg     = 0
    ndim    = grids(auxbas_grid) % pw % pw_grid %npts
    ALLOCATE (xdat(ndim(1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (ydat(ndim(2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zdat(ndim(3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (xder(ndim(1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (yder(ndim(2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zder(ndim(3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (wrk(MAX(ndim(1),ndim(2),ndim(3))), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
!
    Radius: DO IRadTyp = 1, SIZE(pgfs) 
       pgf => pgfs(IRadTyp)%pgf
       pot => potentials(IRadTyp)%pot
       Gaussian: DO IGauss = 1, pgf%Number_of_Gaussians
          alpha     = 1.0_dp / pgf%Gk(IGauss)
          alpha     = alpha * alpha
          height    = pgf%Ak(IGauss)
          ilevel    = pgf%grid_level(IGauss)
          dvol      = grids(ilevel)%pw%pw_grid%dvol
          Atoms: DO Imm = 1, SIZE(pot%mm_atom_index)
             Icg       =   Icg + 1
             LIndMM    =   pot%mm_atom_index(Imm)
             IndMM     =   mm_atom_index(LIndMM)
             ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm,mm_cell)+dOmmOqm
             W         =   mm_charges(LIndMM) * height
             force     =   0.0_dp
             IF (W == 0.0_dp) CYCLE Atoms
             CALL integrate_gf_rspace_NoPBC(zetp=alpha,&
                                            rp=ra,&
                                            scale=-1.0_dp,&
                                            W=W,&
                                            pwgrid=grids(ilevel)%pw,&
                                            cube_info=cube_info(ilevel),&
                                            eps_mm_rspace=eps_mm_rspace,&
                                            ithread=ithread, &
                                            xdat=xdat,&
                                            ydat=ydat,&
                                            zdat=zdat,&
                                            xder=xder,&
                                            yder=yder,&
                                            zder=zder,&
                                            wrk=wrk,&
                                            force=force)
             force = force * dvol
             Forces(LIndMM,:) = Forces(LIndMM,:) + force(:)
             !
             ! Debug Statement
             !
             IF (debug_this_module) THEN
                CALL debug_integrate_gf_rspace_NoPBC(ilevel=ilevel,&
                                                     zetp=alpha,&
                                                     rp=ra,&
                                                     scale=1.0_dp,&
                                                     W=W,&
                                                     pwgrid=grids(ilevel)%pw,&
                                                     cube_info=cube_info(ilevel),&
                                                     eps_mm_rspace=eps_mm_rspace,&
                                                     ithread=ithread,&
                                                     aug_pools=aug_pools,&
                                                     debug_force=force,&
                                                     mm_cell=mm_cell,&
                                                     error=error)
             END IF
          END DO Atoms
       END DO Gaussian
    END DO Radius

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
    IF (ASSOCIATED(xder)) THEN
       DEALLOCATE (xder, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(yder)) THEN
       DEALLOCATE (yder, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(zder)) THEN
       DEALLOCATE (zder, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    IF (ASSOCIATED(wrk)) THEN
       DEALLOCATE (wrk, STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
    
    CALL qmmm_forces_with_gaussian_LR  (pgfs=pgfs,&
                                        grid=grids(coarser_grid)%pw,&
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
                                        error=error)

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_forces_with_gaussian_low


!!****f* qmmm_methods/qmmm_forces_with_gaussian_LR [1.0] *
!!
!!   NAME
!!     qmmm_forces_with_gaussian_LR
!!
!!   FUNCTION
!!     Evaluates the contribution to the forces due to the Long Range
!!     part of the QM/MM potential computed collocating the Electrostatic
!!     Gaussian Potential.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_forces_with_gaussian_LR  (pgfs,grid,num_mm_atoms,mm_charges,mm_atom_index,&
                                            mm_particles,mm_el_pot_radius,para_env,&
                                            coarser_grid_level,Forces, potentials,&
                                            aug_pools, mm_cell, dOmmOqm, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_type), POINTER                   :: grid
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_forces_with_gaussian_LR', &
      routineP = moduleN//':'//routineN    

    INTEGER                                  :: handle, i, Imm, IndMM, &
                                                IRadTyp, ix, j, k, LIndMM, &
                                                my_i, my_j, my_k, n1, n2, n3, &
                                                stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dr1, dr2, dr3, dvol, dx, fac, &
                                                ft1, ft2, ft3, qt, r, r2, &
                                                rt1, rt2, rt3, rv1, rv2, rv3, &
                                                rx, rx2, Term, xs1, xs2, xs3
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: LForces
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pot0_2
    TYPE(qmmm_pot_type), POINTER             :: pot

! Statements

    CALL timeset(routineN,"I","",handle)

    ALLOCATE(LForces(num_mm_atoms,3), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    n1   = grid%pw_grid%npts(1)
    n2   = grid%pw_grid%npts(2)
    n3   = grid%pw_grid%npts(3)
    dr1  = grid%pw_grid%dr(1) 
    dr2  = grid%pw_grid%dr(2) 
    dr3  = grid%pw_grid%dr(3) 
    dvol = grid%pw_grid%dvol
    bo   = grid%pw_grid%bounds

    Radius: DO IRadTyp = 1, SIZE(pgfs) 
       pot => potentials(IRadTyp)%pot
       dx     =  Pot%dx 
       pot0_2 => Pot%pot0_2
       Atoms: DO Imm = 1, SIZE(pot%mm_atom_index)
          LIndMM    =   pot%mm_atom_index(Imm)
          IndMM     =   mm_atom_index(LIndMM)
          ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm,mm_cell)+dOmmOqm
          qt        =   mm_charges(LIndMM)
          IF (qt == 0.0_dp) CYCLE Atoms
          rt1 = ra(1)
          rt2 = ra(2)
          rt3 = ra(3)
          ft1 = 0.0_dp
          ft2 = 0.0_dp
          ft3 = 0.0_dp
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k = (k+n3/2)
             xs3  = REAL(my_k,dp)*dr3
             my_j = bo(1,2)+n2/2
             xs2 = REAL(my_j,dp)*dr2
             rv3 = rt3 - xs3
             DO j =  bo(1,2), bo(2,2)
                my_i= bo(1,1)+n1/2
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

                   fac = grid%cr3d(i,j,k) * Term / r 
                   ft1 = ft1 + fac * rv1
                   ft2 = ft2 + fac * rv2
                   ft3 = ft3 + fac * rv3
                   xs1 = xs1 + dr1
                END DO
                xs2 = xs2 + dr2
             END DO
          END DO LoopOnGrid
          qt = - qt * dvol / dx
          LForces(LindMM,1) = ft1 * qt
          LForces(LindMM,2) = ft2 * qt
          LForces(LindMM,3) = ft3 * qt

          Forces(LIndMM,1) = Forces(LIndMM,1) + LForces(LindMM,1) 
          Forces(LIndMM,2) = Forces(LIndMM,2) + LForces(LindMM,2)
          Forces(LIndMM,3) = Forces(LIndMM,3) + LForces(LindMM,3)
       END DO Atoms
    END DO Radius
    !
    ! Debug Statement
    !
    IF (debug_this_module) THEN
       CALL debug_qmmm_forces_with_gauss_LR(pgfs=pgfs,&
                                            aug_pools=aug_pools,&
                                            rho=grid,&
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
                                            error=error)
    END IF

    DEALLOCATE(LForces, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)    
  END SUBROUTINE qmmm_forces_with_gaussian_LR

!!****f* qmmm_methods/qmmm_debug_forces [1.0] *
!!
!!   NAME
!!     qmmm_debug_forces
!!
!!   FUNCTION
!!     Evaluates numerically QM/MM forces and compares them with
!!     the analytically computed ones.
!!     It is evaluated only when debug_this_module is set to .TRUE.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE qmmm_debug_forces(rho,qs_env,qmmm_env,Analytical_Forces,&
                               mm_particles,mm_atom_index,num_mm_atoms,&
                               Potentials,interp_section,mm_cell,globenv,error)      
    TYPE(pw_type), POINTER                   :: rho
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Analytical_Forces
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: Potentials
    TYPE(section_vals_type), POINTER         :: interp_section
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(global_environment_type), &
      INTENT(in), OPTIONAL, TARGET           :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_debug_forces', &
      routineP = moduleN//':'//routineN      

    INTEGER                                  :: handle, I, IndMM, J, K, stat, iw
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: Coord_save
    REAL(KIND=dp), DIMENSION(2)              :: energy
    REAL(KIND=dp), DIMENSION(3)              :: Err
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: Num_Forces
    TYPE(coeff_type)                         :: v_qmmm_rspace
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc

! Statements
    iw = globenv%scr
    CALL timeset(routineN,"I","",handle)
    NULLIFY( Num_Forces )
    CALL get_qs_env(qs_env=qs_env,&
                    pw_env=pw_env,&
                    para_env=para_env,&
                    error=error)

    CALL pw_env_get(pw_env=pw_env, pw_pools=pw_pools, error=error)
    CALL pw_pool_init_coeff(pw_pools(1)%pool, v_qmmm_rspace,&
         use_data=REALDATA3D, in_space=REALSPACE)
    ALLOCATE(Num_Forces(num_mm_atoms,3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ks_qmmm_env_loc =>  qs_env%ks_qmmm_env
    WRITE(iw,'(/A)')"DEBUG SECTION:"
    Atoms: DO I = 1, num_mm_atoms
       IndMM = mm_atom_index(I)
       Coords: DO J = 1, 3
          Coord_save = mm_particles(IndMM)%r(J)
          energy = 0.0_dp
          Diff: DO K = 1, 2             
             mm_particles(IndMM)%r(J) = Coord_save + (-1)**K * Dx
             CALL coeff_zero(coeff=v_qmmm_rspace)
             SELECT CASE(qmmm_env%qmmm_coupl_type)
             CASE(do_qmmm_grid)
                CALL qmmm_electrostatic_on_grid(qmmm_env=qmmm_env,&
                                                v_qmmm=v_qmmm_rspace,&
                                                mm_particles=mm_particles,&
                                                mm_cell=mm_cell,&
                                                error=error)
             CASE(do_qmmm_gauss)
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
                WRITE(iw,'(T3,A)')"Unknown Coupling..."
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             END SELECT
             energy(K) =  pw_integral_ab ( rho, v_qmmm_rspace%pw )
          END DO Diff
          WRITE(iw,'(A,I6,A,I3,A,2F15.9)')&
          "DEBUG :: MM Atom = ",IndMM," Coord = ",J," Energies (+/-) :: ",energy(2), energy(1)
          Num_Forces(I,J)  = ( energy(2) - energy(1) ) / (2.0_dp * Dx)
          mm_particles(IndMM)%r(J) = Coord_save
       END DO Coords       
    END DO Atoms
    
    SELECT CASE(qmmm_env%qmmm_coupl_type)
    CASE(do_qmmm_grid,do_qmmm_gauss)
       WRITE(iw,'(/A/)')"CHECKING NUMERICAL Vs ANALYTICAL FORCES (Err%):"
       DO I =  1, num_mm_atoms
          IndMM = mm_atom_index(I)
          Err = 0.0_dp
          DO K = 1, 3
             IF (ABS(Num_Forces(I,K)) >= 1.0E-6_dp) THEN 
                Err(K) = (Analytical_Forces(I,K)-Num_Forces(I,K))/Num_Forces(I,K)*100.0_dp
             END IF
          END DO
          WRITE(iw,100)IndMM,Analytical_Forces(I,1),Err(1),&
                             Analytical_Forces(I,2),Err(2),&
                             Analytical_Forces(I,3),Err(3) 
          CPPostcondition(ABS(Err(1))<=MaxErr,cp_failure_level,routineP,error,failure)
          CPPostcondition(ABS(Err(2))<=MaxErr,cp_failure_level,routineP,error,failure)
          CPPostcondition(ABS(Err(3))<=MaxErr,cp_failure_level,routineP,error,failure)
       END DO
    CASE(do_qmmm_none)
       WRITE(iw,'(T3,A)')"No QM/MM Derivatives to debug. Just Mechanical Coupling!"
    CASE DEFAULT
       WRITE(iw,'(T3,A)')"Unknown Coupling..."
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL pw_pool_give_back_coeff ( pw_pools(1)%pool, v_qmmm_rspace)
    DEALLOCATE(Num_Forces,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)    
100 FORMAT(I5,F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ")
  END SUBROUTINE qmmm_debug_forces


!!****f* qmmm_methods/debug_integrate_gf_rspace_NoPBC [1.0] *
!!
!!   NAME
!!     debug_integrate_gf_rspace_NoPBC
!!
!!   FUNCTION
!!     Debugs the integrate_gf_rspace_NoPBC.. It may helps ;-P
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE debug_integrate_gf_rspace_NoPBC(ilevel, zetp, rp, scale, W, pwgrid, cube_info,&
                                             eps_mm_rspace, ithread, aug_pools, debug_force,&
                                             mm_cell,error)
    INTEGER, INTENT(IN)                      :: ilevel
    REAL(KIND=dp), INTENT(IN)                :: zetp
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rp
    REAL(KIND=dp), INTENT(IN)                :: scale, W
    TYPE(pw_type), POINTER                   :: pwgrid
    TYPE(cube_info_type), INTENT(IN)         :: cube_info
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    INTEGER, INTENT(IN), OPTIONAL            :: ithread
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: debug_force
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'debug_integrate_gf_rspace_NoPBC', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, igrid, k, ngrids, &
                                                stat
    INTEGER, DIMENSION(3)                    :: ndim
    INTEGER, SAVE                            :: Icount
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(2)              :: energy
    REAL(KIND=dp), DIMENSION(3)              :: Err, force, myrp, r_pbc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdat, ydat, zdat
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids

    DATA Icount /0/
    ! Statements
    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
    NULLIFY(grids)
    !Statements
    ngrids = SIZE(aug_pools)
    ndim   = 0
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw)
    END DO
    DO igrid = 1, ngrids
       ndim = MAX(grids(igrid)%pw%pw_grid%npts,ndim)
    END DO
    ALLOCATE (xdat(ndim(1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (ydat(ndim(2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zdat(ndim(3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    

    Icount = Icount + 1
    DO i = 1, 3
       DO k = 1, 2
          myrp = rp
          myrp(i) = myrp(i) + (-1.0_dp)**k * Dx
          r_pbc = myrp
          CALL pw_zero(grids(ilevel)%pw)
          CALL collocate_gf_rspace_NoPBC(zetp=zetp,&
                                         rp=myrp,&
                                         scale=-1.0_dp,&
                                         W=W,&
                                         pwgrid=grids(ilevel)%pw,&
                                         cube_info=cube_info,&
                                         eps_mm_rspace=eps_mm_rspace,&
                                         ithread=ithread, &
                                         xdat=xdat,&
                                         ydat=ydat,&
                                         zdat=zdat)

          energy(k) = pw_integral_ab(pwgrid, grids(ilevel)%pw)          
       END DO
       force(i) = ( energy(2) - energy(1) ) / (2.0_dp * Dx)
    END DO
    Err = 0.0_dp
    IF (ALL(force /= 0.0_dp)) THEN 
       Err(1) = (debug_force(1)-force(1))/force(1)*100.0_dp
       Err(2) = (debug_force(2)-force(2))/force(2)*100.0_dp
       Err(3) = (debug_force(3)-force(3))/force(3)*100.0_dp
    END IF
    WRITE(*,100)Icount, debug_force(1), Err(1),&
                        debug_force(2), Err(2),&
                        debug_force(3), Err(3)
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
    CALL timestop(0.0_dp,handle) 
100 FORMAT("Collocation   : ",I5,F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ")
  END SUBROUTINE debug_integrate_gf_rspace_NoPBC



!!****f* qmmm_methods/debug_qmmm_forces_with_gauss_LR [1.0] *
!!
!!   NAME
!!     debug_qmmm_forces_with_gauss_LR
!!
!!   FUNCTION
!!     Debugs qmmm_forces_with_gaussian_LR.. It may helps too ... ;-]
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [tlaino]
!!
!!*************************************************************************
  SUBROUTINE debug_qmmm_forces_with_gauss_LR(pgfs,aug_pools, rho, mm_charges, mm_atom_index,&
                                             mm_particles, mm_el_pot_radius, num_mm_atoms,&
                                             coarser_grid_level, potentials,&
                                             debug_force, para_env, mm_cell,dOmmOqm,error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(pw_type), POINTER                   :: rho
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

    ALLOCATE(Num_Forces(num_mm_atoms,3),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(grids)
    CALL timeset(routineN,"I","",handle)
    ngrids = SIZE(aug_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw)
    END DO
    Atoms: DO I = 1, num_mm_atoms
       IndMM = mm_atom_index(I)
       Coords: DO J = 1, 3
          Coord_save = mm_particles(IndMM)%r(J)
          energy = 0.0_dp
          Diff: DO K = 1, 2             
             mm_particles(IndMM)%r(J) = Coord_save + (-1)**K * Dx
             CALL pw_zero(grids(coarser_grid_level)%pw)

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
                                                 error=error)             

             energy(K) =  pw_integral_ab ( rho, grids(coarser_grid_level)%pw )
          END DO Diff
          WRITE(*,'(A,I6,A,I3,A,2F15.9)')&
          "DEBUG LR:: MM Atom = ",IndMM," Coord = ",J," Energies (+/-) :: ",energy(2), energy(1)
          Num_Forces(I,J)  = ( energy(2) - energy(1) ) / (2.0_dp * Dx)
          mm_particles(IndMM)%r(J) = Coord_save
       END DO Coords       
    END DO Atoms
    
    DO I =  1, num_mm_atoms
       IndMM = mm_atom_index(I)
       Err = 0.0_dp
       IF (ALL(Num_Forces /= 0.0_dp)) THEN 
          Err(1) = (debug_force(I,1)-Num_Forces(I,1))/Num_Forces(I,1)*100.0_dp
          Err(2) = (debug_force(I,2)-Num_Forces(I,2))/Num_Forces(I,2)*100.0_dp
          Err(3) = (debug_force(I,3)-Num_Forces(I,3))/Num_Forces(I,3)*100.0_dp
       END IF
       WRITE(*,100)IndMM,debug_force(I,1),Err(1),&
                         debug_force(I,2),Err(2),&
                         debug_force(I,3),Err(3) 
       CPPostcondition(ABS(Err(1))<=MaxErr,cp_failure_level,routineP,error,failure)
       CPPostcondition(ABS(Err(2))<=MaxErr,cp_failure_level,routineP,error,failure)
       CPPostcondition(ABS(Err(3))<=MaxErr,cp_failure_level,routineP,error,failure)
    END DO

    DEALLOCATE(Num_Forces,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)
    CALL timestop(0.0_dp,handle)
100 FORMAT("MM Atom LR    : ",I5,F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ",F15.9," ( ",F7.2," ) ")
  END SUBROUTINE debug_qmmm_forces_with_gauss_LR


END MODULE qmmm_methods
