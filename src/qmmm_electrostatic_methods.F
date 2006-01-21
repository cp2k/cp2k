!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qmmm_electrostatic_methods [1.0] *
!!
!!   NAME MODULE
!!     qmmm_electrostatic_methods
!!
!!   FUNCTION
!!     A collection of methods to treat the QM/MM electrostatic coupling
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     5.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_electrostatic_methods
  USE coefficient_types,               ONLY: coeff_type,&
                                             coeff_zero
  USE mathconstants,                   ONLY: pi
  USE cp_output_handling, ONLY: cp_print_key_should_output,&
       cp_p_file, cp_print_key_unit_nr, cp_print_key_finished_output,&
       cp_iter_string, cp_iterate,&
       cp_add_iter_level, cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cube_utils,                      ONLY: cube_info_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_qmmm_gauss,&
                                             do_qmmm_grid,&
                                             do_qmmm_none,&
                                             linear_interp,&
                                             spline3_nopbc_interp,&
                                             spline3_pbc_interp
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_val_get,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mm_collocate_potential,          ONLY: collocate_gf_rspace_NoPBC
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pools_create_pws,&
                                             pw_pools_give_back_pws,&
                                             pw_pool_create_pw,&
                                             pw_pool_give_back_pw
  USE pw_spline_utils,                 ONLY: pw_prolongate_s3,&
                                             spl3_nopbct,&
                                             spl3_1d_transf_coeffs,&
                                             spl3_1d_transf_border1,&
                                             add_fine2coarse,&
                                             Eval_Interp_Spl3_pbc
  USE pw_types,                        ONLY: REALDATA3D,RECIPROCALSPACE,&
                                             REALSPACE,COMPLEXDATA1D,&
                                             pw_p_type,&
                                             pw_prolongate_l_NoPBC,&
                                             pw_type,&
                                             pw_zero,&
                                             pw_compare_debug,&
                                             pw_transfer,pw_sumup
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_pot_p_type,&
                                             qmmm_pot_type,&
                                             qmmm_per_pot_p_type,&
                                             qmmm_per_pot_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE cell_types,                 ONLY: cell_type,&
                                             pbc
  USE message_passing, ONLY: mp_stop, mp_sync
  USE machine, only: m_flush
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_electrostatic_methods'

  PUBLIC :: qmmm_el_coupling
  PUBLIC :: qmmm_electrostatic_on_grid
  PUBLIC :: qmmm_elec_with_gaussian,qmmm_elec_with_gaussian_low,&
            qmmm_elec_with_gaussian_LR, qmmm_elec_with_gaussian_LG
!***
!****************************************************************************
CONTAINS
!!****f* qmmm_electrostatic_methods/qmmm_el_coupling [1.0] *
!!
!!   NAME
!!     qmmm_el_coupling
!!
!!   FUNCTION
!!     Main Driver to compute the QM/MM Electrostatic Coupling
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE  qmmm_el_coupling(globenv,qs_env,qmmm_env,mm_particles,mm_cell,error)
    TYPE(global_environment_type), &
      pointer, OPTIONAL           :: globenv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cell_type), POINTER                 :: mm_cell

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_el_coupling', &
      routineP = moduleN//':'//routineN    

    INTEGER                                  :: handle, iw
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(section_vals_type), POINTER         :: input_section, interp_section, print_section
    TYPE(cp_logger_type), POINTER :: logger
    
    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    IF (.NOT.failure) THEN
       NULLIFY(ks_qmmm_env_loc, pw_pools, pw_env,input_section)
       CALL get_qs_env(qs_env=qs_env,&
                       pw_env=pw_env,&
                       para_env=para_env,&
                       input=input_section,&
                       ks_qmmm_env=ks_qmmm_env_loc,&
                       error=error)

       CALL pw_env_get(pw_env=pw_env, pw_pools=pw_pools, error=error)
       print_section => section_vals_get_subs_vals(input_section,"QMMM%PRINT",error=error)
       iw = cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",&
            extension=".qmmmLog",error=error)
       IF (iw>0) &
            WRITE(iw,'(T3,A)')"Information on the QM/MM Electrostatic Potential:"
       !
       ! Initializing vectors:
       !        Zeroing v_qmmm_rspace 
       CALL coeff_zero(coeff=ks_qmmm_env_loc%v_qmmm_rspace)
       !
       ! Electrostatic Interaction type...
       !
       SELECT CASE(qmmm_env%qmmm_coupl_type)
       CASE(do_qmmm_grid)
          IF (iw>0) &
               WRITE(iw,'(T3,A)')"QM/MM Coupling computed on GRID."
          ! Compute the MM electrostatic potential on GRID
          CALL qmmm_electrostatic_on_grid(qmmm_env=qmmm_env,&
                                          v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
                                          mm_particles=mm_particles,&
                                          mm_cell=mm_cell,&
                                          error=error)
       CASE(do_qmmm_gauss)
          IF (iw>0) &
               WRITE(iw,'(T3,A)')"QM/MM Coupling computed collocating the Gaussian Potential Functions."
          interp_section => section_vals_get_subs_vals(input_section,&
               "QMMM%INTERPOLATOR",error=error)
          CALL    qmmm_elec_with_gaussian(qmmm_env=qmmm_env,&
                                          v_qmmm=ks_qmmm_env_loc%v_qmmm_rspace,&
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
          IF (iw>0)&
               WRITE(iw,'(T3,A)')"No QM/MM Electrostatic coupling. Just Mechanical Coupling!"
       CASE DEFAULT
          IF (iw>0) &
               WRITE(iw,'(T3,A)')"Unknown Coupling..."
          CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

       CALL cp_print_key_finished_output(iw,logger,print_section,&
            "PROGRAM_RUN_INFO", error=error)

       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_section,&
            "POTENTIAL",error=error),cp_p_file)) then
          iw = cp_print_key_unit_nr(logger,print_section,"POTENTIAL",&
               extension=".qmmmLog",error=error)
          CALL rs_pw_to_cube(ks_qmmm_env_loc%v_qmmm_rspace%pw,iw,iw>0,&
               "QM/MM: MM ELECTROSTATIC POTENTIAL ",force_env_section=qs_env%input,&
               error=error)
          CALL cp_print_key_finished_output(iw,logger,print_section,&
               "POTENTIAL", error=error)
       END IF
       CALL timestop(0.0_dp,handle)
    END IF
  END SUBROUTINE qmmm_el_coupling

!!****f* qmmm_electrostatic_methods/qmmm_electrostatic_on_grid [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_on_grid
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction on Grid Points
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE  qmmm_electrostatic_on_grid(qmmm_env, v_qmmm, mm_particles, mm_cell, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(coeff_type), INTENT(INOUT)          :: v_qmmm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(cell_type), POINTER                 :: mm_cell
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_electrostatic_on_grid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,"I","",handle)

    CALL qmmm_electrostatic_on_grid_low(&
         v_qmmm,&
         v_qmmm%pw%cr3d,&
         qmmm_env%Potentials,&
         mm_particles,&
         qmmm_env%mm_atom_chrg,&
         qmmm_env%mm_el_pot_radius,&
         qmmm_env%mm_atom_index,&
         qmmm_env%num_mm_atoms,&
         mm_cell,&
         qmmm_env%dOmmOqm,&
         error)

    IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN

       CALL qmmm_electrostatic_on_grid_low(&
            v_qmmm,&
            v_qmmm%pw%cr3d,&
            qmmm_env%added_charges%potentials,&
            qmmm_env%added_charges%added_particles,&
            qmmm_env%added_charges%mm_atom_chrg,&
            qmmm_env%added_charges%mm_el_pot_radius,&
            qmmm_env%added_charges%mm_atom_index,&
            qmmm_env%added_charges%num_mm_atoms,&
            mm_cell,&
            qmmm_env%dOmmOqm,&
            error)       

    END IF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_electrostatic_on_grid

!!****f* qmmm_electrostatic_methods/qmmm_electrostatic_on_grid_low [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_on_grid_low
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction on Grid Points: Low Level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE qmmm_electrostatic_on_grid_low(v_qmmm, Varray, potentials,&
       mm_particles, mm_charges, mm_el_pot_radius, mm_atom_index, &
       num_mm_atoms, mm_cell, dOmmOqm, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(coeff_type), INTENT(INOUT)          :: v_qmmm
    REAL(KIND=dp), DIMENSION(*)              :: Varray
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: &
      routineN = 'qmmm_electrostatic_on_grid_low', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: handle, i, Imm, Imp, IndMM, &
                                                Ipot, ix, j, jj, k, my_j, &
                                                my_k, n1, n2, n3
    INTEGER, DIMENSION(2, 3)                 :: bo,gbo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dr1, dr2, dr3, dx, qt, r, rc, &
                                                rt1, rt2, rt3, rv1, rv2, rv3, &
                                                rx, rx2, rx3, rx4, rx5, xs1, &
                                                xs2, xs3
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pot0_2
    TYPE(qmmm_pot_type), POINTER             :: Pot

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_charges),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mm_atom_index),cp_failure_level,routineP,error,failure)

    bo = v_qmmm%pw%pw_grid%bounds_local  
    gbo = v_qmmm%pw%pw_grid%bounds
    n1 = v_qmmm%pw%pw_grid%npts_local (1)   
    n2 = v_qmmm%pw%pw_grid%npts_local (2)   
    n3 = v_qmmm%pw%pw_grid%npts_local (3)   
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
          r_pbc=pbc(mm_particles(IndMM)%r-dOmmOqm, mm_cell)+dOmmOqm
          rt1= r_pbc(1) 
          rt2= r_pbc(2) 
          rt3= r_pbc(3) 
          qt=mm_charges(Imm)
          !$omp parallel do default (none) private( i, j, k, Imm, IndMM,jj,my_k,my_j,&
          !$omp                            r, rc, qt, ix, rx, rx2, rx3, rx4, rx5,&
          !$omp                   dr1,dr2,dr3, xs1,xs2,xs3, rv1,rv2,rv3, rt1,rt2,rt3)&
          !$omp                            shared(Varray,pot0_2,dx,bo,n1,n2,n3,gbo)
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k=k-gbo(1,3)
             xs3 = REAL(my_k,dp)*dr3
             my_j=bo(1,2)-gbo(1,2)
             jj=(k-bo(1,3))*n1*n2
             xs2 = REAL(my_j,dp)*dr2
             rv3 = rt3 - xs3
             DO j = bo(1,2), bo(2,2)
                xs1 = (bo(1,1)-gbo(1,1))*dr1
                rv2 = rt2 - xs2
                DO i = bo(1,1), bo(2,1)
                   jj=jj+1
                   rv1 = rt1 - xs1
                   r   = SQRT(rv1*rv1 + rv2*rv2 + rv3*rv3)
                   ix  = FLOOR(r/dx)+1
                   rx  = (r-REAL(ix-1,dp)*dx)/dx
                   rx2 = rx*rx
                   rx3 = rx2*rx
                   rx4 = rx3*rx
                   rx5 = rx4*rx

                   Varray(jj) = Varray(jj) + &
                        (pot0_2(1,ix  )*(1_dp-10.0_dp*rx3+15.0_dp*rx4-6.0_dp*rx5) &
                        +pot0_2(2,ix  )*(rx-6.0_dp*rx3+8.0_dp*rx4-3.0_dp*rx5)      &
                        +pot0_2(3,ix  )*(0.5_dp*(rx2-rx5)+1.5_dp*(-rx3+rx4))     &
                        +pot0_2(1,ix+1)*(10.0_dp*rx3-15.0_dp*rx4+6.0_dp*rx5)       &
                        +pot0_2(2,ix+1)*(-4.0_dp*(rx3-rx5)+7.0_dp*(rx4-rx5))       &
                        +pot0_2(3,ix+1)*(0.5_dp*(rx3+rx5)-rx4))                   &
                        *(-1.0_dp) * & ! Sign of electrons
                        qt
                   xs1 = xs1 + dr1
                END DO
                xs2 = xs2 + dr2
             END DO
          END DO LoopOnGrid
       END DO LoopMM
    END DO MainLoopPot
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_electrostatic_on_grid_low

!!****f* qmmm_electrostatic_methods/qmmm_elec_with_gaussian [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_with gaussian
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction collocating the gaussian
!!     Electrostatic Potential
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian(qmmm_env, v_qmmm, mm_particles,&
       aug_pools, cube_info, para_env, eps_mm_rspace, pw_pools,&
       auxbas_grid, coarser_grid, interp_section, mm_cell, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(coeff_type), INTENT(INOUT)          :: v_qmmm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: aug_pools
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    INTEGER, INTENT(IN)                      :: auxbas_grid, coarser_grid
    TYPE(section_vals_type), POINTER         :: interp_section
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, igrid, ilevel, lb(3), &
                                                ngrids, ub(3), kind_interp, handle2
    LOGICAL                                  :: failure
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids

    ! Statements
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(mm_particles),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qmmm_env%mm_atom_chrg),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qmmm_env%mm_atom_index),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(aug_pools),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pools),cp_failure_level,routineP,error,failure)
    !Statements
    CALL timeset(routineN,"I","",handle)
    ngrids=SIZE(pw_pools)
    CALL pw_pools_create_pws(aug_pools,grids,use_data=REALDATA3D,in_space=REALSPACE,error=error)
    DO igrid=1,ngrids
       CALL pw_zero(grids(igrid)%pw)
    END DO
    
    CALL qmmm_elec_with_gaussian_low( grids, mm_particles,&
         qmmm_env%mm_atom_chrg, qmmm_env%mm_el_pot_radius, qmmm_env%mm_atom_index,  &
         qmmm_env%num_mm_atoms, cube_info, para_env, eps_mm_rspace, qmmm_env%pgfs,  &
         auxbas_grid, coarser_grid, pw_pools, aug_pools, qmmm_env%potentials,       &
         interp_section=interp_section, mm_cell=mm_cell, dOmmOqm=qmmm_env%dOmmOqm,  &
         periodic=qmmm_env%periodic, per_potentials=qmmm_env%per_potentials, error=error)

    IF (qmmm_env%move_mm_charges.OR.qmmm_env%add_mm_charges) THEN

       CALL qmmm_elec_with_gaussian_low( grids, qmmm_env%added_charges%added_particles,   &
         qmmm_env%added_charges%mm_atom_chrg, qmmm_env%added_charges%mm_el_pot_radius,    &
         qmmm_env%added_charges%mm_atom_index, qmmm_env%added_charges%num_mm_atoms,       & 
         cube_info, para_env, eps_mm_rspace, qmmm_env%added_charges%pgfs,auxbas_grid,     &
         coarser_grid, pw_pools, aug_pools, qmmm_env%added_charges%potentials,            &
         interp_section=interp_section, mm_cell=mm_cell, dOmmOqm=qmmm_env%dOmmOqm,        &
         periodic=qmmm_env%periodic, per_potentials=qmmm_env%per_potentials, error=error)

    END IF
    CALL section_vals_val_get(interp_section,"kind", i_val=kind_interp, error=error)
    ! RealSpace Interpolation
    SELECT CASE(kind_interp)
    CASE(linear_interp)
       ! Linear Interpolator
       DO Ilevel = coarser_grid, auxbas_grid+1, -1
          CALL pw_prolongate_l_NoPBC(grids(Ilevel  )%pw,&
                                   grids(Ilevel-1)%pw,&
                                   error)
       END DO
    CASE(spline3_nopbc_interp, spline3_pbc_interp)
       ! Spline Iterpolator
       call mp_sync(para_env%group)
       CALL timeset("spline3Int","I","",handle2)
       DO Ilevel = coarser_grid, auxbas_grid+1, -1
          CALL pw_prolongate_s3(grids(Ilevel  )%pw,&
                            grids(Ilevel-1)%pw,&
                            aug_pools(Ilevel)%pool,&
                            param_section=interp_section,&
                            error=error)
       END DO
       CALL timestop(0._dp,handle2)
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    lb = v_qmmm%pw%pw_grid%bounds_local(1,:)
    ub = v_qmmm%pw%pw_grid%bounds_local(2,:)

    v_qmmm%pw%cr3d = grids(auxbas_grid)%pw%cr3d (lb(1):ub(1),&
                                                 lb(2):ub(2),&
                                                 lb(3):ub(3) )

    CALL pw_pools_give_back_pws(aug_pools,grids,error=error)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_elec_with_gaussian

!!****f* qmmm_electrostatic_methods/qmmm_elec_with_gaussian_low [1.0] *
!!
!!   NAME
!!     qmmm_electrostatic_with gaussian_low
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction collocating the gaussian
!!     Electrostatic Potential - Low Level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian_low( tmp_grid, mm_particles,&
         mm_charges, mm_el_pot_radius, mm_atom_index, num_mm_atoms,&
         cube_info, para_env, eps_mm_rspace, pgfs, auxbas_grid, coarser_grid,&
         pw_pools, aug_pools, potentials, interp_section, mm_cell, dOmmOqm,&
         periodic, per_potentials, error)
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: tmp_grid
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges, mm_el_pot_radius
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    INTEGER, INTENT(IN)                      :: num_mm_atoms
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    INTEGER, INTENT(IN)                      :: auxbas_grid, coarser_grid
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools, aug_pools
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(qmmm_per_pot_p_type), DIMENSION(:), &
      POINTER                                :: per_potentials
    TYPE(section_vals_type), POINTER         :: interp_section
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL, INTENT(IN)                      :: periodic

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, Icg, IGauss, &
                                                ilevel, Imm, IndMM, &
                                                IRadTyp, ithread, LIndMM, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, height, W
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:), POINTER     :: xdat, ydat, zdat
    TYPE(qmmm_gaussian_type), POINTER        :: pgf
    TYPE(qmmm_per_pot_type), POINTER         :: per_pot
    TYPE(qmmm_pot_type), POINTER             :: pot
    INTEGER, DIMENSION(2,3)                  :: bo2
    INTEGER                                  :: n_rep_real(3)
!Statements

    NULLIFY(pgf, pot, per_pot, xdat, ydat, zdat)
    CALL timeset(routineN,"I","",handle)

    ithread            = 0
    Icg                = 0
    bo2 = tmp_grid(auxbas_grid)%pw%pw_grid%bounds
    ALLOCATE (xdat(bo2(1,1):bo2(2,1)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (ydat(bo2(1,2):bo2(2,2)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (zdat(bo2(1,3):bo2(2,3)), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

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
          Atoms: DO Imm = 1, SIZE(pot%mm_atom_index)
             Icg       =   Icg + 1
             LIndMM    =   pot%mm_atom_index(Imm)
             IndMM     =   mm_atom_index(LIndMM)
             ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm, mm_cell)+dOmmOqm
             W         =   mm_charges(LIndMM) * height
             IF (W == 0.0_dp) CYCLE Atoms
             CALL collocate_gf_rspace_NoPBC(zetp=alpha,&
                                            rp=ra,&
                                            scale=-1.0_dp,&
                                            W=W,&
                                            pwgrid=tmp_grid(ilevel)%pw,&
                                            cube_info=cube_info(ilevel),&
                                            eps_mm_rspace=eps_mm_rspace,&
                                            ithread=ithread, &
                                            xdat=xdat,&
                                            ydat=ydat,&
                                            zdat=zdat,&
                                            bo2=bo2,&
                                            n_rep_real=n_rep_real,&
                                            mm_cell=mm_cell)
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

    IF (periodic) THEN
       ! Long Range Part of the QM/MM Potential with Gaussians With Periodic Boundary Conditions
       CALL qmmm_elec_with_gaussian_LG    (pgfs=pgfs,&
                                           cgrid=tmp_grid(coarser_grid)%pw,&
                                           mm_charges=mm_charges,&
                                           mm_atom_index=mm_atom_index,&
                                           mm_particles=mm_particles,&
                                           mm_el_pot_radius=mm_el_pot_radius,&
                                           para_env=para_env,&
                                           coarser_grid_level=coarser_grid,&
                                           per_potentials=per_potentials,&
                                           mm_cell=mm_cell,&
                                           dOmmOqm=dOmmOqm,&
                                           error=error)       
    ELSE
       ! Long Range Part of the QM/MM Potential with Gaussians
       CALL qmmm_elec_with_gaussian_LR    (pgfs=pgfs,&
                                           grid=tmp_grid(coarser_grid)%pw,&
                                           mm_charges=mm_charges,&
                                           mm_atom_index=mm_atom_index,&
                                           mm_particles=mm_particles,&
                                           mm_el_pot_radius=mm_el_pot_radius,&
                                           para_env=para_env,&
                                           coarser_grid_level=coarser_grid,&
                                           potentials=potentials,&
                                           mm_cell=mm_cell,&
                                           dOmmOqm=dOmmOqm,&
                                           error=error)
    END IF
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qmmm_elec_with_gaussian_low

!!****f* qmmm_electrostatic_methods/qmmm_elec_with_gaussian_LG [1.0] *
!!
!!   NAME
!!     qmmm_elec_with gaussian_LG
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction collocating 
!!     (1/R - Sum_NG Gaussians) on the coarser grid level in G-SPACE
!!     Long Range QM/MM Electrostatic Potential with Gaussian - Low Level
!!     PERIODIC BOUNDARY CONDITION VERSION
!!
!!   NOTES
!!     This version includes the explicit code of Eval_Interp_Spl3_pbc
!!     in order to achieve better performance
!! 
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian_LG(pgfs, cgrid, mm_charges,&
       mm_atom_index, mm_particles, mm_el_pot_radius, para_env,&
       coarser_grid_level, per_potentials, mm_cell, dOmmOqm, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_type), POINTER                   :: cgrid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: coarser_grid_level
    TYPE(qmmm_per_pot_p_type), DIMENSION(:), &
      POINTER                                :: per_potentials
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_LG', &
      routineP = moduleN//':'//routineN    

    INTEGER                                  :: handle, Imm, IndMM, &
                                                IRadTyp, i, j, k, LIndMM, &
                                                my_j, my_k
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    REAL(KIND=dp)                            :: dr1c, dr2c, dr3c, qt,         &
                                                rt1, rt2, rt3, rv1, rv2, rv3, &
                                                xs1, xs2, xs3
    REAL(KIND=dp), DIMENSION(3)              :: ra, vec
    TYPE(qmmm_per_pot_type), POINTER         :: per_pot
    INTEGER       :: ivec(3), npts(3)
    INTEGER       :: ik1, ik2, ik3, ik4, ij1, ij2, ij3, ij4, ii1, ii2, ii3, ii4
    REAL(KIND=dp) :: dr1, dr2, dr3, xd1, xd2, xd3, val
    REAL(KIND=dp) :: aaa, aba, aca, ada, aab, abb, acb, adb, aac, abc, acc, adc, aad, abd, acd, add   
    REAL(KIND=dp) :: baa, bba, bca, bda, bab, bbb, bcb, bdb, bac, bbc, bcc, bdc, bad, bbd, bcd, bdd 
    REAL(KIND=dp) :: caa, cba, cca, cda, cab, cbb, ccb, cdb, cac, cbc, ccc, cdc, cad, cbd, ccd, cdd 
    REAL(KIND=dp) :: daa, dba, dca, dda, dab, dbb, dcb, ddb, dac, dbc, dcc, ddc, dad, dbd, dcd, ddd 
    REAL(KIND=dp) :: a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, t1, t2, t3, t4
    REAL(KIND=dp) :: e1, e2, e3, f1, f2, f3, g1, g2, g3, h1, h2, h3, s1, s2, s3, s4
    REAL(KIND=dp) :: p1, p2, p3, q1, q2, q3, r1, r2, r3, u1, u2, u3, v1, v2, v3, v4
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER    :: grid
    TYPE(pw_type), POINTER                      :: pw
! Statements

    CALL timeset(routineN,"I","",handle)
    NULLIFY(grid, pw)
    dr1c = cgrid%pw_grid%dr(1) 
    dr2c = cgrid%pw_grid%dr(2) 
    dr3c = cgrid%pw_grid%dr(3) 
    bo   = cgrid%pw_grid%bounds_local
    gbo  = cgrid%pw_grid%bounds

    Radius: DO IRadTyp = 1, SIZE(pgfs) 
       per_pot => per_potentials(IRadTyp)%pot
       pw => per_pot%TabLR
       npts = pw%pw_grid%npts
       dr1  = pw%pw_grid%dr(1)
       dr2  = pw%pw_grid%dr(2)
       dr3  = pw%pw_grid%dr(3)
       grid => pw%cr3d(:,:,:)
       Atoms: DO Imm = 1, SIZE(per_pot%mm_atom_index)
          LIndMM    =   per_pot%mm_atom_index(Imm)
          IndMM     =   mm_atom_index(LIndMM)
          ra(:)     =   pbc(mm_particles(IndMM)%r-dOmmOqm,mm_cell)+dOmmOqm
          qt        =   mm_charges(LIndMM)
          IF (qt == 0.0_dp) CYCLE Atoms
          rt1 = ra(1)
          rt2 = ra(2)
          rt3 = ra(3)
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k=k-gbo(1,3)
             xs3  = REAL(my_k,dp)*dr3c
             my_j=bo(1,2)-gbo(1,2)
             xs2 = REAL(my_j,dp)*dr2c
             rv3 = rt3 - xs3
             DO j =  bo(1,2), bo(2,2)
                xs1 = (bo(1,1)-gbo(1,1))*dr1c
                rv2 = rt2 - xs2
                DO i =  bo(1,1), bo(2,1)
                   rv1  = rt1 - xs1
                   !
                   ! Spline Interpolation 
                   !
                   vec = (/rv1,rv2,rv3/)
                   ivec = FLOOR(vec/pw%pw_grid%dr)
                   xd1  = (vec(1)/dr1)-REAL(ivec(1),kind=dp)
                   xd2  = (vec(2)/dr2)-REAL(ivec(2),kind=dp)
                   xd3  = (vec(3)/dr3)-REAL(ivec(3),kind=dp)    
                   ik1 = MODULO(ivec(3)-1,npts(3))+1
                   ik2 = MODULO(ivec(3)  ,npts(3))+1
                   ik3 = MODULO(ivec(3)+1,npts(3))+1
                   ik4 = MODULO(ivec(3)+2,npts(3))+1                   
                   ij1 = MODULO(ivec(2)-1,npts(2))+1
                   ij2 = MODULO(ivec(2)  ,npts(2))+1
                   ij3 = MODULO(ivec(2)+1,npts(2))+1
                   ij4 = MODULO(ivec(2)+2,npts(2))+1       
                   ii1 = MODULO(ivec(1)-1,npts(1))+1
                   ii2 = MODULO(ivec(1)  ,npts(1))+1
                   ii3 = MODULO(ivec(1)+1,npts(1))+1
                   ii4 = MODULO(ivec(1)+2,npts(1))+1
 
                   aaa   = grid(ii1,ij1,ik1)
                   baa   = grid(ii2,ij1,ik1)
                   caa   = grid(ii3,ij1,ik1)
                   daa   = grid(ii4,ij1,ik1)
                   aba   = grid(ii1,ij2,ik1)
                   bba   = grid(ii2,ij2,ik1)
                   cba   = grid(ii3,ij2,ik1)
                   dba   = grid(ii4,ij2,ik1)
                   aca   = grid(ii1,ij3,ik1)
                   bca   = grid(ii2,ij3,ik1)
                   cca   = grid(ii3,ij3,ik1)
                   dca   = grid(ii4,ij3,ik1)
                   ada   = grid(ii1,ij4,ik1)
                   bda   = grid(ii2,ij4,ik1)
                   cda   = grid(ii3,ij4,ik1)
                   dda   = grid(ii4,ij4,ik1)
                   aab   = grid(ii1,ij1,ik2)
                   bab   = grid(ii2,ij1,ik2)
                   cab   = grid(ii3,ij1,ik2)
                   dab   = grid(ii4,ij1,ik2)
                   abb   = grid(ii1,ij2,ik2)
                   bbb   = grid(ii2,ij2,ik2)
                   cbb   = grid(ii3,ij2,ik2)
                   dbb   = grid(ii4,ij2,ik2)
                   acb   = grid(ii1,ij3,ik2)
                   bcb   = grid(ii2,ij3,ik2)
                   ccb   = grid(ii3,ij3,ik2)
                   dcb   = grid(ii4,ij3,ik2)
                   adb   = grid(ii1,ij4,ik2)
                   bdb   = grid(ii2,ij4,ik2)
                   cdb   = grid(ii3,ij4,ik2)
                   ddb   = grid(ii4,ij4,ik2)
                   aac   = grid(ii1,ij1,ik3)
                   bac   = grid(ii2,ij1,ik3)
                   cac   = grid(ii3,ij1,ik3)
                   dac   = grid(ii4,ij1,ik3)
                   abc   = grid(ii1,ij2,ik3)
                   bbc   = grid(ii2,ij2,ik3)
                   cbc   = grid(ii3,ij2,ik3)
                   dbc   = grid(ii4,ij2,ik3)
                   acc   = grid(ii1,ij3,ik3)
                   bcc   = grid(ii2,ij3,ik3)
                   ccc   = grid(ii3,ij3,ik3)
                   dcc   = grid(ii4,ij3,ik3)
                   adc   = grid(ii1,ij4,ik3)
                   bdc   = grid(ii2,ij4,ik3)
                   cdc   = grid(ii3,ij4,ik3)
                   ddc   = grid(ii4,ij4,ik3)
                   aad   = grid(ii1,ij1,ik4)
                   bad   = grid(ii2,ij1,ik4)
                   cad   = grid(ii3,ij1,ik4)
                   dad   = grid(ii4,ij1,ik4)
                   abd   = grid(ii1,ij2,ik4)
                   bbd   = grid(ii2,ij2,ik4)
                   cbd   = grid(ii3,ij2,ik4)
                   dbd   = grid(ii4,ij2,ik4)
                   acd   = grid(ii1,ij3,ik4)
                   bcd   = grid(ii2,ij3,ik4)
                   ccd   = grid(ii3,ij3,ik4)
                   dcd   = grid(ii4,ij3,ik4)
                   add   = grid(ii1,ij4,ik4)
                   bdd   = grid(ii2,ij4,ik4)
                   cdd   = grid(ii3,ij4,ik4)
                   ddd   = grid(ii4,ij4,ik4)
                   
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
       
                   t1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*a1 + 12.0_dp*a2 - a3)
                   t2 = -22.0_dp/3.0_dp + 10.0_dp*b1 - 4.0_dp*b2 + 0.5_dp*b3
                   t3 =   2.0_dp/3.0_dp -  2.0_dp*c1 + 2.0_dp*c2 - 0.5_dp*c3
                   t4 =   1.0_dp/6.0_dp*d3
                   s1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*e1 + 12.0_dp*e2 - e3)
                   s2 = -22.0_dp/3.0_dp + 10.0_dp*f1 - 4.0_dp*f2 + 0.5_dp*f3
                   s3 =   2.0_dp/3.0_dp -  2.0_dp*g1 + 2.0_dp*g2 - 0.5_dp*g3
                   s4 =   1.0_dp/6.0_dp*h3
                   v1 =   1.0_dp/6.0_dp * (64.0_dp - 48.0_dp*p1 + 12.0_dp*p2 - p3)
                   v2 = -22.0_dp/3.0_dp + 10.0_dp*q1 - 4.0_dp*q2 + 0.5_dp*q3
                   v3 =   2.0_dp/3.0_dp -  2.0_dp*r1 + 2.0_dp*r2 - 0.5_dp*r3
                   v4 =   1.0_dp/6.0_dp*u3
       
                   val = (( aaa * t1 + baa * t2 + caa * t3 + daa * t4 ) * s1  +&
                          ( aba * t1 + bba * t2 + cba * t3 + dba * t4 ) * s2  +&
                          ( aca * t1 + bca * t2 + cca * t3 + dca * t4 ) * s3  +&
                          ( ada * t1 + bda * t2 + cda * t3 + dda * t4 ) * s4  ) * v1 +&
                         (( aab * t1 + bab * t2 + cab * t3 + dab * t4 ) * s1  +&
                          ( abb * t1 + bbb * t2 + cbb * t3 + dbb * t4 ) * s2  +&
                          ( acb * t1 + bcb * t2 + ccb * t3 + dcb * t4 ) * s3  +&
                          ( adb * t1 + bdb * t2 + cdb * t3 + ddb * t4 ) * s4  ) * v2 +&
                         (( aac * t1 + bac * t2 + cac * t3 + dac * t4 ) * s1  +&
                          ( abc * t1 + bbc * t2 + cbc * t3 + dbc * t4 ) * s2  +&
                          ( acc * t1 + bcc * t2 + ccc * t3 + dcc * t4 ) * s3  +&
                          ( adc * t1 + bdc * t2 + cdc * t3 + ddc * t4 ) * s4  ) * v3 +&
                         (( aad * t1 + bad * t2 + cad * t3 + dad * t4 ) * s1  +&
                          ( abd * t1 + bbd * t2 + cbd * t3 + dbd * t4 ) * s2  +&
                          ( acd * t1 + bcd * t2 + ccd * t3 + dcd * t4 ) * s3  +&
                          ( add * t1 + bdd * t2 + cdd * t3 + ddd * t4 ) * s4  ) * v4 

                   cgrid%cr3d(i,j,k) = cgrid%cr3d(i,j,k) - val * qt 
                   xs1 = xs1 + dr1c
                END DO
                xs2 = xs2 + dr2c
             END DO
          END DO LoopOnGrid
       END DO Atoms
    END DO Radius

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_elec_with_gaussian_LG


!!****f* qmmm_electrostatic_methods/qmmm_elec_with_gaussian_LR [1.0] *
!!
!!   NAME
!!     qmmm_elec_with gaussian_LR
!!
!!   FUNCTION
!!     Compute the QM/MM electrostatic Interaction collocating 
!!     (1/R - Sum_NG Gaussians) on the coarser grid level.
!!     Long Range QM/MM Electrostatic Potential with Gaussian - Low Level
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE qmmm_elec_with_gaussian_LR(pgfs, grid, mm_charges,&
       mm_atom_index, mm_particles, mm_el_pot_radius, para_env,&
       coarser_grid_level, potentials, mm_cell, dOmmOqm, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(pw_type), POINTER                   :: grid
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_charges
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: mm_particles
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: dOmmOqm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: coarser_grid_level
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_LR', &
      routineP = moduleN//':'//routineN    

    INTEGER                                  :: handle, i, Imm, IndMM, &
                                                IRadTyp, ix, j, k, LIndMM, &
                                                my_j, my_k, n1, n2, n3
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    REAL(KIND=dp)                            :: dr1, dr2, dr3, dx, qt, r, r2, &
                                                rt1, rt2, rt3, rv1, rv2, rv3, &
                                                rx, rx2, rx3, Term, xs1, xs2, &
                                                xs3
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pot0_2
    TYPE(qmmm_pot_type), POINTER             :: pot

! Statements

    CALL timeset(routineN,"I","",handle)

    n1  = grid%pw_grid%npts(1)
    n2  = grid%pw_grid%npts(2)
    n3  = grid%pw_grid%npts(3)
    dr1 = grid%pw_grid%dr(1) 
    dr2 = grid%pw_grid%dr(2) 
    dr3 = grid%pw_grid%dr(3) 
    bo  = grid%pw_grid%bounds_local
    gbo = grid%pw_grid%bounds

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
          LoopOnGrid: DO k = bo(1,3), bo(2,3)
             my_k=k-gbo(1,3)
             xs3  = REAL(my_k,dp)*dr3
             my_j=bo(1,2)-gbo(1,2)
             xs2 = REAL(my_j,dp)*dr2
             rv3 = rt3 - xs3
             DO j =  bo(1,2), bo(2,2)
                xs1 = (bo(1,1)-gbo(1,1))*dr1
                rv2 = rt2 - xs2
                DO i =  bo(1,1), bo(2,1)
                   rv1  = rt1 - xs1
                   r2   = rv1*rv1 + rv2*rv2 + rv3*rv3
                   r    = SQRT(r2)
                   ix  = FLOOR(r/dx)+1
                   rx  = (r-REAL(ix-1,dp)*dx)/dx
                   rx2 = rx*rx
                   rx3 = rx2*rx                   
                   Term = pot0_2(1,ix  )*(1.0_dp-3.0_dp*rx2+2.0_dp*rx3)  &
                         +pot0_2(2,ix  )*(rx-2.0_dp*rx2+rx3)             &
                         +pot0_2(1,ix+1)*(3.0_dp*rx2-2.0_dp*rx3)         &
                         +pot0_2(2,ix+1)*(-rx2+rx3)
                   grid%cr3d(i,j,k) = grid%cr3d(i,j,k) - Term * qt 
                   xs1 = xs1 + dr1
                END DO
                xs2 = xs2 + dr2
             END DO
          END DO LoopOnGrid
       END DO Atoms
    END DO Radius

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_elec_with_gaussian_LR


END MODULE qmmm_electrostatic_methods
