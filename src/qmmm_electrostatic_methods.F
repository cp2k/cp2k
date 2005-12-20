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
  USE cp_output_handling, ONLY: cp_out_flags_should_output, cp_out_flags_type,&
       cp_p_file, cp_out_flags_unit_nr, cp_out_flags_finished_output,&
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
  USE pw_spline_utils,                 ONLY: pw_prolongate_s3, spl3_nopbct,&
       spl3_1d_transf_coeffs,spl3_1d_transf_border1, add_fine2coarse
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
  USE simulation_cell,                 ONLY: cell_type,&
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

    CHARACTER(LEN=80)                        :: FileName, Title
    INTEGER                                  :: handle, iw, unit_nr
    LOGICAL                                  :: failure
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env_loc
    TYPE(section_vals_type), POINTER         :: input_section, interp_section, print_section
    TYPE(cp_logger_type), POINTER :: logger
    TYPE(section_vals_type), POINTER :: print_pot
    TYPE(cp_out_flags_type), POINTER :: output_control
    
    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(ks_qmmm_env_loc, pw_pools, pw_env,input_section, output_control)
       CALL get_qs_env(qs_env=qs_env,&
                       pw_env=pw_env,&
                       para_env=para_env,&
                       input=input_section,&
                       ks_qmmm_env=ks_qmmm_env_loc,&
                       output_control=output_control,&
                       error=error)

       CALL pw_env_get(pw_env=pw_env, pw_pools=pw_pools, error=error)
       print_section => section_vals_get_subs_vals(input_section,"QMMM%PRINT",error=error)
       iw = cp_out_flags_unit_nr(output_control,"QMMM%PRINT%PROGRAM_RUN_INFO",&
            extension=".qmmmLog",error=error)
       logger => cp_error_get_logger(error)
       print_pot => section_vals_get_subs_vals(print_section,"POTENTIAL",error=error)
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

       CALL cp_out_flags_finished_output(iw,output_control,&
            "QMMM%PRINT%PROGRAM_RUN_INFO", error=error)

       IF (BTEST(cp_out_flags_should_output(output_control,&
            "QMMM%PRINT%POTENTIAL",error=error),cp_p_file)) then
          iw = cp_out_flags_unit_nr(output_control,"QMMM%PRINT%POTENTIAL",&
               extension=".qmmmLog",error=error)
          CALL rs_pw_to_cube(ks_qmmm_env_loc%v_qmmm_rspace%pw,iw,iw>0,&
               "QM/MM: MM ELECTROSTATIC POTENTIAL ")
          CALL cp_out_flags_finished_output(iw,output_control,&
               "QMMM%PRINT%POTENTIAL", error=error)
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
                                                ngrids, ub(3), kind_interp,gbo(2,3),&
                                                jgrid,kgrid,handle2
    LOGICAL                                  :: failure
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: grids
    TYPE(pw_type), POINTER   :: tmp_grid, tmp_grid2,tmp_grid3
    real(dp) :: maxdiff
    TYPE(cp_logger_type), pointer :: logger

! Statements

    logger => cp_error_get_logger(error)
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

    INTEGER                                  :: handle, i, Icg, IGauss, &
                                                ilevel, Imm, IndMM, &
                                                IRadTyp, ithread, LIndMM, stat
    LOGICAL                                  :: failure,ll
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
                                           grid=tmp_grid(coarser_grid)%pw,&
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
  SUBROUTINE qmmm_elec_with_gaussian_LG(pgfs, grid, mm_charges,&
       mm_atom_index, mm_particles, mm_el_pot_radius, para_env,&
       coarser_grid_level, per_potentials, mm_cell, dOmmOqm, error)
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
    TYPE(qmmm_per_pot_p_type), DIMENSION(:), &
      POINTER                                :: per_potentials
    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_elec_with_gaussian_LG', &
      routineP = moduleN//':'//routineN    

    INTEGER                                  :: handle, i, Imm, IndMM, &
                                                IRadTyp, ix, j, k, LIndMM, &
                                                my_j, my_k, n1, n2, n3, kg
    INTEGER, DIMENSION(2, 3)                 :: bo, gbo
    REAL(KIND=dp)                            :: dr1, dr2, dr3, qt, r, r2, &
                                                rt1, rt2, rt3, rv1, rv2, rv3, &
                                                rx, rx2, rx3, Term, xs1, xs2, &
                                                xs3
    REAL(KIND=dp), DIMENSION(3)              :: ra, vec
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Lg, gx, gy, gz
    TYPE(qmmm_per_pot_type), POINTER         :: per_pot

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
       per_pot => per_potentials(IRadTyp)%pot
       Lg => per_pot%Lg
       gx => per_pot%gx
       gy => per_pot%gy
       gz => per_pot%gz
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
             xs3  = REAL(my_k,dp)*dr3
             my_j=bo(1,2)-gbo(1,2)
             xs2 = REAL(my_j,dp)*dr2
             rv3 = rt3 - xs3
             DO j =  bo(1,2), bo(2,2)
                xs1 = (bo(1,1)-gbo(1,1))*dr1
                rv2 = rt2 - xs2
                DO i =  bo(1,1), bo(2,1)
                   rv1  = rt1 - xs1
                   ! Sum over k vectors
                   Term = 0.0_dp
                   DO kg = 1, SIZE(Lg)
                      vec  = (/REAL(gx(kg)),REAL(gy(kg)),REAL(gz(kg))/)
                      Term = Term + lg(kg)*COS(vec(1)*rv1+vec(2)*rv2+vec(3)*rv3)
                   END DO
                   grid%cr3d(i,j,k) = grid%cr3d(i,j,k) - Term * qt 
                   xs1 = xs1 + dr1
                END DO
                xs2 = xs2 + dr2
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
