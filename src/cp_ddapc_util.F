!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE cp_ddapc_util
!!****h* cp2k/cp_ddapc_util *
!!
!!   NAME
!!    cp_ddapc_util
!!
!!   FUNCTION
!!     Density Derived atomic point charges from a QM calculation
!!     (see Bloechl, J. Chem. Phys. Vol. 103 pp. 7422-7428)
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!   SOURCE
!****************************************************************************  
  USE atomic_charges,                  ONLY: print_atomic_charges
  USE cell_types,                      ONLY: cell_type
  USE coefficient_types,               ONLY: coeff_type
  USE cp_control_types,                ONLY: ddapc_restraint_type,&
                                             dft_control_type
  USE cp_ddapc_forces,                 ONLY: evaluate_restraint_functional
  USE cp_ddapc_methods,                ONLY: build_A_matrix,&
                                             build_b_vector,&
                                             build_der_A_matrix,&
                                             build_der_b_vector
  USE cp_ddapc_types,                  ONLY: cp_ddapc_create,&
                                             cp_ddapc_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_constants,                 ONLY: do_full_density,&
                                             do_spin_density
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_subtract,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type
  USE qs_charges_types,                ONLY: qs_charges_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc_util'
  PUBLIC :: get_ddapc,&
            restraint_functional_potential,&
            modify_hartree_pot,&
            cp_ddapc_init
  
CONTAINS

!!****f* cp_ddapc_types/cp_ddapc_init [1.0] *
!!
!!   NAME
!!     cp_ddapc_init
!!
!!   FUNCTION
!!     Initialize the cp_ddapc_environment
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
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE cp_ddapc_init(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_ddapc_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw, iw2, &
                                                n_rep_val, num_gauss, stat
    LOGICAL                                  :: allocate_ddapc_env, failure
    REAL(KIND=dp)                            :: gcut, pfact, rcmin, Vol
    REAL(KIND=dp), DIMENSION(:), POINTER     :: inp_radii, radii
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(coeff_type), POINTER                :: rho0_s_gs, rho_core
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pool
    TYPE(pw_type), POINTER                   :: rho_tot_g
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: density_fit_section

    CALL timeset(routineN,"I","",handle)
    failure         = .FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(dft_control, rho, rho_tot_g, rho_core, rho0_s_gs, pw_env,&
            radii, inp_radii, particle_set, qs_charges, cp_ddapc_env)
    allocate_ddapc_env = qs_env%cp_ddapc_ewald%do_solvation.OR.&
                         qs_env%cp_ddapc_ewald%do_qmmm_periodic.OR.&
                         qs_env%cp_ddapc_ewald%do_decoupling.OR.&
                         qs_env%cp_ddapc_ewald%do_restraint.OR.&
                         qs_env%cp_ddapc_ewald%do_property
    IF (.NOT.failure.AND.allocate_ddapc_env) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       dft_control=dft_control,&
                       rho=rho,&
                       rho_core=rho_core,&
                       rho0_s_gs=rho0_s_gs,&
                       pw_env=pw_env,&
                       qs_charges=qs_charges,&
                       particle_set=particle_set,&
                       cell=cell,&
                       super_cell=super_cell,&
                       error=error)
       density_fit_section => section_vals_get_subs_vals(qs_env%input,"DFT%DENSITY_FITTING",&
            error=error)
       iw=cp_print_key_unit_nr(logger,density_fit_section,&
            "PROGRAM_RUN_INFO",".FitCharge",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(/,A)')" Initializing the DDAPC Environment"
       END IF
       CALL pw_env_get(pw_env=pw_env,auxbas_pw_pool=auxbas_pool,error=error)
       CALL pw_pool_create_pw(auxbas_pool,rho_tot_g, in_space=RECIPROCALSPACE,&
            use_data=COMPLEXDATA1D, error=error)
       Vol = rho_tot_g%pw_grid%vol
       !
       ! Get Input Parameters
       !
       CALL section_vals_val_get(density_fit_section,"RADII",n_rep_val=n_rep_val,error=error)
       IF (n_rep_val /= 0) THEN
          CALL section_vals_val_get(density_fit_section,"RADII",r_vals=inp_radii,error=error)
          num_gauss = SIZE(inp_radii)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          radii = inp_radii
       ELSE
          CALL section_vals_val_get(density_fit_section,"NUM_GAUSS",i_val=num_gauss,error=error)
          CALL section_vals_val_get(density_fit_section,"MIN_RADIUS",r_val=rcmin,error=error)
          CALL section_vals_val_get(density_fit_section,"PFACTOR",r_val=pfact,error=error)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)          
          DO i = 1, num_gauss
             radii(i) = rcmin*pfact**(i-1)
          END DO
       END IF
       CALL section_vals_val_get(density_fit_section,"GCUT",r_val=gcut,error=error)
       ! Create DDAPC environment
       iw2=cp_print_key_unit_nr(logger,density_fit_section,&
            "PROGRAM_RUN_INFO/CONDITION_NUMBER",".FitCharge",error=error)
       ! Initialization of the cp_ddapc_env and of the cp_ddapc_ewald environment
       CALL cp_ddapc_create(qs_env%cp_ddapc_env,&
                            qs_env%cp_ddapc_ewald,&
                            particle_set,&
                            radii,&
                            cell,&
                            super_cell,&
                            rho_tot_g,&
                            gcut,&
                            iw,&
                            iw2,&
                            Vol, &
                            qs_env%input,&
                            error)
       CALL cp_print_key_finished_output(iw2,logger,density_fit_section,&
            "PROGRAM_RUN_INFO/CONDITION_NUMBER",error=error)
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)         
       CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_g,&
            accept_non_compatible=.TRUE.,error=error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_ddapc_init

!!****f* cp_ddapc_util/get_ddapc [1.0] *
!!
!!   NAME
!!     get_ddapc
!!
!!   FUNCTION
!!     Computes the Density Derived Atomic Point Charges 
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
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE get_ddapc(qs_env, calc_force, density_fit_section,&
       qout1, qout2, out_radii, dq_out, ext_rho_tot_g, Itype_of_density, &
       iwc, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    TYPE(section_vals_type), POINTER         :: density_fit_section
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: qout1, qout2, out_radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dq_out
    TYPE(pw_type), OPTIONAL, POINTER         :: ext_rho_tot_g
    CHARACTER(LEN=*), OPTIONAL               :: Itype_of_density
    INTEGER, INTENT(IN), OPTIONAL            :: iwc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_ddapc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: type_of_density
    INTEGER                                  :: handle, handle2, handle3, i, &
                                                ii, iparticle, ispin, iw, j, &
                                                myid, n_rep_val, ndim, &
                                                num_gauss, stat
    LOGICAL                                  :: failure, need_f
    REAL(KIND=dp)                            :: c1, c3, c4, ch_dens, gcut, &
                                                pfact, rcmin, Vol
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: bv, cv, dqv0, qtot, qv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dbv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: dAm, dqv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: inp_radii, radii
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(coeff_type), POINTER                :: rho0_s_gs, rho_core
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pool
    TYPE(pw_type), POINTER                   :: rho_tot_g
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,"I","",handle)
    failure         = .FALSE.
    need_f          = .FALSE.
    IF (PRESENT(calc_force))  need_f      = calc_force    
    logger => cp_error_get_logger(error)
    NULLIFY(dft_control, rho, rho_tot_g, rho_core, rho0_s_gs, pw_env,&
            radii, inp_radii, particle_set, qs_charges, cp_ddapc_env)
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       dft_control=dft_control,&
                       rho=rho,&
                       rho_core=rho_core,&
                       rho0_s_gs=rho0_s_gs,&
                       pw_env=pw_env,&
                       qs_charges=qs_charges,&
                       particle_set=particle_set,&
                       cell=cell,&
                       super_cell=super_cell,&
                       error=error)
       IF (PRESENT(iwc)) THEN 
          iw = iwc
       ELSE
          iw=cp_print_key_unit_nr(logger,density_fit_section,&
               "PROGRAM_RUN_INFO",".FitCharge",error=error)
       END IF
       CALL pw_env_get(pw_env=pw_env,&
                       auxbas_pw_pool=auxbas_pool,&
                       error=error)
       CALL pw_pool_create_pw(auxbas_pool,rho_tot_g, in_space=RECIPROCALSPACE,&
            use_data=COMPLEXDATA1D, error=error)
       IF (PRESENT(ext_rho_tot_g)) THEN
          ! If provided use the input density in g-space
          CALL pw_transfer(ext_rho_tot_g,rho_tot_g)
          type_of_density=Itype_of_density
       ELSE
          CALL section_vals_val_get(density_fit_section,"TYPE_OF_DENSITY",i_val=myid,error=error)
          SELECT CASE(myid)
          CASE(do_full_density)
             ! Otherwise build the total QS density (electron+nuclei) in G-space
             IF(dft_control%qs_control%gapw) THEN
                CALL pw_transfer(rho0_s_gs%pw,rho_tot_g)       
             ELSE
                CALL pw_transfer(rho_core%pw,rho_tot_g)       
             END IF
             DO ispin=1,SIZE(rho%rho_g)
                CALL pw_sumup(rho%rho_g(ispin)%pw,rho_tot_g)
             END DO
             type_of_density = "FULL DENSITY"
          CASE(do_spin_density)
             CALL pw_subtract(rho%rho_g(1)%pw,rho%rho_g(2)%pw,rho_tot_g)
             type_of_density = "SPIN DENSITY"
          END SELECT
       END IF
       Vol = rho%rho_r(1)%pw%pw_grid%vol
       ch_dens = 0.0_dp
       IF (rho_tot_g%pw_grid%have_g0) ch_dens = rho_tot_g%cc(1)
       CALL mp_sum(ch_dens,logger%para_env%group)
       !
       ! Get Input Parameters
       !
       CALL section_vals_val_get(density_fit_section,"RADII",n_rep_val=n_rep_val,error=error)
       IF (n_rep_val /= 0) THEN
          CALL section_vals_val_get(density_fit_section,"RADII",r_vals=inp_radii,error=error)
          num_gauss = SIZE(inp_radii)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          radii = inp_radii
       ELSE
          CALL section_vals_val_get(density_fit_section,"NUM_GAUSS",i_val=num_gauss,error=error)
          CALL section_vals_val_get(density_fit_section,"MIN_RADIUS",r_val=rcmin,error=error)
          CALL section_vals_val_get(density_fit_section,"PFACTOR",r_val=pfact,error=error)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)          
          DO i = 1, num_gauss
             radii(i) = rcmin*pfact**(i-1)
          END DO
       END IF
       IF (PRESENT(out_radii)) THEN
          IF (ASSOCIATED(out_radii)) THEN
             DEALLOCATE(out_radii, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          ALLOCATE(out_radii(SIZE(radii)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          out_radii = radii
       END IF
       CALL section_vals_val_get(density_fit_section,"GCUT",r_val=gcut,error=error)
       cp_ddapc_env => qs_env%cp_ddapc_env 
       ! 
       ! Start with the linear system
       !
       ndim = SIZE(particle_set)*SIZE(radii) 
       ALLOCATE(bv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(qv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       ALLOCATE(qtot(SIZE(particle_set)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)    
       ALLOCATE(cv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL timeset(routineN//"-charges","I","",handle2)
       bv  = 0.0_dp
       cv  = 1.0_dp/Vol
       CALL build_b_vector(bv, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
            particle_set, radii, rho_tot_g, gcut, error);  bv = bv / Vol
       CALL mp_sum(bv,rho_tot_g%pw_grid%para%group)
       c1 = DOT_PRODUCT(cv,MATMUL(cp_ddapc_env%AmI,bv))-ch_dens
       c1 = c1 / cp_ddapc_env%c0
       qv = - MATMUL(cp_ddapc_env%AmI,(bv-c1*cv))
       j  = 0
       qtot = 0.0_dp
       DO i = 1, ndim, num_gauss
          j = j+1
          DO ii = 1, num_gauss
             qtot(j) = qtot(j) + qv((i-1)+ii)
          END DO
       END DO
       IF (PRESENT(qout1)) THEN
          IF (ASSOCIATED(qout1)) THEN
             CPPostcondition(SIZE(qout1)==SIZE(qv),cp_failure_level,routineP,error,failure)
          ELSE
             ALLOCATE(qout1(SIZE(qv)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          qout1 = qv
       END IF
       IF (PRESENT(qout2)) THEN
          IF (ASSOCIATED(qout2)) THEN
             CPPostcondition(SIZE(qout2)==SIZE(qtot),cp_failure_level,routineP,error,failure)
          ELSE
             ALLOCATE(qout2(SIZE(qtot)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          qout2 = qtot
       END IF
       CALL print_atomic_charges(particle_set,iw,title=" DDAP "//&
            TRIM(type_of_density)//" charges:",atomic_charges=qtot)
       CALL timestop(0.0_dp,handle2)
       !
       ! If requested evaluate also the correction to derivatives due to Pulay Forces
       !
       IF (need_f) THEN 
          CALL timeset(routineN//"-forces","I","",handle3)
          IF (iw>0) THEN
             WRITE(iw,'(T3,A)')" Evaluating DDAPC atomic derivatives .."
          END IF
          ALLOCATE(dAm(ndim, ndim,3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dbv(ndim,3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dqv0(ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dqv(ndim,SIZE(particle_set),3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)     
          !
          DO iparticle = 1, SIZE(particle_set)
             dAm = 0.0_dp
             dbv = 0.0_dp
             CALL build_der_b_vector(dbv, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
                  particle_set, radii, rho_tot_g, gcut, iparticle, error);  dbv = dbv / Vol
             CALL mp_sum(dbv,rho_tot_g%pw_grid%para%group)
             IF (debug_this_module) THEN
                CALL debug_der_b_vector(dbv, particle_set, radii, rho_tot_g,&
                     gcut, iparticle, Vol, qs_env, error)
                cp_ddapc_env => qs_env%cp_ddapc_env 
             END IF
             CALL build_der_A_matrix(dAm, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
                  particle_set, radii, rho_tot_g, gcut, iparticle, error);  dAm = dAm /(Vol*Vol)
             CALL mp_sum(dAm,rho_tot_g%pw_grid%para%group)             
             IF (debug_this_module) THEN
                CALL debug_der_A_matrix(dAm, particle_set, radii, rho_tot_g,&
                     gcut, iparticle, Vol, qs_env, error)
                cp_ddapc_env => qs_env%cp_ddapc_env 
             END IF
             DO j = 1, 3
                dqv0 = - MATMUL(cp_ddapc_env%AmI,MATMUL(dAm(:,:,j),-qv))
                c3 = - DOT_PRODUCT(cv,MATMUL(cp_ddapc_env%AmI,MATMUL(dAm(:,:,j),MATMUL(cp_ddapc_env%AmI,bv)))) + &
                       DOT_PRODUCT(cv,MATMUL(cp_ddapc_env%AmI,dbv(:,j))) 
                c3   = c3 / cp_ddapc_env%c0
                c4 = - DOT_PRODUCT(cv,MATMUL(cp_ddapc_env%AmI,MATMUL(dAm(:,:,j),MATMUL(cp_ddapc_env%AmI,cv)))) / &
                     cp_ddapc_env%c0
                dqv(:,iparticle,j) = - ( dqv0 + MATMUL(cp_ddapc_env%AmI,dbv(:,j)) - &
                     c3 * MATMUL(cp_ddapc_env%AmI,cv) + c4 * c1 * MATMUL(cp_ddapc_env%AmI,cv))
             END DO
          END DO
          CPPostcondition(PRESENT(dq_out),cp_failure_level,routineP,error,failure)
          IF (.NOT.ASSOCIATED(dq_out)) THEN
             ALLOCATE(dq_out(SIZE(dqv,1),SIZE(dqv,2),SIZE(dqv,3)), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ELSE
             CPPostcondition(SIZE(dqv,1)==SIZE(dq_out,1),cp_failure_level,routineP,error,failure)
             CPPostcondition(SIZE(dqv,2)==SIZE(dq_out,2),cp_failure_level,routineP,error,failure)
             CPPostcondition(SIZE(dqv,3)==SIZE(dq_out,3),cp_failure_level,routineP,error,failure)
          END IF
          dq_out = dqv
          IF (debug_this_module) THEN
             CALL debug_charge(dqv, qs_env, density_fit_section, &
                  particle_set, radii, rho_tot_g, type_of_density, error)
             cp_ddapc_env => qs_env%cp_ddapc_env 
          END IF
          DEALLOCATE(dqv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dqv0, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dAm, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(dbv, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          CALL timestop(0.0_dp,handle3)
       END IF
       !
       ! End of charge fit
       !
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(bv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(cv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(qv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(qtot, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (.NOT.PRESENT(iwc)) THEN
          CALL cp_print_key_finished_output(iw,logger,density_fit_section,&
               "PROGRAM_RUN_INFO",error=error)
       END IF
       CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_g,&
            accept_non_compatible=.TRUE.,error=error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE get_ddapc

!!****f* cp_ddapc_util/restraint_functional_potential *
!!
!!   NAME
!!      restraint_functional_potential
!!
!!   FUNCTION
!!      modify hartree potential to handle restraints in DDAPC scheme
!!
!!   NOTES
!!      
!!
!!   INPUTS
!!        
!!
!!   MODIFICATION HISTORY
!!     02.2006  modified [Teo]
!!
!!   SOURCE
!!*** **********************************************************************
SUBROUTINE restraint_functional_potential(v_hartree_gspace,&
     density_fit_section, particle_set, AmI, radii, charges,&
     ddapc_restraint_control,energy_res, error)
    TYPE(coeff_type)                         :: v_hartree_gspace
    TYPE(section_vals_type), POINTER         :: density_fit_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: AmI
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii, charges
    TYPE(ddapc_restraint_type), &
      INTENT(INOUT)                          :: ddapc_restraint_control
    REAL(KIND=dp), INTENT(INOUT)             :: energy_res
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'restraint_functional_potential', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: g_corr, phase
    INTEGER                                  :: handle, I, idim, ig, igauss, &
                                                iparticle, n_gauss, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, fac, fac2, g2, gcut, &
                                                gcut2, gfunc, gvec(3), rc, &
                                                rc2, rvec(3), sfac, Vol, w
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cv, uv
    TYPE(pw_type), POINTER                   :: g_hartree

    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
    NULLIFY(g_hartree)
    IF (.NOT.failure) THEN
       n_gauss = SIZE(radii)
       ALLOCATE(cv(n_gauss*SIZE(particle_set)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(uv(n_gauss*SIZE(particle_set)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       uv      = 0.0_dp
       CALL evaluate_restraint_functional(ddapc_restraint_control, n_gauss, uv,&
            charges, energy_res, error)
       !
       CALL section_vals_val_get(density_fit_section,"GCUT",r_val=gcut,error=error)
       gcut2 = gcut*gcut
       g_hartree =>  v_hartree_gspace%pw
       Vol = g_hartree%pw_grid%vol
       cv   =  1.0_dp/Vol
       sfac = -1.0_dp/Vol
       fac  = DOT_PRODUCT(cv, MATMUL(AmI,cv))
       fac2 = DOT_PRODUCT(cv, MATMUL(AmI,uv))
       cv   = uv - cv * fac2/fac
       cv   = MATMUL(AmI,cv)
       IF (g_hartree%pw_grid%have_g0) g_hartree%cc(1) = g_hartree%cc(1) + sfac*fac2/fac
       DO ig = g_hartree % pw_grid % first_gne0, g_hartree % pw_grid % ngpts_cut_local
          g2 = g_hartree % pw_grid % gsq ( ig )
          w  = 4.0_dp * pi * (g2 - gcut2)**2.0_dp/(g2*gcut2)
          IF (g2 > gcut2) EXIT
          gvec  = g_hartree % pw_grid % g (:,ig)
          g_corr = 0.0_dp
          idim = 0
          DO iparticle = 1, SIZE(particle_set)
             DO igauss = 1, SIZE(radii)
                idim   = idim + 1
                rc     = radii(igauss)
                rc2    = rc*rc                
                rvec   = particle_set(iparticle)%r
                arg    = DOT_PRODUCT(gvec,rvec)
                phase  = CMPLX(COS(arg),-SIN(arg),KIND=dp)
                gfunc  = EXP(-g2*rc2/4.0_dp)
                g_corr = g_corr + gfunc * cv(idim) * phase 
             END DO
          END DO
          g_corr = g_corr * w
          g_hartree%cc(ig) = g_hartree%cc(ig) + sfac*g_corr/Vol
       END DO
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE restraint_functional_potential

!!****f* cp_ddapc_util/modify_hartree_pot [1.0] *
!!
!!   NAME
!!     modify_hartree_pot
!!
!!   FUNCTION
!!     Modify the Hartree potential
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
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE modify_hartree_pot(v_hartree_gspace, density_fit_section,&
       particle_set, M, AmI, radii, charges, error)
    TYPE(coeff_type)                         :: v_hartree_gspace
    TYPE(section_vals_type), POINTER         :: density_fit_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: M, AmI
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii, charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'modify_hartree_pot', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: g_corr, phase
    INTEGER                                  :: handle, idim, ig, igauss, &
                                                iparticle, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: arg, fac, fac2, g2, gcut, &
                                                gcut2, gfunc, gvec(3), rc, &
                                                rc2, rvec(3), sfac, Vol, w
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cv, uv
    TYPE(pw_type), POINTER                   :: g_hartree

    failure=.FALSE.
    CALL timeset(routineN,"I","",handle)
    NULLIFY(g_hartree)
    IF (.NOT.failure) THEN
       CALL section_vals_val_get(density_fit_section,"GCUT",r_val=gcut,error=error)
       gcut2 = gcut*gcut
       g_hartree =>  v_hartree_gspace%pw
       Vol = g_hartree%pw_grid%vol
       ALLOCATE(cv(SIZE(M,1)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(uv(SIZE(M,1)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       cv   = 1.0_dp/Vol
       uv   = MATMUL(M,charges)
       sfac = -1.0_dp/Vol
       fac  = DOT_PRODUCT(cv, MATMUL(AmI,cv))
       fac2 = DOT_PRODUCT(cv, MATMUL(AmI,uv))
       cv   = uv - cv * fac2/fac
       cv   = MATMUL(AmI,cv)
       IF (g_hartree%pw_grid%have_g0) g_hartree%cc(1) = g_hartree%cc(1) + sfac*fac2/fac
       DO ig = g_hartree % pw_grid % first_gne0, g_hartree % pw_grid % ngpts_cut_local
          g2 = g_hartree % pw_grid % gsq ( ig )
          w  = 4.0_dp * pi * (g2 - gcut2)**2.0_dp/(g2*gcut2)
          IF (g2 > gcut2) EXIT
          gvec  = g_hartree % pw_grid % g (:,ig)
          g_corr = 0.0_dp
          idim = 0
          DO iparticle = 1, SIZE(particle_set)
             DO igauss = 1, SIZE(radii)
                idim   = idim + 1
                rc     = radii(igauss)
                rc2    = rc*rc                
                rvec   = particle_set(iparticle)%r
                arg    = DOT_PRODUCT(gvec,rvec)
                phase  = CMPLX(COS(arg),-SIN(arg),KIND=dp)
                gfunc  = EXP(-g2*rc2/4.0_dp)
                g_corr = g_corr + gfunc * cv(idim) * phase 
             END DO
          END DO
          g_corr = g_corr * w
          g_hartree%cc(ig) = g_hartree%cc(ig) + sfac*g_corr/Vol
       END DO
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE modify_hartree_pot

!!****f* cp_ddapc_util/debug_der_b_vector [1.0] *
!!
!!   NAME
!!     debug_der_b_vector
!!
!!   FUNCTION
!!     To Debug the derivative of the B vector for the solution of the 
!!     linear system
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
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE debug_der_b_vector(dbv,  particle_set, radii,&
       rho_tot_g, gcut, iparticle, Vol, qs_env, error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: dbv
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(in)                      :: iparticle
    REAL(KIND=dp), INTENT(IN)                :: Vol
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'debug_der_b_vector', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, kk, ndim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, rvec(3), v0
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: bv1, bv2, ddbv
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env

    NULLIFY(cp_ddapc_env)
    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    IF (.NOT.failure) THEN
       dx   = 0.01_dp
       ndim = SIZE(particle_set)*SIZE(radii)
       ALLOCATE(bv1(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(bv2(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ddbv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       rvec = particle_set(iparticle)%r                
       cp_ddapc_env => qs_env%cp_ddapc_env
       DO i = 1, 3
          bv1 = 0.0_dp
          bv2 = 0.0_dp
          particle_set(iparticle)%r(i) = rvec(i) + dx
          CALL build_b_vector(bv1, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
               particle_set, radii, rho_tot_g, gcut, error);  bv1 = bv1 / Vol
          CALL mp_sum(bv1,rho_tot_g%pw_grid%para%group)
          particle_set(iparticle)%r(i) = rvec(i) - dx
          CALL build_b_vector(bv2, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
               particle_set, radii, rho_tot_g, gcut, error);  bv2 = bv2 / Vol
          CALL mp_sum(bv2,rho_tot_g%pw_grid%para%group)
          ddbv = (bv1 - bv2)/(2.0_dp*dx)
          DO kk = 1, SIZE(ddbv)
             IF (ddbv(kk).GT.1.0E-8_dp) THEN
                v0 = ABS(dbv(kk,i)-ddbv(kk))/ddbv(kk)*100.0_dp
                WRITE(*,*)"Error % on B ::",v0
                IF (v0.GT.0.1_dp) THEN
                   WRITE(*,'(A,2I5,2F15.9)')"ERROR IN DERIVATIVE OF B VECTOR, IPARTICLE, ICOORD:",iparticle,i,&
                        dbv(kk,i),ddbv(kk)
                   CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                END IF
             END IF
          END DO
          particle_set(iparticle)%r = rvec
       END DO
       DEALLOCATE(bv1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(bv2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ddbv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE debug_der_b_vector

!!****f* cp_ddapc_util/debug_der_A_matrix [1.0] *
!!
!!   NAME
!!     debug_der_A_matrix
!!
!!   FUNCTION
!!     To Debug the derivative of the A matrix for the solution of the 
!!     linear system
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
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE debug_der_A_matrix(dAm, particle_set, radii,&
       rho_tot_g, gcut, iparticle, Vol, qs_env, error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: dAm
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(in)                      :: iparticle
    REAL(KIND=dp), INTENT(IN)                :: Vol
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'debug_der_A_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, kk, ll, ndim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, rvec(3), v0
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Am1, Am2, ddAm
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env

    failure = .FALSE.
    NULLIFY(cp_ddapc_env)
    CALL timeset(routineN,"I","",handle)
    IF (.NOT.failure) THEN
       dx   = 0.01_dp
       ndim = SIZE(particle_set)*SIZE(radii)
       ALLOCATE(Am1(ndim,ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(Am2(ndim,ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ddAm(ndim,ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       rvec = particle_set(iparticle)%r                
       cp_ddapc_env => qs_env%cp_ddapc_env
       DO i = 1, 3
          Am1 = 0.0_dp
          Am2 = 0.0_dp
          particle_set(iparticle)%r(i) = rvec(i) + dx
          CALL build_A_matrix(Am1, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
               particle_set, radii, rho_tot_g, gcut, error);  Am1 = Am1 /(Vol*Vol)
          CALL mp_sum(Am1,rho_tot_g%pw_grid%para%group)
          particle_set(iparticle)%r(i) = rvec(i) - dx
          CALL build_A_matrix(Am2, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
               particle_set, radii, rho_tot_g, gcut, error);  Am2 = Am2 /(Vol*Vol)
          CALL mp_sum(Am2,rho_tot_g%pw_grid%para%group)
          ddAm = (Am1 - Am2)/(2.0_dp*dx)
          DO kk = 1, SIZE(ddAm,1)
             DO ll = 1, SIZE(ddAm,2)
                IF (ddAm(kk,ll).GT.1.0E-8_dp) THEN
                   v0 = ABS(dAm(kk,ll,i)-ddAm(kk,ll))/ddAm(kk,ll)*100.0_dp
                   WRITE(*,*)"Error % on A ::",v0,Am1(kk,ll), Am2(kk,ll),iparticle,i,kk,ll
                   IF (v0.GT.0.1_dp) THEN
                      WRITE(*,'(A,4I5,2F15.9)')"ERROR IN DERIVATIVE OF A MATRIX, IPARTICLE, ICOORD:",iparticle,i,kk,ll,&
                           dAm(kk,ll,i),ddAm(kk,ll)
                      CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                   END IF
                END IF
             END DO
          END DO
          particle_set(iparticle)%r = rvec
       END DO
       DEALLOCATE(Am1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(Am2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ddAm, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE debug_der_A_matrix
  
!!****f* cp_ddapc_util/debug_charge [1.0] *
!!
!!   NAME
!!     debug_charge
!!
!!   FUNCTION
!!     To Debug the fitted charges
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
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE debug_charge(dqv, qs_env, density_fit_section,&
       particle_set, radii, rho_tot_g, type_of_density, error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: dqv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: density_fit_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    CHARACTER(LEN=*)                         :: type_of_density
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'debug_charge', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iparticle, kk, &
                                                ndim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dx, rvec(3)
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ddqv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: qtot1, qtot2

    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    IF (.NOT.failure) THEN
       WRITE(*,*)"DEBUG_CHARGE_ROUTINE"
       ndim = SIZE(particle_set)*SIZE(radii)
       NULLIFY(qtot1, qtot2)
       ALLOCATE(qtot1(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)  
       ALLOCATE(qtot2(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(ddqv(ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !
       dx = 0.001_dp
       DO iparticle = 1, SIZE(particle_set)
          rvec = particle_set(iparticle)%r
          DO i = 1, 3
             particle_set(iparticle)%r(i) = rvec(i) + dx
             CALL get_ddapc(qs_env,.FALSE.,density_fit_section, qout1=qtot1,&
                  ext_rho_tot_g=rho_tot_g,Itype_of_density=type_of_density,error=error)
             particle_set(iparticle)%r(i) = rvec(i) - dx
             CALL get_ddapc(qs_env,.FALSE.,density_fit_section, qout1=qtot2,&
                  ext_rho_tot_g=rho_tot_g,Itype_of_density=type_of_density,error=error)
             ddqv = (qtot1 - qtot2)/(2.0_dp*dx)
             DO kk = 1, SIZE(qtot1)-1, SIZE(radii)
                IF (ANY(ddqv(kk:kk+2).GT.1.0E-8_dp)) THEN
                   WRITE(*,'(A,2F12.6,F12.2)')"Error :",SUM(dqv(kk:kk+2,iparticle,i)),SUM(ddqv(kk:kk+2)),&
                     ABS((SUM(ddqv(kk:kk+2))-SUM(dqv(kk:kk+2,iparticle,i)))/SUM(ddqv(kk:kk+2))*100.0_dp)
                END IF
             END DO
             particle_set(iparticle)%r = rvec
          END DO
       END DO
       !
       DEALLOCATE(qtot1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(qtot2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(ddqv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE debug_charge
  
END MODULE cp_ddapc_util
