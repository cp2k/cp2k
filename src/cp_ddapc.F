!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE cp_ddapc
!!****h* cp2k/cp_ddapc *
!!
!!   NAME
!!    cp_ddapc
!!
!!   FUNCTION
!!     Density Derived atomic point charges from a QM calculation
!!     (see J. Chem. Phys. Vol. 103 pp. 7422-7428)
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
  USE cell_types,                      ONLY: cell_type
  USE coefficient_types,               ONLY: coeff_type
  USE cp_ddapc_forces,                 ONLY: ewald_ddapc_force,&
                                             reset_ch_pulay
  USE cp_ddapc_methods,                ONLY: build_A_matrix,&
                                             build_b_vector,&
                                             build_der_A_matrix,&
                                             build_der_b_vector
  USE cp_ddapc_types,                  ONLY: cp_ddapc_create,&
                                             cp_ddapc_release,&
                                             cp_ddapc_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_constants,                 ONLY: use_multipole
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
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
                                             pw_integral_ab,&
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
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc'
  PUBLIC :: get_ddapc, cp_decpl_ddapc
  
CONTAINS

!!****f* cp_ddapc/get_ddapc [1.0] *
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
  RECURSIVE SUBROUTINE get_ddapc(qs_env,qmmm_periodic,calc_force,fitcharge_section,&
       qout1,qout2,out_radii,dq_out,just_charge,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: qmmm_periodic, calc_force
    TYPE(section_vals_type), POINTER         :: fitcharge_section
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: qout1, qout2, out_radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: dq_out
    LOGICAL, INTENT(in), OPTIONAL            :: just_charge
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_ddapc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, handle2, handle3, i, &
                                                ii, iparticle, ispin, iw, j, &
                                                n_rep_val, ndim, num_gauss, &
                                                stat
    LOGICAL                                  :: failure, ionode, need_f, &
                                                only_charge
    REAL(KIND=dp)                            :: c1, c3, c4, ch_dens, core_ch, &
                                                elec_ch, gcut, pfact, rcmin, &
                                                Vol
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: bv, cv, dqv0, qtot, qv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Awrk, dbv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: dAm, dqv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: inp_radii, radii
    TYPE(cell_type), POINTER                 :: cell, super_cell
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(cp_ddapc_type), POINTER             :: cp_ddapc_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pool
    TYPE(pw_type), POINTER                   :: rho_tot_g
    TYPE(qs_charges_type), POINTER           :: qs_charges
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    need_f = .FALSE.
    only_charge = .FALSE.
    IF (PRESENT(just_charge)) only_charge=just_charge
    IF (PRESENT(calc_force)) need_f = calc_force    
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) iw = cp_logger_get_default_unit_nr(logger)
    NULLIFY(rho, rho_tot_g, rho_core, pw_env, radii, inp_radii, particle_set, &
         qs_charges, cp_ddapc_env)
    IF (.NOT.failure) THEN
       !
       ! Build the total QS density (electron+nuclei) in G-space
       !
       CALL get_qs_env(qs_env=qs_env,&
                       rho=rho,&
                       rho_core=rho_core,&
                       pw_env=pw_env,&
                       qs_charges=qs_charges,&
                       particle_set=particle_set,&
                       cell=cell,&
                       super_cell=super_cell,&
                       error=error)
       iw=cp_print_key_unit_nr(logger,fitcharge_section,&
            "/PROGRAM_RUN_INFO",".totCharge",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(T3,A,F15.9)')" Evaluating the charges to decouple QM images. " 
       END IF       
       CALL pw_env_get(pw_env=pw_env,&
                       auxbas_pw_pool=auxbas_pool,&
                       error=error)
       CALL pw_pool_create_pw(auxbas_pool,rho_tot_g,&
            in_space=RECIPROCALSPACE, use_data=COMPLEXDATA1D, error=error)
       CALL pw_transfer(rho_core%pw,rho_tot_g)       
       DO ispin=1,SIZE(rho%rho_g)
          CALL pw_sumup(rho%rho_g(ispin)%pw,rho_tot_g)
       END DO
       Vol = rho%rho_r(1)%pw%pw_grid%vol
       core_ch = qs_charges%total_rho_core_rspace
       elec_ch = SUM(rho%tot_rho_r)
       ch_dens = - (elec_ch + core_ch)/Vol
       !
       ! Get Input Parameters
       !
       CALL section_vals_val_get(fitcharge_section,"RADII",n_rep_val=n_rep_val,error=error)
       IF (n_rep_val /= 0) THEN
          CALL section_vals_val_get(fitcharge_section,"RADII",r_vals=inp_radii,error=error)
          num_gauss = SIZE(inp_radii)
          ALLOCATE(radii(num_gauss), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          radii = inp_radii
       ELSE
          CALL section_vals_val_get(fitcharge_section,"NUM_GAUSS",i_val=num_gauss,error=error)
          CALL section_vals_val_get(fitcharge_section,"MIN_RADIUS",r_val=rcmin,error=error)
          CALL section_vals_val_get(fitcharge_section,"PFACTOR",r_val=pfact,error=error)
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
       CALL section_vals_val_get(fitcharge_section,"GCUT",r_val=gcut,error=error)
       ! 
       ! Start with the linear system
       !
       CALL cp_ddapc_create(qs_env%cp_ddapc_env,qs_env%cp_ddapc_ewald, qmmm_periodic,&
            particle_set, radii, cell, super_cell,rho_tot_g, gcut, iw, Vol, &
            fitcharge_section, only_charge, error)
       cp_ddapc_env => qs_env%cp_ddapc_env 

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
       IF (iw>0) THEN
          WRITE(iw,'(/T3,A)')" Atomic point Charges on QM atoms:"
          DO i = 1, SIZE(particle_set)
             WRITE(iw,'(T6,A,I5,F12.6)')"Atom: ",i,qtot(i)
          END DO
          WRITE(iw,'(T5,A,5X,F12.6)')"Total: ",SUM(qtot)
          WRITE(iw,'(/)')
       END IF
       CALL timestop(0.0_dp,handle2)
       !
       ! If requested evaluate also the correction to derivatives due to Pulay Forces
       !
       IF (need_f.AND.(.NOT.only_charge)) THEN          
          CALL timeset(routineN//"-forces","I","",handle3)
          IF (iw>0) THEN
             WRITE(iw,'(T3,A)')" Evaluating Pulay Forces.."
          END IF
          ALLOCATE(dAm(ndim, ndim,3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dbv(ndim,3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dqv0(ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(dqv(ndim,SIZE(particle_set),3), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)     
          ALLOCATE(Awrk(ndim, ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          !
          DO iparticle = 1, SIZE(particle_set)
             dAm = 0.0_dp
             dbv = 0.0_dp
             CALL build_der_b_vector(dbv, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
                  particle_set, radii, rho_tot_g, gcut, iparticle, error);  dbv = dbv / Vol
             CALL mp_sum(dbv,rho_tot_g%pw_grid%para%group)
             IF (debug_this_module) THEN
                CALL debug_der_b_vector(dbv, qmmm_periodic, particle_set, radii, rho_tot_g,&
                     gcut, iparticle, Vol, fitcharge_section, qs_env, error)
                cp_ddapc_env => qs_env%cp_ddapc_env 
             END IF
             CALL build_der_A_matrix(dAm, cp_ddapc_env%gfunc, cp_ddapc_env%w,&
                  particle_set, radii, rho_tot_g, gcut, iparticle, error);  dAm = dAm /(Vol*Vol)
             CALL mp_sum(dAm,rho_tot_g%pw_grid%para%group)             
             IF (debug_this_module) THEN
                CALL debug_der_A_matrix(dAm, qmmm_periodic, particle_set, radii, rho_tot_g,&
                     gcut, iparticle, Vol, fitcharge_section, qs_env, error)
                cp_ddapc_env => qs_env%cp_ddapc_env 
             END IF
             DO j = 1, 3
                dqv0 = - MATMUL(cp_ddapc_env%AmI,MATMUL(dAm(:,:,j),-qv))
                Awrk = - MATMUL(cp_ddapc_env%AmI,MATMUL(dAm(:,:,j),cp_ddapc_env%AmI))
                c3   = DOT_PRODUCT(cv,MATMUL(Awrk,bv))+DOT_PRODUCT(cv,MATMUL(cp_ddapc_env%AmI,dbv(:,j))) ; c3 = c3 / cp_ddapc_env%c0
                c4   = DOT_PRODUCT(cv,MATMUL(Awrk,cv)) / cp_ddapc_env%c0
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
             CALL debug_charge(dqv, qs_env, qmmm_periodic, fitcharge_section, particle_set, radii, error)
             cp_ddapc_env => qs_env%cp_ddapc_env 
          END IF
          !
          DEALLOCATE(Awrk, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)          
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
       CALL cp_print_key_finished_output(iw,logger,fitcharge_section,&
            "/PROGRAM_RUN_INFO",error=error)
       CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_g,&
            accept_non_compatible=.TRUE.,error=error)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE get_ddapc


!!****f* cp_ddapc/cp_decpl_ddapc [1.0] *
!!
!!   NAME
!!     cp_decpl_ddapc
!!
!!   FUNCTION
!!     Routine to couple/decouple periodic images with the Bloechl scheme
!!
!!     The coupling/decoupling is obtaines evaluating terms E2 and E3 in
!!     J. Chem. Phys. Vol. 103 pp. 7422-7428.. The E2 terms is just a 
!!     Ewald summation, and for performance reason I'm writing a specific 
!!     driver instead of using and setting-up the environment of the already 
!!     available routines
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
  SUBROUTINE cp_decpl_ddapc(qs_env,rho_tot_gspace,energy,v_hartree_gspace,calculate_forces,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(coeff_type)                         :: rho_tot_gspace
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    TYPE(coeff_type)                         :: v_hartree_gspace
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_decpl_ddapc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, my_val, stat
    LOGICAL                                  :: failure, ionode, need_f, &
                                                qmmm_periodic
    REAL(KINd=dp)                            :: e_decpl, e_recpl
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges, radii
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: dq
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                multipole_section, &
                                                poisson_section, &
                                                qmmm_periodic_section

    CALL timeset(routineN,"I","",handle)
    failure = .FALSE.
    need_f = .FALSE. 
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) iw = cp_logger_get_default_unit_nr(logger)
    IF (PRESENT(calculate_forces)) need_f = calculate_forces
    IF (.NOT.failure) THEN
       NULLIFY(multipole_section, poisson_section, force_env_section, particle_set,&
               qmmm_periodic_section, charges, radii, dq)
       CALL get_qs_env(qs_env=qs_env,&
                       input=force_env_section,&
                       particle_set=particle_set,&
                       error=error)
       iw=cp_print_key_unit_nr(logger,force_env_section,"PROPERTIES%FIT_CHARGE/PROGRAM_RUN_INFO",&
            extension=".fitChargeLog",error=error)
       poisson_section => section_vals_get_subs_vals(force_env_section,"DFT%POISSON_DFT", error=error)
       qmmm_periodic_section   => section_vals_get_subs_vals(force_env_section,"QMMM%PERIODIC", error=error)
       CALL section_vals_get(qmmm_periodic_section,explicit=qmmm_periodic, error=error)
       qmmm_periodic = qmmm_periodic.AND.qs_env%qmmm
       IF (.NOT.qmmm_periodic) THEN
          CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",i_val=my_val,error=error)
          SELECT CASE (my_val)
          CASE (use_multipole)       
             multipole_section => section_vals_get_subs_vals(poisson_section,"MULTIPOLE", error=error)
          CASE DEFAULT
             CALL timestop(0.0_dp, handle)
             RETURN
          END SELECT
       ELSE
          multipole_section => section_vals_get_subs_vals(qmmm_periodic_section,"MULTIPOLE", error=error)
       END IF
       ! First we evaluate the charges at the corresponding SCF STEP
       IF (need_f) THEN
          CALL reset_ch_pulay(qs_env, error=error)
          CALL get_ddapc(qs_env, qmmm_periodic, need_f, multipole_section, qout1=charges, out_radii=radii,&
               dq_out=dq, error=error)
       ELSE
          CALL get_ddapc(qs_env, qmmm_periodic, need_f, multipole_section, qout1=charges, out_radii=radii,&
               error=error)
       END IF
       ! Evaluate the Ewald contribution to the decoupling/coupling E2 and E3
       IF (iw>0) THEN
          e_decpl = 0.5_dp*DOT_PRODUCT(charges,MATMUL(qs_env%cp_ddapc_env%Md,charges))
          WRITE(iw,'(T3,A,F15.9)')" Decoupling Energy: ", e_decpl
       END IF
       IF (qmmm_periodic) THEN
          IF (iw>0) THEN
             e_recpl = 0.5_dp*DOT_PRODUCT(charges,MATMUL(qs_env%cp_ddapc_env%Mr,charges))
             WRITE(iw,'(T3,A,F15.9)')" Recoupling Energy: ", e_recpl
          END IF
       END IF
       CALL modify_hartree_pot(v_hartree_gspace, rho_tot_gspace, multipole_section, particle_set,&
            qs_env%cp_ddapc_env%Mt, qs_env%cp_ddapc_env%AmI, radii, charges, error)
       ! Re-evaluate Hartree Energy after the Hartree potential has been modified
       energy=0.5_dp * pw_integral_ab ( rho_tot_gspace % pw, v_hartree_gspace % pw )
       IF (need_f) THEN
          CALL ewald_ddapc_force(qs_env, qs_env%cp_ddapc_ewald%coeff_qm,&
               .FALSE.,       multipole_section, particle_set,&
               radii, dq, charges, error=error)
          IF (qmmm_periodic) &
               CALL ewald_ddapc_force(qs_env, qs_env%cp_ddapc_ewald%coeff_mm,&
               qmmm_periodic, multipole_section, particle_set, radii, dq, charges,&
               error=error)
          CALL cp_ddapc_release(qs_env%cp_ddapc_env,error)
       END IF
       !
       DEALLOCATE(charges, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (ASSOCIATED(dq)) THEN
          DEALLOCATE(dq, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL cp_print_key_finished_output(iw,logger,force_env_section,&
            "PROPERTIES%FIT_CHARGE/PROGRAM_RUN_INFO",error=error)
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE cp_decpl_ddapc

!!****f* cp_ddapc/modify_hartree_pot [1.0] *
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
  SUBROUTINE modify_hartree_pot(v_hartree_gspace, rho_tot_gspace, multipole_section, particle_set,&
       M, AmI, radii, charges, error)
    TYPE(coeff_type)                         :: v_hartree_gspace, &
                                                rho_tot_gspace
    TYPE(section_vals_type), POINTER         :: multipole_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: M, AmI
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii, charges
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
       CALL section_vals_val_get(multipole_section,"GCUT",r_val=gcut,error=error)
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

!!****f* cp_ddapc/debug_der_b_vector [1.0] *
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
  SUBROUTINE debug_der_b_vector(dbv, qmmm_periodic, particle_set, radii,&
       rho_tot_g, gcut, iparticle, Vol, fitcharge_section, qs_env, error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: dbv
    LOGICAL, INTENT(in)                      :: qmmm_periodic
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(in)                      :: iparticle
    REAL(KIND=dp), INTENT(IN)                :: Vol
    TYPE(section_vals_type), POINTER         :: fitcharge_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

!!****f* cp_ddapc/debug_der_A_matrix [1.0] *
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
  SUBROUTINE debug_der_A_matrix(dAm, qmmm_periodic, particle_set, radii,&
       rho_tot_g, gcut, iparticle, Vol, fitcharge_section, qs_env, error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: dAm
    LOGICAL, INTENT(in)                      :: qmmm_periodic
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(in)                      :: iparticle
    REAL(KIND=dp), INTENT(IN)                :: Vol
    TYPE(section_vals_type), POINTER         :: fitcharge_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
  
!!****f* cp_ddapc/debug_charge [1.0] *
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
  SUBROUTINE debug_charge(dqv, qs_env, qmmm_periodic, fitcharge_section, particle_set, radii, error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: dqv
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in)                      :: qmmm_periodic
    TYPE(section_vals_type), POINTER         :: fitcharge_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
             CALL get_ddapc(qs_env,qmmm_periodic,.FALSE.,fitcharge_section,qout1=qtot1,error=error)
             particle_set(iparticle)%r(i) = rvec(i) - dx
             CALL get_ddapc(qs_env,qmmm_periodic,.FALSE.,fitcharge_section,qout1=qtot2,error=error)
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
  
END MODULE cp_ddapc
