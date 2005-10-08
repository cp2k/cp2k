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
  USE cp_output_handling,              ONLY: cp_should_output
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE global_types,                    ONLY: global_environment_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE mathconstants,                   ONLY: rootpi, pi
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
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pws,                             ONLY: pw_integral_ab
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALSPACE,&
                                             REALDATA3D,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type,&
                                             pw_zero
  USE particle_types,                  ONLY: particle_type 
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE coefficient_types,               ONLY: coeff_type
  USE mathlib,                         ONLY: invert_matrix
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get
  USE message_passing,                 ONLY: mp_sum
  USE input_constants,                 ONLY: use_multipole
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE qs_force_types,                  ONLY: qs_force_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE erf_fn, ONLY: erfc
#include "cp_common_uses.h"
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_ddapc'
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
  RECURSIVE SUBROUTINE get_ddapc(qs_env,calc_force,fitcharge_section,qout1,qout2,out_radii,dq_out,out_AmI,error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    TYPE(section_vals_type), POINTER         :: fitcharge_section
    REAL(KIND=dp), DIMENSION(:), POINTER, OPTIONAL   :: qout1, qout2, out_radii
    REAL(KINd=dp), DIMENSION(:,:), POINTER, OPTIONAL :: out_AmI
    REAL(KIND=dp), DIMENSION(:,:,:), POINTER, OPTIONAL :: dq_out
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'get_ddapc', &
      routineP = moduleN//':'//routineN  
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(pw_type), POINTER                   :: rho_tot_g
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(cp_logger_type), POINTER :: logger
    LOGICAL :: failure, need_f, ionode
    INTEGER :: iw, handle, ispin, n_rep_val, num_gauss, stat, i, ndim, ii, j, iparticle, kk, ll
    TYPE(pw_pool_type), POINTER              :: auxbas_pool
    REAL(KIND=dp), DIMENSION(:), POINTER     :: inp_radii, radii
    REAL(KIND=dp) :: rcmin, pfact, gcut, Vol, inv_error, c0, c1, c3, c4, dx, v0, ch_dens
    REAL(KIND=dp), DIMENSION(3) :: rvec
    TYPE(particle_type), DIMENSION(:), &
         POINTER                                 :: particle_set
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE     :: bv, cv, qv, qtot, ddqv, bv1, bv2, ddbv, dqv0
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE   :: Am, AmI, dbv, Awrk, Am1, Am2, ddAm
    REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE   :: dAm, dqv
    TYPE(section_vals_type), POINTER             :: print_gen
    REAL(KIND=dp), DIMENSION(:), POINTER         :: qtot1, qtot2

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    need_f = .FALSE.
    IF (PRESENT(calc_force)) need_f = calc_force    
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) iw = cp_logger_get_default_unit_nr(logger)
    NULLIFY(rho, rho_tot_g, rho_core, pw_env, radii, inp_radii, particle_set)
    IF (.NOT.failure) THEN
       print_gen => section_vals_get_subs_vals(fitcharge_section,"PRINT%PROGRAM_RUN_INFO",error=error)
       IF (cp_should_output(logger%iter_info,print_gen,error=error).AND.&
            logger%para_env%mepos==logger%para_env%source) THEN
          WRITE(iw,'(T3,A,F15.9)')" Evaluating the charges to couple QM images. " 
       END IF       
       !
       ! Build the total QS density (electron+nuclei) in G-space
       !
       CALL get_qs_env(qs_env=qs_env,&
                       rho=rho,&
                       rho_core=rho_core,&
                       pw_env=pw_env,&
                       particle_set=particle_set,&
                       error=error)

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
       ch_dens = - REAL(qs_env%dft_control%charge,KIND=dp)/Vol
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
          CPPostcondition(.NOT.ASSOCIATED(out_radii),cp_failure_level,routineP,error,failure) 
          ALLOCATE(out_radii(SIZE(radii)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          out_radii = radii          
       END IF
       CALL section_vals_val_get(fitcharge_section,"GCUT",r_val=gcut,error=error)
       ! 
       ! Start with the linear system
       !
       ndim = SIZE(particle_set)*SIZE(radii)
       ALLOCATE(Am(ndim, ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(AmI(ndim, ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(bv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(cv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(qv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       ALLOCATE(qtot(SIZE(particle_set)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       Am  = 0.0_dp
       AmI = 0.0_dp
       bv  = 0.0_dp
       cv  = 1.0_dp/Vol
       CALL build_b_vector(bv, particle_set, radii, rho_tot_g, gcut, error);  bv = bv / Vol
       CALL build_A_matrix(Am, particle_set, radii, rho_tot_g, gcut, error);  Am = Am /(Vol*Vol)       
       CALL mp_sum(bv,rho_tot_g%pw_grid%para%group)
       CALL mp_sum(Am,rho_tot_g%pw_grid%para%group)
       CALL invert_matrix(Am, AmI, inv_error, "N")
       IF (cp_should_output(logger%iter_info,print_gen,error=error).AND.&
            logger%para_env%mepos==logger%para_env%source) THEN
          WRITE(iw,'(T3,A,F15.9)')" Error inverting the A matrix: ", inv_error 
       END IF
       c0 = DOT_PRODUCT(cv,MATMUL(AmI,cv))
       c1 = DOT_PRODUCT(cv,MATMUL(AmI,bv))-ch_dens
       c1 = c1 / c0
       qv = - MATMUL(AmI,(bv-c1*cv))
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
       IF (cp_should_output(logger%iter_info,print_gen,error=error).AND.&
            logger%para_env%mepos==logger%para_env%source) THEN
          WRITE(iw,'(/T3,A)')" Charges on QM atoms used to couple/decouple the periodic images:"
          DO i = 1, SIZE(particle_set)
             WRITE(iw,'(T6,A,I5,F12.6)')"Atom: ",i,qtot(i)
          END DO
          WRITE(iw,'(T5,A,5X,F12.6)')"Total: ",SUM(qtot)
          WRITE(iw,'(/)')
       END IF
       !
       ! If requested evaluate also the correction to derivatives due to Pulay Forces
       !
       IF (need_f) THEN          
          IF (cp_should_output(logger%iter_info,print_gen,error=error).AND.&
               logger%para_env%mepos==logger%para_env%source) THEN
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
             CALL build_der_b_vector(dbv, particle_set, radii, rho_tot_g, gcut, iparticle, error);  dbv = dbv / Vol
             CALL mp_sum(dbv,rho_tot_g%pw_grid%para%group)
             IF (debug_this_module) THEN
                dx = 0.01_dp
                ALLOCATE(bv1(ndim), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(bv2(ndim), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(ddbv(ndim), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                rvec = particle_set(iparticle)%r                
                DO i = 1, 3
                   bv1 = 0.0_dp
                   bv2 = 0.0_dp
                   particle_set(iparticle)%r = rvec
                   particle_set(iparticle)%r(i) = rvec(i) + dx
                   CALL build_b_vector(bv1, particle_set, radii, rho_tot_g, gcut, error);  bv1 = bv1 / Vol
                   CALL mp_sum(bv1,rho_tot_g%pw_grid%para%group)
                   particle_set(iparticle)%r(i) = rvec(i) - dx
                   CALL build_b_vector(bv2, particle_set, radii, rho_tot_g, gcut, error);  bv2 = bv2 / Vol
                   CALL mp_sum(bv2,rho_tot_g%pw_grid%para%group)
                   ddbv = (bv1 - bv2)/(2.0_dp*dx)
                   DO kk = 1, SIZE(ddbv)
                      IF (ddbv(kk).GT.0.0_dp) THEN
                         v0 = ABS(dbv(kk,i)-ddbv(kk))/ddbv(kk)*100.0_dp
                         WRITE(*,*)"Error % on B ::",v0
                         IF (v0.GT.0.1_dp) THEN
                            WRITE(*,'(A,I5,I5)')"ERROR IN DERIVATIVE OF B VECTOR, IPARTICLE, ICOORD:",iparticle,i
                            CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                         END IF
                      END IF
                   END DO
                END DO                
                DEALLOCATE(bv1, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(bv2, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(ddbv, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             CALL build_der_A_matrix(dAm, particle_set, radii, rho_tot_g, gcut, iparticle, error);  dAm = dAm /(Vol*Vol)
             CALL mp_sum(dAm,rho_tot_g%pw_grid%para%group)             
             IF (debug_this_module) THEN
                ALLOCATE(Am1(ndim,ndim), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(Am2(ndim,ndim), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(ddAm(ndim,ndim), stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                rvec = particle_set(iparticle)%r                
                DO i = 1, 3
                   Am1 = 0.0_dp
                   Am2 = 0.0_dp
                   particle_set(iparticle)%r = rvec
                   particle_set(iparticle)%r(i) = rvec(i) + dx
                   CALL build_A_matrix(Am1, particle_set, radii, rho_tot_g, gcut, error);  Am1 = Am1 /(Vol*Vol)
                   CALL mp_sum(Am1,rho_tot_g%pw_grid%para%group)
                   particle_set(iparticle)%r(i) = rvec(i) - dx
                   CALL build_A_matrix(Am2, particle_set, radii, rho_tot_g, gcut, error);  Am2 = Am2 /(Vol*Vol)
                   CALL mp_sum(Am2,rho_tot_g%pw_grid%para%group)
                   ddAm = (Am1 - Am2)/(2.0_dp*dx)
                   DO kk = 1, SIZE(ddAm,1)
                      DO ll = 1, SIZE(ddAm,2)
                         IF (ddAm(kk,ll).GT.1.0E-8_dp) THEN
                            v0 = ABS(dAm(kk,ll,i)-ddAm(kk,ll))/ddAm(kk,ll)*100.0_dp
                            WRITE(*,*)"Error % on A ::",v0
                            IF (v0.GT.0.1_dp) THEN
                               WRITE(*,'(A,4I5,2F15.9)')"ERROR IN DERIVATIVE OF B VECTOR, IPARTICLE, ICOORD:",iparticle,i,kk,ll,&
                                    dAm(kk,ll,i),ddAm(kk,ll)
                               CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO                
                DEALLOCATE(Am1, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(Am2, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(ddAm, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
             DO j = 1, 3
                dqv0 = - MATMUL(AmI,MATMUL(dAm(:,:,j),-qv))
                Awrk = - MATMUL(AmI,MATMUL(dAm(:,:,j),AmI))
                c3   = DOT_PRODUCT(cv,MATMUL(Awrk,bv))+DOT_PRODUCT(cv,MATMUL(AmI,dbv(:,j))) ; c3 = c3 / c0
                c4   = DOT_PRODUCT(cv,MATMUL(Awrk,cv)) / c0
                dqv(:,iparticle,j) = - ( dqv0 + MATMUL(AmI,dbv(:,j)) - c3 * MATMUL(AmI,cv) + c4 * c1 * MATMUL(AmI,cv))
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
             NULLIFY(qtot1, qtot2)
             ALLOCATE(qtot1(SIZE(qv)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)  
             ALLOCATE(qtot2(SIZE(qv)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(ddqv(ndim),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             !
             dx = 0.001_dp
             DO iparticle = 1, SIZE(particle_set)
                rvec = particle_set(iparticle)%r
                DO i = 1, 3
                   particle_set(iparticle)%r = rvec
                   particle_set(iparticle)%r(i) = rvec(i) + dx
                   CALL get_ddapc(qs_env,.FALSE.,fitcharge_section,qout1=qtot1,error=error)
                   particle_set(iparticle)%r(i) = rvec(i) - dx
                   CALL get_ddapc(qs_env,.FALSE.,fitcharge_section,qout1=qtot2,error=error)
                   ddqv = (qtot1 - qtot2)/(2.0_dp*dx)
                   DO kk = 1, SIZE(qtot1)-1, num_gauss
                      WRITE(*,'(A,2F12.6,F12.2)')"Error :",SUM(dqv(kk:kk+2,iparticle,i)),SUM(ddqv(kk:kk+2)),&
                           ABS((SUM(ddqv(kk:kk+2))-SUM(dqv(kk:kk+2,iparticle,i)))/SUM(ddqv(kk:kk+2))*100.0_dp)
                   END DO
                END DO
                particle_set(iparticle)%r = rvec
             END DO
             !
             DEALLOCATE(qtot1, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(qtot2, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(ddqv, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
       END IF
       !
       ! End of charge fit
       !
       IF (PRESENT(out_AmI)) THEN
          IF (ASSOCIATED(out_AmI)) THEN
             CPPostcondition(SIZE(AmI,1)==SIZE(out_AmI,1),cp_failure_level,routineP,error,failure)
             CPPostcondition(SIZE(AmI,2)==SIZE(out_AmI,2),cp_failure_level,routineP,error,failure)
          ELSE
             ALLOCATE(out_AmI(SIZE(AmI,1),SIZE(AmI,2)), stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          out_AmI = AmI
       END IF
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(Am, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(AmI, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(bv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(cv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(qv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(qtot, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL pw_pool_give_back_pw(auxbas_pool,rho_tot_g,&
         accept_non_compatible=.TRUE.,error=error)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE get_ddapc

!!****f* cp_ddapc/build_b_vector [1.0] *
!!
!!   NAME
!!     build_b_vector
!!
!!   FUNCTION
!!     Computes the B vector for the solution of the linear system
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
  SUBROUTINE build_b_vector(bv, particle_set, radii, rho_tot_g, gcut, error)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: bv
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: radii
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'build_b_vector', &
      routineP = moduleN//':'//routineN      
    COMPLEX(KIND=dp) :: phase
    REAL(KIND=dp)    :: gfunc, rc, rc2, gcut2, g2, w, arg
    INTEGER          :: idim, iparticle, igauss, ig
    REAL(KIND=dp), DIMENSION(3) :: gvec, rvec

    gcut2 = gcut * gcut
    idim  = 0
    DO iparticle = 1, SIZE(particle_set)
       DO igauss = 1, SIZE(radii)
          rc   = radii(igauss)
          rc2  = rc*rc
          idim = idim + 1
          DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
             g2 = rho_tot_g % pw_grid % gsq ( ig )
             w  = 4.0_dp * pi * (g2 - gcut2)**2.0_dp/(g2*gcut2)
             IF (g2 > gcut2) EXIT
             gvec  = rho_tot_g % pw_grid % g (:,ig)
             rvec  = particle_set(iparticle)%r
             arg   = DOT_PRODUCT(gvec,rvec)
             phase = CMPLX(COS(arg),-SIN(arg))
             gfunc = EXP(-g2*rc2/4.0_dp)
             bv (idim) = bv(idim) + w * REAL(CONJG(rho_tot_g%cc(ig))*gfunc*phase,KIND=dp)
          END DO
       END DO
    END DO
        
  END SUBROUTINE build_b_vector
    
!!****f* cp_ddapc/build_A_matrix [1.0] *
!!
!!   NAME
!!     build_A_matrix
!!
!!   FUNCTION
!!     Computes the A matrix for the solution of the linear system
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
  SUBROUTINE build_A_matrix(Am, particle_set, radii, rho_tot_g, gcut, error)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: Am
    REAL(KIND=dp), DIMENSION(:), POINTER       :: radii
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'build_A_vector', &
      routineP = moduleN//':'//routineN      
    COMPLEX(KIND=dp) :: phase1, phase2
    REAL(KIND=dp)    :: gfunc1, gfunc2, rc, rc1, rc2, rc12, gcut2, g2, w, arg1, arg2
    INTEGER          :: idim1, idim2, iparticle1, iparticle2, igauss1, igauss2, ig, istart_g
    REAL(KIND=dp), DIMENSION(3) :: gvec, rvec1, rvec2

    gcut2 = gcut * gcut
    idim1  = 0
    DO iparticle1 = 1, SIZE(particle_set)
       DO igauss1 = 1, SIZE(radii)
          rc  = radii(igauss1)
          rc2 = rc*rc
          idim2 = idim1
          idim1 = idim1 + 1
          DO iparticle2 = iparticle1, SIZE(particle_set)
             istart_g = 1
             IF (iparticle2==iparticle1) istart_g = igauss1
             DO igauss2 = istart_g, SIZE(radii)
                rc1   = radii(igauss2)
                rc12  = rc1*rc1          
                idim2 = idim2 + 1
                DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
                   g2 = rho_tot_g % pw_grid % gsq ( ig )
                   w  = 4.0_dp * pi * (g2 -gcut2)**2/(g2*gcut2)
                   IF (g2 > gcut2) EXIT
                   gvec  = rho_tot_g % pw_grid % g (:,ig)
                   rvec1 = particle_set(iparticle1)%r
                   rvec2 = particle_set(iparticle2)%r
                   arg1   = DOT_PRODUCT(gvec,rvec1)
                   arg2   = DOT_PRODUCT(gvec,rvec2)
                   phase1 = CMPLX(COS(arg1),-SIN(arg1))
                   phase2 = CMPLX(COS(arg2),-SIN(arg2))
                   gfunc1 = EXP(-g2*rc2/4.0_dp)
                   gfunc2 = EXP(-g2*rc12/4.0_dp)
                   Am (idim1,idim2) = Am (idim1,idim2) + w * REAL(CONJG(gfunc1*phase1)*gfunc2*phase2,KIND=dp)
                END DO
                Am (idim2,idim1) = Am (idim1,idim2)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE build_A_matrix

!!****f* cp_ddapc/build_der_b_vector [1.0] *
!!
!!   NAME
!!     build_der_b_vector
!!
!!   FUNCTION
!!     Computes the derivative of B vector for the evaluation of the Pulay forces
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
  SUBROUTINE build_der_b_vector(dbv, particle_set, radii, rho_tot_g, gcut, iparticle0, error)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: dbv
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: radii
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(IN)                      :: iparticle0
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'build_der_b_vector', &
      routineP = moduleN//':'//routineN      
    COMPLEX(KIND=dp) :: dphase
    REAL(KIND=dp)    :: gfunc, rc, rc2, gcut2, g2, w, arg
    INTEGER          :: idim, iparticle, igauss, ig
    REAL(KIND=dp), DIMENSION(3) :: gvec, rvec

    gcut2 = gcut * gcut
    idim  = 0
    DO iparticle = 1, SIZE(particle_set)
       DO igauss = 1, SIZE(radii)
          rc   = radii(igauss)
          rc2  = rc*rc
          idim = idim + 1
          IF (iparticle /= iparticle0) CYCLE
          DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
             g2 = rho_tot_g % pw_grid % gsq ( ig )
             w  = 4.0_dp * pi * (g2 - gcut2)**2.0_dp/(g2*gcut2)
             IF (g2 > gcut2) EXIT
             gvec  = rho_tot_g % pw_grid % g (:,ig)
             rvec  = particle_set(iparticle)%r
             arg   = DOT_PRODUCT(gvec,rvec)
             dphase = - CMPLX(SIN(arg),COS(arg))
             gfunc = EXP(-g2*rc2/4.0_dp)
             dbv (idim,:) = dbv(idim,:) + w * REAL(CONJG(rho_tot_g%cc(ig))*gfunc*dphase,KIND=dp) * gvec
          END DO
       END DO
    END DO
        
  END SUBROUTINE build_der_b_vector
    
!!****f* cp_ddapc/build_der_A_matrix [1.0] *
!!
!!   NAME
!!     build_der_A_matrix
!!
!!   FUNCTION
!!     Computes the derivative of the A matrix for the evaluation of the 
!!     Pulay forces
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
  SUBROUTINE build_der_A_matrix(dAm, particle_set, radii, rho_tot_g, gcut, iparticle0, error)
    IMPLICIT NONE
    ! Arguments
    REAL(KIND=dp), DIMENSION(:,:,:), INTENT(INOUT) :: dAm
    REAL(KIND=dp), DIMENSION(:), POINTER       :: radii
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(IN)                      :: iparticle0
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'build_der_A_vector', &
      routineP = moduleN//':'//routineN      
    COMPLEX(KIND=dp) :: phase1, phase2, dphase1, dphase2
    REAL(KIND=dp)    :: gfunc1, gfunc2, rc, rc1, rc2, rc12, gcut2, g2, w, arg1, arg2
    INTEGER          :: idim1, idim2, iparticle1, iparticle2, igauss1, igauss2, ig, istart_g
    REAL(KIND=dp), DIMENSION(3) :: gvec, rvec1, rvec2

    gcut2 = gcut * gcut
    idim1  = 0
    DO iparticle1 = 1, SIZE(particle_set)
       DO igauss1 = 1, SIZE(radii)
          rc  = radii(igauss1)
          rc2 = rc*rc
          idim2 = idim1
          idim1 = idim1 + 1
          DO iparticle2 = iparticle1, SIZE(particle_set)
             istart_g = 1
             IF (iparticle2==iparticle1) istart_g = igauss1
             DO igauss2 = istart_g, SIZE(radii)
                rc1   = radii(igauss2)
                rc12  = rc1*rc1          
                idim2 = idim2 + 1
                IF ((iparticle1 /= iparticle0).AND.(iparticle2 /= iparticle0)) CYCLE
                DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
                   g2 = rho_tot_g % pw_grid % gsq ( ig )
                   w  = 4.0_dp * pi * (g2 -gcut2)**2/(g2*gcut2)
                   IF (g2 > gcut2) EXIT
                   gvec  = rho_tot_g % pw_grid % g (:,ig)
                   rvec1 = particle_set(iparticle1)%r
                   rvec2 = particle_set(iparticle2)%r
                   arg1   = DOT_PRODUCT(gvec,rvec1)
                   arg2   = DOT_PRODUCT(gvec,rvec2)
                   phase1 = CMPLX(COS(arg1),-SIN(arg1)); dphase1 = - CMPLX(SIN(arg1),COS(arg1))
                   phase2 = CMPLX(COS(arg2),-SIN(arg2)); dphase2 = - CMPLX(SIN(arg2),COS(arg2))
                   gfunc1 = EXP(-g2*rc2/4.0_dp)
                   gfunc2 = EXP(-g2*rc12/4.0_dp)
                   IF    ((iparticle1 == iparticle0).AND.(iparticle2==iparticle0))  THEN 
                      w = w * 2.0_dp
                      dAm (idim1,idim2,1:3) = dAm (idim1,idim2,1:3) + w * REAL(CONJG(gfunc1*dphase1)*gfunc2*phase2,KIND=dp) * gvec
                   ELSEIF (iparticle1 == iparticle0) THEN
                      dAm (idim1,idim2,1:3) = dAm (idim1,idim2,1:3) + w * REAL(CONJG(gfunc1*dphase1)*gfunc2*phase2,KIND=dp) * gvec
                   ELSEIF (iparticle2 == iparticle0) THEN
                      dAm (idim1,idim2,1:3) = dAm (idim1,idim2,1:3) + w * REAL(CONJG(gfunc1*phase1)*gfunc2*dphase2,KIND=dp) * gvec
                   END IF
                END DO
                dAm (idim2,idim1,1:3) = dAm (idim1,idim2,1:3)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE build_der_A_matrix


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
    IMPLICIT NONE
    ! Arguments
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(coeff_type)                         :: v_hartree_gspace, rho_tot_gspace
    LOGICAL, INTENT(IN), OPTIONAL            :: calculate_forces
    REAL(KIND=dp), INTENT(INOUT)             :: energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_decpl_ddapc', &
      routineP = moduleN//':'//routineN      
    TYPE(section_vals_type), POINTER         :: multipole_section,&
                                                force_env_section,&
                                                poisson_section,&
                                                qmmm_periodic_section,&
                                                print_gen
    TYPE(particle_type), DIMENSION(:), &
         POINTER                                 :: particle_set
    LOGICAL :: failure, qmmm_periodic, need_f, ionode
    INTEGER :: stat, my_val, iw
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: radii
    REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: M, AmI
    REAL(KIND=dp), POINTER, DIMENSION(:,:,:) :: dq
    TYPE(cp_logger_type), POINTER :: logger
    REAL(KINd=dp) :: e_decpl, e_recpl

    failure = .FALSE.
    need_f = .FALSE.
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) iw = cp_logger_get_default_unit_nr(logger)
    IF (PRESENT(calculate_forces)) need_f = calculate_forces
    IF (.NOT.failure) THEN
       NULLIFY(multipole_section, poisson_section, force_env_section, particle_set,&
               qmmm_periodic_section, charges, radii, M, AmI, dq)
       CALL get_qs_env(qs_env=qs_env,&
                       input=force_env_section,&
                       particle_set=particle_set,&
                       error=error)
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
             RETURN
          END SELECT
       ELSE
          multipole_section => section_vals_get_subs_vals(qmmm_periodic_section,"MULTIPOLE", error=error)
       END IF
       print_gen => section_vals_get_subs_vals(multipole_section,"PRINT%PROGRAM_RUN_INFO",error=error)
       ! First we evaluate the charges at the corresponding SCF STEP
       IF (need_f) THEN
          CALL get_ddapc(qs_env, need_f, multipole_section, qout1=charges, out_radii=radii, out_AmI=AmI,&
               dq_out=dq, error=error)
          CALL reset_ch_pulay(qs_env, error=error)
       ELSE
          CALL get_ddapc(qs_env, need_f, multipole_section, qout1=charges, out_radii=radii, out_AmI=AmI,&
               error=error)
       END IF
       ! Evaluate the Ewald contribution to the decoupling/coupling E2 and E3
       ALLOCATE(M(SIZE(charges),SIZE(charges)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       M       = 0.0_dp
       e_decpl = 0.0_dp
       IF (qmmm_periodic) THEN
          CALL ewald_ddapc_pot(qs_env, .FALSE., multipole_section, particle_set,&
               M, radii, need_f, dq, charges, error=error)
          IF (cp_should_output(logger%iter_info,print_gen,error=error).AND.&
               logger%para_env%mepos==logger%para_env%source) THEN
             e_decpl = 0.5_dp*DOT_PRODUCT(charges,MATMUL(M,charges))
             WRITE(iw,'(T3,A,F15.9)')" Decoupling Energy: ", e_decpl
          END IF
       END IF
       CALL ewald_ddapc_pot(qs_env, qmmm_periodic, multipole_section, particle_set,&
            M, radii, need_f, dq, charges, error=error)
       IF (cp_should_output(logger%iter_info,print_gen,error=error).AND.&
            logger%para_env%mepos==logger%para_env%source) THEN
          e_recpl = 0.5_dp*DOT_PRODUCT(charges,MATMUL(M,charges))
          IF (qmmm_periodic) THEN
             WRITE(iw,'(T3,A,F15.9)')" Recoupling Energy: ", e_recpl
          ELSE
             WRITE(iw,'(T3,A,F15.9)')" Decoupling Energy: ", e_recpl
          END IF
       END IF
       CALL modify_hartree_pot(v_hartree_gspace, rho_tot_gspace, multipole_section, particle_set, M, &
            AmI, radii, charges, error)
       ! Re-evaluate Hartree Energy after the Hartree potential has been modified
       energy=0.5_dp * pw_integral_ab ( rho_tot_gspace % pw, v_hartree_gspace % pw )
       !
       DEALLOCATE(AmI, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(charges, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       DEALLOCATE(radii, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(M, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (ASSOCIATED(dq)) THEN
          DEALLOCATE(dq, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE cp_decpl_ddapc

!!****f* cp_ddapc/ewald_ddapc_pot [1.0] *
!!
!!   NAME
!!     ewald_ddapc_pot
!!
!!   FUNCTION
!!     Evaluates the Ewald term E2 and E3 energy term for the decoupling/coupling
!!     of periodic images
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
  RECURSIVE SUBROUTINE ewald_ddapc_pot(qs_env, qmmm_periodic, multipole_section, &
       particle_set, M, radii, need_f, dq, charges, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: multipole_section    
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    LOGICAL, INTENT(IN)                      :: qmmm_periodic
    REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: M
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: radii
    LOGICAL, INTENT(IN)                      :: need_f
    REAL(KIND=dp), POINTER, DIMENSION(:,:,:), OPTIONAL :: dq
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: charges
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                            :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_ddapc_pot', &
         routineP = moduleN//':'//routineN      
    LOGICAL :: failure
    INTEGER :: nmax1, nmax2, nmax3, k1, k2, k3, n_rep, rmax1, rmax2, rmax3, ich1, ich2, stat
    INTEGER :: iparticle1, iparticle2, idim1, idim2, igauss1, igauss2, istart_g,r1,r2,r3, sfact
    REAL(KIND=dp) :: eps, rcut, tol, alpha, tol1, fac, rc, rc2, rc1, rc12, galpha, gsq, dfact
    REAL(KIND=dp) :: g_ewald, gsqi, fac3, r_ewald, r, fs, eps0, ew_neut, fac2, q1, q2, t1, t2
    REAL(KIND=dp), DIMENSION(3) :: rvec, gvec, fvec, ra, rb, drvec
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: d_el
    TYPE(cell_type), POINTER                 :: cell
    TYPE(unit_convert_type), POINTER         :: units

    failure = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(ASSOCIATED(M),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(radii),cp_failure_level,routineP,error,failure)
       CALL get_cp2k_units ( units )
       sfact =  1.0_dp
       IF (qmmm_periodic) THEN 
          sfact = -1.0_dp
          CALL get_qs_env(qs_env, super_cell=cell)
       ELSE
          CALL get_qs_env(qs_env, cell=cell)
       END IF
       rcut = MIN(cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3))/2.0_dp
       CALL section_vals_val_get(multipole_section,"RCUT",n_rep_val=n_rep,error=error)
       IF (n_rep==1) CALL section_vals_val_get(multipole_section,"RCUT",r_val=rcut,error=error)
       CALL section_vals_val_get(multipole_section,"EWALD_PRECISION",r_val=eps,error=error)
       !
       ! Setting-up parameters for Ewald summation
       !
       eps0   = units%eps0
       eps    = MIN(ABS(eps),0.5_dp)
       tol    = SQRT(ABS(LOG(eps*rcut)))
       alpha  = SQRT(ABS(LOG(eps*rcut*tol)))/rcut
       galpha = 1.0_dp/(4.0_dp*alpha*alpha)
       tol1   = SQRT(-LOG(eps*rcut*(2.0_dp*tol*alpha)**2))
       nmax1  = NINT(0.25_dp + cell%hmat(1,1)*alpha*tol1/pi)
       nmax2  = NINT(0.25_dp + cell%hmat(2,2)*alpha*tol1/pi)
       nmax3  = NINT(0.25_dp + cell%hmat(3,3)*alpha*tol1/pi)
       
       rmax1  = NINT(rcut/cell%hmat(1,1))
       rmax2  = NINT(rcut/cell%hmat(2,2))
       rmax3  = NINT(rcut/cell%hmat(3,3))
       IF (need_f) THEN 
          ALLOCATE(d_el(3,SIZE(particle_set)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          d_el = 0.0_dp
       END IF       
       fac    = 1.d0/cell%deth
       fac2   = 1.0_dp / ( 4.0_dp * pi * eps0 )
       fac3   = fac/8.0_dp
       fvec   = 2.0_dp * pi / (/cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3)/)
       idim1  = 0
       ich1   = 0
       DO iparticle1 = 1, SIZE(particle_set)
          DO igauss1 = 1, SIZE(radii)
             ich2  = ich1
             ich1  = ich1+1
             q1    = charges(ich1)
             rc    = radii(igauss1)
             rc2   = rc*rc
             idim2 = idim1
             idim1 = idim1 + 1
             DO iparticle2 = iparticle1, SIZE(particle_set)
                istart_g = 1
                IF (iparticle2==iparticle1) istart_g = igauss1
                DO igauss2 = istart_g, SIZE(radii)
                   ich2  = ich2+1
                   q2    = charges(ich2)
                   rc1   = radii(igauss2)
                   rc12  = rc1*rc1          
                   idim2 = idim2 + 1
                   rvec   = particle_set(iparticle1)%r - particle_set(iparticle2)%r
                   !
                   ! Real-Space Contribution
                   !
                   r_ewald = 0.0_dp
                   dfact = 1.0_dp
                   IF (iparticle1 == iparticle2) dfact = 0.5_dp
                   IF (iparticle1 /= iparticle2) THEN
                      ra = pbc(rvec,cell)
                      r  = SQRT(DOT_PRODUCT(ra,ra))
                      IF (r<=rcut) THEN 
                         t1 = erfc(alpha*r) / r
                         r_ewald = t1
                         IF (need_f) THEN
                            drvec = ra / r * q1 * q2 * fac2 * sfact
                            t2 = -2.0_dp * alpha * EXP(-alpha*alpha*r*r) / (r*rootpi) - t1/r
                            d_el(1:3,iparticle1) = d_el(1:3,iparticle1) - t2 * drvec
                            d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + t2 * drvec
                         END IF
                      END IF
                   END IF
                   DO r1 = -rmax1, rmax1
                      DO r2 = -rmax2, rmax2
                         DO r3 = -rmax3, rmax3
                            IF ((r1==0).AND.(r2==0).AND.(r3==0)) CYCLE
                            rb = (/REAL(r1,KIND=dp),REAL(r2,KIND=dp),REAL(r3,KIND=dp)/)
                            ra = pbc(rvec,cell) + MATMUL(cell%hmat,rb)
                            r  = SQRT(DOT_PRODUCT(ra,ra))
                            IF (r<=rcut) THEN 
                               t1 = erfc(alpha*r) / r
                               r_ewald = r_ewald + t1                           
                               IF (need_f) THEN
                                  drvec = ra / r * q1 * q2 * fac2 * dfact * sfact
                                  t2 = -2.0_dp * alpha * EXP(-alpha*alpha*r*r) / (r*rootpi) - t1/r
                                  d_el(1:3,iparticle1) = d_el(1:3,iparticle1) - t2 * drvec
                                  d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + t2 * drvec
                               END IF
                            END IF
                         END DO
                      END DO
                   END DO
                   !
                   ! G-space Contribution
                   !
                   g_ewald = 0.0_dp
                   DO k1 = 0, nmax1
                      DO k2 = -nmax2, nmax2
                         DO k3 = -nmax3, nmax3
                            IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE
                            fs = 2.0_dp; IF (k1==0) fs = 1.0_dp
                            gvec = fvec * (/REAL(k1,KIND=dp),REAL(k2,KIND=dp),REAL(k3,KIND=dp)/)
                            gsq  = DOT_PRODUCT(gvec,gvec)
                            gsqi = fs/gsq
                            t1   = fac * gsqi * EXP(-galpha*gsq)
                            g_ewald = g_ewald + t1 * COS(DOT_PRODUCT(gvec,rvec))
                            IF (need_f) THEN
                               t2 = -SIN(DOT_PRODUCT(gvec,rvec)) * t1 * q1 * q2 / eps0 * dfact * sfact
                               d_el(1:3,iparticle1) = d_el(1:3,iparticle1) - t2 * gvec
                               d_el(1:3,iparticle2) = d_el(1:3,iparticle2) + t2 * gvec
                            END IF
                         END DO
                      END DO
                   END DO
                   !
                   ! G-EWALD, R-EWALD and BACKGROUND CHARGE
                   !
                   g_ewald = r_ewald*fac2 + g_ewald/eps0
                   ew_neut = - fac / ( 4.0_dp * eps0 * alpha ** 2 )
                   g_ewald = g_ewald + ew_neut
                   !
                   ! Self Contribution
                   !
                   IF (iparticle1 == iparticle2) THEN
                      g_ewald = g_ewald - alpha/(2.0_dp*eps0*pi*SQRT(pi))
                   END IF
                   !
                   IF (iparticle1/=iparticle2) THEN 
                      ra = pbc(rvec,cell)
                      r  = SQRT(DOT_PRODUCT(ra,ra))
                      M(idim1,idim2) = M(idim1,idim2) + sfact / r
                      IF (need_f) THEN
                         t2    = - 1.0_dp /(r*r) * sfact
                         drvec = ra / r * q1 * q2 
                         d_el(1:3,iparticle1) = d_el(1:3,iparticle1) + t2 * drvec
                         d_el(1:3,iparticle2) = d_el(1:3,iparticle2) - t2 * drvec
                      END IF
                   END IF
                   M(idim1,idim2) = M(idim1,idim2) - sfact * (g_ewald + fac3*(rc2+rc12))
                   M(idim2,idim1) = M(idim1,idim2)
                END DO
             END DO
          END DO
       END DO
       IF (need_f) THEN
          IF (debug_this_module) CALL debug_ewald_ddapc_pot(qs_env, qmmm_periodic, multipole_section,&
               particle_set, M, radii, charges, d_el, error=error)
          CALL cp_decpl_ddapc_forces(qs_env, M, charges, dq, d_el, particle_set, error)
          DEALLOCATE(d_el, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)          
       END IF
    END IF
  END SUBROUTINE ewald_ddapc_pot

!!****f* cp_ddapc/debug_ewald_ddapc_pot [1.0] *
!!
!!   NAME
!!     debug_ewald_ddapc_pot
!!
!!   FUNCTION
!!     Routine to Debug the ewald contributions to the coupling/decoupling technique
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
!!     09.2005 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE debug_ewald_ddapc_pot(qs_env, qmmm_periodic, multipole_section, particle_set,&
            M, radii, charges, d_el, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: multipole_section    
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    LOGICAL, INTENT(IN)                      :: qmmm_periodic
    REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: M
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: radii
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: charges
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: d_el
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                            :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_ddapc_pot', &
         routineP = moduleN//':'//routineN      
    LOGICAL :: failure    
    INTEGER :: ip, k, j, stat
    REAL(KIND=dp) :: Dx, std_value, numer_energy(2), numer_forces, my_error
    REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: Mloc

    failure = .FALSE.
    IF (.NOT.failure) THEN
       Dx = 0.001_dp
       ALLOCATE(Mloc(SIZE(charges),SIZE(charges)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
       Atom: DO ip = 1, SIZE(particle_set)
          Coord: DO k = 1, 3
             numer_energy = 0.0_dp
             std_value = particle_set(ip)%r(k)             
             DO j = 1, 2 
                particle_set(ip)%r(k) = std_value - (-1.0_dp)**j * Dx
                Mloc = 0.0_dp
                CALL ewald_ddapc_pot(qs_env, qmmm_periodic, multipole_section, particle_set, Mloc, radii, &
                     .FALSE., charges=charges, error=error)
                numer_energy(j) = 0.5_dp * DOT_PRODUCT(charges,MATMUL(Mloc,charges))
             END DO
             particle_set(ip)%r(k) = std_value
             numer_forces = (numer_energy(1) - numer_energy(2) ) / (2.0_dp*Dx)
             my_error = 0.0_dp
             IF (ABS(numer_forces) >= 1.0E-7_dp) & 
                  my_error = ABS((numer_forces - d_el(k,ip))/numer_forces*100.0_dp)
             WRITE(*,'(A,I5,A,I5,3F15.9)')"ATOM NUMBER ::",ip," COORD:",k,numer_forces,d_el(k,ip),my_error
             CPPostcondition(my_error<=5.0_dp,cp_failure_level,routineP,error,failure)
          END DO Coord
       END DO Atom
       DEALLOCATE(Mloc, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)         
    END IF

  END SUBROUTINE debug_ewald_ddapc_pot

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
    IMPLICIT NONE
    ! Arguments
    TYPE(coeff_type)                         :: v_hartree_gspace, rho_tot_gspace 
    TYPE(section_vals_type), POINTER         :: multipole_section
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: M, AmI
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: radii, charges
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'modify_hartree_pot', &
      routineP = moduleN//':'//routineN          
    LOGICAL :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: uv, cv
    TYPE(pw_type), POINTER                   :: g_hartree
    REAL(kind=dp) :: gcut, Vol, fac, fac2, g2, gcut2, gfunc, w
    INTEGER       :: ig, idim, iparticle, igauss, stat
    REAL(KINd=dp) :: rc, rc2, sfac, gvec(3), rvec(3), arg
    COMPLEX(KIND=dp) :: phase, g_corr

    failure=.FALSE.
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
                phase  = CMPLX(COS(arg),-SIN(arg))
                gfunc  = EXP(-g2*rc2/4.0_dp)
                g_corr = g_corr + gfunc * cv(idim) * phase 
             END DO
          END DO
          g_corr = g_corr * w
          g_hartree%cc(ig) = g_hartree%cc(ig) + sfac*g_corr/Vol
       END DO
    END IF

  END SUBROUTINE modify_hartree_pot

!!****f* cp_ddapc/cp_decpl_ddapc_forces [1.0] *
!!
!!   NAME
!!     cp_decpl_ddapc_forces
!!
!!   FUNCTION
!!     Evaluation of the pulay forces due to the fitted charge density
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
  SUBROUTINE cp_decpl_ddapc_forces(qs_env, M, charges, dq, d_el, particle_set, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(particle_type), DIMENSION(:), &
         POINTER                             :: particle_set
    REAL(KIND=dp), POINTER, DIMENSION(:,:)   :: M, d_el
    REAL(KIND=dp), POINTER, DIMENSION(:)     :: charges    
    REAL(KIND=dp), POINTER, DIMENSION(:,:,:) :: dq
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_decpl_ddapc_forces', &
      routineP = moduleN//':'//routineN          
    LOGICAL :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)   :: uv
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: chf
    INTEGER :: iatom, natom, i, stat, ikind, k, j
    TYPE(qs_force_type), DIMENSION(:), &
         POINTER                                     :: force
    TYPE(atomic_kind_type), DIMENSION(:), &
         POINTER                                     :: atomic_kind_set
    INTEGER, ALLOCATABLE, DIMENSION(:)               :: atom_of_kind, kind_of
    
    failure = .FALSE.
    IF (.NOT.failure) THEN
       natom = SIZE(particle_set)
       CALL get_qs_env(qs_env=qs_env,&
                       atomic_kind_set=atomic_kind_set,&
                       force=force)
       ALLOCATE (atom_of_kind(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (kind_of(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE (chf(3,natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                                atom_of_kind=atom_of_kind,&
                                kind_of=kind_of)

       ALLOCATE(uv(SIZE(M,1)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       uv   = MATMUL(M,charges)
       DO k = 1, natom
          DO j = 1, 3
             chf(j,k) =  DOT_PRODUCT(uv,dq(:,k,j))
          END DO
       END DO
       DO iatom=1,natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          force(ikind)%ch_pulay(1:3,i) = force(ikind)%ch_pulay(1:3,i) + chf(1:3,iatom) + d_el(1:3,iatom) 
       END DO
       DEALLOCATE(atom_of_kind, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(kind_of, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(chf, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(uv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE cp_decpl_ddapc_forces

!!****f* cp_ddapc/reset_ch_pulay [1.0] *
!!
!!   NAME
!!     reset_ch_pulay
!!
!!   FUNCTION
!!     Evaluation of the pulay forces due to the fitted charge density
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
  SUBROUTINE reset_ch_pulay(qs_env, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'reset_ch_pulay', &
      routineP = moduleN//':'//routineN          
    LOGICAL :: failure
    INTEGER :: ind
    TYPE(qs_force_type), DIMENSION(:), &
         POINTER                                     :: force
    failure = .FALSE.
    IF (.NOT.failure) THEN
       CALL get_qs_env(qs_env=qs_env,&
                       force=force)
       DO ind = 1, SIZE(force)
          force(ind)%ch_pulay = 0.0_dp
       END DO
    END IF

  END SUBROUTINE reset_ch_pulay

END MODULE cp_ddapc
