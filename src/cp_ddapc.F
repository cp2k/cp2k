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
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALSPACE,&
                                             REALDATA3D,&
                                             RECIPROCALSPACE,&
                                             pw_inject_NoPBC,&
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
  RECURSIVE SUBROUTINE get_ddapc(qs_env,calc_force,fitcharge_section,qout1,qout2,error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: calc_force
    TYPE(section_vals_type), POINTER         :: fitcharge_section
    REAL(KIND=dp), DIMENSION(:), POINTER, OPTIONAL :: qout1, qout2
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
       ch_dens = REAL(qs_env%dft_control%charge,KIND=dp)/Vol
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
          CPPostcondition(SIZE(qout1)==SIZE(qv),cp_failure_level,routineP,error,failure)
          qout1 = qv
       END IF
       IF (PRESENT(qout2)) THEN
          CPPostcondition(SIZE(qout2)==SIZE(qtot),cp_failure_level,routineP,error,failure)
          qout2 = qtot
       END IF
       IF (cp_should_output(logger%iter_info,print_gen,error=error).AND.&
            logger%para_env%mepos==logger%para_env%source) THEN
          WRITE(iw,'(/T3,A)')" Charges on QM atoms used to couple/decouple the periodic images:"
          DO i = 1, SIZE(particle_set)
             WRITE(iw,'(T6,A,I5,F12.6)')"Atom: ",i,qtot(i)
          END DO
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
  SUBROUTINE cp_decpl_ddapc(qs_env,v_hartree_gspace,error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(coeff_type)                         :: v_hartree_gspace    
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'cp_decpl_ddapc', &
      routineP = moduleN//':'//routineN      
    TYPE(section_vals_type), POINTER         :: fitcharge_section,&
                                                force_env_section,&
                                                poisson_section,&
                                                qmmm_periodic_section
    TYPE(particle_type), DIMENSION(:), &
         POINTER                                 :: particle_set
    LOGICAL :: failure, qmmm_periodic
    INTEGER :: stat, handle, my_val
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges

    failure = .FALSE.
    IF (.NOT.failure) THEN
       NULLIFY(fitcharge_section, poisson_section, force_env_section, particle_set,&
               qmmm_periodic_section)
       CALL get_qs_env(qs_env=qs_env,&
                       input=force_env_section,&
                       particle_set=particle_set,&
                       error=error)
       poisson_section => section_vals_get_subs_vals(force_env_section,"DFT%POISSON_DFT", error=error)
       qmmm_periodic_section   => section_vals_get_subs_vals(force_env_section,"QMMM%PERIODIC", error=error)
       CALL section_vals_get(qmmm_periodic_section,explicit=qmmm_periodic, error=error)
       IF (.NOT.qmmm_periodic) THEN
          CALL section_vals_val_get(poisson_section,"POISSON_SOLVER",i_val=my_val,error=error)
          SELECT CASE (my_val)
          CASE (use_multipole)       
             fitcharge_section => section_vals_get_subs_vals(force_env_section,"PROPERTIES%FIT_CHARGE", error=error)
          CASE DEFAULT
             RETURN
          END SELECT
       ELSE
          fitcharge_section => section_vals_get_subs_vals(qmmm_periodic_section,"FIT_CHARGE", error=error)
       END IF
       ! First we evaluate the charges at the corresponding SCF STEP
       CALL get_ddapc(qs_env,.FALSE.,fitcharge_section,qout2=charges,error=error)
    END IF
  END SUBROUTINE cp_decpl_ddapc

END MODULE cp_ddapc
