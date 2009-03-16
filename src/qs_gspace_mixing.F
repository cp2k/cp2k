!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE qs_gspace_mixing

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_trace,&
                                             cp_fm_upper_to_full
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: copy_fm_to_sm,&
                                             copy_sm_to_fm,&
                                             cp_sm_fm_multiply
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all,&
                                             invert_matrix
  USE message_passing,                 ONLY: mp_sum
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_scale,&
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
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_methods,                  ONLY: cp_sm_mix
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE sparse_matrix_types,             ONLY: cp_sm_scale_and_add,&
                                             cp_sm_sm_trace,&
                                             get_matrix_info,&
                                             real_matrix_p_type,&
                                             scale_matrix,&
                                             set_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gspace_mixing'

  PUBLIC :: broyden_mixing, gmix_potential_only, gspace_density_mixing, pulay_mixing, self_consistency_check

CONTAINS

  SUBROUTINE self_consistency_check(qs_env,scf_env,p_out,delta,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_out
    REAL(KIND=dp), INTENT(INOUT)             :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'self_consistency_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_q, p_in

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(p_out),cp_failure_level,routineP,error,failure)
    NULLIFY(rho, matrix_q, p_in, para_env)

    CALL get_qs_env(qs_env=qs_env,  &
                    para_env=para_env,&
                    rho=rho,&
                    error=error)

    p_in => rho%rho_ao
    matrix_q => scf_env%p_delta
    nspin = SIZE(p_in)

    ! Compute the difference (p_out - p_in)and check convergence
    delta = 0.0_dp
    DO ispin = 1,nspin
       CALL set_matrix(matrix_q(ispin)%matrix, 0.0_dp)
 !      CALL cp_sm_mix(m1=p_in(ispin)%matrix,m2=p_out(ispin)%matrix,&
       CALL cp_sm_mix(m1=p_out(ispin)%matrix,m2=p_in(ispin)%matrix,&
            p_mix=1.0_dp,delta=tmp, para_env=para_env,&
            m3=matrix_q(ispin)%matrix,error=error)
      delta = MAX(tmp,delta)
    END DO

   CALL timestop(handle)

  END SUBROUTINE self_consistency_check

! *****************************************************************************
!> \brief 
!> \param qs_env the qs_environment where to perform the scf procedure
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note At the moment the mixing of density matrix in g-space does not work, 
!>       The g-space mixing is then performed only on the grid and affects the potential only
!>       This is done in gmix_potential_only
!> \par History
!>      02.2009 
!> \author MI
! *****************************************************************************
      
  SUBROUTINE  gspace_density_mixing(qs_env,scf_env,p_out,delta,eps_scf,diis,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: p_out
    REAL(KIND=dp), INTENT(INOUT)             :: delta
    REAL(KIND=dp), INTENT(IN)                :: eps_scf
    LOGICAL, INTENT(in)                      :: diis
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gspace_density_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ig, ispin, nao, ng, &
                                                nspin
    LOGICAL                                  :: failure, gapw
    REAL(KIND=dp)                            :: a_mix, b_mix, t1, t2, tmp, &
                                                trace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: g2
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_p_type)                          :: kerker_pot_rspace, rhog
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_q, matrix_s, p_in

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(p_out),cp_failure_level,routineP,error,failure)
    NULLIFY(auxbas_pw_pool, matrix_s, rho, matrix_q, p_in, para_env, pw_env)
    NULLIFY(g2, pw_grid)

    CALL get_qs_env(qs_env=qs_env,  &
                    dft_control=dft_control,&
                    matrix_s=matrix_s,&
                    pw_env=pw_env,&
                    rho=rho,& 
                    para_env=para_env,&
                    error=error)
    p_in => rho%rho_ao
    matrix_q => scf_env%p_delta
    nspin = SIZE(p_in)
    CALL get_matrix_info(matrix= rho%rho_ao(1)%matrix, nrow=nao)

    ! Compute the difference (p_out - p_in)and check convergence
    delta = 0.0_dp
    DO ispin = 1,nspin
       CALL set_matrix(matrix_q(ispin)%matrix, 0.0_dp)
 !      CALL cp_sm_mix(m1=p_in(ispin)%matrix,m2=p_out(ispin)%matrix,&
       CALL cp_sm_mix(m1=p_out(ispin)%matrix,m2=p_in(ispin)%matrix,&
            p_mix=1.0_dp,delta=tmp, para_env=para_env,&
            m3=matrix_q(ispin)%matrix,error=error)
      delta = MAX(tmp,delta)
    END DO
    ! For the moment skip the density matrix mixing in g-space
    ! In its place, Matthias method of mixing only the potential is restored
      scf_env%skip_mixing = .TRUE.
    IF(delta<eps_scf .OR. diis .OR. scf_env%skip_mixing) THEN
      scf_env%skip_mixing = .FALSE.
      CALL timestop(handle)
      RETURN
    END IF

    a_mix = scf_env%p_mix_alpha
    b_mix = scf_env%p_mix_beta

    !Collocate the density difference on the grid in realspace
    gapw = dft_control%qs_control%gapw
    CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                           rhog%pw,&
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE,error=error)
    CALL pw_pool_create_pw(auxbas_pw_pool,&
                           kerker_pot_rspace%pw,&
                           use_data=REALDATA3D, &
                           in_space=REALSPACE,error=error)
    DO ispin=1,nspin
       !collocate the real space density on the real space grid 
       CALL calculate_rho_elec(matrix_p=matrix_q(ispin)%matrix,&
            rho=rho%rho_r(ispin),&
            rho_gspace=rho%rho_g(ispin),&
            total_rho=rho%tot_rho_r(ispin),&
            qs_env=qs_env,soft_valid=gapw,error=error)

       CALL pw_zero(rhog%pw)
       ! point pw
       pw_grid => auxbas_pw_pool % pw_grid
       ng = SIZE ( pw_grid % gsq )
       g2 => pw_grid %gsq
       CALL pw_transfer ( rho%rho_g(ispin)%pw, rhog%pw )
!       rhog%pw =>  rho%rho_g(ispin)%pw
       DO ig=1,ng
!          rhog%pw%cc(ig) = rhog%pw%cc(ig) *  (1.0_dp-a_mix)*g2(ig)/(g2(ig)+b_mix*b_mix)
          rhog%pw%cc(ig) = rhog%pw%cc(ig) *  g2(ig)/(g2(ig)+b_mix*b_mix)
!dbg
 !         rhog%pw%cc(ig) = rhog%pw%cc(ig) *  (1.0_dp-a_mix)
!          rhog%pw%cc(ig) = rhog%pw%cc(ig)
!dbg
       ENDDO  
       CALL pw_transfer ( rhog%pw, kerker_pot_rspace%pw)
       
       CALL pw_scale ( kerker_pot_rspace%pw, pw_grid%dvol)

       CALL set_matrix(matrix_q(ispin)%matrix, 0.0_dp)
       CALL integrate_v_rspace(v_rspace=kerker_pot_rspace, &
            h=matrix_q(ispin), qs_env=qs_env,&
            calculate_forces =.FALSE., gapw=gapw, error=error)
       
       ! Impose Tr[S^-1 Q S^-1 S] = 0, i.e. the correction to P_out
       ! should not change the total number of electrons 
    IF(.FALSE.) THEN
       CALL copy_sm_to_fm(matrix_q(ispin)%matrix,scf_env%scf_work2,error=error)
       CALL cp_fm_upper_to_full(scf_env%scf_work2,scf_env%scf_work1(ispin)%matrix,error=error) 
       CALL cp_fm_trace(scf_env%s_minus_one,scf_env%scf_work2,trace,error=error)
       trace = -trace / REAL(nao)
       CALL cp_sm_scale_and_add(matrix_a=matrix_q(ispin)%matrix,alpha=1.0_dp,&
                                matrix_b=matrix_s(1)%matrix,beta=trace,error=error)
    END IF

!       CALL cp_sm_sm_trace(p_out(ispin)%matrix,matrix_s(1)%matrix,trace, para_env, error)
       CALL cp_sm_sm_trace(p_in(ispin)%matrix,matrix_s(1)%matrix,trace, para_env, error)
       t1 = trace
 
       CALL cp_sm_fm_multiply(matrix_q(ispin)%matrix,scf_env%s_minus_one,&
            scf_env%scf_work2, nao,error=error)
       CALL cp_fm_gemm("N","N",nao,nao,nao,1.0_dp,scf_env%s_minus_one,&
            scf_env%scf_work2,0.0_dp,scf_env%scf_work1(ispin)%matrix,error=error )

!       a_mix = 1.0_dp-scf_env%p_mix_alpha
!       CALL copy_fm_to_sm(scf_env%scf_work1(ispin)%matrix,&
!            p_out(ispin)%matrix,a_mix,1.0_dp)

       a_mix = scf_env%p_mix_alpha
       CALL copy_fm_to_sm(scf_env%scf_work1(ispin)%matrix,&
            p_out(ispin)%matrix,a_mix,0.0_dp)
       CALL cp_sm_scale_and_add(p_out(ispin)%matrix,1.0_dp,p_in(ispin)%matrix,1.0_dp,error)
       CALL cp_sm_sm_trace(p_out(ispin)%matrix,matrix_s(1)%matrix,trace, para_env, error)
       t2 =t1/ trace
       CALL scale_matrix(p_out(ispin)%matrix,t2)
       
    END DO  !ispin
    CALL pw_pool_give_back_pw(auxbas_pw_pool,rhog%pw,error=error)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,kerker_pot_rspace%pw,error=error)

    CALL timestop(handle)

  END SUBROUTINE gspace_density_mixing

! *****************************************************************************
!> \brief G-space mixing performed via the Kerker damping on the density on the grid
!>        thus affecting only the caluclation of the potential, but not the denisity matrix 
!> \par History
!>      02.2009 
!> \author MI
! *****************************************************************************

  SUBROUTINE gmix_potential_only (scf_env, rho, error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gmix_potential_only', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), DIMENSION(:), POINTER       :: cc_new
    INTEGER                                  :: ig, ispin, ng, nspin
    LOGICAL                                  :: failure
    REAL(dp)                                 :: f_mix
    REAL(dp), DIMENSION(:), POINTER          :: g2

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(scf_env%rhoin_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_env%kerker_factor),cp_failure_level,routineP,error,failure)
    NULLIFY(cc_new, g2)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(rho%rho_g(1)%pw%pw_grid%gsq)
    g2 => rho%rho_g(1)%pw%pw_grid%gsq

    DO ispin = 1,nspin
      cc_new =>  rho%rho_g(ispin)%pw%cc
      DO ig = 1,ng
         f_mix = scf_env%p_mix_alpha*scf_env%kerker_factor(ig)
         cc_new(ig) = (1.0_dp-f_mix)*scf_env%rhoin_buffer(1,ispin)%cc(ig)+f_mix*cc_new(ig) 
         scf_env%rhoin_buffer(1,ispin)%cc(ig) = cc_new(ig)
      END DO
    END DO

  END SUBROUTINE  gmix_potential_only

! *****************************************************************************
!> \brief Pulay Mixing using as metrics for the residual the multiplication by the Kerer damping factor
!>        The mixing is applied directly on the density expansion in reciprocal space, therefore it affects the potentials
!>        on the grid but not the density matrix
!> \par History
!>      03.2009 
!> \author MI
! *****************************************************************************

  SUBROUTINE pulay_mixing(scf_env, rho, para_env, error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pulay_mixing', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), ALLOCATABLE, DIMENSION(:)   :: cc_mix
    INTEGER                                  :: handle, ib, ibb, ig, ispin, &
                                                istat, jb, nb, nb1, nbuffer, &
                                                ng, ng_tot, ngpts, nspin
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha_kerker, alpha_pulay, &
                                                beta, f_mix, norm, vol
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: ev
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a, b
    REAL(dp), DIMENSION(:), POINTER          :: g2

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(scf_env%res_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_env%rhoin_buffer),cp_failure_level,routineP,error,failure)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(scf_env%res_buffer(1,1)%cc)
    ng_tot = SIZE(rho%rho_g(1)%pw%pw_grid%gsq)
    ngpts = rho%rho_g(1)%pw%pw_grid%ngpts
    g2 => rho%rho_g(1)%pw%pw_grid%gsq
    vol = rho%rho_g(1)%pw%pw_grid%vol

    alpha_kerker=scf_env%p_mix_alpha
    beta=scf_env%p_mix_beta
    alpha_pulay = scf_env%pulay_alpha
    nbuffer = scf_env%npulay

    ib = MODULO(scf_env%mixing_ncall,nbuffer) + 1
    scf_env%mixing_ncall = scf_env%mixing_ncall + 1
    nb = MIN(scf_env%mixing_ncall,nbuffer)
    ibb = MODULO(scf_env%mixing_ncall,nbuffer) + 1

    nb1 = nb + 1
    ALLOCATE(a(nb1,nb1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    a=0.0_dp
    ALLOCATE(b(nb1,nb1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    b=0.0_dp
    ALLOCATE(ev(nb1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ev=0.0_dp
    ALLOCATE(cc_mix(ng), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
  
    DO ispin = 1,nspin
      scf_env%res_buffer(ib,ispin)%cc(:) = CMPLX(0._dp,0._dp,KIND=dp)
      norm = 0.0_dp
      DO ig = 1,ng
          f_mix = scf_env%kerker_factor(ig)
!           f_mix = (g2(ig)+beta*beta)/(g2(ig)+beta*beta*0.01_dp*0.01_dp)
          scf_env%res_buffer(ib,ispin)%cc(ig) = f_mix*(rho%rho_g(ispin)%pw%cc(ig)-scf_env%rhoin_buffer(ib,ispin)%cc(ig))
      END DO

      DO jb = 1,nb
         scf_env%pulay_matrix(jb,ib) = 0.0_dp
         IF(rho%rho_g(1)%pw%pw_grid%have_g0) THEN
           scf_env%pulay_matrix(jb,ib) = scf_env%pulay_matrix(jb,ib) +  (&
                             REAL(scf_env%res_buffer(jb,ispin)%cc(1),dp)*REAL(scf_env%res_buffer(ib,ispin)%cc(1),dp) +&
                             AIMAG(scf_env%res_buffer(jb,ispin)%cc(1))*AIMAG(scf_env%res_buffer(ib,ispin)%cc(1)))
         ELSE
           f_mix = 1.0_dp*scf_env%kerker_factor(1)
           scf_env%pulay_matrix(jb,ib) = scf_env%pulay_matrix(jb,ib) + f_mix* (&
                             REAL(scf_env%res_buffer(jb,ispin)%cc(1),dp)*REAL(scf_env%res_buffer(ib,ispin)%cc(1),dp) +&
                             AIMAG(scf_env%res_buffer(jb,ispin)%cc(1))*AIMAG(scf_env%res_buffer(ib,ispin)%cc(1)))
         END IF
         DO ig = 2,ng
           f_mix = 1.0_dp*scf_env%kerker_factor(ig)
           scf_env%pulay_matrix(jb,ib) = scf_env%pulay_matrix(jb,ib) + f_mix* (&
                             REAL(scf_env%res_buffer(jb,ispin)%cc(ig),dp)*REAL(scf_env%res_buffer(ib,ispin)%cc(ig),dp) +&
                             AIMAG(scf_env%res_buffer(jb,ispin)%cc(ig))*AIMAG(scf_env%res_buffer(ib,ispin)%cc(ig)))
         END DO
         CALL mp_sum(scf_env%pulay_matrix(jb,ib),para_env%group)
         scf_env%pulay_matrix(jb,ib)=scf_env%pulay_matrix(jb,ib)* vol *vol*100.0_dp
         scf_env%pulay_matrix(ib,jb) = scf_env%pulay_matrix(jb,ib)
      END DO

      IF( nb==1) THEN
        DO ig = 1,ng
           f_mix = alpha_kerker*scf_env%kerker_factor(ig)
           cc_mix(ig) = rho%rho_g(ispin)%pw%cc(ig) - scf_env%rhoin_buffer(ib,ispin)%cc(ig)
           rho%rho_g(ispin)%pw%cc(ig) = f_mix*cc_mix(ig)+scf_env%rhoin_buffer(ib,ispin)%cc(ig)
           scf_env%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) 
        END DO
      ELSE

        b(1:nb,1:nb) = scf_env%pulay_matrix(1:nb, 1:nb)
        b(nb1,1:nb)  = 1.0_dp
        b(1:nb,nb1)  = 1.0_dp
        b(nb1,nb1)   = 0.0_dp

        CALL diamat_all(b(1:nb1,1:nb1),ev(1:nb1),error=error)

        a(1:nb1,1:nb1) = b(1:nb1,1:nb1)
        DO jb = 1,nb1
          IF(ABS(ev(jb)) < 1.E-5_dp ) THEN
            a(1:nb1,jb) = 0.0_dp
          ELSE
            a(1:nb1,jb) = a(1:nb1,jb) / ev(jb)
          END IF
        END DO
        ev(1:nb) = MATMUL(a(1:nb,1:nb1),b(nb1,1:nb1))

        cc_mix = CMPLX(0._dp,0._dp,KIND=dp) 
        DO jb = 1,nb
          DO ig = 1,ng
          !   cc_mix(ig) = cc_mix(ig) +  ev(jb) * scf_env%rhoin_buffer(jb,ispin)%cc(ig)
           cc_mix(ig) = cc_mix(ig) +  ev(jb) * (scf_env%rhoin_buffer(jb,ispin)%cc(ig) +scf_env%res_buffer(jb,ispin)%cc(ig))
          END DO
        END DO
        scf_env%rhoin_buffer(ibb,ispin)%cc = CMPLX(0._dp,0._dp,KIND=dp)
        IF(alpha_pulay > 0.0_dp) THEN
          DO ig = 1,ng
            f_mix = alpha_pulay*scf_env%kerker_factor(ig)
            rho%rho_g(ispin)%pw%cc(ig) =  f_mix*rho%rho_g(ispin)%pw%cc(ig) + (1.0_dp-f_mix) *cc_mix(ig)
            scf_env%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) 
          END DO
        ELSE
          DO ig = 1,ng
            rho%rho_g(ispin)%pw%cc(ig) =  cc_mix(ig)
            scf_env%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) 
          END DO
        END IF

      END IF  ! nb==1
    END DO ! ispin

    DEALLOCATE(a, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(b, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ev, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(cc_mix, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE pulay_mixing

! *****************************************************************************
!> \brief Broyden Mixing using as metrics for the residual the multiplication by the Kerer damping factor
!>        The mixing is applied directly on the density expansion in reciprocal space, therefore it affects the potentials
!>        on the grid but not the density matrix
!> \par History
!>      03.2009 
!> \author MI
! *****************************************************************************

  SUBROUTINE broyden_mixing(scf_env, rho, para_env, error)

    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'broyden_mixing', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp)                              :: cc_mix
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:)   :: res_rho
    INTEGER                                  :: handle, ib, ibb, ig, ispin, &
                                                istat, jb, kb, nb, nbuffer, &
                                                ng, ng_tot, ngpts, ngpts_tot, &
                                                nspin
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha, beta, delta_norm, &
                                                f_mix, inv_err, res_norm, &
                                                vol, w0
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: c, g
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a, b
    REAL(dp), DIMENSION(:), POINTER          :: g2

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(scf_env%res_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_env%rhoin_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_env%drho_buffer),cp_failure_level,routineP,error,failure)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(scf_env%res_buffer(1,1)%cc)
    ng_tot = SIZE(rho%rho_g(1)%pw%pw_grid%gsq)
    ngpts_tot = rho%rho_g(1)%pw%pw_grid%ngpts
    ngpts = ng
    CALL mp_sum(ngpts,para_env%group)
    g2 => rho%rho_g(1)%pw%pw_grid%gsq
    vol = rho%rho_g(1)%pw%pw_grid%vol

    alpha=scf_env%p_mix_alpha
    beta=scf_env%p_mix_beta
    w0 = scf_env%broy_w0
    nbuffer = scf_env%nbroy

    ib = MODULO(scf_env%mixing_ncall,nbuffer) + 1
    scf_env%mixing_ncall = scf_env%mixing_ncall + 1
    nb = MIN(scf_env%mixing_ncall,nbuffer)
    ibb = MODULO(scf_env%mixing_ncall,nbuffer) + 1

    ALLOCATE(a(ib-1,ib-1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    a=0.0_dp
    ALLOCATE(b(ib-1,ib-1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    b=0.0_dp
    ALLOCATE(c(ib-1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    c=0.0_dp
    ALLOCATE(g(ib-1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    g=0.0_dp
    ALLOCATE(res_rho(ng), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    DO ispin = 1,nspin
      DO ig = 1,ng
         res_rho(ig) = rho%rho_g(ispin)%pw%cc(ig) - scf_env%rhoin_buffer(ib,ispin)%cc(ig)
      END DO

      IF(ib==1) THEN
        DO ig = 1,ng
          scf_env%broy_last_res(ispin)%cc(ig) = res_rho(ig)
          f_mix = alpha*scf_env%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = scf_env%rhoin_buffer(ib,ispin)%cc(ig) + f_mix*res_rho(ig) 
          scf_env%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO
      ELSE

        delta_norm = 0.0_dp
        res_norm = 0.0_dp
        DO ig = 1,ng
          scf_env%res_buffer(ib-1,ispin)%cc(ig) = res_rho(ig) - scf_env%broy_last_res(ispin)%cc(ig)
          scf_env%broy_last_res(ispin)%cc(ig) = res_rho(ig)
          res_norm = res_norm + &
                     REAL(res_rho(ig),dp)*REAL(res_rho(ig),dp) + AIMAG(res_rho(ig))*AIMAG(res_rho(ig)) 
          delta_norm = delta_norm + &
                       REAL(scf_env%res_buffer(ib-1,ispin)%cc(ig),dp)*REAL(scf_env%res_buffer(ib-1,ispin)%cc(ig),dp) +&
                       AIMAG(scf_env%res_buffer(ib-1,ispin)%cc(ig))*AIMAG(scf_env%res_buffer(ib-1,ispin)%cc(ig))
          scf_env%drho_buffer(ib-1,ispin)%cc(ig) = scf_env%rhoin_buffer(ib,ispin)%cc(ig) - &
                                                   scf_env%rhoin_buffer(ib-1,ispin)%cc(ig)
        END DO 
        CALL mp_sum(delta_norm, para_env%group)
        delta_norm = SQRT(delta_norm)
        CALL mp_sum(res_norm,para_env%group)
        res_norm = SQRT(res_norm)

        scf_env%res_buffer(ib-1,ispin)%cc(:) = scf_env%res_buffer(ib-1,ispin)%cc(:)/delta_norm
        scf_env%drho_buffer(ib-1,ispin)%cc(:) = scf_env%drho_buffer(ib-1,ispin)%cc(:)/delta_norm

        a(:,:) = 0.0_dp
        DO ig = 1,ng
          f_mix = alpha*scf_env%kerker_factor(ig)
          scf_env%drho_buffer(ib-1,ispin)%cc(ig) = f_mix*scf_env%res_buffer(ib-1,ispin)%cc(ig)+&
                                                    scf_env%drho_buffer(ib-1,ispin)%cc(ig)
        END DO 
        DO jb = 1,ib-1
           DO kb = jb,ib-1
              DO ig = 1,ng
                  a(kb,jb) =  a(kb,jb) +&
                          REAL(scf_env%res_buffer(jb,ispin)%cc(ig),dp)*REAL(scf_env%res_buffer(kb,ispin)%cc(ig),dp) +&
                          AIMAG(scf_env%res_buffer(jb,ispin)%cc(ig))*AIMAG(scf_env%res_buffer(kb,ispin)%cc(ig))
              END DO
              a(jb,kb) = a(kb,jb)
           END DO
        END DO
        CALL mp_sum(a,para_env%group)

        C = 0.0_dp
        DO jb = 1,ib-1
           a(jb,jb) = w0 + a(jb,jb)
           DO ig = 1,ng
             c(jb) = c(jb) + REAL(scf_env%res_buffer(jb,ispin)%cc(ig),dp)*REAL(res_rho(ig),dp) + &
                     AIMAG(scf_env%res_buffer(jb,ispin)%cc(ig))*AIMAG(res_rho(ig))
           END DO
        END DO
        CALL mp_sum(c,para_env%group)
        CALL invert_matrix(a,b,inv_err,error=error)

        CALL dgemv('T',IB-1,IB-1,1.0_dp,B,IB-1,C,1,0.0_dp,G,1)

        DO ig = 1,ng
          cc_mix = CMPLX(0.0_dp,0.0_dp,kind=dp)
          DO jb =  1,ib-1
            cc_mix = cc_mix - G(jb)*scf_env%drho_buffer(jb,ispin)%cc(ig)
          END DO
          f_mix = alpha*scf_env%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = scf_env%rhoin_buffer(ib,ispin)%cc(ig) + f_mix * res_rho(ig) + cc_mix 
          scf_env%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO

      END IF

    END DO

    DEALLOCATE(a,b,c,g,res_rho, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

 END SUBROUTINE broyden_mixing

END MODULE qs_gspace_mixing
