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
  USE pw_methods,                      ONLY: pw_integrate_function,&
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
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_methods,                  ONLY: cp_sm_mix
  USE qs_scf_types,                    ONLY: broyden_mixing_new_nr,&
                                             broyden_mixing_nr,&
                                             cp_1d_z_p_type,&
                                             gspace_mixing_nr,&
                                             mixing_storage_type,&
                                             multisecant_mixing_nr,&
                                             pulay_mixing_nr,&
                                             qs_scf_env_type
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

  PUBLIC :: gspace_mixing, broyden_mixing, broyden_mixing_new, gmix_potential_only, &
            gspace_density_mixing, pulay_mixing, self_consistency_check

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

    a_mix = scf_env%mixing_store%alpha
    b_mix = scf_env%mixing_store%beta

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
!> \brief  Driver for the g-space mixing, calls the proper routine given the
!requested method
!> \par History
!>      05.2009 
!> \author MI
! *****************************************************************************
  SUBROUTINE gspace_mixing(scf_env, rho, para_env, error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gspace_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ig, ispin, ng, nspin
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha, f_mix

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(scf_env%mixing_store),cp_failure_level,routineP,error,failure)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(rho%rho_g(1)%pw%pw_grid%gsq)
    alpha = scf_env%mixing_store%alpha

    IF(scf_env%iter_count+1 <= scf_env%nskip_mixing) THEN
      ! skip mixing
      DO ispin = 1,nspin
        DO ig = 1,ng
          scf_env%mixing_store%rhoin_buffer(1,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO
      END DO
       scf_env%mixing_store%iter_method= "NoMix"
      RETURN
    END IF
    IF((scf_env%iter_count+1 - scf_env%nskip_mixing) <= scf_env%mixing_store%n_simple_mix) THEN
      ! Simple Kerker damping : rho(G)=rho_in(G)-alpha k(G)*(rho_out(G)-rho_in(G))
      DO ispin = 1,nspin
        DO ig = 1,ng
          f_mix = alpha*scf_env%mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = scf_env%mixing_store%rhoin_buffer(1,ispin)%cc(ig) + &
            f_mix*(rho%rho_g(ispin)%pw%cc(ig) - scf_env%mixing_store%rhoin_buffer(1,ispin)%cc(ig))
          scf_env%mixing_store%rhoin_buffer(1,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO
      END DO
       scf_env%mixing_store%iter_method= "Kerker"
      RETURN
    END IF

    CALL timeset(routineN,handle)

      IF(scf_env%mixing_method==gspace_mixing_nr) THEN
            CALL gmix_potential_only(scf_env%mixing_store, rho, error)
       scf_env%mixing_store%iter_method= "Kerker"
      ELSEIF(scf_env%mixing_method==pulay_mixing_nr) THEN
            CALL pulay_mixing(scf_env%mixing_store, rho, para_env, error)
       scf_env%mixing_store%iter_method= "Pulay"
      ELSEIF(scf_env%mixing_method==broyden_mixing_nr) THEN
            CALL broyden_mixing(scf_env%mixing_store, rho, para_env, error)
       scf_env%mixing_store%iter_method= "Broy"
      ELSEIF(scf_env%mixing_method==broyden_mixing_new_nr) THEN
            CALL broyden_mixing_new(scf_env%mixing_store, rho, para_env, error)
       scf_env%mixing_store%iter_method= "Broy"
      ELSEIF(scf_env%mixing_method==multisecant_mixing_nr) THEN
            CALL multisecant_mixing(scf_env%mixing_store, rho, para_env, error)
       scf_env%mixing_store%iter_method= "MSec"
      END IF

    DO ispin=1,nspin
      CALL pw_transfer(rho%rho_g(ispin)%pw,rho%rho_r(ispin)%pw)
      rho%tot_rho_r(ispin) = pw_integrate_function(rho%rho_r(ispin)%pw,isign=-1)
    END DO

    CALL timestop(handle)

  END SUBROUTINE gspace_mixing


! *****************************************************************************
!> \brief G-space mixing performed via the Kerker damping on the density on the grid
!>        thus affecting only the caluclation of the potential, but not the denisity matrix 
!> \par History
!>      02.2009 
!> \author MI
! *****************************************************************************

  SUBROUTINE gmix_potential_only (mixing_store, rho, error)

    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gmix_potential_only', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), DIMENSION(:), POINTER       :: cc_new
    INTEGER                                  :: handle, ig, ispin, ng, nspin
    LOGICAL                                  :: failure
    REAL(dp)                                 :: f_mix

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store%rhoin_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store%kerker_factor),cp_failure_level,routineP,error,failure)
    CPPrecondition(mixing_store%nbuffer==1,cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    NULLIFY(cc_new)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(rho%rho_g(1)%pw%pw_grid%gsq)

    DO ispin = 1,nspin
      cc_new =>  rho%rho_g(ispin)%pw%cc
      DO ig = 1,ng
         f_mix = mixing_store%alpha*mixing_store%kerker_factor(ig)
         cc_new(ig) = (1.0_dp-f_mix)*mixing_store%rhoin_buffer(1,ispin)%cc(ig)+f_mix*cc_new(ig) 
         mixing_store%rhoin_buffer(1,ispin)%cc(ig) = cc_new(ig)
      END DO
    END DO

    CALL timestop(handle)

  END SUBROUTINE  gmix_potential_only

! *****************************************************************************
!> \brief Pulay Mixing using as metrics for the residual the Kerer damping factor
!>        The mixing is applied directly on the density in reciprocal space,
!>        therefore it affects the potentials
!>        on the grid but not the density matrix
!> \par History
!>      03.2009 
!> \author MI
! *****************************************************************************

  SUBROUTINE pulay_mixing(mixing_store, rho, para_env, error)

    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pulay_mixing', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), ALLOCATABLE, DIMENSION(:)   :: cc_mix
    INTEGER                                  :: handle, i, ib, ibb, ig, &
                                                ispin, istat, jb, nb, nb1, &
                                                nbuffer, ng, nspin
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha_kerker, alpha_pulay, &
                                                beta, f_mix, inv_err, norm, &
                                                norm_c_inv, vol
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: alpha_c, ev
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a, b, c, c_inv

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store%res_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store%rhoin_buffer),cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(mixing_store%res_buffer(1,1)%cc)
    vol = rho%rho_g(1)%pw%pw_grid%vol

    alpha_kerker=mixing_store%alpha
    beta=mixing_store%beta
    alpha_pulay = mixing_store%pulay_alpha
    nbuffer = mixing_store%nbuffer

    ib = MODULO(mixing_store%ncall,nbuffer) + 1
    mixing_store%ncall = mixing_store%ncall + 1
    nb = MIN(mixing_store%ncall,nbuffer)
    ibb = MODULO(mixing_store%ncall,nbuffer) + 1

    nb1 = nb + 1
    ALLOCATE(a(nb1,nb1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    a=0.0_dp
    ALLOCATE(b(nb1,nb1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    b=0.0_dp
    ALLOCATE(c(nb,nb), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    c=0.0_dp
    ALLOCATE(c_inv(nb,nb), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    c_inv=0.0_dp
    ALLOCATE(alpha_c(nb), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    alpha_c=0.0_dp
    norm_c_inv = 0.0_dp
    ALLOCATE(ev(nb1), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ev=0.0_dp
    ALLOCATE(cc_mix(ng), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
  
    DO ispin = 1,nspin
      mixing_store%res_buffer(ib,ispin)%cc(:) = CMPLX(0._dp,0._dp,KIND=dp)
      norm = 0.0_dp
      DO ig = 1,ng
          f_mix = mixing_store%kerker_factor(ig)
          mixing_store%res_buffer(ib,ispin)%cc(ig) = f_mix*(rho%rho_g(ispin)%pw%cc(ig)-&
                              mixing_store%rhoin_buffer(ib,ispin)%cc(ig))
      END DO

      DO jb = 1,nb
         mixing_store%pulay_matrix(jb,ib) = 0.0_dp
         IF(rho%rho_g(1)%pw%pw_grid%have_g0) THEN
           mixing_store%pulay_matrix(jb,ib) = mixing_store%pulay_matrix(jb,ib) +  (&
                       REAL(mixing_store%res_buffer(jb,ispin)%cc(1),dp)*&
                       REAL(mixing_store%res_buffer(ib,ispin)%cc(1),dp) +&
                       AIMAG(mixing_store%res_buffer(jb,ispin)%cc(1))*&
                       AIMAG(mixing_store%res_buffer(ib,ispin)%cc(1)))
         ELSE
           f_mix = 1.0_dp*mixing_store%kerker_factor(1)
           mixing_store%pulay_matrix(jb,ib) = mixing_store%pulay_matrix(jb,ib) + f_mix* (&
                      REAL(mixing_store%res_buffer(jb,ispin)%cc(1),dp)*&
                      REAL(mixing_store%res_buffer(ib,ispin)%cc(1),dp) +&
                      AIMAG(mixing_store%res_buffer(jb,ispin)%cc(1))*&
                      AIMAG(mixing_store%res_buffer(ib,ispin)%cc(1)))
         END IF
         DO ig = 2,ng
           f_mix = 1.0_dp*mixing_store%kerker_factor(ig)
           mixing_store%pulay_matrix(jb,ib) = mixing_store%pulay_matrix(jb,ib) + &
                  f_mix*(REAL(mixing_store%res_buffer(jb,ispin)%cc(ig),dp)&
                         *REAL(mixing_store%res_buffer(ib,ispin)%cc(ig),dp) +&
                          AIMAG(mixing_store%res_buffer(jb,ispin)%cc(ig))*&
                          AIMAG(mixing_store%res_buffer(ib,ispin)%cc(ig)))
         END DO
         CALL mp_sum(mixing_store%pulay_matrix(jb,ib),para_env%group)
         mixing_store%pulay_matrix(jb,ib)=mixing_store%pulay_matrix(jb,ib)*vol*vol*100.0_dp
         mixing_store%pulay_matrix(ib,jb) = mixing_store%pulay_matrix(jb,ib)
      END DO

      IF( nb==1) THEN
        DO ig = 1,ng
           f_mix = alpha_kerker*mixing_store%kerker_factor(ig)
           cc_mix(ig) = rho%rho_g(ispin)%pw%cc(ig) -&
                         mixing_store%rhoin_buffer(ib,ispin)%cc(ig)
           rho%rho_g(ispin)%pw%cc(ig) = f_mix*cc_mix(ig)+ &
                                        mixing_store%rhoin_buffer(ib,ispin)%cc(ig)
           mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) 
        END DO
      ELSE

        b(1:nb,1:nb) = mixing_store%pulay_matrix(1:nb, 1:nb)
        c(1:nb,1:nb) = b(1:nb,1:nb)
        b(nb1,1:nb)  = 1.0_dp
        b(1:nb,nb1)  = 1.0_dp
        b(nb1,nb1)   = 0.0_dp

        CALL diamat_all(b(1:nb1,1:nb1),ev(1:nb1),error=error)
        CALL invert_matrix(c,c_inv,inv_err,error=error)
        alpha_c = 0.0_dp
        DO i = 1,nb
          DO jb=1,nb
            alpha_c(i) =  alpha_c(i) + c_inv(jb,i)
            norm_c_inv = norm_c_inv + c_inv(jb,i)
          END DO
        END DO
        alpha_c(1:nb) = alpha_c(1:nb)/norm_c_inv

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
           cc_mix(ig) = cc_mix(ig) +&
                        alpha_c(jb)*(mixing_store%rhoin_buffer(jb,ispin)%cc(ig) +&
                        mixing_store%res_buffer(jb,ispin)%cc(ig))
          END DO
        END DO
        mixing_store%rhoin_buffer(ibb,ispin)%cc = CMPLX(0._dp,0._dp,KIND=dp)
        IF(alpha_pulay > 0.0_dp) THEN
          DO ig = 1,ng
            f_mix = alpha_pulay*mixing_store%kerker_factor(ig)
            rho%rho_g(ispin)%pw%cc(ig) =  f_mix*rho%rho_g(ispin)%pw%cc(ig) +&
                      (1.0_dp-f_mix) *cc_mix(ig)
            mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) 
          END DO
        ELSE
          DO ig = 1,ng
            rho%rho_g(ispin)%pw%cc(ig) =  cc_mix(ig)
            mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) 
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
    DEALLOCATE(c,c_inv,alpha_c, STAT=istat)
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

  SUBROUTINE broyden_mixing(mixing_store, rho, para_env, error)

    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'broyden_mixing', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp)                              :: cc_mix
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:)   :: res_rho
    INTEGER                                  :: handle, ib, ibb, ig, ispin, &
                                                istat, jb, kb, nb, nbuffer, &
                                                ng, nspin
    LOGICAL                                  :: failure
    REAL(dp)                                 :: alpha, beta, delta_norm, &
                                                f_mix, inv_err, res_norm, w0
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: c, g
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a, b

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store%res_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store%rhoin_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store%drho_buffer),cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(mixing_store%res_buffer(1,1)%cc)

    alpha=mixing_store%alpha
    beta=mixing_store%beta
    w0 = mixing_store%broy_w0
    nbuffer = mixing_store%nbuffer

    ib = MODULO(mixing_store%ncall,nbuffer) + 1
    mixing_store%ncall = mixing_store%ncall + 1
 !   ib = MIN(mixing_store%ncall,nbuffer)
    nb = MIN(mixing_store%ncall,nbuffer)
    ibb = MODULO(mixing_store%ncall,nbuffer) + 1
  !  ibb = MIN(mixing_store%ncall+1,nbuffer)

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
         res_rho(ig) = rho%rho_g(ispin)%pw%cc(ig) - mixing_store%rhoin_buffer(ib,ispin)%cc(ig)
      END DO

      IF(ib==1) THEN
        DO ig = 1,ng
          mixing_store%last_res(ispin)%cc(ig) = res_rho(ig)
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(ib,ispin)%cc(ig) + f_mix*res_rho(ig) 
          mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO
      ELSE

        delta_norm = 0.0_dp
        res_norm = 0.0_dp
        DO ig = 1,ng
          mixing_store%res_buffer(ib-1,ispin)%cc(ig) = res_rho(ig) - mixing_store%last_res(ispin)%cc(ig)
          mixing_store%last_res(ispin)%cc(ig) = res_rho(ig)
          res_norm = res_norm + &
               REAL(res_rho(ig),dp)*REAL(res_rho(ig),dp) +&
               AIMAG(res_rho(ig))*AIMAG(res_rho(ig)) 
          delta_norm = delta_norm + &
               REAL(mixing_store%res_buffer(ib-1,ispin)%cc(ig),dp)*&
               REAL(mixing_store%res_buffer(ib-1,ispin)%cc(ig),dp) +&
               AIMAG(mixing_store%res_buffer(ib-1,ispin)%cc(ig))*&
               AIMAG(mixing_store%res_buffer(ib-1,ispin)%cc(ig))
          mixing_store%drho_buffer(ib-1,ispin)%cc(ig) = &
                     mixing_store%rhoin_buffer(ib,ispin)%cc(ig) - &
                     mixing_store%rhoin_buffer(ib-1,ispin)%cc(ig)
        END DO 
        CALL mp_sum(delta_norm, para_env%group)
        delta_norm = SQRT(delta_norm)
        CALL mp_sum(res_norm,para_env%group)
        res_norm = SQRT(res_norm)

        mixing_store%res_buffer(ib-1,ispin)%cc(:) = mixing_store%res_buffer(ib-1,ispin)%cc(:)/delta_norm
        mixing_store%drho_buffer(ib-1,ispin)%cc(:) = mixing_store%drho_buffer(ib-1,ispin)%cc(:)/delta_norm

        a(:,:) = 0.0_dp
        DO ig = 1,ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
          mixing_store%drho_buffer(ib-1,ispin)%cc(ig) = &
                 f_mix*mixing_store%res_buffer(ib-1,ispin)%cc(ig)+&
                 mixing_store%drho_buffer(ib-1,ispin)%cc(ig)
        END DO 
        DO jb = 1,ib-1
           DO kb = jb,ib-1
              DO ig = 1,ng
                  a(kb,jb) =  a(kb,jb) +&
                    REAL(mixing_store%res_buffer(jb,ispin)%cc(ig),dp)*&
                    REAL(mixing_store%res_buffer(kb,ispin)%cc(ig),dp) +&
                    AIMAG(mixing_store%res_buffer(jb,ispin)%cc(ig))*&
                    AIMAG(mixing_store%res_buffer(kb,ispin)%cc(ig))
              END DO
              a(jb,kb) = a(kb,jb)
           END DO
        END DO
        CALL mp_sum(a,para_env%group)

        C = 0.0_dp
        DO jb = 1,ib-1
           a(jb,jb) = w0 + a(jb,jb)
           DO ig = 1,ng
             c(jb) = c(jb) + REAL(mixing_store%res_buffer(jb,ispin)%cc(ig),dp)*REAL(res_rho(ig),dp) + &
                     AIMAG(mixing_store%res_buffer(jb,ispin)%cc(ig))*AIMAG(res_rho(ig))
           END DO
        END DO
        CALL mp_sum(c,para_env%group)
        CALL invert_matrix(a,b,inv_err,error=error)

        CALL dgemv('T',IB-1,IB-1,1.0_dp,B,IB-1,C,1,0.0_dp,G,1)

        DO ig = 1,ng
          cc_mix = CMPLX(0.0_dp,0.0_dp,kind=dp)
          DO jb =  1,ib-1
            cc_mix = cc_mix - G(jb)*mixing_store%drho_buffer(jb,ispin)%cc(ig)
          END DO
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(ib,ispin)%cc(ig) +&
                   f_mix * res_rho(ig) + cc_mix 
          mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO

!        IF(mixing_store%ncall>=nbuffer) THEN
!          DO jb = 1,nbuffer - 1
!             mixing_store%rhoin_buffer(jb,ispin)%cc(:) = mixing_store%rhoin_buffer(jb+1,ispin)%cc(:)
!             mixing_store%drho_buffer(jb,ispin)%cc(:) = mixing_store%drho_buffer(jb+1,ispin)%cc(:)
!             mixing_store%res_buffer(jb,ispin)%cc(:) = mixing_store%res_buffer(jb+1,ispin)%cc(:)
!          END DO 
!        END IF

      END IF

    END DO

    DEALLOCATE(a,b,c,g,res_rho, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

 END SUBROUTINE broyden_mixing

 SUBROUTINE broyden_mixing_new(mixing_store, rho, para_env,  error)

    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'broyden_mixing_new', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), ALLOCATABLE, DIMENSION(:)   :: delta_res_p, res_rho, &
                                                res_rho_p, tmp
    INTEGER                                  :: handle, ib, ibb, ig, info, &
                                                ispin, istat, jb, kb, kkb, &
                                                lwork, nb, nbuffer, ng, nspin
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: IWORK
    LOGICAL                                  :: failure, skip_bq
    REAL(dp) :: alpha, beta, delta, delta_p, delta_rhog, delta_rhog_p, f_mix, &
      imp, imp_j, inv_err, norm, norm_ig, rep, rep_j, vol, wc, wmax
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: aval, work_dgesdd
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a, au, av, b, bq, work
    REAL(dp), DIMENSION(:), POINTER          :: p_metric
    REAL(dp), DIMENSION(:, :), POINTER       :: fmat, gmat, smat, weight
    TYPE(cp_1d_z_p_type), DIMENSION(:), &
      POINTER                                :: tmp_z
    TYPE(cp_1d_z_p_type), DIMENSION(:, :), &
      POINTER                                :: delta_res, u_vec, z_vec

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store%rhoin_buffer),cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    NULLIFY(delta_res, u_vec, z_vec)
    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(mixing_store%rhoin_buffer(1,1)%cc)
    vol = rho%rho_g(1)%pw%pw_grid%vol
    alpha=mixing_store%alpha
    beta=mixing_store%beta


    wc = mixing_store%wc
    wmax = mixing_store%wmax
    nbuffer = mixing_store%nbuffer

    mixing_store%ncall = mixing_store%ncall + 1
    ib = MIN(mixing_store%ncall,nbuffer) 
    nb = MIN(mixing_store%ncall,nbuffer)
    ibb = MIN(mixing_store%ncall+1,nbuffer)

    ALLOCATE(res_rho(ng), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(res_rho_p(ng), STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(ib>1) THEN
      ALLOCATE(a(ib-1,ib-1), STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      a=0.0_dp
      ALLOCATE(b(ib-1,ib-1), STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      b=0.0_dp
      ALLOCATE(bq(ib-1,ib-1), STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      bq=0.0_dp
      ALLOCATE(tmp(ng), STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(delta_res_p(ng),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    p_metric => mixing_store%p_metric
    weight => mixing_store%weight
    CPPrecondition(ASSOCIATED(mixing_store%delta_res),cp_failure_level,routineP,error,failure)
    delta_res => mixing_store%delta_res
    CPPrecondition(ASSOCIATED(mixing_store%u_vec),cp_failure_level,routineP,error,failure)
    u_vec => mixing_store%u_vec
    CPPrecondition(ASSOCIATED(mixing_store%z_vec),cp_failure_level,routineP,error,failure)
    z_vec => mixing_store%z_vec


    delta_rhog = 0.0_dp
    delta_rhog_p = 0.0_dp

    DO ispin = 1,nspin
 
      fmat => mixing_store%fmat(:,:,ispin)
      smat => mixing_store%smat(:,:,ispin)
      gmat => mixing_store%gmat(:,:,ispin)

      delta = 0.0_dp
      delta_p = 0.0_dp
      ! Residual at this step R_i(G) (rho_out(G)-rho_in(G))
      ! Residual multiplied by the metrics RP_i(G) = (rho_out(G)-rho_in(G)) * P(G)
      ! Delta is the norm of the residual, measures how far we are from convergence
      DO ig = 1,ng
         res_rho(ig) = rho%rho_g(ispin)%pw%cc(ig) - mixing_store%rhoin_buffer(ib,ispin)%cc(ig)
         res_rho_p(ig) = res_rho(ig)*p_metric(ig)
         norm_ig = REAL(res_rho(ig),dp)*REAL(res_rho(ig),dp)+AIMAG(res_rho(ig))*AIMAG(res_rho(ig))
         delta = delta + norm_ig
         delta_p = delta_p +norm_ig * p_metric(ig)
      END DO
      CALL mp_sum(delta,para_env%group)
      delta = SQRT(delta)
      delta = delta
      CALL mp_sum(delta_p,para_env%group)
      delta_p = SQRT(delta_p)
      delta_p = delta_p

      delta_rhog = delta_rhog + delta
      delta_rhog_p = delta_rhog_p + delta_p

      weight(ib,ispin) = mixing_store%wc
      IF(wc< 0.0_dp) weight(ib,ispin) = 0.01_dp*ABS(wc)/(delta_p*delta_p)
      IF(weight(ib,ispin)==0.0_dp) weight(ib,ispin) = 100.0_dp
      IF(weight(ib,ispin)<1.0_dp) weight(ib,ispin) = 1.0_dp

      IF(ib==1) THEN
        ! Simple Kerker damping : linear mixing rho(G) = rho_in(G) - alpha k(G)*(rho_out(G)-rho_in(G))
        DO ig = 1,ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(ib,ispin)%cc(ig) + f_mix*res_rho(ig)
          mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO
      ELSE
        norm = 0.0_dp
        ! Difference of residuals DR_{i-1)} (G) = R_i(G) - R_{i-1}(G)
        DO ig = 1,ng
          delta_res(ib-1, ispin)%cc(ig) = res_rho(ig) - mixing_store%last_res(ispin)%cc(ig)
          delta_res_p(ig) = p_metric(ig) * (res_rho(ig) - mixing_store%last_res(ispin)%cc(ig))
          norm_ig = REAL(delta_res(ib-1, ispin)%cc(ig),dp)*REAL(delta_res(ib-1, ispin)%cc(ig),dp) + &
                    AIMAG(delta_res(ib-1, ispin)%cc(ig)) * AIMAG(delta_res(ib-1, ispin)%cc(ig)) 
          norm = norm + norm_ig*p_metric(ig) 
        END DO
        CALL mp_sum(norm, para_env%group)
        norm = 1.0_dp/SQRT(norm)
        delta_res(ib-1, ispin)%cc(:) = delta_res(ib-1, ispin)%cc(:)*norm 
        delta_res_p(:) = delta_res_p(:)*norm 

        ! Vector U_{i-1}(G) =  Drho_{i-1} + k(G)  * DR_{i-1}(G)
        DO ig = 1,ng
          tmp(ig) = (mixing_store%rhoin_buffer(ib,ispin)%cc(ig) - mixing_store%rhoin_buffer(ib-1,ispin)%cc(ig))*norm 
          u_vec(ib-1,ispin)%cc(ig) = (tmp(ig) + mixing_store%kerker_factor(ig) * delta_res(ib-1, ispin)%cc(ig))
        END DO

        DO jb = 1,ib-1
          fmat(jb,ib-1) = 0.0_dp 
          smat(jb,ib-1) = 0.0_dp 
          gmat(jb,ib-1) = 0.0_dp 
          gmat(ib-1,jb) = 0.0_dp 

          DO ig = 1,ng
            rep_j = REAL(delta_res(jb, ispin)%cc(ig),dp) 
            imp_j = AIMAG(delta_res(jb, ispin)%cc(ig)) 
           ! < DR_{j} | DR_{i-1} >
            rep = REAL(delta_res_p(ig),dp)
            imp = AIMAG(delta_res_p(ig))
            fmat(jb,ib-1) = fmat(jb,ib-1) + rep_j*rep + imp_j*imp

           ! < DR_{j} | K | DR_{i-1} >
            rep = REAL(delta_res(ib-1,ispin)%cc(ig),dp)*mixing_store%kerker_factor(ig)
            imp = AIMAG(delta_res(ib-1,ispin)%cc(ig))*mixing_store%kerker_factor(ig)
            smat(jb,ib-1) = smat(jb,ib-1) + rep_j*rep + imp_j*imp

           ! < DR_{j} | U_{i-1} > = < DR_{j} | J_{i-1} - K | DR_{i-i} >
            rep = REAL(u_vec(ib-1,ispin)%cc(ig),dp)
            imp = AIMAG(u_vec(ib-1,ispin)%cc(ig))
            gmat(jb,ib-1) = gmat(jb,ib-1) - rep_j*rep + imp_j*imp

         END DO

         gmat(ib-1,jb) = 0.0_dp
         ! < U_{j} | DR_{i-1} >
         DO ig = 1,ng
            rep_j = REAL(u_vec(jb,ispin)%cc(ig),dp)
            imp_j = AIMAG(u_vec(jb,ispin)%cc(ig))
            rep = REAL(delta_res(ib-1,ispin)%cc(ig),dp)
            imp = AIMAG(delta_res(ib-1,ispin)%cc(ig))
            gmat(ib-1,jb) = gmat(ib-1,jb) - rep_j*rep + imp_j*imp
          END DO

        END DO
        CALL mp_sum(fmat(1:ib-1,ib-1),para_env%group) 
        CALL mp_sum(smat(1:ib-1,ib-1),para_env%group) 
        CALL mp_sum(gmat(1:ib-1,ib-1),para_env%group) 
        CALL mp_sum(gmat(ib-1,1:ib-2),para_env%group)

        fmat(ib-1,ib-1) = 1.0_dp

        DO jb = 1,ib-2
          fmat(ib-1,jb) = fmat(jb,ib-1) 
          smat(ib-1,jb) = smat(jb,ib-1)
        ENDDO

        DO jb = 1,ib-1
          a(jb,jb) = 1.0_dp + weight(jb,ispin)*weight(jb,ispin)*fmat(jb,jb)          
          DO kb = 1,jb-1
             a(jb,kb) = weight(jb,ispin)*weight(kb,ispin)*fmat(jb,kb)
             a(kb,jb) = weight(jb,ispin)*weight(kb,ispin)*fmat(kb,jb)
          ENDDO
        END DO

        IF(.TRUE.) THEN
          b=0.0_dp
          CALL invert_matrix(a,b,inv_err,error=error)
        ELSE
          b=0.0_dp
          ALLOCATE(work(ib-1,ib-1),aval(ib-1),av(ib-1,ib-1),au(ib-1,ib-1))
          work=a
          ALLOCATE(iwork(8*(ib-1)),work_dgesdd(1))
             lwork=-1
             CALL DGESDD('S',ib-1,ib-1,work,ib-1,aval,au,ib-1,av,ib-1,work_dgesdd,lwork,iwork,info)
             lwork=work_dgesdd(1)
             DEALLOCATE(work_dgesdd) ; ALLOCATE(work_dgesdd(lwork))
             CALL DGESDD('S',ib-1,ib-1,work,ib-1,aval,au,ib-1,av,ib-1,work_dgesdd,lwork,iwork,info)
             ! construct the inverse
             DO kb=1,ib-1
                ! invert SV
                IF (aval(kb)<1.E-6_dp) THEN
                   aval(kb)=0.0_dp
                ELSE
                   aval(kb)=1.0_dp/aval(kb)
                ENDIF
                av(kb,:)=av(kb,:)*aval(kb)
             ENDDO
             CALL DGEMM('T','T',ib-1,ib-1,ib-1,1.0_dp,av,ib-1,au,ib-1,0.0_dp,b,ib-1)
             DEALLOCATE(iwork,aval,au,av,work_dgesdd, work)
        END IF


        bq = 0.0_dp
        ! Broyden second method requires also bq (NYI) 
        skip_bq = .TRUE.
        DO jb = 1,ib-1
          skip_bq = skip_bq .AND. (weight(jb,ispin)>mixing_store%wmax) 
          DO kb = 1, ib - 2
            bq(jb,kb) = 0.0_dp
            DO kkb = 1, ib-1
              bq(jb,kb) = bq(jb,kb) - weight(jb,ispin)*weight(kkb,ispin)*b(jb,kkb)*fmat(kkb,kb)
            END DO
          END DO
          bq(jb,jb) = 1.0_dp + bq(jb,jb)
        END DO

        IF(.NOT. skip_bq) THEN
           ! in this case the old z_vec is needed 
           ! a temporary array is needed to store the new one
           ALLOCATE(tmp_z(ib-1),STAT=istat)
           DO jb = 1,ib-1
             ALLOCATE(tmp_z(jb)%cc(ng),STAT=istat)
           END DO 
        END IF
        DO jb = 1,ib-1
          tmp(:) = CMPLX(0.0_dp,0.0_dp)
          IF(.NOT. skip_bq) THEN
            ! sum_{kb} bq(jb,kb) * z_vec_{kb,iter-2}
            ! added only for small weights
            DO kb = 1,ib-2
               IF(weight(kb,ispin) >=(10._dp*mixing_store%wmax))  CYCLE
               DO ig = 1,ng
                 tmp(ig) = tmp(ig) + bq(jb,kb)*z_vec(kb,ispin)%cc(ig)
               END DO
            END DO  ! kb
          END IF

           ! sum_{kb} w(jb)*w(kb)*b(jb,kb) * u_vec_{kb}
          DO kb = 1,ib-1
             DO ig = 1,ng
                tmp(ig) = tmp(ig) + weight(kb,ispin)*weight(jb,ispin)*b(jb,kb)*u_vec(kb,ispin)%cc(ig) 
             END DO
          END DO

          ! store the new z_vec(jb)
          IF(skip_bq .OR. (weight(jb,ispin) >=(10._dp*mixing_store%wmax))) THEN
             z_vec(jb,ispin)%cc(:) = tmp(:)
          ELSE
             ! temporary array: old z_vec may still be needed
             tmp_z(jb)%cc(:) = tmp(:) 
          END IF
        END DO !jb

        IF(.NOT. skip_bq ) THEN
          DO jb = 1,ib-1
             IF (weight(jb,ispin) < (10._dp*mixing_store%wmax))  z_vec(jb,ispin)%cc(:) = tmp_z(jb)%cc(:)
             DEALLOCATE(tmp_z(jb)%cc,STAT=istat)
          END DO
          DEALLOCATE(tmp_z,STAT=istat)
        END IF

        ! Overwrite the density i reciprocal space
        rho%rho_g(ispin)%pw%cc(:) = CMPLX(0.0_dp,0.0_dp)
        DO jb =1,ib-1
          norm = 0.0_dp
          DO ig= 1,ng
            rep_j = REAL(delta_res(jb,ispin)%cc(ig),dp)
            imp_j = AIMAG(delta_res(jb,ispin)%cc(ig))
            rep = REAL(res_rho_p(ig),dp)
            imp = AIMAG(res_rho_p(ig))
            norm = norm + rep_j*rep + imp_j*imp
          END DO
          CALL mp_sum(norm,para_env%group)
          ! Subtract |Z_jb)><DR_jb|P|R_{iter}> 
          DO ig = 1,ng
            rho%rho_g(ispin)%pw%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) - norm * z_vec(jb,ispin)%cc(ig) 
          END DO
        END DO

        DO ig = 1,ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) =rho%rho_g(ispin)%pw%cc(ig) +&
                   mixing_store%rhoin_buffer(ib,ispin)%cc(ig) + f_mix*res_rho(ig)
          mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO
      END IF  ! ib

      mixing_store%last_res(ispin)%cc(:) = res_rho(:)
      IF(mixing_store%ncall>=nbuffer) THEN
          DO jb = 1,nbuffer - 2
             mixing_store%rhoin_buffer(jb,ispin)%cc(:) = mixing_store%rhoin_buffer(jb+1,ispin)%cc(:)
             mixing_store%delta_res(jb,ispin)%cc(:) = mixing_store%delta_res(jb+1,ispin)%cc(:)
             mixing_store%u_vec(jb,ispin)%cc(:) = mixing_store%u_vec(jb+1,ispin)%cc(:)
             mixing_store%z_vec(jb,ispin)%cc(:) = mixing_store%z_vec(jb+1,ispin)%cc(:)
             DO kb = 1, nbuffer - 2
                mixing_store%fmat(kb,jb,ispin) = mixing_store%fmat(kb+1,jb+1,ispin) 
                mixing_store%smat(kb,jb,ispin) = mixing_store%smat(kb+1,jb+1,ispin) 
                mixing_store%gmat(kb,jb,ispin) = mixing_store%gmat(kb+1,jb+1,ispin) 
             END DO
          END DO 
          mixing_store%rhoin_buffer(nbuffer-1,ispin)%cc(:) =&
               mixing_store%rhoin_buffer(nbuffer,ispin)%cc(:)
          mixing_store%weight(nbuffer-1,ispin) = mixing_store%weight(nbuffer,ispin)
       END IF

    END DO  ! ispin
    IF(ib>1) THEN
      DEALLOCATE(a,b,bq, STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(delta_res_p,tmp, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(res_rho,res_rho_p, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

 END SUBROUTINE broyden_mixing_new

! *****************************************************************************
!> \brief Multisecant scheme to perform charge density Mixing 
!>        as preconditioner we use the Kerer damping factor
!>        The mixing is applied directly on the density in reciprocal space,
!>        therefore it affects the potentials
!>        on the grid but not the density matrix
!> \par History
!>      03.2009 
!> \author MI
! *****************************************************************************
 SUBROUTINE multisecant_mixing(mixing_store,rho,para_env,error)

    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'multisecant_mixing', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: cmone = (-1.0_dp,0.0_dp), &
                                                cone = (1.0_dp,0.0_dp), &
                                                czero = (0.0_dp,0.0_dp)

    COMPLEX(dp)                              :: saa, yaa
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:)   :: gn_global, pgn, &
                                                res_matrix_global, tmp_vec, &
                                                ugn
    COMPLEX(dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a_matrix, res_matrix, sa, &
                                                step_matrix, ya
    COMPLEX(dp), DIMENSION(:), POINTER       :: gn
    INTEGER                                  :: handle, ib, ib_next, ib_prev, &
                                                ig, ig_global, iig, ispin, &
                                                istat, jb, kb, nb, nbuffer, &
                                                ng, ng_global, nspin
    LOGICAL                                  :: failure, use_zgemm, &
                                                use_zgemm_rev
    REAL(dp) :: alpha, f_mix, gn_norm, gn_norm_old, inv_err, pgn_norm, prec, &
      r_step, reg_par, sigma_max, sigma_tilde, step_size
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: norm_res
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: b_matrix, binv_matrix
    REAL(dp), SAVE                           :: sigma_old = 1.0_dp

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store),cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    NULLIFY(gn)

    use_zgemm = .FALSE.
    use_zgemm_rev = .TRUE.

    ! prepare the parameters
    nspin = SIZE(rho%rho_g,1)
    ng_global = rho%rho_g(1)%pw%pw_grid%ngpts
    ng = SIZE(mixing_store%rhoin_buffer(1,1)%cc)
    alpha=mixing_store%alpha
  
    sigma_max = mixing_store%sigma_max
    reg_par = mixing_store%reg_par
    r_step  =  mixing_store%r_step
    nbuffer = mixing_store%nbuffer

    ! determine the step number, and multisecant iteration 
    nb  = MIN(mixing_store%ncall,nbuffer-1)
    ib  = MODULO(mixing_store%ncall,nbuffer) +1 
    IF(mixing_store%ncall>0) THEN
      ib_prev = MODULO(mixing_store%ncall-1,nbuffer) +1 
    ELSE
      ib_prev = 0 
    END IF
    mixing_store%ncall = mixing_store%ncall + 1
    ib_next = MODULO(mixing_store%ncall,nbuffer) +1
   
    ! compute the residual gn and its norm gn_norm
    DO ispin = 1,nspin
      gn =>  mixing_store%res_buffer(ib,ispin)%cc
      gn_norm = 0.0_dp
      DO ig = 1,ng 
        gn(ig) = (rho%rho_g(ispin)%pw%cc(ig)-mixing_store%rhoin_buffer(ib,ispin)%cc(ig))
        gn_norm = gn_norm + &
           REAL(gn(ig),dp)*REAL(gn(ig),dp) + AIMAG( gn(ig))*AIMAG( gn(ig))
      END DO
      CALL mp_sum(gn_norm,para_env%group)
      gn_norm = SQRT(gn_norm)
      mixing_store%norm_res_buffer(ib,ispin) = gn_norm
    END DO 

    IF( nb == 0) THEN
      DO ispin = 1,nspin
        DO ig = 1,ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(1,ispin)%cc(ig) + &
              f_mix*  mixing_store%res_buffer(1,ispin)%cc(ig)
          mixing_store%rhoin_buffer(ib_next,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO
      END DO
      CALL timestop(handle)
      RETURN
    END IF

    ! allocate temporary arrays
    ! step_matrix  S ngxnb
    ALLOCATE(step_matrix(ng,nb),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ! res_matrix Y  ngxnb
    ALLOCATE(res_matrix(ng,nb),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ! matrix A  nbxnb
    ALLOCATE(a_matrix(nb,ng_global),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ! PSI nb vector of norms
    ALLOCATE(norm_res(nb),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    ! matrix B   nbxnb
    ALLOCATE(b_matrix(nb,nb),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ! matrix B_inv   nbxnb
    ALLOCATE(binv_matrix(nb,nb),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(gn_global(ng_global),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(res_matrix_global(ng_global),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(use_zgemm) THEN
      ALLOCATE(sa(ng,ng_global),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(ya(ng,ng_global),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(use_zgemm_rev) THEN
      ALLOCATE(tmp_vec(nb),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    ALLOCATE(pgn(ng),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(ugn(ng),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    
    DO ispin = 1,nspin
     ! generate the global vector with the present residual
      gn =>  mixing_store%res_buffer(ib,ispin)%cc
      gn_global = CMPLX(0.0_dp, 0.0_dp)
      DO ig = 1,ng
         ig_global = mixing_store%ig_global_index(ig) 
         gn_global(ig_global) = gn(ig)
      END DO
      CALL mp_sum(gn_global,para_env%group)

    ! compute steps (matrix S) and residual differences (matrix Y) 
    ! with respect to the present
    ! step and the present residual (use stored rho_in and res_buffer)

    ! These quantities are pre-conditioned by means of Kerker factor multipliccation
      DO jb = 1,nb+1
        IF(jb<ib) THEN
           kb = jb
        ELSEIF(jb>ib) THEN
           kb = jb-1
        ELSE
           CYCLE
        END IF
        norm_res(kb) = 0.0_dp
        DO ig = 1,ng

          prec = mixing_store%kerker_factor(ig)

          step_matrix(ig,kb) = prec*(mixing_store%rhoin_buffer(jb,ispin)%cc(ig) - &
                               mixing_store%rhoin_buffer(ib,ispin)%cc(ig))
          res_matrix(ig,kb) = (mixing_store%res_buffer(jb,ispin)%cc(ig) -&
                             mixing_store%res_buffer(ib,ispin)%cc(ig))
          norm_res(kb)= norm_res(kb) + REAL(res_matrix(ig,kb),dp)*REAL(res_matrix(ig,kb),dp) +&
                        AIMAG(res_matrix(ig,kb))* AIMAG(res_matrix(ig,kb))
          res_matrix(ig,kb) = prec * res_matrix(ig,kb)
        END DO
      END DO  !jb

    ! normalize each column of S and Y => Snorm Ynorm
      CALL mp_sum(norm_res,para_env%group)
      norm_res(1:nb) = 1.0_dp/SQRT(norm_res(1:nb))
      DO kb = 1,nb
         step_matrix(1:ng,kb) = step_matrix(1:ng,kb)*norm_res(kb)
         res_matrix(1:ng,kb) = res_matrix(1:ng,kb)*norm_res(kb)
      END DO

    ! compute A as [(Ynorm^t Ynorm) + (alpha I)]^(-1) Ynorm^t
      ! compute B
      DO jb = 1,nb
           DO kb = 1,nb
             b_matrix(kb,jb) = 0.0_dp
             DO ig = 1,ng
               ! it is assumed that summing over all G vector gives a real, because
               !  y(-G,kb) = (y(G,kb))*
               b_matrix(kb,jb) = b_matrix(kb,jb) + REAL(res_matrix(ig,kb)*res_matrix(ig,jb),dp) 
             END DO
           END DO
      END DO

      CALL mp_sum(b_matrix,para_env%group)
      DO jb= 1,nb
         b_matrix(jb,jb) = b_matrix(jb,jb) + reg_par
      END DO
      ! invert B
      CALL invert_matrix(b_matrix,binv_matrix,inv_err,error=error)

      ! A = Binv Ynorm^t 
      a_matrix = CMPLX(0.0_dp,0.0_dp)
      DO kb = 1,nb
         res_matrix_global = CMPLX(0.0_dp, 0.0_dp)
         DO ig = 1,ng
            ig_global = mixing_store%ig_global_index(ig) 
            res_matrix_global(ig_global) = res_matrix(ig,kb)
         END DO
         CALL mp_sum(res_matrix_global,para_env%group)

         DO jb = 1,nb
            DO ig = 1,ng 
               ig_global = mixing_store%ig_global_index(ig)
               prec = mixing_store%kerker_factor(ig)
               a_matrix(jb,ig_global) = a_matrix(jb,ig_global) + &
                  prec*binv_matrix(jb,kb)*res_matrix_global(ig_global)
            END DO
         END DO
      END DO
      CALL mp_sum(a_matrix,para_env%group)

      ! compute the two components of gn that will be used to update rho
      gn =>  mixing_store%res_buffer(ib,ispin)%cc
      pgn_norm = 0.0_dp

      IF(use_zgemm) THEN

        CALL zgemm("N","N",ng,ng_global,nb,cmone,step_matrix(1,1),ng,&
                   a_matrix(1,1),nb,czero,sa(1,1),ng)
        CALL zgemm("N","N",ng,ng_global,nb,cmone,res_matrix(1,1),ng,&
                   a_matrix(1,1),nb,czero,ya(1,1),ng)
        DO ig = 1,ng
           ig_global = mixing_store%ig_global_index(ig) 
           ya(ig,ig_global) =  ya(ig,ig_global) + CMPLX(1.0_dp,0.0_dp)
        END DO

        CALL zgemv("N",ng,ng_global,cone,sa(1,1),&
                    ng,gn_global(1),1,czero,pgn(1),1)
        CALL zgemv("N",ng,ng_global,cone,ya(1,1),&
                    ng,gn_global(1),1,czero,ugn(1),1)

        DO ig = 1,ng
          pgn_norm = pgn_norm + REAL(pgn(ig),dp)*REAL(pgn(ig),dp) + &
                  AIMAG(pgn(ig))*AIMAG(pgn(ig))
        END DO
        CALL mp_sum(pgn_norm,para_env%group)
      ELSEIF(use_zgemm_rev)  THEN

        CALL zgemv("N",nb,ng_global,cone,a_matrix(1,1),&
                     nb,gn_global(1),1,czero,tmp_vec(1),1)

        CALL zgemv("N",ng,nb,cmone,step_matrix(1,1),ng,&
                     tmp_vec(1),1,czero,pgn(1),1)

        CALL zgemv("N",ng,nb,cmone,res_matrix(1,1),ng,&
                     tmp_vec(1),1,czero,ugn(1),1)

        DO ig = 1,ng
          pgn_norm = pgn_norm + REAL(pgn(ig),dp)*REAL(pgn(ig),dp) + &
                  AIMAG(pgn(ig))*AIMAG(pgn(ig))
          ugn(ig) = ugn(ig) + gn(ig)
        END DO
        CALL mp_sum(pgn_norm,para_env%group)

      ELSE
        DO ig = 1,ng
          pgn(ig) =   CMPLX(0.0_dp,0.0_dp)
          ugn(ig) =   CMPLX(0.0_dp,0.0_dp)
          ig_global =  mixing_store%ig_global_index(ig)
          DO iig = 1,ng_global
             saa = CMPLX(0.0_dp,0.0_dp)
             yaa = CMPLX(0.0_dp,0.0_dp)

            IF(ig_global==iig) yaa =  CMPLX(1.0_dp,0.0_dp)
  
            DO jb = 1,nb
               saa = saa - step_matrix(ig,jb)*a_matrix(jb,iig)
               yaa = yaa - res_matrix(ig,jb)*a_matrix(jb,iig)
            END DO
            pgn(ig) = pgn(ig) + saa * gn_global(iig) 
            ugn(ig) = ugn(ig) + yaa * gn_global(iig) 
          END DO
        END DO
        DO ig = 1,ng
          pgn_norm = pgn_norm + REAL(pgn(ig),dp)*REAL(pgn(ig),dp) + &
                  AIMAG(pgn(ig))*AIMAG(pgn(ig))
        END DO
        CALL mp_sum(pgn_norm,para_env%group)
      END IF

      gn_norm =  mixing_store%norm_res_buffer(ib,ispin)
      gn_norm_old =  mixing_store%norm_res_buffer(ib_prev,ispin)

      IF(ib_prev/=0) THEN
        sigma_tilde = sigma_old * MAX(0.5_dp,MIN(2.0_dp,gn_norm_old/gn_norm))
      ELSE
        sigma_tilde = 0.5_dp
      END IF
      ! Step size for the unpredicted component
      step_size = MIN(sigma_tilde,r_step*pgn_norm/gn_norm,sigma_max)

      sigma_old = 0.1_dp! step_size

    ! update the density
      DO ig = 1,ng
         rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(ib,ispin)%cc(ig) -&
              step_size*ugn(ig) + alpha*pgn(ig)
         mixing_store%rhoin_buffer(ib_next,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) 
      END DO

    END DO ! ispin

    ! Deallocate  temporary arrays
    DEALLOCATE(step_matrix,res_matrix, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(norm_res, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(a_matrix,b_matrix,binv_matrix, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(ugn,pgn, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    IF(use_zgemm) THEN
      DEALLOCATE(sa,ya, STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(use_zgemm_rev) THEN
      DEALLOCATE(tmp_vec, STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(gn_global,res_matrix_global, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

 END SUBROUTINE multisecant_mixing

END MODULE qs_gspace_mixing
