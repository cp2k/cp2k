!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE qs_gspace_mixing

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_create, &
       cp_dbcsr_distribution, cp_dbcsr_get_info, cp_dbcsr_init, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_set, cp_dbcsr_trace
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_alloc_block_from_nbl,&
                                             cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_trace,&
                                             cp_fm_upper_to_full
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_types,                     ONLY: dbcsr_type_symmetric
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all,&
                                             invert_matrix
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_min,&
                                             mp_sum
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_integrate_function,&
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
  USE qs_density_mixing_types,         ONLY: broyden_mixing_new_nr,&
                                             broyden_mixing_nr,&
                                             cp_1d_z_p_type,&
                                             gspace_mixing_nr,&
                                             mixing_storage_type,&
                                             multisecant_mixing_nr,&
                                             pulay_mixing_nr
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_methods,                  ONLY: cp_sm_mix
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gspace_mixing'

  PUBLIC :: gspace_mixing, broyden_mixing, broyden_mixing_new, gmix_potential_only, &
            gspace_density_mixing, mixing_allocate, mixing_init,&
            pulay_mixing, self_consistency_check

CONTAINS

  SUBROUTINE self_consistency_check(rho_ao,p_delta,para_env,p_out,delta,error)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao, p_delta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: p_out
    REAL(KIND=dp), INTENT(INOUT)             :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'self_consistency_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspin
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_q, p_in

    CALL timeset(routineN,handle)

    NULLIFY(matrix_q, p_in)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(p_out),cp_failure_level,routineP,error,failure)
    NULLIFY(matrix_q, p_in)
    p_in => rho_ao
    matrix_q => p_delta
    nspin = SIZE(p_in)

    ! Compute the difference (p_out - p_in)and check convergence
    delta = 0.0_dp
    DO ispin = 1,nspin
       CALL cp_dbcsr_set(matrix_q(ispin)%matrix, 0.0_dp, error=error)
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
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
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
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_q, matrix_s, p_in
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_p_type)                          :: kerker_pot_rspace, rhog
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)

    NULLIFY(matrix_q,p_in)
    NULLIFY(auxbas_pw_pool, matrix_s, rho, para_env, pw_env)
    NULLIFY(g2, pw_grid)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(p_out),cp_failure_level,routineP,error,failure)

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

    CALL cp_dbcsr_get_info(matrix= p_in(1)%matrix, nfullrows_total=nao)

    ! Compute the difference (p_out - p_in)and check convergence
    delta = 0.0_dp
    DO ispin = 1,nspin
       CALL cp_dbcsr_set(matrix_q(ispin)%matrix, 0.0_dp,error=error)
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

       CALL pw_zero(rhog%pw, error=error)
       ! point pw
       pw_grid => auxbas_pw_pool % pw_grid
       ng = SIZE ( pw_grid % gsq )
       g2 => pw_grid %gsq
       CALL pw_transfer ( rho%rho_g(ispin)%pw, rhog%pw, error=error)
!       rhog%pw =>  rho%rho_g(ispin)%pw
       DO ig=1,ng
!          rhog%pw%cc(ig) = rhog%pw%cc(ig) *  (1.0_dp-a_mix)*g2(ig)/(g2(ig)+b_mix*b_mix)
          rhog%pw%cc(ig) = rhog%pw%cc(ig) *  g2(ig)/(g2(ig)+b_mix*b_mix)
       ENDDO
       CALL pw_transfer ( rhog%pw, kerker_pot_rspace%pw, error=error)

       CALL pw_scale ( kerker_pot_rspace%pw, pw_grid%dvol, error=error)

       CALL cp_dbcsr_set(matrix_q(ispin)%matrix, 0.0_dp,error=error)
       CALL integrate_v_rspace(v_rspace=kerker_pot_rspace, &
            h=matrix_q(ispin), qs_env=qs_env,&
            calculate_forces =.FALSE., gapw=gapw, error=error)

       ! Impose Tr[S^-1 Q S^-1 S] = 0, i.e. the correction to P_out
       ! should not change the total number of electrons
    IF(.FALSE.) THEN
       CALL copy_dbcsr_to_fm(matrix_q(ispin)%matrix,scf_env%scf_work2,error=error)
       CALL cp_fm_upper_to_full(scf_env%scf_work2,scf_env%scf_work1(ispin)%matrix,error=error)
       CALL cp_fm_trace(scf_env%s_minus_one,scf_env%scf_work2,trace,error=error)
       trace = -trace / REAL(nao,KIND=dp)
       CALL cp_dbcsr_add(matrix_a=matrix_q(ispin)%matrix,matrix_b=matrix_s(1)%matrix,&
            alpha_scalar=1.0_dp,beta_scalar=trace,error=error)
    END IF

       CALL cp_dbcsr_trace(p_in(ispin)%matrix,matrix_s(1)%matrix,trace,error=error)
       t1 = trace

       CALL cp_dbcsr_sm_fm_multiply(matrix_q(ispin)%matrix,scf_env%s_minus_one,&
            scf_env%scf_work2, nao,error=error)
       CALL cp_fm_gemm("N","N",nao,nao,nao,1.0_dp,scf_env%s_minus_one,&
            scf_env%scf_work2,0.0_dp,scf_env%scf_work1(ispin)%matrix,error=error )

!       a_mix = 1.0_dp-scf_env%p_mix_alpha
!       CALL copy_fm_to_sm(scf_env%scf_work1(ispin)%matrix,&
!            p_out(ispin)%matrix,a_mix,1.0_dp)

       a_mix = scf_env%p_mix_alpha
       CALL copy_fm_to_dbcsr(scf_env%scf_work1(ispin)%matrix,&
            p_out(ispin)%matrix,a_mix,0.0_dp,keep_sparsity=.TRUE., error=error)
       CALL cp_dbcsr_add(p_out(ispin)%matrix,p_in(ispin)%matrix,&
            alpha_scalar=1.0_dp,beta_scalar=1.0_dp,error=error)
       CALL cp_dbcsr_trace(p_out(ispin)%matrix,matrix_s(1)%matrix,trace, error=error)
       t2 =t1/ trace
       CALL cp_dbcsr_scale(p_out(ispin)%matrix,t2,error=error)

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
  SUBROUTINE gspace_mixing(qs_env, scf_env, mixing_store, rho, para_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gspace_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, ig, ispin, &
                                                natom, ng, nspin
    LOGICAL                                  :: failure, gapw
    REAL(dp)                                 :: alpha
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_tmp
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store),cp_failure_level,routineP,error,failure)

    NULLIFY(auxbas_pw_pool, dft_control, pw_env, rho_atom)
    CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env,error=error)
    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(rho%rho_g(1)%pw%pw_grid%gsq)
    CPPrecondition((ng == SIZE(mixing_store%rhoin(1)%cc)),cp_failure_level,routineP,error,failure)

    alpha = mixing_store%alpha
    gapw = dft_control%qs_control%gapw

    IF(nspin==2) THEN
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool,&
                           rho_tmp%pw,&
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE,error=error)
      CALL pw_zero(rho_tmp%pw, error=error)
      CALL pw_copy( rho%rho_g(1)%pw, rho_tmp%pw , error=error)
      CALL pw_axpy(rho%rho_g(2)%pw,rho%rho_g(1)%pw, 1.0_dp, error=error)
      CALL pw_axpy(rho_tmp%pw,rho%rho_g(2)%pw, -1.0_dp, error=error)
      CALL pw_scale(rho%rho_g(2)%pw, -1.0_dp, error=error)
    END IF

    IF(scf_env%iter_count+1 <= mixing_store%nskip_mixing) THEN
      ! skip mixing
      DO ispin = 1,nspin
        DO ig = 1,ng
!          mixing_store%rhoin_buffer(1,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
          mixing_store%rhoin(ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO
      END DO
      IF(gapw) THEN
        CALL get_qs_env(qs_env=qs_env,&
              rho_atom_set=rho_atom,error=error)
        natom = SIZE(rho_atom)
        DO ispin = 1, nspin
          DO iatom = 1,natom
            IF(mixing_store%paw(iatom)) THEN
               mixing_store%cpc_h_in(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
               mixing_store%cpc_s_in(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
            END IF
          END DO
        END DO
      END IF
      mixing_store%iter_method= "NoMix"
      IF(nspin==2) THEN
        CALL pw_axpy(rho%rho_g(2)%pw,rho%rho_g(1)%pw, 1.0_dp, error=error)
        CALL pw_scale(rho%rho_g(1)%pw, 0.5_dp, error=error)
        CALL pw_axpy(rho%rho_g(1)%pw,rho%rho_g(2)%pw, -1.0_dp, error=error)
        CALL pw_scale(rho%rho_g(2)%pw,-1.0_dp, error=error)
        CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tmp%pw,error=error)
      END IF
      RETURN
    END IF

    CALL timeset(routineN,handle)
    IF((scf_env%iter_count+1 - mixing_store%nskip_mixing) <= mixing_store%n_simple_mix) THEN
       CALL gmix_potential_only(qs_env,mixing_store, rho, error)
       mixing_store%iter_method= "Kerker"
    ELSE

       IF(scf_env%mixing_method==gspace_mixing_nr) THEN
          CALL gmix_potential_only(qs_env,mixing_store, rho, error)
          mixing_store%iter_method= "Kerker"
       ELSEIF(scf_env%mixing_method==pulay_mixing_nr) THEN
          CPPrecondition(.NOT. gapw,cp_failure_level,routineP,error,failure)
          CALL pulay_mixing(mixing_store, rho, para_env, error)
          mixing_store%iter_method= "Pulay"
       ELSEIF(scf_env%mixing_method==broyden_mixing_nr) THEN
          CALL broyden_mixing(qs_env,mixing_store, rho, para_env, error)
          mixing_store%iter_method= "Broy."
       ELSEIF(scf_env%mixing_method==broyden_mixing_new_nr) THEN
          CPPrecondition(.NOT. gapw,cp_failure_level,routineP,error,failure)
          CALL broyden_mixing_new(mixing_store, rho, para_env, error)
          mixing_store%iter_method= "Broy."
       ELSEIF(scf_env%mixing_method==multisecant_mixing_nr) THEN
          CPPrecondition(.NOT. gapw,cp_failure_level,routineP,error,failure)
          CALL multisecant_mixing(mixing_store, rho, para_env, error)
          mixing_store%iter_method= "MSec."
       END IF
    END IF

    IF(nspin==2) THEN
        CALL pw_axpy(rho%rho_g(2)%pw,rho%rho_g(1)%pw, 1.0_dp, error=error)
        CALL pw_scale(rho%rho_g(1)%pw, 0.5_dp, error=error)
        CALL pw_axpy(rho%rho_g(1)%pw,rho%rho_g(2)%pw, -1.0_dp, error=error)
        CALL pw_scale(rho%rho_g(2)%pw,-1.0_dp, error=error)
      CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tmp%pw,error=error)
    END IF


    DO ispin=1,nspin
      CALL pw_transfer(rho%rho_g(ispin)%pw,rho%rho_r(ispin)%pw, error=error)
      rho%tot_rho_r(ispin) = pw_integrate_function(rho%rho_r(ispin)%pw,isign=-1,error=error)
!dbg
!  write(*,*) 'rho int 4', ispin , rho%tot_rho_r(ispin)
!dbg
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

  SUBROUTINE gmix_potential_only (qs_env,mixing_store, rho, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gmix_potential_only', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), DIMENSION(:), POINTER       :: cc_new
    INTEGER                                  :: handle, iatom, ig, ispin, &
                                                natom, ng, nspin
    LOGICAL                                  :: failure, gapw
    REAL(dp)                                 :: alpha, f_mix
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store%rhoin),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store%kerker_factor),cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    NULLIFY(cc_new)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(rho%rho_g(1)%pw%pw_grid%gsq)

    gapw = qs_env%dft_control%qs_control%gapw
    alpha = mixing_store%alpha

    DO ispin = 1,nspin
      cc_new =>  rho%rho_g(ispin)%pw%cc
      DO ig = 1,mixing_store%ig_max ! ng
         f_mix = mixing_store%alpha*mixing_store%kerker_factor(ig)
         cc_new(ig) = (1.0_dp-f_mix)*mixing_store%rhoin(ispin)%cc(ig)+f_mix*cc_new(ig)
         mixing_store%rhoin(ispin)%cc(ig) = cc_new(ig)
      END DO
      DO ig = mixing_store%ig_max +1 , ng
         f_mix = mixing_store%alpha
         cc_new(ig) = (1.0_dp-f_mix)*mixing_store%rhoin(ispin)%cc(ig)+f_mix*cc_new(ig)
         mixing_store%rhoin(ispin)%cc(ig) = cc_new(ig)
      END DO
    END DO

    IF(gapw) THEN
      CALL get_qs_env(qs_env=qs_env,&
            rho_atom_set=rho_atom,error=error)
      natom = SIZE(rho_atom)
      DO ispin = 1, nspin
        DO iatom = 1,natom
          IF(mixing_store%paw(iatom)) THEN
            rho_atom(iatom)%cpc_h(ispin)%r_coef = alpha*rho_atom(iatom)%cpc_h(ispin)%r_coef + &
                   mixing_store%cpc_h_in(iatom,ispin)%r_coef*(1._dp-alpha)
            rho_atom(iatom)%cpc_s(ispin)%r_coef = alpha*rho_atom(iatom)%cpc_s(ispin)%r_coef + &
                   mixing_store%cpc_s_in(iatom,ispin)%r_coef*(1._dp-alpha)
            mixing_store%cpc_h_in(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
            mixing_store%cpc_s_in(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
          END IF
        END DO
      END DO
    END IF

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
      IF(nb==1) mixing_store%rhoin_buffer(1,ispin)%cc = mixing_store%rhoin(ispin)%cc
      DO ig = 1,ng
          f_mix = mixing_store%kerker_factor(ig)
 !         f_mix = 1.0_dp
          mixing_store%res_buffer(ib,ispin)%cc(ig) = f_mix*(rho%rho_g(ispin)%pw%cc(ig)-&
                              mixing_store%rhoin_buffer(ib,ispin)%cc(ig))
      END DO

      DO jb = 1,nb
         mixing_store%pulay_matrix(jb,ib) = 0.0_dp
         DO ig = 1,ng
!           f_mix = 1.0_dp*mixing_store%kerker_factor(ig)
           f_mix = mixing_store%special_metric(ig)
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
                        mixing_store%pulay_beta*mixing_store%res_buffer(jb,ispin)%cc(ig))
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
!> \brief Broyden Mixing using as metrics for the residual the Kerer damping factor
!>        The mixing is applied directly on the density expansion in reciprocal space
!> \par History
!>      03.2009
!> \author MI
! *****************************************************************************

  SUBROUTINE broyden_mixing(qs_env,mixing_store, rho, para_env, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'broyden_mixing', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp)                              :: cc_mix
    COMPLEX(dp), ALLOCATABLE, DIMENSION(:)   :: res_rho
    INTEGER                                  :: handle, i, iatom, ib, ibb, &
                                                ig, ispin, istat, j, jb, kb, &
                                                n1, n2, natom, nb, nbuffer, &
                                                ng, nspin
    LOGICAL                                  :: failure, gapw
    REAL(dp)                                 :: alpha, beta, dcpc_h_res, &
                                                dcpc_s_res, delta_norm, &
                                                f_mix, inv_err, res_norm, &
                                                valh, vals, w0
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: c, g
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a, b
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store%res_buffer),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store%rhoin),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store%rhoin_old),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mixing_store%drho_buffer),cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(mixing_store%res_buffer(1,1)%cc)

    alpha=mixing_store%alpha
    beta=mixing_store%beta
    w0 = mixing_store%broy_w0
    nbuffer = mixing_store%nbuffer
    gapw = qs_env%dft_control%qs_control%gapw

    ib = MODULO(mixing_store%ncall,nbuffer) + 1
    mixing_store%ncall = mixing_store%ncall + 1
!dbg
!    ib = MIN(mixing_store%ncall,nbuffer)
!dbg
    nb = MIN(mixing_store%ncall,nbuffer)
    ibb = MODULO(mixing_store%ncall,nbuffer) + 1
!dbg
!    ibb = MIN(mixing_store%ncall+1,nbuffer)
!dbg

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

    IF(gapw) THEN
      CALL get_qs_env(qs_env=qs_env,&
            rho_atom_set=rho_atom,error=error)
      natom = SIZE(rho_atom)
    END IF

    DO ispin = 1,nspin
      res_rho  = CMPLX(0.0_dp,0.0_dp,KIND=dp)
      DO ig = 1,ng
         res_rho(ig) = rho%rho_g(ispin)%pw%cc(ig) - mixing_store%rhoin(ispin)%cc(ig)
      END DO

      IF(ib==1) THEN
        DO ig = 1,ng
          mixing_store%last_res(ispin)%cc(ig) = res_rho(ig)
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin(ispin)%cc(ig) + f_mix*res_rho(ig)
          mixing_store%rhoin_old(ispin)%cc(ig) =  mixing_store%rhoin(ispin)%cc(ig)
          mixing_store%rhoin(ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO

        IF(gapw) THEN
          DO iatom = 1,natom
            IF(mixing_store%paw(iatom)) THEN
              mixing_store%cpc_h_lastres(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef - &
                     mixing_store%cpc_h_in(iatom,ispin)%r_coef
              mixing_store%cpc_s_lastres(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef - &
                     mixing_store%cpc_s_in(iatom,ispin)%r_coef

              rho_atom(iatom)%cpc_h(ispin)%r_coef = alpha*rho_atom(iatom)%cpc_h(ispin)%r_coef + &
                     mixing_store%cpc_h_in(iatom,ispin)%r_coef*(1._dp-alpha)
              rho_atom(iatom)%cpc_s(ispin)%r_coef = alpha*rho_atom(iatom)%cpc_s(ispin)%r_coef + &
                     mixing_store%cpc_s_in(iatom,ispin)%r_coef*(1._dp-alpha)

              mixing_store%cpc_h_old(iatom,ispin)%r_coef = mixing_store%cpc_h_in(iatom,ispin)%r_coef
              mixing_store%cpc_s_old(iatom,ispin)%r_coef = mixing_store%cpc_s_in(iatom,ispin)%r_coef
              mixing_store%cpc_h_in(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
              mixing_store%cpc_s_in(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
            END IF
          END DO
        END IF
      ELSE

        delta_norm = 0.0_dp
        res_norm = 0.0_dp
        DO ig = 1, ng
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
        END DO
        DO ig = 1, ng
          mixing_store%drho_buffer(ib-1,ispin)%cc(ig) = &
                     mixing_store%rhoin(ispin)%cc(ig) - &
                     mixing_store%rhoin_old(ispin)%cc(ig)
        END DO
        CALL mp_sum(delta_norm, para_env%group)
        delta_norm = SQRT(delta_norm)
        CALL mp_sum(res_norm,para_env%group)
        res_norm = SQRT(res_norm)

        mixing_store%res_buffer(ib-1,ispin)%cc(:) = mixing_store%res_buffer(ib-1,ispin)%cc(:)/delta_norm
        mixing_store%drho_buffer(ib-1,ispin)%cc(:) = mixing_store%drho_buffer(ib-1,ispin)%cc(:)/delta_norm

        IF(gapw) THEN

          DO iatom = 1,natom
            IF(mixing_store%paw(iatom)) THEN
              n1 = SIZE(mixing_store%cpc_s_in(iatom,ispin)%r_coef,1)
              n2 = SIZE(mixing_store%cpc_s_in(iatom,ispin)%r_coef,2)
              DO i = 1,n2
              DO j = 1,n1
              mixing_store%dcpc_h_in(ib-1,iatom,ispin)%r_coef(j,i) = &
                    (mixing_store%cpc_h_in(iatom,ispin)%r_coef(j,i)-&
                     mixing_store%cpc_h_old(iatom,ispin)%r_coef(j,i))/delta_norm
              dcpc_h_res = ((rho_atom(iatom)%cpc_h(ispin)%r_coef(j,i) - &
                     mixing_store%cpc_h_in(iatom,ispin)%r_coef(j,i)) -  &
                     mixing_store%cpc_h_lastres(iatom,ispin)%r_coef(j,i))/delta_norm
              mixing_store%cpc_h_lastres(iatom,ispin)%r_coef(j,i) = rho_atom(iatom)%cpc_h(ispin)%r_coef(j,i) - &
                     mixing_store%cpc_h_in(iatom,ispin)%r_coef(j,i)

              mixing_store%dcpc_s_in(ib-1,iatom,ispin)%r_coef(j,i) = &
                    (mixing_store%cpc_s_in(iatom,ispin)%r_coef(j,i)-&
                     mixing_store%cpc_s_old(iatom,ispin)%r_coef(j,i))/delta_norm
               dcpc_s_res = ((rho_atom(iatom)%cpc_s(ispin)%r_coef(j,i) - &
                     mixing_store%cpc_s_in(iatom,ispin)%r_coef(j,i)) -  &
                     mixing_store%cpc_s_lastres(iatom,ispin)%r_coef(j,i))/delta_norm
              mixing_store%cpc_s_lastres(iatom,ispin)%r_coef(j,i) = rho_atom(iatom)%cpc_s(ispin)%r_coef(j,i) - &
                     mixing_store%cpc_s_in(iatom,ispin)%r_coef(j,i)

                 mixing_store%dcpc_h_in(ib-1,iatom,ispin)%r_coef(j,i) = &
                    alpha*dcpc_h_res +&
                     mixing_store%dcpc_h_in(ib-1,iatom,ispin)%r_coef(j,i)
                 mixing_store%dcpc_s_in(ib-1,iatom,ispin)%r_coef(j,i) = &
                    alpha*dcpc_s_res +&
                     mixing_store%dcpc_s_in(ib-1,iatom,ispin)%r_coef(j,i)
              END DO
              END DO
            END IF
          END DO

        END IF

        a(:,:) = 0.0_dp
        DO ig = 1, ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
          mixing_store%drho_buffer(ib-1,ispin)%cc(ig) = &
                 f_mix*mixing_store%res_buffer(ib-1,ispin)%cc(ig)+&
                 mixing_store%drho_buffer(ib-1,ispin)%cc(ig)
        END DO
        DO jb = 1,ib-1
           DO kb = jb,ib-1
              DO ig = 1, mixing_store%ig_max !ng
                  a(kb,jb) =  a(kb,jb) + mixing_store%p_metric(ig)*(&
                    REAL(mixing_store%res_buffer(jb,ispin)%cc(ig),dp)*&
                    REAL(mixing_store%res_buffer(kb,ispin)%cc(ig),dp) +&
                    AIMAG(mixing_store%res_buffer(jb,ispin)%cc(ig))*&
                    AIMAG(mixing_store%res_buffer(kb,ispin)%cc(ig)))
              END DO
              a(jb,kb) = a(kb,jb)
           END DO
        END DO
        CALL mp_sum(a,para_env%group)

        C = 0.0_dp
        DO jb = 1,ib-1
           a(jb,jb) = w0 + a(jb,jb)
           DO ig = 1, mixing_store%ig_max !ng
             c(jb) = c(jb) +  mixing_store%p_metric(ig)*(&
                     REAL(mixing_store%res_buffer(jb,ispin)%cc(ig),dp)*REAL(res_rho(ig),dp) + &
                     AIMAG(mixing_store%res_buffer(jb,ispin)%cc(ig))*AIMAG(res_rho(ig)))
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
          rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin(ispin)%cc(ig) +&
                   f_mix * res_rho(ig) + cc_mix
          mixing_store%rhoin_old(ispin)%cc(ig) =  mixing_store%rhoin(ispin)%cc(ig)
          mixing_store%rhoin(ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
        END DO

        IF(gapw) THEN
          DO iatom = 1,natom
            IF(mixing_store%paw(iatom)) THEN
              n1 = SIZE(mixing_store%cpc_s_in(iatom,ispin)%r_coef,1)
              n2 = SIZE(mixing_store%cpc_s_in(iatom,ispin)%r_coef,2)
              DO i = 1,n2
              DO j = 1,n1
                valh = 0.0_dp
                vals = 0.0_dp
                DO jb = 1,ib-1
                  valh = valh - G(jb)*mixing_store%dcpc_h_in(jb,iatom,ispin)%r_coef(j,i)
                  vals = vals - G(jb)*mixing_store%dcpc_s_in(jb,iatom,ispin)%r_coef(j,i)
                END DO
                rho_atom(iatom)%cpc_h(ispin)%r_coef(j,i) = alpha*rho_atom(iatom)%cpc_h(ispin)%r_coef(j,i) + &
                     mixing_store%cpc_h_in(iatom,ispin)%r_coef(j,i)*(1._dp-alpha) + valh
                rho_atom(iatom)%cpc_s(ispin)%r_coef(j,i) = alpha*rho_atom(iatom)%cpc_s(ispin)%r_coef(j,i) + &
                     mixing_store%cpc_s_in(iatom,ispin)%r_coef(j,i)*(1._dp-alpha) + vals
              END DO
              END DO

              mixing_store%cpc_h_old(iatom,ispin)%r_coef = mixing_store%cpc_h_in(iatom,ispin)%r_coef
              mixing_store%cpc_s_old(iatom,ispin)%r_coef = mixing_store%cpc_s_in(iatom,ispin)%r_coef
              mixing_store%cpc_h_in(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_h(ispin)%r_coef
              mixing_store%cpc_s_in(iatom,ispin)%r_coef = rho_atom(iatom)%cpc_s(ispin)%r_coef
            END IF
          END DO
        END IF

!  MI: this type of propagation makes the convergence worse, but it is not crystal-clear why
!      possibly if is still affected by some bug. Re-check!
!        IF(mixing_store%ncall>=nbuffer) THEN
!         DO jb = 1,nbuffer - 1
!            mixing_store%drho_buffer(jb,ispin)%cc(:) = mixing_store%drho_buffer(jb+1,ispin)%cc(:)
!            mixing_store%res_buffer(jb,ispin)%cc(:) = mixing_store%res_buffer(jb+1,ispin)%cc(:)
!         END DO
!       END IF

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
      imp, imp_j, inv_err, norm, norm_ig, rep, rep_j, sqt_uvol, sqt_vol, &
      uvol, vol, wc, wmax
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: aval, work_dgesdd
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: a, a2, au, av, b, bq, fm2, &
                                                work
    REAL(dp), DIMENSION(:), POINTER          :: p_metric
    REAL(dp), DIMENSION(:, :), POINTER       :: fmat, weight
    TYPE(cp_1d_z_p_type), DIMENSION(:), &
      POINTER                                :: tmp_z
    TYPE(cp_1d_z_p_type), DIMENSION(:, :), &
      POINTER                                :: delta_res, u_vec, z_vec

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store%rhoin_buffer),cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    NULLIFY(delta_res, u_vec, z_vec)
    NULLIFY(fmat)
    nspin = SIZE(rho%rho_g,1)
    ng = SIZE(mixing_store%rhoin_buffer(1,1)%cc)
    vol = rho%rho_g(1)%pw%pw_grid%vol
    sqt_vol =  SQRT(vol)
    uvol = 1.0_dp/vol
    sqt_uvol =  SQRT(uvol)
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

      delta = 0.0_dp
      delta_p = 0.0_dp
      ! Residual at this step R_i(G) (rho_out(G)-rho_in(G))
      ! Residual multiplied by the metrics RP_i(G) = (rho_out(G)-rho_in(G)) * P(G)
      ! Delta is the norm of the residual, measures how far we are from convergence
      DO ig = 1,ng
         res_rho(ig) = (rho%rho_g(ispin)%pw%cc(ig) - mixing_store%rhoin_buffer(ib,ispin)%cc(ig))*sqt_uvol
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

      weight(ib,ispin) = 1.0_dp ! wc
      IF(wc< 0.0_dp) weight(ib,ispin) = 0.01_dp*ABS(wc)/(delta_p*delta_p)
      IF(weight(ib,ispin)==0.0_dp) weight(ib,ispin) = 100.0_dp
      IF(weight(ib,ispin)<1.0_dp) weight(ib,ispin) = 1.0_dp

      IF(ib==1) THEN
        ! Simple Kerker damping : linear mixing rho(G) = rho_in(G) - alpha k(G)*(rho_out(G)-rho_in(G))
        DO ig = 1,ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
!          mixing_store%rhoin_buffer(ib,ispin)%cc(ig) = mixing_store%rhoin_buffer(ib,ispin)%cc(ig)*sqt_vol
          rho%rho_g(ispin)%pw%cc(ig) = (mixing_store%rhoin_buffer(ib,ispin)%cc(ig) + f_mix*res_rho(ig))
          mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
!          rho%rho_g(ispin)%pw%cc(ig) =  rho%rho_g(ispin)%pw%cc(ig)*sqt_vol
        END DO
      ELSE
        norm = 0.0_dp
        ! Difference of residuals DR_{i-1)} (G) = R_i(G) - R_{i-1}(G)
        DO ig = 1,ng
          delta_res(ib-1, ispin)%cc(ig) = res_rho(ig) - mixing_store%last_res(ispin)%cc(ig)
          delta_res_p(ig) = p_metric(ig) * (res_rho(ig) - mixing_store%last_res(ispin)%cc(ig))
          norm_ig = REAL(delta_res(ib-1, ispin)%cc(ig),dp)*REAL(delta_res_p(ig),dp)+&
                    AIMAG(delta_res(ib-1, ispin)%cc(ig))*AIMAG(delta_res_p(ig))
          norm = norm + norm_ig
        END DO
        CALL mp_sum(norm,para_env%group)
        norm = 1._dp/SQRT(norm)
        delta_res(ib-1, ispin)%cc(:) = delta_res(ib-1, ispin)%cc(:)*norm
        delta_res_p(:) = delta_res_p(:)*norm

!dbg
!    write(*,*) 'ibroy ', ib, weight(ib,ispin)
    IF(para_env%ionode) WRITE(*,*) 'norm ', norm
!dbg
        ! Vector U_{i-1}(G) =  Drho_{i-1} + k(G)  * DR_{i-1}(G)
        DO ig = 1,ng
          tmp(ig) = (mixing_store%rhoin_buffer(ib,ispin)%cc(ig) - &
                     mixing_store%rhoin_buffer(ib-1,ispin)%cc(ig))*norm
          u_vec(ib-1,ispin)%cc(ig) = (tmp(ig) + &
                     mixing_store%kerker_factor(ig) * delta_res(ib-1,ispin)%cc(ig))
        END DO

        DO jb = 1,ib-1
          fmat(jb,ib-1) = 0.0_dp

          DO ig = 1,ng
            rep_j = REAL(delta_res(jb, ispin)%cc(ig),dp)
            imp_j = AIMAG(delta_res(jb, ispin)%cc(ig))
           ! < DR_{j} | DR_{i-1} >
            rep = REAL(delta_res_p(ig),dp)
            imp = AIMAG(delta_res_p(ig))
            fmat(jb,ib-1) = fmat(jb,ib-1) + rep_j*rep + imp_j*imp
         END DO

        END DO
        CALL mp_sum(fmat(1:ib-1,ib-1),para_env%group)

        fmat(ib-1,ib-1) = 1.0_dp

        DO jb = 1,ib-2
          fmat(ib-1,jb) = fmat(jb,ib-1)
        ENDDO
!dbg
   IF(para_env%ionode) THEN
      WRITE(*,*) 'fmat '
        DO jb = 1,ib-1
          WRITE(*,*)  jb, (fmat(jb,kb),kb=1,ib-1)
        END DO
   END IF
!dbg

!dbg
        ALLOCATE (fm2(ib-1,ib-1),STAT=istat)
        fm2 = 0.0_dp
        DO kb=1,ib-1
          DO jb=1,kb
            DO ig = 1,ng
               rep_j = REAL(delta_res(jb, ispin)%cc(ig),dp)
               imp_j = AIMAG(delta_res(jb, ispin)%cc(ig))
               rep =   REAL(delta_res(kb, ispin)%cc(ig),dp)
               imp =   AIMAG(delta_res(kb, ispin)%cc(ig))
               fm2(jb,kb) = fm2(jb,kb) +( rep_j*rep + imp_j*imp)* p_metric(ig)
            END DO ! ig
          END DO
        END DO
        CALL mp_sum(fm2,para_env%group)
        DO jb = 1,ib-1
          fm2(jb,jb) = 1.0_dp
          DO kb = 1,jb-1
            fmat(jb,kb) = fmat(kb,jb)
          ENDDO
        END DO
   IF(para_env%ionode) THEN
      WRITE(*,*) 'fm2 '
        DO jb = 1,ib-1
          WRITE(*,*)  jb, (fm2(jb,kb),kb=1,ib-1)
        END DO
   END IF
!dbg

        DO jb = 1,ib-1
!          a(jb,jb) = 1.0_dp + weight(jb,ispin)*weight(jb,ispin)*fmat(jb,jb)
          a(jb,jb) =  weight(jb,ispin)*weight(jb,ispin)*fmat(jb,jb)
          DO kb = 1,jb-1
             a(jb,kb) = weight(jb,ispin)*weight(kb,ispin)*fmat(jb,kb)
             a(kb,jb) = weight(jb,ispin)*weight(kb,ispin)*fmat(kb,jb)
          ENDDO
        END DO
!dbg
   IF(para_env%ionode) THEN
      WRITE(*,*) 'a mat'
      DO kb = 1,ib-1
         WRITE(*,*)  'a ', kb, (a(kb,jb),jb=1,ib-1)
      END DO
   END IF
!dbg
!dbg
        ALLOCATE(a2(ib-1,ib-1),STAT=istat)
        DO jb = 1,ib-1
        DO kb = 1,ib-1
          a2(kb,jb) =  weight(jb,ispin)*weight(kb,ispin)*fm2(kb,jb)
        END DO
!          a2(jb,jb) = 1.0_dp + a2(jb,jb)
        END DO

        IF(.TRUE.) THEN
          b=0.0_dp
          CALL invert_matrix(a,b,inv_err,error=error)
!dbg
   IF(para_env%ionode) THEN
      WRITE(*,*) 'invert version 1'
      DO kb = 1,ib-1
         WRITE(*,*)  'b ', kb, (b(kb,jb),jb=1,ib-1)
      END DO
   END IF
!dbg
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

!          CALL diamat_all(work(1:ib-1,1:ib-1),aval(1:ib-1),error=error)
!
!          avec(1:ib-1,1:ib-1) = work(1:ib-1,1:ib-1)
!          DO jb = 1,ib-1
!            IF(ABS(aval(jb)) < 1.E-5_dp ) THEN
!              avec(1:ib-1,jb) = 0.0_dp
!            ELSE
!              avec(1:ib-1,jb) = work(1:ib-1,jb) / aval(jb)
!            END IF
!          END DO
!          b(1:ib-1,1:ib-1) = MATMUL(avec(1:ib-1,1:ib-1),work(1:ib-1,1:ib-1))

!          DEALLOCATE(work,avec,aval)
        END IF


        bq = 0.0_dp
        ! Broyden second method requires also bq (NYI)
        skip_bq = .TRUE.
        DO jb = 1,ib-1
!          skip_bq = skip_bq .AND. (weight(jb,ispin) > wmax)
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
          tmp(:) = CMPLX(0.0_dp,0.0_dp,KIND=dp)
          IF(.NOT. skip_bq) THEN
            ! sum_{kb} bq(jb,kb) * z_vec_{kb,iter-2}
            ! added only for small weights
            DO kb = 1,ib-2
               IF (weight(kb,ispin) >= (10.0_dp*wmax)) CYCLE
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
          IF(skip_bq .OR. (weight(jb,ispin) >=(10._dp*wmax))) THEN
             z_vec(jb,ispin)%cc(:) = tmp(:)
          ELSE
             ! temporary array: old z_vec may still be needed
             tmp_z(jb)%cc(:) = tmp(:)
          END IF
        END DO !jb

        IF(.NOT. skip_bq ) THEN
          DO jb = 1,ib-1
             IF (weight(jb,ispin) < (10._dp*wmax))  z_vec(jb,ispin)%cc(:) = tmp_z(jb)%cc(:)
             DEALLOCATE(tmp_z(jb)%cc,STAT=istat)
          END DO
          DEALLOCATE(tmp_z,STAT=istat)
        END IF

        ! Overwrite the density i reciprocal space
        rho%rho_g(ispin)%pw%cc(:) = CMPLX(0.0_dp,0.0_dp,KIND=dp)
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
!dbg
         IF(para_env%ionode) WRITE(*,*) 'norm ', norm
!dbg
          ! Subtract |Z_jb)><DR_jb|P|R_{iter}>
          DO ig = 1,ng
            rho%rho_g(ispin)%pw%cc(ig) = rho%rho_g(ispin)%pw%cc(ig) - norm * z_vec(jb,ispin)%cc(ig)* sqt_vol

          END DO
        END DO

        DO ig = 1,ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) =rho%rho_g(ispin)%pw%cc(ig) +&
                   mixing_store%rhoin_buffer(ib,ispin)%cc(ig) + f_mix*res_rho(ig)
          mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
!          rho%rho_g(ispin)%pw%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)*sqt_vol
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
!>        as preconditioner we use the Kerker damping factor
!>        The mixing is applied directly on the density in reciprocal space,
!>        therefore it affects the potentials
!>        on the grid but not the density matrix
!> \par History
!>      05.2009
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
    REAL(dp) :: alpha, f_mix, gn_norm, gn_norm_old, inv_err, n_low, n_up, &
      pgn_norm, prec, r_step, reg_par, sigma_max, sigma_tilde, step_size
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: norm_res, norm_res_low, &
                                                norm_res_up
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: b_matrix, binv_matrix
    REAL(dp), DIMENSION(:), POINTER          :: g2
    REAL(dp), SAVE                           :: sigma_old = 1.0_dp

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mixing_store),cp_failure_level,routineP,error,failure)

    CALL timeset(routineN,handle)

    NULLIFY(gn)

    use_zgemm = .FALSE.
    use_zgemm_rev = .TRUE.
!   use_zgemm_rev = .FALSE.

    ! prepare the parameters
    nspin = SIZE(rho%rho_g,1)
    ng_global = rho%rho_g(1)%pw%pw_grid%ngpts
    ng = SIZE(mixing_store%rhoin_buffer(1,1)%cc)
    alpha=mixing_store%alpha

    sigma_max = mixing_store%sigma_max
    reg_par = mixing_store%reg_par
    r_step = mixing_store%r_step
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

!dbg
    IF(para_env%ionode) WRITE(*,*) mixing_store%ncall, ib, nb, ib_prev, ib_next
!dbg

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
      !simple mixing
      DO ispin = 1,nspin
        DO ig = 1,ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(1,ispin)%cc(ig) + &
              f_mix*mixing_store%res_buffer(1,ispin)%cc(ig)
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
    ALLOCATE(norm_res_low(nb),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(norm_res_up(nb),STAT=istat)
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
      gn_global = CMPLX(0.0_dp, 0.0_dp,KIND=dp)
      DO ig = 1,ng
         ig_global = mixing_store%ig_global_index(ig)
         gn_global(ig_global) = gn(ig)
      END DO
      CALL mp_sum(gn_global,para_env%group)

    ! compute steps (matrix S) and residual differences (matrix Y)
    ! with respect to the present
    ! step and the present residual (use stored rho_in and res_buffer)

    ! These quantities are pre-conditioned by means of Kerker factor multipliccation
      g2 =>rho%rho_g(1)%pw%pw_grid%gsq

      DO jb = 1,nb+1
        IF(jb<ib) THEN
           kb = jb
        ELSEIF(jb>ib) THEN
           kb = jb-1
        ELSE
           CYCLE
        END IF
        norm_res(kb) = 0.0_dp
        norm_res_low(kb) = 0.0_dp
        norm_res_up(kb) = 0.0_dp
        DO ig = 1,ng

!          prec = mixing_store%kerker_factor(ig)
          prec = 1.0_dp

          step_matrix(ig,kb) = prec*(mixing_store%rhoin_buffer(jb,ispin)%cc(ig) - &
                               mixing_store%rhoin_buffer(ib,ispin)%cc(ig))
          res_matrix(ig,kb) = (mixing_store%res_buffer(jb,ispin)%cc(ig) -&
                             mixing_store%res_buffer(ib,ispin)%cc(ig))
          norm_res(kb)= norm_res(kb) + REAL(res_matrix(ig,kb),dp)*REAL(res_matrix(ig,kb),dp) +&
                        AIMAG(res_matrix(ig,kb))* AIMAG(res_matrix(ig,kb))
          IF(g2(ig)<4.0_dp) THEN
             norm_res_low(kb) = norm_res_low(kb) +&
                        REAL(res_matrix(ig,kb),dp)*REAL(res_matrix(ig,kb),dp) +&
                        AIMAG(res_matrix(ig,kb))* AIMAG(res_matrix(ig,kb))
          ELSE
             norm_res_up(kb) = norm_res_up(kb) +&
                        REAL(res_matrix(ig,kb),dp)*REAL(res_matrix(ig,kb),dp) +&
                        AIMAG(res_matrix(ig,kb))* AIMAG(res_matrix(ig,kb))
          END IF
          res_matrix(ig,kb) = prec * res_matrix(ig,kb)
        END DO
      END DO  !jb

    ! normalize each column of S and Y => Snorm Ynorm
      CALL mp_sum(norm_res,para_env%group)
      CALL mp_sum(norm_res_up,para_env%group)
      CALL mp_sum(norm_res_low,para_env%group)
!dbg
      IF(para_env%ionode)  WRITE(*,*) 'norm res   ' , norm_res(1:nb)
      IF(para_env%ionode)  WRITE(*,*) 'norm res  L' , norm_res_low(1:nb)
      IF(para_env%ionode)  WRITE(*,*) 'norm res  U' , norm_res_up(1:nb)
!dbg
      norm_res(1:nb) = 1.0_dp/SQRT(norm_res(1:nb))
      n_low = 0.0_dp
      n_up = 0.0_dp
      DO jb=1,nb
        n_low = n_low + norm_res_low(jb)/norm_res(jb)
        n_up = n_up + norm_res_up(jb)/norm_res(jb)
      END DO
      DO ig=1,ng
         IF(g2(ig)>4.0_dp) THEN
            step_matrix(ig,1:nb) = step_matrix(ig,1:nb) * SQRT(n_low/n_up)
            res_matrix(ig,1:nb) = res_matrix(ig,1:nb)  * SQRT(n_low/n_up)
         END IF
      END DO
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
      a_matrix = CMPLX(0.0_dp,0.0_dp,KIND=dp)
      DO kb = 1,nb
         res_matrix_global = CMPLX(0.0_dp, 0.0_dp,KIND=dp)
         DO ig = 1,ng
            ig_global = mixing_store%ig_global_index(ig)
            res_matrix_global(ig_global) = res_matrix(ig,kb)
         END DO
         CALL mp_sum(res_matrix_global,para_env%group)

         DO jb = 1,nb
            DO ig = 1,ng
               ig_global = mixing_store%ig_global_index(ig)
               a_matrix(jb,ig_global) = a_matrix(jb,ig_global) + &
                  binv_matrix(jb,kb)*res_matrix_global(ig_global)
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
           ya(ig,ig_global) =  ya(ig,ig_global) + CMPLX(1.0_dp,0.0_dp,KIND=dp)
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
          pgn(ig) =   CMPLX(0.0_dp,0.0_dp,KIND=dp)
          ugn(ig) =   CMPLX(0.0_dp,0.0_dp,KIND=dp)
          ig_global =  mixing_store%ig_global_index(ig)
          DO iig = 1,ng_global
             saa = CMPLX(0.0_dp,0.0_dp,KIND=dp)
             yaa = CMPLX(0.0_dp,0.0_dp,KIND=dp)

!            IF(ig_global==iig) ya(ig,ig_global) =  CMPLX(1.0_dp,0.0_dp,KIND=dp)
            IF(ig_global==iig) yaa =  CMPLX(1.0_dp,0.0_dp,KIND=dp)

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
        sigma_tilde = 0.1_dp
      ! Step size for the unpredicted component
      step_size = MIN(sigma_tilde,r_step*pgn_norm/gn_norm,sigma_max)
      sigma_old =  step_size

    ! update the density
!dbg
      IF(para_env%ionode)  WRITE(*,*) 'step ', step_size, pgn_norm,gn_norm,sigma_tilde
      DO ig = 1,ng
         prec = mixing_store%kerker_factor(ig)
         rho%rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(ib,ispin)%cc(ig)  &
              -prec*step_size*ugn(ig) +  prec*pgn(ig)! - 0.1_dp * prec* gn(ig)
         mixing_store%rhoin_buffer(ib_next,ispin)%cc(ig) = rho%rho_g(ispin)%pw%cc(ig)
      END DO

    END DO ! ispin

    ! Deallocate  temporary arrays
    DEALLOCATE(step_matrix,res_matrix, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(norm_res, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(norm_res_up, STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(norm_res_low, STAT=istat)
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

! *****************************************************************************
!> \brief  allocation needed when density mixing is used
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2009 created [MI]
!> \author fawzi
! *****************************************************************************
 SUBROUTINE mixing_allocate(qs_env,mixing_method,mixing_store,&
      nspins,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: mixing_method
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mixing_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iat, ispin, natom, &
                                                nbuffer, stat
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    failure = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY(sab_orb, scf_env, matrix_s)
    CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb, &
         scf_env=scf_env, &
         matrix_s=matrix_s, &
         error=error)

!   *** allocate p_mix_new ***
    IF (.NOT.ASSOCIATED(scf_env%p_mix_new)) THEN
       CALL cp_dbcsr_allocate_matrix_set(scf_env%p_mix_new,nspins,error=error)
       DO i=1,nspins
          ALLOCATE(scf_env%p_mix_new(i)%matrix)
          CALL cp_dbcsr_init(scf_env%p_mix_new(i)%matrix, error=error)
          CALL cp_dbcsr_create(matrix=scf_env%p_mix_new(i)%matrix, &
               name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//"DENSITY", &
               dist=cp_dbcsr_distribution(matrix_s(1)%matrix), matrix_type=dbcsr_type_symmetric,&
               row_blk_size=cp_dbcsr_row_block_sizes(matrix_s(1)%matrix), &
               col_blk_size=cp_dbcsr_col_block_sizes(matrix_s(1)%matrix), &
               nblks=0, nze=0, error=error)
          CALL cp_dbcsr_alloc_block_from_nbl(scf_env%p_mix_new(i)%matrix,sab_orb,error=error)
          CALL cp_dbcsr_set(scf_env%p_mix_new(i)%matrix,0.0_dp,error=error)
       ENDDO
    END IF

!   *** allocate p_delta ***
    IF (mixing_method>=gspace_mixing_nr) THEN
       IF(.NOT.ASSOCIATED(scf_env%p_delta)) THEN
          CALL cp_dbcsr_allocate_matrix_set(scf_env%p_delta,nspins,error=error)
          DO i=1,nspins
             ALLOCATE(scf_env%p_delta(i)%matrix)
             CALL cp_dbcsr_init(scf_env%p_delta(i)%matrix, error=error)
             CALL cp_dbcsr_create(matrix=scf_env%p_delta(i)%matrix, &
                  name="SCF"//TRIM(ADJUSTL(cp_to_string(scf_env%id_nr)))//"DENSITY", &
                  dist=cp_dbcsr_distribution(matrix_s(1)%matrix), matrix_type=dbcsr_type_symmetric,&
                  row_blk_size=cp_dbcsr_row_block_sizes(matrix_s(1)%matrix), &
                  col_blk_size=cp_dbcsr_col_block_sizes(matrix_s(1)%matrix), &
                  nblks=0, nze=0, error=error)
             CALL cp_dbcsr_alloc_block_from_nbl(scf_env%p_delta(i)%matrix,sab_orb,error=error)
             CALL cp_dbcsr_set(scf_env%p_delta(i)%matrix,0.0_dp,error=error)
          ENDDO
       END IF
       CPPrecondition(ASSOCIATED(mixing_store),cp_failure_level,routineP,error,failure)
    END IF

!   *** allocate  buffer for gspace mixing ***
    IF (mixing_method>=gspace_mixing_nr) THEN
      nbuffer = mixing_store%nbuffer
      mixing_store%ncall = 0
      IF(.NOT. ASSOCIATED(mixing_store%rhoin)) THEN
        ALLOCATE(mixing_store%rhoin(nspins), STAT=stat)
        DO ispin = 1,nspins
          NULLIFY(mixing_store%rhoin(ispin)%cc)
        END DO
      END IF
      IF(qs_env % dft_control % qs_control%gapw) THEN
        CALL get_qs_env(qs_env=qs_env,&
               rho_atom_set=rho_atom,error=error)
        natom = SIZE(rho_atom)
        IF(.NOT. ASSOCIATED(mixing_store%paw)) THEN
          ALLOCATE(mixing_store%paw(natom),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          mixing_store%paw = .FALSE.
          ALLOCATE(mixing_store%cpc_h_in(natom,nspins),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mixing_store%cpc_s_in(natom,nspins),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ispin = 1,nspins
            DO iat = 1,natom
               NULLIFY(mixing_store%cpc_h_in(iat,ispin)%r_coef)
               NULLIFY(mixing_store%cpc_s_in(iat,ispin)%r_coef)
            END DO
          END DO
        END IF
      END IF
    END IF

!   *** allocate rhoin_buffer if needed
    IF (mixing_method==pulay_mixing_nr .OR.mixing_method==broyden_mixing_new_nr&
           .OR. mixing_method==multisecant_mixing_nr ) THEN
      IF(.NOT. ASSOCIATED(mixing_store%rhoin_buffer)) THEN
        ALLOCATE(mixing_store%rhoin_buffer(nbuffer,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin = 1,nspins
          DO i = 1,nbuffer
            NULLIFY(mixing_store%rhoin_buffer(i,ispin)%cc)
          END DO
        END DO
      END IF
    END IF

!   *** allocare res_buffer if needed
    IF (mixing_method>=pulay_mixing_nr) THEN
      IF(.NOT. ASSOCIATED(mixing_store%res_buffer)) THEN
        ALLOCATE(mixing_store%res_buffer(nbuffer,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin = 1,nspins
          DO i = 1,nbuffer
            NULLIFY(mixing_store%res_buffer(i,ispin)%cc)
          END DO
        END DO
      END IF
    END IF

!   *** allocate pulay buffer ***
    IF (mixing_method==pulay_mixing_nr) THEN
      IF(.NOT. ASSOCIATED(mixing_store%pulay_matrix)) THEN
        ALLOCATE(mixing_store%pulay_matrix(nbuffer,nbuffer),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF

    END IF
!   *** allocate broyden buffer ***
    IF (mixing_method==broyden_mixing_nr) THEN
      IF(.NOT. ASSOCIATED(mixing_store%rhoin_old)) THEN
        ALLOCATE(mixing_store%rhoin_old(nspins), STAT=stat)
        DO ispin = 1,nspins
          NULLIFY(mixing_store%rhoin_old(ispin)%cc)
        END DO
      END IF
      IF(.NOT. ASSOCIATED(mixing_store%drho_buffer)) THEN
        ALLOCATE(mixing_store%drho_buffer(nbuffer,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(mixing_store%last_res(nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin = 1,nspins
          DO i = 1,nbuffer
           NULLIFY(mixing_store%drho_buffer(i,ispin)%cc)
          END DO
          NULLIFY(mixing_store%last_res(ispin)%cc)
        END DO
      END IF
      IF(qs_env % dft_control % qs_control%gapw) THEN
        IF(.NOT. ASSOCIATED(mixing_store%cpc_h_old)) THEN
          ALLOCATE(mixing_store%cpc_h_old(natom,nspins),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mixing_store%cpc_s_old(natom,nspins),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ispin = 1,nspins
            DO iat = 1,natom
               NULLIFY(mixing_store%cpc_h_old(iat,ispin)%r_coef)
               NULLIFY(mixing_store%cpc_s_old(iat,ispin)%r_coef)
            END DO
          END DO
        END IF
        IF(.NOT. ASSOCIATED(mixing_store%dcpc_h_in)) THEN
          ALLOCATE(mixing_store%dcpc_h_in(nbuffer,natom,nspins),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mixing_store%dcpc_s_in(nbuffer,natom,nspins),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mixing_store%cpc_h_lastres(natom,nspins),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(mixing_store%cpc_s_lastres(natom,nspins),STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ispin = 1,nspins
            DO iat = 1,natom
              DO i = 1,nbuffer
                NULLIFY(mixing_store%dcpc_h_in(i,iat,ispin)%r_coef)
                NULLIFY(mixing_store%dcpc_s_in(i,iat,ispin)%r_coef)
              END DO
              NULLIFY(mixing_store%cpc_h_lastres(iat,ispin)%r_coef)
              NULLIFY(mixing_store%cpc_s_lastres(iat,ispin)%r_coef)
            END DO
          END DO
        END IF
      END IF
    END IF
!   *** allocate broyden buffer ***
    IF (mixing_method==broyden_mixing_new_nr) THEN
      IF(.NOT. ASSOCIATED(mixing_store%u_vec)) THEN
        ALLOCATE(mixing_store%last_res(nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(mixing_store%delta_res(nbuffer-1,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(mixing_store%u_vec(nbuffer-1,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(mixing_store%z_vec(nbuffer-1,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(mixing_store%weight(nbuffer,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ALLOCATE(mixing_store%fmat(nbuffer-1,nbuffer-1,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        DO ispin = 1,nspins
        DO i = 1,nbuffer-1
           NULLIFY(mixing_store%delta_res(i,ispin)%cc)
           NULLIFY(mixing_store%u_vec(i,ispin)%cc)
           NULLIFY(mixing_store%z_vec(i,ispin)%cc)
        END DO
        NULLIFY(mixing_store%last_res(ispin)%cc)
        END DO
      END IF
    END IF

!   *** allocate multisecant buffer ***
    IF (mixing_method==multisecant_mixing_nr) THEN
      IF(.NOT. ASSOCIATED(mixing_store%norm_res_buffer)) THEN
        ALLOCATE(mixing_store%norm_res_buffer(nbuffer,nspins), STAT=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF
    END IF

    CALL timestop(handle)

 END SUBROUTINE mixing_allocate

! *****************************************************************************
!> \brief  initialiation needed when gspace mixing is used
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2009 created [MI]
!> \author MI
! *****************************************************************************
 SUBROUTINE mixing_init(mixing_method,rho,mixing_store,para_env,rho_atom,error)
    INTEGER, INTENT(IN)                      :: mixing_method
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_atom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'mixing_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iat, ib, ig, ig1, &
                                                ig_count, iproc, ispin, &
                                                istat, n1, n2, natom, &
                                                nbuffer, ng, nspin
    LOGICAL                                  :: failure
    REAL(dp)                                 :: bconst, beta, fdamp, g2max, &
                                                g2min, kmin
    REAL(dp), DIMENSION(:), POINTER          :: g2
    REAL(dp), DIMENSION(:, :), POINTER       :: g_vec

    CALL timeset(routineN,handle)

    NULLIFY(g2,g_vec)
    failure = .FALSE.
    nspin = SIZE(rho%rho_g)
    ng = SIZE(rho%rho_g(1)%pw%pw_grid%gsq,1)
    mixing_store%ig_max = ng
    g2 =>rho%rho_g(1)%pw%pw_grid%gsq
    g_vec => rho%rho_g(1)%pw%pw_grid%g

    IF(mixing_store%max_gvec_exp > 0._dp) THEN
      DO ig =1,ng
        IF(g2(ig) > mixing_store%max_g2) THEN
           mixing_store%ig_max = ig
          EXIT
        END IF
      END DO
    END IF

    IF(.NOT. ASSOCIATED(mixing_store%kerker_factor)) THEN
      ALLOCATE(mixing_store%kerker_factor(ng),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF(.NOT. ASSOCIATED(mixing_store%special_metric)) THEN
      ALLOCATE(mixing_store%special_metric(ng),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF
    beta = mixing_store%beta
    kmin = 0.1_dp
    mixing_store%kerker_factor = 1.0_dp
    mixing_store%special_metric = 1.0_dp
    ig1 = 1
    IF (rho%rho_g(1)%pw%pw_grid%have_g0) ig1 = 2
    DO ig =ig1,mixing_store%ig_max
      mixing_store%kerker_factor(ig) = MAX(g2(ig)/(g2(ig)+beta*beta),kmin)
      mixing_store%special_metric(ig) = 1.0_dp+50.0_dp/8.0_dp*( 1.0_dp+&
         COS(g_vec(1,ig))+COS(g_vec(2,ig))+COS(g_vec(3,ig)) + COS(g_vec(1,ig))*COS(g_vec(2,ig)) +&
         COS(g_vec(2,ig))*COS(g_vec(3,ig)) + COS(g_vec(1,ig))*COS(g_vec(3,ig)) + &
         COS(g_vec(1,ig))*COS(g_vec(2,ig))*COS(g_vec(3,ig)))
    END DO

    nbuffer = mixing_store%nbuffer
    DO ispin = 1,nspin
      IF(.NOT. ASSOCIATED(mixing_store%rhoin(ispin)%cc)) THEN
        ALLOCATE(mixing_store%rhoin(ispin)%cc(ng),STAT=istat)
      END IF
      mixing_store%rhoin(ispin)%cc = rho%rho_g(ispin)%pw%cc

      IF(ASSOCIATED(mixing_store%rhoin_buffer)) THEN
        IF(.NOT. ASSOCIATED(mixing_store%rhoin_buffer(1,ispin)%cc)) THEN
           DO ib = 1,nbuffer
              ALLOCATE(mixing_store%rhoin_buffer(ib,ispin)%cc(ng),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
            END DO
         END IF
         mixing_store%rhoin_buffer(1,ispin)%cc(1:ng) = &
                       rho%rho_g(ispin)%pw%cc(1:ng)
      END IF
      IF(ASSOCIATED(mixing_store%res_buffer)) THEN
        IF(.NOT. ASSOCIATED(mixing_store%res_buffer(1,ispin)%cc)) THEN
          DO ib = 1,nbuffer
            ALLOCATE(mixing_store%res_buffer(ib,ispin)%cc(ng),STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          END DO
        END IF
      END IF
    END DO

    IF(nspin==2) THEN
      mixing_store%rhoin(1)%cc =  rho%rho_g(1)%pw%cc +  rho%rho_g(2)%pw%cc
      mixing_store%rhoin(2)%cc =  rho%rho_g(1)%pw%cc -  rho%rho_g(2)%pw%cc
      IF(ASSOCIATED(mixing_store%rhoin_buffer)) THEN
        mixing_store%rhoin_buffer(1,1)%cc =  rho%rho_g(1)%pw%cc +  rho%rho_g(2)%pw%cc
        mixing_store%rhoin_buffer(1,2)%cc =  rho%rho_g(1)%pw%cc -  rho%rho_g(2)%pw%cc
      END IF
    END IF


    IF(PRESENT(rho_atom)) THEN
     natom = SIZE(rho_atom)
     DO ispin = 1,nspin
     DO iat = 1,natom
       IF(ASSOCIATED(rho_atom(iat)%cpc_s(ispin)%r_coef)) THEN
         mixing_store%paw(iat) = .TRUE.
         n1 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef,1)
         n2 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef,2)
         IF(ASSOCIATED(mixing_store%cpc_s_in)) THEN
           IF(.NOT. ASSOCIATED(mixing_store%cpc_s_in(iat,ispin)%r_coef)) THEN
              ALLOCATE(mixing_store%cpc_s_in(iat,ispin)%r_coef(n1,n2),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
              ALLOCATE(mixing_store%cpc_h_in(iat,ispin)%r_coef(n1,n2),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
           END IF
           mixing_store%cpc_h_in(iat,ispin)%r_coef = rho_atom(iat)%cpc_h(ispin)%r_coef
           mixing_store%cpc_s_in(iat,ispin)%r_coef = rho_atom(iat)%cpc_s(ispin)%r_coef
         END IF
       END IF
     END DO
     END DO
    END IF

    IF (mixing_method==gspace_mixing_nr) THEN
    ELSEIF(mixing_method==pulay_mixing_nr) THEN
    ELSEIF(mixing_method==broyden_mixing_nr) THEN
       DO ispin = 1,nspin
         IF(.NOT. ASSOCIATED(mixing_store%rhoin_old(ispin)%cc)) THEN
            ALLOCATE(mixing_store%rhoin_old(ispin)%cc(ng),STAT=istat)
         END IF
         IF(.NOT. ASSOCIATED(mixing_store%drho_buffer(1,ispin)%cc)) THEN
           DO ib = 1,nbuffer
             ALLOCATE(mixing_store%drho_buffer(ib,ispin)%cc(ng),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
           END DO
           ALLOCATE(mixing_store%last_res(ispin)%cc(ng),STAT=istat)
           CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         END IF
         DO ib = 1,nbuffer
            mixing_store%drho_buffer(ib,ispin)%cc = CMPLX(0.0_dp,0.0_dp,kind=dp)
         END DO
         mixing_store%last_res(ispin)%cc = CMPLX(0.0_dp,0.0_dp,kind=dp)
         mixing_store%rhoin_old(ispin)%cc = CMPLX(0.0_dp,0.0_dp,kind=dp)
       END DO
       IF(PRESENT(rho_atom)) THEN
       DO ispin = 1,nspin
       DO iat = 1,natom
         IF(mixing_store%paw(iat)) THEN
           n1 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef,1)
           n2 = SIZE(rho_atom(iat)%cpc_s(ispin)%r_coef,2)
           IF(.NOT. ASSOCIATED(mixing_store%cpc_s_old(iat,ispin)%r_coef)) THEN
             ALLOCATE(mixing_store%cpc_s_old(iat,ispin)%r_coef(n1,n2),STAT=istat)
             ALLOCATE(mixing_store%cpc_h_old(iat,ispin)%r_coef(n1,n2),STAT=istat)
           END IF
           mixing_store%cpc_h_old(iat,ispin)%r_coef = 0.0_dp
           mixing_store%cpc_s_old(iat,ispin)%r_coef = 0.0_dp
           IF(.NOT. ASSOCIATED(mixing_store%dcpc_s_in(1,iat,ispin)%r_coef)) THEN
             DO ib = 1,nbuffer
               ALLOCATE(mixing_store%dcpc_h_in(ib,iat,ispin)%r_coef(n1,n2),STAT=istat)
               CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
               ALLOCATE(mixing_store%dcpc_s_in(ib,iat,ispin)%r_coef(n1,n2),STAT=istat)
               CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             END DO
             ALLOCATE(mixing_store%cpc_h_lastres(iat,ispin)%r_coef(n1,n2),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(mixing_store%cpc_s_lastres(iat,ispin)%r_coef(n1,n2),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
           END IF
           DO ib = 1,nbuffer
             mixing_store%dcpc_h_in(ib,iat,ispin)%r_coef = 0.0_dp
             mixing_store%dcpc_s_in(ib,iat,ispin)%r_coef = 0.0_dp
           END DO
           mixing_store%cpc_h_lastres(iat,ispin)%r_coef=0.0_dp
           mixing_store%cpc_s_lastres(iat,ispin)%r_coef=0.0_dp
         END IF
       END DO
       END DO
       END IF

       IF(.NOT. ASSOCIATED(mixing_store%p_metric)) THEN
         ALLOCATE(mixing_store%p_metric(ng),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
              bconst = mixing_store%bconst
              g2min = 1.E30_dp
              DO ig = 1,ng
                IF(g2(ig)>1.E-10_dp) g2min =  MIN(g2min,g2(ig))
              END DO
         g2max = -1.E30_dp
         DO ig = 1,ng
           g2max =  MAX(g2max,g2(ig))
         END DO
         CALL mp_min(g2min,para_env%group)
         CALL mp_max(g2max,para_env%group)
!      fdamp/g2 varies between (bconst-1) and 0
!      i.e. p_metric varies between bconst and 1
!         fdamp = (bconst-1.0_dp)*g2min
         fdamp = (bconst-1.0_dp)*g2min*g2max/(g2max-g2min*bconst)
         DO ig = 1,ng
             mixing_store%p_metric(ig) = (g2(ig)+fdamp)/MAX(g2(ig),1.E-10_dp)
         END DO
         IF(rho%rho_g(1)%pw%pw_grid%have_g0) mixing_store%p_metric(1) = bconst
       END IF

    ELSEIF(mixing_method==broyden_mixing_new_nr) THEN
       DO ispin = 1,nspin
         IF(.NOT. ASSOCIATED(mixing_store%u_vec(1,ispin)%cc)) THEN
            DO ib = 1,nbuffer-1
              ALLOCATE(mixing_store%delta_res(ib,ispin)%cc(ng),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
              ALLOCATE(mixing_store%u_vec(ib,ispin)%cc(ng),STAT=istat)
               CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
              ALLOCATE(mixing_store%z_vec(ib,ispin)%cc(ng),STAT=istat)
              CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
            END DO
            ALLOCATE(mixing_store%last_res(ispin)%cc(ng),STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         END IF
       END DO
       IF(.NOT. ASSOCIATED(mixing_store%p_metric)) THEN
         ALLOCATE(mixing_store%p_metric(ng),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         bconst = mixing_store%bconst
         g2min = 1.E30_dp
         DO ig = 1,ng
            IF(g2(ig)>1.E-10_dp) g2min =  MIN(g2min,g2(ig))
         END DO
         CALL mp_min(g2min,para_env%group)
         fdamp = (bconst-1.0_dp)*g2min
         DO ig = 1,ng
             mixing_store%p_metric(ig) = (g2(ig)+fdamp)/MAX(g2(ig),1.E-10_dp)
         END DO
         IF(rho%rho_g(1)%pw%pw_grid%have_g0) mixing_store%p_metric(1) = bconst
!dbg
!        mixing_store%p_metric = 1.0_dp
!dbg
      END IF
    ELSEIF(mixing_method==multisecant_mixing_nr) THEN
      IF(.NOT. ASSOCIATED(mixing_store%ig_global_index)) THEN
         ALLOCATE(mixing_store%ig_global_index(ng),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      mixing_store%ig_global_index = 0
      ig_count = 0
      DO iproc = 0,para_env%num_pe-1
        IF(para_env%mepos == iproc) THEN
          DO ig = 1,ng
             ig_count = ig_count+1
             mixing_store%ig_global_index(ig) = ig_count
          END DO
        END IF
        CALL mp_bcast(ig_count,iproc,para_env%group)
      END DO
    END IF

    CALL timestop(handle)

 END SUBROUTINE mixing_init


END MODULE qs_gspace_mixing
