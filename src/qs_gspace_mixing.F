!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE qs_gspace_mixing

  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_allocate_matrix_set, cp_dbcsr_col_block_sizes, &
       cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_row_block_sizes, cp_dbcsr_set, &
       cp_dbcsr_type, dbcsr_type_symmetric
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all,&
                                             invert_matrix
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_min,&
                                             mp_sum
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
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
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_density_mixing_types,         ONLY: broyden_mixing_new_nr,&
                                             broyden_mixing_nr,&
                                             cp_1d_z_p_type,&
                                             gspace_mixing_nr,&
                                             mixing_storage_type,&
                                             multisecant_mixing_nr,&
                                             pulay_mixing_nr
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_scf_methods,                  ONLY: cp_sm_mix
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gspace_mixing'

  PUBLIC :: gspace_mixing,&
            mixing_allocate, mixing_init,&
            self_consistency_check

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param rho_ao ...
!> \param p_delta ...
!> \param para_env ...
!> \param p_out ...
!> \param delta ...
! *****************************************************************************
  SUBROUTINE self_consistency_check(rho_ao,p_delta,para_env,p_out,delta)
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: rho_ao, p_delta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: p_out
    REAL(KIND=dp), INTENT(INOUT)             :: delta

    CHARACTER(len=*), PARAMETER :: routineN = 'self_consistency_check', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ic, ispin, nimg, &
                                                nspins
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_q, p_in

    CALL timeset(routineN,handle)

    NULLIFY(matrix_q, p_in)

    failure = .FALSE.
    CPASSERT(ASSOCIATED(p_out))
    NULLIFY(matrix_q, p_in)
    p_in => rho_ao
    matrix_q => p_delta
    nspins = SIZE(p_in,1)
    nimg = SIZE(p_in,2)

    ! Compute the difference (p_out - p_in)and check convergence
    delta = 0.0_dp
    DO ispin=1, nspins
       DO ic=1, nimg
          CALL cp_dbcsr_set(matrix_q(ispin,ic)%matrix, 0.0_dp)
          CALL cp_sm_mix(m1=p_out(ispin,ic)%matrix,m2=p_in(ispin,ic)%matrix,&
               p_mix=1.0_dp,delta=tmp, para_env=para_env,&
               m3=matrix_q(ispin,ic)%matrix)
         delta = MAX(tmp,delta)
       END DO
    END DO
    CALL timestop(handle)

  END SUBROUTINE self_consistency_check

! *****************************************************************************
!> \brief  Driver for the g-space mixing, calls the proper routine given the
!>requested method
!> \param qs_env ...
!> \param mixing_method ...
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
!> \param iter_count ...
!> \par History
!>      05.2009
!>      02.2015 changed input to make g-space mixing available in linear scaling 
!>              SCF [Patrick Seewald]
!> \author MI
! *****************************************************************************
  SUBROUTINE gspace_mixing(qs_env, mixing_method, mixing_store, rho, para_env, iter_count)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER                                  :: mixing_method
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: iter_count

    CHARACTER(len=*), PARAMETER :: routineN = 'gspace_mixing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, ig, ispin, &
                                                natom, ng, nspin
    LOGICAL                                  :: failure, gapw
    REAL(dp)                                 :: alpha
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_tmp
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPASSERT(ASSOCIATED(mixing_store))
    NULLIFY(auxbas_pw_pool, dft_control, pw_env, rho_atom, rho_g, rho_r, tot_rho_r)

    CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env)
    CALL qs_rho_get(rho, rho_g=rho_g, rho_r=rho_r, tot_rho_r= tot_rho_r)

    nspin = SIZE(rho_g,1)
    ng = SIZE(rho_g(1)%pw%pw_grid%gsq)
    CPASSERT((ng == SIZE(mixing_store%rhoin(1)%cc)))

    alpha = mixing_store%alpha
    gapw = dft_control%qs_control%gapw

    IF(nspin==2) THEN
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool)
      CALL pw_pool_create_pw(auxbas_pw_pool,&
                           rho_tmp%pw,&
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE)
      CALL pw_zero(rho_tmp%pw)
      CALL pw_copy(rho_g(1)%pw, rho_tmp%pw)
      CALL pw_axpy(rho_g(2)%pw,rho_g(1)%pw, 1.0_dp)
      CALL pw_axpy(rho_tmp%pw,rho_g(2)%pw, -1.0_dp)
      CALL pw_scale(rho_g(2)%pw, -1.0_dp)
    END IF

    IF(iter_count+1 <= mixing_store%nskip_mixing) THEN
      ! skip mixing
      DO ispin = 1,nspin
        DO ig = 1,ng
!          mixing_store%rhoin_buffer(1,ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
          mixing_store%rhoin(ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
        END DO
      END DO
      IF(gapw) THEN
        CALL get_qs_env(qs_env=qs_env,&
              rho_atom_set=rho_atom)
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
        CALL pw_axpy(rho_g(2)%pw,rho_g(1)%pw, 1.0_dp)
        CALL pw_scale(rho_g(1)%pw, 0.5_dp)
        CALL pw_axpy(rho_g(1)%pw,rho_g(2)%pw, -1.0_dp)
        CALL pw_scale(rho_g(2)%pw,-1.0_dp)
        CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tmp%pw)
      END IF
      CALL timestop(handle)
      RETURN
    END IF

    IF((iter_count+1 - mixing_store%nskip_mixing) <= mixing_store%n_simple_mix) THEN
       CALL gmix_potential_only(qs_env,mixing_store, rho)
       mixing_store%iter_method= "Kerker"
    ELSE

       IF(mixing_method==gspace_mixing_nr) THEN
          CALL gmix_potential_only(qs_env,mixing_store, rho)
          mixing_store%iter_method= "Kerker"
       ELSEIF(mixing_method==pulay_mixing_nr) THEN
          CPASSERT(.NOT. gapw)
          CALL pulay_mixing(mixing_store, rho, para_env)
          mixing_store%iter_method= "Pulay"
       ELSEIF(mixing_method==broyden_mixing_nr) THEN
          CALL broyden_mixing(qs_env,mixing_store, rho, para_env)
          mixing_store%iter_method= "Broy."
       ELSEIF(mixing_method==broyden_mixing_new_nr) THEN
          CPASSERT(.NOT. gapw)
          CALL broyden_mixing_new(mixing_store, rho, para_env)
          mixing_store%iter_method= "Broy."
       ELSEIF(mixing_method==multisecant_mixing_nr) THEN
          CPASSERT(.NOT. gapw)
          CALL multisecant_mixing(mixing_store, rho, para_env)
          mixing_store%iter_method= "MSec."
       END IF
    END IF

    IF(nspin==2) THEN
        CALL pw_axpy(rho_g(2)%pw,rho_g(1)%pw, 1.0_dp)
        CALL pw_scale(rho_g(1)%pw, 0.5_dp)
        CALL pw_axpy(rho_g(1)%pw,rho_g(2)%pw, -1.0_dp)
        CALL pw_scale(rho_g(2)%pw,-1.0_dp)
      CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_tmp%pw)
    END IF


    DO ispin=1,nspin
      CALL pw_transfer(rho_g(ispin)%pw,rho_r(ispin)%pw)
      tot_rho_r(ispin) = pw_integrate_function(rho_r(ispin)%pw,isign=-1)
!dbg
!  write(*,*) 'rho int 4', ispin , tot_rho_r(ispin)
!dbg
    END DO

    CALL timestop(handle)

  END SUBROUTINE gspace_mixing


! *****************************************************************************
!> \brief G-space mixing performed via the Kerker damping on the density on the grid
!>        thus affecting only the caluclation of the potential, but not the denisity matrix
!> \param qs_env ...
!> \param mixing_store ...
!> \param rho ...
!> \par History
!>      02.2009
!> \author MI
! *****************************************************************************

  SUBROUTINE gmix_potential_only (qs_env,mixing_store, rho)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho

    CHARACTER(len=*), PARAMETER :: routineN = 'gmix_potential_only', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), DIMENSION(:), POINTER       :: cc_new
    INTEGER                                  :: handle, iatom, ig, ispin, &
                                                natom, ng, nspin
    LOGICAL                                  :: failure, gapw
    REAL(dp)                                 :: alpha, f_mix
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    failure = .FALSE.
    CPASSERT(ASSOCIATED(mixing_store%rhoin))
    CPASSERT(ASSOCIATED(mixing_store%kerker_factor))

    CALL timeset(routineN,handle)
    NULLIFY(cc_new, dft_control, rho_g)

    CALL get_qs_env(qs_env, dft_control=dft_control)
    CALL qs_rho_get(rho, rho_g=rho_g)

    nspin = SIZE(rho_g,1)
    ng = SIZE(rho_g(1)%pw%pw_grid%gsq)

    gapw = dft_control%qs_control%gapw
    alpha = mixing_store%alpha

    DO ispin = 1,nspin
      cc_new =>  rho_g(ispin)%pw%cc
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
            rho_atom_set=rho_atom)
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
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
!> \par History
!>      03.2009
!> \author MI
! *****************************************************************************

  SUBROUTINE pulay_mixing(mixing_store, rho, para_env)

    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env

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
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g

    failure = .FALSE.
    CPASSERT(ASSOCIATED(mixing_store%res_buffer))
    CPASSERT(ASSOCIATED(mixing_store%rhoin_buffer))
    NULLIFY(rho_g)
    CALL timeset(routineN,handle)

    CALL qs_rho_get(rho, rho_g=rho_g)
    nspin = SIZE(rho_g,1)
    ng = SIZE(mixing_store%res_buffer(1,1)%cc)
    vol = rho_g(1)%pw%pw_grid%vol

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
    CPASSERT(istat==0)
    a=0.0_dp
    ALLOCATE(b(nb1,nb1), STAT=istat)
    CPASSERT(istat==0)
    b=0.0_dp
    ALLOCATE(c(nb,nb), STAT=istat)
    CPASSERT(istat==0)
    c=0.0_dp
    ALLOCATE(c_inv(nb,nb), STAT=istat)
    CPASSERT(istat==0)
    c_inv=0.0_dp
    ALLOCATE(alpha_c(nb), STAT=istat)
    CPASSERT(istat==0)
    alpha_c=0.0_dp
    norm_c_inv = 0.0_dp
    ALLOCATE(ev(nb1), STAT=istat)
    CPASSERT(istat==0)
    ev=0.0_dp
    ALLOCATE(cc_mix(ng), STAT=istat)
    CPASSERT(istat==0)

    DO ispin = 1,nspin
      mixing_store%res_buffer(ib,ispin)%cc(:) = CMPLX(0._dp,0._dp,KIND=dp)
      norm = 0.0_dp
      IF(nb==1) mixing_store%rhoin_buffer(1,ispin)%cc = mixing_store%rhoin(ispin)%cc
      DO ig = 1,ng
          f_mix = mixing_store%kerker_factor(ig)
 !         f_mix = 1.0_dp
          mixing_store%res_buffer(ib,ispin)%cc(ig) = f_mix*(rho_g(ispin)%pw%cc(ig)-&
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
           cc_mix(ig) = rho_g(ispin)%pw%cc(ig) -&
                         mixing_store%rhoin_buffer(ib,ispin)%cc(ig)
           rho_g(ispin)%pw%cc(ig) = f_mix*cc_mix(ig)+ &
                                        mixing_store%rhoin_buffer(ib,ispin)%cc(ig)
           mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
        END DO
      ELSE

        b(1:nb,1:nb) = mixing_store%pulay_matrix(1:nb, 1:nb)
        c(1:nb,1:nb) = b(1:nb,1:nb)
        b(nb1,1:nb)  = 1.0_dp
        b(1:nb,nb1)  = 1.0_dp
        b(nb1,nb1)   = 0.0_dp

        CALL diamat_all(b(1:nb1,1:nb1),ev(1:nb1))
        CALL invert_matrix(c,c_inv,inv_err)
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
            rho_g(ispin)%pw%cc(ig) =  f_mix*rho_g(ispin)%pw%cc(ig) +&
                      (1.0_dp-f_mix) *cc_mix(ig)
            mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
          END DO
        ELSE
          DO ig = 1,ng
            rho_g(ispin)%pw%cc(ig) =  cc_mix(ig)
            mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
          END DO
        END IF

      END IF  ! nb==1
    END DO ! ispin

    DEALLOCATE(a, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(b, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(ev, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(cc_mix, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(c,c_inv,alpha_c, STAT=istat)
    CPASSERT(istat==0)

    CALL timestop(handle)

  END SUBROUTINE pulay_mixing

! *****************************************************************************
!> \brief Broyden Mixing using as metrics for the residual the Kerer damping factor
!>        The mixing is applied directly on the density expansion in reciprocal space
!> \param qs_env ...
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
!> \par History
!>      03.2009
!> \author MI
! *****************************************************************************

  SUBROUTINE broyden_mixing(qs_env,mixing_store, rho, para_env)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env

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
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    failure = .FALSE.
    CPASSERT(ASSOCIATED(mixing_store%res_buffer))
    CPASSERT(ASSOCIATED(mixing_store%rhoin))
    CPASSERT(ASSOCIATED(mixing_store%rhoin_old))
    CPASSERT(ASSOCIATED(mixing_store%drho_buffer))
    NULLIFY(dft_control, rho_g)
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env, dft_control=dft_control)
    CALL qs_rho_get(rho, rho_g=rho_g)

    nspin = SIZE(rho_g,1)
    ng = SIZE(mixing_store%res_buffer(1,1)%cc)

    alpha=mixing_store%alpha
    beta=mixing_store%beta
    w0 = mixing_store%broy_w0
    nbuffer = mixing_store%nbuffer
    gapw = dft_control%qs_control%gapw

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
    CPASSERT(istat==0)
    a=0.0_dp
    ALLOCATE(b(ib-1,ib-1), STAT=istat)
    CPASSERT(istat==0)
    b=0.0_dp
    ALLOCATE(c(ib-1), STAT=istat)
    CPASSERT(istat==0)
    c=0.0_dp
    ALLOCATE(g(ib-1), STAT=istat)
    CPASSERT(istat==0)
    g=0.0_dp
    ALLOCATE(res_rho(ng), STAT=istat)
    CPASSERT(istat==0)

    IF(gapw) THEN
      CALL get_qs_env(qs_env=qs_env,&
            rho_atom_set=rho_atom)
      natom = SIZE(rho_atom)
    ELSE
      natom = 0
    END IF

    DO ispin = 1,nspin
      res_rho  = CMPLX(0.0_dp,0.0_dp,KIND=dp)
      DO ig = 1,ng
         res_rho(ig) = rho_g(ispin)%pw%cc(ig) - mixing_store%rhoin(ispin)%cc(ig)
      END DO

      IF(ib==1) THEN
        DO ig = 1,ng
          mixing_store%last_res(ispin)%cc(ig) = res_rho(ig)
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin(ispin)%cc(ig) + f_mix*res_rho(ig)
          mixing_store%rhoin_old(ispin)%cc(ig) =  mixing_store%rhoin(ispin)%cc(ig)
          mixing_store%rhoin(ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
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
        CALL invert_matrix(a,b,inv_err)

        CALL dgemv('T',IB-1,IB-1,1.0_dp,B,IB-1,C,1,0.0_dp,G,1)

        DO ig = 1,ng
          cc_mix = CMPLX(0.0_dp,0.0_dp,kind=dp)
          DO jb =  1,ib-1
            cc_mix = cc_mix - G(jb)*mixing_store%drho_buffer(jb,ispin)%cc(ig)
          END DO
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin(ispin)%cc(ig) +&
                   f_mix * res_rho(ig) + cc_mix
          mixing_store%rhoin_old(ispin)%cc(ig) =  mixing_store%rhoin(ispin)%cc(ig)
          mixing_store%rhoin(ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
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
    CPASSERT(istat==0)

    CALL timestop(handle)

 END SUBROUTINE broyden_mixing

! *****************************************************************************
!> \brief ...
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
! *****************************************************************************
 SUBROUTINE broyden_mixing_new(mixing_store, rho, para_env)

    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env

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
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g

    failure = .FALSE.
    CPASSERT(ASSOCIATED(mixing_store%rhoin_buffer))

    CALL timeset(routineN,handle)

    NULLIFY(delta_res, u_vec, z_vec)
    NULLIFY(fmat, rho_g)
    CALL qs_rho_get(rho, rho_g=rho_g)

    nspin = SIZE(rho_g,1)
    ng = SIZE(mixing_store%rhoin_buffer(1,1)%cc)
    vol = rho_g(1)%pw%pw_grid%vol
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
    CPASSERT(istat==0)
    ALLOCATE(res_rho_p(ng), STAT=istat)
    CPASSERT(istat==0)
    IF(ib>1) THEN
      ALLOCATE(a(ib-1,ib-1), STAT=istat)
      CPASSERT(istat==0)
      a=0.0_dp
      ALLOCATE(b(ib-1,ib-1), STAT=istat)
      CPASSERT(istat==0)
      b=0.0_dp
      ALLOCATE(bq(ib-1,ib-1), STAT=istat)
      CPASSERT(istat==0)
      bq=0.0_dp
      ALLOCATE(tmp(ng), STAT=istat)
      CPASSERT(istat==0)
      ALLOCATE(delta_res_p(ng),STAT=istat)
      CPASSERT(istat==0)
    END IF

    p_metric => mixing_store%p_metric
    weight => mixing_store%weight
    CPASSERT(ASSOCIATED(mixing_store%delta_res))
    delta_res => mixing_store%delta_res
    CPASSERT(ASSOCIATED(mixing_store%u_vec))
    u_vec => mixing_store%u_vec
    CPASSERT(ASSOCIATED(mixing_store%z_vec))
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
         res_rho(ig) = (rho_g(ispin)%pw%cc(ig) - mixing_store%rhoin_buffer(ib,ispin)%cc(ig))*sqt_uvol
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
          rho_g(ispin)%pw%cc(ig) = (mixing_store%rhoin_buffer(ib,ispin)%cc(ig) + f_mix*res_rho(ig))
          mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
!          rho_g(ispin)%pw%cc(ig) =  rho_g(ispin)%pw%cc(ig)*sqt_vol
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
          CALL invert_matrix(a,b,inv_err)
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
          work(:,:)=a
          ALLOCATE(iwork(8*(ib-1)),work_dgesdd(1))
             lwork=-1
             CALL DGESDD('S',ib-1,ib-1,work,ib-1,aval,au,ib-1,av,ib-1,work_dgesdd,lwork,iwork,info)
             lwork=INT(work_dgesdd(1))
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

!          CALL diamat_all(work(1:ib-1,1:ib-1),aval(1:ib-1))
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
        rho_g(ispin)%pw%cc(:) = CMPLX(0.0_dp,0.0_dp,KIND=dp)
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
            rho_g(ispin)%pw%cc(ig) = rho_g(ispin)%pw%cc(ig) - norm * z_vec(jb,ispin)%cc(ig)* sqt_vol

          END DO
        END DO

        DO ig = 1,ng
          f_mix = alpha*mixing_store%kerker_factor(ig)
          rho_g(ispin)%pw%cc(ig) =rho_g(ispin)%pw%cc(ig) +&
                   mixing_store%rhoin_buffer(ib,ispin)%cc(ig) + f_mix*res_rho(ig)
          mixing_store%rhoin_buffer(ibb,ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
!          rho_g(ispin)%pw%cc(ig) = rho_g(ispin)%pw%cc(ig)*sqt_vol
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
      CPASSERT(istat==0)
      DEALLOCATE(delta_res_p,tmp, STAT=istat)
      CPASSERT(istat==0)
    END IF
    DEALLOCATE(res_rho,res_rho_p, STAT=istat)
    CPASSERT(istat==0)

    CALL timestop(handle)

 END SUBROUTINE broyden_mixing_new

! *****************************************************************************
!> \brief Multisecant scheme to perform charge density Mixing
!>        as preconditioner we use the Kerker damping factor
!>        The mixing is applied directly on the density in reciprocal space,
!>        therefore it affects the potentials
!>        on the grid but not the density matrix
!> \param mixing_store ...
!> \param rho ...
!> \param para_env ...
!> \par History
!>      05.2009
!> \author MI
! *****************************************************************************
 SUBROUTINE multisecant_mixing(mixing_store,rho,para_env)

    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_para_env_type), POINTER          :: para_env

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
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g

    failure = .FALSE.
    CPASSERT(ASSOCIATED(mixing_store))

    CALL timeset(routineN,handle)

    NULLIFY(gn, rho_g)

    use_zgemm = .FALSE.
    use_zgemm_rev = .TRUE.
!   use_zgemm_rev = .FALSE.

    ! prepare the parameters
    CALL qs_rho_get(rho, rho_g=rho_g)

    nspin = SIZE(rho_g,1)
    ! not implemented for large grids.
    CPASSERT(rho_g(1)%pw%pw_grid%ngpts<HUGE(ng_global))
    ng_global = INT(rho_g(1)%pw%pw_grid%ngpts)
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
        gn(ig) = (rho_g(ispin)%pw%cc(ig)-mixing_store%rhoin_buffer(ib,ispin)%cc(ig))
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
          rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(1,ispin)%cc(ig) + &
              f_mix*mixing_store%res_buffer(1,ispin)%cc(ig)
          mixing_store%rhoin_buffer(ib_next,ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
        END DO
      END DO
      CALL timestop(handle)
      RETURN
    END IF

    ! allocate temporary arrays
    ! step_matrix  S ngxnb
    ALLOCATE(step_matrix(ng,nb),STAT=istat)
    CPASSERT(istat==0)
    ! res_matrix Y  ngxnb
    ALLOCATE(res_matrix(ng,nb),STAT=istat)
    CPASSERT(istat==0)
    ! matrix A  nbxnb
    ALLOCATE(a_matrix(nb,ng_global),STAT=istat)
    CPASSERT(istat==0)
    ! PSI nb vector of norms
    ALLOCATE(norm_res(nb),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(norm_res_low(nb),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(norm_res_up(nb),STAT=istat)
    CPASSERT(istat==0)

    ! matrix B   nbxnb
    ALLOCATE(b_matrix(nb,nb),STAT=istat)
    CPASSERT(istat==0)
    ! matrix B_inv   nbxnb
    ALLOCATE(binv_matrix(nb,nb),STAT=istat)
    CPASSERT(istat==0)

    ALLOCATE(gn_global(ng_global),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(res_matrix_global(ng_global),STAT=istat)
    CPASSERT(istat==0)
    IF(use_zgemm) THEN
      ALLOCATE(sa(ng,ng_global),STAT=istat)
      CPASSERT(istat==0)
      ALLOCATE(ya(ng,ng_global),STAT=istat)
      CPASSERT(istat==0)
    END IF
    IF(use_zgemm_rev) THEN
      ALLOCATE(tmp_vec(nb),STAT=istat)
      CPASSERT(istat==0)
    END IF
    ALLOCATE(pgn(ng),STAT=istat)
    CPASSERT(istat==0)
    ALLOCATE(ugn(ng),STAT=istat)
    CPASSERT(istat==0)


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
      g2 =>rho_g(1)%pw%pw_grid%gsq

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
      CALL invert_matrix(b_matrix,binv_matrix,inv_err)

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
         rho_g(ispin)%pw%cc(ig) = mixing_store%rhoin_buffer(ib,ispin)%cc(ig)  &
              -prec*step_size*ugn(ig) +  prec*pgn(ig)! - 0.1_dp * prec* gn(ig)
         mixing_store%rhoin_buffer(ib_next,ispin)%cc(ig) = rho_g(ispin)%pw%cc(ig)
      END DO

    END DO ! ispin

    ! Deallocate  temporary arrays
    DEALLOCATE(step_matrix,res_matrix, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(norm_res, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(norm_res_up, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(norm_res_low, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(a_matrix,b_matrix,binv_matrix, STAT=istat)
    CPASSERT(istat==0)
    DEALLOCATE(ugn,pgn, STAT=istat)
    CPASSERT(istat==0)
    IF(use_zgemm) THEN
      DEALLOCATE(sa,ya, STAT=istat)
      CPASSERT(istat==0)
    END IF
    IF(use_zgemm_rev) THEN
      DEALLOCATE(tmp_vec, STAT=istat)
      CPASSERT(istat==0)
    END IF
    DEALLOCATE(gn_global,res_matrix_global, STAT=istat)
    CPASSERT(istat==0)

    CALL timestop(handle)

 END SUBROUTINE multisecant_mixing

! *****************************************************************************
!> \brief  allocation needed when density mixing is used
!> \param qs_env ...
!> \param mixing_method ...
!> \param p_mix_new ...
!> \param p_delta ...
!> \param nspins ...
!> \param mixing_store ...
!> \par History
!>      05.2009 created [MI]
!>      08.2014 kpoints [JGH]
!>      02.2015 changed input to make g-space mixing available in linear scaling 
!>              SCF [Patrick Seewald]
!> \author fawzi
! *****************************************************************************
 SUBROUTINE mixing_allocate(qs_env, mixing_method, p_mix_new, p_delta, nspins, mixing_store)
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER                                  :: mixing_method
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: p_mix_new, p_delta
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(mixing_storage_type), POINTER       :: mixing_store

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mixing_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iat, ic, ispin, &
                                                natom, nbuffer, nimg, stat
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: refmatrix
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    failure = .FALSE.

    CALL timeset(routineN,handle)

    NULLIFY( matrix_s, dft_control, sab_orb, refmatrix, rho_atom)
    CALL get_qs_env(qs_env,&
                    sab_orb=sab_orb, &
                    matrix_s_kp=matrix_s, &
                    dft_control=dft_control)

    refmatrix => matrix_s(1,1)%matrix
    nimg = dft_control%nimages

!   *** allocate p_mix_new ***
    IF (PRESENT(p_mix_new)) THEN
      IF (.NOT.ASSOCIATED(p_mix_new)) THEN
         CALL cp_dbcsr_allocate_matrix_set(p_mix_new,nspins,nimg)
         DO i=1,nspins
            DO ic=1,nimg
               ALLOCATE(p_mix_new(i,ic)%matrix)
               CALL cp_dbcsr_init(p_mix_new(i,ic)%matrix)
               CALL cp_dbcsr_create(matrix=p_mix_new(i,ic)%matrix, &
                    name="SCF DENSITY", &
                    dist=cp_dbcsr_distribution(refmatrix), matrix_type=dbcsr_type_symmetric,&
                    row_blk_size=cp_dbcsr_row_block_sizes(refmatrix), &
                    col_blk_size=cp_dbcsr_col_block_sizes(refmatrix), &
                    nze=0)
               CALL cp_dbcsr_alloc_block_from_nbl(p_mix_new(i,ic)%matrix,sab_orb)
               CALL cp_dbcsr_set(p_mix_new(i,ic)%matrix,0.0_dp)
            ENDDO
         ENDDO
      END IF
    END IF

!   *** allocate p_delta ***
    IF (PRESENT(p_delta)) THEN
      IF (mixing_method>=gspace_mixing_nr) THEN
         IF(.NOT.ASSOCIATED(p_delta)) THEN
            CALL cp_dbcsr_allocate_matrix_set(p_delta,nspins,nimg)
            DO i=1,nspins
               DO ic=1,nimg
                  ALLOCATE(p_delta(i,ic)%matrix)
                  CALL cp_dbcsr_init(p_delta(i,ic)%matrix)
                  CALL cp_dbcsr_create(matrix=p_delta(i,ic)%matrix, &
                       name="SCF DENSITY", &
                       dist=cp_dbcsr_distribution(refmatrix), matrix_type=dbcsr_type_symmetric,&
                       row_blk_size=cp_dbcsr_row_block_sizes(refmatrix), &
                       col_blk_size=cp_dbcsr_col_block_sizes(refmatrix), &
                       nze=0)
                  CALL cp_dbcsr_alloc_block_from_nbl(p_delta(i,ic)%matrix,sab_orb)
                  CALL cp_dbcsr_set(p_delta(i,ic)%matrix,0.0_dp)
               ENDDO
            ENDDO
         END IF
         CPASSERT(ASSOCIATED(mixing_store))
      END IF
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
      IF(dft_control%qs_control%gapw) THEN
        CALL get_qs_env(qs_env=qs_env,&
               rho_atom_set=rho_atom)
        natom = SIZE(rho_atom)
        IF(.NOT. ASSOCIATED(mixing_store%paw)) THEN
          ALLOCATE(mixing_store%paw(natom),STAT=stat)
          CPASSERT(stat==0)
          mixing_store%paw = .FALSE.
          ALLOCATE(mixing_store%cpc_h_in(natom,nspins),STAT=stat)
          CPASSERT(stat==0)
          ALLOCATE(mixing_store%cpc_s_in(natom,nspins),STAT=stat)
          CPASSERT(stat==0)
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
        CPASSERT(stat==0)
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
        CPASSERT(stat==0)
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
        CPASSERT(stat==0)
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
        CPASSERT(stat==0)
        ALLOCATE(mixing_store%last_res(nspins), STAT=stat)
        CPASSERT(stat==0)
        DO ispin = 1,nspins
          DO i = 1,nbuffer
           NULLIFY(mixing_store%drho_buffer(i,ispin)%cc)
          END DO
          NULLIFY(mixing_store%last_res(ispin)%cc)
        END DO
      END IF
      IF(dft_control%qs_control%gapw) THEN
        IF(.NOT. ASSOCIATED(mixing_store%cpc_h_old)) THEN
          ALLOCATE(mixing_store%cpc_h_old(natom,nspins),STAT=stat)
          CPASSERT(stat==0)
          ALLOCATE(mixing_store%cpc_s_old(natom,nspins),STAT=stat)
          CPASSERT(stat==0)
          DO ispin = 1,nspins
            DO iat = 1,natom
               NULLIFY(mixing_store%cpc_h_old(iat,ispin)%r_coef)
               NULLIFY(mixing_store%cpc_s_old(iat,ispin)%r_coef)
            END DO
          END DO
        END IF
        IF(.NOT. ASSOCIATED(mixing_store%dcpc_h_in)) THEN
          ALLOCATE(mixing_store%dcpc_h_in(nbuffer,natom,nspins),STAT=stat)
        CPASSERT(stat==0)
          ALLOCATE(mixing_store%dcpc_s_in(nbuffer,natom,nspins),STAT=stat)
        CPASSERT(stat==0)
          ALLOCATE(mixing_store%cpc_h_lastres(natom,nspins),STAT=stat)
        CPASSERT(stat==0)
          ALLOCATE(mixing_store%cpc_s_lastres(natom,nspins),STAT=stat)
        CPASSERT(stat==0)
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
        CPASSERT(stat==0)
        ALLOCATE(mixing_store%delta_res(nbuffer-1,nspins), STAT=stat)
        CPASSERT(stat==0)
        ALLOCATE(mixing_store%u_vec(nbuffer-1,nspins), STAT=stat)
        CPASSERT(stat==0)
        ALLOCATE(mixing_store%z_vec(nbuffer-1,nspins), STAT=stat)
        CPASSERT(stat==0)
        ALLOCATE(mixing_store%weight(nbuffer,nspins), STAT=stat)
        CPASSERT(stat==0)
        ALLOCATE(mixing_store%fmat(nbuffer-1,nbuffer-1,nspins), STAT=stat)
        CPASSERT(stat==0)
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
        CPASSERT(stat==0)
      END IF
    END IF

    CALL timestop(handle)

 END SUBROUTINE mixing_allocate

! *****************************************************************************
!> \brief  initialiation needed when gspace mixing is used
!> \param mixing_method ...
!> \param rho ...
!> \param mixing_store ...
!> \param para_env ...
!> \param rho_atom ...
!> \par History
!>      05.2009 created [MI]
!> \author MI
! *****************************************************************************
 SUBROUTINE mixing_init(mixing_method,rho,mixing_store,para_env,rho_atom)
    INTEGER, INTENT(IN)                      :: mixing_method
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(mixing_storage_type), POINTER       :: mixing_store
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_atom

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
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g

    CALL timeset(routineN,handle)

    NULLIFY(g2,g_vec,rho_g)
    failure = .FALSE.
    CALL qs_rho_get(rho, rho_g=rho_g)

    nspin = SIZE(rho_g)
    ng = SIZE(rho_g(1)%pw%pw_grid%gsq,1)
    mixing_store%ig_max = ng
    g2 => rho_g(1)%pw%pw_grid%gsq
    g_vec => rho_g(1)%pw%pw_grid%g

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
      CPASSERT(istat==0)
    END IF
    IF(.NOT. ASSOCIATED(mixing_store%special_metric)) THEN
      ALLOCATE(mixing_store%special_metric(ng),STAT=istat)
      CPASSERT(istat==0)
    END IF
    beta = mixing_store%beta
    kmin = 0.1_dp
    mixing_store%kerker_factor = 1.0_dp
    mixing_store%special_metric = 1.0_dp
    ig1 = 1
    IF (rho_g(1)%pw%pw_grid%have_g0) ig1 = 2
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
      mixing_store%rhoin(ispin)%cc = rho_g(ispin)%pw%cc

      IF(ASSOCIATED(mixing_store%rhoin_buffer)) THEN
        IF(.NOT. ASSOCIATED(mixing_store%rhoin_buffer(1,ispin)%cc)) THEN
           DO ib = 1,nbuffer
              ALLOCATE(mixing_store%rhoin_buffer(ib,ispin)%cc(ng),STAT=istat)
              CPASSERT(istat==0)
            END DO
         END IF
         mixing_store%rhoin_buffer(1,ispin)%cc(1:ng) = &
                       rho_g(ispin)%pw%cc(1:ng)
      END IF
      IF(ASSOCIATED(mixing_store%res_buffer)) THEN
        IF(.NOT. ASSOCIATED(mixing_store%res_buffer(1,ispin)%cc)) THEN
          DO ib = 1,nbuffer
            ALLOCATE(mixing_store%res_buffer(ib,ispin)%cc(ng),STAT=istat)
            CPASSERT(istat==0)
          END DO
        END IF
      END IF
    END DO

    IF(nspin==2) THEN
      mixing_store%rhoin(1)%cc =  rho_g(1)%pw%cc +  rho_g(2)%pw%cc
      mixing_store%rhoin(2)%cc =  rho_g(1)%pw%cc -  rho_g(2)%pw%cc
      IF(ASSOCIATED(mixing_store%rhoin_buffer)) THEN
        mixing_store%rhoin_buffer(1,1)%cc =  rho_g(1)%pw%cc +  rho_g(2)%pw%cc
        mixing_store%rhoin_buffer(1,2)%cc =  rho_g(1)%pw%cc -  rho_g(2)%pw%cc
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
              CPASSERT(istat==0)
              ALLOCATE(mixing_store%cpc_h_in(iat,ispin)%r_coef(n1,n2),STAT=istat)
              CPASSERT(istat==0)
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
             CPASSERT(istat==0)
           END DO
           ALLOCATE(mixing_store%last_res(ispin)%cc(ng),STAT=istat)
           CPASSERT(istat==0)
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
               CPASSERT(istat==0)
               ALLOCATE(mixing_store%dcpc_s_in(ib,iat,ispin)%r_coef(n1,n2),STAT=istat)
               CPASSERT(istat==0)
             END DO
             ALLOCATE(mixing_store%cpc_h_lastres(iat,ispin)%r_coef(n1,n2),STAT=istat)
             CPASSERT(istat==0)
             ALLOCATE(mixing_store%cpc_s_lastres(iat,ispin)%r_coef(n1,n2),STAT=istat)
             CPASSERT(istat==0)
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
         CPASSERT(istat==0)
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
         IF(rho_g(1)%pw%pw_grid%have_g0) mixing_store%p_metric(1) = bconst
       END IF

    ELSEIF(mixing_method==broyden_mixing_new_nr) THEN
       DO ispin = 1,nspin
         IF(.NOT. ASSOCIATED(mixing_store%u_vec(1,ispin)%cc)) THEN
            DO ib = 1,nbuffer-1
              ALLOCATE(mixing_store%delta_res(ib,ispin)%cc(ng),STAT=istat)
              CPASSERT(istat==0)
              ALLOCATE(mixing_store%u_vec(ib,ispin)%cc(ng),STAT=istat)
               CPASSERT(istat==0)
              ALLOCATE(mixing_store%z_vec(ib,ispin)%cc(ng),STAT=istat)
              CPASSERT(istat==0)
            END DO
            ALLOCATE(mixing_store%last_res(ispin)%cc(ng),STAT=istat)
            CPASSERT(istat==0)
         END IF
       END DO
       IF(.NOT. ASSOCIATED(mixing_store%p_metric)) THEN
         ALLOCATE(mixing_store%p_metric(ng),STAT=istat)
         CPASSERT(istat==0)
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
         IF(rho_g(1)%pw%pw_grid%have_g0) mixing_store%p_metric(1) = bconst
!dbg
!        mixing_store%p_metric = 1.0_dp
!dbg
      END IF
    ELSEIF(mixing_method==multisecant_mixing_nr) THEN
      IF(.NOT. ASSOCIATED(mixing_store%ig_global_index)) THEN
         ALLOCATE(mixing_store%ig_global_index(ng),STAT=istat)
         CPASSERT(istat==0)
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
