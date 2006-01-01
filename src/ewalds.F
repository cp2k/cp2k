!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/ewalds [1.0] *
!!
!!   NAME
!!     ewalds
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     JGH (15-Mar-2001) : New routine ewald_setup (former pme_setup)
!!     JGH (23-Mar-2001) : Get rid of global variable ewald_grp
!!
!!   SOURCE
!******************************************************************************

MODULE ewalds

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dg_types,                        ONLY: dg_type, dg_get
  USE dg_rho0_types,                   ONLY: dg_rho0_type
  USE ewald_environment_types,         ONLY: ewald_environment_type, &
                                             ewald_env_get
  USE ewald_pw_types,                  ONLY: ewald_pw_type, &
                                             ewald_pw_get
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp, default_string_length
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_types,                        ONLY: pw_type
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE simulation_cell,                 ONLY: cell_type
  USE structure_factor_types,          ONLY: structure_factor_type
  USE structure_factors,               ONLY: structure_factor_allocate,&
                                             structure_factor_deallocate,&
                                             structure_factor_evaluate
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: matvec_3x3
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ewald_evaluate, ewald_correction, ewald_self
  PUBLIC :: ewald_print
!******************************************************************************

  CONTAINS

!******************************************************************************
!!****** ewald_pw_methods/ewald_evaluate [1.0] *
!!
!!   NAME
!!     ewald_evaluate
!!
!!   FUNCTION
!!     computes the potential and the force from the g-space part of
!!     the 1/r potential
!!     Ref.: J.-P. Hansen, Enrico Fermi School, 1985         
!!     Note: Only the positive G-vectors are used in the sum.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (28-Feb-2001) : added fft_usage flag to pw_grid_setup and calls to
!!                         pw_grid_construct
!!     JGH (21-Feb-2001) : changed name
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_evaluate (ewald_env, ewald_pw, atomic_kind_set, particle_set,  &
                           local_particles, fg_coulomb, vg_coulomb, pv_g )
    TYPE ( ewald_environment_type ), POINTER   :: ewald_env
    TYPE ( atomic_kind_type ), POINTER        :: atomic_kind_set ( : )
    TYPE ( ewald_pw_type ), POINTER :: ewald_pw
    TYPE ( particle_type ), POINTER :: particle_set ( : )
    TYPE (distribution_1d_type), POINTER :: local_particles
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: fg_coulomb
    REAL(KIND=dp), INTENT(OUT)                   :: vg_coulomb
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)  :: pv_g
! locals
    TYPE ( cell_type ), POINTER :: cell
    COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: summe
    INTEGER                                  :: gpt, handle, isos, lp, mp, &
                                                node, np, iparticle_kind, &
                                                nparticle_kind, iparticle, &
                                                iparticle_local, nparticle_local, &
                                                group, nnodes
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(KIND=dp)                                :: common, denom, &
                                                e_igdotr, flops, &
                                                four_alpha_sq, gauss, gsq, &
                                                pref, alpha, eps0, q
    REAL(KIND=dp), DIMENSION(3)                  :: vec
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER   :: rho0
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(structure_factor_type)              :: exp_igr
    TYPE ( pw_pool_type ), POINTER           :: pw_pool
    TYPE ( pw_type ), POINTER                :: pw
    TYPE ( cp_error_type ), POINTER          :: error
    TYPE ( dg_type ), POINTER                :: dg
    TYPE ( dg_rho0_type ), POINTER           :: dg_rho0
    REAL (KIND=dp), DIMENSION ( : ), POINTER   :: charge

!------------------------------------------------------------------------------

  CALL timeset ( 'EWALD', 'E', 'Mflops', handle )
  flops = 0.0_dp

! pointing
  CALL ewald_env_get (ewald_env, alpha=alpha, eps0 = eps0, group = group )
  CALL ewald_pw_get (ewald_pw, pw_big_pool=pw_pool, dg = dg, &
                     cell = cell )
  CALL dg_get ( dg, dg_rho0=dg_rho0 )
  rho0 => dg_rho0 % density % pw % cr3d
  pw_grid => pw_pool % pw_grid
  bds =>  pw_grid % bounds

! allocating
  nparticle_kind = SIZE ( atomic_kind_set )
  nnodes = 0
  DO iparticle_kind = 1, nparticle_kind
     nnodes = nnodes + local_particles%n_el(iparticle_kind) 
  ENDDO 

  CALL structure_factor_allocate ( pw_grid % bounds, nnodes, exp_igr)

  ALLOCATE ( summe ( 1:pw_grid % ngpts_cut ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald', 'summe', 0 )

  ALLOCATE ( charge ( 1:nnodes ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald', 'charge', 0 )

! Initializing vg_coulomb and fg_coulomb
  vg_coulomb = 0.0_dp
  fg_coulomb = 0.0_dp
  pv_g = 0.0_dp
! defining four_alpha_sq
  four_alpha_sq = 4.0_dp * alpha ** 2
! zero node count
  node = 0
  DO iparticle_kind = 1, nparticle_kind
     atomic_kind => atomic_kind_set(iparticle_kind)
     CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q)
     nparticle_local = local_particles%n_el(iparticle_kind)
     DO iparticle_local=1,nparticle_local
        node = node + 1
        charge ( node ) = q
        iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
        vec = matvec_3x3 ( cell % h_inv, particle_set ( iparticle ) % r )
        CALL structure_factor_evaluate ( vec, pw_grid % npts, exp_igr % lb, &
        exp_igr % ex(:,node), exp_igr % ey(:,node), exp_igr % ez(:,node) )
     END DO
  END DO
  
  summe ( : ) = CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
! looping over the positive g-vectors
  DO gpt = 1, pw_grid % ngpts_cut

     lp = pw_grid % mapl % pos ( pw_grid % g_hat ( 1, gpt ) )
     mp = pw_grid % mapm % pos ( pw_grid % g_hat ( 2, gpt ) )
     np = pw_grid % mapn % pos ( pw_grid % g_hat ( 3, gpt ) )

     lp = lp + bds ( 1, 1 )
     mp = mp + bds ( 1, 2 )
     np = np + bds ( 1, 3 )

! initializing sum to be used in the energy and force
     DO node = 1, nnodes
       summe ( gpt ) = summe ( gpt ) + charge ( node ) * &
            ( exp_igr % ex ( lp, node ) &
            * exp_igr % ey ( mp, node ) &
            * exp_igr % ez ( np, node ) )
     END DO
  END DO
  flops = flops &
       + pw_grid % ngpts_cut * 10.0_dp &
       + pw_grid % ngpts_cut * node * 20.0_dp

#if defined(__parallel)
  CALL mp_sum ( summe, group )
#endif

! looping over the positive g-vectors
  DO gpt = 1, pw_grid % ngpts_cut
! computing the potential energy
     lp = pw_grid % mapl % pos ( pw_grid % g_hat ( 1, gpt ) )
     mp = pw_grid % mapm % pos ( pw_grid % g_hat ( 2, gpt ) )
     np = pw_grid % mapn % pos ( pw_grid % g_hat ( 3, gpt ) )

     lp = lp + bds ( 1, 1 )
     mp = mp + bds ( 1, 2 )
     np = np + bds ( 1, 3 )

     IF ( pw_grid % gsq ( gpt ) <= 1.0E-10_dp ) CYCLE

     gauss = ( rho0 ( lp, mp, np ) * pw_grid % vol ) ** 2 &
          / pw_grid % gsq ( gpt )
     common = gauss * REAL ( summe ( gpt ) * CONJG ( summe ( gpt ) ) )
     vg_coulomb = vg_coulomb + common

! computing the force
    node = 0
    DO node = 1, nnodes
      e_igdotr = AIMAG ( summe(gpt) * CONJG &
         ( exp_igr % ex ( lp, node ) &
         * exp_igr % ey ( mp, node ) &
         * exp_igr % ez ( np, node ) ) )
      fg_coulomb ( :, node ) = fg_coulomb ( :, node ) &
      + charge ( node ) * gauss * e_igdotr * pw_grid % g ( :, gpt )
    END DO

! compute the virial P*V

     denom = 1.0_dp / four_alpha_sq + 1.0_dp / pw_grid % gsq ( gpt )
     pv_g(1,1) = pv_g(1,1) + common * ( 1.0_dp &
          - 2.0_dp * pw_grid % g ( 1, gpt ) * pw_grid % g ( 1, gpt ) * denom )
     pv_g(1,2) = pv_g(1,2) - common * &
          ( 2.0_dp * pw_grid % g ( 1, gpt ) * pw_grid % g ( 2, gpt ) * denom )
     pv_g(1,3) = pv_g(1,3) - common * &
          ( 2.0_dp * pw_grid % g ( 1, gpt ) * pw_grid % g ( 3, gpt ) * denom )
     pv_g(2,1) = pv_g(2,1) - common * &
          ( 2.0_dp * pw_grid % g ( 2, gpt ) * pw_grid % g ( 1, gpt ) * denom )
     pv_g(2,2) = pv_g(2,2) + common * ( 1.0_dp &
          - 2.0_dp * pw_grid % g ( 2, gpt ) * pw_grid % g ( 2, gpt ) * denom )
     pv_g(2,3) = pv_g(2,3) - common * &
          ( 2.0_dp * pw_grid % g ( 2, gpt ) * pw_grid % g ( 3, gpt ) * denom )
     pv_g(3,1) = pv_g(3,1) - common * &
          ( 2.0_dp * pw_grid % g ( 3, gpt ) * pw_grid % g ( 1, gpt ) * denom )
     pv_g(3,2) = pv_g(3,2) - common * &
          ( 2.0_dp * pw_grid % g ( 3, gpt ) * pw_grid % g ( 2, gpt ) * denom )
     pv_g(3,3) = pv_g(3,3) + common * ( 1.0_dp &
          - 2.0_dp * pw_grid % g ( 3, gpt ) * pw_grid % g ( 3, gpt ) * denom )
  END DO

  flops = flops + pw_grid % ngpts_cut * 55.0_dp &
       + pw_grid % ngpts_cut * node * 23.0_dp

  pref = 1.0_dp / eps0 / pw_grid % vol 
  vg_coulomb = vg_coulomb * pref
  pv_g = pv_g * pref

  pref = 2.0_dp * pref
  fg_coulomb = fg_coulomb * pref

  CALL structure_factor_deallocate ( exp_igr )

  DEALLOCATE ( charge, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald', 'charge' )
  NULLIFY ( charge )

  DEALLOCATE ( summe, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald', 'summe' )

  flops = flops * 1.0E-6_dp
  CALL timestop ( flops, handle )

END SUBROUTINE ewald_evaluate

!!*****
!******************************************************************************
!!****** ewalds/ewald_self [1.0] *
!!
!!   NAME
!!     ewald_self
!!
!!   FUNCTION
!     Computes the self interaction from g-space 
!!     and the neutralizing background
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_self ( ewald_env, atomic_kind_set,  &
                        local_particles, e_self, e_neut  )

    TYPE ( ewald_environment_type ), POINTER  :: ewald_env
    TYPE ( atomic_kind_type ), POINTER        :: atomic_kind_set ( : )
    TYPE ( distribution_1d_type ), POINTER    :: local_particles 
    REAL (KIND=dp), INTENT ( OUT )              :: e_self, e_neut

    INTEGER                                  :: nparticle_kind, group
    INTEGER                                  :: iparticle_kind, nparticle_local
    REAL(KIND=dp)                                :: q, q_self, q_sum, eps0, alpha
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    INTEGER     :: ewald_type

!------------------------------------------------------------------------------


  CALL ewald_env_get ( ewald_env, ewald_type = ewald_type,  &
                       eps0 = eps0, alpha = alpha, group = group ) 
  q_self = 0.0_dp
  q_sum = 0.0_dp
  nparticle_kind = SIZE ( atomic_kind_set )
  DO iparticle_kind=1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q)
        nparticle_local = local_particles%n_el(iparticle_kind)
        q_self = q_self + q*q*nparticle_local
        q_sum = q_sum + q*nparticle_local
  END DO

#if defined(__parallel)
  CALL mp_sum ( q_self, group )
  CALL mp_sum ( q_sum, group )
#endif

  e_neut = 0.0_dp
  e_self = 0.0_dp
  IF ( ewald_type /= do_ewald_none) THEN

     e_self = -q_self *  alpha &
          / ( 4.0_dp *  eps0 * pi * SQRT ( pi ) )

     e_neut = -q_sum ** 2 / ( 8.0_dp * eps0 * alpha ** 2 )
  END IF

END SUBROUTINE ewald_self

!!*****
!******************************************************************************
!!****** ewalds/ewald_correction [1.0] *
!!
!!   NAME
!!     ewald_correction
!!
!!   FUNCTION
!!     computing the shift in the potential energy 
!!     for the real space coulomb potential
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

FUNCTION ewald_correction ( ewald_env, qi, qj, rcutsq ) &
     RESULT ( e_cutoff_coul )


    TYPE(ewald_environment_type), POINTER    :: ewald_env
    REAL(KIND=dp), INTENT(IN)                    :: qi, qj, rcutsq
    REAL(KIND=dp)                                :: e_cutoff_coul

    REAL(KIND=dp), PARAMETER :: ac1 = 0.254829592_dp, ac2 = -0.284496736_dp, &
      ac3 = 1.421413741_dp, ac4 = -1.453152027_dp, ac5 = 1.061405429_dp, &
      pc = 0.3275911_dp

    REAL(KIND=dp)                                :: arg, e_arg_arg, &
                                                erfc_over_rij, erfcf, &
                                                inv_fourpieps0, tc, eps0, alpha

!------------------------------------------------------------------------------
! compute the value of the real-space coulomb at the cut-off

  CALL ewald_env_get ( ewald_env, eps0 = eps0, alpha = alpha )
  inv_fourpieps0 = 1.0_dp / ( 4.0_dp * pi * eps0 )
  arg = alpha * SQRT ( rcutsq )
  e_arg_arg = EXP ( - arg ** 2 )
  tc = 1.0_dp / ( 1.0_dp + pc * arg )
  erfcf = ((((ac5*tc+ac4)*tc+ac3)*tc+ac2)*tc+ac1)*tc * e_arg_arg
  erfc_over_rij = erfcf / SQRT ( rcutsq )
  e_cutoff_coul = inv_fourpieps0 * qi * qj * erfc_over_rij

END FUNCTION ewald_correction

!!*****
!******************************************************************************
!!****** ewalds/ewald_print [1.0] *
!!
!!   NAME
!!     ewald_print
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_print ( iw, cell, e_gspace, e_self, e_neut, e_bonded )


    INTEGER, INTENT(IN)                      :: iw
    TYPE(cell_type), INTENT(IN)              :: cell
    REAL(KIND=dp), INTENT(IN)                    :: e_gspace, e_self, e_neut, &
                                                e_bonded

    TYPE(unit_convert_type), POINTER         :: units

!------------------------------------------------------------------------------

  CALL get_cp2k_units ( units )

  WRITE ( iw, '( A, A )' ) ' *********************************', &
       '**********************************************'
  WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' INITIAL GSPACE ENERGY', &
       units % e_label, '= ', e_gspace
  WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' SELF ENERGY CORRECTION', &
       units % e_label, '= ', e_self
  WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' NEUT. BACKGROUND', units % e_label, &
       '= ', e_neut / cell % deth
  WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' BONDED CORRECTION', units % e_label, &
       '= ', e_bonded
  WRITE ( iw, '( A, A )' ) ' *********************************', &
       '**********************************************'

END SUBROUTINE ewald_print

!******************************************************************************

END MODULE ewalds

!******************************************************************************
