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
  USE coefficient_types,               ONLY: coeff_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE dg_rho0s,                        ONLY: dg_rho0_setup
  USE dg_types,                        ONLY: dg_type
  USE dgs,                             ONLY: dg_pme_grid_setup
  USE ewald_parameters_types,          ONLY: ewald_parameters_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE molecule_types,                  ONLY: particle_node_type
  USE particle_lists,                  ONLY: particle_list_select
  USE particle_types,                  ONLY: particle_type
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_construct,&
                                             pw_grid_setup
  USE qs_interactions,                 ONLY: exp_radius
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_grid_setup
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

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ewald_evaluate, ewald_correction, ewald_self, ewald_setup
  PUBLIC :: ewald_print, ewald_initialize

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** ewalds/ewald_evaluate [1.0] *
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

SUBROUTINE ewald_evaluate (dg, ewald_param, fg_coulomb, vg_coulomb, &
     pv_g, pnode, box )
  

    TYPE(dg_type), INTENT(IN), TARGET        :: dg
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: fg_coulomb
    REAL(dbl), INTENT(OUT)                   :: vg_coulomb
    REAL(dbl), DIMENSION(:, :), INTENT(OUT)  :: pv_g
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(cell_type), INTENT(IN)              :: box

    COMPLEX(dbl), ALLOCATABLE, DIMENSION(:)  :: summe
    INTEGER                                  :: gpt, handle, isos, lp, mp, &
                                                nnodes, node, np
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(dbl)                                :: charge, common, denom, &
                                                e_igdotr, flops, &
                                                four_alpha_sq, gauss, gsq, &
                                                pref
    REAL(dbl), DIMENSION(3)                  :: vec
    REAL(dbl), DIMENSION(:, :, :), POINTER   :: rho0
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(coeff_type), POINTER                :: dg_rho0
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(structure_factor_type)              :: exp_igr

!------------------------------------------------------------------------------

  CALL timeset ( 'EWALD', 'E', 'Mflops', handle )
  flops = 0.0_dbl

! pointing
  dg_rho0 => dg % dg_rho0 % density
  rho0 => dg % dg_rho0 % density % pw % cr3d
  pw_grid => dg % dg_rho0 % density % pw % pw_grid
  bds => dg % dg_rho0 % density % pw % pw_grid % bounds

! allocating
  nnodes = SIZE ( pnode )

  CALL structure_factor_allocate ( pw_grid % bounds, nnodes, exp_igr)

  ALLOCATE ( summe ( 1:pw_grid % ngpts_cut ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald', 'summe', 0 )

! Initializing vg_coulomb and fg_coulomb
  vg_coulomb = 0.0_dbl
  fg_coulomb = 0.0_dbl
  pv_g = 0.0_dbl
! defining four_alpha_sq
  four_alpha_sq = 4.0_dbl * ewald_param % alpha ** 2

  DO node = 1, nnodes

     vec = matvec_3x3 ( box % h_inv, pnode ( node ) % p % r )

     CALL structure_factor_evaluate ( vec, pw_grid % npts, exp_igr % lb, &
          exp_igr % ex(:,node), exp_igr % ey(:,node), exp_igr % ez(:,node) )

  END DO
  
  summe ( : ) = CMPLX ( 0.0_dbl, 0.0_dbl, dbl )
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
        atomic_kind => pnode (node) % p % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=charge)
        summe ( gpt ) = summe ( gpt ) + charge * &
             ( exp_igr % ex ( lp, node ) &
             * exp_igr % ey ( mp, node ) &
             * exp_igr % ez ( np, node ) )
     END DO
  END DO
  flops = flops &
       + pw_grid % ngpts_cut * 10.0_dbl &
       + pw_grid % ngpts_cut * nnodes * 20.0_dbl

#if defined(__parallel)
  CALL mp_sum ( summe, ewald_param % group )
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

     IF ( pw_grid % gsq ( gpt ) <= 1.0E-10_dbl ) CYCLE

     gauss = ( rho0 ( lp, mp, np ) * pw_grid % vol ) ** 2 &
          / pw_grid % gsq ( gpt )
     common = gauss * REAL ( summe ( gpt ) * CONJG ( summe ( gpt ) ) )
     vg_coulomb = vg_coulomb + common

! computing the force
     DO node = 1, nnodes
        e_igdotr = AIMAG ( summe(gpt) * CONJG &
             ( exp_igr % ex ( lp, node ) &
             * exp_igr % ey ( mp, node ) &
             * exp_igr % ez ( np, node ) ) )
        atomic_kind => pnode (node) % p % atomic_kind
        CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=charge)
        charge =  charge * gauss * e_igdotr
        fg_coulomb ( :, node ) = fg_coulomb ( :, node ) &
             + charge * pw_grid % g ( :, gpt )
     END DO

! compute the virial P*V

     denom = 1.0_dbl / four_alpha_sq + 1.0_dbl / pw_grid % gsq ( gpt )
     pv_g(1,1) = pv_g(1,1) + common * ( 1.0_dbl &
          - 2.0_dbl * pw_grid % g ( 1, gpt ) * pw_grid % g ( 1, gpt ) * denom )
     pv_g(1,2) = pv_g(1,2) - common * &
          ( 2.0_dbl * pw_grid % g ( 1, gpt ) * pw_grid % g ( 2, gpt ) * denom )
     pv_g(1,3) = pv_g(1,3) - common * &
          ( 2.0_dbl * pw_grid % g ( 1, gpt ) * pw_grid % g ( 3, gpt ) * denom )
     pv_g(2,1) = pv_g(2,1) - common * &
          ( 2.0_dbl * pw_grid % g ( 2, gpt ) * pw_grid % g ( 1, gpt ) * denom )
     pv_g(2,2) = pv_g(2,2) + common * ( 1.0_dbl &
          - 2.0_dbl * pw_grid % g ( 2, gpt ) * pw_grid % g ( 2, gpt ) * denom )
     pv_g(2,3) = pv_g(2,3) - common * &
          ( 2.0_dbl * pw_grid % g ( 2, gpt ) * pw_grid % g ( 3, gpt ) * denom )
     pv_g(3,1) = pv_g(3,1) - common * &
          ( 2.0_dbl * pw_grid % g ( 3, gpt ) * pw_grid % g ( 1, gpt ) * denom )
     pv_g(3,2) = pv_g(3,2) - common * &
          ( 2.0_dbl * pw_grid % g ( 3, gpt ) * pw_grid % g ( 2, gpt ) * denom )
     pv_g(3,3) = pv_g(3,3) + common * ( 1.0_dbl &
          - 2.0_dbl * pw_grid % g ( 3, gpt ) * pw_grid % g ( 3, gpt ) * denom )
  END DO

  flops = flops + pw_grid % ngpts_cut * 55.0_dbl &
       + pw_grid % ngpts_cut * nnodes * 23.0_dbl

  pref = 1.0_dbl / ( ewald_param % eps0 * pw_grid % vol )
  vg_coulomb = vg_coulomb * pref
  pv_g = pv_g * pref

  pref = 2.0_dbl * pref
  fg_coulomb = fg_coulomb * pref

  CALL structure_factor_deallocate ( exp_igr )

  DEALLOCATE ( summe, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald', 'summe' )

  flops = flops * 1.0E-6_dbl
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
!!     Computes the self interaction from g-space 
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

SUBROUTINE ewald_self ( ewald_param, pnode, e_self, e_neut  )

    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    REAL(dbl), INTENT(OUT)                   :: e_self, e_neut

    INTEGER                                  :: i, nnodes
    REAL(dbl)                                :: q, q_self, q_sum
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!------------------------------------------------------------------------------

  q_self = 0.0_dbl
  q_sum = 0.0_dbl
  nnodes = SIZE ( pnode )
  DO i = 1, nnodes
     atomic_kind => pnode (i) % p % atomic_kind
     CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q)
     q_self = q_self + q*q
     q_sum = q_sum + q
  END DO

#if defined(__parallel)
  CALL mp_sum(q_self,ewald_param%group)
  CALL mp_sum(q_sum,ewald_param%group)
#endif

  IF ( ewald_param % ewald_type == 'EWALD' .OR. &
       ewald_param % ewald_type == 'PME' .OR. &
       ewald_param % ewald_type == 'SPME' ) THEN

     e_self = -q_self * ewald_param % alpha &
          / ( 4.0_dbl * ewald_param % eps0 * pi * SQRT ( pi ) )

     e_neut = -q_sum ** 2 / ( 8.0_dbl * ewald_param % eps0 &
          * ewald_param % alpha ** 2 )
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

FUNCTION ewald_correction ( ewald_param, qi, qj, rcutsq ) &
     RESULT ( e_cutoff_coul )


    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    REAL(dbl), INTENT(IN)                    :: qi, qj, rcutsq
    REAL(dbl)                                :: e_cutoff_coul

    REAL(dbl), PARAMETER :: ac1 = 0.254829592_dbl, ac2 = -0.284496736_dbl, &
      ac3 = 1.421413741_dbl, ac4 = -1.453152027_dbl, ac5 = 1.061405429_dbl, &
      pc = 0.3275911_dbl

    REAL(dbl)                                :: arg, e_arg_arg, &
                                                erfc_over_rij, erfcf, &
                                                inv_fourpieps0, tc

!------------------------------------------------------------------------------
! compute the value of the real-space coulomb at the cut-off

  inv_fourpieps0 = 1.0_dbl / ( 4.0_dbl * pi * ewald_param % eps0 )
  arg = ewald_param % alpha * SQRT ( rcutsq )
  e_arg_arg = EXP ( - arg ** 2 )
  tc = 1.0_dbl / ( 1.0_dbl + pc * arg )
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

SUBROUTINE ewald_print ( iw, box, e_gspace, e_self, e_neut, e_bonded )


    INTEGER, INTENT(IN)                      :: iw
    TYPE(cell_type), INTENT(IN)              :: box
    REAL(dbl), INTENT(IN)                    :: e_gspace, e_self, e_neut, &
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
       '= ', e_neut / box % deth
  WRITE ( iw, '( A, A, T35, A, T66, E15.7 )' ) ' BONDED CORRECTION', units % e_label, &
       '= ', e_bonded
  WRITE ( iw, '( A, A )' ) ' *********************************', &
       '**********************************************'

END SUBROUTINE ewald_print

!!*****
!******************************************************************************
!!****** ewalds/ewald_initialize [1.0] *
!!
!!   NAME
!!     ewald_initialize
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (12-Jan-2001): Added SPME part
!!     JGH (15-Mar-2001): Work newly distributed between initialize, setup,
!!                        and force routine
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_initialize ( dg, part, pnode, pnode_grp, ewald_param, box, &
          e_self, e_neut, ewald_grid, pme_small_grid, pme_big_grid, spme_grid )
  

    TYPE(dg_type), INTENT(OUT)               :: dg
    TYPE(particle_type), DIMENSION(:), &
      INTENT(IN)                             :: part
    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN)                             :: pnode
    TYPE(global_environment_type), &
      INTENT(IN)                             :: pnode_grp
    TYPE(ewald_parameters_type), &
      INTENT(INOUT)                          :: ewald_param
    TYPE(cell_type), INTENT(INOUT)           :: box
    REAL(dbl), INTENT(OUT)                   :: e_self, e_neut
    TYPE(pw_grid_type), INTENT(OUT), &
      OPTIONAL                               :: ewald_grid, pme_small_grid, &
                                                pme_big_grid, spme_grid

    INTEGER                                  :: gmax, iw, npts_s( 3 )
    REAL(dbl)                                :: alphasq, cutoff, &
                                                cutoff_radius, ecut, norm
    TYPE(realspace_grid_type)                :: rs

!------------------------------------------------------------------------------
! parallelisation is over atoms (pnodes), so the group of processors
! has to be the same as the group for the pnodes
! This holds only for normal Ewald sums!
! For the grid based methods, G-Space/Real Space are distributed differently

  ewald_param % group = pnode_grp % group

! fire up the reciprocal space and compute self interaction and
! term from the neutralizing background.

  IF ( ewald_param % ewald_type /= 'NONE' ) THEN

     CALL ewald_self ( ewald_param, pnode, e_self, e_neut )
     
     IF ( PRESENT ( ewald_grid ) ) THEN
        ! set up EWALD "Classic"
        CALL pw_grid_construct( ewald_grid )
        gmax = ewald_param % gmax
        IF ( gmax == 2 * ( gmax / 2 ) ) THEN
           IF ( pnode_grp % ionode ) &
             CALL stop_program ( "initialize_ewalds", "gmax has to be odd" )
        END IF
        ewald_grid % bounds ( 1, : ) = -gmax / 2
        ewald_grid % bounds ( 2, : ) = +gmax / 2
        npts_s = (/ gmax, gmax, gmax /)
        ewald_grid % grid_span = HALFSPACE

        CALL pw_find_cutoff ( npts_s, box, cutoff )

        ecut = 0.5_dbl * cutoff * cutoff
        CALL pw_grid_setup( box, ewald_grid, cutoff = ecut, &
                            info = iw, fft_usage = .FALSE. )
     
     ELSEIF ( PRESENT ( pme_small_grid ) .AND. PRESENT ( pme_big_grid ) ) THEN
        ! set up PME, the double grid method
        CALL pw_grid_construct( pme_small_grid )
        CALL pw_grid_construct( pme_big_grid )

        npts_s ( : ) = ewald_param % ns_max
! compute cut-off radius
        alphasq = ( ewald_param % alpha ) ** 2
        norm = ( 2.0_dbl * alphasq / pi ) ** ( 1.5_dbl )
        alphasq = ( ewald_param % alpha ) ** 2
        norm = ( 2.0_dbl * alphasq / pi ) ** ( 1.5_dbl )
        cutoff_radius = exp_radius ( 0, 2._dbl * alphasq, ewald_param % epsilon , norm )

        CALL dg_pme_grid_setup ( box, npts_s, cutoff_radius, &
                        pme_small_grid, pme_big_grid, pnode_grp )

        CALL rs_grid_setup ( rs, pme_big_grid, MAXVAL ( pme_small_grid % npts ), &
             pnode_grp % scr )

        ! set up a list with all particles having a charge
        CALL particle_list_select ( pnode, dg % plist, dg % nparts )

        IF ( pnode_grp % ionode ) THEN
          WRITE ( iw, '( A,T71,E10.4 )' ) &
             ' EWALD| Gaussian tolerance (effective) ', ewald_param % epsilon
          WRITE ( iw, '( A,T63,3I6 )' ) &
             ' EWALD| Small box grid ', pme_small_grid % npts
          WRITE ( iw, '( A,T63,3I6 )' ) &
             ' EWALD| Full box grid ', pme_big_grid % npts
        END IF

     ELSEIF ( PRESENT ( spme_grid ) ) THEN
        ! set up SPME
        CALL pw_grid_construct( spme_grid )
        gmax = ewald_param % gmax
        npts_s = (/ gmax, gmax, gmax /)
        CALL pw_find_cutoff ( npts_s, box, cutoff )
        spme_grid % grid_span = HALFSPACE
        spme_grid % bounds ( 1, : ) = -gmax/2
        spme_grid % bounds ( 2, : ) = -gmax/2 + gmax - 1

        ecut = 0.5_dbl * cutoff * cutoff
        CALL pw_grid_setup( box, spme_grid, cutoff = ecut, &
             pe_group = pnode_grp % group, info = pnode_grp % scr, &
             fft_usage = .TRUE. )

        npts_s ( : ) = ewald_param % o_spline
        CALL rs_grid_setup ( rs, spme_grid, MAXVAL ( npts_s ), pnode_grp % scr )

        ! set up a list with all particles having a charge
        CALL particle_list_select ( pnode, dg % plist, dg % nparts )

     END IF

  END IF

END SUBROUTINE ewald_initialize

!!*****
!******************************************************************************
!!****** ewalds/ewald_setup [1.0] *
!!
!!   NAME
!!     ewald_setup
!!
!!   FUNCTION
!!     Calculates the Fourier transform of the "Ewald function"
!!
!!   AUTHOR
!!     JGH (15-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ewald_setup ( pw_grid, ewald_param, dg )


    TYPE(pw_grid_type), INTENT(IN)           :: pw_grid
    TYPE(ewald_parameters_type), INTENT(IN)  :: ewald_param
    TYPE(dg_type), INTENT(OUT)               :: dg

    INTEGER                                  :: isos

!------------------------------------------------------------------------------
! This is the first (and only) double grid

  dg % grid_index = 1

  ALLOCATE ( dg % dg_rho0 % zet (1) , STAT = isos)
  IF ( isos /= 0 ) CALL stop_memory ('ewald_setup', 'zet', 0 )

! No contracted Gaussians are used here
  NULLIFY ( dg % dg_rho0 % gcc )

  dg % dg_rho0 % type = ewald_param % ewald_type
  dg % dg_rho0 % zet ( 1 ) = ewald_param % alpha

  CALL dg_rho0_setup ( dg % dg_rho0, pw_grid )

END SUBROUTINE ewald_setup

!!*****
!******************************************************************************

END MODULE ewalds

!******************************************************************************
