!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pws [1.0] *
!!
!!   NAME
!!     pws
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (4-Jan-2001) : fft_wrap moved to coefficient_types and pw_types
!!
!!   SOURCE
!******************************************************************************

MODULE pws

  USE coefficient_types, ONLY : coeff_type, coeff_allocate, coeff_deallocate, &
       SQUARE, SQUAREROOT, PLANEWAVES
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : fourpi
  USE message_passing, ONLY : mp_sum
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE, FULLSPACE
  USE pw_types, ONLY : pw_type, pw_allocate, pw_deallocate, &
       pw_gather, pw_scatter, pw_fft_wrap, &
       REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D, &
       NOSPACE, REALSPACE, RECIPROCALSPACE
  USE termination, ONLY : stop_memory, stop_program
  USE util, ONLY : dotprod_3d

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_poisson_solver, pw_kinetic_evaluate
  PUBLIC :: pw_integral_a, pw_integral_ab, pw_h_psi

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** pws/pw_poisson_solver [1.0] *
!!
!!   NAME
!!     pw_poisson_solver
!!
!!   FUNCTION
!!     Solve Poisson equation in a plane wave basis set
!!     Obtains electrostatic potential and its derivatives with respect to r
!!     from the density
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_poisson_solver ( density, vhartree, dvhartree, h_stress )

  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ), INTENT ( IN ), TARGET :: density
  TYPE ( coeff_type ), INTENT ( OUT ), OPTIONAL :: vhartree
  TYPE ( coeff_type ), DIMENSION ( 3 ), INTENT ( OUT ), OPTIONAL :: dvhartree
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( OUT ), OPTIONAL :: h_stress

! Locals
  INTEGER :: gpt, dir
  REAL ( dbl ) :: gf
  COMPLEX ( dbl ) :: i_cmplx
  INTEGER :: lp, mp, np, ln, mn, nn
  INTEGER, DIMENSION ( : ), POINTER :: ghat

  INTEGER, DIMENSION ( :, : ), POINTER :: g_hat
  TYPE ( pw_grid_type ), POINTER :: pw_grid
  TYPE ( pw_type ) :: rho_g, cfft

!------------------------------------------------------------------------------

  i_cmplx = ( 0._dbl, 1._dbl )

  pw_grid => density % pw % pw_grid
  g_hat => pw_grid % g_hat

  CALL pw_allocate ( rho_g, pw_grid, use_data = COMPLEXDATA3D )

  CALL pw_allocate ( cfft, pw_grid, use_data = COMPLEXDATA3D )
  cfft % in_space = RECIPROCALSPACE

!..FWFFT ( real to g space )
  CALL pw_fft_wrap ( density%pw, rho_g )

  IF ( PRESENT ( vhartree ) ) THEN
     cfft % cc3d = 0._dbl

     DO gpt = 1, pw_grid % ngpts_cut

        ! compute Green's function  '''
        CALL pw_greens_function ( gf, pw_grid % g ( :, gpt ) )

        ghat => g_hat ( :, gpt )

        lp = pw_grid % mapl % pos ( ghat ( 1 ) )
        ln = pw_grid % mapl % neg ( ghat ( 1 ) )
        mp = pw_grid % mapm % pos ( ghat ( 2 ) )
        mn = pw_grid % mapm % neg ( ghat ( 2 ) )
        np = pw_grid % mapn % pos ( ghat ( 3 ) )
        nn = pw_grid % mapn % neg ( ghat ( 3 ) )

        lp = lp + LBOUND ( cfft % cc3d, 1 )
        mp = mp + LBOUND ( cfft % cc3d, 2 )
        np = np + LBOUND ( cfft % cc3d, 3 )
        ln = ln + LBOUND ( cfft % cc3d, 1 )
        mn = mn + LBOUND ( cfft % cc3d, 2 )
        nn = nn + LBOUND ( cfft % cc3d, 3 )

        cfft % cc3d ( lp, mp, np ) = gf * rho_g % cc3d ( lp, mp, np )
        cfft % cc3d ( ln, mn, nn ) = gf * rho_g % cc3d ( ln, mn, nn )
     END DO

! get phi_r
!..BWFFT ( g to real space )
     CALL pw_fft_wrap ( cfft, vhartree%pw )
  END IF

! multiply by the Green's function   '''

  IF ( PRESENT ( dvhartree ) ) THEN
     DO dir = 1, 3
        cfft % cc3d = 0.0_dbl

        DO gpt = 1, pw_grid % ngpts_cut

           ! compute Green's function   '''
           CALL pw_greens_function ( gf, pw_grid % g ( :, gpt ) )

           ghat => g_hat ( :, gpt )

           lp = pw_grid % mapl % pos ( ghat ( 1 ) )
           ln = pw_grid % mapl % neg ( ghat ( 1 ) )
           mp = pw_grid % mapm % pos ( ghat ( 2 ) )
           mn = pw_grid % mapm % neg ( ghat ( 2 ) )
           np = pw_grid % mapn % pos ( ghat ( 3 ) )
           nn = pw_grid % mapn % neg ( ghat ( 3 ) )

           lp = lp + LBOUND ( cfft % cc3d, 1 )
           mp = mp + LBOUND ( cfft % cc3d, 2 )
           np = np + LBOUND ( cfft % cc3d, 3 )
           ln = ln + LBOUND ( cfft % cc3d, 1 )
           mn = mn + LBOUND ( cfft % cc3d, 2 )
           nn = nn + LBOUND ( cfft % cc3d, 3 )

           cfft % cc3d ( lp, mp, np ) = &
                  i_cmplx * pw_grid % g ( dir, gpt ) * gf &
                * rho_g % cc3d ( lp, mp, np )
           cfft % cc3d ( ln, mn, nn ) = &
                - i_cmplx * pw_grid % g ( dir, gpt ) * gf &
                * rho_g % cc3d ( ln, mn, nn )
        END DO

! get derivative of phi_r
!..BWFFT ( g to real space )
        CALL pw_fft_wrap ( cfft, dvhartree ( dir ) %pw)
     END DO
  END IF

  CALL pw_deallocate ( cfft )

  CALL pw_deallocate ( rho_g )

END SUBROUTINE pw_poisson_solver

!!*****
!******************************************************************************

SUBROUTINE pw_greens_function ( gf, g )
  IMPLICIT NONE

! Arguments
  REAL ( dbl ), DIMENSION ( : ), INTENT ( IN ) :: g
  REAL ( dbl ), INTENT ( OUT ) :: gf

! Locals
  INTEGER :: dimensionality  !TMPTMPTMP
  REAL ( dbl ) :: g2

!------------------------------------------------------------------------------

  dimensionality = SIZE ( g )

  SELECT CASE ( dimensionality )
    CASE ( 3 )                   ! 3-dim
       g2 = dotprod_3d ( g, g )
       IF ( g2 > 1.0E-10_dbl ) THEN
          gf = fourpi / g2
       ELSE
          gf = 0.0_dbl
       END IF
    CASE ( 2 )                   ! 2-dim
       stop "greens_function, no 2-dim function implemented yet"
    CASE ( 1 )                   ! 1-dim
       stop "greens_function, no 1-dim function implemented yet"
    CASE ( 0 )                   ! 0-dim = cluster
       stop "greens_function, no 0-dim function implemented yet"
  END SELECT

END SUBROUTINE pw_greens_function

!******************************************************************************
!!****** pws/pw_integral [1.0] *
!!
!!   NAME
!!     pw_integral
!!
!!   FUNCTION
!!     Calculate integral over unit cell for functions in plane wave basis
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

FUNCTION pw_integral_ab ( pw1, pw2, flag ) RESULT ( integral_value )

  TYPE ( pw_type ), INTENT ( IN ) :: pw1, pw2
  INTEGER, INTENT ( IN ), OPTIONAL :: flag

  REAL ( dbl ) :: integral_value

!------------------------------------------------------------------------------

  IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
     integral_value = SUM ( pw1 % cr3d ( :, :, : ) &
          * pw2 % cr3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == REALDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     integral_value = SUM ( pw1 % cr3d ( :, :, : ) &
          * pw2 % cc3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == REALDATA3D ) THEN
     integral_value = SUM ( pw1 % cc3d ( :, :, : ) &
          * pw2 % cr3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     integral_value = SUM ( CONJG ( pw1 % cc3d ( :, :, : ) ) &
          * pw2 % cc3d ( :, :, : ) )

  ELSE IF ( pw1 % in_use == REALDATA1D &
       .AND. pw2 % in_use == REALDATA1D ) THEN
     integral_value = SUM ( pw1 % cr ( : ) * pw2 % cr ( : ) )
  ELSE IF ( pw1 % in_use == REALDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = SUM ( pw1 % cr ( : ) * pw2 % cc ( : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == REALDATA1D ) THEN
     integral_value = SUM ( pw1 % cc ( : ) * pw2 % cr ( : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = SUM ( CONJG ( pw1 % cc ( : ) ) * pw2 % cc ( : ) )
  ELSE
     CALL stop_program ( "integral_ab", "no possible data" )
  END IF

  IF ( pw1 % in_use == REALDATA3D .OR. pw1 % in_use == COMPLEXDATA3D ) THEN
     integral_value = integral_value * pw1 % pw_grid % dvol
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
     integral_value = integral_value
  ELSE
     CALL stop_program ( "pw_integral_ab", "volume element not implemented" )
  END IF

!tmptmp  IF ( pw1 % group_id >= 0 ) CALL mp_sum ( integral_value, pw1 % group_id )

END FUNCTION pw_integral_ab

!******************************************************************************

FUNCTION pw_integral_a ( pw1, flag ) RESULT ( integral_value )

  TYPE ( pw_type ), INTENT ( IN ) :: pw1
  INTEGER, INTENT ( IN ), OPTIONAL :: flag

  REAL ( dbl ) :: integral_value

!------------------------------------------------------------------------------

  IF ( PRESENT ( flag ) ) THEN
     IF ( flag == SQUARE ) THEN
        IF ( pw1 % in_use == REALDATA3D ) THEN
           integral_value = SUM ( pw1 % cr3d ( :, :, : ) ** 2 )
        ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
           integral_value = SUM ( CONJG ( pw1 % cc3d ( :, :, : ) ) &
                * pw1 % cc3d ( :, :, : ) )
        ELSE IF ( pw1 % in_use == REALDATA1D ) THEN
           integral_value = SUM ( pw1 % cr ( : ) ** 2 )
        ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
           integral_value = SUM ( CONJG ( pw1 % cc ( : ) ) &
                * pw1 % cc ( : ) )
        ELSE
           CALL stop_program ( "integral_aa", "no possible SQUARE data" )
        END IF

     ELSE IF ( flag == SQUAREROOT ) THEN
        CALL stop_program ( "integral_aa", "no SQUAREROOT defined" )
     END IF

  ELSE
     IF ( pw1 % in_use == REALDATA3D ) THEN
        integral_value = SUM ( pw1 % cr3d ( :, :, : ) )
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
        integral_value = SUM ( pw1 % cc3d ( :, :, : ) )
     ELSE IF ( pw1 % in_use == REALDATA1D ) THEN
        integral_value = SUM ( pw1 % cr ( : ) )
     ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
        integral_value = SUM ( pw1 % cc ( : ) )
     ELSE
        CALL stop_program ( "integral_aa", "no possible data" )
     END IF
  END IF

  IF ( pw1 % in_use == REALDATA3D ) &
       integral_value = integral_value * pw1 % pw_grid % dvol

!tmptmp  IF ( pw1 % group_id >= 0 ) CALL mp_sum ( integral_value, pw1 % group_id )

END FUNCTION pw_integral_a

!!*****
!******************************************************************************
!!****** pws/pw_kinetic_evaluate [1.0] *
!!
!!   NAME
!!     pw_kinetic_evaluate
!!
!!   FUNCTION
!!     Calculate kinetic energy, potential, derivative of potential and
!!     stress tensor
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_kinetic_evaluate ( psi, xkpt, e_kinetic, kinetic_oper_psi, &
     stress_kinetic )
  
  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: psi
  REAL ( dbl ), INTENT ( IN ) :: xkpt ( 3 )
  
  REAL ( dbl ), INTENT ( OUT ), OPTIONAL :: e_kinetic
  TYPE ( coeff_type ), INTENT ( INOUT ), OPTIONAL :: kinetic_oper_psi
  REAL ( dbl ), INTENT ( INOUT ), OPTIONAL :: stress_kinetic ( 3, 3 )
  
! Local
  INTEGER :: gpt

!------------------------------------------------------------------------------
  
  IF ( .NOT. PRESENT ( kinetic_oper_psi ) ) THEN
     CALL stop_program ( "pw_kinetic_evaluate", &
          "currently kinetic_oper_psi must be present" )
  END IF
  
  IF ( psi % pw % in_use /= COMPLEXDATA1D .OR. &
       kinetic_oper_psi % pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_program ( "pw_kinetic_evaluate", "wrong data in use" )
  END IF

  !OpenMP: parallel do
  DO gpt = 1, psi % pw % pw_grid % ngpts_cut
     kinetic_oper_psi % pw % cc ( gpt ) = kinetic_oper_psi % pw % cc ( gpt ) &
          + 0.5_dbl &
          * SUM ( ( psi % pw % pw_grid % g ( :, gpt ) + xkpt ( : ) ) ** 2 ) &
          * psi % pw % cc ( gpt )
  END DO
  
  IF ( PRESENT ( e_kinetic ) ) &
       e_kinetic = pw_integral_ab ( psi % pw, kinetic_oper_psi % pw )
  
  IF ( PRESENT ( stress_kinetic ) ) THEN
     CALL stop_program ( "pw_kinetic_evaluate", "stress not yet implemented" )
  END IF
  
END SUBROUTINE pw_kinetic_evaluate

!!*****
!******************************************************************************
!!****** pws/pw_h_psi [1.0] *
!!
!!   NAME
!!     pw_h_psi
!!
!!   FUNCTION
!!     Apply local potential to a wavefunction in plane wave basis
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_h_psi ( psi_iks, xkpt, vloc, h_psi, psi_iks_r )

  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: psi_iks
  TYPE ( coeff_type ), INTENT ( IN ), OPTIONAL :: psi_iks_r
  TYPE ( coeff_type ), INTENT ( IN ) :: vloc
  REAL ( dbl ), INTENT ( IN ) :: xkpt ( 3 )

  TYPE ( coeff_type ), INTENT ( INOUT ) :: h_psi

! Locals
  INTEGER :: gpt, gpt2, i1, i2, i3, nhits !TMPTMPTMP
  REAL ( dbl ) :: ekinetic !TMPTMPTMP
!  COMPLEX ( dbl ) :: foo !TMPTMPTMP
  TYPE ( pw_type ) :: cfft

!------------------------------------------------------------------------------

  CALL pw_allocate ( cfft, psi_iks % pw % pw_grid, &
       use_data = COMPLEXDATA3D )

#if defined ( DEBUG )
  nhits = 0
  OPEN ( UNIT = 99, FILE = "psik", FORM = "FORMATTED" )
  DO gpt = 1, psi_iks % grid % ngpts
     READ ( UNIT = 99, FMT = * ) i1, i2, i3, foo
     DO gpt2 = 1, psi_iks % grid % ngpts
        IF ( psi_iks % grid % g_hat ( 1, gpt2 ) == i1 .AND. &
             psi_iks % grid % g_hat ( 2, gpt2 ) == i2 .AND. &
             psi_iks % grid % g_hat ( 3, gpt2 ) == i3 ) THEN
           psi_iks % cc ( gpt2 ) = foo
           nhits = nhits + 1
           EXIT
        END IF
     END DO
  END DO
  CLOSE ( UNIT = 99 )
  WRITE ( 6, * ) "DDD: #hits = ", nhits
#endif

  IF ( PRESENT ( psi_iks_r ) ) THEN
     cfft % cr3d = vloc % pw % cr3d * psi_iks_r % pw % cc3d
  ELSE
     CALL pw_fft_wrap ( psi_iks%pw, cfft )

     cfft % cc3d = cfft % cc3d * vloc % pw % cr3d
  END IF

  CALL pw_fft_wrap ( cfft, h_psi%pw )

  !! This must be done after the local potential
  CALL pw_kinetic_evaluate ( psi_iks, xkpt, ekinetic, h_psi )

  CALL pw_deallocate ( cfft )

END SUBROUTINE pw_h_psi

!******************************************************************************

SUBROUTINE pw_stateset_density ( psi_iks, n_s )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), DIMENSION ( : ), INTENT ( IN ) :: psi_iks
  
  TYPE ( coeff_type ), INTENT ( INOUT ) :: n_s
  
! Locals
  INTEGER :: state, nstates
  REAL ( dbl ) :: wkptf
  TYPE ( pw_type ) :: cfft
  
!------------------------------------------------------------------------------
  
  CALL pw_allocate ( cfft, psi_iks ( 1 ) % pw % pw_grid, &
       use_data = COMPLEXDATA3D )
  
  nstates = SIZE ( psi_iks )
  
  !OpenMP parallel do ?
  DO state = 1, nstates
     CALL pw_fft_wrap ( psi_iks ( state )%pw , cfft )
     
     wkptf = 1.0_dbl / psi_iks ( state ) % pw % pw_grid % vol ! TMPTMPTMP
     
     n_s % pw % cr3d = n_s % pw % cr3d + wkptf &
          * ( REAL ( cfft % cc3d ) ** 2 + AIMAG ( cfft % cc3d ) ** 2 )
  END DO
  
  CALL pw_deallocate ( cfft )
  
END SUBROUTINE pw_stateset_density

!!*****
!******************************************************************************

END MODULE pws

!******************************************************************************
