!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/dg_rho0s [1.0] *
!!
!!   NAME
!!     dg_rho0s
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE dg_rho0s

  USE coefficient_types,               ONLY: PW_COMPLEXDATA3D,&
                                             PW_REALDATA3D,&
                                             coeff_allocate,&
                                             coeff_type,&
                                             coeff_zero
  USE kinds,                           ONLY: dbl
  USE mathconstants,                   ONLY: pi
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC:: dg_rho0_type, dg_rho0_setup

  ! type = type of gaussian (PME)
  ! grid = grid number
  ! gcc = Gaussian contraction coefficient
  ! zet = Gaussian exponent
  TYPE dg_rho0_type
     CHARACTER ( LEN = 12 ) :: type
     INTEGER :: grid
     INTEGER :: kind
     REAL ( dbl ) :: cutoff_radius
     REAL ( dbl ), DIMENSION ( : ), POINTER :: gcc
     REAL ( dbl ), DIMENSION ( : ), POINTER :: zet
     TYPE ( coeff_type ) :: density
  END TYPE dg_rho0_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE dg_rho0_setup ( dg_rho0, pw_grid )


    TYPE(dg_rho0_type), INTENT(INOUT)        :: dg_rho0
    TYPE(pw_grid_type), INTENT(IN)           :: pw_grid

!------------------------------------------------------------------------------

  SELECT CASE ( dg_rho0 % type )
  CASE ( "EWALD" )
    CALL coeff_allocate ( dg_rho0 % density, pw_grid, PW_REALDATA3D )
    CALL dg_rho0_pme_gauss ( dg_rho0 % density, dg_rho0 % zet ( 1 ) )
  CASE ( "PME" )
    CALL coeff_allocate ( dg_rho0 % density, pw_grid, PW_REALDATA3D )
    CALL dg_rho0_pme_gauss ( dg_rho0 % density, dg_rho0 % zet ( 1 ) )
  CASE ( "SPME" )
    CALL stop_program ( 'dg_rho0_setup', 'type not implemented' )
  CASE ( "S" )
    CALL coeff_allocate ( dg_rho0 % density, pw_grid, PW_REALDATA3D )
     CALL dg_rho0_s_cgf ( dg_rho0 % density, dg_rho0 % zet , &
                          dg_rho0 % gcc)
  CASE ( "PX", "PY", "PZ" )
    CALL coeff_allocate ( dg_rho0 % density, pw_grid, PW_COMPLEXDATA3D )
    CALL dg_rho0_p_cgf ( dg_rho0 % density, dg_rho0 % zet, & 
                         dg_rho0 % gcc, dg_rho0 % type)
  CASE DEFAULT
    CALL stop_program ( 'dg_rho0_setup', 'type not implemented' )
  END SELECT

END SUBROUTINE dg_rho0_setup

!******************************************************************************

SUBROUTINE dg_rho0_pme_gauss ( dg_rho0, alpha )


    TYPE(coeff_type), INTENT(INOUT), TARGET  :: dg_rho0
    REAL(dbl), INTENT(IN)                    :: alpha

    INTEGER, PARAMETER                       :: IMPOSSIBLE = 10000

    INTEGER                                  :: gpt, l0, ln, lp, m0, mn, mp, &
                                                n0, nn, np
    INTEGER, DIMENSION(:), POINTER           :: ghat
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(dbl)                                :: const, e_gsq
    REAL(dbl), DIMENSION(:, :, :), POINTER   :: rho0
    TYPE(pw_grid_type), POINTER              :: pw_grid

!------------------------------------------------------------------------------

  const = 1.0_dbl / ( 8.0_dbl * alpha ** 2 )

  pw_grid => dg_rho0 % pw % pw_grid
  bds => pw_grid % bounds

  IF ( -bds ( 1, 1 ) == bds ( 2, 1 ) ) THEN
    l0 = IMPOSSIBLE
  ELSE
    l0 = bds ( 1, 1 )
  END IF

  IF ( -bds ( 1, 2 ) == bds ( 2, 2 ) ) THEN
    m0 = IMPOSSIBLE
  ELSE
    m0 = bds ( 1, 2 )
  END IF

  IF ( -bds ( 1, 3 ) == bds ( 2, 3 ) ) THEN
    n0 = IMPOSSIBLE
  ELSE
    n0 = bds ( 1, 3 )
  END IF

  CALL coeff_zero ( dg_rho0 )

  rho0 => dg_rho0 % pw % cr3d

  DO gpt = 1, pw_grid % ngpts_cut
    ghat => pw_grid % g_hat ( :, gpt )

    lp = pw_grid % mapl % pos ( ghat ( 1 ) )
    ln = pw_grid % mapl % neg ( ghat ( 1 ) )
    mp = pw_grid % mapm % pos ( ghat ( 2 ) )
    mn = pw_grid % mapm % neg ( ghat ( 2 ) )
    np = pw_grid % mapn % pos ( ghat ( 3 ) )
    nn = pw_grid % mapn % neg ( ghat ( 3 ) )

    e_gsq = EXP ( -const * pw_grid % gsq ( gpt ) ) / pw_grid % vol

    !*apsi
    lp = lp + bds ( 1, 1 )
    mp = mp + bds ( 1, 2 )
    np = np + bds ( 1, 3 )
    ln = ln + bds ( 1, 1 )
    mn = mn + bds ( 1, 2 )
    nn = nn + bds ( 1, 3 )

    rho0 ( lp, mp, np ) = e_gsq
    rho0 ( ln, mn, nn ) = e_gsq

    IF ( ghat ( 1 ) == l0 .OR. ghat ( 2 ) == m0 .OR. ghat ( 3 ) == n0 ) THEN
      rho0 ( lp, mp, np ) = 0.0_dbl
      rho0 ( ln, mn, nn ) = 0.0_dbl
    END IF

  END DO

END SUBROUTINE dg_rho0_pme_gauss
!******************************************************************************

SUBROUTINE dg_rho0_s_cgf ( dg_rho0, zet, gcc )

!******************************************************************************
!   NAME
!!     dg_rho0_s_cgf
!!
!!   FUNCTION
!!     puts a contracted gaussian function of s-type                
!!     in the g-space                            
!!
!!   AUTHOR
!!     gloria      
!!
!!   SOURCE
!******************************************************************************




    TYPE(coeff_type), INTENT(INOUT), TARGET  :: dg_rho0
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: zet, gcc

    INTEGER, PARAMETER                       :: IMPOSSIBLE = 10000

    INTEGER                                  :: gpt, icg, l0, ln, lp, m0, mn, &
                                                mp, n0, ncg, nn, np
    INTEGER, DIMENSION(:), POINTER           :: ghat
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(dbl)                                :: charge, coef, const, e_gsq
    REAL(dbl), DIMENSION(:, :, :), POINTER   :: rho0
    TYPE(pw_grid_type), POINTER              :: pw_grid

!------------------------------------------------------------------------------

  ncg = SIZE ( zet  )

  pw_grid => dg_rho0 % pw % pw_grid
  bds => pw_grid % bounds

  IF ( -bds ( 1, 1 ) == bds ( 2, 1 ) ) THEN
    l0 = IMPOSSIBLE
  ELSE
    l0 = bds ( 1, 1 )
  END IF

  IF ( -bds ( 1, 2 ) == bds ( 2, 2 ) ) THEN
    m0 = IMPOSSIBLE
  ELSE
    m0 = bds ( 1, 2 )
  END IF

  IF ( -bds ( 1, 3 ) == bds ( 2, 3 ) ) THEN
    n0 = IMPOSSIBLE
  ELSE
    n0 = bds ( 1, 3 )
  END IF

  CALL coeff_zero ( dg_rho0 )

  rho0 => dg_rho0 % pw % cr3d

  charge = 0._dbl
  DO icg = 1, ncg

     const = 1.0_dbl / ( 4.0_dbl * zet ( icg ) )
     coef  = ( gcc ( icg ) * pi ** 1.5_dbl ) / zet ( icg )**1.5_dbl

     rho0(l0,m0,n0)  =  rho0(l0,m0,n0)  + coef * 1._dbl / pw_grid % vol

     DO gpt = 1, pw_grid % ngpts_cut
        ghat => pw_grid % g_hat ( :, gpt )

        lp = pw_grid % mapl % pos ( ghat ( 1 ) )
        ln = pw_grid % mapl % neg ( ghat ( 1 ) )
        mp = pw_grid % mapm % pos ( ghat ( 2 ) )
        mn = pw_grid % mapm % neg ( ghat ( 2 ) )
        np = pw_grid % mapn % pos ( ghat ( 3 ) )
        nn = pw_grid % mapn % neg ( ghat ( 3 ) )

        e_gsq = EXP ( -const * pw_grid % gsq ( gpt ) ) / pw_grid % vol

    !*apsi
        lp = lp + bds ( 1, 1 )
        mp = mp + bds ( 1, 2 )
        np = np + bds ( 1, 3 )
        ln = ln + bds ( 1, 1 )
        mn = mn + bds ( 1, 2 )
        nn = nn + bds ( 1, 3 )

        IF ( lp == l0 .AND. mp == m0 .AND. np == n0 ) CYCLE
        rho0 ( lp, mp, np ) = rho0 ( lp, mp, np ) + e_gsq * coef
        rho0 ( ln, mn, nn ) = rho0 ( ln, mn, nn ) + e_gsq * coef

        IF ( ghat ( 1 ) == l0 .OR. ghat ( 2 ) == m0 .OR. ghat ( 3 ) == n0 ) THEN
          rho0 ( lp, mp, np ) = 0.0_dbl
          rho0 ( ln, mn, nn ) = 0.0_dbl
        END IF

      END DO

   END DO

END SUBROUTINE dg_rho0_s_cgf

!******************************************************************************

SUBROUTINE dg_rho0_p_cgf ( dg_rho0, zet, gcc, type )

!******************************************************************************
!   NAME
!!     dg_rho0_p_cgf
!!
!!   FUNCTION
!!     puts a contracted gaussian function of p-type                
!!     in the g-space. The p is imaginary, so
!!     the p-function has then to be multiplied by (0.0, -1.0)
!!     in order to get the true p.
!!
!!   AUTHOR
!!     gloria      
!!
!!   SOURCE
!******************************************************************************




    TYPE(coeff_type), INTENT(INOUT), TARGET  :: dg_rho0
    REAL(dbl), DIMENSION(:), INTENT(IN)      :: zet, gcc
    CHARACTER(LEN=12), INTENT(IN)            :: type

    INTEGER, PARAMETER                       :: IMPOSSIBLE = 10000

    COMPLEX(dbl)                             :: minus_i
    COMPLEX(dbl), DIMENSION(:, :, :), &
      POINTER                                :: rho0
    INTEGER                                  :: gpt, icg, l0, ln, lp, m0, mn, &
                                                mp, n0, ncg, nn, np
    INTEGER, DIMENSION(:), POINTER           :: ghat
    INTEGER, DIMENSION(:, :), POINTER        :: bds
    REAL(dbl)                                :: charge, coef, const, e_gsq, &
                                                p_on_grid
    TYPE(pw_grid_type), POINTER              :: pw_grid

!------------------------------------------------------------------------------

  minus_i = CMPLX( 0._dbl, -1.0_dbl, dbl )
  ncg = size(zet)

  pw_grid => dg_rho0 % pw % pw_grid
  bds => pw_grid % bounds

  IF ( -bds ( 1, 1 ) == bds ( 2, 1 ) ) THEN
    l0 = IMPOSSIBLE
  ELSE
    l0 = bds ( 1, 1 )
  END IF

  IF ( -bds ( 1, 2 ) == bds ( 2, 2 ) ) THEN
    m0 = IMPOSSIBLE
  ELSE
    m0 = bds ( 1, 2 )
  END IF

  IF ( -bds ( 1, 3 ) == bds ( 2, 3 ) ) THEN
    n0 = IMPOSSIBLE
  ELSE
    n0 = bds ( 1, 3 )
  END IF

  CALL coeff_zero ( dg_rho0 )

  rho0 => dg_rho0 % pw % cc3d
  charge = 0._dbl
  DO icg = 1, ncg

     IF ( zet (icg) == 0 ) CYCLE
     
     const = 1.0_dbl / ( 4.0_dbl * zet ( icg ) )
     coef  = ( gcc ( icg ) * pi**1.5 ) / ( 2.0_dbl * zet ( icg )**2.5 )

     DO gpt = 1, pw_grid % ngpts_cut
        ghat => pw_grid % g_hat ( :, gpt )

        lp = pw_grid % mapl % pos ( ghat ( 1 ) )
        ln = pw_grid % mapl % neg ( ghat ( 1 ) )
        mp = pw_grid % mapm % pos ( ghat ( 2 ) )
        mn = pw_grid % mapm % neg ( ghat ( 2 ) )
        np = pw_grid % mapn % pos ( ghat ( 3 ) )
        nn = pw_grid % mapn % neg ( ghat ( 3 ) )

        e_gsq = EXP ( -const * pw_grid % gsq ( gpt ) ) / pw_grid % vol

    !*apsi
        lp = lp + bds ( 1, 1 )
        mp = mp + bds ( 1, 2 )
        np = np + bds ( 1, 3 )
        ln = ln + bds ( 1, 1 )
        mn = mn + bds ( 1, 2 )
        nn = nn + bds ( 1, 3 )

        IF ( lp == l0 .AND. mp == m0 .AND. np == n0 ) CYCLE

        SELECT CASE (type)
 
        CASE ("PX")

           p_on_grid = pw_grid % g (1, gpt) * e_gsq * coef
           
        CASE ("PY")

           p_on_grid = pw_grid % g (2, gpt) * e_gsq * coef

        CASE ("PZ")

           p_on_grid = pw_grid % g (3, gpt) * e_gsq * coef

        END SELECT

        rho0 ( lp, mp, np ) = rho0 ( lp, mp, np ) + p_on_grid * minus_i  
        rho0 ( ln, mn, nn ) = rho0 ( ln, mn, nn ) - p_on_grid * minus_i

        IF ( ghat ( 1 ) == l0 .OR. ghat ( 2 ) == m0 .OR. ghat ( 3 ) == n0 ) THEN
          rho0 ( lp, mp, np ) = 0.0_dbl
          rho0 ( ln, mn, nn ) = 0.0_dbl
        END IF

      END DO

   END DO

END SUBROUTINE dg_rho0_p_cgf

!****************************************************************************

!******************************************************************************

END MODULE dg_rho0s
