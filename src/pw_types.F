!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/pw_types [1.0] *
!!
!!   NAME
!!     pw_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (29-Dec-2000) : Changes for parallel use
!!     JGH (13-Mar-2001) : added timing calls
!!     JGH (26-Feb-2003) : OpenMP enabled
!!     JGH (17-Nov-2007) : Removed mass arrays
!!     JGH (01-Dec-2007) : Removed and renamed routines
!!
!!   NOTES
!!     If parallel mode is distributed certain combination of
!!     "in_use" and "in_space" can not be used.
!!     For performance reasons it would be better to have the loops
!!     over g-vectros in the gather/scatter routines in new subprograms
!!     with the actual arrays (also the adressing) in the parameter list
!!
!!
!!   SOURCE
!******************************************************************************

MODULE pw_types
  USE fft_tools,                       ONLY: BWFFT,&
                                             FWFFT,&
                                             fft3d
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_loc_c,&
                                             m_loc_r
  USE message_passing,                 ONLY: mp_sum
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             PW_MODE_DISTRIBUTED,&
                                             PW_MODE_LOCAL,&
                                             pw_grid_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: print_stack,&
                                             timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_type, pw_p_type
  PUBLIC :: pw_retain, pw_release, pw_create
  PUBLIC :: REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D, NODATA
  PUBLIC :: NOSPACE, REALSPACE, RECIPROCALSPACE
  PUBLIC :: pw_zero, pw_structure_factor, pw_smoothing
  PUBLIC :: pw_copy, pw_axpy, pw_transfer, pw_scale
  PUBLIC :: pw_derive, pw_dr2, pw_fft_wrap, pw_write
  PUBLIC :: pw_compare_debug
  PUBLIC :: pw_integral_aa, pw_integral_ab, pw_integral_a2b
  PUBLIC :: pw_dr2_gg, pw_integrate_function

  TYPE pw_type
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: cr
     REAL (KIND=dp), DIMENSION ( :, :, : ), POINTER :: cr3d
     COMPLEX (KIND=dp), DIMENSION ( : ), POINTER :: cc
     COMPLEX (KIND=dp), DIMENSION ( :, :, : ), POINTER :: cc3d

     INTEGER :: in_use              ! Which data is used [r1d/c1d/r3d/c3d]
     INTEGER :: in_space            ! Real/Reciprocal space
     INTEGER :: id_nr               ! unique identifier
     INTEGER :: ref_count           ! reference count

     TYPE ( pw_grid_type ), POINTER :: pw_grid
  END TYPE pw_type
  
  TYPE pw_p_type
     TYPE(pw_type), POINTER :: pw
  END TYPE pw_p_type

  ! Flags for the structure member 'in_use'
  INTEGER, PARAMETER :: REALDATA1D = 301, COMPLEXDATA1D = 302
  INTEGER, PARAMETER :: REALDATA3D = 303, COMPLEXDATA3D = 304, NODATA = 305

  ! Flags for the structure member 'in_space'
  INTEGER, PARAMETER :: NOSPACE = 371, REALSPACE = 372, RECIPROCALSPACE = 373

  ! to generate unique id_nr
  INTEGER, SAVE, PRIVATE :: last_pw_id_nr=0
  INTEGER, SAVE, PRIVATE :: allocated_pw_count=0

  INTEGER, PUBLIC, PARAMETER :: SQUARE = 391, SQUAREROOT = 392

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_types'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

  INTERFACE pw_gather
     MODULE PROCEDURE pw_gather_s, pw_gather_p
  END INTERFACE

  INTERFACE pw_scatter
     MODULE PROCEDURE pw_scatter_s, pw_scatter_p
  END INTERFACE

  INTERFACE pw_fft_wrap
     MODULE PROCEDURE fft_wrap_pw1, fft_wrap_pw1pw2
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

!!****** pw_types/pw_zero [1.0] *
!!
!!   NAME
!!     pw_zero
!!
!!   SYNOPSIS
!!     Subroutine pw_zero(pw)
!!       Type(pw_type), Intent (INOUT):: pw
!!     End Subroutine pw_zero
!!
!!   FUNCTION
!!     Set values of a pw type to zero
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_zero ( pw )


    TYPE(pw_type), INTENT(INOUT)             :: pw

    INTEGER                                  :: handle, i, ns
    REAL(KIND=dp)                            :: zr

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_zero', 'I', 'Mzero', handle )
  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_zero")
  IF ( pw % in_use == REALDATA1D ) THEN
     ns = SIZE(pw%cr)
     pw%cr(:) = 0._dp
  ELSE IF ( pw % in_use == COMPLEXDATA1D ) THEN
     ns = SIZE(pw%cc)
     pw%cc(:) = CMPLX(0._dp,0._dp,KIND=dp)
  ELSE IF ( pw % in_use == REALDATA3D ) THEN
     ns = SIZE(pw%cr3d)
     pw%cr3d(:,:,:) = 0._dp
  ELSE IF ( pw % in_use == COMPLEXDATA3D ) THEN
     ns = SIZE(pw%cc3d)
     pw%cc3d(:,:,:) = CMPLX(0._dp,0._dp,KIND=dp)
  ELSE
     CALL stop_program ( "pw_zero", "no possible data field" )
  END IF

  zr = REAL ( ns,KIND=dp) * 1.e-6_dp
  CALL timestop ( zr, handle )

END SUBROUTINE pw_zero

!!*****
!******************************************************************************
!!****** pw_types/pw_copy [1.0] *
!!
!!   NAME
!!     pw_copy
!!
!!   SYNOPSIS
!!     Subroutine pw_copy(pw1, pw2)
!!       Type(pw_type), Intent (IN):: pw1
!!       Type(pw_type), Intent (INOUT):: pw2
!!     End Subroutine pw_copy
!!
!!   FUNCTION
!!     copy a pw type variable
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (7-Mar-2001) : check for pw_grid %id_nr, allow copy if 
!!       in_use == COMPLEXDATA1D and in_space == RECIPROCALSPACE
!!     JGH (21-Feb-2003) : Code for generalized reference grids
!!
!!   NOTES
!!     Currently only copying of respective types allowed, 
!!     in order to avoid errors
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_copy ( pw1, pw2 )
  
  
    TYPE(pw_type), INTENT(IN)                :: pw1
    TYPE(pw_type), INTENT(INOUT)             :: pw2

    INTEGER                                  :: handle, i, j, ng, ng1, ng2, ns
    REAL(KIND=dp)                            :: zc

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_copy', 'I', 'Mcopy', handle )
  CPPreconditionNoErr(pw1%ref_count>0,cp_failure_level,"pw_copy")
  CPPreconditionNoErr(pw2%ref_count>0,cp_failure_level,"pw_copy")
  IF ( pw1 % pw_grid %id_nr /= pw2 % pw_grid %id_nr ) THEN

    IF ( pw1 % pw_grid % spherical .AND. pw2 % pw_grid % spherical ) THEN

      IF ( pw_compatible ( pw1 % pw_grid, pw2 % pw_grid ) ) THEN

        IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
             pw2 % in_use == COMPLEXDATA1D .AND. &
             pw1 % in_space == RECIPROCALSPACE ) THEN
           ng1 = SIZE ( pw1 % cc )
           ng2 = SIZE ( pw2 % cc )
           ng = MIN ( ng1, ng2 )
!$omp parallel do private(i)
           DO i = 1, ng
              pw2%cc(i) = pw1%cc(i)
           END DO
           IF ( ng2 > ng ) THEN
!$omp parallel do private(i)
              DO i = ng+1, ng2
                 pw2%cc(i) = CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
              END DO
           END IF
           ns = 2 * MAX ( ng1, ng2 )
        ELSE
           CALL stop_program ( "pw_copy", "no suitable data field" )
        END IF

      ELSE
        WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                 " grid 1 :",pw1 % pw_grid %id_nr, &
                 " sperical :",pw1 % pw_grid % spherical, &
                 " reference :",pw1 % pw_grid % reference 
        WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                 " grid 2 :",pw2 % pw_grid %id_nr, &
                 " sperical :",pw2 % pw_grid % spherical, &
                 " reference :",pw2 % pw_grid % reference 
        CALL stop_program ( "pw_copy", "Incompatible grids" )
      END IF

    ELSE IF ( .NOT. ( pw1 % pw_grid % spherical .OR. &
                      pw2 % pw_grid % spherical ) ) THEN

      ng1 = SIZE ( pw1 % cc )
      ng2 = SIZE ( pw2 % cc )
      ns = 2 * MAX ( ng1, ng2 )

      IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
           pw2 % in_use == COMPLEXDATA1D .AND. &
           pw1 % in_space == RECIPROCALSPACE ) THEN

        IF ( ( pw1 % pw_grid %id_nr == pw2 % pw_grid % reference ) ) THEN
          IF( ng1 >= ng2 ) THEN
!$omp parallel do private(i,j)
            DO i = 1, ng2
              j = pw2 % pw_grid % gidx ( i )
              pw2 % cc ( i ) = pw1 % cc ( j )
            END DO
          ELSE
            CALL pw_zero ( pw2 )
!$omp parallel do private(i,j)
            DO i = 1, ng1
              j = pw2 % pw_grid % gidx ( i )
              pw2 % cc ( j ) = pw1 % cc ( i )
            END DO
          END IF
        ELSE IF ( ( pw2 % pw_grid %id_nr == pw1 % pw_grid % reference ) ) THEN
          IF( ng1 >= ng2 ) THEN
!$omp parallel do private(i,j)
            DO i = 1, ng2
              j = pw1 % pw_grid % gidx ( i )
              pw2 % cc ( i ) = pw1 % cc ( j )
            END DO
          ELSE
            CALL pw_zero ( pw2 )
!$omp parallel do private(i,j)
            DO i = 1, ng1
              j = pw1 % pw_grid % gidx ( i )
              pw2 % cc ( j ) = pw1 % cc ( i )
            END DO
          END IF
        ELSE
          WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                 " grid 1 :",pw1 % pw_grid %id_nr, &
                 " sperical :",pw1 % pw_grid % spherical, &
                 " reference :",pw1 % pw_grid % reference 
          WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                 " grid 2 :",pw2 % pw_grid %id_nr, &
                 " sperical :",pw2 % pw_grid % spherical, &
                 " reference :",pw2 % pw_grid % reference 
          CALL stop_program ( "pw_copy", "Incompatible grids" )
        END IF

      ELSE
         CALL stop_program ( "pw_copy", "no suitable data field" )
      END IF

      pw2 % in_space = RECIPROCALSPACE

    ELSE
      WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 1 :",pw1 % pw_grid %id_nr, &
               " sperical :",pw1 % pw_grid % spherical, &
               " reference :",pw1 % pw_grid % reference 
      WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 2 :",pw2 % pw_grid %id_nr, &
               " sperical :",pw2 % pw_grid % spherical, &
               " reference :",pw2 % pw_grid % reference 
      CALL stop_program ( "pw_copy", "Incompatible grids" )
    END IF

  ELSE

    IF ( pw1 % in_use == REALDATA1D .AND. pw2 % in_use == REALDATA1D ) THEN
       ns = SIZE ( pw1 % cr )
!$omp parallel do private(i)
       DO i = 1, ns
          pw2 % cr(i) = pw1 % cr(i)
       END DO
    ELSE IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
         pw2 % in_use == COMPLEXDATA1D ) THEN
       ns = SIZE ( pw1 % cc )
!$omp parallel do private(i)
       DO i = 1, ns
          pw2 % cc(i) = pw1 % cc(i)
       END DO
    ELSE IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
       ns = SIZE ( pw1 % cr3d )
       pw2 % cr3d(:,:,:) = pw1 % cr3d(:,:,:)
    ELSE IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
         pw2 % in_use == COMPLEXDATA3D ) THEN
       ns = SIZE ( pw1 % cc3d )
       pw2 % cc3d(:,:,:) = pw1 % cc3d(:,:,:)
    ELSE
       CALL stop_program ( "pw_copy", "no suitable data field" )
    END IF

  END IF

  pw2 % in_space = pw1 % in_space

  zc = REAL ( ns,KIND=dp) * 1.e-6_dp
  CALL timestop ( zc, handle )
  
END SUBROUTINE pw_copy

!!****f* pw_types/pw_scale *
!!
!!   NAME
!!     pw_scale
!!   FUNCTION
!!     multiplies pw coeffs with a number
!!   NOTES
!!
!!   INPUTS
!!    -
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
!!*****
SUBROUTINE pw_scale( pw, a )

    TYPE(pw_type), INTENT(INOUT)             :: pw
    REAL(KIND=dp), INTENT(IN)                :: a

    INTEGER                                  :: handle, ns
    REAL(KIND=dp)                            :: flop

    CALL timeset ( 'pw_scale', 'I', 'Mflops', handle )
    CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_scale")

    SELECT CASE ( pw % in_use )
    CASE ( REALDATA1D )
       ns = SIZE(pw%cr)
       pw%cr(:) = a*pw%cr(:)
    CASE ( COMPLEXDATA1D )
       ns = 2*SIZE(pw%cc)
       pw%cc(:) = a*pw%cc(:)
    CASE ( REALDATA3D)
       ns = SIZE(pw%cr3d)
       pw%cr3d(:,:,:) = a*pw%cr3d(:,:,:)
    CASE (COMPLEXDATA3D )
       ns = 2*SIZE(pw%cc3d)
       pw%cc3d(:,:,:) = a*pw%cc3d(:,:,:)
    CASE DEFAULT  
       CALL stop_program ( "pw_scale", "no suitable data field" )
    END SELECT

    flop = REAL ( ns, KIND=dp) * 1.e-6_dp
    CALL timestop ( flop, handle )

END SUBROUTINE pw_scale

!******************************************************************************
!!****** pw_types/pw_derive [1.0] *
!!
!!   NAME
!!     pw_derive
!!
!!   SYNOPSIS
!!     Subroutine pw_derive(pw, n)
!!       Type(pw_type), Intent (INOUT):: pw
!!       Integer, Dimension(3), Intent (IN):: n
!!     End Subroutine pw_derive
!!
!!   FUNCTION
!!     Calculate the derivative of a plane wave vector
!!
!!   AUTHOR
!!     JGH (25-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (06-10-2002) allow only for inplace derivatives
!!
!!   NOTES
!!     Calculate the derivative dx^n(1) dy^n(2) dz^n(3) PW
!!     PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_derive ( pw, n )

    TYPE(pw_type), INTENT(INOUT)             :: pw
    INTEGER, DIMENSION(3), INTENT(IN)        :: n

    COMPLEX(KIND=dp)                         :: im
    INTEGER                                  :: cnt, handle, i, m
    REAL(KIND=dp)                            :: flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_derive', 'I', 'Mflops', handle )

  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_derive")
  IF ( ANY ( n < 0 ) ) CALL stop_program ( "pw_derive", "n < 0" )

  m = SUM ( n )
  im = CMPLX ( 0.0_dp, 1.0_dp,KIND=dp) ** m

  flop = 0.0_dp

  IF ( pw % in_space == RECIPROCALSPACE .AND. &
       pw % in_use == COMPLEXDATA1D ) THEN

    cnt = SIZE ( pw % cc )

    IF ( n ( 1 ) == 1 ) THEN
!$omp parallel do private (i)
      DO i = 1, cnt
         pw % cc ( i ) = pw % cc ( i ) * pw % pw_grid % g ( 1, i )
      END DO
      flop = flop + 6 * cnt
    ELSE IF ( n ( 1 ) > 1 ) THEN
!$omp parallel do private (i)
      DO i = 1, cnt
         pw % cc ( i ) = pw % cc ( i ) * ( pw % pw_grid % g ( 1, i ) ** n ( 1 ) )
      END DO
      flop = flop + 7 * cnt
    END IF
    IF ( n ( 2 ) == 1 ) THEN
!$omp parallel do private (i)
      DO i = 1, cnt
         pw % cc ( i ) = pw % cc ( i ) * pw % pw_grid % g ( 2, i )
      END DO
      flop = flop + 6 * cnt
    ELSE IF ( n ( 2 ) > 1 ) THEN
!$omp parallel do private (i)
      DO i = 1, cnt
         pw % cc ( i ) = pw % cc ( i ) * ( pw % pw_grid % g ( 2, i ) ** n ( 2 ) )
      END DO
      flop = flop + 7 * cnt
    END IF
    IF ( n ( 3 ) == 1 ) THEN
!$omp parallel do private (i)
      DO i = 1, cnt
         pw % cc ( i ) = pw % cc ( i ) * pw % pw_grid % g ( 3, i )
      END DO
      flop = flop + 6 * cnt
    ELSE IF ( n ( 3 ) > 1 ) THEN
!$omp parallel do private (i)
      DO i = 1, cnt
         pw % cc ( i ) = pw % cc ( i ) * ( pw % pw_grid % g ( 3, i ) ** n ( 3 ) )
      END DO
      flop = flop + 7 * cnt
    END IF

    ! im can take the values 1, -1, i, -i
    ! skip this if im == 1
    IF ( ABS ( REAL ( im,KIND=dp) - 1.0_dp ) > 1.e-10 ) THEN
!$omp parallel do private (i)
      DO i = 1, cnt
         pw % cc ( i ) = im * pw % cc ( i )
      END DO
      flop = flop + 6 * cnt
    END IF

  ELSE

    CALL stop_program ( "pw_derive", "no suitable data field" )

  END IF

  flop = flop * 1.e-6_dp
  CALL timestop ( flop, handle )

END SUBROUTINE pw_derive

!******************************************************************************
!!****** pw_types/pw_dr2 [1.0] *
!!
!!   NAME
!!     pw_dr2
!!
!!   SYNOPSIS
!!     Subroutine pw_dr2(pw, pwdr2)
!!       Type(pw_type), Intent (IN)               :: pw
!!       Type(pw_type), DIMENSION(:), Intent (OUT):: pwdr2
!!     End Subroutine pw_dr2
!!
!!   FUNCTION
!!     Calculate the tensorial 2nd derivative of a plane wave vector
!!
!!   AUTHOR
!!     JGH (05-May-2006)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
!!
!!*****
!******************************************************************************

SUBROUTINE pw_dr2 ( pw, pwdr2, i, j )

    TYPE(pw_type), INTENT(IN)                :: pw
    TYPE(pw_type), INTENT(INOUT)             :: pwdr2
    INTEGER, INTENT(IN)                      :: i, j

    INTEGER                                  :: cnt, handle, ig
    REAL(KIND=dp)                            :: flop, gg, o3

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_dr2', 'I', 'Mflops', handle )

  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_dr2")

  flop = 0.0_dp
  o3 = 1._dp/3._dp

  IF ( pw % in_space == RECIPROCALSPACE .AND. &
       pw % in_use == COMPLEXDATA1D ) THEN

    cnt = SIZE ( pw % cc )

    IF ( i == j ) THEN
!$omp parallel do private (ig)
      DO ig = 1, cnt
        gg = pw%pw_grid%g(i,ig)**2 - o3*pw%pw_grid%gsq(ig)
        pwdr2%cc(ig) = gg * pw%cc(ig)
      END DO
      flop = flop + 5 * cnt
    ELSE
!$omp parallel do private (ig)
      DO ig = 1, cnt
        pwdr2%cc(ig) = pw%cc(ig) * (pw%pw_grid%g(i,ig)*pw%pw_grid%g(j,ig))
      END DO
      flop = flop + 4 * cnt
    END IF

  ELSE

    CALL stop_program ( "pw_dr2", "no suitable data field" )

  END IF

  flop = flop * 1.e-6_dp
  CALL timestop ( flop, handle )

END SUBROUTINE pw_dr2

!******************************************************************************
!!****** pw_types/pw_dr2_gg [1.0] *
!!
!!   NAME
!!     pw_dr2_gg
!!
!!   SYNOPSIS
!!     Subroutine pw_dr2(pw, pwdr2)
!!       Type(pw_type), Intent (IN)               :: pw
!!       Type(pw_type), DIMENSION(:), Intent (OUT):: pwdr2_gg
!!     End Subroutine pw_dr2
!!
!!   FUNCTION
!!     Calculate the tensorial 2nd derivative of a plane wave vector
!!     and divides by |G|^2. pwdr2_gg(G=0) is put to zero.
!!
!!   AUTHOR
!!     RD (20-Nov-2006)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
!!     Adapted from pw_dr2
!!
!!*****
!******************************************************************************

SUBROUTINE pw_dr2_gg ( pw, pwdr2_gg, i, j )

    TYPE(pw_type), INTENT(IN)                :: pw
    TYPE(pw_type), INTENT(INOUT)             :: pwdr2_gg
    INTEGER, INTENT(IN)                      :: i, j

    INTEGER                                  :: cnt, handle, ig
    REAL(KIND=dp)                            :: flop, gg, o3

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_dr2_gg', 'I', 'Mflops', handle )

  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_dr2_gg")

  flop = 0.0_dp
  o3 = 1._dp/3._dp

  IF ( pw % in_space == RECIPROCALSPACE .AND. &
       pw % in_use == COMPLEXDATA1D ) THEN

    cnt = SIZE ( pw % cc )

    IF ( i == j ) THEN
!$omp parallel do private (ig)
      DO ig = pw%pw_grid % first_gne0, cnt
        gg = pw%pw_grid%g(i,ig)**2 - o3*pw%pw_grid%gsq(ig)
        pwdr2_gg%cc(ig) = gg * pw%cc(ig) / pw%pw_grid%gsq(ig)
      END DO
      flop = flop + 6 * cnt
    ELSE
!$omp parallel do private (ig)
      DO ig = pw%pw_grid % first_gne0, cnt
        pwdr2_gg%cc(ig) = pw%cc(ig) * (pw%pw_grid%g(i,ig)*pw%pw_grid%g(j,ig)) &
                          / pw%pw_grid%gsq(ig) 
      END DO
      flop = flop + 5 * cnt
    END IF

    IF ( pw%pw_grid % have_g0 ) pwdr2_gg%cc ( 1 ) = 0.0_dp

  ELSE

    CALL stop_program ( "pw_dr2_gg", "no suitable data field" )

  END IF

  flop = flop * 1.e-6_dp
  CALL timestop ( flop, handle )

END SUBROUTINE pw_dr2_gg

!******************************************************************************
!!****** pw_types/pw_smoothing [1.0] *
!!
!!   NAME
!!     pw_smoothing
!!
!!   SYNOPSIS
!!     Subroutine pw_smoothing(pw, ecut, sigma)
!!       Type(pw_type), Intent (INOUT)            :: pw
!!       REAL, Intent(IN)                         :: ecut, sigma
!!     End Subroutine pw_smoothing
!!
!!   FUNCTION
!!     Multiplies a G-space function with a smoothing factor of the form
!!     f(|G|) = exp((ecut - G^2)/sigma)/(1+exp((ecut - G^2)/sigma))
!!
!!   AUTHOR
!!     JGH (09-June-2006)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
!!
!!*****
!******************************************************************************

SUBROUTINE pw_smoothing ( pw, ecut, sigma )

    TYPE(pw_type), INTENT(INOUT)             :: pw
    REAL(KIND=dp), INTENT(IN)                :: ecut, sigma

    INTEGER                                  :: cnt, handle, ig
    REAL(KIND=dp)                            :: arg, f, flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_smooting', 'I', 'Mflops', handle )

  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_dr2")

  flop = 0.0_dp

  IF ( pw % in_space == RECIPROCALSPACE .AND. &
       pw % in_use == COMPLEXDATA1D ) THEN

    cnt = SIZE ( pw % cc )

!$omp parallel do private (ig)
    DO ig = 1, cnt
      arg = (ecut - pw%pw_grid%gsq(ig))/sigma
      f = EXP(arg)/(1+EXP(arg))
      pw%cc(ig) = f * pw%cc(ig)
    END DO
    flop = flop + 6 * cnt

  ELSE

    CALL stop_program ( "pw_smoothing", "no suitable data field" )

  END IF

  flop = flop * 1.e-6_dp
  CALL timestop ( flop, handle )

END SUBROUTINE pw_smoothing

!******************************************************************************
!!****** pw_types/pw_transfer [1.0] *
!!
!!   NAME
!!     pw_transfer
!!
!!   SYNOPSIS
!!     Subroutine pw_transfer(pw1, pw2)
!!       Type(pw_type), Intent (IN):: pw1
!!       Type(pw_type), Intent (INOUT):: pw2
!!     End Subroutine pw_transfer
!!
!!   FUNCTION
!!     Generalize copy of pw types
!!
!!   AUTHOR
!!     JGH (25-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (13-Mar-2001) : added gather/scatter cases
!!
!!   NOTES
!!     Copy routine that allows for in_space changes
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_transfer ( pw1, pw2, debug )

    TYPE(pw_type), INTENT(IN), TARGET        :: pw1
    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw2
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

!------------------------------------------------------------------------------

  CPPreconditionNoErr(pw1%ref_count>0,cp_failure_level,"pw_transfer")
  CPPreconditionNoErr(pw2%ref_count>0,cp_failure_level,"pw_transfer")

  IF ( pw1 % in_space == REALSPACE .AND. pw2 % in_space == REALSPACE ) THEN

     ! simple copy should do
     CALL pw_copy ( pw1, pw2 )

  ELSEIF ( pw1 % in_space == RECIPROCALSPACE .AND. &
         pw2 % in_space == RECIPROCALSPACE ) THEN

     IF ( pw1 % in_use == pw2 % in_use ) THEN

        ! simple copy should do
        CALL pw_copy ( pw1, pw2 )
     
     ELSE

       ! we have to gather/scatter the data
       IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
         CALL pw_scatter ( pw1, pw2 % cc3d )
       ELSEIF ( pw2 % in_use == COMPLEXDATA1D ) THEN
         CALL pw_gather ( pw2, pw1 % cc3d )
       ELSE
         CALL stop_program ( "pw_transfer", "Do not know what to do" )
       END IF

     END IF

  ELSE

     ! FFT needed, all further tests done in pw_fft_wrap
     CALL pw_fft_wrap ( pw1, pw2, debug )

  END IF

END SUBROUTINE pw_transfer

!!*****
!******************************************************************************
!!****** pw_types/pw_axpy [1.0] *
!!
!!   NAME
!!     pw_axpy
!!
!!   SYNOPSIS
!!     Subroutine pw_axpy(pw1, pw2, alpha)
!!       Type(pw_type), Intent (IN)           :: pw1
!!       Type(pw_type), Intent (INOUT)        :: pw2
!!       REAL(KIND=dp), Intent (IN), OPTIONAL :: alpha
!!     End Subroutine pw_axpy
!!
!!   FUNCTION
!!     pw2 = alpha*pw1 + pw2
!!     alpha defaults to 1
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (21-Feb-2003) : added reference grid functionality
!!     JGH (01-Dec-2007) : rename and remove complex alpha
!!
!!   NOTES
!!     Currently only summing up of respective types allowed,
!!     in order to avoid errors
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_axpy ( pw1, pw2, alpha )


    TYPE(pw_type), INTENT(IN)                :: pw1
    TYPE(pw_type), INTENT(INOUT)             :: pw2
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha

    INTEGER                                  :: handle, i, j, ng, ng1, ng2
    REAL(KIND=dp)                            :: flop, my_alpha

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_axpy', 'I', 'Mflops', handle )
  CPPreconditionNoErr(pw1%ref_count>0,cp_failure_level,"pw_axpy")
  CPPreconditionNoErr(pw2%ref_count>0,cp_failure_level,"pw_axpy")

  my_alpha=1.0_dp
  IF (PRESENT(alpha)) my_alpha=alpha

  IF ( pw1 % pw_grid %id_nr == pw2 % pw_grid %id_nr ) THEN

     IF ( pw1 % in_use == REALDATA1D .AND. pw2 % in_use == REALDATA1D ) THEN
        IF (my_alpha==1.0_dp) THEN
           !$omp parallel do private(i)
           DO i = 1, SIZE ( pw2 % cr )
              pw2 % cr(i) = pw2 % cr(i) + pw1 % cr(i)
           END DO
           flop = REAL ( SIZE ( pw2 % cr ),KIND=dp)
        ELSE
           !$omp parallel do private(i)
           DO i = 1, SIZE ( pw2 % cr )
              pw2 % cr(i) = pw2 % cr(i) + my_alpha*pw1 % cr(i)
           END DO
           flop = REAL ( 2*SIZE ( pw2 % cr ),KIND=dp)
        END IF
     ELSE IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
          pw2 % in_use == COMPLEXDATA1D ) THEN
        IF (my_alpha==1.0_dp) THEN
           !$omp parallel do private(i)
           DO i = 1, SIZE ( pw2 % cc )
              pw2 % cc(i) = pw2 % cc(i) + pw1 % cc(i)
           END DO
           flop = REAL ( 2 * SIZE ( pw2 % cc ),KIND=dp)
        ELSE
           !$omp parallel do private(i)
           DO i = 1, SIZE ( pw2 % cc )
              pw2 % cc(i) = pw2 % cc(i) + my_alpha*pw1 % cc(i)
           END DO
           flop = REAL ( 4 * SIZE ( pw2 % cc ),KIND=dp)
        END IF
     ELSE IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
        IF (my_alpha==1.0_dp) THEN
           pw2 % cr3d = pw2 % cr3d + pw1 % cr3d
           flop = REAL ( SIZE ( pw2 % cr3d ),KIND=dp)
        ELSE
           pw2 % cr3d = pw2 % cr3d + my_alpha * pw1 % cr3d
           flop = REAL ( 2*SIZE ( pw2 % cr3d ),KIND=dp)
        END IF
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
          pw2 % in_use == COMPLEXDATA3D ) THEN
        IF (my_alpha==1.0_dp) THEN
           pw2 % cc3d = pw2 % cc3d + pw1 % cc3d
           flop = REAL ( 2 * SIZE ( pw2 % cc3d ),KIND=dp)
        ELSE
           pw2 % cc3d = pw2 % cc3d + my_alpha*pw1 % cc3d
           flop = REAL ( 4 * SIZE ( pw2 % cc3d ),KIND=dp)
        END IF
     ELSE
        CALL stop_program ( "pw_axpy", "no suitable data field" )
     END IF

  ELSE IF ( pw_compatible ( pw1 % pw_grid, pw2 % pw_grid ) ) THEN

     IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
          pw2 % in_use == COMPLEXDATA1D .AND. &
          pw1 % in_space == RECIPROCALSPACE .AND. &
          pw1 % in_space == RECIPROCALSPACE ) THEN

        ng1 = SIZE ( pw1 % cc )
        ng2 = SIZE ( pw2 % cc )
        ng  = MIN ( ng1, ng2 )
        flop = REAL ( 2 * ng,KIND=dp)

        IF ( pw1 % pw_grid % spherical ) THEN

           IF (my_alpha==1.0_dp) THEN           
              !$omp parallel do private(i)
              DO i = 1, ng
                 pw2 % cc (i) = pw2 % cc (i) + pw1 % cc (i)
              END DO
           ELSE
              !$omp parallel do private(i)
              DO i = 1, ng
                 pw2 % cc (i) = pw2 % cc (i) + my_alpha*pw1 % cc (i)
              END DO
           END IF

        ELSEIF ( ( pw1 % pw_grid %id_nr == pw2 % pw_grid % reference ) ) THEN
  
           IF (my_alpha==1.0_dp) THEN
              IF( ng1 >= ng2 ) THEN
                 !$omp parallel do private(i,j)
                 DO i = 1, ng2
                    j = pw2 % pw_grid % gidx ( i )
                    pw2 % cc ( i ) = pw2 % cc ( i ) + pw1 % cc ( j )
                 END DO
              ELSE
                 !$omp parallel do private(i,j)
                 DO i = 1, ng1
                    j = pw2 % pw_grid % gidx ( i )
                    pw2 % cc ( j ) = pw2 % cc ( j ) + pw1 % cc ( i )
                 END DO
              END IF
           ELSE
              IF( ng1 >= ng2 ) THEN
                 !$omp parallel do private(i,j)
                 DO i = 1, ng2
                    j = pw2 % pw_grid % gidx ( i )
                    pw2 % cc ( i ) = pw2 % cc ( i ) + my_alpha*pw1 % cc ( j )
                 END DO
              ELSE
                 !$omp parallel do private(i,j)
                 DO i = 1, ng1
                    j = pw2 % pw_grid % gidx ( i )
                    pw2 % cc ( j ) = pw2 % cc ( j ) + my_alpha*pw1 % cc ( i )
                 END DO
              END IF
           END IF
        ELSEIF ( ( pw2 % pw_grid %id_nr == pw1 % pw_grid % reference ) ) THEN

           IF (my_alpha==1.0_dp) THEN
              IF( ng1 >= ng2 ) THEN
                 !$omp parallel do private(i,j)
                 DO i = 1, ng2
                    j = pw1 % pw_grid % gidx ( i )
                    pw2 % cc ( i ) = pw2 % cc ( i ) + pw1 % cc ( j )
                 END DO
              ELSE
                 !$omp parallel do private(i,j)
                 DO i = 1, ng1
                    j = pw1 % pw_grid % gidx ( i )
                    pw2 % cc ( j ) = pw2 % cc ( j ) + pw1 % cc ( i )
                 END DO
              END IF
           ELSE
              IF( ng1 >= ng2 ) THEN
                 !$omp parallel do private(i,j)
                 DO i = 1, ng2
                    j = pw1 % pw_grid % gidx ( i )
                    pw2 % cc ( i ) = pw2 % cc ( i ) + my_alpha* pw1 % cc ( j )
                 END DO
              ELSE
                 !$omp parallel do private(i,j)
                 DO i = 1, ng1
                    j = pw1 % pw_grid % gidx ( i )
                    pw2 % cc ( j ) = pw2 % cc ( j ) + my_alpha*pw1 % cc ( i )
                 END DO
              END IF
           END IF
        ELSE

           WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                " grid 1 :",pw1 % pw_grid %id_nr, &
                " sperical :",pw1 % pw_grid % spherical, &
                " reference :",pw1 % pw_grid % reference
           WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                " grid 2 :",pw2 % pw_grid %id_nr, &
                " sperical :",pw2 % pw_grid % spherical, &
                " reference :",pw2 % pw_grid % reference
           CALL stop_program ( "pw_axpy", "grids not compatible" )

        END IF

     ELSE
        CALL stop_program ( "pw_axpy", "no suitable data field" )
     END IF

  ELSE

     WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
          " grid 1 :",pw1 % pw_grid %id_nr, &
          " sperical :",pw1 % pw_grid % spherical, &
          " reference :",pw1 % pw_grid % reference
     WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
          " grid 2 :",pw2 % pw_grid %id_nr, &
          " sperical :",pw2 % pw_grid % spherical, &
          " reference :",pw2 % pw_grid % reference
     CALL stop_program ( "pw_axpy", "grids not compatible" )

  END IF

  flop = flop * 1.e-6_dp
  CALL timestop ( flop, handle )

END SUBROUTINE pw_axpy

!!*****
!******************************************************************************
!!****** pw_types/pw_gather [1.0] *
!!
!!   NAME
!!     pw_gather
!!
!!   FUNCTION
!!     Gathers the pw vector from a 3d data field
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_gather_s ( pw, c, scale )
 
    TYPE(pw_type), INTENT(INOUT)             :: pw
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale

    INTEGER                                  :: gpt, handle, l, m, n, ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat
    REAL(KIND=dp)                            :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_gather', 'I', 'Mcopy', handle )
  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_gather_s")

  IF ( pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_program ( "pw_gather_s", "Data field has to be COMPLEXDATA1D" )
  ENDIF

! after the gather we are in g-space
  pw % in_space = RECIPROCALSPACE

  mapl => pw % pw_grid % mapl % pos
  mapm => pw % pw_grid % mapm % pos
  mapn => pw % pw_grid % mapn % pos

  ngpts = SIZE ( pw % pw_grid % gsq  )

  ghat => pw % pw_grid % g_hat

  IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      pw % cc ( gpt ) = scale * c ( l, m, n )

    END DO

  ELSE

!$omp parallel do private(gpt,l,m,n)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      pw % cc ( gpt ) = c ( l, m, n )

    END DO

  END IF

  cpy = REAL ( ngpts,KIND=dp) * 1.e-6_dp
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_gather_s

!******************************************************************************

SUBROUTINE pw_gather_p ( pw, c, scale )


    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale

    INTEGER                                  :: gpt, handle, l, m, mn, n, &
                                                ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat, yzq
    REAL(KIND=dp)                            :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_gather', 'I', 'Mcopy', handle )
  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_gather_p")

  IF ( pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_program ( "pw_gather_p", "Data field has to be COMPLEXDATA1D" )
  ENDIF

  IF ( pw % pw_grid % para % mode /= PW_MODE_DISTRIBUTED ) THEN
     CALL stop_program ( "pw_gather_p", "This grid type is not distributed")
  ENDIF

! after the gather we are in g-space
  pw % in_space = RECIPROCALSPACE

  mapl => pw % pw_grid % mapl % pos
  mapm => pw % pw_grid % mapm % pos
  mapn => pw % pw_grid % mapn % pos

  ngpts = SIZE ( pw % pw_grid % gsq  )

  ghat => pw % pw_grid % g_hat
  yzq => pw % pw_grid % para % yzq

  IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n,mn)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      mn = yzq ( m, n )
      pw % cc ( gpt ) = scale * c ( l, mn )

    END DO

  ELSE

!$omp parallel do private(gpt,l,m,n,mn)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      mn = yzq ( m, n )
      pw % cc ( gpt ) = c ( l, mn )

    END DO

  END IF

  cpy = REAL ( ngpts,KIND=dp) * 1.e-6_dp
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_gather_p

!!*****
!******************************************************************************
!!****** pw_types/pw_scatter [1.0] *
!!
!!   NAME
!!     pw_scatter
!!
!!   FUNCTION
!!     Scatters a pw vector to a 3d data field
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_scatter_s ( pw, c, scale )


    TYPE(pw_type), INTENT(IN)                :: pw
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale

    INTEGER                                  :: gpt, handle, l, m, n, ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat
    REAL(KIND=dp)                            :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_scatter', 'I', 'Mcopy', handle )
  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_scatter_s")

  IF ( pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_program ( "pw_scatter_s", "Data field has to be COMPLEXDATA1D" )
  ENDIF

  IF ( pw % in_space /= RECIPROCALSPACE ) THEN
     CALL stop_program ( "pw_scatter_s", "Data has to be in RECIPROCALSPACE" )
  ENDIF

  mapl => pw % pw_grid % mapl % pos
  mapm => pw % pw_grid % mapm % pos
  mapn => pw % pw_grid % mapn % pos

  ghat => pw % pw_grid % g_hat

  ngpts = SIZE ( pw % pw_grid % gsq  )

  ! should only zero the unused bits (but the zero is needed)
  IF ( .NOT. PRESENT ( scale ) ) c = 0.0_dp

  IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      c ( l, m, n ) = scale * pw % cc ( gpt )

    END DO

  ELSE

!$omp parallel do private(gpt,l,m,n)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      c ( l, m, n ) = pw % cc ( gpt )

    END DO

  END IF

  IF ( pw % pw_grid % grid_span == HALFSPACE ) THEN

    mapl => pw % pw_grid % mapl % neg
    mapm => pw % pw_grid % mapm % neg
    mapn => pw % pw_grid % mapn % neg

    IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n)
      DO gpt = 1, ngpts

        l = mapl ( ghat ( 1, gpt ) ) + 1
        m = mapm ( ghat ( 2, gpt ) ) + 1
        n = mapn ( ghat ( 3, gpt ) ) + 1
        c ( l, m, n ) = scale * CONJG ( pw % cc ( gpt ) )

      END DO

    ELSE

!$omp parallel do private(gpt,l,m,n)
      DO gpt = 1, ngpts

        l = mapl ( ghat ( 1, gpt ) ) + 1
        m = mapm ( ghat ( 2, gpt ) ) + 1
        n = mapn ( ghat ( 3, gpt ) ) + 1
        c ( l, m, n ) = CONJG ( pw % cc ( gpt ) )

      END DO

    END IF

  END IF

  cpy = REAL ( ngpts,KIND=dp) * 1.e-6_dp
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_scatter_s

!******************************************************************************

SUBROUTINE pw_scatter_p ( pw, c, scale )


    TYPE(pw_type), INTENT(IN), TARGET        :: pw
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale

    INTEGER                                  :: gpt, handle, l, m, mn, n, &
                                                ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat, yzq
    REAL(KIND=dp)                            :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_scatter', 'I', 'Mcopy', handle )
  CPPreconditionNoErr(pw%ref_count>0,cp_failure_level,"pw_gather_p")

  IF ( pw % in_use /= COMPLEXDATA1D ) THEN
     CALL stop_program ( "pw_scatter_p", "Data field has to be COMPLEXDATA1D" )
  ENDIF

  IF ( pw % in_space /= RECIPROCALSPACE ) THEN
     CALL stop_program ( "pw_scatter_p", "Data has to be in RECIPROCALSPACE" )
  ENDIF

  IF ( pw % pw_grid % para % mode /= PW_MODE_DISTRIBUTED ) THEN
     CALL stop_program ( "pw_scatter_p", "This grid type is not distributed")
  ENDIF

  mapl => pw % pw_grid % mapl % pos
  mapm => pw % pw_grid % mapm % pos
  mapn => pw % pw_grid % mapn % pos

  ghat => pw % pw_grid % g_hat
  yzq => pw % pw_grid % para % yzq

  ngpts = SIZE ( pw % pw_grid % gsq  )

  IF ( .NOT. PRESENT ( scale ) ) c = 0.0_dp

  IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n,mn)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      mn = yzq ( m, n )
      c ( l, mn ) = scale * pw % cc ( gpt )

    END DO

  ELSE

!$omp parallel do private(gpt,l,m,n,mn)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      mn = yzq ( m, n )
      c ( l, mn ) = pw % cc ( gpt )

    END DO

  END IF

  IF ( pw % pw_grid % grid_span == HALFSPACE ) THEN

    mapm => pw % pw_grid % mapm % neg
    mapn => pw % pw_grid % mapn % neg
    mapl => pw % pw_grid % mapl % neg

    IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n,mn)
      DO gpt = 1, ngpts

        l = mapl ( ghat ( 1, gpt ) ) + 1
        m = mapm ( ghat ( 2, gpt ) ) + 1
        n = mapn ( ghat ( 3, gpt ) ) + 1
        mn = yzq ( m, n )
        c ( l, mn ) = scale * CONJG ( pw % cc ( gpt ) )

      END DO

    ELSE

!$omp parallel do private(gpt,l,m,n,mn)
      DO gpt = 1, ngpts

        l = mapl ( ghat ( 1, gpt ) ) + 1
        m = mapm ( ghat ( 2, gpt ) ) + 1
        n = mapn ( ghat ( 3, gpt ) ) + 1
        mn = yzq ( m, n )
        c ( l, mn ) = CONJG ( pw % cc ( gpt ) )

      END DO

    END IF

  END IF

  cpy = REAL ( ngpts,KIND=dp) * 1.e-6_dp
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_scatter_p

!!*****
!******************************************************************************
!!****** pw_types/pw_fft_wrap [1.1] *
!!
!!   NAME
!!     pw_fft_wrap
!!
!!   FUNCTION
!!     Generic function for 3d FFT of a coefficient_type or pw_type
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (30-12-2000): New setup of functions and adaptation to parallelism
!!     JGH (04-01-2001): Moved routine from pws to this module, only covers
!!                       pw_types, no more coefficient types
!!
!!   NOTES
!!     The following functions are covered
!!      fft_wrap_pw1, fft_wrap_pw1pw2
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE fft_wrap_pw1 ( pw1, debug )


    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw1
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    INTEGER                                  :: dir, handle, out_space
    LOGICAL                                  :: test
    REAL(KIND=dp)                            :: norm

!------------------------------------------------------------------------------

  CALL timeset("fft_wrap_pw1","I"," ",handle)
  CPPreconditionNoErr(pw1%ref_count>0,cp_failure_level,"fft_wrap_pw1")

  IF ( PRESENT ( debug ) ) THEN
    test = debug
  ELSE
    test = .FALSE.
  END IF

  IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
       pw1 % pw_grid % para % mode == PW_MODE_LOCAL ) THEN
!..what dirction will the transform be?
    IF ( pw1 % in_space == REALSPACE ) THEN
       dir = FWFFT
       norm = pw1 % pw_grid % dvol
       out_space = RECIPROCALSPACE
    ELSEIF ( pw1 % in_space == RECIPROCALSPACE ) THEN
       dir = BWFFT
       norm = 1.0_dp
       out_space = REALSPACE
    ELSE
       CALL stop_program ( "fft_wrap_pw1", "PW structure is missing a "//&
                       "proper tag to identidy its space" )
    END IF

    CALL fft3d ( dir, pw1 % pw_grid % npts, pw1 % cc3d, &
                 scale = norm, debug=test )

!..tag new data with correct space
    pw1 % in_space = out_space
  ELSE
    CALL stop_program ( "fft_wrap_pw1", "In place FFT only possible for "//&
                    "replicated data with COMPLEXDATA3D structure" )
  END IF
  CALL timestop(0.0_dp,handle) 

END SUBROUTINE fft_wrap_pw1


!******************************************************************************

SUBROUTINE fft_wrap_pw1pw2 ( pw1, pw2, debug )
  

    TYPE(pw_type), INTENT(IN), TARGET        :: pw1
    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw2
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    CHARACTER(LEN=9)                         :: mode
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: grays
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: c_in, c_out
    INTEGER                                  :: dir, handle, handle2, ierr, &
                                                l1, l2, l3, my_pos, &
                                                nloc( 3 ), nrays, nsize, &
                                                out_space
    INTEGER, DIMENSION(:), POINTER           :: n
    LOGICAL                                  :: test
    REAL(KIND=dp)                            :: norm

!------------------------------------------------------------------------------

  CALL timeset("fft_wrap_pw1pw2_all","I"," ",handle2)
  CALL timeset("fft_wrap_pw1pw2_"//TRIM(ADJUSTL(cp_to_string( &
               CEILING (pw1%pw_grid%cutoff/10)*10))),"I"," ",handle)

  CPPreconditionNoErr(pw1%ref_count>0,cp_failure_level,"fft_wrap_pw1pw2")
  CPPreconditionNoErr(pw2%ref_count>0,cp_failure_level,"fft_wrap_pw1pw2")
  NULLIFY ( c_in )
  NULLIFY ( c_out )

  IF ( PRESENT ( debug ) ) THEN
    test = debug
  ELSE
    test = .FALSE.
  END IF

!..check if grids are compatible
  IF ( pw1 % pw_grid %id_nr /= pw2 % pw_grid %id_nr ) THEN
    IF ( pw1 % pw_grid % dvol /= pw2 % pw_grid % dvol ) THEN
      CALL stop_program ( "fft_wrap_pw1pw2", "PW grids not compatible" )
    END IF
    IF ( pw1 % pw_grid % para %group /= pw2 % pw_grid % para % group) THEN
      CALL stop_program ( "fft_wrap_pw1pw2", &
                      "PW grids have not compatible MPI groups" )
    END IF
  END IF

!..prepare input
  IF ( pw1 % in_space == REALSPACE ) THEN
    dir = FWFFT
    norm = 1.0_dp / pw1 % pw_grid % ngpts
    out_space = RECIPROCALSPACE
  ELSE IF ( pw1 % in_space == RECIPROCALSPACE ) THEN
    dir = BWFFT
    norm = 1.0_dp
    out_space = REALSPACE
  ELSE
    CALL stop_program ( "fft_wrap_pw1pw2", "Error in space tag" )
  END IF

  n => pw1 % pw_grid % npts

  mode = fftselect ( pw1 % in_use, pw2 % in_use, pw1 % in_space )

  IF ( pw1 % pw_grid % para % mode == PW_MODE_LOCAL ) THEN

!
!..replicated data, use local FFT
!

    IF ( test ) THEN
      WRITE ( *,'(A)') " FFT Protocol "
      IF ( dir == FWFFT ) WRITE ( *,'(A,T76,A)') "  Transform direction ","FWFFT"
      IF ( dir == BWFFT ) WRITE ( *,'(A,T76,A)') "  Transform direction ","BWFFT"
      IF ( pw1 % in_space == REALSPACE ) &
         WRITE ( *,'(A,T72,A)') "  in space ","REALSPACE"
      IF ( pw1 % in_space == RECIPROCALSPACE ) &
         WRITE ( *,'(A,T66,A)') "  in space ","RECIPROCALSPACE"
      IF ( out_space == REALSPACE ) &
         WRITE ( *,'(A,T72,A)') "  out space ","REALSPACE"
      IF ( out_space == RECIPROCALSPACE ) &
         WRITE ( *,'(A,T66,A)') "  out space ","RECIPROCALSPACE"
      WRITE ( *,'(A,T66,E15.6)') "  scale factor",norm
    END IF

    SELECT CASE ( mode )
      CASE DEFAULT
        CALL stop_program ( "fft_wrap_pw1pw2", &
                        "Illegal combination of in_use and in_space" )
      CASE ( "FW_C3DC3D" )
        c_in => pw1 % cc3d
        c_out => pw2 % cc3d
        CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug=test )
      CASE ( "FW_R3DC3D" )
        pw2 % cc3d = CMPLX ( pw1 % cr3d, 0.0_dp, KIND=dp)
        c_out => pw2 % cc3d
        CALL fft3d ( dir, n, c_out, scale = norm, debug=test )
      CASE ( "FW_C3DC1D" )
        c_in => pw1 % cc3d
        IF (.NOT. ASSOCIATED (c_out)) &
            ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
        ! transform
        CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug=test )
        ! gather results
        IF ( test ) WRITE ( *,'(A)') "  PW_GATHER : 3d -> 1d "
        CALL pw_gather ( pw2, c_out )
        DEALLOCATE ( c_out, STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
      CASE ( "FW_R3DC1D" )
        IF (.NOT. ASSOCIATED (c_out)) &
            ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
        nsize=SIZE(pw1 % cr3d,1)*SIZE(pw1 % cr3d,2)*SIZE(pw1 % cr3d,3)
        l1 = LBOUND(pw1%cr3d,1)
        l2 = LBOUND(pw1%cr3d,2)
        l3 = LBOUND(pw1%cr3d,3)
        CALL copy_rc(nsize,pw1%cr3d(l1,l2,l3),c_out(1,1,1))
        CALL fft3d ( dir, n, c_out, scale = 1._dp, debug=test )
        IF ( test ) WRITE ( *,'(A)') "  PW_GATHER : 3d -> 1d "
        CALL pw_gather ( pw2, c_out, scale = norm )
        DEALLOCATE ( c_out, STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
      CASE ( "BW_C3DC3D" )
        c_in => pw1 % cc3d
        c_out => pw2 % cc3d
        CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug=test )
      CASE ( "BW_C3DR3D" )
        c_in => pw1 % cc3d
        IF (.NOT. ASSOCIATED (c_out)) &
            ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
        CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug=test )
        ! use real part only
        IF ( test ) WRITE ( *,'(A)') "  REAL part "
        pw2 % cr3d = REAL ( c_out,KIND=dp)
        DEALLOCATE ( c_out, STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
      CASE ( "BW_C1DC3D" )
        c_out => pw2 % cc3d
        IF ( test ) WRITE ( *,'(A)') "  PW_SCATTER : 3d -> 1d "
        CALL pw_scatter ( pw1, c_out )
        CALL fft3d ( dir, n, c_out, scale = norm, debug=test )
      CASE ( "BW_C1DR3D" )
        IF (.NOT. ASSOCIATED (c_out)) &
            ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
        IF ( test ) WRITE ( *,'(A)') "  PW_SCATTER : 3d -> 1d "
        CALL pw_scatter ( pw1, c_out )
        ! transform
        CALL fft3d ( dir, n, c_out, scale = norm, debug=test )
        ! use real part only
        IF ( test ) WRITE ( *,'(A)') "  REAL part "
        ! pw2 % cr3d = REAL ( c_out,KIND=dp)
        nsize=SIZE(pw2 % cr3d,1)*SIZE(pw2 % cr3d,2)*SIZE(pw2 % cr3d,3)
        l1 = LBOUND(pw2%cr3d,1)
        l2 = LBOUND(pw2%cr3d,2)
        l3 = LBOUND(pw2%cr3d,3)
        CALL copy_cr(nsize,c_out(1,1,1),pw2%cr3d(l1,l2,l3))
        DEALLOCATE ( c_out, STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
    END SELECT

    IF ( test ) WRITE ( *,'(A)') " End of FFT Protocol "

  ELSE

!
!..parallel FFT
!

    IF ( test .AND. pw1 % pw_grid % para % group_head ) THEN
      WRITE ( *,'(A)') " FFT Protocol "
      IF ( dir == FWFFT ) WRITE ( *,'(A,T76,A)') "  Transform direction ","FWFFT"
      IF ( pw1 % in_space == REALSPACE ) &
         WRITE ( *,'(A,T72,A)') "  in space ","REALSPACE"
      IF ( pw1 % in_space == RECIPROCALSPACE ) &
         WRITE ( *,'(A,T66,A)') "  in space ","RECIPROCALSPACE"
      IF ( out_space == REALSPACE ) &
         WRITE ( *,'(A,T72,A)') "  out space ","REALSPACE"
      IF ( out_space == RECIPROCALSPACE ) &
         WRITE ( *,'(A,T66,A)') "  out space ","RECIPROCALSPACE"
      WRITE ( *,'(A,T66,E15.6)') "  scale factor",norm
    END IF

    my_pos = pw1 % pw_grid % para % my_pos
    nrays = pw1 % pw_grid % para % nyzray ( my_pos )
    ALLOCATE ( grays ( n(1), nrays ), STAT = ierr )
    CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
    grays = 0.0_dp

    SELECT CASE ( mode )
      CASE DEFAULT
        CALL stop_program ( "fft_wrap_pw1pw2", &
                        "Illegal combination of in_use and in_space "//&
                        "in parallel 3d FFT" )
      CASE ( "FW_C3DC1D" )
        c_in => pw1 % cc3d
      CASE ( "FW_R3DC1D" )
        nloc = pw1 % pw_grid % npts_local
        ALLOCATE ( c_in( nloc(1), nloc(2), nloc(3) ), STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
        c_in = CMPLX ( pw1 % cr3d, 0.0_dp,KIND=dp)
      CASE ( "BW_C1DC3D" )
        IF ( test .AND. pw1%pw_grid%para%group_head ) &
            WRITE ( *,'(A)') "  PW_SCATTER : 2d -> 1d "
        CALL pw_scatter ( pw1, grays )
        c_in => pw2 % cc3d
      CASE ( "BW_C1DR3D" )
        IF ( test .AND. pw1%pw_grid%para%group_head ) &
            WRITE ( *,'(A)') "  PW_SCATTER : 2d -> 1d "
        CALL pw_scatter ( pw1, grays )
        nloc = pw2 % pw_grid % npts_local
        ALLOCATE ( c_in( nloc(1), nloc(2), nloc(3) ), STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
    END SELECT

!..transform
    IF ( pw1 % pw_grid % para % ray_distribution ) THEN
      CALL fft3d ( dir, n, c_in, grays, pw1 % pw_grid % para % group, &
                   pw1 % pw_grid % para % rs_group, &
                   pw1 % pw_grid % para % yzp, pw1 % pw_grid % para % nyzray, &
                   pw1 % pw_grid % para % bo, scale = norm, debug=test )
    ELSE
      CALL fft3d ( dir, n, c_in, grays, pw1 % pw_grid % para % rs_group, &
                   pw1 % pw_grid % para % bo, scale = norm, debug=test )
    END IF

!..prepare output
    SELECT CASE ( mode )
      CASE DEFAULT
        CALL stop_program ( "fft_wrap_pw1pw2", &
                        "Illegal combination of in_use and in_space "//&
                        "in parallel 3d FFT" )
      CASE ( "FW_C3DC1D" )
        IF ( test .AND. pw1%pw_grid%para%group_head ) &
            WRITE ( *,'(A)') "  PW_GATHER : 2d -> 1d "
        CALL pw_gather ( pw2, grays )
      CASE ( "FW_R3DC1D" )
        IF ( test .AND. pw1%pw_grid%para%group_head ) &
            WRITE ( *,'(A)') "  PW_GATHER : 2d -> 1d "
        CALL pw_gather ( pw2, grays )
        DEALLOCATE ( c_in, STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
      CASE ( "BW_C1DC3D" )
        ! nothing to do
      CASE ( "BW_C1DR3D" )
        IF ( test .AND. pw1%pw_grid%para%group_head ) &
            WRITE ( *,'(A)') "  Real part "
        pw2 % cr3d = REAL ( c_in,KIND=dp)
        DEALLOCATE ( c_in, STAT = ierr )
        CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")
    END SELECT

    DEALLOCATE ( grays, STAT = ierr )
    CPPreconditionNoErr(ierr == 0,cp_failure_level,"pw_zero")

  END IF

!..update the space tag for pw2
  pw2 % in_space = out_space

  IF ( test .AND. pw1 % pw_grid % para % group_head ) THEN
      WRITE ( *,'(A)') " End of FFT Protocol "
  END IF
  CALL timestop(0.0_dp,handle) 
  CALL timestop(0.0_dp,handle2) 

END SUBROUTINE fft_wrap_pw1pw2

!******************************************************************************

FUNCTION fftselect ( use1, use2, space1 ) RESULT ( mode )


    INTEGER, INTENT(IN)                      :: use1, use2, space1
    CHARACTER(LEN=9)                         :: mode

!------------------------------------------------------------------------------

  IF ( space1 == REALSPACE ) THEN
    mode ( 1 : 3 ) = "FW_"
  ELSE IF ( space1 == RECIPROCALSPACE ) THEN
    mode ( 1 : 3 ) = "BW_"
  ELSE
    CALL stop_program ( "fftselect", "Error in space tag" )
  END IF

  SELECT CASE ( use1 )
    CASE ( COMPLEXDATA3D )
      mode ( 4 : 6 ) = "C3D"
    CASE ( REALDATA3D )
      mode ( 4 : 6 ) = "R3D"
    CASE ( COMPLEXDATA1D )
      mode ( 4 : 6 ) = "C1D"
    CASE ( REALDATA1D )
      mode ( 4 : 6 ) = "R1D"
    CASE DEFAULT
      CALL stop_program ( "fftselect", "Error in use1 tag" )
  END SELECT

  SELECT CASE ( use2 )
    CASE ( COMPLEXDATA3D )
      mode ( 7 : 9 ) = "C3D"
    CASE ( REALDATA3D )
      mode ( 7 : 9 ) = "R3D"
    CASE ( COMPLEXDATA1D )
      mode ( 7 : 9 ) = "C1D"
    CASE ( REALDATA1D )
      mode ( 7 : 9 ) = "R1D"
    CASE DEFAULT
      CALL stop_program ( "fftselect", "Error in use1 tag" )
  END SELECT

END FUNCTION fftselect

!!*****
!******************************************************************************

!!****f* pw_types/pw_write [1.0] *
!!
!!   NAME
!!     pw_write
!!
!!   SYNOPSIS
!!     Subroutine pw_write(pw, unit_nr, error)
!!       Type(pw_type), Intent (IN):: pw
!!       Integer, Intent (IN):: unit_nr
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_write
!!
!!   FUNCTION
!!     writes a small description of the actual grid
!!     (change to output the data as cube file, maybe with an
!!     optional long_description arg?)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw: the pw data to output
!!     - unit_nr: the unit to output to
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_write(pw, unit_nr, error)
    TYPE(pw_type), INTENT(in)                :: pw
    INTEGER, INTENT(in)                      :: unit_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iostat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  WRITE (unit=unit_nr, fmt="('<pw>:{ id_nr=',i8,',')",iostat=iostat)&
       pw%id_nr

  SELECT CASE(pw%in_use)
  CASE (REALDATA1D)
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_use=REALDATA1D"
     IF (ASSOCIATED(pw%cr)) THEN
        WRITE (unit=unit_nr, fmt="(' cr=<real(',i8,':',i8,')at 0x',z16.16,'>')")&
             LBOUND(pw%cr,1),UBOUND(pw%cr,1),m_loc_r(pw%cr(LBOUND(pw%cr)))
     ELSE
        WRITE (unit=unit_nr, fmt="(' cr=*null*')")
     END IF
  CASE (REALDATA3D)
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_use=REALDATA3D"
     IF (ASSOCIATED(pw%cr3d)) THEN
        WRITE (unit=unit_nr, fmt="(' cr3d=<real(',i8,':',i8,',',i8,':',i8,',',i8,':',i8,')at 0x',z16.16,'>')")&
             LBOUND(pw%cr3d,1),UBOUND(pw%cr3d,1),LBOUND(pw%cr3d,2),UBOUND(pw%cr3d,2),&
             LBOUND(pw%cr3d,3),UBOUND(pw%cr3d,3),&
             m_loc_r(pw%cr3d)
     ELSE
        WRITE (unit=unit_nr, fmt="(' cr3d=*null*')")
     END IF
  CASE (COMPLEXDATA1D)
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_use=COMPLEXDATA1D"
     IF (ASSOCIATED(pw%cc)) THEN
        WRITE (unit=unit_nr, fmt="(' cc=<real(',i8,':',i8,') at 0x',z16.16,'>')")&
             LBOUND(pw%cc,1),UBOUND(pw%cc,1),m_loc_c(pw%cc)
     ELSE
       WRITE (unit=unit_nr, fmt="(' cc=*null*')")
    END IF
  CASE (COMPLEXDATA3D)
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_use=COMPLEXDATA3D"
     IF (ASSOCIATED(pw%cc3d)) THEN
        WRITE (unit=unit_nr, fmt="(' cc3d=<real(',i8,':',i8,',',i8,':',i8,',',i8,':',i8,') at 0x',z16.16,'>')")&
             LBOUND(pw%cc3d,1),UBOUND(pw%cc3d,1),LBOUND(pw%cc3d,2),UBOUND(pw%cc3d,2),&
             LBOUND(pw%cc3d,3),UBOUND(pw%cc3d,3),&
             m_loc_c(pw%cc3d)
     ELSE
        WRITE (unit=unit_nr, fmt="(' cr3d=*null*')")
     END IF
  CASE default
     WRITE (unit=unit_nr, fmt="(' in_use=',i8,',')",iostat=iostat)&
          pw%in_use
  END SELECT

  SELECT CASE(pw%in_space)
  CASE (NOSPACE)
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_space=NOSPACE"
  CASE (REALSPACE)
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_space=REALSPACE"
  CASE (RECIPROCALSPACE)
     WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_space=RECIPROCALSPACE"
  CASE default
     WRITE (unit=unit_nr, fmt="(' in_space=',i8,',')",iostat=iostat)&
          pw%in_space
  END SELECT

  WRITE (unit=unit_nr, fmt="(' pw_grid%id_nr=',i8,/,' }')",iostat=iostat)&
       pw%pw_grid%id_nr

END SUBROUTINE pw_write
!***************************************************************************

FUNCTION pw_compatible ( grida, gridb ) RESULT ( compat )
    TYPE(pw_grid_type), INTENT(IN)           :: grida, gridb
    LOGICAL                                  :: compat

   compat = .FALSE.
   IF ( grida%id_nr == gridb%id_nr ) THEN
      compat = .TRUE.
   ELSE IF ( grida%reference == gridb%id_nr ) THEN
      compat = .TRUE.
   ELSE IF ( gridb%reference == grida%id_nr ) THEN
      compat = .TRUE.
   END IF
   
END FUNCTION pw_compatible

!***************************************************************************

!!****f* pw_types/pw_retain [1.0] *
!!
!!   NAME
!!     pw_retain
!!
!!   SYNOPSIS
!!     Subroutine pw_retain(pw, error)
!!       Type(pw_type), Pointer:: pw
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_retain
!!
!!   FUNCTION
!!     retains a pw type
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - pw: the pw to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_retain(pw, error)
    TYPE(pw_type), POINTER                   :: pw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pw),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(pw%ref_count>0,cp_failure_level,routineP,error)
     pw%ref_count=pw%ref_count+1
  END IF
END SUBROUTINE pw_retain
!***************************************************************************

!!****f* pw_types/pw_release [1.0] *
!!
!!   NAME
!!     pw_release
!!
!!   SYNOPSIS
!!     Subroutine pw_release(pw, error)
!!       Type(pw_type), Pointer:: pw
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_release
!!
!!   FUNCTION
!!     releases the given pw
!!
!!   NOTES
!!     see doc/ReferenceCounting.html
!!
!!   ARGUMENTS
!!     - pw: the pw to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     04.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_release(pw, error)
    TYPE(pw_type), POINTER                   :: pw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(pw)) THEN
     CPPreconditionNoFail(pw%ref_count>0,cp_failure_level,routineP,error)
     pw%ref_count=pw%ref_count-1
     IF (pw%ref_count==0) THEN
        pw%ref_count=1

        allocated_pw_count = allocated_pw_count - 1
        SELECT CASE(pw % in_use)
        CASE (REALDATA1D)
           DEALLOCATE ( pw % cr, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CASE(COMPLEXDATA1D)
           DEALLOCATE ( pw % cc, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CASE(REALDATA3D)
           IF (ASSOCIATED(pw%cr3d)) THEN
              !FM optimizations of pools might have removed the 3d field to cache it
              DEALLOCATE ( pw % cr3d, STAT = stat )
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
        CASE(COMPLEXDATA3D)
           DEALLOCATE ( pw % cc3d, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CASE(NODATA)
        CASE default
           CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                "unknown data type "//cp_to_string(pw%in_use),error,failure)
        END SELECT
        pw%ref_count=0
        DEALLOCATE(pw, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(pw)
END SUBROUTINE pw_release
!***************************************************************************

!!****f* pw_types/pw_create [1.0] *
!!
!!   NAME
!!     pw_create
!!
!!   SYNOPSIS
!!     Subroutine pw_create(pw, pw_grid, use_data, in_space,&
!!         cr3d_ptr, error)
!!       Type(pw_type), Pointer:: pw
!!       Type(pw_grid_type), Pointer:: pw_grid
!!       Integer, Intent (IN):: use_data
!!       Integer, Intent (IN), Optional:: in_space
!!       REAL(KIND=dp), Dimension(:,:,:), Optional, Pointer:: cr3d_ptr
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_create
!!
!!   FUNCTION
!!     allocates and initializes pw_type
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pw: the type that will bw allocated and initialized
!!     - use_data: which kind of data will be used
!!     - in_space: in which space the pw is (real or reciprocal)
!!     - cr3d_ptr: pointer with the cr3d data (make sense only if 
!!       use_data==REALDATA3D)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_create(pw, pw_grid, use_data, in_space, cr3d_ptr, error)
    TYPE(pw_type), POINTER                   :: pw
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, INTENT(in)                      :: use_data
    INTEGER, INTENT(in), OPTIONAL            :: in_space
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: cr3d_ptr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bounds
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

  failure=.FALSE.
  
  CALL timeset(routineN,"I","",handle)
  CPPrecondition(.NOT.ASSOCIATED(pw),cp_failure_level,routineP,error,failure)
  ALLOCATE(pw,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
  logger => cp_error_get_logger(error)
  IF (debug_this_module) THEN
     WRITE (cp_logger_get_default_unit_nr(logger),"('*** allocated pw ***')")
     IF (PRESENT(cr3d_ptr)) THEN
        IF (ASSOCIATED(cr3d_ptr)) THEN
           WRITE (cp_logger_get_default_unit_nr(logger),"('*** cr3d associated ***')")
        END IF
     END IF
     CALL print_stack(cp_logger_get_default_unit_nr(logger))
  END IF

  IF (PRESENT(cr3d_ptr)) THEN
     IF (ASSOCIATED(cr3d_ptr)) THEN
        CPAssertNoFail(use_data==REALDATA3D,cp_failure_level,routineP,error)
     END IF
  END IF

  last_pw_id_nr=last_pw_id_nr+1
  pw % id_nr = last_pw_id_nr
  pw % ref_count = 1
  NULLIFY ( pw % pw_grid )
  pw % in_use = use_data
  pw % pw_grid => pw_grid
  pw % in_space = NOSPACE
  bounds => pw % pw_grid % bounds_local

  allocated_pw_count = allocated_pw_count + 1
  ! write(6,*) "allocated_pw_count ",allocated_pw_count

  NULLIFY ( pw % cr, pw % cc, pw % cr3d, pw % cc3d )

  SELECT CASE(use_data)
  CASE(REALDATA1D)
     ALLOCATE ( pw % cr ( pw % pw_grid % ngpts_cut_local ), STAT = stat )
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

  CASE(COMPLEXDATA1D)
     ALLOCATE ( pw % cc ( pw % pw_grid % ngpts_cut_local ), STAT = stat )
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

  CASE (REALDATA3D)
     IF (PRESENT(cr3d_ptr)) THEN
        IF (ASSOCIATED(cr3d_ptr)) THEN
           IF (ALL(bounds(1,:) <= bounds(2,:))) THEN
             CPPreconditionNoFail(ALL(LBOUND(cr3d_ptr)==bounds(1,:)),cp_failure_level,routineP,error)
             CPPreconditionNoFail(ALL(UBOUND(cr3d_ptr)==bounds(2,:)),cp_failure_level,routineP,error)
           END IF
           pw%cr3d => cr3d_ptr
        END IF
     END IF
     IF (.NOT.ASSOCIATED(pw%cr3d)) THEN
        ALLOCATE ( pw % cr3d ( &
             bounds ( 1, 1 ) : bounds ( 2, 1 ), &
             bounds ( 1, 2 ) : bounds ( 2, 2 ), &
             bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = stat )
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     END IF

  CASE(COMPLEXDATA3D)
     ALLOCATE ( pw % cc3d ( &
          bounds ( 1, 1 ) : bounds ( 2, 1 ), &
          bounds ( 1, 2 ) : bounds ( 2, 2 ), &
          bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = stat )
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  CASE(NODATA)
  CASE default
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "unknown data type",error,failure)
  END SELECT
  IF (PRESENT(in_space)) pw%in_space=in_space
  END IF
  CALL timestop(0.0_dp,handle)
END SUBROUTINE pw_create
!***************************************************************************

!!****f* pw_types/pw_compare_debug [1.0] *
!!
!!   NAME
!!     pw_compare_debug
!!
!!   SYNOPSIS
!!     Subroutine pw_compare_debug(pw1, pw2, maxdiff, error)
!!       Type(pw_type), Pointer:: pw1, pw2
!!       REAL(KIND=dp), Intent (OUT), Optional:: maxdiff
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine pw_compare_debug
!!
!!   FUNCTION
!!     compares two pw, only for debugging purposes
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pw1,pw2: the two pw to compare
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_compare_debug(pw1,pw2,maxdiff,error)
    TYPE(pw_type), POINTER                   :: pw1, pw2
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: maxdiff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_compare_debug', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: diff, mdiff
    TYPE(cp_logger_type), POINTER            :: logger

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pw1),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pw2),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw1%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(pw2%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(pw1%in_space==pw2%in_space,cp_warning_level,routineP,error,failure)
  CPPrecondition(pw1%in_use==pw2%in_use,cp_warning_level,routineP,error,failure)
  CALL cp_assert(ALL(pw1%pw_grid%bounds_local==pw2%pw_grid%bounds_local),&
       cp_failure_level,cp_assertion_failed,routineP,&
       "wrong pw distribution",error,failure)
  IF (.NOT. failure) THEN
     logger => cp_error_get_logger(error)
     unit_nr=-1
     mdiff=0.0_dp
     SELECT CASE(pw1%in_use)
     CASE(REALDATA3D)
        DO k=pw1%pw_grid%bounds_local(1,3),pw1%pw_grid%bounds_local(2,3)
           DO j=pw1%pw_grid%bounds_local(1,2),pw1%pw_grid%bounds_local(2,2)
              DO i=pw1%pw_grid%bounds_local(1,1),pw1%pw_grid%bounds_local(2,1)
                 diff=ABS(pw1%cr3d(i,j,k)-pw2%cr3d(i,j,k))
                 IF (mdiff<diff) THEN
                    WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                         fmt="(' diff=',e10.4,'at(',i4,i4,i4,'),',e10.4,'vs',e10.4)")&
                         diff,i,j,k,pw1%cr3d(i,j,k),pw2%cr3d(i,j,k)
                    mdiff=diff
                 END IF
              END DO
           END DO
        END DO
     CASE(COMPLEXDATA1D)
        DO i=1,pw1%pw_grid%ngpts_local
           diff=ABS(pw1%cc(i)-pw2%cc(i))
           IF (mdiff<diff) THEN
              WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                   fmt="(' diff=',e10.4,'at(',i4,i4,i4,'),',e9.3,e9.3,'vs',e9.3,e9.3)")&
                   diff,i,pw1%cc(i),pw2%cc(i)
              mdiff=diff
           END IF
        END DO
     CASE default
        CPPrecondition(.FALSE.,cp_warning_level,routineP,error,failure)
     END SELECT
     WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
          fmt="(' maxdiff=',e10.4)") mdiff
     IF (PRESENT(maxdiff)) maxdiff=mdiff
  ELSE
     WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
          fmt="(' incompatible pws')")
     IF (PRESENT(maxdiff)) maxdiff=HUGE(0.0_dp)
  END IF
END SUBROUTINE pw_compare_debug

!***************************************************************************
!!****** pws/pw_integral [1.0] *
!!
!!   NAME
!!     pw_integral
!!
!!   FUNCTION
!!     Calculate integral over unit cell for functions in plane wave basis
!!     only returns the real part of it ......
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (14-Mar-2001) : Parallel sum and some tests, HALFSPACE case
!!
!******************************************************************************

FUNCTION pw_integral_ab ( pw1, pw2, flag ) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    INTEGER, INTENT(IN), OPTIONAL            :: flag
    REAL(KIND=dp)                            :: integral_value

!------------------------------------------------------------------------------

  IF ( pw1 % pw_grid %id_nr /= pw2 % pw_grid %id_nr ) THEN
     CALL stop_program ( "integral_ab", "grids incompatible" )
  END IF

  ! since the return value is real, only do accurate sum on the real bit ?
  IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
     integral_value = accurate_sum ( pw1 % cr3d ( :, :, : )  &
          * pw2 % cr3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == REALDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     integral_value = accurate_sum (  pw1 % cr3d ( :, :, : )  &
          * pw2 % cc3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == REALDATA3D ) THEN
     integral_value = accurate_sum ( pw1 % cc3d ( :, :, : ) &
          *  pw2 % cr3d ( :, :, : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D &
       .AND. pw2 % in_use == COMPLEXDATA3D ) THEN
     integral_value = accurate_sum ( CONJG ( pw1 % cc3d ( :, :, : ) ) &
          * pw2 % cc3d ( :, :, : ) )

  ELSE IF ( pw1 % in_use == REALDATA1D &
       .AND. pw2 % in_use == REALDATA1D ) THEN
     integral_value = accurate_sum ( pw1 % cr ( : ) *  pw2 % cr ( : ) )
  ELSE IF ( pw1 % in_use == REALDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = accurate_sum ( pw1 % cr ( : ) * pw2 % cc ( : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == REALDATA1D ) THEN
     integral_value = accurate_sum ( pw1 % cc ( : ) * pw2 % cr ( : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D &
       .AND. pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = accurate_sum ( CONJG ( pw1 % cc ( : ) ) * pw2 % cc ( : ) )
  ELSE
     CALL stop_program ( "integral_ab", "no possible DATA" )
  END IF

  IF ( pw1 % in_use == REALDATA3D .OR. pw1 % in_use == COMPLEXDATA3D ) THEN
     integral_value = integral_value * pw1 % pw_grid % dvol
  ELSE
     integral_value = integral_value * pw1 % pw_grid % vol
  ENDIF
  IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
     IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
        integral_value = 2.0_dp * integral_value
        IF ( pw1 % pw_grid % have_g0 ) integral_value = integral_value - &
             CONJG ( pw1 % cc ( 1 ) ) * pw2 % cc ( 1 )
     END IF
  END IF

  IF ( pw1 % pw_grid % para % mode == PW_MODE_DISTRIBUTED ) &
       CALL mp_sum ( integral_value, pw1 % pw_grid % para % group )

END FUNCTION pw_integral_ab

!******************************************************************************

FUNCTION pw_integral_aa ( pw1, flag ) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1
    INTEGER, INTENT(IN), OPTIONAL            :: flag
    REAL(KIND=dp)                            :: integral_value

!------------------------------------------------------------------------------

  IF ( PRESENT ( flag ) ) THEN
     IF ( flag == SQUARE ) THEN
        IF ( pw1 % in_use == REALDATA3D ) THEN
           integral_value = accurate_sum ( pw1 % cr3d ( :, :, : )** 2 )
        ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
           integral_value = accurate_sum ( CONJG ( pw1 % cc3d ( :, :, : ) ) &
                * pw1 % cc3d ( :, :, : )  )
        ELSE IF ( pw1 % in_use == REALDATA1D ) THEN
           integral_value = accurate_sum ( pw1 % cr ( : ) ** 2 )
        ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
           integral_value = accurate_sum ( CONJG ( pw1 % cc ( : ) ) &
                * pw1 % cc ( : ) )
           IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
              integral_value = 2.0_dp * integral_value
              IF ( pw1 % pw_grid % have_g0 ) integral_value = integral_value - &
                   CONJG ( pw1 % cc ( 1 ) ) * pw1 % cc ( 1 )
           END IF
        ELSE
           CALL stop_program ( "integral_aa", "no possible SQUARE DATA" )
        END IF

     ELSE IF ( flag == SQUAREROOT ) THEN
        CALL stop_program ( "integral_aa", "no SQUAREROOT defined" )
     END IF

  ELSE
     IF ( pw1 % in_use == REALDATA3D ) THEN
        integral_value = accurate_sum ( pw1 % cr3d ( :, :, : ) )
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
        integral_value = accurate_sum ( pw1 % cc3d ( :, :, : ) )
     ELSE IF ( pw1 % in_use == REALDATA1D ) THEN
        integral_value = accurate_sum ( pw1 % cr ( : ) )
     ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
        integral_value = accurate_sum ( pw1 % cc ( : ) )
     ELSE
        CALL stop_program ( "integral_aa", "no possible DATA" )
     END IF
  END IF

  IF ( pw1 % in_use == REALDATA3D .OR. pw1 % in_use == COMPLEXDATA3D ) THEN
     integral_value = integral_value * pw1 % pw_grid % dvol
  ELSE
     integral_value = integral_value * pw1 % pw_grid % vol
  END IF

  IF ( pw1 % pw_grid % para % mode == PW_MODE_DISTRIBUTED ) &
       CALL mp_sum ( integral_value, pw1 % pw_grid % para % group )

END FUNCTION pw_integral_aa

!******************************************************************************

FUNCTION pw_integral_a2b ( pw1, pw2 ) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    REAL(KIND=dp)                            :: integral_value

!------------------------------------------------------------------------------

  IF ( pw1 % pw_grid %id_nr /= pw2 % pw_grid %id_nr ) THEN
     CALL stop_program ( "integral_a2b", "grids incompatible" )
  END IF
  IF ( pw1 % in_use == REALDATA1D .AND. &
       pw2 % in_use == REALDATA1D ) THEN
     integral_value = accurate_sum ( pw1 % cr ( : ) * pw2 % cr ( : ) &
          * pw1 % pw_grid % gsq ( : ) )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
       pw2 % in_use == COMPLEXDATA1D ) THEN
     integral_value = accurate_sum ( REAL ( CONJG ( pw1 % cc ( : ) ) &
          *  pw2 % cc ( : ) ,KIND=dp) * pw1 % pw_grid % gsq ( : ) )
     IF ( pw1 % pw_grid % grid_span == HALFSPACE ) THEN
        integral_value = 2.0_dp * integral_value
     END IF
  ELSE
     CALL stop_program ( "integral_a2b", "no possible DATA" )
  END IF

  IF ( pw1 % in_use == REALDATA3D .OR. pw1 % in_use == COMPLEXDATA3D ) THEN
     integral_value = integral_value * pw1 % pw_grid % dvol
  ELSE
     integral_value = integral_value * pw1 % pw_grid % vol
  END IF

  IF ( pw1 % pw_grid % para % mode == PW_MODE_DISTRIBUTED ) &
       CALL mp_sum ( integral_value, pw1 % pw_grid % para % group )

END FUNCTION pw_integral_a2b

!******************************************************************************
!!****** pw_types/pw_structure_factor [1.0] *
!!
!!   NAME
!!     pw_structure_factor
!!
!!   FUNCTION
!!     Calculate the structure factor for point r
!!
!!   AUTHOR
!!     JGH (05-May-2006)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
!!
!****
!******************************************************************************

SUBROUTINE pw_structure_factor ( sf, r )

    TYPE(pw_type), INTENT(INOUT)             :: sf
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: r

    INTEGER                                  :: cnt, handle, ig
    REAL(KIND=dp)                            :: arg, flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_structure_factor', 'I', 'Mflops', handle )

  CPPreconditionNoErr(sf%ref_count>0,cp_failure_level,"sf")

  flop = 0.0_dp

  IF ( sf % in_space == RECIPROCALSPACE .AND. &
       sf % in_use == COMPLEXDATA1D ) THEN

    cnt = SIZE ( sf % cc )

!$omp parallel do private (ig)
    DO ig = 1, cnt
      arg = DOT_PRODUCT(sf%pw_grid%g(:,ig),r)
      sf%cc(ig)=CMPLX(COS(arg),-SIN(arg),KIND=dp)
    END DO
    flop = flop + 7 * cnt
  ELSE

    CALL stop_program ( "pw_structure_factor", "no suitable data field" )

  END IF

  flop = flop * 1.e-6_dp
  CALL timestop ( flop, handle )

END SUBROUTINE pw_structure_factor

!******************************************************************************

FUNCTION pw_integrate_function(fun,isign,OPERATOR) RESULT(total_fun)

    TYPE(pw_type), INTENT(IN)                :: fun
    INTEGER, INTENT(IN), OPTIONAL            :: isign
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: OPERATOR
    REAL(KIND=dp)                            :: total_fun

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_integrate_function', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iop

!   ---------------------------------------------------------------------------

    iop = 0
    IF ( PRESENT(OPERATOR) ) THEN
      SELECT CASE (OPERATOR)
        CASE ("ABS","abs")
          iop = 1
        CASE DEFAULT
          CALL stop_program(routineN,"Unknown operator")
      END SELECT
    END IF

    total_fun = 0._dp

    IF (fun%in_space == REALSPACE) THEN
      IF (fun%in_use == REALDATA3D) THEN
        ! do reduction using maximum accuracy
        IF (iop==1) THEN
          total_fun = fun%pw_grid%dvol*accurate_sum(ABS(fun%cr3d))
        ELSE
          total_fun = fun%pw_grid%dvol*accurate_sum(fun%cr3d)
        END IF
      ELSE
        CALL stop_program(routineN,"in_space/in_use combination not implemented")
      END IF
    ELSEIF (fun%in_space == RECIPROCALSPACE) THEN
      IF (iop==1) &
        CALL stop_program(routineN,"Operator ABS not implemented")
      IF (fun%in_use == COMPLEXDATA3D) THEN
        total_fun = fun%pw_grid%vol*fun%cc3d(      &
                            fun%pw_grid%bounds(1,1), &
                            fun%pw_grid%bounds(1,2), &
                            fun%pw_grid%bounds(1,3))
      ELSEIF (fun%in_use == COMPLEXDATA1D) THEN
        IF ( fun%pw_grid%have_g0 ) total_fun = fun%pw_grid%vol*fun%cc(1)
      ELSE
        CALL stop_program(routineN,"in_space/in_use combination not implemented")
      END IF
    ELSE
      CALL stop_program(routineN,"No space defined")
    END IF
    IF (fun%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
       CALL mp_sum(total_fun,fun%pw_grid%para%group)
    END IF
    IF ( PRESENT(isign) ) THEN
       total_fun = total_fun * SIGN(1._dp,REAL(isign,dp))
    END IF

END FUNCTION pw_integrate_function

!******************************************************************************

END MODULE pw_types

!******************************************************************************
