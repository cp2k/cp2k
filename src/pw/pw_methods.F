
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \note
!>      If parallel mode is distributed certain combination of
!>      "in_use" and "in_space" can not be used.
!>      For performance reasons it would be better to have the loops
!>      over g-vectros in the gather/scatter routines in new subprograms
!>      with the actual arrays (also the adressing) in the parameter list
!> \par History
!>      JGH (29-Dec-2000) : Changes for parallel use
!>      JGH (13-Mar-2001) : added timing calls
!>      JGH (26-Feb-2003) : OpenMP enabled
!>      JGH (17-Nov-2007) : Removed mass arrays
!>      JGH (01-Dec-2007) : Removed and renamed routines
!>      03.2008 [tlaino] : Splitting pw_types into pw_types and pw_methods
!> \author apsi
! *****************************************************************************
MODULE pw_methods
  USE f77_blas
  USE fft_tools,                       ONLY: BWFFT,&
                                             FWFFT,&
                                             fft3d
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE message_passing,                 ONLY: mp_sum
  USE pw_cuda,                         ONLY: pw_cuda_c1dr3d_3d,&
                                             pw_cuda_c1dr3d_3d_ps,&
                                             pw_cuda_r3dc1d_3d,&
                                             pw_cuda_r3dc1d_3d_ps
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             PW_MODE_DISTRIBUTED,&
                                             PW_MODE_LOCAL,&
                                             pw_grid_type
  USE pw_types,                        ONLY: &
       COMPLEXDATA1D, COMPLEXDATA3D, NOSPACE, REALDATA1D, REALDATA3D, &
       REALSPACE, RECIPROCALSPACE, SQUARE, SQUAREROOT, pw_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pw_zero, pw_structure_factor, pw_smoothing
  PUBLIC :: pw_copy, pw_axpy, pw_transfer, pw_scale
  PUBLIC :: pw_derive, pw_dr2, pw_write
  PUBLIC :: pw_integral_ab, pw_integral_a2b
  PUBLIC :: pw_integral_aa_fast, pw_integral_ab_fast
  PUBLIC :: pw_dr2_gg, pw_integrate_function

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_methods'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

  INTERFACE pw_gather
     MODULE PROCEDURE pw_gather_s, pw_gather_p
  END INTERFACE

  INTERFACE pw_scatter
     MODULE PROCEDURE pw_scatter_s, pw_scatter_p
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief Set values of a pw type to zero
!> \param pw ...
!> \param error ...
!> \par History
!>      none
!> \author apsi
! *****************************************************************************
  SUBROUTINE pw_zero ( pw, error )

    TYPE(pw_type), INTENT(INOUT)             :: pw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_zero', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ns
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: zr

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( pw%in_use == REALDATA1D ) THEN
       ns = SIZE(pw%cr)
       pw%cr(:) = 0._dp
    ELSE IF ( pw%in_use == COMPLEXDATA1D ) THEN
       ns = SIZE(pw%cc)
       pw%cc(:) = CMPLX(0._dp,0._dp,KIND=dp)
    ELSE IF ( pw%in_use == REALDATA3D ) THEN
       ns = SIZE(pw%cr3d)
       pw%cr3d(:,:,:) = 0._dp
    ELSE IF ( pw%in_use == COMPLEXDATA3D ) THEN
       ns = SIZE(pw%cc3d)
       pw%cc3d(:,:,:) = CMPLX(0._dp,0._dp,KIND=dp)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"No possible data field!")
    END IF

    zr = REAL ( ns,KIND=dp) * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_zero

! *****************************************************************************
!> \brief copy a pw type variable
!> \param pw1 ...
!> \param pw2 ...
!> \param error ...
!> \par History
!>      JGH (7-Mar-2001) : check for pw_grid %id_nr, allow copy if
!>        in_use == COMPLEXDATA1D and in_space == RECIPROCALSPACE
!>      JGH (21-Feb-2003) : Code for generalized reference grids
!> \author apsi
!> \note
!>      Currently only copying of respective types allowed,
!>      in order to avoid errors
! *****************************************************************************
  SUBROUTINE pw_copy ( pw1, pw2, error)
    TYPE(pw_type), INTENT(IN)                :: pw1
    TYPE(pw_type), INTENT(INOUT)             :: pw2
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ng, ng1, ng2, ns
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: zc

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw1%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(pw2%ref_count>0,cp_failure_level,routineP,error,failure)
    IF ( pw1%pw_grid %id_nr /= pw2%pw_grid %id_nr ) THEN

       IF ( pw1%pw_grid%spherical .AND. pw2%pw_grid%spherical ) THEN

          IF ( pw_compatible (pw1%pw_grid, pw2%pw_grid, error=error) ) THEN

             IF ( pw1%in_use == COMPLEXDATA1D .AND. &
                  pw2%in_use == COMPLEXDATA1D .AND. &
                  pw1%in_space == RECIPROCALSPACE ) THEN
                ng1 = SIZE ( pw1%cc )
                ng2 = SIZE ( pw2%cc )
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
                CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")
             END IF

          ELSE
             WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                        " grid 1 :",pw1%pw_grid %id_nr, &
                        " spherical :",pw1%pw_grid%spherical, &
                        " reference :",pw1%pw_grid%reference
             WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                        " grid 2 :",pw2%pw_grid %id_nr, &
                        " spherical :",pw2%pw_grid%spherical, &
                        " reference :",pw2%pw_grid%reference
             CALL stop_program(routineN,moduleN,__LINE__,"Incompatible grids")
          END IF

       ELSE IF ( .NOT. ( pw1%pw_grid%spherical .OR. &
                         pw2%pw_grid%spherical ) ) THEN

          ng1 = SIZE ( pw1%cc )
          ng2 = SIZE ( pw2%cc )
          ns = 2 * MAX ( ng1, ng2 )

          IF ( pw1%in_use == COMPLEXDATA1D .AND. &
               pw2%in_use == COMPLEXDATA1D .AND. &
               pw1%in_space == RECIPROCALSPACE ) THEN

             IF ( ( pw1%pw_grid %id_nr == pw2%pw_grid%reference ) ) THEN
                IF( ng1 >= ng2 ) THEN
!$omp parallel do private(i,j)
                   DO i = 1, ng2
                      j = pw2%pw_grid%gidx ( i )
                      pw2%cc ( i ) = pw1%cc ( j )
                   END DO
                ELSE
                   CALL pw_zero ( pw2, error=error)
!$omp parallel do private(i,j)
                   DO i = 1, ng1
                      j = pw2%pw_grid%gidx ( i )
                      pw2%cc ( j ) = pw1%cc ( i )
                   END DO
                END IF
             ELSE IF ( ( pw2%pw_grid %id_nr == pw1%pw_grid%reference ) ) THEN
                IF( ng1 >= ng2 ) THEN
!$omp parallel do private(i,j)
                   DO i = 1, ng2
                      j = pw1%pw_grid%gidx ( i )
                      pw2%cc ( i ) = pw1%cc ( j )
                   END DO
                ELSE
                   CALL pw_zero ( pw2, error=error)
!$omp parallel do private(i,j)
                   DO i = 1, ng1
                      j = pw1%pw_grid%gidx ( i )
                      pw2%cc ( j ) = pw1%cc ( i )
                   END DO
                END IF
             ELSE
                WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                           " grid 1 :",pw1%pw_grid %id_nr, &
                           " spherical :",pw1%pw_grid%spherical, &
                           " reference :",pw1%pw_grid%reference
                WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                           " grid 2 :",pw2%pw_grid %id_nr, &
                           " spherical :",pw2%pw_grid%spherical, &
                           " reference :",pw2%pw_grid%reference
                CALL stop_program(routineN,moduleN,__LINE__,"Incompatible grids")
             END IF

          ELSE
             CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")
          END IF

          pw2%in_space = RECIPROCALSPACE

       ELSE
          WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                     " grid 1 :",pw1%pw_grid %id_nr, &
                     " spherical :",pw1%pw_grid%spherical, &
                     " reference :",pw1%pw_grid%reference
          WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                     " grid 2 :",pw2%pw_grid %id_nr, &
                     " spherical :",pw2%pw_grid%spherical, &
                     " reference :",pw2%pw_grid%reference
          CALL stop_program(routineN,moduleN,__LINE__,"Incompatible grids")
       END IF

    ELSE

       IF ( pw1%in_use == REALDATA1D .AND. pw2%in_use == REALDATA1D ) THEN
          ns = SIZE ( pw1%cr )
!$omp parallel do private(i)
          DO i = 1, ns
             pw2%cr(i) = pw1%cr(i)
          END DO
       ELSE IF ( pw1%in_use == COMPLEXDATA1D .AND. &
            pw2%in_use == COMPLEXDATA1D ) THEN
          ns = SIZE ( pw1%cc )
!$omp parallel do private(i)
          DO i = 1, ns
             pw2%cc(i) = pw1%cc(i)
          END DO
       ELSE IF ( pw1%in_use == REALDATA3D .AND. pw2%in_use == REALDATA3D ) THEN
          ns = SIZE ( pw1%cr3d )
          pw2%cr3d(:,:,:) = pw1%cr3d(:,:,:)
       ELSE IF ( pw1%in_use == COMPLEXDATA3D .AND. &
            pw2%in_use == COMPLEXDATA3D ) THEN
          ns = SIZE ( pw1%cc3d )
          pw2%cc3d(:,:,:) = pw1%cc3d(:,:,:)
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")
       END IF

    END IF

    pw2%in_space = pw1%in_space
    zc = REAL ( ns,KIND=dp) * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_copy

! *****************************************************************************
!> \brief multiplies pw coeffs with a number
!> \param pw ...
!> \param a ...
!> \param error ...
!> \par History
!>      11.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE pw_scale( pw, a , error)

    TYPE(pw_type), INTENT(INOUT)             :: pw
    REAL(KIND=dp), INTENT(IN)                :: a
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_scale', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ns
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: flop

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)

    SELECT CASE ( pw%in_use )
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
       CALL stop_program(routineN,moduleN,__LINE__, "No suitable data field" )
    END SELECT

    flop = REAL ( ns, KIND=dp) * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_scale

! *****************************************************************************
!> \brief Calculate the derivative of a plane wave vector
!> \param pw ...
!> \param n ...
!> \param error ...
!> \par History
!>      JGH (06-10-2002) allow only for inplace derivatives
!> \author JGH (25-Feb-2001)
!> \note
!>      Calculate the derivative dx^n(1) dy^n(2) dz^n(3) PW
!>      PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
! *****************************************************************************
  SUBROUTINE pw_derive ( pw, n, error )

    TYPE(pw_type), INTENT(INOUT)             :: pw
    INTEGER, DIMENSION(3), INTENT(IN)        :: n
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_derive', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: im
    INTEGER                                  :: cnt, handle, i, m
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: flop

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ALL(n>=0),cp_failure_level,routineP,error,failure)

    m = SUM ( n )
    im = CMPLX ( 0.0_dp, 1.0_dp,KIND=dp) ** m

    flop = 0.0_dp

    IF ( pw%in_space == RECIPROCALSPACE .AND. &
         pw%in_use == COMPLEXDATA1D ) THEN

       cnt = SIZE ( pw%cc )

       IF ( n ( 1 ) == 1 ) THEN
!$omp parallel do private (i)
          DO i = 1, cnt
             pw%cc ( i ) = pw%cc ( i ) * pw%pw_grid%g ( 1, i )
          END DO
          flop = flop + 6 * cnt
       ELSE IF ( n ( 1 ) > 1 ) THEN
!$omp parallel do private (i)
          DO i = 1, cnt
             pw%cc ( i ) = pw%cc ( i ) * ( pw%pw_grid%g ( 1, i ) ** n ( 1 ) )
          END DO
          flop = flop + 7 * cnt
       END IF
       IF ( n ( 2 ) == 1 ) THEN
!$omp parallel do private (i)
          DO i = 1, cnt
             pw%cc ( i ) = pw%cc ( i ) * pw%pw_grid%g ( 2, i )
          END DO
          flop = flop + 6 * cnt
       ELSE IF ( n ( 2 ) > 1 ) THEN
!$omp parallel do private (i)
          DO i = 1, cnt
             pw%cc ( i ) = pw%cc ( i ) * ( pw%pw_grid%g ( 2, i ) ** n ( 2 ) )
          END DO
          flop = flop + 7 * cnt
       END IF
       IF ( n ( 3 ) == 1 ) THEN
!$omp parallel do private (i)
          DO i = 1, cnt
             pw%cc ( i ) = pw%cc ( i ) * pw%pw_grid%g ( 3, i )
          END DO
          flop = flop + 6 * cnt
       ELSE IF ( n ( 3 ) > 1 ) THEN
!$omp parallel do private (i)
          DO i = 1, cnt
             pw%cc ( i ) = pw%cc ( i ) * ( pw%pw_grid%g ( 3, i ) ** n ( 3 ) )
          END DO
          flop = flop + 7 * cnt
       END IF

       ! im can take the values 1, -1, i, -i
       ! skip this if im == 1
       IF ( ABS ( REAL ( im,KIND=dp) - 1.0_dp ) > 1.e-10 ) THEN
!$omp parallel do private (i)
          DO i = 1, cnt
             pw%cc ( i ) = im * pw%cc ( i )
          END DO
          flop = flop + 6 * cnt
       END IF

    ELSE

       CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")

    END IF

    flop = flop * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_derive

! *****************************************************************************
!> \brief Calculate the tensorial 2nd derivative of a plane wave vector
!> \param pw ...
!> \param pwdr2 ...
!> \param i ...
!> \param j ...
!> \param error ...
!> \par History
!>      none
!> \author JGH (05-May-2006)
!> \note
!>      PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
! *****************************************************************************
  SUBROUTINE pw_dr2 ( pw, pwdr2, i, j, error)

    TYPE(pw_type), INTENT(IN)                :: pw
    TYPE(pw_type), INTENT(INOUT)             :: pwdr2
    INTEGER, INTENT(IN)                      :: i, j
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_dr2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cnt, handle, ig
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: flop, gg, o3

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)

    flop = 0.0_dp
    o3 = 1._dp/3._dp

    IF ( pw%in_space == RECIPROCALSPACE .AND. &
         pw%in_use == COMPLEXDATA1D ) THEN

       cnt = SIZE ( pw%cc )

       IF ( i == j ) THEN
!$omp parallel do default(none) private(ig,gg) shared(i,o3,pw,pwdr2,cnt)
          DO ig = 1, cnt
             gg = pw%pw_grid%g(i,ig)**2 - o3*pw%pw_grid%gsq(ig)
             pwdr2%cc(ig) = gg * pw%cc(ig)
          END DO
          flop = flop + 5 * cnt
       ELSE
!$omp parallel do private (ig) default(none) shared(pwdr2,pw,i,j,cnt)
          DO ig = 1, cnt
             pwdr2%cc(ig) = pw%cc(ig) * (pw%pw_grid%g(i,ig)*pw%pw_grid%g(j,ig))
          END DO
          flop = flop + 4 * cnt
       END IF

    ELSE

       CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")

    END IF

    flop = flop * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_dr2

! *****************************************************************************
!> \brief Calculate the tensorial 2nd derivative of a plane wave vector
!>      and divides by |G|^2. pwdr2_gg(G=0) is put to zero.
!> \param pw ...
!> \param pwdr2_gg ...
!> \param i ...
!> \param j ...
!> \param error ...
!> \par History
!>      none
!> \author RD (20-Nov-2006)
!> \note
!>      PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
!>      Adapted from pw_dr2
! *****************************************************************************
  SUBROUTINE pw_dr2_gg ( pw, pwdr2_gg, i, j, error)

    TYPE(pw_type), INTENT(IN)                :: pw
    TYPE(pw_type), INTENT(INOUT)             :: pwdr2_gg
    INTEGER, INTENT(IN)                      :: i, j
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: cnt, handle, ig
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: flop, gg, o3
    CHARACTER(len=*), PARAMETER :: routineN = 'pw_dr2_gg', &
      routineP = moduleN//':'//routineN

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)

    flop = 0.0_dp
    o3 = 1._dp/3._dp

    IF ( pw%in_space == RECIPROCALSPACE .AND. &
         pw%in_use == COMPLEXDATA1D ) THEN

       cnt = SIZE ( pw%cc )

       IF ( i == j ) THEN
!$omp parallel do private (ig) default(none) private(gg) shared(cnt,pw,o3,pwdr2_gg,i)
          DO ig = pw%pw_grid%first_gne0, cnt
             gg = pw%pw_grid%g(i,ig)**2 - o3*pw%pw_grid%gsq(ig)
             pwdr2_gg%cc(ig) = gg * pw%cc(ig) / pw%pw_grid%gsq(ig)
          END DO
          flop = flop + 6 * cnt
       ELSE
!$omp parallel do private (ig) default(none) shared(pwdr2_gg,pw,i,j,cnt)
          DO ig = pw%pw_grid%first_gne0, cnt
             pwdr2_gg%cc(ig) = pw%cc(ig) * (pw%pw_grid%g(i,ig)*pw%pw_grid%g(j,ig)) &
                  / pw%pw_grid%gsq(ig)
          END DO
          flop = flop + 5 * cnt
       END IF

       IF ( pw%pw_grid%have_g0 ) pwdr2_gg%cc ( 1 ) = 0.0_dp

    ELSE

       CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")

    END IF

    flop = flop * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_dr2_gg

! *****************************************************************************
!> \brief Multiplies a G-space function with a smoothing factor of the form
!>      f(|G|) = exp((ecut - G^2)/sigma)/(1+exp((ecut - G^2)/sigma))
!> \param pw ...
!> \param ecut ...
!> \param sigma ...
!> \param error ...
!> \par History
!>      none
!> \author JGH (09-June-2006)
!> \note
!>      PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
! *****************************************************************************
  SUBROUTINE pw_smoothing ( pw, ecut, sigma, error)

    TYPE(pw_type), INTENT(INOUT)             :: pw
    REAL(KIND=dp), INTENT(IN)                :: ecut, sigma
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_smoothing', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cnt, handle, ig
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, f, flop

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)

    flop = 0.0_dp

    IF ( pw%in_space == RECIPROCALSPACE .AND. &
         pw%in_use == COMPLEXDATA1D ) THEN

       cnt = SIZE ( pw%cc )

!$omp parallel do default(none) private(ig,arg,f) shared(sigma,ecut,pw,cnt)
       DO ig = 1, cnt
          arg = (ecut - pw%pw_grid%gsq(ig))/sigma
          f = EXP(arg)/(1+EXP(arg))
          pw%cc(ig) = f * pw%cc(ig)
       END DO
       flop = flop + 6 * cnt

    ELSE

       CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")

    END IF

    flop = flop * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_smoothing

! *****************************************************************************
!> \brief Generalize copy of pw types
!> \param pw1 ...
!> \param pw2 ...
!> \param debug ...
!> \param error ...
!> \par History
!>      JGH (13-Mar-2001) : added gather/scatter cases
!> \author JGH (25-Feb-2001)
!> \note
!>      Copy routine that allows for in_space changes
! *****************************************************************************
  SUBROUTINE pw_transfer ( pw1, pw2, debug, error)

    TYPE(pw_type), INTENT(IN), TARGET        :: pw1
    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw2
    LOGICAL, INTENT(IN), OPTIONAL            :: debug
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_transfer', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPrecondition(pw1%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(pw2%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( pw1%in_space == REALSPACE .AND. pw2%in_space == REALSPACE ) THEN

       ! simple copy should do
       CALL pw_copy ( pw1, pw2, error=error)

    ELSEIF ( pw1%in_space == RECIPROCALSPACE .AND. &
         pw2%in_space == RECIPROCALSPACE ) THEN

       IF ( pw1%in_use == pw2%in_use ) THEN

          ! simple copy should do
          CALL pw_copy ( pw1, pw2, error=error)

       ELSE

          ! we have to gather/scatter the data
          IF ( pw1%in_use == COMPLEXDATA1D ) THEN
             CALL pw_scatter ( pw1, pw2%cc3d, error=error)
          ELSEIF ( pw2%in_use == COMPLEXDATA1D ) THEN
             CALL pw_gather ( pw2, pw1%cc3d, error=error)
          ELSE
             CALL stop_program(routineN,moduleN,__LINE__, "Do not know what to do" )
          END IF

       END IF

    ELSE

       ! FFT needed, all further tests done in fft_wrap_pw1pw2
       CALL fft_wrap_pw1pw2 ( pw1, pw2, debug, error=error)

    END IF

  END SUBROUTINE pw_transfer

! *****************************************************************************
!> \brief pw2 = alpha*pw1 + pw2
!>      alpha defaults to 1
!> \param pw1 ...
!> \param pw2 ...
!> \param alpha ...
!> \param error ...
!> \par History
!>      JGH (21-Feb-2003) : added reference grid functionality
!>      JGH (01-Dec-2007) : rename and remove complex alpha
!> \author apsi
!> \note
!>      Currently only summing up of respective types allowed,
!>      in order to avoid errors
! *****************************************************************************
  SUBROUTINE pw_axpy ( pw1, pw2, alpha, error)

    TYPE(pw_type), INTENT(IN)                :: pw1
    TYPE(pw_type), INTENT(INOUT)             :: pw2
    REAL(KIND=dp), INTENT(in), OPTIONAL      :: alpha
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_axpy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ng, ng1, ng2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: flop, my_alpha

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw1%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(pw2%ref_count>0,cp_failure_level,routineP,error,failure)

    my_alpha=1.0_dp
    IF (PRESENT(alpha)) my_alpha=alpha

    IF ( pw1%pw_grid %id_nr == pw2%pw_grid %id_nr ) THEN

       IF ( pw1%in_use == REALDATA1D .AND. pw2%in_use == REALDATA1D ) THEN
          IF (my_alpha==1.0_dp) THEN
!$omp parallel do private(i) default(none) shared(pw2,pw1)
             DO i = 1, SIZE ( pw2%cr )
                pw2%cr(i) = pw2%cr(i) + pw1%cr(i)
             END DO
             flop = REAL ( SIZE ( pw2%cr ),KIND=dp)
          ELSE
!$omp parallel do private(i) default(none) shared(pw2,pw1,my_alpha)
             DO i = 1, SIZE ( pw2%cr )
                pw2%cr(i) = pw2%cr(i) + my_alpha*pw1%cr(i)
             END DO
             flop = REAL ( 2*SIZE ( pw2%cr ),KIND=dp)
          END IF
       ELSE IF ( pw1%in_use == COMPLEXDATA1D .AND. &
            pw2%in_use == COMPLEXDATA1D ) THEN
          IF (my_alpha==1.0_dp) THEN
!$omp parallel do private(i) default(none) shared(pw2,pw1)
             DO i = 1, SIZE ( pw2%cc )
                pw2%cc(i) = pw2%cc(i) + pw1%cc(i)
             END DO
             flop = REAL ( 2 * SIZE ( pw2%cc ),KIND=dp)
          ELSE
!$omp parallel do private(i) default(none) shared(pw2,my_alpha,pw1)
             DO i = 1, SIZE ( pw2%cc )
                pw2%cc(i) = pw2%cc(i) + my_alpha*pw1%cc(i)
             END DO
             flop = REAL ( 4 * SIZE ( pw2%cc ),KIND=dp)
          END IF
       ELSE IF ( pw1%in_use == REALDATA3D .AND. pw2%in_use == REALDATA3D ) THEN
          IF (my_alpha==1.0_dp) THEN
             pw2%cr3d = pw2%cr3d + pw1%cr3d
             flop = REAL ( SIZE ( pw2%cr3d ),KIND=dp)
          ELSE
             pw2%cr3d = pw2%cr3d + my_alpha * pw1%cr3d
             flop = REAL ( 2*SIZE ( pw2%cr3d ),KIND=dp)
          END IF
       ELSE IF ( pw1%in_use == COMPLEXDATA3D .AND. &
            pw2%in_use == COMPLEXDATA3D ) THEN
          IF (my_alpha==1.0_dp) THEN
             pw2%cc3d = pw2%cc3d + pw1%cc3d
             flop = REAL ( 2 * SIZE ( pw2%cc3d ),KIND=dp)
          ELSE
             pw2%cc3d = pw2%cc3d + my_alpha*pw1%cc3d
             flop = REAL ( 4 * SIZE ( pw2%cc3d ),KIND=dp)
          END IF
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")
       END IF

    ELSE IF ( pw_compatible ( pw1%pw_grid, pw2%pw_grid, error=error) ) THEN

       IF ( pw1%in_use == COMPLEXDATA1D .AND. &
            pw2%in_use == COMPLEXDATA1D .AND. &
            pw1%in_space == RECIPROCALSPACE .AND. &
            pw1%in_space == RECIPROCALSPACE ) THEN

          ng1 = SIZE ( pw1%cc )
          ng2 = SIZE ( pw2%cc )
          ng  = MIN ( ng1, ng2 )
          flop = REAL ( 2 * ng,KIND=dp)

          IF ( pw1%pw_grid%spherical ) THEN

             IF (my_alpha==1.0_dp) THEN
!$omp parallel do private(i) default(none) shared(pw2,pw1,ng)
                DO i = 1, ng
                   pw2%cc (i) = pw2%cc (i) + pw1%cc (i)
                END DO
             ELSE
!$omp parallel do private(i) default(none) shared(pw2,pw1,my_alpha,ng)
                DO i = 1, ng
                   pw2%cc (i) = pw2%cc (i) + my_alpha*pw1%cc (i)
                END DO
             END IF

          ELSEIF ( ( pw1%pw_grid %id_nr == pw2%pw_grid%reference ) ) THEN

             IF (my_alpha==1.0_dp) THEN
                IF( ng1 >= ng2 ) THEN
!$omp parallel do private(i,j) default(none) shared(pw2,pw1,ng2)
                   DO i = 1, ng2
                      j = pw2%pw_grid%gidx ( i )
                      pw2%cc ( i ) = pw2%cc ( i ) + pw1%cc ( j )
                   END DO
                ELSE
!$omp parallel do private(i,j) default(none) shared(pw2,pw1,ng1)
                   DO i = 1, ng1
                      j = pw2%pw_grid%gidx ( i )
                      pw2%cc ( j ) = pw2%cc ( j ) + pw1%cc ( i )
                   END DO
                END IF
             ELSE
                IF( ng1 >= ng2 ) THEN
!$omp parallel do private(i,j) default(none) shared(pw2,pw1,my_alpha,ng2)
                   DO i = 1, ng2
                      j = pw2%pw_grid%gidx ( i )
                      pw2%cc ( i ) = pw2%cc ( i ) + my_alpha*pw1%cc ( j )
                   END DO
                ELSE
!$omp parallel do private(i,j) default(none) shared(pw2,pw1,my_alpha,ng1)
                   DO i = 1, ng1
                      j = pw2%pw_grid%gidx ( i )
                      pw2%cc ( j ) = pw2%cc ( j ) + my_alpha*pw1%cc ( i )
                   END DO
                END IF
             END IF
          ELSEIF ( ( pw2%pw_grid %id_nr == pw1%pw_grid%reference ) ) THEN

             IF (my_alpha==1.0_dp) THEN
                IF( ng1 >= ng2 ) THEN
!$omp parallel do private(i,j) default(none) shared(pw2,pw1,ng2)
                   DO i = 1, ng2
                      j = pw1%pw_grid%gidx ( i )
                      pw2%cc ( i ) = pw2%cc ( i ) + pw1%cc ( j )
                   END DO
                ELSE
!$omp parallel do private(i,j) default(none) shared(pw2,pw1,ng1)
                   DO i = 1, ng1
                      j = pw1%pw_grid%gidx ( i )
                      pw2%cc ( j ) = pw2%cc ( j ) + pw1%cc ( i )
                   END DO
                END IF
             ELSE
                IF( ng1 >= ng2 ) THEN
!$omp parallel do private(i,j) default(none) shared(pw2,pw1,my_alpha,ng2)
                   DO i = 1, ng2
                      j = pw1%pw_grid%gidx ( i )
                      pw2%cc ( i ) = pw2%cc ( i ) + my_alpha* pw1%cc ( j )
                   END DO
                ELSE
!$omp parallel do private(i,j) default(none) shared(pw2,pw1,my_alpha,ng1)
                   DO i = 1, ng1
                      j = pw1%pw_grid%gidx ( i )
                      pw2%cc ( j ) = pw2%cc ( j ) + my_alpha*pw1%cc ( i )
                   END DO
                END IF
             END IF
          ELSE

             WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                        " grid 1 :",pw1%pw_grid %id_nr, &
                        " sperical :",pw1%pw_grid%spherical, &
                        " reference :",pw1%pw_grid%reference
             WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                        " grid 2 :",pw2%pw_grid %id_nr, &
                        " sperical :",pw2%pw_grid%spherical, &
                        " reference :",pw2%pw_grid%reference
             CALL stop_program(routineN,moduleN,__LINE__,"Grids not compatible")

          END IF

       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")
       END IF

    ELSE

       WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                  " grid 1 :",pw1%pw_grid %id_nr, &
                  " sperical :",pw1%pw_grid%spherical, &
                  " reference :",pw1%pw_grid%reference
       WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                  " grid 2 :",pw2%pw_grid %id_nr, &
                  " sperical :",pw2%pw_grid%spherical, &
                  " reference :",pw2%pw_grid%reference
       CALL stop_program(routineN,moduleN,__LINE__,"Grids not compatible")

    END IF

    flop = flop * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_axpy

! *****************************************************************************
!> \brief Gathers the pw vector from a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \param error ...
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE pw_gather_s ( pw, c, scale, error)

    TYPE(pw_type), INTENT(INOUT)             :: pw
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_gather_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: gpt, handle, l, m, n, ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cpy

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( pw%in_use /= COMPLEXDATA1D ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Data field has to be COMPLEXDATA1D")
    ENDIF

    ! after the gather we are in g-space
    pw%in_space = RECIPROCALSPACE

    mapl => pw%pw_grid%mapl%pos
    mapm => pw%pw_grid%mapm%pos
    mapn => pw%pw_grid%mapn%pos

    ngpts = SIZE ( pw%pw_grid%gsq  )

    ghat => pw%pw_grid%g_hat

    IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n) default(none) shared(ngpts, mapl,mapm,mapn,ghat,scale,pw,c)
       DO gpt = 1, ngpts

          l = mapl ( ghat ( 1, gpt ) ) + 1
          m = mapm ( ghat ( 2, gpt ) ) + 1
          n = mapn ( ghat ( 3, gpt ) ) + 1
          pw%cc ( gpt ) = scale * c ( l, m, n )

       END DO

    ELSE

!$omp parallel do private(gpt,l,m,n) default(none) shared(ngpts, mapl,mapm,mapn,ghat,pw,c)
       DO gpt = 1, ngpts

          l = mapl ( ghat ( 1, gpt ) ) + 1
          m = mapm ( ghat ( 2, gpt ) ) + 1
          n = mapn ( ghat ( 3, gpt ) ) + 1
          pw%cc ( gpt ) = c ( l, m, n )

       END DO

    END IF

    cpy = REAL ( ngpts,KIND=dp) * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_gather_s

! *****************************************************************************
!> \brief ...
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pw_gather_p ( pw, c, scale, error)

    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_gather_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: gpt, handle, l, m, mn, n, &
                                                ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat, yzq
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cpy

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( pw%in_use /= COMPLEXDATA1D ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Data field has to be COMPLEXDATA1D")
    ENDIF

    IF ( pw%pw_grid%para%mode /= PW_MODE_DISTRIBUTED ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"This grid type is not distributed")
    ENDIF

    ! after the gather we are in g-space
    pw%in_space = RECIPROCALSPACE

    mapl => pw%pw_grid%mapl%pos
    mapm => pw%pw_grid%mapm%pos
    mapn => pw%pw_grid%mapn%pos

    ngpts = SIZE ( pw%pw_grid%gsq  )

    ghat => pw%pw_grid%g_hat
    yzq => pw%pw_grid%para%yzq

    IF ( PRESENT ( scale ) ) THEN
!$omp parallel do default(none), &
!$omp             private(l,m,n,mn), &
!$omp             shared(ngpts,mapl,mapm,mapn,ghat,yzq,pw,c,scale)
       DO gpt = 1, ngpts

          l = mapl ( ghat ( 1, gpt ) ) + 1
          m = mapm ( ghat ( 2, gpt ) ) + 1
          n = mapn ( ghat ( 3, gpt ) ) + 1
          mn = yzq ( m, n )
          pw%cc ( gpt ) = scale * c ( l, mn )

       END DO
!$omp end parallel do
    ELSE
!$omp parallel do default(none), &
!$omp             private(l,m,n,mn), &
!$omp             shared(ngpts,mapl,mapm,mapn,ghat,yzq,pw,c)
       DO gpt = 1, ngpts

          l = mapl ( ghat ( 1, gpt ) ) + 1
          m = mapm ( ghat ( 2, gpt ) ) + 1
          n = mapn ( ghat ( 3, gpt ) ) + 1
          mn = yzq ( m, n )
          pw%cc ( gpt ) = c ( l, mn )

       END DO
!$omp end parallel do
    END IF

    cpy = REAL ( ngpts,KIND=dp) * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_gather_p

! *****************************************************************************
!> \brief Scatters a pw vector to a 3d data field
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \param error ...
!> \par History
!>      none
!> \author JGH
! *****************************************************************************
  SUBROUTINE pw_scatter_s ( pw, c, scale, error)

    TYPE(pw_type), INTENT(IN)                :: pw
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_scatter_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: gpt, handle, l, m, n, ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cpy

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( pw%in_use /= COMPLEXDATA1D ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Data field has to be COMPLEXDATA1D")
    ENDIF

    IF ( pw%in_space /= RECIPROCALSPACE ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Data has to be in RECIPROCALSPACE")
    ENDIF

    mapl => pw%pw_grid%mapl%pos
    mapm => pw%pw_grid%mapm%pos
    mapn => pw%pw_grid%mapn%pos

    ghat => pw%pw_grid%g_hat

    ngpts = SIZE ( pw%pw_grid%gsq  )

    ! should only zero the unused bits (but the zero is needed)
    IF ( .NOT. PRESENT ( scale ) ) c = 0.0_dp

    IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n) default(none) shared(ngpts, mapl,mapm,mapn,ghat,scale,pw,c)
       DO gpt = 1, ngpts

          l = mapl ( ghat ( 1, gpt ) ) + 1
          m = mapm ( ghat ( 2, gpt ) ) + 1
          n = mapn ( ghat ( 3, gpt ) ) + 1
          c ( l, m, n ) = scale * pw%cc ( gpt )

       END DO

    ELSE

!$omp parallel do private(gpt,l,m,n) default(none) shared(ngpts, mapl,mapm,mapn,ghat,pw,c)
       DO gpt = 1, ngpts

          l = mapl ( ghat ( 1, gpt ) ) + 1
          m = mapm ( ghat ( 2, gpt ) ) + 1
          n = mapn ( ghat ( 3, gpt ) ) + 1
          c ( l, m, n ) = pw%cc ( gpt )

       END DO

    END IF

    IF ( pw%pw_grid%grid_span == HALFSPACE ) THEN

       mapl => pw%pw_grid%mapl%neg
       mapm => pw%pw_grid%mapm%neg
       mapn => pw%pw_grid%mapn%neg

       IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n) default(none) shared(ngpts, mapl,mapm,mapn,ghat,scale,pw,c)
          DO gpt = 1, ngpts

             l = mapl ( ghat ( 1, gpt ) ) + 1
             m = mapm ( ghat ( 2, gpt ) ) + 1
             n = mapn ( ghat ( 3, gpt ) ) + 1
             c ( l, m, n ) = scale * CONJG ( pw%cc ( gpt ) )

          END DO

       ELSE

!$omp parallel do private(gpt,l,m,n) default(none) shared(ngpts, mapl,mapm,mapn,ghat,pw,c)
          DO gpt = 1, ngpts

             l = mapl ( ghat ( 1, gpt ) ) + 1
             m = mapm ( ghat ( 2, gpt ) ) + 1
             n = mapn ( ghat ( 3, gpt ) ) + 1
             c ( l, m, n ) = CONJG ( pw%cc ( gpt ) )

          END DO

       END IF

    END IF

    cpy = REAL ( ngpts,KIND=dp) * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_scatter_s

! *****************************************************************************
!> \brief ...
!> \param pw ...
!> \param c ...
!> \param scale ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pw_scatter_p ( pw, c, scale, error)

    TYPE(pw_type), INTENT(IN), TARGET        :: pw
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: c
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_scatter_p', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: gpt, handle, l, m, mn, n, &
                                                ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat, yzq
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: cpy

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(pw%ref_count>0,cp_failure_level,routineP,error,failure)

    IF ( pw%in_use /= COMPLEXDATA1D ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Data field has to be COMPLEXDATA1D")
    ENDIF

    IF ( pw%in_space /= RECIPROCALSPACE ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Data has to be in RECIPROCALSPACE")
    ENDIF

    IF ( pw%pw_grid%para%mode /= PW_MODE_DISTRIBUTED ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"This grid type is not distributed")
    ENDIF

    mapl => pw%pw_grid%mapl%pos
    mapm => pw%pw_grid%mapm%pos
    mapn => pw%pw_grid%mapn%pos

    ghat => pw%pw_grid%g_hat
    yzq => pw%pw_grid%para%yzq

    ngpts = SIZE ( pw%pw_grid%gsq  )

    IF ( .NOT. PRESENT ( scale ) ) CALL zero_c(SIZE(c),c)

    IF ( PRESENT ( scale ) ) THEN
!$omp parallel do default(none), &
!$omp             private(l,m,n,mn), &
!$omp             shared(ngpts,mapl,mapm,mapn,ghat,yzq,pw,c,scale)
       DO gpt = 1, ngpts

          l = mapl ( ghat ( 1, gpt ) ) + 1
          m = mapm ( ghat ( 2, gpt ) ) + 1
          n = mapn ( ghat ( 3, gpt ) ) + 1
          mn = yzq ( m, n )
          c ( l, mn ) = scale * pw%cc ( gpt )

       END DO
!$omp end parallel do
    ELSE
!$omp parallel do default(none), &
!$omp             private(l,m,n,mn), &
!$omp             shared(ngpts,mapl,mapm,mapn,ghat,yzq,pw,c)
       DO gpt = 1, ngpts

          l = mapl ( ghat ( 1, gpt ) ) + 1
          m = mapm ( ghat ( 2, gpt ) ) + 1
          n = mapn ( ghat ( 3, gpt ) ) + 1
          mn = yzq ( m, n )
          c ( l, mn ) = pw%cc ( gpt )

       END DO
!$omp end parallel do
    END IF

    IF ( pw%pw_grid%grid_span == HALFSPACE ) THEN

       mapm => pw%pw_grid%mapm%neg
       mapn => pw%pw_grid%mapn%neg
       mapl => pw%pw_grid%mapl%neg

       IF ( PRESENT ( scale ) ) THEN
!$omp parallel do default(none), &
!$omp             private(l,m,n,mn), &
!$omp             shared(ngpts,mapl,mapm,mapn,ghat,yzq,pw,c,scale)
          DO gpt = 1, ngpts

             l = mapl ( ghat ( 1, gpt ) ) + 1
             m = mapm ( ghat ( 2, gpt ) ) + 1
             n = mapn ( ghat ( 3, gpt ) ) + 1
             mn = yzq ( m, n )
             c ( l, mn ) = scale * CONJG ( pw%cc ( gpt ) )

          END DO
!$omp end parallel do
       ELSE
!$omp parallel do default(none), &
!$omp             private(l,m,n,mn), &
!$omp             shared(ngpts,mapl,mapm,mapn,ghat,yzq,pw,c)
          DO gpt = 1, ngpts

             l = mapl ( ghat ( 1, gpt ) ) + 1
             m = mapm ( ghat ( 2, gpt ) ) + 1
             n = mapn ( ghat ( 3, gpt ) ) + 1
             mn = yzq ( m, n )
             c ( l, mn ) = CONJG ( pw%cc ( gpt ) )

          END DO
!$omp end parallel do
       END IF

    END IF

    cpy = REAL ( ngpts,KIND=dp) * 1.e-6_dp
    CALL timestop(handle)


  END SUBROUTINE pw_scatter_p

! *****************************************************************************
!> \brief Generic function for 3d FFT of a coefficient_type or pw_type
!> \param pw1 ...
!> \param pw2 ...
!> \param debug ...
!> \param error ...
!> \par History
!>      JGH (30-12-2000): New setup of functions and adaptation to parallelism
!>      JGH (04-01-2001): Moved routine from pws to this module, only covers
!>                        pw_types, no more coefficient types
!> \author apsi
!> \note
!>       fft_wrap_pw1pw2
! *****************************************************************************
  SUBROUTINE fft_wrap_pw1pw2 ( pw1, pw2, debug, error)

    TYPE(pw_type), INTENT(IN), TARGET        :: pw1
    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw2
    LOGICAL, INTENT(IN), OPTIONAL            :: debug
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fft_wrap_pw1pw2', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=9)                         :: mode
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      POINTER                                :: grays
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: c_in, c_out
    INTEGER                                  :: dir, handle, handle2, my_pos, &
                                                nrays, nsize, out_space, stat
    INTEGER, DIMENSION(3)                    :: nloc
    INTEGER, DIMENSION(:), POINTER           :: n
    LOGICAL                                  :: failure, test, use_pw_cuda
    REAL(KIND=dp)                            :: norm

    failure = .FALSE.
    CALL timeset(routineN,handle2)
    CALL timeset(routineN//"_"//TRIM(ADJUSTL(cp_to_string( &
         CEILING (pw1%pw_grid%cutoff/10)*10))),handle)

    CPPrecondition(pw1%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(pw2%ref_count>0,cp_failure_level,routineP,error,failure)
    NULLIFY ( c_in )
    NULLIFY ( c_out )

    IF ( PRESENT ( debug ) ) THEN
       test = debug
    ELSE
       test = .FALSE.
    END IF

    !..check if grids are compatible
    IF ( pw1%pw_grid %id_nr /= pw2%pw_grid %id_nr ) THEN
       IF ( pw1%pw_grid%dvol /= pw2%pw_grid%dvol ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "PW grids not compatible")
       END IF
       IF ( pw1%pw_grid%para%group /= pw2%pw_grid%para%group) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "PW grids have not compatible MPI groups")
       END IF
    END IF

    !..prepare input
    IF ( pw1%in_space == REALSPACE ) THEN
       dir = FWFFT
       norm = 1.0_dp / pw1%pw_grid%ngpts
       out_space = RECIPROCALSPACE
    ELSE IF ( pw1%in_space == RECIPROCALSPACE ) THEN
       dir = BWFFT
       norm = 1.0_dp
       out_space = REALSPACE
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"Error in space tag")
    END IF

    n => pw1%pw_grid%npts

    mode = fftselect ( pw1%in_use, pw2%in_use, pw1%in_space, error=error)

    IF ( pw1%pw_grid%para%mode == PW_MODE_LOCAL ) THEN

       !
       !..replicated data, use local FFT
       !

       IF ( test ) THEN
          WRITE ( *,'(A)') " FFT Protocol "
          IF ( dir == FWFFT ) WRITE ( *,'(A,T76,A)') "  Transform direction ","FWFFT"
          IF ( dir == BWFFT ) WRITE ( *,'(A,T76,A)') "  Transform direction ","BWFFT"
          IF ( pw1%in_space == REALSPACE ) &
               WRITE ( *,'(A,T72,A)') "  in space ","REALSPACE"
          IF ( pw1%in_space == RECIPROCALSPACE ) &
               WRITE ( *,'(A,T66,A)') "  in space ","RECIPROCALSPACE"
          IF ( out_space == REALSPACE ) &
               WRITE ( *,'(A,T72,A)') "  out space ","REALSPACE"
          IF ( out_space == RECIPROCALSPACE ) &
               WRITE ( *,'(A,T66,A)') "  out space ","RECIPROCALSPACE"
          WRITE ( *,'(A,T66,E15.6)') "  scale factor",norm
       END IF

       SELECT CASE ( mode )
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Illegal combination of in_use and in_space")
       CASE ( "FW_C3DC3D" )
          c_in => pw1%cc3d
          c_out => pw2%cc3d
          CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug = test )
       CASE ( "FW_R3DC3D" )
          pw2%cc3d = CMPLX ( pw1%cr3d, 0.0_dp, KIND=dp)
          c_out => pw2%cc3d
          CALL fft3d ( dir, n, c_out, scale = norm, debug = test )
       CASE ( "FW_C3DC1D" )
          c_in => pw1%cc3d
          ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = stat)
          IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "c_out",2*dp_size*PRODUCT(n(1:3)))
          ! transform
          CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug = test )
          ! gather results
          IF ( test ) WRITE ( *,'(A)') "  PW_GATHER : 3d -> 1d "
          CALL pw_gather ( pw2, c_out, error = error)
          DEALLOCATE ( c_out, STAT = stat )
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CASE ( "FW_R3DC1D" )
#if defined (__PW_CUDA)
          CALL pw_cuda_r3dc1d_3d(pw1, pw2, scale = norm, error = error)
#else
          ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = stat )
          IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "c_out",2*dp_size*PRODUCT(n(1:3)))
          nsize=SIZE(pw1%cr3d,1)*SIZE(pw1%cr3d,2)*SIZE(pw1%cr3d,3)
          CALL copy_rc(nsize,pw1%cr3d,c_out)
          CALL fft3d ( dir, n, c_out, scale = norm, debug = test )
          CALL pw_gather ( pw2, c_out, error = error)
          DEALLOCATE ( c_out, STAT = stat )
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
#endif
       CASE ( "BW_C3DC3D" )
          c_in => pw1%cc3d
          c_out => pw2%cc3d
          CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug=test )
       CASE ( "BW_C3DR3D" )
          c_in => pw1%cc3d
          ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = stat )
          IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "c_out",2*dp_size*PRODUCT(n(1:3)))
          CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug = test )
          ! use real part only
          IF ( test ) WRITE ( *,'(A)') "  REAL part "
          pw2%cr3d = REAL ( c_out,KIND=dp)
          DEALLOCATE ( c_out, STAT = stat )
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CASE ( "BW_C1DC3D" )
          c_out => pw2%cc3d
          IF ( test ) WRITE ( *,'(A)') "  PW_SCATTER : 3d -> 1d "
          CALL pw_scatter ( pw1, c_out, error=error)
          CALL fft3d ( dir, n, c_out, scale = norm, debug = test )
       CASE ( "BW_C1DR3D" )
#if defined (__PW_CUDA)
          CALL pw_cuda_c1dr3d_3d(pw1, pw2, scale = norm, error = error)
#else
          ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = stat)
          IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "c_out",2*dp_size*PRODUCT(n(1:3)))
          IF ( test ) WRITE ( *,'(A)') "  PW_SCATTER : 3d -> 1d "
          CALL pw_scatter ( pw1, c_out, error = error)
          ! transform
          CALL fft3d ( dir, n, c_out, scale = norm, debug = test )
          ! use real part only
          IF ( test ) WRITE ( *,'(A)') "  REAL part "
          nsize=SIZE(pw2%cr3d,1)*SIZE(pw2%cr3d,2)*SIZE(pw2%cr3d,3)
          CALL copy_cr(nsize,c_out,pw2%cr3d)
          DEALLOCATE ( c_out, STAT = stat )
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
#endif
       END SELECT

       IF ( test ) WRITE ( *,'(A)') " End of FFT Protocol "

    ELSE

       !
       !..parallel FFT
       !

       IF ( test .AND. pw1%pw_grid%para%group_head ) THEN
          WRITE ( *,'(A)') " FFT Protocol "
          IF ( dir == FWFFT ) WRITE ( *,'(A,T76,A)') "  Transform direction ","FWFFT"
          IF ( pw1%in_space == REALSPACE ) &
               WRITE ( *,'(A,T72,A)') "  in space ","REALSPACE"
          IF ( pw1%in_space == RECIPROCALSPACE ) &
               WRITE ( *,'(A,T66,A)') "  in space ","RECIPROCALSPACE"
          IF ( out_space == REALSPACE ) &
               WRITE ( *,'(A,T72,A)') "  out space ","REALSPACE"
          IF ( out_space == RECIPROCALSPACE ) &
               WRITE ( *,'(A,T66,A)') "  out space ","RECIPROCALSPACE"
          WRITE ( *,'(A,T66,E15.6)') "  scale factor",norm
       END IF

       my_pos = pw1%pw_grid%para%my_pos
       nrays = pw1%pw_grid%para%nyzray ( my_pos )
       grays => pw1%pw_grid%grays
       CPPostcondition(SIZE(grays,1)==n(1),cp_failure_level,routineP,error,failure)
       CPPostcondition(SIZE(grays,2)==nrays,cp_failure_level,routineP,error,failure)

       SELECT CASE ( mode )
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Illegal combination of in_use and in_space "//&
                            "in parallel 3d FFT")
       CASE ( "FW_C3DC1D" )
          !..prepare input
          c_in => pw1%cc3d
          CALL zero_c(SIZE(grays),grays)
          !..transform
          IF ( pw1%pw_grid%para%ray_distribution ) THEN
             CALL fft3d ( dir, n, c_in, grays, pw1%pw_grid%para%group, &
                  pw1%pw_grid%para%rs_group, &
                  pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                  pw1%pw_grid%para%bo, scale = norm, debug=test )
          ELSE
             CALL fft3d ( dir, n, c_in, grays, pw1%pw_grid%para%rs_group, &
                  pw1%pw_grid%para%bo, scale = norm, debug=test )
          END IF
          !..prepare output
          IF ( test .AND. pw1%pw_grid%para%group_head ) &
               WRITE ( *,'(A)') "  PW_GATHER : 2d -> 1d "
          CALL pw_gather ( pw2, grays, error=error)
       CASE ( "FW_R3DC1D" )
#if defined (__PW_CUDA)
         ! (no ray dist. is not efficient in CUDA)
          use_pw_cuda = pw1%pw_grid%para%ray_distribution
#else
          use_pw_cuda = .FALSE.
#endif
          IF (use_pw_cuda) THEN
             CALL pw_cuda_r3dc1d_3d_ps(pw1, pw2, scale = norm, error = error)
          ELSE
             !..prepare input
             nloc = pw1%pw_grid%npts_local
             ALLOCATE ( c_in( nloc(1), nloc(2), nloc(3) ), STAT = stat)
             IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                             "c_in",2*dp_size*PRODUCT(nloc(1:3)))
             nsize=nloc(1)*nloc(2)*nloc(3)
             CALL copy_rc(nsize,pw1%cr3d,c_in)
             CALL zero_c(SIZE(grays),grays)
             !..transform
             IF ( pw1%pw_grid%para%ray_distribution ) THEN
                CALL fft3d ( dir, n, c_in, grays, pw1%pw_grid%para%group, &
                     pw1%pw_grid%para%rs_group, &
                     pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                     pw1%pw_grid%para%bo, scale = norm, debug=test )
             ELSE
                CALL fft3d ( dir, n, c_in, grays, pw1%pw_grid%para%rs_group, &
                     pw1%pw_grid%para%bo, scale = norm, debug=test )
             END IF
             !..prepare output
             IF ( test .AND. pw1%pw_grid%para%group_head ) &
                  WRITE ( *,'(A)') "  PW_GATHER : 2d -> 1d "
             CALL pw_gather ( pw2, grays, error=error)
             DEALLOCATE ( c_in, STAT = stat )
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       CASE ( "BW_C1DC3D" )
          !..prepare input
          IF ( test .AND. pw1%pw_grid%para%group_head ) &
               WRITE ( *,'(A)') "  PW_SCATTER : 2d -> 1d "
          CALL zero_c(SIZE(grays),grays)
          CALL pw_scatter ( pw1, grays, error=error)
          c_in => pw2%cc3d
          !..transform
          IF ( pw1%pw_grid%para%ray_distribution ) THEN
             CALL fft3d ( dir, n, c_in, grays, pw1%pw_grid%para%group, &
                  pw1%pw_grid%para%rs_group, &
                  pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                  pw1%pw_grid%para%bo, scale = norm, debug=test )
          ELSE
             CALL fft3d ( dir, n, c_in, grays, pw1%pw_grid%para%rs_group, &
                  pw1%pw_grid%para%bo, scale = norm, debug=test )
          END IF
          !..prepare output (nothing to do)
       CASE ( "BW_C1DR3D" )
#if defined (__PW_CUDA)
         ! (no ray dist. is not efficient in CUDA)
          use_pw_cuda = pw1%pw_grid%para%ray_distribution
#else
          use_pw_cuda = .FALSE.
#endif
          IF (use_pw_cuda) THEN
             CALL pw_cuda_c1dr3d_3d_ps(pw1, pw2, scale = norm, error = error)
          ELSE
             !..prepare input
             IF ( test .AND. pw1%pw_grid%para%group_head ) &
                  WRITE ( *,'(A)') "  PW_SCATTER : 2d -> 1d "
             CALL zero_c(SIZE(grays),grays)
             CALL pw_scatter ( pw1, grays, error=error)
             nloc = pw2%pw_grid%npts_local
             ALLOCATE ( c_in( nloc(1), nloc(2), nloc(3) ), STAT = stat )
             IF (stat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                             "c_in",2*dp_size*PRODUCT(nloc(1:3)))
             nsize=nloc(1)*nloc(2)*nloc(3)
             !..transform
             IF ( pw1%pw_grid%para%ray_distribution ) THEN
                CALL fft3d ( dir, n, c_in, grays, pw1%pw_grid%para%group, &
                     pw1%pw_grid%para%rs_group, &
                     pw1%pw_grid%para%yzp, pw1%pw_grid%para%nyzray, &
                     pw1%pw_grid%para%bo, scale = norm, debug=test )
             ELSE
                CALL fft3d ( dir, n, c_in, grays, pw1%pw_grid%para%rs_group, &
                     pw1%pw_grid%para%bo, scale = norm, debug=test )
             END IF
             !..prepare output
             IF ( test .AND. pw1%pw_grid%para%group_head ) &
                  WRITE ( *,'(A)') "  Real part "
             CALL copy_cr(nsize,c_in,pw2%cr3d)
             DEALLOCATE ( c_in, STAT = stat )
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
       END SELECT

    END IF

    !..update the space tag for pw2
    pw2%in_space = out_space

    IF ( test .AND. pw1%pw_grid%para%group_head ) THEN
       WRITE ( *,'(A)') " End of FFT Protocol "
    END IF
    CALL timestop(handle)
    CALL timestop(handle2)


  END SUBROUTINE fft_wrap_pw1pw2

! *****************************************************************************
!> \brief ...
!> \param use1 ...
!> \param use2 ...
!> \param space1 ...
!> \param error ...
! *****************************************************************************
  FUNCTION fftselect ( use1, use2, space1, error) RESULT ( mode )
    INTEGER, INTENT(IN)                      :: use1, use2, space1
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=9)                         :: mode

    CHARACTER(len=*), PARAMETER :: routineN = 'fftselect', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF ( space1 == REALSPACE ) THEN
       mode ( 1 : 3 ) = "FW_"
    ELSE IF ( space1 == RECIPROCALSPACE ) THEN
       mode ( 1 : 3 ) = "BW_"
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"Error in space tag")
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
       CALL stop_program(routineN,moduleN,__LINE__,"Error in use1 tag")
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
       CALL stop_program(routineN,moduleN,__LINE__,"Error in use1 tag")
    END SELECT

  END FUNCTION fftselect

! *****************************************************************************
!> \brief writes a small description of the actual grid
!>      (change to output the data as cube file, maybe with an
!>      optional long_description arg?)
!> \param pw the pw data to output
!> \param unit_nr the unit to output to
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
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
          WRITE (unit=unit_nr, fmt="(' cr=<real(',i8,':',i8,')>')")&
               LBOUND(pw%cr,1),UBOUND(pw%cr,1)
       ELSE
          WRITE (unit=unit_nr, fmt="(' cr=*null*')")
       END IF
    CASE (REALDATA3D)
       WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_use=REALDATA3D"
       IF (ASSOCIATED(pw%cr3d)) THEN
          WRITE (unit=unit_nr, fmt="(' cr3d=<real(',i8,':',i8,',',i8,':',i8,',',i8,':',i8,')>')")&
               LBOUND(pw%cr3d,1),UBOUND(pw%cr3d,1),LBOUND(pw%cr3d,2),UBOUND(pw%cr3d,2),&
               LBOUND(pw%cr3d,3),UBOUND(pw%cr3d,3)
       ELSE
          WRITE (unit=unit_nr, fmt="(' cr3d=*null*')")
       END IF
    CASE (COMPLEXDATA1D)
       WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_use=COMPLEXDATA1D"
       IF (ASSOCIATED(pw%cc)) THEN
          WRITE (unit=unit_nr, fmt="(' cc=<real(',i8,':',i8,')>')")&
               LBOUND(pw%cc,1),UBOUND(pw%cc,1)
       ELSE
          WRITE (unit=unit_nr, fmt="(' cc=*null*')")
       END IF
    CASE (COMPLEXDATA3D)
       WRITE (unit=unit_nr, fmt="(a)",iostat=iostat) " in_use=COMPLEXDATA3D"
       IF (ASSOCIATED(pw%cc3d)) THEN
          WRITE (unit=unit_nr, fmt="(' cc3d=<real(',i8,':',i8,',',i8,':',i8,',',i8,':',i8,')>')")&
               LBOUND(pw%cc3d,1),UBOUND(pw%cc3d,1),LBOUND(pw%cc3d,2),UBOUND(pw%cc3d,2),&
               LBOUND(pw%cc3d,3),UBOUND(pw%cc3d,3)
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

! *****************************************************************************
!> \brief ...
!> \param grida ...
!> \param gridb ...
!> \param error ...
! *****************************************************************************
  FUNCTION pw_compatible ( grida, gridb, error) RESULT ( compat )
    TYPE(pw_grid_type), INTENT(IN)           :: grida, gridb
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL                                  :: compat

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_compatible', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    compat = .FALSE.
    IF ( grida%id_nr == gridb%id_nr ) THEN
       compat = .TRUE.
    ELSE IF ( grida%reference == gridb%id_nr ) THEN
       compat = .TRUE.
    ELSE IF ( gridb%reference == grida%id_nr ) THEN
       compat = .TRUE.
    END IF

  END FUNCTION pw_compatible

! *****************************************************************************
!> \brief Calculate integral over unit cell for functions in plane wave basis
!>      only returns the real part of it ......
!> \param pw1 ...
!> \param pw2 ...
!> \param error ...
!> \par History
!>      JGH (14-Mar-2001) : Parallel sum and some tests, HALFSPACE case
!> \author apsi
! *****************************************************************************
  FUNCTION pw_integral_ab ( pw1, pw2, error) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: integral_value

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_integral_ab', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF ( pw1%pw_grid %id_nr /= pw2%pw_grid %id_nr ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Grids incompatible")
    END IF

    ! since the return value is real, only do accurate sum on the real bit ?
    IF ( pw1%in_use == REALDATA3D .AND. pw2%in_use == REALDATA3D ) THEN
       integral_value = accurate_sum ( pw1%cr3d ( :, :, : )  &
            * pw2%cr3d ( :, :, : ) )
    ELSE IF ( pw1%in_use == REALDATA3D &
         .AND. pw2%in_use == COMPLEXDATA3D ) THEN
       integral_value = accurate_sum (  pw1%cr3d ( :, :, : )  &
            * pw2%cc3d ( :, :, : ) )
    ELSE IF ( pw1%in_use == COMPLEXDATA3D &
         .AND. pw2%in_use == REALDATA3D ) THEN
       integral_value = accurate_sum ( pw1%cc3d ( :, :, : ) &
            *  pw2%cr3d ( :, :, : ) )
    ELSE IF ( pw1%in_use == COMPLEXDATA3D &
         .AND. pw2%in_use == COMPLEXDATA3D ) THEN
       integral_value = accurate_sum ( CONJG ( pw1%cc3d ( :, :, : ) ) &
            * pw2%cc3d ( :, :, : ) )

    ELSE IF ( pw1%in_use == REALDATA1D &
         .AND. pw2%in_use == REALDATA1D ) THEN
       integral_value = accurate_sum ( pw1%cr (:) *  pw2%cr (:) )
    ELSE IF ( pw1%in_use == REALDATA1D &
         .AND. pw2%in_use == COMPLEXDATA1D ) THEN
       integral_value = accurate_sum ( pw1%cr (:) * pw2%cc (:) )
    ELSE IF ( pw1%in_use == COMPLEXDATA1D &
         .AND. pw2%in_use == REALDATA1D ) THEN
       integral_value = accurate_sum ( pw1%cc (:) * pw2%cr (:) )
    ELSE IF ( pw1%in_use == COMPLEXDATA1D &
         .AND. pw2%in_use == COMPLEXDATA1D ) THEN
       integral_value = accurate_sum ( CONJG ( pw1%cc (:) ) * pw2%cc (:) )
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"No possible DATA")
    END IF

    IF ( pw1%in_use == REALDATA3D .OR. pw1%in_use == COMPLEXDATA3D ) THEN
       integral_value = integral_value * pw1%pw_grid%dvol
    ELSE
       integral_value = integral_value * pw1%pw_grid%vol
    ENDIF
    IF ( pw1%in_use == COMPLEXDATA1D ) THEN
       IF ( pw1%pw_grid%grid_span == HALFSPACE ) THEN
          integral_value = 2.0_dp * integral_value
          IF ( pw1%pw_grid%have_g0 ) integral_value = integral_value - &
               CONJG ( pw1%cc ( 1 ) ) * pw2%cc ( 1 )
       END IF
    END IF

    IF ( pw1%pw_grid%para%mode == PW_MODE_DISTRIBUTED ) &
         CALL mp_sum ( integral_value, pw1%pw_grid%para%group )

  END FUNCTION pw_integral_ab

! *****************************************************************************
!> \brief Calculate integral over unit cell for functions in plane wave basis
!>      only returns the real part of it ......
!> \param pw1 ...
!> \param pw2 ...
!> \param error ...
!> \par History
!>      JGH (14-Mar-2001) : Parallel sum and some tests, HALFSPACE case
!> \author apsi
! *****************************************************************************
  FUNCTION pw_integral_ab_fast ( pw1, pw2, error) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: integral_value

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_integral_ab_fast', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF ( pw1%pw_grid %id_nr /= pw2%pw_grid %id_nr ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Grids incompatible")
    END IF

    ! since the return value is real, only do accurate sum on the real bit ?
    IF ( pw1%in_use == REALDATA3D .AND. pw2%in_use == REALDATA3D ) THEN
       integral_value = SUM ( pw1%cr3d ( :, :, : )  &
            * pw2%cr3d ( :, :, : ) )
    ELSE IF ( pw1%in_use == REALDATA3D &
         .AND. pw2%in_use == COMPLEXDATA3D ) THEN
       integral_value = SUM (  pw1%cr3d ( :, :, : )  &
            * pw2%cc3d ( :, :, : ) )
    ELSE IF ( pw1%in_use == COMPLEXDATA3D &
         .AND. pw2%in_use == REALDATA3D ) THEN
       integral_value = SUM ( pw1%cc3d ( :, :, : ) &
            *  pw2%cr3d ( :, :, : ) )
    ELSE IF ( pw1%in_use == COMPLEXDATA3D &
         .AND. pw2%in_use == COMPLEXDATA3D ) THEN
       integral_value = SUM ( CONJG ( pw1%cc3d ( :, :, : ) ) &
            * pw2%cc3d ( :, :, : ) )

    ELSE IF ( pw1%in_use == REALDATA1D &
         .AND. pw2%in_use == REALDATA1D ) THEN
       integral_value = SUM ( pw1%cr (:) *  pw2%cr (:) )
    ELSE IF ( pw1%in_use == REALDATA1D &
         .AND. pw2%in_use == COMPLEXDATA1D ) THEN
       integral_value = SUM ( pw1%cr (:) * pw2%cc (:) )
    ELSE IF ( pw1%in_use == COMPLEXDATA1D &
         .AND. pw2%in_use == REALDATA1D ) THEN
       integral_value = SUM ( pw1%cc (:) * pw2%cr (:) )
    ELSE IF ( pw1%in_use == COMPLEXDATA1D &
         .AND. pw2%in_use == COMPLEXDATA1D ) THEN
       integral_value = SUM ( CONJG ( pw1%cc (:) ) * pw2%cc (:) )
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"No possible DATA")
    END IF

    IF ( pw1%in_use == REALDATA3D .OR. pw1%in_use == COMPLEXDATA3D ) THEN
       integral_value = integral_value * pw1%pw_grid%dvol
    ELSE
       integral_value = integral_value * pw1%pw_grid%vol
    ENDIF
    IF ( pw1%in_use == COMPLEXDATA1D ) THEN
       IF ( pw1%pw_grid%grid_span == HALFSPACE ) THEN
          integral_value = 2.0_dp * integral_value
          IF ( pw1%pw_grid%have_g0 ) integral_value = integral_value - &
               CONJG ( pw1%cc ( 1 ) ) * pw2%cc ( 1 )
       END IF
    END IF

    IF ( pw1%pw_grid%para%mode == PW_MODE_DISTRIBUTED ) &
         CALL mp_sum ( integral_value, pw1%pw_grid%para%group )

  END FUNCTION pw_integral_ab_fast

! *****************************************************************************
!> \brief ...
!> \param pw1 ...
!> \param flag ...
!> \param error ...
! *****************************************************************************
  FUNCTION pw_integral_aa_fast ( pw1, flag, error) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1
    INTEGER, INTENT(IN), OPTIONAL            :: flag
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: integral_value

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_integral_aa_fast', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF ( PRESENT ( flag ) ) THEN
       IF ( flag == SQUARE ) THEN
          IF ( pw1%in_use == REALDATA3D ) THEN
             integral_value = SUM ( pw1%cr3d ( :, :, : )** 2 )
          ELSE IF ( pw1%in_use == COMPLEXDATA3D ) THEN
             integral_value = SUM ( CONJG ( pw1%cc3d ( :, :, : ) ) &
                  * pw1%cc3d ( :, :, : )  )
          ELSE IF ( pw1%in_use == REALDATA1D ) THEN
             integral_value = SUM ( pw1%cr (:) ** 2 )
          ELSE IF ( pw1%in_use == COMPLEXDATA1D ) THEN
             integral_value = SUM ( CONJG ( pw1%cc (:) ) &
                  * pw1%cc (:) )
             IF ( pw1%pw_grid%grid_span == HALFSPACE ) THEN
                integral_value = 2.0_dp * integral_value
                IF ( pw1%pw_grid%have_g0 ) integral_value = integral_value - &
                     CONJG ( pw1%cc ( 1 ) ) * pw1%cc ( 1 )
             END IF
          ELSE
             CALL stop_program(routineN,moduleN,__LINE__,"No possible SQUARE DATA")
          END IF

       ELSE IF ( flag == SQUAREROOT ) THEN
          CALL stop_program(routineN,moduleN,__LINE__,"No SQUAREROOT defined")
       END IF

    ELSE
       IF ( pw1%in_use == REALDATA3D ) THEN
          integral_value = SUM ( pw1%cr3d ( :, :, : ) )
       ELSE IF ( pw1%in_use == COMPLEXDATA3D ) THEN
          integral_value = SUM ( pw1%cc3d ( :, :, : ) )
       ELSE IF ( pw1%in_use == REALDATA1D ) THEN
          integral_value = SUM ( pw1%cr (:) )
       ELSE IF ( pw1%in_use == COMPLEXDATA1D ) THEN
          integral_value = SUM ( pw1%cc (:) )
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,"No possible DATA")
       END IF
    END IF

    IF ( pw1%in_use == REALDATA3D .OR. pw1%in_use == COMPLEXDATA3D ) THEN
       integral_value = integral_value * pw1%pw_grid%dvol
    ELSE
       integral_value = integral_value * pw1%pw_grid%vol
    END IF

    IF ( pw1%pw_grid%para%mode == PW_MODE_DISTRIBUTED ) &
         CALL mp_sum ( integral_value, pw1%pw_grid%para%group )

  END FUNCTION pw_integral_aa_fast

! *****************************************************************************
!> \brief ...
!> \param pw1 ...
!> \param pw2 ...
!> \param error ...
! *****************************************************************************
  FUNCTION pw_integral_a2b ( pw1, pw2, error) RESULT ( integral_value )

    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: integral_value

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_integral_a2b', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF ( pw1%pw_grid %id_nr /= pw2%pw_grid %id_nr ) THEN
       CALL stop_program(routineN,moduleN,__LINE__,"Grids incompatible")
    END IF
    IF ( pw1%in_use == REALDATA1D .AND. &
         pw2%in_use == REALDATA1D ) THEN
       integral_value = accurate_sum ( pw1%cr (:) * pw2%cr (:) &
            * pw1%pw_grid%gsq (:) )
    ELSE IF ( pw1%in_use == COMPLEXDATA1D .AND. &
         pw2%in_use == COMPLEXDATA1D ) THEN
       integral_value = accurate_sum ( REAL ( CONJG ( pw1%cc (:) ) &
            *  pw2%cc (:) ,KIND=dp) * pw1%pw_grid%gsq (:))
       IF ( pw1%pw_grid%grid_span == HALFSPACE ) THEN
          integral_value = 2.0_dp * integral_value
       END IF
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"No possible DATA")
    END IF

    IF ( pw1%in_use == REALDATA3D .OR. pw1%in_use == COMPLEXDATA3D ) THEN
       integral_value = integral_value * pw1%pw_grid%dvol
    ELSE
       integral_value = integral_value * pw1%pw_grid%vol
    END IF

    IF ( pw1%pw_grid%para%mode == PW_MODE_DISTRIBUTED ) &
         CALL mp_sum ( integral_value, pw1%pw_grid%para%group )

  END FUNCTION pw_integral_a2b

! *****************************************************************************
!> \brief Calculate the structure factor for point r
!> \param sf ...
!> \param r ...
!> \param error ...
!> \par History
!>      none
!> \author JGH (05-May-2006)
!> \note
!>      PW has to be in RECIPROCALSPACE and data in use is COMPLEXDATA1D
! *****************************************************************************
  SUBROUTINE pw_structure_factor ( sf, r, error)

    TYPE(pw_type), INTENT(INOUT)             :: sf
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: r
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_structure_factor', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cnt, handle, ig
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, flop

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CPPrecondition(sf%ref_count>0,cp_failure_level,routineP,error,failure)
    flop = 0.0_dp

    IF ( sf%in_space == RECIPROCALSPACE .AND. &
         sf%in_use == COMPLEXDATA1D ) THEN

       cnt = SIZE ( sf%cc )

!$omp parallel do private (ig,arg) default(none) shared(sf,r,cnt)
       DO ig = 1, cnt
          arg = DOT_PRODUCT(sf%pw_grid%g(:,ig),r)
          sf%cc(ig)=CMPLX(COS(arg),-SIN(arg),KIND=dp)
       END DO
       flop = flop + 7 * cnt
    ELSE

       CALL stop_program(routineN,moduleN,__LINE__,"No suitable data field")

    END IF

    flop = flop * 1.e-6_dp
    CALL timestop(handle)

  END SUBROUTINE pw_structure_factor

! *****************************************************************************
!> \brief ...
!> \param fun ...
!> \param isign ...
!> \param oprt ...
!> \param error ...
!> \retval total_fun ...
! *****************************************************************************
  FUNCTION pw_integrate_function(fun,isign,oprt,error) RESULT(total_fun)

    TYPE(pw_type), INTENT(IN)                :: fun
    INTEGER, INTENT(IN), OPTIONAL            :: isign
    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: oprt
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: total_fun

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_integrate_function', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iop
    LOGICAL                                  :: failure

    failure = .FALSE.
    iop = 0
    IF ( PRESENT(oprt) ) THEN
       SELECT CASE (oprt)
       CASE ("ABS","abs")
          iop = 1
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,"Unknown operator")
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
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "In_space/in_use combination not implemented")
       END IF
    ELSEIF (fun%in_space == RECIPROCALSPACE) THEN
       IF (iop==1) &
            CALL stop_program(routineN,moduleN,__LINE__,&
                              "Operator ABS not implemented")
       IF (fun%in_use == COMPLEXDATA3D) THEN
          total_fun = fun%pw_grid%vol*fun%cc3d(      &
               fun%pw_grid%bounds(1,1), &
               fun%pw_grid%bounds(1,2), &
               fun%pw_grid%bounds(1,3))
       ELSEIF (fun%in_use == COMPLEXDATA1D) THEN
          IF ( fun%pw_grid%have_g0 ) total_fun = fun%pw_grid%vol*fun%cc(1)
       ELSE
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "In_space/in_use combination not implemented")
       END IF
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,"No space defined")
    END IF
    IF (fun%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
       CALL mp_sum(total_fun,fun%pw_grid%para%group)
    END IF
    IF ( PRESENT(isign) ) THEN
       total_fun = total_fun * SIGN(1._dp,REAL(isign,dp))
    END IF

  END FUNCTION pw_integrate_function

END MODULE pw_methods

