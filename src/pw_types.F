!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
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
!!     JGH (26-feb-2003) : OpenMP enabled
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
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE fft_tools,                       ONLY: BWFFT,&
                                             FWFFT,&
                                             fft3d
  USE kinds,                           ONLY: dbl
  USE machine,                         ONLY: m_loc_c,&
                                             m_loc_r
  USE message_passing,                 ONLY: mp_sendrecv
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             PW_MODE_DISTRIBUTED,&
                                             PW_MODE_LOCAL,&
                                             pw_grid_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop,&
                                             print_stack

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_type, pw_p_type
  PUBLIC :: pw_allocate, pw_deallocate, pw_retain, pw_release, pw_create
  PUBLIC :: pw_zero, pw_gather, pw_scatter
  PUBLIC :: pw_copy, pw_add, pw_sumup, pw_subtract, pw_reduce, pw_transfer
  PUBLIC :: pw_derive, pw_derive_fd, pw_sumup_square, pw_fft_wrap, pw_write
  PUBLIC :: pw_prolongate, pw_inject, pw_compare_debug
  PUBLIC :: pw_prolongate_NoPBC, pw_inject_NoPBC
  PUBLIC :: REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D, NODATA
  PUBLIC :: NOSPACE, REALSPACE, RECIPROCALSPACE

  TYPE pw_type
     REAL ( dbl ), DIMENSION ( : ), POINTER :: cr, mass_cr
     REAL ( dbl ), DIMENSION ( :, :, : ), POINTER :: cr3d
     COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: cc, mass_cc
     COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: cc3d

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

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN='pw_types'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.false.

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

!******************************************************************************
!!****** pw_types/pw_allocate [1.0] *
!!
!!   NAME
!!     pw_allocate
!!
!!   SYNOPSIS
!!     Subroutine pw_allocate(pw, pw_grid, use_data, mass_flag, cr3d_ptr,&
!!         error)
!!       Type(pw_type), Intent (INOUT), Target:: pw
!!       Type(pw_grid_type), Intent (IN), Target:: pw_grid
!!       Integer, Intent (IN):: use_data
!!       Logical, Intent (IN), Optional:: mass_flag
!!       Real(Kind=dbl), Dimension(:,:,:), Optional, Pointer:: cr3d_ptr
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine pw_allocate
!!
!!   FUNCTION
!!     Allocate a plane wave type variable
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_allocate ( pw, pw_grid, use_data, mass_flag, cr3d_ptr, error )


    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw
    TYPE(pw_grid_type), INTENT(IN), TARGET   :: pw_grid
    INTEGER, INTENT(IN)                      :: use_data
    LOGICAL, INTENT(IN), OPTIONAL            :: mass_flag
    REAL(kind=dbl), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: cr3d_ptr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_allocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: allocstat, handle
    INTEGER, DIMENSION(:, :), POINTER        :: bounds
    LOGICAL                                  :: mass_flag_local
    TYPE(cp_logger_type), pointer :: logger

!------------------------------------------------------------------------------

  CALL timeset(routineN,"I","",handle)
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
  mass_flag_local = .FALSE.
  IF ( PRESENT ( mass_flag ) ) mass_flag_local = mass_flag

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

  NULLIFY ( pw % cr, pw % cc, pw % cr3d, pw % cc3d ,&
       pw % mass_cr, pw % mass_cc)

  IF ( use_data == REALDATA1D ) THEN
     ALLOCATE ( pw % cr ( pw % pw_grid % ngpts_cut_local ), STAT = allocstat )
     IF ( allocstat /= 0 ) CALL stop_memory ( "pw_allocate", &
        "cr", pw % pw_grid % ngpts_cut_local )
     IF ( mass_flag_local ) THEN
       ALLOCATE ( pw % mass_cr ( pw % pw_grid % ngpts_cut_local ), STAT = allocstat )
       IF ( allocstat /= 0 ) CALL stop_memory ( "pw_allocate", &
          "mass_cr", pw % pw_grid % ngpts_cut_local )
     ENDIF

  ELSE IF ( use_data == COMPLEXDATA1D ) THEN
     ALLOCATE ( pw % cc ( pw % pw_grid % ngpts_cut_local ), STAT = allocstat )
     IF ( allocstat /= 0 ) CALL stop_memory ( "pw_allocate", &
        "cc", pw % pw_grid % ngpts_cut_local )
     IF ( mass_flag_local ) THEN
       ALLOCATE ( pw % mass_cc ( pw % pw_grid % ngpts_cut_local ), STAT = allocstat )
       IF ( allocstat /= 0 ) CALL stop_memory ( "pw_allocate", &
          "mass_cc", pw % pw_grid % ngpts_cut_local )
     ENDIF

  ELSE IF ( use_data == REALDATA3D ) THEN
     IF (PRESENT(cr3d_ptr)) THEN
        IF (ASSOCIATED(cr3d_ptr)) THEN
           CPPreconditionNoFail(ALL(LBOUND(cr3d_ptr)==bounds(1,:)),cp_failure_level,routineP,error)
           CPPreconditionNoFail(ALL(UBOUND(cr3d_ptr)==bounds(2,:)),cp_failure_level,routineP,error)
           pw%cr3d => cr3d_ptr
        END IF
     END IF
     IF (.NOT.ASSOCIATED(pw%cr3d)) THEN
        ALLOCATE ( pw % cr3d ( &
             bounds ( 1, 1 ) : bounds ( 2, 1 ), &
             bounds ( 1, 2 ) : bounds ( 2, 2 ), &
             bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = allocstat )
        IF ( allocstat /= 0 ) CALL stop_memory ( "pw_allocate", "cr3d", 0 )
     END IF

  ELSE IF ( use_data == COMPLEXDATA3D ) THEN
     ALLOCATE ( pw % cc3d ( &
          bounds ( 1, 1 ) : bounds ( 2, 1 ), &
          bounds ( 1, 2 ) : bounds ( 2, 2 ), &
          bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = allocstat )
     IF ( allocstat /= 0 ) CALL stop_memory ( "pw_allocate", "cc3d", 0 )

  ELSE
     CALL stop_program ( "pw_allocate", "no possible data field" )
  END IF
  CALL timestop(0.0_dbl,handle)

END SUBROUTINE pw_allocate

!!*****
!******************************************************************************
!!****** pw_types/pw_deallocate [1.0] *
!!
!!   NAME
!!     pw_deallocate
!!
!!   SYNOPSIS
!!     Subroutine pw_deallocate(pw)
!!       Type(pw_type), Intent (INOUT):: pw
!!     End Subroutine pw_deallocate
!!
!!   FUNCTION
!!     Deallocate a plane wave type variable
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_deallocate ( pw )


    TYPE(pw_type), INTENT(INOUT)             :: pw

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_deallocate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: allocstat, handle

!------------------------------------------------------------------------------

  CALL timeset(routineN,"I","",handle)
  allocated_pw_count = allocated_pw_count - 1
  IF ( pw % in_use == REALDATA1D ) THEN
     DEALLOCATE ( pw % cr, STAT = allocstat )
     IF ( allocstat /= 0 ) CALL stop_memory ( "pw_deallocate", "cr" )
  ELSE IF ( pw % in_use == COMPLEXDATA1D ) THEN
     DEALLOCATE ( pw % cc, STAT = allocstat )
     IF ( allocstat /= 0 ) CALL stop_memory ( "pw_deallocate", "cc" )
  ELSE IF ( pw % in_use == REALDATA3D ) THEN
     IF ( ASSOCIATED(pw % cr3d) ) THEN
        DEALLOCATE ( pw % cr3d, STAT = allocstat )
        IF ( allocstat /= 0 ) CALL stop_memory ( "pw_deallocate", "cr3d" )
     END IF
  ELSE IF ( pw % in_use == COMPLEXDATA3D ) THEN
     DEALLOCATE ( pw % cc3d, STAT = allocstat )
     IF ( allocstat /= 0 ) CALL stop_memory ( "pw_deallocate", "cc3d" )
  ELSE IF ( pw % in_use /= NODATA ) THEN
     CALL stop_program ( "pw_deallocate", "no possible data field" )
  END IF
  CALL timestop(0.0_dbl,handle)

END SUBROUTINE pw_deallocate

!!*****
!******************************************************************************
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
    REAL(dbl)                                :: zr

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_zero', 'I', 'Mzero', handle )

  IF ( pw % in_use == REALDATA1D ) THEN
     ns = SIZE ( pw % cr  )
!$omp parallel do private(i)
     DO i = 1, ns
        pw % cr ( i ) = 0.0_dbl
     END DO
  ELSE IF ( pw % in_use == COMPLEXDATA1D ) THEN
     ns = 2 * SIZE ( pw % cc  )
!$omp parallel do private(i)
     DO i = 1, ns/2
        pw % cc ( i ) = CMPLX ( 0.0_dbl, 0.0_dbl, dbl )
     END DO
  ELSE IF ( pw % in_use == REALDATA3D ) THEN
     ns = SIZE ( pw % cr3d  )
     CALL DCOPY(ns,0.0_dbl,0,pw % cr3d,1)
  ELSE IF ( pw % in_use == COMPLEXDATA3D ) THEN
     ns = 2 * SIZE ( pw % cc3d  )
     CALL DCOPY(ns,0.0_dbl,0,pw % cc3d,1) 
  ELSE
     CALL stop_program ( "pw_zero", "no possible data field" )
  END IF

  zr = REAL ( ns, dbl ) * 1.e-6_dbl
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
!!     JGH (7-Mar-2001) : check for pw_grid % identifier, allow copy if 
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
    REAL(dbl)                                :: zc

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_copy', 'I', 'Mcopy', handle )

  IF ( pw1 % pw_grid % identifier /= pw2 % pw_grid % identifier ) THEN

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
                 pw2%cc(i) = CMPLX ( 0.0_dbl, 0.0_dbl, dbl )
              END DO
           END IF
           ns = 2 * MAX ( ng1, ng2 )
        ELSE
           CALL stop_program ( "pw_copy", "no suitable data field" )
        END IF

      ELSE
        WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                 " grid 1 :",pw1 % pw_grid % identifier, &
                 " sperical :",pw1 % pw_grid % spherical, &
                 " reference :",pw1 % pw_grid % reference 
        WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                 " grid 2 :",pw2 % pw_grid % identifier, &
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

        IF ( ( pw1 % pw_grid % identifier == pw2 % pw_grid % reference ) ) THEN
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
        ELSE IF ( ( pw2 % pw_grid % identifier == pw1 % pw_grid % reference ) ) THEN
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
                 " grid 1 :",pw1 % pw_grid % identifier, &
                 " sperical :",pw1 % pw_grid % spherical, &
                 " reference :",pw1 % pw_grid % reference 
          WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
                 " grid 2 :",pw2 % pw_grid % identifier, &
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
               " grid 1 :",pw1 % pw_grid % identifier, &
               " sperical :",pw1 % pw_grid % spherical, &
               " reference :",pw1 % pw_grid % reference 
      WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 2 :",pw2 % pw_grid % identifier, &
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
       ns = 2 * SIZE ( pw1 % cc )
!$omp parallel do private(i)
       DO i = 1, ns/2
          pw2 % cc(i) = pw1 % cc(i)
       END DO
    ELSE IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
       ns = SIZE ( pw1 % cr3d )
       CALL DCOPY(ns, pw1 % cr3d, 1, pw2 % cr3d, 1)
    ELSE IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
         pw2 % in_use == COMPLEXDATA3D ) THEN
       ns = 2 * SIZE ( pw1 % cc3d )
       CALL DCOPY(ns, pw1 % cc3d, 1, pw2 % cc3d, 1)
    ELSE
       CALL stop_program ( "pw_copy", "no suitable data field" )
    END IF

  END IF

  pw2 % in_space = pw1 % in_space

  zc = REAL ( ns, dbl ) * 1.e-6_dbl
  CALL timestop ( zc, handle )
  
END SUBROUTINE pw_copy

!!*****
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

    COMPLEX(dbl)                             :: im
    INTEGER                                  :: cnt, handle, i, m
    REAL(dbl)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_derive', 'I', 'Mflops', handle )

  IF ( ANY ( n < 0 ) ) CALL stop_program ( "pw_derive", "n < 0" )

  m = SUM ( n )
  im = CMPLX ( 0._dbl, 1._dbl, dbl ) ** m

  flop = 0._dbl

  IF ( pw % in_space == RECIPROCALSPACE .AND. &
       pw % in_use == COMPLEXDATA1D ) THEN

    cnt = REAL ( SIZE ( pw % cc ), dbl )

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
    IF ( ABS ( REAL ( im, dbl ) - 1._dbl ) > 1.e-10 ) THEN
!$omp parallel do private (i)
      DO i = 1, cnt
         pw % cc ( i ) = im * pw % cc ( i )
      END DO
      flop = flop + 6 * cnt
    END IF

  ELSE

    CALL stop_program ( "pw_derive", "no suitable data field" )

  END IF

  flop = flop * 1.e-6_dbl
  CALL timestop ( flop, handle )

END SUBROUTINE pw_derive

!******************************************************************************
!!****** pw_types/pw_get_left_right [1.0] *
!!
!!   NAME
!!     pw_get_left_right
!!
!!   SYNOPSIS
!!     Subroutine pw_get_left_right(pw, left, right )
!!       Implicit None
!!       Type(pw_type), Intent (IN):: pw
!!       Real(dbl), DIMENSION(:,:), POINTER :: left,right
!!     End Subroutine pw_derive
!!
!!   FUNCTION
!!     Gets the slice on the left and on the right of the pw_grid owned by the current cpu
!!     user allocates left & right 
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE pw_get_left_right ( pw, left, right )

    TYPE(pw_type), INTENT(IN)                :: pw
    REAL(KIND=dbl), DIMENSION(:, :), POINTER :: left, right

    INTEGER                                  :: left_cpu, left_x, right_cpu, &
                                                right_x
    INTEGER, DIMENSION(2, 3)                 :: bl
    REAL(KIND=dbl), DIMENSION(:, :), POINTER :: send

!!*****

  bl = pw%pw_grid%bounds_local
  IF ( pw % in_space /= REALSPACE .OR.  pw % in_use /= REALDATA3D ) THEN
     CALL stop_program("pw_get_left_right","wrong pw type")
  ENDIF
  IF ( .NOT. ASSOCIATED(left) .OR. .NOT. ASSOCIATED(right)) THEN
     CALL stop_program("pw_get_left_right","left right not associated")
  ENDIF
  IF ( .NOT. ASSOCIATED(left) .OR. .NOT. ASSOCIATED(right)) THEN
     CALL stop_program("pw_get_left_right","left right not associated")
  ENDIF

  ! if we have a part of x we try to get the neighbor (trying to avoid mistakes if CPU>npoints_x)
  ALLOCATE(  send(bl(1,2):bl(2,2) , bl(1,3):bl(2,3)) )
  left_x=bl(1,1)-1
  IF (left_x  .LT. pw%pw_grid%bounds(1,1)) left_x =pw%pw_grid%bounds(2,1)
  left_cpu=pw%pw_grid%para%pos_of_x(left_x)
  right_x=bl(2,1)+1
  IF (right_x .GT. pw%pw_grid%bounds(2,1)) right_x=pw%pw_grid%bounds(1,1)
  right_cpu=pw%pw_grid%para%pos_of_x(right_x)
  ! send the left data (receive the right)
  send = pw%cr3d(bl(1,1),:,:)
  CALL mp_sendrecv(send,left_cpu,right,right_cpu,pw%pw_grid%para%group)
  ! send the right data (receive the left)
  send = pw%cr3d(bl(1,2),:,:)
  CALL mp_sendrecv(send,right_cpu,left,left_cpu,pw%pw_grid%para%group)
  DEALLOCATE(send)

END SUBROUTINE pw_get_left_right

!******************************************************************************
!!****** pw_types/pw_derive_fd [1.0] *
!!
!!   NAME
!!     pw_derive_fd
!!
!!   SYNOPSIS
!!     Subroutine pw_derive_fd(pw, pw_der, n)
!!       Type(pw_type), Intent (IN):: pw
!!       Type(pw_type), Intent (INOUT):: pw_der
!!       Integer, Dimension(3), Intent (IN):: n
!!     End Subroutine pw_derive_fd
!!
!!   FUNCTION
!!     Calculate the derivative of a plane wave vector using finite
!!     differences. Only for cr3d and order 1 right now.
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_derive_fd ( pw, pw_der, n )

    TYPE(pw_type), INTENT(IN)                :: pw
    TYPE(pw_type), INTENT(INOUT)             :: pw_der
    INTEGER, DIMENSION(3), INTENT(IN)        :: n

    INTEGER                                  :: handle, i, j, k
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wrap
    INTEGER, DIMENSION(2, 3)                 :: bl
    REAL(dbl)                                :: itwodelta
    REAL(KIND=dbl), DIMENSION(:, :), POINTER :: left, right

!!*****

  CALL timeset ( 'pw_derive_fd', 'I', '', handle )
 
  ! just one component one for the time being
  IF (ANY(n.LT.0) .OR. SUM(n).NE.1) THEN
     CALL stop_program("pw_derive_fd","wrong derivative order") 
  ENDIF
  ! only cr3d supported
  IF ( pw % in_space /= REALSPACE .OR.  pw % in_use /= REALDATA3D ) THEN
     CALL stop_program("pw_derive_fd","wrong pw type") 
  ENDIF
  ! only distribution over x is allowed
  IF (ANY(pw%pw_grid%bounds_local(:,2:3) .NE. pw%pw_grid%bounds(:,2:3))) THEN
     CALL stop_program("pw_derive_fd","wrong pw distribution") 
  ENDIF
  ! pw and pw_der should be fully equivalent. Only one test here.
  IF ( pw_der % in_space /= REALSPACE .OR.  pw_der % in_use /= REALDATA3D ) THEN
     CALL stop_program("pw_derive_fd","wrong pw_der type") 
  ENDIF

  bl = pw%pw_grid%bounds_local
 
  ! fully local derivs, z
  IF (n(3) .EQ. 1) THEN
    itwodelta = 1.0_dbl/( 2.0_dbl * pw%pw_grid%dr(3) )

    ! precomputed wrapping for PBC
    ALLOCATE(wrap(bl(1,3)-1:bl(2,3)+1))
    DO k=bl(1,3),bl(2,3)
       wrap(k)=k
    ENDDO
    wrap(bl(1,3)-1)=bl(2,3)
    wrap(bl(2,3)+1)=bl(1,3)

    ! deriv
    DO k=bl(1,3),bl(2,3)
    DO j=bl(1,2),bl(2,2)
    DO i=bl(1,1),bl(2,1)
       pw_der % cr3d(i,j,k) = (pw % cr3d(i,j,wrap(k+1)) - pw % cr3d(i,j,wrap(k-1)))*itwodelta
    ENDDO
    ENDDO
    ENDDO

    DEALLOCATE(wrap)
  ENDIF

  ! fully local derivs, y
  IF (n(2) .EQ. 1) THEN
    itwodelta = 1.0_dbl/( 2.0_dbl * pw%pw_grid%dr(2) )

    ! precomputed wrapping for PBC
    ALLOCATE(wrap(bl(1,2)-1:bl(2,2)+1))
    DO j=bl(1,2),bl(2,2)
       wrap(j)=j
    ENDDO
    wrap(bl(1,2)-1)=bl(2,2)
    wrap(bl(2,2)+1)=bl(1,2)

    ! deriv
    DO k=bl(1,3),bl(2,3)
    DO j=bl(1,2),bl(2,2)
    DO i=bl(1,1),bl(2,1)
       pw_der % cr3d(i,j,k) = (pw % cr3d(i,wrap(j+1),k) - pw % cr3d(i,wrap(j-1),k))*itwodelta
    ENDDO
    ENDDO
    ENDDO

    DEALLOCATE(wrap)
  ENDIF

  ! possibly non-local deriv, x
  IF (n(1) .EQ. 1) THEN
     ! treat local case first
     IF ( ALL(pw%pw_grid%bounds(1:2,1) .EQ. pw%pw_grid%bounds_local(1:2,1)) ) THEN

        itwodelta = 1.0_dbl/( 2.0_dbl * pw%pw_grid%dr(1) )
        ! precomputed wrapping for PBC
        ALLOCATE(wrap(bl(1,1)-1:bl(2,1)+1))
        DO i=bl(1,1),bl(2,1)
          wrap(i)=i
        ENDDO
        wrap(bl(1,1)-1)=bl(2,1)
        wrap(bl(2,1)+1)=bl(1,1)

        ! deriv
        DO k=bl(1,3),bl(2,3)
        DO j=bl(1,2),bl(2,2)
        DO i=bl(1,1),bl(2,1)
           pw_der % cr3d(i,j,k) = (pw % cr3d(wrap(i+1),j,k) - pw % cr3d(wrap(i-1),j,k))*itwodelta
        ENDDO
        ENDDO
        ENDDO

        DEALLOCATE(wrap)
     ELSE ! the non-local case, involves exchanging boundaries with other CPUs.

        IF (bl(2,1)-bl(1,1).GT.0) THEN

           ALLOCATE(  left(bl(1,2):bl(2,2) , bl(1,3):bl(2,3)) )
           ALLOCATE( right(bl(1,2):bl(2,2) , bl(1,3):bl(2,3)) )

           CALL pw_get_left_right(pw,left,right)

           itwodelta = 1.0_dbl/( 2.0_dbl * pw%pw_grid%dr(1) )
           ! deriv the inner part
           DO k=bl(1,3),bl(2,3)
           DO j=bl(1,2),bl(2,2)
           DO i=bl(1,1)+1,bl(2,1)-1
              pw_der % cr3d(i,j,k) = (pw % cr3d(i+1,j,k) - pw % cr3d(i-1,j,k))*itwodelta
           ENDDO
           ENDDO
           ENDDO

           ! deriv the left part
           i=bl(1,1)
           DO k=bl(1,3),bl(2,3)
           DO j=bl(1,2),bl(2,2)
              pw_der % cr3d(i,j,k) = ( pw % cr3d(i+1,j,k) - left(j,k))*itwodelta
           ENDDO
           ENDDO

           ! deriv the right part
           i=bl(2,1)
           DO k=bl(1,3),bl(2,3)
           DO j=bl(1,2),bl(2,2)
              pw_der % cr3d(i,j,k) = (right(j,k) - pw % cr3d(i-1,j,k))*itwodelta
           ENDDO
           ENDDO

           DEALLOCATE(left,right)
        ENDIF
     ENDIF
  ENDIF

  CALL timestop(0.0_dbl,handle)

END SUBROUTINE pw_derive_fd

! ******************************************************************************
!          prolongates the values of a coarse grid to a fine grid
!          uses simple linear interpolation
! ******************************************************************************
SUBROUTINE pw_prolongate(pw_coarse, pw_fine )

    TYPE(pw_type), INTENT(IN)                :: pw_coarse
    TYPE(pw_type), INTENT(INOUT)             :: pw_fine

    INTEGER                                  :: handle, i, j, k
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wrapx, wrapy, wrapz
    INTEGER, DIMENSION(3)                    :: lb, n, ub

  CALL timeset ( 'pw_prolongate', 'I', '', handle )

  ! only cr3d supported
  IF ( pw_coarse % in_space /= REALSPACE .OR.  pw_coarse % in_use /= REALDATA3D ) THEN
     CALL stop_program("pw_prolongate","wrong pw type (I)")
  ENDIF
  ! only cr3d supported
  IF ( pw_fine % in_space /= REALSPACE .OR.  pw_fine % in_use /= REALDATA3D ) THEN
     CALL stop_program("pw_prolongate","wrong pw type (II)")
  ENDIF
  ! only local for the time being
  IF ( pw_coarse % pw_grid % para % mode /= PW_MODE_LOCAL ) THEN
     CALL stop_program("pw_prolongate","sorry local only")
  ENDIF
  ! only for the easiest case 
  IF ( ANY(pw_coarse % pw_grid % npts * 2 .NE. pw_fine % pw_grid % npts )) THEN
     CALL stop_program("pw_prolongate","sorry too difficult ")
  ENDIF
  ! we use the bounds of the smaller grid
  n   = pw_coarse % pw_grid % npts
  lb  = pw_coarse % pw_grid % bounds_local(1,:)
  ub  = pw_coarse % pw_grid % bounds_local(2,:)

  ! we'll need a wrap for PBC
  ALLOCATE(wrapx(lb(1)*2-1:ub(1)*2+2))
  DO i=lb(1)*2,ub(1)*2+1
     wrapx(i)=i
  ENDDO
  wrapx(lb(1)*2-1)=ub(1)*2+1
  wrapx(ub(1)*2+2)=lb(1)*2
  ALLOCATE(wrapy(lb(2)*2-1:ub(2)*2+2))
  DO i=lb(2)*2,ub(2)*2+1
     wrapy(i)=i
  ENDDO
  wrapy(lb(2)*2-1)=ub(2)*2+1
  wrapy(ub(2)*2+2)=lb(2)*2
  ALLOCATE(wrapz(lb(3)*2-1:ub(3)*2+2))
  DO i=lb(3)*2,ub(3)*2+1
     wrapz(i)=i
  ENDDO
  wrapz(lb(3)*2-1)=ub(3)*2+1
  wrapz(ub(3)*2+2)=lb(3)*2

  ! copy the indentical points
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_fine % cr3d(2*i,2*j,2*k ) = pw_coarse % cr3d(i,j,k)
    ENDDO
   ENDDO
  ENDDO
  ! other points are averages of the points we have already. The points we have are the corners
  ! of a cube, we fill in edges, faces and center, in this order
  ! put values on the x edges
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_fine % cr3d(2*i+1,2*j,2*k ) = 0.5_dbl* (pw_fine % cr3d(2*i,2*j,2*k ) + &
                                                  pw_fine % cr3d(wrapx(2*i+2),2*j,2*k ))
    ENDDO
   ENDDO
  ENDDO
  ! put values on the y edges
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_fine % cr3d(2*i,2*j+1,2*k ) = 0.5_dbl* (pw_fine % cr3d(2*i,2*j,2*k ) +  &
                                                  pw_fine % cr3d(2*i,wrapy(2*j+2),2*k ))
    ENDDO
   ENDDO
  ENDDO
  ! put values on the z edges
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_fine % cr3d(2*i,2*j,2*k+1 ) = 0.5_dbl* (pw_fine % cr3d(2*i,2*j,2*k ) +  &
                                                  pw_fine % cr3d(2*i,2*j,wrapz(2*k+2)))
    ENDDO
   ENDDO
  ENDDO
  ! put values on the xy faces
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_fine % cr3d(2*i+1,2*j+1,2*k ) = 0.5_dbl* (pw_fine % cr3d(2*i,2*j+1,2*k ) + &  
                                                    pw_fine % cr3d(wrapx(2*i+2),2*j+1,2*k ))
    ENDDO
   ENDDO
  ENDDO
  ! put values on the xz faces
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_fine % cr3d(2*i+1,2*j,2*k+1 ) = 0.5_dbl* (pw_fine % cr3d(2*i,2*j,2*k+1 ) + &
                                                    pw_fine % cr3d(wrapx(2*i+2),2*j,2*k+1 ))
    ENDDO
   ENDDO
  ENDDO
  ! put values on the yz faces
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_fine % cr3d(2*i,2*j+1,2*k+1 ) = 0.5_dbl* (pw_fine % cr3d(2*i,2*j,2*k+1 ) + &
                                                    pw_fine % cr3d(2*i,wrapy(2*j+2),2*k+1 ))
    ENDDO
   ENDDO
  ENDDO
  ! put values on the center
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_fine % cr3d(2*i+1,2*j+1,2*k+1 ) = 0.5_dbl* (pw_fine % cr3d(2*i ,2*j+1,2*k+1 ) + &
                                                      pw_fine % cr3d(wrapx(2*i+2),2*j+1,2*k+1 ))
    ENDDO
   ENDDO
  ENDDO
  DEALLOCATE(wrapx,wrapy,wrapz)
  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE pw_prolongate


! ******************************************************************************
!          injectes the values of a fine grid to a coarse grid
!          uses simple linear interpolation, and is the 'inverse' of the pw_prolongate
!          i.e. read prolongate backwards
!          xx WARNING xx for efficiency reasons, pw_fine gets destroyed ! sorry...
! ******************************************************************************
SUBROUTINE pw_inject(pw_fine, pw_coarse )

    TYPE(pw_type), INTENT(INOUT)             :: pw_fine, pw_coarse

    INTEGER                                  :: handle, i, j, k
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: wrapx, wrapy, wrapz
    INTEGER, DIMENSION(3)                    :: lb, n, ub
    REAL(KIND=dbl)                           :: dum

  CALL timeset ( 'pw_inject', 'I', '', handle )

  ! only cr3d supported
  IF ( pw_coarse % in_space /= REALSPACE .OR.  pw_coarse % in_use /= REALDATA3D ) THEN
     CALL stop_program("pw_inject","wrong pw type (I)")
  ENDIF
  ! only cr3d supported
  IF ( pw_fine % in_space /= REALSPACE .OR.  pw_fine % in_use /= REALDATA3D ) THEN
     CALL stop_program("pw_inject","wrong pw type (II)")
  ENDIF
  ! only local for the time being
  IF ( pw_coarse % pw_grid % para % mode /= PW_MODE_LOCAL ) THEN
     CALL stop_program("pw_inject","sorry local only")
  ENDIF
  ! only for the easiest case
  IF ( ANY(pw_coarse % pw_grid % npts * 2 .NE. pw_fine % pw_grid % npts )) THEN
     CALL stop_program("pw_inject","sorry too difficult ")
  ENDIF
  n   = pw_coarse % pw_grid % npts
  lb  = pw_coarse % pw_grid % bounds_local(1,:)
  ub  = pw_coarse % pw_grid % bounds_local(2,:)

  ! we'll need a wrap for PBC
  ALLOCATE(wrapx(lb(1)*2-1:ub(1)*2+2))
  DO i=lb(1)*2,ub(1)*2+1
     wrapx(i)=i
  ENDDO
  wrapx(lb(1)*2-1)=ub(1)*2+1
  wrapx(ub(1)*2+2)=lb(1)*2
  ALLOCATE(wrapy(lb(2)*2-1:ub(2)*2+2))
  DO i=lb(2)*2,ub(2)*2+1
     wrapy(i)=i
  ENDDO
  wrapy(lb(2)*2-1)=ub(2)*2+1
  wrapy(ub(2)*2+2)=lb(2)*2
  ALLOCATE(wrapz(lb(3)*2-1:ub(3)*2+2))
  DO i=lb(3)*2,ub(3)*2+1
     wrapz(i)=i
  ENDDO
  wrapz(lb(3)*2-1)=ub(3)*2+1
  wrapz(ub(3)*2+2)=lb(3)*2

  ! values from the center
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       dum = 0.5_dbl * pw_fine % cr3d(2*i+1,2*j+1,2*k+1 )
       pw_fine % cr3d(2*i ,2*j+1,2*k+1 )         = pw_fine % cr3d(2*i ,2*j+1,2*k+1 )+ dum
       pw_fine % cr3d(wrapx(2*i+2),2*j+1,2*k+1 ) = pw_fine % cr3d(wrapx(2*i+2),2*j+1,2*k+1 )+ dum
    ENDDO
   ENDDO
  ENDDO
  ! values from the yz faces
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       dum = 0.5_dbl * pw_fine % cr3d(2*i,2*j+1,2*k+1 )
       pw_fine % cr3d(2*i,2*j,2*k+1 )            = pw_fine % cr3d(2*i,2*j,2*k+1 ) + dum
       pw_fine % cr3d(2*i,wrapy(2*j+2),2*k+1 )   = pw_fine % cr3d(2*i,wrapy(2*j+2),2*k+1 ) + dum
    ENDDO
   ENDDO
  ENDDO
  ! values from the xz faces
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       dum = 0.5_dbl * pw_fine % cr3d(2*i+1,2*j,2*k+1 )
       pw_fine % cr3d(2*i,2*j,2*k+1 )          = pw_fine % cr3d(2*i,2*j,2*k+1 ) + dum
       pw_fine % cr3d(wrapx(2*i+2),2*j,2*k+1 ) = pw_fine % cr3d(wrapx(2*i+2),2*j,2*k+1 ) + dum
    ENDDO
   ENDDO
  ENDDO
  ! values from the xy faces
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       dum = 0.5_dbl * pw_fine % cr3d(2*i+1,2*j+1,2*k )
       pw_fine % cr3d(2*i,2*j+1,2*k )          = pw_fine % cr3d(2*i,2*j+1,2*k ) + dum
       pw_fine % cr3d(wrapx(2*i+2),2*j+1,2*k ) = pw_fine % cr3d(wrapx(2*i+2),2*j+1,2*k ) + dum
    ENDDO
   ENDDO
  ENDDO
  ! values from the z edges
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       dum = 0.5_dbl * pw_fine % cr3d(2*i,2*j,2*k+1 )
       pw_fine % cr3d(2*i,2*j,2*k )         = pw_fine % cr3d(2*i,2*j,2*k ) + dum
       pw_fine % cr3d(2*i,2*j,wrapz(2*k+2)) = pw_fine % cr3d(2*i,2*j,wrapz(2*k+2)) + dum
    ENDDO
   ENDDO
  ENDDO
  ! values from the y edges
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       dum = 0.5_dbl * pw_fine % cr3d(2*i,2*j+1,2*k )
       pw_fine % cr3d(2*i,2*j,2*k )          = pw_fine % cr3d(2*i,2*j,2*k ) + dum
       pw_fine % cr3d(2*i,wrapy(2*j+2),2*k ) = pw_fine % cr3d(2*i,wrapy(2*j+2),2*k ) + dum
    ENDDO
   ENDDO
  ENDDO
  ! values from the x edges
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       dum = 0.5_dbl * pw_fine % cr3d(2*i+1,2*j,2*k )
       pw_fine % cr3d(2*i,2*j,2*k )          = pw_fine % cr3d(2*i,2*j,2*k ) + dum
       pw_fine % cr3d(wrapx(2*i+2),2*j,2*k ) = pw_fine % cr3d(wrapx(2*i+2),2*j,2*k ) + dum
    ENDDO
   ENDDO
  ENDDO
  ! copy the indentical points
  DO k=lb(3),ub(3)
   DO j=lb(2),ub(2)
    DO i=lb(1),ub(1)
       pw_coarse % cr3d(i,j,k) = pw_fine % cr3d(2*i,2*j,2*k )
    ENDDO
   ENDDO
  ENDDO

  DEALLOCATE(wrapx,wrapy,wrapz)
  CALL timestop ( 0.0_dbl, handle )

END SUBROUTINE pw_inject

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

SUBROUTINE pw_transfer ( pw1, pw2 )


    TYPE(pw_type), INTENT(IN)                :: pw1
    TYPE(pw_type), INTENT(INOUT)             :: pw2

!------------------------------------------------------------------------------

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
     CALL pw_fft_wrap ( pw1, pw2 )

  END IF

END SUBROUTINE pw_transfer

!!*****
!******************************************************************************
!!****** pw_types/pw_add [1.0] *
!!
!!   NAME
!!     pw_add
!!
!!   SYNOPSIS
!!     Subroutine pw_add(pw1, pw2, pw3)
!!       Type(pw_type), Intent (IN):: pw1, pw2
!!       Type(pw_type), Intent (INOUT):: pw3
!!     End Subroutine pw_add
!!
!!   FUNCTION
!!     pw3 = pw1 + pw2
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     CJM (23-May-01)  : check for pw_grid % identifier, allow copy if 
!!       in_use == COMPLEXDATA1D and in_space == RECIPROCALSPACE
!!     JGH (21-Feb-03)  : generalization to all reference grid types
!!
!!   NOTES
!!     Currently only adding of respective types allowed,
!!     in order to avoid errors
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_add ( pw1, pw2, pw3 )
  
  
    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    TYPE(pw_type), INTENT(INOUT)             :: pw3

    INTEGER                                  :: handle, i, j, ng, ng1, ng2, &
                                                ng3, ngm
    REAL(dbl)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_add', 'I', 'Mflops', handle )

  IF ( pw3 % pw_grid % identifier /= pw2 % pw_grid % identifier .OR. &
       pw3 % pw_grid % identifier /= pw1 % pw_grid % identifier  ) THEN

    ! different grids: this works only with reference grids

    IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
         pw2 % in_use == COMPLEXDATA1D .AND. &
         pw3 % in_use == COMPLEXDATA1D .AND. &
         pw1 % in_space == RECIPROCALSPACE .AND. &
         pw2 % in_space == RECIPROCALSPACE ) THEN

       ng1 = SIZE ( pw1 % cc )
       ng2 = SIZE ( pw2 % cc )
       ng3 = SIZE ( pw3 % cc )
       
       IF ( pw1 % pw_grid % spherical .AND. pw2 % pw_grid % spherical &
                            .AND. pw3 % pw_grid % spherical ) THEN

         IF ( pw_compatible ( pw1 % pw_grid, pw3 % pw_grid ) .AND. &
              pw_compatible ( pw2 % pw_grid, pw3 % pw_grid ) ) THEN

            ng = MIN ( ng1, ng2, ng3 )
!$omp parallel do private(i)
            DO i = 1,ng
               pw3 % cc ( i ) = pw2 % cc ( i ) + pw1 % cc ( i )
            END DO
            IF ( ng3 > ng ) THEN
!$omp parallel do private(i)
               DO i = ng+1, ng3
                  pw3 % cc ( i ) = CMPLX ( 0._dbl, 0._dbl, dbl )
               END DO
               IF ( ng1 > ng ) THEN
                  ngm = MIN ( ng3, ng1 )
!$omp parallel do private(i)
                  DO i = ng+1, ngm
                     pw3 % cc ( i ) = pw1 % cc ( i )
                  END DO
               END IF
               IF ( ng2 > ng ) THEN
                  ngm = MIN ( ng3, ng2 )
!$omp parallel do private(i)
                  DO i = ng+1, ngm
                     pw3 % cc ( i ) = pw2 % cc ( i )
                  END DO
               END IF
            END IF

         ELSE
            WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 1 :",pw1 % pw_grid % identifier, &
               " sperical :",pw1 % pw_grid % spherical, &
               " reference :",pw1 % pw_grid % reference 
            WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 2 :",pw2 % pw_grid % identifier, &
               " sperical :",pw2 % pw_grid % spherical, &
               " reference :",pw2 % pw_grid % reference 
            WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 3 :",pw3 % pw_grid % identifier, &
               " sperical :",pw3 % pw_grid % spherical, &
               " reference :",pw3 % pw_grid % reference 
            CALL stop_program ( "pw_add", "incompatible grids" )
         END IF

       ELSE IF ( .NOT. (pw1 % pw_grid % spherical .OR. &
                        pw2 % pw_grid % spherical .OR. &
                        pw3 % pw_grid % spherical) ) THEN
         
         CALL pw_copy ( pw1, pw3 )

         IF ( pw_compatible ( pw2 % pw_grid, pw3 % pw_grid ) ) THEN

           IF ( pw2 % pw_grid % identifier == pw3 % pw_grid % identifier ) THEN

!$omp parallel do private(i)
             DO i = 1, ng3
                pw3 % cc ( i ) = pw3 % cc ( i ) + pw2 % cc ( i )
             END DO

           ELSEIF ( ( pw2 % pw_grid % identifier == pw3 % pw_grid % reference ) ) THEN

             IF( ng2 >= ng3 ) THEN
!$omp parallel do private(i,j)
               DO i = 1, ng3
                 j = pw3 % pw_grid % gidx ( i )
                 pw3 % cc ( i ) = pw3 % cc ( i ) + pw2 % cc ( j )
               END DO
             ELSE
!$omp parallel do private(i,j)
               DO i = 1, ng2
                 j = pw3 % pw_grid % gidx ( i )
                 pw3 % cc ( j ) = pw3 % cc ( j ) + pw2 % cc ( i )
               END DO
             END IF

           ELSEIF ( ( pw3 % pw_grid % identifier == pw2 % pw_grid % reference ) ) THEN

             IF( ng2 >= ng3 ) THEN
!$omp parallel do private(i,j)
               DO i = 1, ng3
                 j = pw2 % pw_grid % gidx ( i )
                 pw3 % cc ( i ) = pw3 % cc ( i ) + pw2 % cc ( j )
               END DO
             ELSE
!$omp parallel do private(i,j)
               DO i = 1, ng2
                 j = pw2 % pw_grid % gidx ( i )
                 pw3 % cc ( j ) = pw3 % cc ( j ) + pw2 % cc ( i )
               END DO
             END IF

           END IF

         ELSE
            WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 1 :",pw1 % pw_grid % identifier, &
               " sperical :",pw1 % pw_grid % spherical, &
               " reference :",pw1 % pw_grid % reference 
            WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 2 :",pw2 % pw_grid % identifier, &
               " sperical :",pw2 % pw_grid % spherical, &
               " reference :",pw2 % pw_grid % reference 
            WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
               " grid 3 :",pw3 % pw_grid % identifier, &
               " sperical :",pw3 % pw_grid % spherical, &
               " reference :",pw3 % pw_grid % reference 
            CALL stop_program ( "pw_add", "Incompatible grids" )
         END IF

       ELSE

          WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
             " grid 1 :",pw1 % pw_grid % identifier, &
             " sperical :",pw1 % pw_grid % spherical, &
             " reference :",pw1 % pw_grid % reference 
          WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
             " grid 2 :",pw2 % pw_grid % identifier, &
             " sperical :",pw2 % pw_grid % spherical, &
             " reference :",pw2 % pw_grid % reference 
          WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
             " grid 3 :",pw3 % pw_grid % identifier, &
             " sperical :",pw3 % pw_grid % spherical, &
             " reference :",pw3 % pw_grid % reference 
          CALL stop_program ( "pw_add", "incompatible grids" )

       END IF

       flop = 2 * ng3
       pw3 % in_space = RECIPROCALSPACE

    ELSE
       CALL stop_program ( "pw_add", "no suitable data field" )
    END IF

  ELSE 

    IF ( pw1 % in_use == REALDATA1D .AND. pw2 % in_use == REALDATA1D .AND. &
         pw3 % in_use == REALDATA1D ) THEN
!$omp parallel do private(i)
       DO i = 1, SIZE ( pw3 % cr )
          pw3 % cr ( i ) = pw1 % cr ( i ) + pw2 % cr ( i )
       END DO
       flop = REAL ( SIZE ( pw2 % cr ), dbl )
    ELSE IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
         pw2 % in_use == COMPLEXDATA1D .AND. pw3 % in_use == COMPLEXDATA1D ) THEN
!$omp parallel do private(i)
       DO i = 1, SIZE ( pw3 % cc )
          pw3 % cc ( i ) = pw1 % cc ( i ) + pw2 % cc ( i )
       END DO
       flop = REAL ( 2 * SIZE ( pw2 % cc ), dbl )
    ELSE IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D .AND. &
         pw3 % in_use == REALDATA3D ) THEN
       pw3 % cr3d = pw1 % cr3d + pw2 % cr3d
       flop = REAL ( SIZE ( pw2 % cr3d ), dbl )
    ELSE IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
         pw2 % in_use == COMPLEXDATA3D .AND. pw3 % in_use == COMPLEXDATA3D ) THEN
       pw3 % cc3d = pw1 % cc3d + pw2 % cc3d
       flop = REAL ( 2 * SIZE ( pw2 % cc3d ), dbl )
    ELSE
       CALL stop_program ( "pw_add", "no suitable data field" )
    END IF

  END IF

  flop = flop * 1.e-6_dbl
  CALL timestop ( flop, handle )
  
END SUBROUTINE pw_add

!!*****
!******************************************************************************
!!****** pw_types/pw_sumup_square [1.0] *
!!
!!   NAME
!!     pw_sumup_square
!!
!!   SYNOPSIS
!!     Subroutine pw_sumup_square(pw1, pw2)
!!       Type(pw_type), Intent (INOUT):: pw1
!!       Type(pw_type), Intent (IN):: pw2
!!     End Subroutine pw_sumup_square
!!
!!   FUNCTION
!!     pw1 = pw1 + pw2*pw2
!!
!!   AUTHOR
!!     JGH ( 25-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_sumup_square ( pw1, pw2 )
  
  
    TYPE(pw_type), INTENT(INOUT)             :: pw1
    TYPE(pw_type), INTENT(IN)                :: pw2

    INTEGER                                  :: handle, i
    REAL(dbl)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_sumup_square', 'I', 'Mflops', handle )
  
  IF ( pw1 % pw_grid % identifier /= pw2 % pw_grid % identifier ) &
     CALL stop_program ( "pw_sumup_square", "grids not compatible" )

  IF ( pw1 % in_use == REALDATA1D ) THEN
     IF ( pw2 % in_use == REALDATA1D ) THEN
!$omp parallel do private(i)
        DO i = 1, SIZE ( pw2 % cr )
           pw1 % cr(i) = pw1 % cr(i) + pw2 % cr(i) * pw2 % cr(i)
        END DO
        flop = REAL ( 2 * SIZE ( pw2 % cr ), dbl )
     ELSE IF ( pw2 % in_use == COMPLEXDATA1D ) THEN
!$omp parallel do private(i)
        DO i = 1, SIZE ( pw2 % cc )
           pw1 % cr(i) = pw1 % cr(i) + REAL ( pw2 % cc(i) * CONJG ( pw2 % cc(i) ), dbl )
        END DO
        flop = REAL ( 7 * SIZE ( pw2 % cc ), dbl )
     ELSE
        CALL stop_program ( "pw_sumup_square", "no suitable data fields" )
     END IF
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D ) THEN
     IF ( pw2 % in_use == COMPLEXDATA1D ) THEN
!$omp parallel do private(i)
        DO i = 1, SIZE ( pw2 % cc )
           pw1 % cc(i) = pw1 % cc(i) + pw2 % cc(i) * CONJG ( pw2 % cc(i) )
        END DO
        flop = REAL ( 8 * SIZE ( pw2 % cc ), dbl )
     ELSE
        CALL stop_program ( "pw_sumup_square", "no suitable data fields" )
     END IF
  ELSE IF ( pw1 % in_use == REALDATA3D ) THEN
     IF ( pw2 % in_use == REALDATA3D ) THEN
        pw1 % cr3d = pw1 % cr3d + pw2 % cr3d * pw2 % cr3d
        flop = REAL ( 2 * SIZE ( pw2 % cr3d ), dbl )
     ELSE IF ( pw2 % in_use == COMPLEXDATA3D ) THEN
        pw1 % cr3d = pw1 % cr3d + REAL ( pw2 % cc3d * CONJG ( pw2 % cc3d ), dbl )
        flop = REAL ( 7 * SIZE ( pw2 % cr3d ), dbl )
     ELSE
        CALL stop_program ( "pw_sumup_square", "no suitable data fields" )
     END IF
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D ) THEN
     IF ( pw2 % in_use == COMPLEXDATA3D ) THEN
        pw1 % cc3d = pw1 % cc3d + pw2 % cc3d * CONJG ( pw2 % cc3d )
        flop = REAL ( 8 * SIZE ( pw2 % cr3d ), dbl )
     ELSE
        CALL stop_program ( "pw_sumup_square", "no suitable data fields" )
     END IF
  END IF

  flop = flop * 1.e-6_dbl
  CALL timestop ( flop, handle )
  
END SUBROUTINE pw_sumup_square

!!*****
!******************************************************************************
!!****** pw_types/pw_sumup [1.0] *
!!
!!   NAME
!!     pw_sumup
!!
!!   SYNOPSIS
!!     Subroutine pw_sumup(pw1, pw2)
!!       Type(pw_type), Intent (IN):: pw1
!!       Type(pw_type), Intent (INOUT):: pw2
!!     End Subroutine pw_sumup
!!
!!   FUNCTION
!!     pw2 = pw1 + pw2
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     JGH (21-Feb-2003) : added reference grid functionality
!!
!!   NOTES
!!     Currently only summing up of respective types allowed,
!!     in order to avoid errors
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_sumup ( pw1, pw2 )
 
 
    TYPE(pw_type), INTENT(IN)                :: pw1
    TYPE(pw_type), INTENT(INOUT)             :: pw2

    INTEGER                                  :: handle, i, j, ng, ng1, ng2
    REAL(dbl)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_sumup', 'I', 'Mflops', handle )

  IF ( pw1 % pw_grid % identifier == pw2 % pw_grid % identifier ) THEN

     IF ( pw1 % in_use == REALDATA1D .AND. pw2 % in_use == REALDATA1D ) THEN
!$omp parallel do private(i)
        DO i = 1, SIZE ( pw2 % cr )
           pw2 % cr(i) = pw2 % cr(i) + pw1 % cr(i)
        END DO
        flop = REAL ( SIZE ( pw2 % cr ), dbl )
     ELSE IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
          pw2 % in_use == COMPLEXDATA1D ) THEN
!$omp parallel do private(i)
        DO i = 1, SIZE ( pw2 % cc )
           pw2 % cc(i) = pw2 % cc(i) + pw1 % cc(i)
        END DO
        flop = REAL ( 2 * SIZE ( pw2 % cc ), dbl )
     ELSE IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
        pw2 % cr3d = pw2 % cr3d + pw1 % cr3d
        flop = REAL ( SIZE ( pw2 % cr3d ), dbl )
     ELSE IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
          pw2 % in_use == COMPLEXDATA3D ) THEN
        pw2 % cc3d = pw2 % cc3d + pw1 % cc3d
        flop = REAL ( 2 * SIZE ( pw2 % cc3d ), dbl )
     ELSE
        CALL stop_program ( "pw_sumup", "no suitable data field" )
     END IF

  ELSE IF ( pw_compatible ( pw1 % pw_grid, pw2 % pw_grid ) ) THEN

     IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
          pw2 % in_use == COMPLEXDATA1D .AND. &
          pw1 % in_space == RECIPROCALSPACE .AND. &
          pw1 % in_space == RECIPROCALSPACE ) THEN

        ng1 = SIZE ( pw1 % cc )
        ng2 = SIZE ( pw2 % cc )
        ng  = MIN ( ng1, ng2 )
        flop = REAL ( 2 * ng, dbl )

        IF ( pw1 % pw_grid % spherical ) THEN

!$omp parallel do private(i)
          DO i = 1, ng
             pw2 % cc (i) = pw2 % cc (i) + pw1 % cc (i)
          END DO

        ELSEIF ( ( pw1 % pw_grid % identifier == pw1 % pw_grid % reference ) ) THEN

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

       ELSEIF ( ( pw2 % pw_grid % identifier == pw1 % pw_grid % reference ) ) THEN

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

         WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
            " grid 1 :",pw1 % pw_grid % identifier, &
            " sperical :",pw1 % pw_grid % spherical, &
            " reference :",pw1 % pw_grid % reference
         WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
            " grid 2 :",pw2 % pw_grid % identifier, &
            " sperical :",pw2 % pw_grid % spherical, &
            " reference :",pw2 % pw_grid % reference
         CALL stop_program ( "pw_sumup_square", "grids not compatible" )

       END IF
    
     ELSE
       CALL stop_program ( "pw_sumup", "no suitable data field" )
     END IF

  ELSE

     WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
        " grid 1 :",pw1 % pw_grid % identifier, &
        " sperical :",pw1 % pw_grid % spherical, &
        " reference :",pw1 % pw_grid % reference
     WRITE ( *, "(A,I5,T30,A,L1,T60,A,I5)" ) &
        " grid 2 :",pw2 % pw_grid % identifier, &
        " sperical :",pw2 % pw_grid % spherical, &
        " reference :",pw2 % pw_grid % reference
     CALL stop_program ( "pw_sumup_square", "grids not compatible" )

  END IF

  flop = flop * 1.e-6_dbl
  CALL timestop ( flop, handle )
 
END SUBROUTINE pw_sumup

!!*****
!******************************************************************************
!!****** pw_types/pw_subtract [1.0] *
!!
!!   NAME
!!     pw_subtract
!!
!!   SYNOPSIS
!!     Subroutine pw_subtract(pw1, pw2, pw3)
!!       Type(pw_type), Intent (IN):: pw1, pw2
!!       Type(pw_type), Intent (INOUT):: pw3
!!     End Subroutine pw_subtract
!!
!!   FUNCTION
!!     pw3 = pw1 - pw2
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     Currently only subtracting of respective types allowed,
!!     in order to avoid errors
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_subtract ( pw1, pw2, pw3 )
  
  
    TYPE(pw_type), INTENT(IN)                :: pw1, pw2
    TYPE(pw_type), INTENT(INOUT)             :: pw3

    INTEGER                                  :: handle
    REAL(dbl)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_subtract', 'I', 'Mflops', handle )

  IF ( pw1 % pw_grid % identifier /= pw2 % pw_grid % identifier ) &
     CALL stop_program ( "pw_sumup_square", "grids not compatible" )

  IF ( pw1 % in_use == REALDATA1D .AND. pw2 % in_use == REALDATA1D .AND. &
       pw3 % in_use == REALDATA1D ) THEN
     pw3 % cr = pw1 % cr - pw2 % cr
     flop = REAL ( SIZE ( pw2 % cr ), dbl )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
       pw2 % in_use == COMPLEXDATA1D .AND. pw3 % in_use == COMPLEXDATA1D ) THEN
     pw3 % cc = pw1 % cc - pw2 % cc
     flop = REAL ( 2 * SIZE ( pw2 % cc ), dbl )
  ELSE IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D .AND. &
       pw3 % in_use == REALDATA3D ) THEN
     pw3 % cr3d = pw1 % cr3d - pw2 % cr3d
     flop = REAL ( SIZE ( pw2 % cr3d ), dbl )
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
       pw2 % in_use == COMPLEXDATA3D .AND. pw3 % in_use == COMPLEXDATA3D ) THEN
     pw3 % cc3d = pw1 % cc3d - pw2 % cc3d
     flop = REAL ( 2 * SIZE ( pw2 % cc3d ), dbl )
  ELSE
     CALL stop_program ( "pw_subtract", "no suitable data field" )
  END IF

  flop = flop * 1.e-6_dbl
  CALL timestop ( flop, handle )
  
END SUBROUTINE pw_subtract

!!*****
!******************************************************************************
!!****** pw_types/pw_reduce [1.0] *
!!
!!   NAME
!!     pw_reduce
!!
!!   SYNOPSIS
!!     Subroutine pw_reduce(pw1, pw2)
!!       Type(pw_type), Intent (IN):: pw1
!!       Type(pw_type), Intent (INOUT):: pw2
!!     End Subroutine pw_reduce
!!
!!   FUNCTION
!!     pw2 = pw2 - pw1
!!
!!   AUTHOR
!!     apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!     Currently only reduction of respective types allowed,
!!     in order to avoid errors
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pw_reduce ( pw1, pw2 )
  
  
    TYPE(pw_type), INTENT(IN)                :: pw1
    TYPE(pw_type), INTENT(INOUT)             :: pw2

    INTEGER                                  :: handle
    REAL(dbl)                                :: flop

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_reduce', 'I', 'Mflops', handle )

  IF ( pw1 % pw_grid % identifier /= pw2 % pw_grid % identifier ) &
     CALL stop_program ( "pw_sumup_square", "grids not compatible" )

  IF ( pw1 % in_use == REALDATA1D .AND. pw2 % in_use == REALDATA1D ) THEN
     pw2 % cr = pw2 % cr - pw1 % cr
     flop = REAL ( SIZE ( pw2 % cr ), dbl )
  ELSE IF ( pw1 % in_use == COMPLEXDATA1D .AND. &
       pw2 % in_use == COMPLEXDATA1D ) THEN
     pw2 % cc = pw2 % cc - pw1 % cc
     flop = REAL ( 2 * SIZE ( pw2 % cc ), dbl )
  ELSE IF ( pw1 % in_use == REALDATA3D .AND. pw2 % in_use == REALDATA3D ) THEN
     pw2 % cr3d = pw2 % cr3d - pw1 % cr3d
     flop = REAL ( SIZE ( pw2 % cr3d ), dbl )
  ELSE IF ( pw1 % in_use == COMPLEXDATA3D .AND. &
       pw2 % in_use == COMPLEXDATA3D ) THEN
     pw2 % cc3d = pw2 % cc3d - pw1 % cc3d
     flop = REAL ( 2 * SIZE ( pw2 % cc3d ), dbl )
  ELSE
     CALL stop_program ( "pw_reduce", "no suitable data field" )
  END IF

  flop = flop * 1.e-6_dbl
  CALL timestop ( flop, handle )
  
END SUBROUTINE pw_reduce

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
    COMPLEX(dbl), DIMENSION(:, :, :), &
      INTENT(IN)                             :: c
    COMPLEX(dbl), INTENT(IN), OPTIONAL       :: scale

    INTEGER                                  :: gpt, handle, l, m, n, ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat
    REAL(dbl)                                :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_gather', 'I', 'Mcopy', handle )

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
      pw % cc ( gpt ) = pw % cc ( gpt ) + scale * c ( l, m, n )

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

  cpy = REAL ( ngpts, dbl ) * 1.e-6_dbl
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_gather_s

!******************************************************************************

SUBROUTINE pw_gather_p ( pw, c, scale )


    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw
    COMPLEX(dbl), DIMENSION(:, :), &
      INTENT(IN)                             :: c
    COMPLEX(dbl), INTENT(IN), OPTIONAL       :: scale

    INTEGER                                  :: gpt, handle, l, m, mn, n, &
                                                ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat, yzq
    REAL(dbl)                                :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_gather', 'I', 'Mcopy', handle )

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
      pw % cc ( gpt ) = pw % cc ( gpt ) + scale * c ( l, mn )

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

  cpy = REAL ( ngpts, dbl ) * 1.e-6_dbl
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
    COMPLEX(dbl), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: c
    COMPLEX(dbl), INTENT(IN), OPTIONAL       :: scale

    INTEGER                                  :: gpt, handle, l, m, n, ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat
    REAL(dbl)                                :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_scatter', 'I', 'Mcopy', handle )

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

  IF ( .NOT. PRESENT ( scale ) ) c = 0._dbl

  IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      c ( l, m, n ) = c ( l, m, n ) + scale * pw % cc ( gpt )

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
        c ( l, m, n ) = c ( l, m, n ) + scale * CONJG ( pw % cc ( gpt ) )

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

  cpy = REAL ( ngpts, dbl ) * 1.e-6_dbl
  CALL timestop ( cpy, handle )

END SUBROUTINE pw_scatter_s

!******************************************************************************

SUBROUTINE pw_scatter_p ( pw, c, scale )


    TYPE(pw_type), INTENT(IN), TARGET        :: pw
    COMPLEX(dbl), DIMENSION(:, :), &
      INTENT(INOUT)                          :: c
    COMPLEX(dbl), INTENT(IN), OPTIONAL       :: scale

    INTEGER                                  :: gpt, handle, l, m, mn, n, &
                                                ngpts
    INTEGER, DIMENSION(:), POINTER           :: mapl, mapm, mapn
    INTEGER, DIMENSION(:, :), POINTER        :: ghat, yzq
    REAL(dbl)                                :: cpy

!------------------------------------------------------------------------------

  CALL timeset ( 'pw_scatter', 'I', 'Mcopy', handle )

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

  IF ( .NOT. PRESENT ( scale ) ) c = 0._dbl

  IF ( PRESENT ( scale ) ) THEN

!$omp parallel do private(gpt,l,m,n,mn)
    DO gpt = 1, ngpts

      l = mapl ( ghat ( 1, gpt ) ) + 1
      m = mapm ( ghat ( 2, gpt ) ) + 1
      n = mapn ( ghat ( 3, gpt ) ) + 1
      mn = yzq ( m, n )
      c ( l, mn ) = c ( l, mn ) + scale * pw % cc ( gpt )

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
        c ( l, mn ) = c ( l, mn ) + scale * CONJG ( pw % cc ( gpt ) )

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

  cpy = REAL ( ngpts, dbl ) * 1.e-6_dbl
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
    REAL(dbl)                                :: norm

!------------------------------------------------------------------------------

  CALL timeset("fft_wrap_pw1","I","",handle)

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
       norm = 1._dbl
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
  CALL timestop(0.0_dbl,handle) 

END SUBROUTINE fft_wrap_pw1


!******************************************************************************

SUBROUTINE fft_wrap_pw1pw2 ( pw1, pw2, debug )
  

    TYPE(pw_type), INTENT(IN), TARGET        :: pw1
    TYPE(pw_type), INTENT(INOUT), TARGET     :: pw2
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    CHARACTER(LEN=9)                         :: mode
    COMPLEX(dbl), ALLOCATABLE, &
      DIMENSION(:, :)                        :: grays
    COMPLEX(dbl), DIMENSION(:, :, :), &
      POINTER                                :: c_in, c_out
    INTEGER                                  :: dir, handle, handle2, ierr, &
                                                l1, l2, l3, my_pos, &
                                                nloc( 3 ), nrays, nsize, &
                                                out_space
    INTEGER, DIMENSION(:), POINTER           :: n
    LOGICAL                                  :: test
    REAL(dbl)                                :: norm

!------------------------------------------------------------------------------

  CALL timeset("fft_wrap_pw1pw2_all","I","",handle2)
  CALL timeset("fft_wrap_pw1pw2_"//TRIM(ADJUSTL(cp_to_string(pw1%pw_grid%identifier))),"I","",handle)

  NULLIFY ( c_in )
  NULLIFY ( c_out )

  IF ( PRESENT ( debug ) ) THEN
    test = debug
  ELSE
    test = .FALSE.
  END IF

!..check if grids are compatible
  IF ( pw1 % pw_grid % identifier /= pw2 % pw_grid % identifier ) THEN
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
    norm = 1._dbl / pw1 % pw_grid % ngpts
    out_space = RECIPROCALSPACE
  ELSE IF ( pw1 % in_space == RECIPROCALSPACE ) THEN
    dir = BWFFT
    norm = 1._dbl
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
      CASE ( "FW_R3DC3D" )
        ! nsize=size(pw1 % cr3d,1)*size(pw1 % cr3d,2)*size(pw1 % cr3d,3)
        ! CALL copy_rc(nsize,pw1%cr3d(1,1,1),pw2%cc3d(1,1,1))
        pw2 % cc3d = CMPLX ( pw1 % cr3d, 0._dbl, dbl )
        c_out => pw2 % cc3d
      CASE ( "FW_C3DC1D" )
        c_in => pw1 % cc3d
        IF (.NOT. ASSOCIATED (c_out)) ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out", &
            pw1 % pw_grid % ngpts )
      CASE ( "FW_R3DC1D" )
        IF (.NOT. ASSOCIATED (c_out)) ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out", &
            pw1 % pw_grid % ngpts )
        nsize=SIZE(pw1 % cr3d,1)*SIZE(pw1 % cr3d,2)*SIZE(pw1 % cr3d,3)
        l1 = LBOUND(pw1%cr3d,1)
        l2 = LBOUND(pw1%cr3d,2)
        l3 = LBOUND(pw1%cr3d,3)
        CALL copy_rc(nsize,pw1%cr3d(l1,l2,l3),c_out(1,1,1))
        ! c_out = CMPLX ( pw1 % cr3d, 0._dbl, dbl )
      CASE ( "BW_C3DC3D" )
        c_in => pw1 % cc3d
        c_out => pw2 % cc3d
      CASE ( "BW_C3DR3D" )
        c_in => pw1 % cc3d
        IF (.NOT. ASSOCIATED (c_out)) ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out", &
            pw1 % pw_grid % ngpts )
      CASE ( "BW_C1DC3D" )
        c_out => pw2 % cc3d
        IF ( test ) WRITE ( *,'(A)') "  PW_SCATTER : 3d -> 1d "
        CALL pw_scatter ( pw1, c_out )
      CASE ( "BW_C1DR3D" )
        IF (.NOT. ASSOCIATED (c_out)) ALLOCATE ( c_out( n(1), n(2), n(3) ), STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out", &
            pw1 % pw_grid % ngpts )
        IF ( test ) WRITE ( *,'(A)') "  PW_SCATTER : 3d -> 1d "
        CALL pw_scatter ( pw1, c_out )
    END SELECT

!..transform
    IF ( ASSOCIATED ( c_in ) ) THEN
      CALL fft3d ( dir, n, c_in, c_out, scale = norm, debug=test )
    ELSE
      CALL fft3d ( dir, n, c_out, scale = norm, debug=test )
    ENDIF

!..prepare output
    SELECT CASE ( mode )
      CASE DEFAULT
        CALL stop_program ( "fft_wrap_pw1pw2", &
                        "Illegal combination of in_use and in_space" )
      CASE ( "FW_C3DC3D" )
        !     nothing to do
      CASE ( "FW_R3DC3D" )
        !     nothing to do
      CASE ( "FW_C3DC1D" )
        !     gather results
        IF ( test ) WRITE ( *,'(A)') "  PW_GATHER : 3d -> 1d "
        CALL pw_gather ( pw2, c_out )
        DEALLOCATE ( c_out, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out" )
      CASE ( "FW_R3DC1D" )
        !     gather results
        IF ( test ) WRITE ( *,'(A)') "  PW_GATHER : 3d -> 1d "
        CALL pw_gather ( pw2, c_out )
        DEALLOCATE ( c_out, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out" )
      CASE ( "BW_C3DC3D" )
        !     nothing to do
      CASE ( "BW_C3DR3D" )
        !     use real part only
        IF ( test ) WRITE ( *,'(A)') "  REAL part "
        pw2 % cr3d = REAL ( c_out, dbl )
        DEALLOCATE ( c_out, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out" )
      CASE ( "BW_C1DC3D" )
        !     nothing to do
      CASE ( "BW_C1DR3D" )
        !     use real part only
        IF ( test ) WRITE ( *,'(A)') "  REAL part "
        ! pw2 % cr3d = REAL ( c_out, dbl )
        nsize=SIZE(pw2 % cr3d,1)*SIZE(pw2 % cr3d,2)*SIZE(pw2 % cr3d,3)
        l1 = LBOUND(pw2%cr3d,1)
        l2 = LBOUND(pw2%cr3d,2)
        l3 = LBOUND(pw2%cr3d,3)
        CALL copy_cr(nsize,c_out(1,1,1),pw2%cr3d(l1,l2,l3))
        
        DEALLOCATE ( c_out, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_out" )
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
    IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", &
      "grays", n(1) * nrays )
    grays = 0._dbl

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
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_in", &
            PRODUCT ( nloc ) )
        c_in = CMPLX ( pw1 % cr3d, 0._dbl, dbl )
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
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_in", &
            PRODUCT ( nloc ) )
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
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_in" )
      CASE ( "BW_C1DC3D" )
        ! nothing to do
      CASE ( "BW_C1DR3D" )
        IF ( test .AND. pw1%pw_grid%para%group_head ) &
            WRITE ( *,'(A)') "  Real part "
        pw2 % cr3d = REAL ( c_in, dbl )
        DEALLOCATE ( c_in, STAT = ierr )
        IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "c_in" )
    END SELECT

    DEALLOCATE ( grays, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "fft_wrap_pw1pw2", "grays" )

  END IF

!..update the space tag for pw2
  pw2 % in_space = out_space

  IF ( test .AND. pw1 % pw_grid % para % group_head ) THEN
      WRITE ( *,'(A)') " End of FFT Protocol "
  END IF
  CALL timestop(0.0_dbl,handle) 
  CALL timestop(0.0_dbl,handle2) 

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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

  WRITE (unit=unit_nr, fmt="(' pw_grid%identifier=',i8,/,' }')",iostat=iostat)&
       pw%pw_grid%identifier

END SUBROUTINE pw_write
!***************************************************************************

FUNCTION pw_compatible ( grida, gridb ) RESULT ( compat )
    TYPE(pw_grid_type), INTENT(IN)           :: grida, gridb
    LOGICAL                                  :: compat

   compat = .FALSE.
   IF ( grida%identifier == gridb%identifier ) THEN
      compat = .TRUE.
   ELSE IF ( grida%reference == gridb%identifier ) THEN
      compat = .TRUE.
   ELSE IF ( gridb%reference == grida%identifier ) THEN
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
        CALL pw_deallocate(pw)
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
!!     Subroutine pw_create(pw, pw_grid, use_data, in_space, mass_flag,&
!!         cr3d_ptr, error)
!!       Type(pw_type), Pointer:: pw
!!       Type(pw_grid_type), Pointer:: pw_grid
!!       Integer, Intent (IN):: use_data
!!       Integer, Intent (IN), Optional:: in_space
!!       Logical, Intent (IN), Optional:: mass_flag
!!       Real(Kind=dbl), Dimension(:,:,:), Optional, Pointer:: cr3d_ptr
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
!!     - mass_flag: if mass should be allocated
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
SUBROUTINE pw_create(pw, pw_grid, use_data, in_space, mass_flag, cr3d_ptr, error)
    TYPE(pw_type), POINTER                   :: pw
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, INTENT(in)                      :: use_data
    INTEGER, INTENT(in), OPTIONAL            :: in_space
    LOGICAL, INTENT(in), OPTIONAL            :: mass_flag
    REAL(kind=dbl), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: cr3d_ptr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  ALLOCATE(pw,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(cr3d_ptr)) THEN
        CPPrecondition(.NOT.ASSOCIATED(cr3d_ptr).OR.(use_data==REALDATA3D),cp_failure_level,routineP,error,failure)
        CALL pw_allocate(pw, pw_grid, use_data=use_data, &
             mass_flag=mass_flag,cr3d_ptr=cr3d_ptr)
     ELSE
        CALL pw_allocate(pw, pw_grid, use_data=use_data, &
             mass_flag=mass_flag)
     END IF

     pw%in_space=0
     IF (PRESENT(in_space)) pw%in_space=in_space
  END IF
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
!!       Real(Kind=dbl), Intent (OUT), Optional:: maxdiff
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
    REAL(kind=dbl), INTENT(out), OPTIONAL    :: maxdiff
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_compare_debug', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, unit_nr
    LOGICAL                                  :: failure
    REAL(kind=dbl)                           :: diff, mdiff
    TYPE(cp_logger_type), POINTER            :: logger

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pw1),cp_failure_level,routineP,error,failure)
  CPPrecondition(ASSOCIATED(pw2),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw1%in_space==pw2%in_space,cp_warning_level,routineP,error,failure)
  CPPrecondition(pw1%in_use==pw2%in_use,cp_warning_level,routineP,error,failure)
  CALL cp_assert(ALL(pw1%pw_grid%bounds_local==pw2%pw_grid%bounds_local),&
       cp_failure_level,cp_assertion_failed,routineP,&
       "wrong pw distribution",error,failure)
  IF (.NOT. failure) THEN
     logger => cp_error_get_logger(error)
     unit_nr=-1
     mdiff=0.0_dbl
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
                   diff,i,j,k,pw1%cc(i),pw2%cc(i)
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
     IF (PRESENT(maxdiff)) maxdiff=HUGE(0.0_dbl)
  END IF
END SUBROUTINE pw_compare_debug
!***************************************************************************

!!****f* pw_types/pw_prolongate_NoPBC [1.0] *
!!
!!   NAME
!!     pw_prolongate_NoPBC
!!
!!   FUNCTION
!!     Prolongates the values of a pw_type coarse grid to a pw_type fine grid
!!     uses simple linear interpolation. First Order.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE pw_prolongate_NoPBC(pw_coarse, pw_fine, error)
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(pw_type), POINTER       :: pw_coarse
    TYPE(pw_type), POINTER       :: pw_fine
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'pw_prolongate_NoPBC', &
         routineP = moduleN//':'//routineN
    LOGICAL                                  :: failure
    INTEGER                                  :: handle, i, j, k, stat
    INTEGER, DIMENSION(3)                    :: lb, n, ub
    REAL(KIND=dbl), DIMENSION(:,:,:), ALLOCATABLE :: tgr 


    CALL timeset ( 'pw_prolongate_NoPBC', 'I', '', handle )

    ! only for the easiest case 
    IF ( ANY((pw_coarse % pw_grid % npts -1)* 2 .NE. (pw_fine % pw_grid % npts -1))) THEN
       CALL stop_program("pw_prolongate_NoPBC","sorry too difficult ")
    ENDIF
    !WRITE(*,*)'coarse grid points  : ',pw_coarse % pw_grid% npts
    !WRITE(*,*)'       upper bounds : ',pw_coarse % pw_grid% bounds(2,:)
    !WRITE(*,*)'       lower bounds : ',pw_coarse % pw_grid% bounds(1,:)
    !WRITE(*,*)'fine   grid points: ',  pw_fine % pw_grid%npts
    !WRITE(*,*)'       upper bounds : ',pw_fine % pw_grid%bounds(2,:)
    !WRITE(*,*)'       lower bounds : ',pw_fine % pw_grid%bounds(1,:)
    ! we use the bounds of the smaller grid
    n   = pw_coarse % pw_grid %npts
    lb  = pw_fine % pw_grid % bounds(1,:)
    ub  = pw_fine % pw_grid % bounds(2,:)

    ALLOCATE( tgr( lb(1):ub(1),&
                   lb(2):ub(2),&
                   lb(3):ub(3)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    tgr = 0.0_dbl
    lb  = pw_coarse % pw_grid % bounds(1,:)
    ub  = pw_coarse % pw_grid % bounds(2,:)

    ! Copy the indentical points
    DO k=lb(3),ub(3)
       DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
              tgr (2*i,2*j,2*k ) = pw_coarse % cr3d (i,j,k)
          ENDDO
       ENDDO
    ENDDO

    ! Put values on the x edges
    DO k=lb(3),ub(3)
       DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)-1
             tgr(2*i+1,2*j,2*k ) = 0.5_dbl* (tgr(2*i,   2*j, 2*k ) + &
                                             tgr(2*i+2, 2*j, 2*k ))
          ENDDO
       ENDDO
    ENDDO
    ! Put values on the y edges
    DO k=lb(3),ub(3)
       DO j=lb(2),ub(2)-1
          DO i=lb(1),ub(1)
             tgr(2*i,2*j+1,2*k ) = 0.5_dbl* (tgr(2*i, 2*j,   2*k ) +  &
                                             tgr(2*i, 2*j+2, 2*k ))
          ENDDO
       ENDDO
    ENDDO
    ! Put values on the z edges
    DO k=lb(3),ub(3)-1
       DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
             tgr(2*i,2*j,2*k+1 ) = 0.5_dbl* (tgr(2*i, 2*j, 2*k  ) +  &
                                             tgr(2*i, 2*j, 2*k+2 ))
          ENDDO
       ENDDO
    ENDDO

    ! Put values on the xy faces
    DO k=lb(3),ub(3)
       DO j=lb(2),ub(2)-1
          DO i=lb(1),ub(1)-1
             tgr(2*i+1,2*j+1,2*k ) = 0.5_dbl* (tgr(2*i,   2*j+1, 2*k ) + &  
                                               tgr(2*i+2, 2*j+1, 2*k ))
          ENDDO
       ENDDO
    ENDDO
    ! Put values on the xz faces
    DO k=lb(3),ub(3)-1
       DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)-1
             tgr(2*i+1,2*j,2*k+1 ) = 0.5_dbl* (tgr(2*i,   2*j, 2*k+1 ) + &
                                               tgr(2*i+2, 2*j, 2*k+1 ))
          ENDDO
       ENDDO
    ENDDO
    ! Put values on the yz faces
    DO k=lb(3),ub(3)-1
       DO j=lb(2),ub(2)-1
          DO i=lb(1),ub(1)
             tgr(2*i,2*j+1,2*k+1 ) = 0.5_dbl* (tgr(2*i, 2*j,   2*k+1 ) + &
                                               tgr(2*i, 2*j+2, 2*k+1 ))
          ENDDO
       ENDDO
    ENDDO
    ! Put values on the center
    DO k=lb(3),ub(3)-1
       DO j=lb(2),ub(2)-1
          DO i=lb(1),ub(1)-1
             tgr(2*i+1,2*j+1,2*k+1 ) = 0.5_dbl* (tgr(2*i ,  2*j+1, 2*k+1 ) + &
                                                 tgr(2*i+2, 2*j+1, 2*k+1 ))
          ENDDO
       ENDDO
    ENDDO
    pw_fine % cr3d = pw_fine % cr3d + tgr
    DEALLOCATE( tgr, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop ( 0.0_dbl, handle )    

  END SUBROUTINE pw_prolongate_NoPBC

!!****f* pw_types/pw_inject_NoPBC [1.0] *
!!
!!   NAME
!!     pw_inject_NoPBC
!!
!!   FUNCTION
!!     Injectes the values of a fine pw_type grid to a coarse pw_type grid
!!     uses simple linear interpolation, and is the 'inverse' of the 
!!     rs_prolongate_NoPBC
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [tlaino]
!!
!!**************************************************************************
  SUBROUTINE pw_inject_NoPBC(pw_fine, pw_coarse, error)
    ! Arguments
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error
    TYPE(pw_type), POINTER        :: pw_fine
    TYPE(pw_type), POINTER        :: pw_coarse
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'pw_inject_NoPBC', &
         routineP = moduleN//':'//routineN
    LOGICAL                                   :: failure
    INTEGER                                   :: handle, i, j, k, stat
    INTEGER, DIMENSION(3)                     :: lb, ub
    REAL(KIND=dbl)                            :: dum
    REAL(KIND=dbl), DIMENSION(:,:,:), ALLOCATABLE :: tgr 

    CALL timeset ( 'pw_inject_NoPBC', 'I', '', handle )

    ! only for the easiest case
    IF ( ANY((pw_coarse % pw_grid % npts - 1 )* 2 .NE. (pw_fine % pw_grid % npts -1))) THEN
       CALL stop_program("pw_inject_NoPBC","sorry too difficult ")
    ENDIF
    lb  = pw_fine % pw_grid % bounds(1,:)
    ub  = pw_fine % pw_grid % bounds(2,:)

    ALLOCATE( tgr( lb(1):ub(1),&
                   lb(2):ub(2),&
                   lb(3):ub(3)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    tgr = pw_fine % cr3d
    lb  = pw_coarse % pw_grid % bounds(1,:)
    ub  = pw_coarse % pw_grid % bounds(2,:)

    ! Values from the center
    DO k=lb(3),ub(3)-1
       DO j=lb(2),ub(2)-1
          DO i=lb(1),ub(1)-1
             dum = 0.5_dbl * tgr(2*i+1, 2*j+1, 2*k+1 )
             tgr(2*i  , 2*j+1, 2*k+1 ) = tgr(2*i ,  2*j+1, 2*k+1 ) + dum
             tgr(2*i+2, 2*j+1, 2*k+1 ) = tgr(2*i+2, 2*j+1, 2*k+1 ) + dum
          ENDDO
       ENDDO
    ENDDO
    ! Values from the yz faces
    DO k=lb(3),ub(3)-1
       DO j=lb(2),ub(2)-1
          DO i=lb(1),ub(1)
             dum = 0.5_dbl * tgr(2*i  , 2*j+1, 2*k+1 )
             tgr(2*i  , 2*j  , 2*k+1 ) = tgr(2*i  ,2*j  , 2*k+1 ) + dum
             tgr(2*i  , 2*j+2, 2*k+1 ) = tgr(2*i  ,2*j+2, 2*k+1 ) + dum
          ENDDO
       ENDDO
    ENDDO
    ! Values from the xz faces
    DO k=lb(3),ub(3)-1
       DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)-1
             dum = 0.5_dbl * tgr(2*i+1, 2*j  , 2*k+1 )
             tgr(2*i  , 2*j  , 2*k+1 ) = tgr(2*i  , 2*j  , 2*k+1 ) + dum
             tgr(2*i+2, 2*j  , 2*k+1 ) = tgr(2*i+2, 2*j  , 2*k+1 ) + dum
          ENDDO
       ENDDO
    ENDDO
    ! Values from the xy faces
    DO k=lb(3),ub(3)
       DO j=lb(2),ub(2)-1
          DO i=lb(1),ub(1)-1
             dum = 0.5_dbl * tgr(2*i+1, 2*j+1, 2*k   )
             tgr(2*i  , 2*j+1, 2*k   ) = tgr(2*i  , 2*j+1, 2*k   ) + dum
             tgr(2*i+2, 2*j+1, 2*k   ) = tgr(2*i+2, 2*j+1, 2*k   ) + dum
          ENDDO
       ENDDO
    ENDDO
    ! Values from the z edges
    DO k=lb(3),ub(3)-1
       DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
             dum = 0.5_dbl * tgr(2*i  , 2*j  , 2*k+1 )
             tgr(2*i  , 2*j  , 2*k   ) = tgr(2*i  , 2*j  , 2*k   ) + dum
             tgr(2*i  , 2*j  , 2*k+2 ) = tgr(2*i  , 2*j  , 2*k+2 ) + dum
          ENDDO
       ENDDO
    ENDDO
    ! Values from the y edges
    DO k=lb(3),ub(3)
       DO j=lb(2),ub(2)-1
          DO i=lb(1),ub(1)
             dum = 0.5_dbl * tgr(2*i  , 2*j+1, 2*k   )
             tgr(2*i  , 2*j  , 2*k   ) = tgr(2*i  , 2*j  , 2*k   ) + dum
             tgr(2*i  , 2*j+2, 2*k   ) = tgr(2*i  , 2*j+2, 2*k   ) + dum
          ENDDO
       ENDDO
    ENDDO
    ! Values from the x edges
    DO k=lb(3),ub(3)
       DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)-1
             dum = 0.5_dbl * tgr(2*i+1, 2*j  , 2*k   )
             tgr(2*i  , 2*j  , 2*k   ) = tgr(2*i  , 2*j  , 2*k   ) + dum
             tgr(2*i+2, 2*j  , 2*k   ) = tgr(2*i+2, 2*j  , 2*k   ) + dum
          ENDDO
       ENDDO
    ENDDO
    ! Copy the indentical points
    DO k=lb(3),ub(3)
       DO j=lb(2),ub(2)
          DO i=lb(1),ub(1)
             pw_coarse % cr3d(i,j,k) = pw_coarse % cr3d(i,j,k) + tgr(2*i,2*j,2*k )
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE( tgr, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop ( 0.0_dbl, handle )
    
  END SUBROUTINE pw_inject_NoPBC

END MODULE pw_types

!******************************************************************************
