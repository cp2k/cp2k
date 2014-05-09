!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines the basic variable types
!> \note
!>      Data type definitions; tested on:
!>          - IBM AIX xlf90
!>          - SGI IRIX  f90
!>          - CRAY T3E  f90
!>          - DEC ALPHA f90
!>          - NAG_F90
!>          - SUN
!>          - HITACHI
!> \par History
!>      Adapted for CP2K by JGH
!> \author Matthias Krack
! *****************************************************************************
MODULE kinds

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sp, dp, print_kind_info, dp_size, sp_size, int_size
  PUBLIC :: int_1, int_2, int_4, int_8, int_1_size, int_2_size, int_4_size, int_8_size
  PUBLIC :: real_4, real_8, real_4_size, real_8_size
  PUBLIC :: default_string_length, default_path_length, max_line_length

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
  ! we rely on this (libraries) but do not check this
  INTEGER, PARAMETER :: dp_size  = 8,&
                        int_size = BIT_SIZE(0)/8,&
                        sp_size  = 4

  INTEGER, PARAMETER :: real_4 = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER :: real_8 = SELECTED_REAL_KIND ( 14, 200 )
  INTEGER, PARAMETER :: real_4_size = 4
  INTEGER, PARAMETER :: real_8_size = 8


  INTEGER, PARAMETER :: int_1  = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: int_1_size = BIT_SIZE(INT(0, int_1))/8

  INTEGER, PARAMETER :: int_2  = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: int_2_size = BIT_SIZE(INT(0, int_2))/8

  INTEGER, PARAMETER :: int_4 = SELECTED_INT_KIND(5)
  INTEGER, PARAMETER :: int_4_size = BIT_SIZE(INT(0, int_4))/8

  INTEGER, PARAMETER :: int_8 = SELECTED_INT_KIND(10)
  INTEGER, PARAMETER :: int_8_size = BIT_SIZE(INT(0, int_8))/8


  INTEGER, PARAMETER :: default_string_length=80
  INTEGER, PARAMETER :: default_path_length=1024
  INTEGER, PARAMETER :: max_line_length = 2*default_path_length
  CHARACTER(LEN=1),PARAMETER, PUBLIC :: default_blank_character(2) = (/" ",CHAR(9)/)

CONTAINS

! *****************************************************************************
!> \brief Print informations about the used data types.
!> \par History
!>      Adapted by JGH for Cp2k
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE print_kind_info ( iw )

    INTEGER, INTENT(IN)                      :: iw

    WRITE ( iw, '( /, T2, A )' ) 'DATA TYPE INFORMATION:'

    WRITE ( iw, '( /,T2,A,T79,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E14.8) )' ) &
         'REAL: Data type name:', 'dp', '      Kind value:', KIND ( 0.0_dp ), &
         '      Precision:', PRECISION ( 0.0_dp ), &
         '      Smallest non-negligible quantity relative to 1:', &
         EPSILON ( 0.0_dp ), &
         '      Smallest positive number:', TINY ( 0.0_dp ), &
         '      Largest representable number:', HUGE ( 0.0_dp )
    WRITE ( iw, '( /,T2,A,T79,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E14.8) )' ) &
         '      Data type name:', 'sp', '      Kind value:', KIND ( 0.0_sp ), &
         '      Precision:', PRECISION ( 0.0_sp ), &
         '      Smallest non-negligible quantity relative to 1:', &
         EPSILON ( 0.0_sp ), &
         '      Smallest positive number:', TINY ( 0.0_sp ), &
         '      Largest representable number:', HUGE ( 0.0_sp )
    WRITE ( iw, '( /,T2,A,T72,A,4(/,T2,A,T61,I20) )' ) &
         'INTEGER: Data type name:', '(default)', '         Kind value:', &
         KIND ( 0 ), &
         '         Bit size:', BIT_SIZE ( 0 ), &
         '         Largest representable number:', HUGE ( 0 )
    WRITE ( iw, '( /,T2,A,T72,A,/,T2,A,T75,I6,/ )' ) &
         'LOGICAL: Data type name:', '(default)', &
         '         Kind value:', KIND ( .TRUE. )
    WRITE ( iw, '( /,T2,A,T72,A,/,T2,A,T75,I6,/ )' ) &
         'CHARACTER: Data type name:', '(default)', &
         '           Kind value:', KIND ( 'C' )

  END SUBROUTINE print_kind_info

END MODULE kinds

