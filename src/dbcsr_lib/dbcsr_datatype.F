!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR data types
!> \author  Urban Borstnik
!> \date    2009-05-12
!> \version 0.95
!>
!> <b>Modification history:</b>
!> - Created 2009-05-12
! *****************************************************************************
MODULE dbcsr_datatype
  PUBLIC :: dbcsr_type_real_4, dbcsr_type_real_8,&
            dbcsr_type_complex_4, dbcsr_type_complex_8,&
            dbcsr_type_real_default, dbcsr_type_complex_default,&
            dbcsr_type_real_4_2d, dbcsr_type_real_8_2d,&
            dbcsr_type_complex_4_2d, dbcsr_type_complex_8_2d,&
            dbcsr_type_int_4, dbcsr_type_int_8


!> Type definitions:
!> * bit 0: always 1
!> * bit 1: single (0: 4) vs. double (1: 1)
!> * bit 2: real (0) vs. complex (1)
!> * bit 3: dimension (0: 1, 1: 2)
!> * bit 4: floating point (0) or integer type (1)
  INTEGER, PARAMETER          :: dbcsr_type_real_4 = 1         !001
  INTEGER, PARAMETER          :: dbcsr_type_real_8 = 3         !011
  INTEGER, PARAMETER          :: dbcsr_type_complex_4 = 5      !101
  INTEGER, PARAMETER          :: dbcsr_type_complex_8 = 7      !111

  INTEGER, PARAMETER          :: dbcsr_type_real_4_2d = 9      !1001
  INTEGER, PARAMETER          :: dbcsr_type_real_8_2d = 11     !1011
  INTEGER, PARAMETER          :: dbcsr_type_complex_4_2d = 13  !1101
  INTEGER, PARAMETER          :: dbcsr_type_complex_8_2d = 15  !1111

  INTEGER, PARAMETER          :: dbcsr_type_int_4 = 17         !10001
  INTEGER, PARAMETER          :: dbcsr_type_int_8 = 19         !10011

  INTEGER, PARAMETER          :: dbcsr_type_real_default = dbcsr_type_real_8
  INTEGER, PARAMETER          :: dbcsr_type_complex_default = dbcsr_type_complex_8

  ! Different method for dbcsr_filter
  INTEGER, PARAMETER            :: dbcsr_filter_frobenius = 1

  ! Different norm for dbcsr_norm
  INTEGER, PARAMETER            :: dbcsr_norm_frobenius = 1
  INTEGER, PARAMETER            :: dbcsr_norm_maxabsnorm = 2
  INTEGER, PARAMETER            :: dbcsr_norm_gershgorin = 3
  INTEGER, PARAMETER            :: dbcsr_norm_column = 4



END MODULE dbcsr_datatype
