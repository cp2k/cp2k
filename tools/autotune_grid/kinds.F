MODULE kinds
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: sp, dp, dp_size, sp_size, int_size, int_8
  PUBLIC :: default_string_length, default_path_length

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND ( 6, 30 )
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND ( 14, 200 )
  ! we rely on this (libraries) but do not check this
  INTEGER, PARAMETER :: dp_size  = 8,&
                        int_size = BIT_SIZE(0)/8,&
                        sp_size  = 4

  ! this int holds more than the normal 4byte ints
  ! on standard machines it ought to be an 8byte int but this is not guaranteed
  INTEGER, PARAMETER :: int_8 = SELECTED_INT_KIND(10)

  INTEGER, PARAMETER :: default_string_length=80
  INTEGER, PARAMETER :: default_path_length=250
END MODULE
