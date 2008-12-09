MODULE mpfr
  USE ISO_C_BINDING

  IMPLICIT NONE
  

  INTEGER, PARAMETER               :: MAX_CHAR = 10000
  INTEGER, PARAMETER, PUBLIC       :: dp=KIND(0.0D0)

  INTEGER, PARAMETER, PUBLIC       :: GMP_RNDN=0,&
                                      GMP_RNDZ=1,&
                                      GMP_RNDU=2,&
                                      GMP_RNDD=3,&
                                      GMP_RND_MAX=4,&
                                      GMP_RNDNA=-1
                                           
  TYPE, BIND(C)    :: mpfr_type
    INTEGER(C_SHORT)            :: mpfr_prec
    INTEGER(C_LONG)             :: mpfr_sign
    INTEGER(C_LONG)             :: mpfr
    TYPE(C_PTR)                 :: mpfr_d
  END TYPE mpfr_type

  INTERFACE
    SUBROUTINE mpfr_init2(value, precision) BIND(C, name="mpfr_init2")
      IMPORT
      TYPE(mpfr_type)             :: value
      INTEGER(C_SHORT), VALUE     :: precision
    END SUBROUTINE mpfr_init2
  
    SUBROUTINE mpfr_init(value) BIND(C, name="mpfr_init")
      IMPORT
      TYPE(mpfr_type)             :: value
    END SUBROUTINE mpfr_init
  
    SUBROUTINE mpfr_set_default_precision(precision) BIND(C, name="mpfr_set_default_prec")
      IMPORT
      INTEGER(C_SHORT), VALUE     :: precision
    END SUBROUTINE mpfr_set_default_precision
  
    FUNCTION mpfr_get_default_precision() BIND(C, name="mpfr_get_default_prec")
      IMPORT
      INTEGER(C_SHORT)            :: mpfr_get_default_precision
    END FUNCTIOn mpfr_get_default_precision
  
    FUNCTION mpfr_get_precision(variable) BIND(C, name="mpfr_get_prec")
      IMPORT
      INTEGER(C_SHORT)            :: mpfr_get_precision
      TYPE(mpfr_type)             :: variable
    END FUNCTION mpfr_get_precision
  
    FUNCTION mpfr_set_d(variable,real_value,rounding) BIND(C, name="mpfr_set_d")
      IMPORT
      INTEGER(C_INT)              :: mpfr_set_d
      TYPE(mpfr_type)             :: variable
      REAL(C_DOUBLE), VALUE       :: real_value
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_set_d

    FUNCTION mpfr_set_str(variable,str,base,rounding) BIND(C, name="mpfr_set_str")
      IMPORT
      INTEGER(C_INT)              :: mpfr_set_str
      TYPE(mpfr_type)             :: variable
      CHARACTER(C_CHAR)           :: str
      INTEGER(C_INT), VALUE       :: base
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_set_str
    
    FUNCTION mpfr_strtofr(variable,str1,str2,base,rounding) BIND(C,name="mpfr_strtofr")
      IMPORT
      INTEGER(C_INT)              :: mpfr_strtofr
      TYPE(mpfr_type)             :: variable
      CHARACTER(C_CHAR)           :: str1
      CHARACTER(C_CHAR),VALUE     :: str2
      INTEGER(C_INT), VALUE       :: base
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_strtofr
  
    PURE FUNCTION mpfr_get_d(variable,rounding) BIND(C, name="mpfr_get_d")
      IMPORT
      REAL(C_DOUBLE)              :: mpfr_get_d
      TYPE(mpfr_type), INTENT(IN) :: variable
      INTEGER(C_INT), VALUE, INTENT(IN) :: rounding
    END FUNCTION mpfr_get_d

    FUNCTION mpfr_cmp(op1,op2) BIND(C,name="mpfr_cmp")
      IMPORT
      INTEGER(C_INT)              :: mpfr_cmp
      TYPE(mpfr_type)             :: op1,op2
    END FUNCTION mpfr_cmp
  
    FUNCTION mpfr_add(result,op1,op2,rounding) BIND(C, name="mpfr_add")
      IMPORT
      INTEGER(C_INT)              :: mpfr_add
      TYPE(mpfr_type)             :: result,op1,op2
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_add
  
    FUNCTION mpfr_sub(result,op1,op2,rounding) BIND(C, name="mpfr_sub")
      IMPORT
      INTEGER(C_INT)              :: mpfr_sub
      TYPE(mpfr_type)             :: result,op1,op2
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_sub 
  
    FUNCTION mpfr_mul_ui(result,op1,int,rounding) BIND(C, name="mpfr_mul_ui")
      IMPORT
      INTEGER(C_INT)              :: mpfr_mul_ui
      TYPE(mpfr_type)             :: result,op1
      INTEGER(C_INT), VALUE       :: int
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_mul_ui
  
    FUNCTION mpfr_mul(result,op1,op2,rounding) BIND(C, name="mpfr_mul")
      IMPORT
      INTEGER(C_INT)              :: mpfr_mul
      TYPE(mpfr_type)             :: result,op1,op2
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_mul
  
    FUNCTION mpfr_div(result,op1,op2,rounding) BIND(C, name="mpfr_div")
      IMPORT
      INTEGER(C_INT)              :: mpfr_div
      TYPE(mpfr_type)             :: result,op1,op2
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_div
  
    FUNCTION mpfr_pow(result,op1,op2,rounding) BIND(C, name="mpfr_pow")
      IMPORT
      INTEGER(C_INT)              :: mpfr_pow
      TYPE(mpfr_type)             :: result,op1,op2
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_pow
  
    FUNCTION mpfr_pow_si(result,op1,op2,rounding) BIND(C, name="mpfr_pow_si")
      IMPORT
      INTEGER(C_INT)              :: mpfr_pow_si
      TYPE(mpfr_type)             :: result,op1
      INTEGER(C_SHORT)            :: op2
      INTEGER(C_INT), VALUE       :: rounding
    END FUNCTION mpfr_pow_si
  
    FUNCTION mpfr_dump(variable) BIND(C, name="mpfr_dump")
      IMPORT
      INTEGER(C_INT)              :: mpfr_dump
      TYPE(mpfr_type)             :: variable
    END FUNCTION mpfr_dump
  
    SUBROUTINE mpfr_clear(variable) BIND(C, name="mpfr_clear")
      IMPORT
      TYPE(mpfr_type)             :: variable
    END SUBROUTINE mpfr_clear
  
    FUNCTION mpfr_get_str(str,exp,base,n,variable,rounding) BIND(C,name="mpfr_get_str")
      IMPORT
      TYPE(C_PTR)                :: mpfr_get_str
      CHARACTER(C_CHAR)          :: str
      INTEGER(C_INT)             :: exp
      INTEGER(C_INT), VALUE      :: base
      INTEGER(C_SHORT), VALUE    :: n
      TYPE(mpfr_type)            :: variable
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_get_str
  
    FUNCTION mpfr_log(result,op,rounding) BIND(C,name="mpfr_log")
      IMPORT
      INTEGER(C_INT)             :: mpfr_log
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_log
  
    FUNCTION mpfr_log2(result,op,rounding) BIND(C,name="mpfr_log2")
      IMPORT
      INTEGER(C_INT)             :: mpfr_log2
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_log2
  
    FUNCTION mpfr_log10(result,op,rounding) BIND(C,name="mpfr_log10")
      IMPORT
      INTEGER(C_INT)             :: mpfr_log10
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_log10
  
    FUNCTION mpfr_exp(result,op,rounding) BIND(C,name="mpfr_exp")
      IMPORT
      INTEGER(C_INT)             :: mpfr_exp
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_exp
  
    FUNCTION mpfr_exp2(result,op,rounding) BIND(C,name="mpfr_exp2")
      IMPORT
      INTEGER(C_INT)             :: mpfr_exp2
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_exp2
  
    FUNCTION mpfr_exp10(result,op,rounding) BIND(C,name="mpfr_exp10")
      IMPORT
      INTEGER(C_INT)             :: mpfr_exp10
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_exp10
  
    FUNCTION mpfr_cos(result,op,rounding) BIND(C,name="mpfr_cos")
      IMPORT
      INTEGER(C_INT)             :: mpfr_cos
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_cos
  
    FUNCTION mpfr_sin(result,op,rounding) BIND(C,name="mpfr_sin")
      IMPORT
      INTEGER(C_INT)             :: mpfr_sin
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_sin
   
    FUNCTION mpfr_tan(result,op,rounding) BIND(C,name="mpfr_tan")
      IMPORT
      INTEGER(C_INT)             :: mpfr_tan
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_tan
  
    FUNCTION mpfr_sec(result,op,rounding) BIND(C,name="mpfr_sec")
      IMPORT
      INTEGER(C_INT)             :: mpfr_sec
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_sec
   
    FUNCTION mpfr_csc(result,op,rounding) BIND(C,name="mpfr_csc")
      IMPORT
      INTEGER(C_INT)             :: mpfr_csc
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_csc
  
    FUNCTION mpfr_cot(result,op,rounding) BIND(C,name="mpfr_cot")
      IMPORT
      INTEGER(C_INT)             :: mpfr_cot
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_cot
  
    FUNCTION mpfr_acos(result,op,rounding) BIND(C,name="mpfr_acos")
      IMPORT
      INTEGER(C_INT)             :: mpfr_acos
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_acos
  
    FUNCTION mpfr_asin(result,op,rounding) BIND(C,name="mpfr_asin")
      IMPORT
      INTEGER(C_INT)             :: mpfr_asin
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_asin
  
    FUNCTION mpfr_atan(result,op,rounding) BIND(C,name="mpfr_atan")
      IMPORT
      INTEGER(C_INT)             :: mpfr_atan
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_atan
  
    FUNCTION mpfr_atan2(result,x,y,rounding) BIND(C,name="mpfr_atan2")
      IMPORT
      INTEGER(C_INT)             :: mpfr_atan2
      TYPE(mpfr_type)            :: result, x,y
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_atan2
  
    FUNCTION mpfr_cosh(result,op,rounding) BIND(C,name="mpfr_cosh")
      IMPORT
      INTEGER(C_INT)             :: mpfr_cosh
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_cosh
  
    FUNCTION mpfr_sinh(result,op,rounding) BIND(C,name="mpfr_sinh")
      IMPORT
      INTEGER(C_INT)             :: mpfr_sinh
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_sinh
  
    FUNCTION mpfr_tanh(result,op,rounding) BIND(C,name="mpfr_tanh")
      IMPORT
      INTEGER(C_INT)             :: mpfr_tanh
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_tanh
  
    FUNCTION mpfr_sech(result,op,rounding) BIND(C,name="mpfr_sech")
      IMPORT
      INTEGER(C_INT)             :: mpfr_sech
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_sech
  
    FUNCTION mpfr_csch(result,op,rounding) BIND(C,name="mpfr_csch")
      IMPORT
      INTEGER(C_INT)             :: mpfr_csch
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_csch
  
    FUNCTION mpfr_coth(result,op,rounding) BIND(C,name="mpfr_coth")
      IMPORT
      INTEGER(C_INT)             :: mpfr_coth
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_coth
  
    FUNCTION mpfr_acosh(result,op,rounding) BIND(C,name="mpfr_acosh")
      IMPORT
      INTEGER(C_INT)             :: mpfr_acosh
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_acosh
  
    FUNCTION mpfr_asinh(result,op,rounding) BIND(C,name="mpfr_asinh")
      IMPORT
      INTEGER(C_INT)             :: mpfr_asinh
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_asinh
  
    FUNCTION mpfr_atanh(result,op,rounding) BIND(C,name="mpfr_atanh")
      IMPORT
      INTEGER(C_INT)             :: mpfr_atanh
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_atanh
  
    FUNCTION mpfr_eint(x,y,rounding) BIND(C,name="mpfr_eint")
      IMPORT
      INTEGER(C_INT)             :: mpfr_eint
      TYPE(mpfr_type)            :: x,y
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_eint
  
    FUNCTION mpfr_gamma(result,op,rounding) BIND(C,name="mpfr_gamma")
      IMPORT
      INTEGER(C_INT)             :: mpfr_gamma
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_gamma
  
    FUNCTION mpfr_lngamma(result,op,rounding) BIND(C,name="mpfr_lngamma")
      IMPORT
      INTEGER(C_INT)             :: mpfr_lngamma
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_lngamma
  
    FUNCTION mpfr_erf(result,op,rounding) BIND(C,name="mpfr_erf")
      IMPORT
      INTEGER(C_INT)             :: mpfr_erf
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_erf
  
    FUNCTION mpfr_erfc(result,op,rounding) BIND(C,name="mpfr_erfc")
      IMPORT
      INTEGER(C_INT)             :: mpfr_erfc
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_erfc
  
    FUNCTION mpfr_bessel_j0(result,op,rounding) BIND(C,name="mpfr_j0")
      IMPORT
      INTEGER(C_INT)             :: mpfr_bessel_j0
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_bessel_j0
  
    FUNCTION mpfr_bessel_j1(result,op,rounding) BIND(C,name="mpfr_j1")
      IMPORT
      INTEGER(C_INT)             :: mpfr_bessel_j1
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_bessel_j1
  
    FUNCTION mpfr_bessel_y0(result,op,rounding) BIND(C,name="mpfr_y0")
      IMPORT
      INTEGER(C_INT)             :: mpfr_bessel_y0
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_bessel_y0
  
    FUNCTION mpfr_bessel_y1(result,op,rounding) BIND(C,name="mpfr_y1")
      IMPORT
      INTEGER(C_INT)             :: mpfr_bessel_y1
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_bessel_y1
  
    FUNCTION mpfr_const_log2(result,rounding) BIND(C,name="mpfr_const_log2")
      IMPORT
      INTEGER(C_INT)             :: mpfr_const_log2
      TYPE(mpfr_type)            :: result
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_const_log2
  
    FUNCTION mpfr_const_pi(result,rounding) BIND(C,name="mpfr_const_pi")
      IMPORT
      INTEGER(C_INT)             :: mpfr_const_pi
      TYPE(mpfr_type)            :: result
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_const_pi
  
    FUNCTION mpfr_const_euler(result,rounding) BIND(C,name="mpfr_const_euler")
      IMPORT
      INTEGER(C_INT)             :: mpfr_const_euler
      TYPE(mpfr_type)            :: result
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_const_euler

    FUNCTION mpfr_sqrt(result,op,rounding) BIND(C,name="mpfr_sqrt")
      IMPORT
      INTEGER(C_INT)             :: mpfr_sqrt
      TYPE(mpfr_type)            :: result, op
      INTEGER(C_INT), VALUE      :: rounding
    END FUNCTION mpfr_sqrt
  END INTERFACE

END MODULE mpfr

MODULE mpfr_ops
  USE mpfr

  INTERFACE OPERATOR (+)
    MODULE PROCEDURE mpfr_addition_mp_mp,&
                     mpfr_addition_mp_real,&
                     mpfr_addition_real_mp,&
                     mpfr_addition_mp_int,&
                     mpfr_addition_int_mp
  END INTERFACE
  INTERFACE OPERATOR (-)
    MODULE PROCEDURE mpfr_subtraction_mp_mp,&
                     mpfr_subtraction_mp_real,&
                     mpfr_subtraction_real_mp,&
                     mpfr_subtraction_int_mp,&
                     mpfr_subtraction_mp_int,&
                     mpfr_minus
  END INTERFACE
  INTERFACE OPERATOR (*)
    MODULE PROCEDURE mpfr_multiplication_mp_mp,&
                     mpfr_multiplication_real_mp,&
                     mpfr_multiplication_mp_real,&
                     mpfr_multiplication_int_mp,&
                     mpfr_multiplication_mp_int
  END INTERFACE
  INTERFACE OPERATOR(/)
    MODULE PROCEDURE mpfr_division_mp_mp,&
                     mpfr_division_real_mp,&
                     mpfr_division_mp_real,&
                     mpfr_division_int_mp,&
                     mpfr_division_mp_int
  END INTERFACE
  INTERFACE OPERATOR(**)
    MODULE PROCEDURE mpfr_power_mp_mp,&
                     mpfr_power_mp_int,&
                     mpfr_power_int_mp,&
                     mpfr_power_real_mp,&
                     mpfr_power_mp_real
  END INTERFACE
  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE mpfr_assign_mp_real,&
                      mpfr_assign_mp_mp,&
                      mpfr_assign_mp_str,&
                      mpfr_assign_mp_int
  END INTERFACE

  INTERFACE OPERATOR(.CONVERT.)
    MODULE PROCEDURE mpfr_convert_str
  END INTERFACE

  INTERFACE OPERATOR(<)
    MODULE PROCEDURE mpfr_lt_mp_mp,&
                     mpfr_lt_mp_real,&
                     mpfr_lt_real_mp,&
                     mpfr_lt_mp_int,&
                     mpfr_lt_int_mp
  END INTERFACE

  INTERFACE OPERATOR(>)
    MODULE PROCEDURE mpfr_gt_mp_mp,&
                     mpfr_gt_mp_real,&
                     mpfr_gt_real_mp,&
                     mpfr_gt_mp_int,&
                     mpfr_gt_int_mp
  END INTERFACE
 
  INTERFACE OPERATOR(<=)
    MODULE PROCEDURE mpfr_lte_mp_mp,&
                     mpfr_lte_mp_real,&
                     mpfr_lte_real_mp,&
                     mpfr_lte_mp_int,&
                     mpfr_lte_int_mp
  END INTERFACE

  INTERFACE OPERATOR(>=)
    MODULE PROCEDURE mpfr_gte_mp_mp,&
                     mpfr_gte_mp_real,&
                     mpfr_gte_real_mp,&
                     mpfr_gte_mp_int,&
                     mpfr_gte_int_mp
  END INTERFACE

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE mpfr_eq_mp_mp,&
                     mpfr_eq_mp_real,&
                     mpfr_eq_real_mp,&
                     mpfr_eq_mp_int,&
                     mpfr_eq_int_mp
  END INTERFACE

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE mpfr_neq_mp_mp,&
                     mpfr_neq_mp_real,&
                     mpfr_neq_real_mp,&
                     mpfr_neq_mp_int,&
                     mpfr_neq_int_mp
  END INTERFACE

  INTERFACE set_value
    MODULE PROCEDURE set_value_real, set_value_int, set_value_str
  END INTERFACE

  INTERFACE log
    MODULE PROCEDURE log_mp
  END INTERFACE
  
  INTERFACE log2
    MODULE PROCEDURE log2_mp
  END INTERFACE
   
  INTERFACE log10
    MODULE PROCEDURE log10_mp
  END INTERFACE

  INTERFACE exp
    MODULE PROCEDURE exp_mp
  END INTERFACE

  INTERFACE exp2
    MODULE PROCEDURE exp2_mp
  END INTERFACE

  INTERFACE exp10
    MODULE PROCEDURE exp10_mp
  END INTERFACE

  INTERFACE cos
    MODULE PROCEDURE cos_mp
  END INTERFACE

  INTERFACE sin
    MODULE PROCEDURE sin_mp
  END INTERFACE

  INTERFACE tan
    MODULE PROCEDURE tan_mp
  END INTERFACE

  INTERFACE sec
    MODULE PROCEDURE sec_mp
  END INTERFACE

  INTERFACE csc
    MODULE PROCEDURE csc_mp
  END INTERFACE

  INTERFACE cot
    MODULE PROCEDURE cot_mp
  END INTERFACE

  INTERFACE acos
    MODULE PROCEDURE acos_mp
  END INTERFACE

  INTERFACE asin
    MODULE PROCEDURE asin_mp
  END INTERFACE

  INTERFACE atan
    MODULE PROCEDURE atan_mp
  END INTERFACE

  INTERFACE atan2
    MODULE PROCEDURE atan2_mp
  END INTERFACE

  INTERFACE cosh
    MODULE PROCEDURE cosh_mp
  END INTERFACE

  INTERFACE sinh
    MODULE PROCEDURE sinh_mp
  END INTERFACE

  INTERFACE tanh
    MODULE PROCEDURE tanh_mp
  END INTERFACE

  INTERFACE sech
    MODULE PROCEDURE sech_mp
  END INTERFACE

  INTERFACE csch
    MODULE PROCEDURE csch_mp
  END INTERFACE

  INTERFACE coth
    MODULE PROCEDURE coth_mp
  END INTERFACE

  INTERFACE acosh
    MODULE PROCEDURE acosh_mp
  END INTERFACE

  INTERFACE asinh
    MODULE PROCEDURE asinh_mp
  END INTERFACE

  INTERFACE atanh
    MODULE PROCEDURE atanh_mp
  END INTERFACE

  INTERFACE ei
    MODULE PROCEDURE ei_mp
  END INTERFACE

  INTERFACE gamma
    MODULE PROCEDURE gamma_mp
  END INTERFACE

  INTERFACE lngamma
    MODULE PROCEDURE lngamma_mp
  END INTERFACE
 
  INTERFACE erf
    MODULE PROCEDURE erf_mp
  END INTERFACE

  INTERFACE erfc
    MODULE PROCEDURE erfc_mp
  END INTERFACE

  INTERFACE bessel_j0
    MODULE PROCEDURE bessel_j0_mp
  END INTERFACE

  INTERFACE bessel_j1
    MODULE PROCEDURE bessel_j1_mp
  END INTERFACE

  INTERFACE bessel_y0
    MODULE PROCEDURE bessel_y0_mp
  END INTERFACE

  INTERFACE bessel_y1
    MODULE PROCEDURE bessel_y1_mp
  END INTERFACE

  INTERFACE sqrt
    MODULE PROCEDURE sqrt_mp
  END INTERFACE

  INTERFACE REAL
    MODULE PROCEDURE mp_to_real
  END INTERFACE

  CONTAINS

  SUBROUTINE initialize(variable,precision)
    TYPE(mpfr_type)            :: variable
    INTEGER*2, OPTIONAL        :: precision

    IF( PRESENT(precision) ) THEN
      CALL mpfr_init2(variable,precision)
    ELSE
      CALL mpfr_init(variable)
    END IF
  END SUBROUTINE initialize

  SUBROUTINE mpfr_assign_mp_real(op1,op2)
    TYPE(mpfr_type),&
        INTENT(INOUT)          :: op1
    REAL(dp),&
        INTENT(IN)             :: op2

    CALL initialize(op1)
    CALL set_value(op1,op2)
  END SUBROUTINE mpfr_assign_mp_real

  SUBROUTINE mpfr_assign_mp_int(op1,op2)
    TYPE(mpfr_type),&
        INTENT(INOUT)          :: op1
    INTEGER,&
        INTENT(IN)             :: op2

    REAL(dp)                   :: op2_real

    CALL initialize(op1)
    op2_real = REAL(op2,dp)
    CALL set_value(op1,op2_real)
  END SUBROUTINE mpfr_assign_mp_int

  SUBROUTINE mpfr_assign_mp_str(op1,op2)
    TYPE(mpfr_type),&
        INTENT(INOUT)          :: op1
    CHARACTER(LEN=*),&
        INTENT(IN)             :: op2
  
    CALL initialize(op1)
    CALL set_value(op1,op2)
  END SUBROUTINE mpfr_assign_mp_str

  SUBROUTINE mpfr_assign_mp_mp(op1,op2)
    TYPE(mpfr_type),&
        INTENT(INOUT)          :: op1
    TYPE(mpfr_type),&
        INTENT(IN)             :: op2
 
    CALL initialize(op1)
    op1%mpfr_prec = op2%mpfr_prec
    op1%mpfr_sign = op2%mpfr_sign
    op1%mpfr = op2%mpfr
    op1%mpfr_d = op2%mpfr_d
  END SUBROUTINE mpfr_assign_mp_mp

  ELEMENTAL FUNCTION mp_to_real(variable)
    REAL(dp)                    :: mp_to_real
    TYPE(mpfr_type), INTENT(IN) :: variable

    mp_to_real = mpfr_get_d(variable,GMP_RNDN)
  END FUNCTION mp_to_real

  SUBROUTINE set_value_real(variable,value)
    TYPE(mpfr_type)            :: variable
    REAL(dp)                   :: value
  
    INTEGER                    :: retval

    retval = mpfr_set_d(variable,value,GMP_RNDN)
  END SUBROUTINE set_value_real

  SUBROUTINE set_value_int(variable,value)
    TYPE(mpfr_type)            :: variable
    INTEGER                    :: value

    REAL(dp)                   :: real_value

    real_value = REAL(value,dp)
    retval = mpfr_set_d(variable,real_value,GMP_RNDN)
  END SUBROUTINE set_value_int

  SUBROUTINE set_value_str(variable,str)
    TYPE(mpfr_type)            :: variable
    CHARACTER(LEN=*)           :: str

    retval = mpfr_set_str(variable,str,10,GMP_RNDN)
  END SUBROUTINE set_value_str

  FUNCTION mpfr_convert_str(a)
    TYPE(mpfr_type)            :: mpfr_convert_str
    CHARACTER(LEN=*),&
           INTENT(IN)          :: a

    CHARACTER(LEN=MAX_CHAR)    :: buffer
    INTEGER                    :: retval

    CALL initialize(mpfr_convert_str)
    buffer = TRIM(a)//C_NULL_CHAR
    DO i=1,LEN_TRIM(a)
       IF(buffer(i:i)=="D" .OR. buffer(i:i)=="d" ) buffer(i:i)="E"
    ENDDO

    retval = mpfr_set_str(mpfr_convert_str,buffer,10,GMP_RNDN)   
  END FUNCTION mpfr_convert_str


  FUNCTION mpfr_addition_mp_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_addition_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: retval
    
    CALL initialize(mpfr_addition_mp_mp)
    retval = mpfr_add(mpfr_addition_mp_mp,a1,a2,GMP_RNDN)
  END FUNCTION mpfr_addition_mp_mp

  FUNCTION mpfr_addition_mp_real(a1,a2)
    TYPE(mpfr_type)            :: mpfr_addition_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2

    INTEGER                    :: retval
    TYPE(mpfr_type)            :: a2_mpfr    
   
    CALL initialize(a2_mpfr)
    CALL set_value(a2_mpfr,a2)
    mpfr_addition_mp_real = mpfr_addition_mp_mp(a1,a2_mpfr)
    CALL mpfr_clear(a2_mpfr)
  END FUNCTION mpfr_addition_mp_real

  FUNCTION mpfr_addition_real_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_addition_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    mpfr_addition_real_mp = mpfr_addition_mp_real(a2,a1)
  END FUNCTION mpfr_addition_real_mp
 
  FUNCTION mpfr_addition_mp_int(a1,a2)
    TYPE(mpfr_type)            :: mpfr_addition_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    INTEGER                    :: retval
    TYPE(mpfr_type)            :: a2_mpfr
    REAL(dp)                   :: a2_real

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mpfr)
    CALL set_value(a2_mpfr,a2_real)
    mpfr_addition_mp_int = mpfr_addition_mp_mp(a1,a2_mpfr)
    CALL mpfr_clear(a2_mpfr)
  END FUNCTION mpfr_addition_mp_int

  FUNCTION mpfr_addition_int_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_addition_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    mpfr_addition_int_mp = mpfr_addition_mp_int(a2,a1)
  END FUNCTION mpfr_addition_int_mp

  FUNCTION mpfr_subtraction_mp_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_subtraction_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: retval
    
    CALL initialize(mpfr_subtraction_mp_mp)
    retval = mpfr_sub(mpfr_subtraction_mp_mp,a1,a2,GMP_RNDN)
  END FUNCTION mpfr_subtraction_mp_mp

  FUNCTION mpfr_minus(a1)
    TYPE(mpfr_type)            :: mpfr_minus
    TYPE(mpfr_type),INTENT(IN) :: a1

    INTEGER                    :: retval
    
    CALL initialize(mpfr_minus)
    mpfr_minus = 0.0_dp - a1
  END FUNCTION mpfr_minus

  FUNCTION mpfr_subtraction_real_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_subtraction_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: retval
    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    mpfr_subtraction_real_mp = mpfr_subtraction_mp_mp(a1_mp,a2)
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_subtraction_real_mp

  FUNCTION mpfr_subtraction_mp_real(a1,a2)
    TYPE(mpfr_type)            :: mpfr_subtraction_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2
    
    TYPE(mpfr_type)            :: a2_mp

    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2)
    mpfr_subtraction_mp_real = mpfr_subtraction_mp_mp(a1,a2_mp)
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_subtraction_mp_real

  FUNCTION mpfr_subtraction_int_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_subtraction_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    REAL(dp)                   :: a1_real
    TYPE(mpfr_type)            :: a1_mp

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    mpfr_subtraction_int_mp = mpfr_subtraction_mp_mp(a1_mp,a2)
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_subtraction_int_mp

  FUNCTION mpfr_subtraction_mp_int(a1,a2)
    TYPE(mpfr_type)            :: mpfr_subtraction_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    REAL(dp)                   :: a2_real
    TYPE(mpfr_type)            :: a2_mp

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    mpfr_subtraction_mp_int = mpfr_subtraction_mp_mp(a1,a2_mp)
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_subtraction_mp_int

  FUNCTION mpfr_multiplication_mp_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_multiplication_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: retval
    
    CALL initialize(mpfr_multiplication_mp_mp)
    retval = mpfr_mul(mpfr_multiplication_mp_mp,a1,a2,GMP_RNDN)
  END FUNCTION mpfr_multiplication_mp_mp

  FUNCTION mpfr_multiplication_real_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_multiplication_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: retval
    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    mpfr_multiplication_real_mp = mpfr_multiplication_mp_mp(a1_mp,a2)
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_multiplication_real_mp

  FUNCTION mpfr_multiplication_mp_real(a1,a2)
    TYPE(mpfr_type)            :: mpfr_multiplication_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2
    
    mpfr_multiplication_mp_real = mpfr_multiplication_real_mp(a2,a1)
  END FUNCTION mpfr_multiplication_mp_real

  FUNCTION mpfr_multiplication_int_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_multiplication_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: retval
    TYPE(mpfr_type)            :: a1_mp
    REAL(dp)                   :: a1_real

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    mpfr_multiplication_int_mp = mpfr_multiplication_mp_mp(a1_mp,a2)
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_multiplication_int_mp

  FUNCTION mpfr_multiplication_mp_int(a1,a2)
    TYPE(mpfr_type)            :: mpfr_multiplication_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    INTEGER                    :: retval

    mpfr_multiplication_mp_int = mpfr_multiplication_int_mp(a2,a1)
  END FUNCTION mpfr_multiplication_mp_int

  FUNCTION mpfr_division_mp_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_division_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: retval
    
    CALL initialize(mpfr_division_mp_mp)
    retval = mpfr_div(mpfr_division_mp_mp,a1,a2,GMP_RNDN)
  END FUNCTION mpfr_division_mp_mp

  FUNCTION mpfr_division_real_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_division_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    mpfr_division_real_mp = mpfr_division_mp_mp(a1_mp,a2)
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_division_real_mp

  FUNCTION mpfr_division_mp_real(a1,a2)
    TYPE(mpfr_type)            :: mpfr_division_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2
   
    TYPE(mpfr_type)            :: a2_mp
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2)
    mpfr_division_mp_real = mpfr_division_mp_mp(a1,a2_mp)
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_division_mp_real

  FUNCTION mpfr_division_int_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_division_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    TYPE(mpfr_type)            :: a1_mp
    REAL(dp)                   :: a1_real

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    mpfr_division_int_mp = mpfr_division_mp_mp(a1_mp,a2)
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_division_int_mp

  FUNCTION mpfr_division_mp_int(a1,a2)
    TYPE(mpfr_type)            :: mpfr_division_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    TYPE(mpfr_type)            :: a2_mp
    REAL(dp)                   :: a2_real

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    mpfr_division_mp_int = mpfr_division_mp_mp(a1,a2_mp)
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_division_mp_int

  FUNCTION mpfr_power_mp_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_power_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: retval
    
    CALL initialize(mpfr_power_mp_mp)
    retval = mpfr_pow(mpfr_power_mp_mp,a1,a2,GMP_RNDN)
  END FUNCTION mpfr_power_mp_mp

  FUNCTION mpfr_power_real_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_power_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    mpfr_power_real_mp = mpfr_power_mp_mp(a1_mp,a2)
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_power_real_mp

  FUNCTION mpfr_power_mp_real(a1,a2)
    TYPE(mpfr_type)            :: mpfr_power_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2

    TYPE(mpfr_type)            :: a2_mp

    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2) 
    mpfr_power_mp_real = mpfr_power_mp_mp(a1,a2_mp)
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_power_mp_real

  FUNCTION mpfr_power_mp_int(a1,a2)
    TYPE(mpfr_type)            :: mpfr_power_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    REAL(dp)                   :: a2_real
    TYPE(mpfr_type)            :: a2_mp

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    mpfr_power_mp_int = mpfr_power_mp_mp(a1,a2_mp)
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_power_mp_int

  FUNCTION mpfr_power_int_mp(a1,a2)
    TYPE(mpfr_type)            :: mpfr_power_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2
 
    REAL(dp)                   :: a1_real
    TYPE(mpfr_type)            :: a1_mp

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    mpfr_power_int_mp = mpfr_power_mp_mp(a1_mp,a2)
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_power_int_mp

  FUNCTION mpfr_lt_mp_mp(a1,a2)
    LOGICAL                    :: mpfr_lt_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: comp

    comp = mpfr_cmp(a1,a2)
    IF( comp < 0) THEN 
      mpfr_lt_mp_mp = .TRUE.
    ELSE IF(comp>=0) THEN
      mpfr_lt_mp_mp = .FALSE.    
    END IF
  END FUNCTION mpfr_lt_mp_mp

  FUNCTION mpfr_lt_mp_real(a1,a2)
    LOGICAL                    :: mpfr_lt_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp

    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp < 0) THEN 
      mpfr_lt_mp_real = .TRUE.
    ELSE IF(comp>=0) THEN
      mpfr_lt_mp_real = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_lt_mp_real

  FUNCTION mpfr_lt_real_mp(a1,a2)
    LOGICAL                    :: mpfr_lt_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp < 0) THEN 
      mpfr_lt_real_mp = .TRUE.
    ELSE IF(comp>=0) THEN
      mpfr_lt_real_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_lt_real_mp

  FUNCTION mpfr_lt_mp_int(a1,a2)
    LOGICAL                    :: mpfr_lt_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp
    REAL(dp)                   :: a2_real

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp < 0) THEN 
      mpfr_lt_mp_int = .TRUE.
    ELSE IF(comp>=0) THEN
      mpfr_lt_mp_int = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_lt_mp_int

  FUNCTION mpfr_lt_int_mp(a1,a2)
    LOGICAL                    :: mpfr_lt_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp
    REAL(dp)                   :: a1_real

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp < 0) THEN 
      mpfr_lt_int_mp = .TRUE.
    ELSE IF(comp>=0) THEN
      mpfr_lt_int_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_lt_int_mp

  FUNCTION mpfr_gt_mp_mp(a1,a2)
    LOGICAL                    :: mpfr_gt_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: comp

    comp = mpfr_cmp(a1,a2)
    IF( comp > 0) THEN 
      mpfr_gt_mp_mp = .TRUE.
    ELSE IF(comp<=0) THEN
      mpfr_gt_mp_mp = .FALSE.    
    END IF
  END FUNCTION mpfr_gt_mp_mp

  FUNCTION mpfr_gt_mp_real(a1,a2)
    LOGICAL                    :: mpfr_gt_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp

    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp > 0) THEN 
      mpfr_gt_mp_real = .TRUE.
    ELSE IF(comp<=0) THEN
      mpfr_gt_mp_real = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_gt_mp_real

  FUNCTION mpfr_gt_real_mp(a1,a2)
    LOGICAL                    :: mpfr_gt_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp > 0) THEN 
      mpfr_gt_real_mp = .TRUE.
    ELSE IF(comp<=0) THEN
      mpfr_gt_real_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_gt_real_mp

  FUNCTION mpfr_gt_mp_int(a1,a2)
    LOGICAL                    :: mpfr_gt_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp
    REAL(dp)                   :: a2_real

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp > 0) THEN 
      mpfr_gt_mp_int = .TRUE.
    ELSE IF(comp<=0) THEN
      mpfr_gt_mp_int = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_gt_mp_int

  FUNCTION mpfr_gt_int_mp(a1,a2)
    LOGICAL                    :: mpfr_gt_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp
    REAL(dp)                   :: a1_real

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp > 0) THEN 
      mpfr_gt_int_mp = .TRUE.
    ELSE IF(comp<=0) THEN
      mpfr_gt_int_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_gt_int_mp

  FUNCTION mpfr_lte_mp_mp(a1,a2)
    LOGICAL                    :: mpfr_lte_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: comp

    comp = mpfr_cmp(a1,a2)
    IF( comp <= 0) THEN 
      mpfr_lte_mp_mp = .TRUE.
    ELSE IF(comp>0) THEN
      mpfr_lte_mp_mp = .FALSE.    
    END IF
  END FUNCTION mpfr_lte_mp_mp

  FUNCTION mpfr_lte_mp_real(a1,a2)
    LOGICAL                    :: mpfr_lte_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp

    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp <= 0) THEN 
      mpfr_lte_mp_real = .TRUE.
    ELSE IF(comp>0) THEN
      mpfr_lte_mp_real = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_lte_mp_real

  FUNCTION mpfr_lte_real_mp(a1,a2)
    LOGICAL                    :: mpfr_lte_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp <= 0) THEN 
      mpfr_lte_real_mp = .TRUE.
    ELSE IF(comp>0) THEN
      mpfr_lte_real_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_lte_real_mp

  FUNCTION mpfr_lte_mp_int(a1,a2)
    LOGICAL                    :: mpfr_lte_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp
    REAL(dp)                   :: a2_real

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp <= 0) THEN 
      mpfr_lte_mp_int = .TRUE.
    ELSE IF(comp>0) THEN
      mpfr_lte_mp_int = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_lte_mp_int

  FUNCTION mpfr_lte_int_mp(a1,a2)
    LOGICAL                    :: mpfr_lte_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp
    REAL(dp)                   :: a1_real

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp <= 0) THEN 
      mpfr_lte_int_mp = .TRUE.
    ELSE IF(comp>0) THEN
      mpfr_lte_int_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_lte_int_mp

  FUNCTION mpfr_gte_mp_mp(a1,a2)
    LOGICAL                    :: mpfr_gte_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: comp

    comp = mpfr_cmp(a1,a2)
    IF( comp >= 0) THEN 
      mpfr_gte_mp_mp = .TRUE.
    ELSE IF(comp<0) THEN
      mpfr_gte_mp_mp = .FALSE.    
    END IF
  END FUNCTION mpfr_gte_mp_mp

  FUNCTION mpfr_gte_mp_real(a1,a2)
    LOGICAL                    :: mpfr_gte_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp

    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp >= 0) THEN 
      mpfr_gte_mp_real = .TRUE.
    ELSE IF(comp>0) THEN
      mpfr_gte_mp_real = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_gte_mp_real

  FUNCTION mpfr_gte_real_mp(a1,a2)
    LOGICAL                    :: mpfr_gte_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp >= 0) THEN 
      mpfr_gte_real_mp = .TRUE.
    ELSE IF(comp<0) THEN
      mpfr_gte_real_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_gte_real_mp

  FUNCTION mpfr_gte_mp_int(a1,a2)
    LOGICAL                    :: mpfr_gte_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp
    REAL(dp)                   :: a2_real

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp >= 0) THEN 
      mpfr_gte_mp_int = .TRUE.
    ELSE IF(comp<0) THEN
      mpfr_gte_mp_int = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_gte_mp_int

  FUNCTION mpfr_gte_int_mp(a1,a2)
    LOGICAL                    :: mpfr_gte_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp
    REAL(dp)                   :: a1_real

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp >= 0) THEN 
      mpfr_gte_int_mp = .TRUE.
    ELSE IF(comp<0) THEN
      mpfr_gte_int_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_gte_int_mp

  FUNCTION mpfr_eq_mp_mp(a1,a2)
    LOGICAL                    :: mpfr_eq_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: comp

    comp = mpfr_cmp(a1,a2)
    IF( comp == 0) THEN 
      mpfr_eq_mp_mp = .TRUE.
    ELSE
      mpfr_eq_mp_mp = .FALSE.    
    END IF
  END FUNCTION mpfr_eq_mp_mp

  FUNCTION mpfr_eq_mp_real(a1,a2)
    LOGICAL                    :: mpfr_eq_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp

    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp == 0) THEN 
      mpfr_eq_mp_real = .TRUE.
    ELSE
      mpfr_eq_mp_real = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_eq_mp_real

  FUNCTION mpfr_eq_real_mp(a1,a2)
    LOGICAL                    :: mpfr_eq_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp == 0) THEN 
      mpfr_eq_real_mp = .TRUE.
    ELSE 
      mpfr_eq_real_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_eq_real_mp

  FUNCTION mpfr_eq_mp_int(a1,a2)
    LOGICAL                    :: mpfr_eq_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp
    REAL(dp)                   :: a2_real

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp == 0) THEN 
      mpfr_eq_mp_int = .TRUE.
    ELSE
      mpfr_eq_mp_int = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_eq_mp_int

  FUNCTION mpfr_eq_int_mp(a1,a2)
    LOGICAL                    :: mpfr_eq_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp
    REAL(dp)                   :: a1_real

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp == 0) THEN 
      mpfr_eq_int_mp = .TRUE.
    ELSE
      mpfr_eq_int_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_eq_int_mp

  FUNCTION mpfr_neq_mp_mp(a1,a2)
    LOGICAL                    :: mpfr_neq_mp_mp
    TYPE(mpfr_type),INTENT(IN) :: a1,a2

    INTEGER                    :: comp

    comp = mpfr_cmp(a1,a2)
    IF( comp /= 0) THEN 
      mpfr_neq_mp_mp = .TRUE.
    ELSE
      mpfr_neq_mp_mp = .FALSE.    
    END IF
  END FUNCTION mpfr_neq_mp_mp

  FUNCTION mpfr_neq_mp_real(a1,a2)
    LOGICAL                    :: mpfr_neq_mp_real
    TYPE(mpfr_type),INTENT(IN) :: a1
    REAL(dp), INTENT(IN)       :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp

    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp /= 0) THEN 
      mpfr_neq_mp_real = .TRUE.
    ELSE
      mpfr_neq_mp_real = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_neq_mp_real

  FUNCTION mpfr_neq_real_mp(a1,a2)
    LOGICAL                    :: mpfr_neq_real_mp
    REAL(dp), INTENT(IN)       :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp

    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp /= 0) THEN 
      mpfr_neq_real_mp = .TRUE.
    ELSE 
      mpfr_neq_real_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_neq_real_mp

  FUNCTION mpfr_neq_mp_int(a1,a2)
    LOGICAL                    :: mpfr_neq_mp_int
    TYPE(mpfr_type),INTENT(IN) :: a1
    INTEGER, INTENT(IN)        :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a2_mp
    REAL(dp)                   :: a2_real

    a2_real = REAL(a2,dp)
    CALL initialize(a2_mp)
    CALL set_value(a2_mp,a2_real)
    comp = mpfr_cmp(a1,a2_mp)
    IF( comp /= 0) THEN 
      mpfr_neq_mp_int = .TRUE.
    ELSE
      mpfr_neq_mp_int = .FALSE.    
    END IF
    CALL mpfr_clear(a2_mp)
  END FUNCTION mpfr_neq_mp_int

  FUNCTION mpfr_neq_int_mp(a1,a2)
    LOGICAL                    :: mpfr_neq_int_mp
    INTEGER, INTENT(IN)        :: a1
    TYPE(mpfr_type),INTENT(IN) :: a2

    INTEGER                    :: comp
    TYPE(mpfr_type)            :: a1_mp
    REAL(dp)                   :: a1_real

    a1_real = REAL(a1,dp)
    CALL initialize(a1_mp)
    CALL set_value(a1_mp,a1_real)
    comp = mpfr_cmp(a1_mp,a2)
    IF( comp /= 0) THEN 
      mpfr_neq_int_mp = .TRUE.
    ELSE
      mpfr_neq_int_mp = .FALSE.    
    END IF
    CALL mpfr_clear(a1_mp)
  END FUNCTION mpfr_neq_int_mp

  FUNCTION log_mp(op)
    TYPE(mpfr_type)            :: log_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(log_mp)
    retval = mpfr_log(log_mp,op,GMP_RNDN)
  END FUNCTION log_mp    

  FUNCTION log2_mp(op)
    TYPE(mpfr_type)            :: log2_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(log2_mp)
    retval = mpfr_log2(log2_mp,op,GMP_RNDN)
  END FUNCTION log2_mp    

  FUNCTION log10_mp(op)
    TYPE(mpfr_type)            :: log10_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(log10_mp)
    retval = mpfr_log10(log10_mp,op,GMP_RNDN)
  END FUNCTION log10_mp    

  FUNCTION exp_mp(op)
    TYPE(mpfr_type)            :: exp_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(exp_mp)
    retval = mpfr_exp(exp_mp,op,GMP_RNDN)
  END FUNCTION exp_mp    

  FUNCTION exp2_mp(op)
    TYPE(mpfr_type)            :: exp2_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(exp2_mp)
    retval = mpfr_exp2(exp2_mp,op,GMP_RNDN)
  END FUNCTION exp2_mp    

  FUNCTION exp10_mp(op)
    TYPE(mpfr_type)            :: exp10_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(exp10_mp)
    retval = mpfr_exp10(exp10_mp,op,GMP_RNDN)
  END FUNCTION exp10_mp    

  FUNCTION sin_mp(op)
    TYPE(mpfr_type)            :: sin_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(sin_mp)
    retval = mpfr_sin(sin_mp,op,GMP_RNDN)
  END FUNCTION sin_mp    

  FUNCTION cos_mp(op)
    TYPE(mpfr_type)            :: cos_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(cos_mp)
    retval = mpfr_cos(cos_mp,op,GMP_RNDN)
  END FUNCTION cos_mp    

  FUNCTION tan_mp(op)
    TYPE(mpfr_type)            :: tan_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(tan_mp)
    retval = mpfr_tan(tan_mp,op,GMP_RNDN)
  END FUNCTION tan_mp    

  FUNCTION sec_mp(op)
    TYPE(mpfr_type)            :: sec_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(sec_mp)
    retval = mpfr_sec(sec_mp,op,GMP_RNDN)
  END FUNCTION sec_mp    

  FUNCTION csc_mp(op)
    TYPE(mpfr_type)            :: csc_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(csc_mp)
    retval = mpfr_csc(csc_mp,op,GMP_RNDN)
  END FUNCTION csc_mp    

  FUNCTION cot_mp(op)
    TYPE(mpfr_type)            :: cot_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(cot_mp)
    retval = mpfr_cot(cot_mp,op,GMP_RNDN)
  END FUNCTION cot_mp    

  FUNCTION acos_mp(op)
    TYPE(mpfr_type)            :: acos_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(acos_mp)
    retval = mpfr_acos(acos_mp,op,GMP_RNDN)
  END FUNCTION acos_mp    

  FUNCTION asin_mp(op)
    TYPE(mpfr_type)            :: asin_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(asin_mp)
    retval = mpfr_asin(asin_mp,op,GMP_RNDN)
  END FUNCTION asin_mp    

  FUNCTION atan_mp(op)
    TYPE(mpfr_type)            :: atan_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(atan_mp)
    retval = mpfr_atan(atan_mp,op,GMP_RNDN)
  END FUNCTION atan_mp    

  FUNCTION atan2_mp(x,y)
    TYPE(mpfr_type)            :: atan2_mp
    TYPE(mpfr_type)            :: x,y

    INTEGER                    :: retval

    CALL initialize(atan2_mp)
    retval = mpfr_atan2(atan2_mp,x,y,GMP_RNDN)
  END FUNCTION atan2_mp    

  FUNCTION cosh_mp(op)
    TYPE(mpfr_type)            :: cosh_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(cosh_mp)
    retval = mpfr_cosh(cosh_mp,op,GMP_RNDN)
  END FUNCTION cosh_mp    

  FUNCTION sinh_mp(op)
    TYPE(mpfr_type)            :: sinh_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(sinh_mp)
    retval = mpfr_sinh(sinh_mp,op,GMP_RNDN)
  END FUNCTION sinh_mp    

  FUNCTION tanh_mp(op)
    TYPE(mpfr_type)            :: tanh_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(tanh_mp)
    retval = mpfr_tanh(tanh_mp,op,GMP_RNDN)
  END FUNCTION tanh_mp    

  FUNCTION sech_mp(op)
    TYPE(mpfr_type)            :: sech_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(sech_mp)
    retval = mpfr_sech(sech_mp,op,GMP_RNDN)
  END FUNCTION sech_mp    

  FUNCTION csch_mp(op)
    TYPE(mpfr_type)            :: csch_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(csch_mp)
    retval = mpfr_csch(csch_mp,op,GMP_RNDN)
  END FUNCTION csch_mp    

  FUNCTION coth_mp(op)
    TYPE(mpfr_type)            :: coth_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(coth_mp)
    retval = mpfr_coth(coth_mp,op,GMP_RNDN)
  END FUNCTION coth_mp    

  FUNCTION acosh_mp(op)
    TYPE(mpfr_type)            :: acosh_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(acosh_mp)
    retval = mpfr_acosh(acosh_mp,op,GMP_RNDN)
  END FUNCTION acosh_mp    

  FUNCTION asinh_mp(op)
    TYPE(mpfr_type)            :: asinh_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(asinh_mp)
    retval = mpfr_asinh(asinh_mp,op,GMP_RNDN)
  END FUNCTION asinh_mp    

  FUNCTION atanh_mp(op)
    TYPE(mpfr_type)            :: atanh_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(atanh_mp)
    retval = mpfr_atanh(atanh_mp,op,GMP_RNDN)
  END FUNCTION atanh_mp    

  FUNCTION ei_mp(op)
    TYPE(mpfr_type)            :: ei_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(ei_mp)
    retval = mpfr_eint(ei_mp,op,GMP_RNDN)
  END FUNCTION ei_mp

  FUNCTION gamma_mp(op)
    TYPE(mpfr_type)            :: gamma_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(gamma_mp)
    retval = mpfr_gamma(gamma_mp,op,GMP_RNDN)
  END FUNCTION gamma_mp    

  FUNCTION lngamma_mp(op)
    TYPE(mpfr_type)            :: lngamma_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(lngamma_mp)
    retval = mpfr_lngamma(lngamma_mp,op,GMP_RNDN)
  END FUNCTION lngamma_mp    

  FUNCTION erf_mp(op)
    TYPE(mpfr_type)            :: erf_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(erf_mp)
    retval = mpfr_erf(erf_mp,op,GMP_RNDN)
  END FUNCTION erf_mp    

  FUNCTION erfc_mp(op)
    TYPE(mpfr_type)            :: erfc_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(erfc_mp)
    retval = mpfr_erfc(erfc_mp,op,GMP_RNDN)
  END FUNCTION erfc_mp    

  FUNCTION bessel_j0_mp(op)
    TYPE(mpfr_type)            :: bessel_j0_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(bessel_j0_mp)
    retval = mpfr_bessel_j0(bessel_j0_mp,op,GMP_RNDN)
  END FUNCTION bessel_j0_mp    

  FUNCTION bessel_j1_mp(op)
    TYPE(mpfr_type)            :: bessel_j1_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(bessel_j1_mp)
    retval = mpfr_bessel_j1(bessel_j1_mp,op,GMP_RNDN)
  END FUNCTION bessel_j1_mp    

  FUNCTION bessel_y0_mp(op)
    TYPE(mpfr_type)            :: bessel_y0_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(bessel_y0_mp)
    retval = mpfr_bessel_y0(bessel_y0_mp,op,GMP_RNDN)
  END FUNCTION bessel_y0_mp    

  FUNCTION bessel_y1_mp(op)
    TYPE(mpfr_type)            :: bessel_y1_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(bessel_y1_mp)
    retval = mpfr_bessel_y1(bessel_y1_mp,op,GMP_RNDN)
  END FUNCTION bessel_y1_mp    

  FUNCTION get_pi()
    TYPE(mpfr_type)            :: get_pi

    INTEGER                    :: retval

    CALL initialize(get_pi)    
    retval = mpfr_const_pi(get_pi,GMP_RNDN)
  END FUNCTION get_pi
  
  FUNCTION get_e()
    TYPE(mpfr_type)            :: get_e

    INTEGER                    :: retval

    CALL initialize(get_e)    
    retval = mpfr_const_euler(get_e,GMP_RNDN)
  END FUNCTION get_e

  FUNCTION get_log2()
    TYPE(mpfr_type)            :: get_log2

    INTEGER                    :: retval

    CALL initialize(get_log2)    
    retval = mpfr_const_log2(get_log2,GMP_RNDN)
  END FUNCTION get_log2

  FUNCTION sqrt_mp(op)
    TYPE(mpfr_type)            :: sqrt_mp
    TYPE(mpfr_type)            :: op

    INTEGER                    :: retval

    CALL initialize(sqrt_mp)
    retval = mpfr_sqrt(sqrt_mp,op,GMP_RNDN)
  END FUNCTION sqrt_mp    

END MODULE mpfr_ops
