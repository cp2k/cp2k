!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Writes information on XC functionals to output
! *****************************************************************************
MODULE xc_write_output
  
  USE input_constants,                 ONLY: xc_none
  USE input_cp2k_check,                ONLY: xc_functionals_expand
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_get_subs_vals2,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length
  USE xc_derivatives,                  ONLY: xc_functional_get_info
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_write_output'

  PUBLIC :: xc_write

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param output_unit ...
!> \param xc_section ...
!> \param lsd ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE xc_write(output_unit, xc_section, lsd, error)
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(section_vals_type), POINTER         :: xc_section
    LOGICAL, INTENT(IN)                      :: lsd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'xc_write', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2*default_string_length)   :: reference, shortform
    CHARACTER(LEN=20)                        :: tmpStr
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: func_name
    INTEGER                                  :: ifun, ifunc_name, il, myfun, &
                                                vmajor, vminor, vmicro
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: xc_fun, xc_fun_section

    failure = .FALSE.
    IF (output_unit>0) THEN

       xc_fun_section => section_vals_get_subs_vals(xc_section,&
            "XC_FUNCTIONAL",error=error)
       CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",i_val=myfun,error=error)
       IF (myfun/=xc_none) THEN
          CALL section_vals_val_get(xc_section,"FUNCTIONAL_ROUTINE",&
               c_val=tmpStr,error=error)
          WRITE (output_unit, fmt="(' FUNCTIONAL| ROUTINE=',a)")&
               TRIM(tmpStr)
          CALL xc_functionals_expand(xc_fun_section,xc_section,error=error)
          ifun=0
          DO
             ifun=ifun+1
             xc_fun => section_vals_get_subs_vals2(xc_fun_section,i_section=ifun,&
                  error=error)
             IF (.NOT.ASSOCIATED(xc_fun)) EXIT
             CALL xc_functional_get_info(xc_fun,&
                                         lsd=lsd,&
                                         reference=reference,&
                                         shortform=shortform,&
                                         ifunc_name=1,&
                                         error=error)
             IF (TRIM(xc_fun%section%name) == "LIBXC") THEN
#if defined(__LIBXC3)
               CALL xc_f90_version(vmajor,vminor,vmicro)
#elif defined(__LIBXC2)
               CALL xc_f90_version(vmajor,vminor)
               vmicro=0
#else
               vmajor=0 ; vminor=0 ; vmicro=0
#endif
               WRITE (output_unit, fmt="(A,i1,A,i1,A,i1,A)") ' FUNCTIONAL| LIBXC Vers. ', vmajor, ".", vminor, &
                      ".", vmicro,' (Marques, Oliveira, Burnus, CPC 183, 2272 (2012))'
               WRITE (output_unit, fmt="(' FUNCTIONAL| ',a,':')") TRIM(shortform)
               il = LEN_TRIM(reference)
               IF (il < 67) THEN
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:il)
               ELSE IF (il < 134) THEN
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(68:il)
               ELSE
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(68:134)
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(134:il)
               END IF
               CALL section_vals_val_get(xc_fun,"functional",c_vals=func_name,error=error)
               DO ifunc_name=2, SIZE(func_name)
                 CALL xc_functional_get_info(xc_fun,&
                                             lsd=lsd,&
                                             reference=reference,&
                                             shortform=shortform,&
                                             ifunc_name=ifunc_name,&
                                             error=error)
               WRITE (output_unit, fmt="(A,i1,A,i1,A)") ' FUNCTIONAL| LIBXC version ', vmajor, ".", vminor, &
                      ' (Marques, Oliveira, Burnus, CPC 183, 2272 (2012))'
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a,':')") TRIM(shortform)
                 il = LEN_TRIM(reference)
                 IF (il < 67) THEN
                   WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:il)
                 ELSE IF (il < 134) THEN
                   WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
                   WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(68:il)
                 ELSE
                   WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
                   WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(68:134)
                   WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(134:il)
                 END IF
               END DO
             ELSE
               WRITE (output_unit, fmt="(' FUNCTIONAL| ',a,':')")&
                    TRIM(xc_fun%section%name)
               il = LEN_TRIM(reference)
               IF (il < 67) THEN
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:il)
               ELSE IF (il < 134) THEN
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(68:il)
               ELSE
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(1:67)
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(68:134)
                 WRITE (output_unit, fmt="(' FUNCTIONAL| ',a)") reference(134:il)
               END IF
             END IF
          END DO
       ELSE
          WRITE (output_unit, fmt="(' FUNCTIONAL| NO EXCHANGE-CORRELATION FUNCTIONAL USED.')")
       END IF
    END IF

  END SUBROUTINE xc_write

END MODULE xc_write_output
