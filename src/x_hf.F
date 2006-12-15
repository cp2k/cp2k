!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!***** cp2k/x_hf [1.0] *
!!
!!   NAME
!!     x_hf
!!
!!   FUNCTION
!!     routines and types for Hartree-Fock-Exchange
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Manuel Guidon 
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!****************************************************************************
MODULE x_hf   
  USE kinds,                           ONLY: dp,dp_size
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_get,&
                                             section_vals_val_get


#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  x_hf_type, x_hf_create, x_hf_release 
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'x_hf'

!***
!****************************************************************************

!!***** x_hf/x_hf_type [1.0] *
!!
!!   NAME
!!     x_hf_type
!!
!!   FUNCTION
!!     stores some data used in construction of Kohn-Sham matrix
!!     holds a certain amount of precomputed four-center-integrals 
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - x_buffer: stores four-center-integrals
!!     - number_of_elements: according to memory usage parameter in input file
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!***************************************************************************

  TYPE x_hf_type
    REAL(dp), DIMENSION(:), POINTER          :: x_buffer
    INTEGER                                  :: number_of_elements
  END TYPE

  CONTAINS

  SUBROUTINE x_hf_create(x_data,hf_section,error)
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(section_vals_type), POINTER         :: hf_section
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    INTEGER                                  :: number_of_precalc_integrals,&
                                                max_memory, stat
    LOGICAL                                  :: failure
    CHARACTER(LEN=*), PARAMETER              :: routineN = 'x_hf_create',&
                                                routineP = moduleN//':'//routineN

    CALL section_vals_val_get(hf_section,"MAX_MEMORY",i_val=max_memory,error=error)
    number_of_precalc_integrals = max_memory * 1024*1024/dp_size
    ALLOCATE(x_data,STAT=stat)
    ALLOCATE(x_data%x_buffer(number_of_precalc_integrals),STAT=stat)
    x_data%number_of_elements = number_of_precalc_integrals

  END SUBROUTINE x_hf_create 

  SUBROUTINE x_hf_release(x_data,error)
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)  &
                                             :: error

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    CHARACTER(LEN=*), PARAMETER              :: routineN = 'x_hf_release',&
                                                routineP = moduleN//':'//routineN

    DEALLOCATE(x_data%x_buffer,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(x_data,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE x_hf_release

END MODULE x_hf


