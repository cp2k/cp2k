!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             dp_size
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  x_hf_type, x_hf_create, x_hf_release, x_hf_screening_type ,&
          x_hf_set_distr_energy, proc_assoc_energy,&
          x_hf_set_distr_forces, proc_assoc_forces

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
!!     - x_hf_screening_type: stores information for screening
!!
!!   AUTHOR
!!     Manuel Guidon
!!
!!   MODIFICATION HISTORY
!!     11.2006 created [Manuel Guidon]
!!
!***************************************************************************
  TYPE x_hf_screening_type
    REAL(dp)                                 :: omega
    INTEGER                                  :: screening_type
  END TYPE
  
  TYPE x_hf_type
    REAL(dp), DIMENSION(:), POINTER          :: x_buffer
    REAL(dp), DIMENSION(:,:,:,:), &
              POINTER                        :: x_screening_buffer
    INTEGER                                  :: number_of_elements
    TYPE(x_hf_screening_type)                :: screening_parameter
    INTEGER*8, DIMENSION(:), POINTER         :: integral_distribution_energy
    INTEGER*8, DIMENSION(:), POINTER         :: integral_distribution_forces
  END TYPE

  CONTAINS

  SUBROUTINE x_hf_create(x_data,hf_section,natom,max_set,error)
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(section_vals_type), POINTER         :: hf_section
    INTEGER, INTENT(IN)                      :: natom, max_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: max_memory, &
                                                number_of_precalc_integrals, &
                                                stat, TYPE
    REAL(dp)                                 :: omega

    CALL section_vals_val_get(hf_section,"MAX_MEMORY",i_val=max_memory,error=error)
    number_of_precalc_integrals = max_memory * 1024*1024/dp_size
    ALLOCATE(x_data,STAT=stat)
    ALLOCATE(x_data%x_buffer(number_of_precalc_integrals),STAT=stat)
    x_data%number_of_elements = number_of_precalc_integrals
    CALL section_vals_val_get(hf_section,"screening_type",i_val=TYPE,error=error)
    CALL section_vals_val_get(hf_section,"omega",r_val=omega,error=error)
    ALLOCATE(x_data%x_screening_buffer(natom,natom,max_set,max_set))
 
    x_data%screening_parameter%screening_type = TYPE
    x_data%screening_parameter%omega = omega
#ifdef __parallel
    ALLOCATE(x_data%integral_distribution_energy(1))
    ALLOCATE(x_data%integral_distribution_forces(1))
#endif
  END SUBROUTINE x_hf_create 

  SUBROUTINE x_hf_release(x_data,error)
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    DEALLOCATE(x_data%x_buffer,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(x_data%x_screening_buffer,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
#ifdef __parallel
    DEALLOCATE(x_data%integral_distribution_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(x_data%integral_distribution_forces,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
#endif
    DEALLOCATE(x_data,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE x_hf_release

  SUBROUTINE x_hf_set_distr_energy(distr, x_data, error)
    INTEGER*8, DIMENSION(:), POINTER         :: distr
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error
    
    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'x_hf_set_dist_energy', &
      routineP = moduleN//':'//routineN

    DEALLOCATE(x_data%integral_distribution_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%integral_distribution_energy(SIZE(distr)))
    x_data%integral_distribution_energy = distr
  END SUBROUTINE x_hf_set_distr_energy

  FUNCTION proc_assoc_energy(counter, distribution)
    LOGICAL                                  :: proc_assoc_energy
    INTEGER*8                                :: counter
    INTEGER*8, DIMENSION(:), POINTER         :: distribution

    INTEGER                                  :: mysize, ipos

#ifndef __parallel
    proc_assoc_energy = .TRUE.
    RETURN
#else
    mysize = SIZE(distribution,1)/2
    
    proc_assoc_energy = .FALSE.
    DO ipos = 1,mysize
      IF(counter>=distribution(2*ipos-1) .AND. counter<=distribution(2*ipos) ) THEN
        proc_assoc_energy = .TRUE.
        RETURN
      END IF
    END DO
#endif
  END FUNCTION proc_assoc_energy

  SUBROUTINE x_hf_set_distr_forces(distr, x_data, error)
    INTEGER*8, DIMENSION(:), POINTER         :: distr
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'x_hf_set_distr_forces', &
      routineP = moduleN//':'//routineN


    DEALLOCATE(x_data%integral_distribution_forces,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%integral_distribution_forces(SIZE(distr)))
    x_data%integral_distribution_forces = distr
  END SUBROUTINE x_hf_set_distr_forces

  FUNCTION proc_assoc_forces(counter, distribution)
    LOGICAL                                  :: proc_assoc_forces
    INTEGER*8                                :: counter
    INTEGER*8, DIMENSION(:), POINTER         :: distribution

    INTEGER                                  :: mysize, ipos

#ifndef __parallel
    proc_assoc_forces = .TRUE.
    RETURN
#else
    mysize = SIZE(distribution,1)/2
    
    proc_assoc_forces = .FALSE.
    DO ipos = 1,mysize
      IF(counter>=distribution(2*ipos-1) .AND. counter<=distribution(2*ipos) ) THEN
        proc_assoc_forces = .TRUE.
        RETURN
      END IF
    END DO
#endif
  END FUNCTION proc_assoc_forces
 
END MODULE x_hf


