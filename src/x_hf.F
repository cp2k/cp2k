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
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_8
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  x_hf_type, x_hf_create, x_hf_release, x_hf_screening_type ,&
          x_hf_set_distr_energy, proc_assoc_energy,&
          x_hf_set_distr_forces, proc_assoc_forces, x_hf_create_neighbor_cells,&
          x_hf_cell_type

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

  TYPE x_hf_cell_type
    INTEGER                                  :: cell(3)
  END TYPE
  
  TYPE x_hf_type
    REAL(dp), DIMENSION(:), POINTER          :: x_buffer
    REAL(dp), DIMENSION(:, :, :, :), &
              POINTER                        :: x_screening_buffer
    REAL(dp), DIMENSION(:, :, :, :, :, :), &
              POINTER                        :: x_screening_buffer_periodic
    INTEGER                                  :: number_of_elements
    TYPE(x_hf_screening_type)                :: screening_parameter
    INTEGER(KIND=int_8), DIMENSION(:),&
                             POINTER         :: integral_distribution_energy
    INTEGER(KIND=int_8), DIMENSION(:),&
                             POINTER         :: integral_distribution_forces
    TYPE(x_hf_cell_type), DIMENSION(:), &
               POINTER                       :: neighbor_cells
    LOGICAL                                  :: do_periodic
    INTEGER                                  :: number_of_k_loops
  END TYPE

  CONTAINS

  SUBROUTINE x_hf_create(x_data,hf_section,natom,max_set,error)
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(section_vals_type), POINTER         :: hf_section
    INTEGER, INTENT(IN)                      :: natom, max_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: max_memory, nkloops, &
                                                number_of_precalc_integrals, &
                                                pbc_shells, stat, TYPE
    LOGICAL                                  :: do_periodic
    REAL(dp)                                 :: omega
    TYPE(section_vals_type), POINTER         :: hf_pbc_section

    CALL section_vals_val_get(hf_section,"MAX_MEMORY",i_val=max_memory,error=error)
    number_of_precalc_integrals = max_memory * 1024*1024/dp_size
    ALLOCATE(x_data,STAT=stat)
    ALLOCATE(x_data%x_buffer(number_of_precalc_integrals),STAT=stat)
    x_data%number_of_elements = number_of_precalc_integrals
    CALL section_vals_val_get(hf_section,"screening_type",i_val=TYPE,error=error)
    CALL section_vals_val_get(hf_section,"omega",r_val=omega,error=error)
 
    x_data%screening_parameter%screening_type = TYPE
    x_data%screening_parameter%omega = omega

    hf_pbc_section => section_vals_get_subs_vals(hf_section,"PERIODIC",error=error)
    CALL section_vals_get(hf_pbc_section,explicit=do_periodic,error=error)
    x_data%do_periodic = do_periodic
    x_data%number_of_k_loops = 1
    IF(do_periodic) THEN
      hf_pbc_section => section_vals_get_subs_vals(hf_section,"PERIODIC",error=error)
      CALL section_vals_val_get(hf_pbc_section,"NUMBER_OF_SHELLS",i_val=pbc_shells,error=error)
      ALLOCATE(x_data%neighbor_cells(1),STAT=stat)
      CALL x_hf_create_neighbor_cells(x_data,pbc_shells, error=error)
      ALLOCATE(x_data%x_screening_buffer_periodic(natom,natom,max_set,max_set,&
                                                  SIZE(x_data%neighbor_cells),SIZE(x_data%neighbor_cells)))
      CALL section_vals_val_get(hf_pbc_section,"NUMBER_OF_K_LOOPS",i_val=nkloops,error=error)
      x_data%number_of_k_loops = nkloops
    ELSE
      ALLOCATE(x_data%x_screening_buffer(natom,natom,max_set,max_set))
    END IF

    ALLOCATE(x_data%integral_distribution_energy(1))
    ALLOCATE(x_data%integral_distribution_forces(1))
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
    
    IF( x_data%do_periodic ) THEN
      DEALLOCATE(x_data%x_screening_buffer_periodic,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(x_data%neighbor_cells,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
      DEALLOCATE(x_data%x_screening_buffer,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    DEALLOCATE(x_data%integral_distribution_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(x_data%integral_distribution_forces,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(x_data,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE x_hf_release

  SUBROUTINE x_hf_set_distr_energy(distr, x_data, error)
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: distr
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_set_distr_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    DEALLOCATE(x_data%integral_distribution_energy,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%integral_distribution_energy(SIZE(distr)))
    x_data%integral_distribution_energy = distr
  END SUBROUTINE x_hf_set_distr_energy

  FUNCTION proc_assoc_energy(counter, distribution)
    INTEGER(KIND=int_8)                      :: counter
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: distribution
    LOGICAL                                  :: proc_assoc_energy

    INTEGER                                  :: ipos, mysize

    mysize = SIZE(distribution,1)/2
    
    proc_assoc_energy = .FALSE.
    DO ipos = 1,mysize
      IF(counter>=distribution(2*ipos-1) .AND. counter<=distribution(2*ipos) ) THEN
        proc_assoc_energy = .TRUE.
        RETURN
      END IF
    END DO
  END FUNCTION proc_assoc_energy

  SUBROUTINE x_hf_set_distr_forces(distr, x_data, error)
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: distr
    TYPE(x_hf_type), POINTER                 :: x_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_set_distr_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    DEALLOCATE(x_data%integral_distribution_forces,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%integral_distribution_forces(SIZE(distr)))
    x_data%integral_distribution_forces = distr
  END SUBROUTINE x_hf_set_distr_forces

  FUNCTION proc_assoc_forces(counter, distribution)
    INTEGER(KIND=int_8)                      :: counter
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: distribution
    LOGICAL                                  :: proc_assoc_forces

    INTEGER                                  :: ipos, mysize

    mysize = SIZE(distribution,1)/2

    proc_assoc_forces = .FALSE.
    DO ipos = 1,mysize
      IF(counter>=distribution(2*ipos-1) .AND. counter<=distribution(2*ipos) ) THEN
        proc_assoc_forces = .TRUE.
        RETURN
      END IF
    END DO
  END FUNCTION proc_assoc_forces
  
  SUBROUTINE x_hf_create_neighbor_cells(x_data, pbc_shells, error)
    TYPE(x_hf_type), POINTER                 :: x_data
    INTEGER, INTENT(IN)                      :: pbc_shells
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'x_hf_create_neighbor_cells', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, m(3), stat, &
                                                total_number_of_cells
    LOGICAL                                  :: failure

    failure = .FALSE.
    total_number_of_cells = 0
    DO i = 0,pbc_shells
      total_number_of_cells = total_number_of_cells + count_cells(i)
    END DO
    
    DEALLOCATE(x_data%neighbor_cells,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(x_data%neighbor_cells(total_number_of_cells), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    
    m = 0
    i = 1
    DO WHILE(SUM(m**2)<=pbc_shells)
       x_data%neighbor_cells(i)%cell = m
       CALL next_image_cell(m)
       i=i+1
    ENDDO
  END SUBROUTINE x_hf_create_neighbor_cells

   FUNCTION count_cells(shell)
    INTEGER, INTENT(IN)                      :: shell
    INTEGER                                  :: count_cells

    INTEGER                                  :: i, j, k

    count_cells=0
    DO i=-shell,shell
      DO j=-shell,shell
        DO k=-shell,shell
          IF ((i**2+j**2+k**2==shell)) count_cells=count_cells+1
        ENDDO
      ENDDO
    ENDDO
  END FUNCTION count_cells

  SUBROUTINE next_image_cell(m)
    INTEGER                                  :: m(3)

    INTEGER                                  :: i, j, k, shell
    LOGICAL                                  :: found

    found=.FALSE.
    shell=SUM(m**2)
    outer: DO
           DO i=-shell,shell
           DO j=-shell,shell
    inner: DO k=-shell,shell
              IF (.NOT.(i**2+j**2+k**2==shell)) CYCLE inner
              IF (found) THEN
                  m=(/i,j,k/)
                  EXIT outer
              ENDIF
              IF (ALL(M.EQ.(/i,j,k/))) found=.TRUE.
           ENDDO inner
           ENDDO
           ENDDO
           shell=shell+1
           ENDDO outer
  END SUBROUTINE next_image_cell

END MODULE x_hf


