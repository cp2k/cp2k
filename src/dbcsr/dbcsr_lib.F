!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Routines that affect the DBCSR library as a whole
!> \author  Ole Schuett
!> \date    10-2014
! **************************************************************************************************

MODULE dbcsr_lib
   USE dbcsr_multiply_api,              ONLY: dbcsr_multiply_clear_mempools,&
                                              dbcsr_multiply_lib_finalize,&
                                              dbcsr_multiply_lib_init
   USE kinds,                           ONLY: int_1_size,&
                                              int_2_size,&
                                              int_4_size,&
                                              int_8_size
#include "../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_lib'

   PUBLIC :: dbcsr_init_lib, dbcsr_finalize_lib, dbcsr_clear_mempools

   LOGICAL, PRIVATE, SAVE :: is_initialized = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Initialize the DBCSR library
!>
!> Prepares the DBCSR library for use.
! **************************************************************************************************
   SUBROUTINE dbcsr_init_lib()

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_init_lib', routineP = moduleN//':'//routineN

      INTEGER                                            :: error_handle

!n_stack_buffers,mem_type, n_stack_mem_regions, stack_size, nstacks,&
!INTEGER, DIMENSION(3)                    :: nxstacks
!   ---------------------------------------------------------------------------
!TODO: problem: init/finalize are called by cp2k_runs AND f77_interface

      IF (is_initialized) RETURN
      CALL timeset(routineN, error_handle)
      !
      IF (int_1_size /= 1) &
         CPABORT("Incorrect assumption of an 8-bit integer size!")
      IF (int_2_size /= 2) &
         CPABORT("Incorrect assumption of a 16-bit integer size!")
      IF (int_4_size /= 4) &
         CPABORT("Incorrect assumption of a 32-bit integer size!")
      IF (int_8_size /= 8) &
         CPABORT("Incorrect assumption of a 64-bit integer size!")

!$OMP     PARALLEL DEFAULT(NONE)
      CALL dbcsr_multiply_lib_init()
!$OMP     END PARALLEL

      is_initialized = .TRUE.
      CALL timestop(error_handle)
   END SUBROUTINE dbcsr_init_lib

! **************************************************************************************************
!> \brief Finalize the DBCSR library
!>
!> Cleans up after the DBCSR library.  Used to deallocate persistent objects.
!> \param group ...
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE dbcsr_finalize_lib(group, output_unit)
      INTEGER, INTENT(IN)                                :: group, output_unit

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_finalize_lib', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: error_handle

!TODO: problem: init/finalize are called by cp2k_runs AND f77_interface

      IF (.NOT. is_initialized) RETURN
      CALL timeset(routineN, error_handle)

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,T2,A)") REPEAT("-", 79)
         WRITE (UNIT=output_unit, FMT="(T2,A,T80,A)") "-", "-"
         WRITE (UNIT=output_unit, FMT="(T2,A,T35,A,T80,A)") "-", "DBCSR STATISTICS", "-"
         WRITE (UNIT=output_unit, FMT="(T2,A,T80,A)") "-", "-"
         WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("-", 79)
      END IF

!$OMP     PARALLEL DEFAULT(NONE) SHARED(output_unit, group)
      CALL dbcsr_multiply_lib_finalize(group, output_unit)
!$OMP     END PARALLEL
      IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T2,A)") REPEAT("-", 79)

      is_initialized = .FALSE.
      CALL timestop(error_handle)
   END SUBROUTINE dbcsr_finalize_lib

! **************************************************************************************************
!> \brief  Deallocate memory contained in mempools
! **************************************************************************************************
   SUBROUTINE dbcsr_clear_mempools()

!$OMP     PARALLEL DEFAULT(NONE)
      CALL dbcsr_multiply_clear_mempools()
!$OMP     END PARALLEL
   END SUBROUTINE dbcsr_clear_mempools

END MODULE dbcsr_lib
