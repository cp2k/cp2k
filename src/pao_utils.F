!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Helper routines shared by multiple PAO modules
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_utils
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_iterator,&
                                             cp_dbcsr_iterator_blocks_left,&
                                             cp_dbcsr_iterator_next_block,&
                                             cp_dbcsr_iterator_start,&
                                             cp_dbcsr_iterator_stop,&
                                             cp_dbcsr_type
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_utils'

  PUBLIC :: fill_randomly

CONTAINS

! *****************************************************************************
!> \brief Fills a DBCSR matrix with random noise.
!> \param matrix ...
! *****************************************************************************
  SUBROUTINE fill_randomly(matrix)
    TYPE(cp_dbcsr_type)                      :: matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'fill_randomly', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, j
    REAL(dp), DIMENSION(:, :), POINTER       :: block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(rng_stream_type), POINTER           :: rng_stream

    NULLIFY(rng_stream)
    CALL timeset(routineN,handle)

    CALL cp_dbcsr_iterator_start(iter, matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block)
       CALL create_rng_stream(rng_stream, name="pao")
       DO i=1, SIZE(block, 1)
          DO j=1, SIZE(block, 2)
             block(i,j) = next_random_number(rng_stream)
          ENDDO
       ENDDO
       CALL delete_rng_stream(rng_stream)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE fill_randomly

END MODULE pao_utils
