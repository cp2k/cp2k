!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Some auxiliary functions and subroutines needed for HFX calculations 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_helpers 
  
  USE kinds,                           ONLY: int_8
  USE mathconstants
  USE orbital_pointers
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC FNV64, count_cells, next_image_cell

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_helpers'

!***
  
  CONTAINS

! *****************************************************************************
!> \brief Calculates random number using hash-table 
!> \note
!>      - FNV-1a hash from
!>          http://www.isthe.com/chongo/tech/comp/fnv/ 
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  FUNCTION FNV64(val) RESULT(hash)
    INTEGER(int_8), INTENT(IN)               :: val
    INTEGER(int_8)                           :: hash

    INTEGER                                  :: i
    INTEGER(int_8)                           :: tmp

    hash=1099511628211_int_8
    tmp=val
    DO i=1,8
      hash=IOR(IAND(hash,NOT(255_int_8)),IAND(IEOR(IAND(tmp,255_int_8),hash),255_int_8))
      tmp=ISHFT(tmp,-8)
      hash = hash + ISHFT(hash,1) + ISHFT(hash,4)+ISHFT(hash,5)+ISHFT(hash,7)+ISHFT(hash,8)+ISHFT(hash,40)      
    ENDDO
  END FUNCTION FNV64


! *****************************************************************************
!> \brief - Auxiliary function for creating periodic neighbor cells
!> \param shell number of shells in each coordinate direction
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
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

! *****************************************************************************
!> \brief - Auxiliary function for creating periodic neighbor cells
!> \param shell shell for which we need the next, stored as number of shifts
!>               in each direction
!> \par History
!>      09.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
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

! *****************************************************************************

END MODULE hfx_helpers 
