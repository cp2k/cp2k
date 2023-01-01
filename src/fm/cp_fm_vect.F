!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routine to handle vectors of full matrixes
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE cp_fm_vect
   USE cp_fm_types, ONLY: cp_fm_p_type, &
                          cp_fm_release, cp_fm_type
#include "../base/base_uses.f90"

   #:mute
      #:set types = [("cp_fm_type", "array", ""), ("cp_fm_p_type", "p_type","%matrix")]
      #:set attributes = [("ALLOCATABLE", "alloc", "ALLOCATED"), ("POINTER", "pointer", "ASSOCIATED")]
   #:endmute

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_vect'

   PUBLIC :: cp_fm_vect_dealloc

   INTERFACE cp_fm_vect_dealloc
      #:for typename, shortname, appendix in types
         #:for attr, shortattr, create in attributes
            MODULE PROCEDURE cp_fm_m1_dealloc_${shortname}$_${shortattr}$
         #:endfor
      #:endfor
   END INTERFACE
!***
CONTAINS

   #:for typename, shortname, appendix in types
      #:for attr, shortattr, create in attributes
! **************************************************************************************************
!> \brief deallocate an array of pointers to blacs matrixes
!> \param matrices the array of matrices to deallocate
!> \par History
!>      07.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
         SUBROUTINE cp_fm_m1_dealloc_${shortname}$_${shortattr}$ (matrices)
            TYPE(${typename}$), DIMENSION(:), ${attr}$          :: matrices

            INTEGER                                            :: i

            IF (${create}$ (matrices)) THEN
               DO i = 1, SIZE(matrices)
                  CALL cp_fm_release(matrices(i) ${appendix}$)
                  #:if typename=="cp_fm_p_type"
                     DEALLOCATE (matrices(i)%matrix)
                  #:endif
               END DO
               DEALLOCATE (matrices)
               #:if attr=="POINTER"
                  NULLIFY (matrices)
               #:endif
            END IF
         END SUBROUTINE cp_fm_m1_dealloc_${shortname}$_${shortattr}$
      #:endfor
   #:endfor

END MODULE cp_fm_vect
