!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility methods to build 3-center integral tensors of various types.
! **************************************************************************************************

MODULE qs_3c_tensors_types
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_iterator_p_type,&
                                              neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_3c_tensors_types'

   PUBLIC :: distribution_3d_type, neighbor_list_3c_type, neighbor_list_3c_iterator_type

   INTEGER, PARAMETER, PUBLIC :: symmetric_none = 0, symmetric_ij = 1, symmetric_jk = 2, symmetrik_ik = 3, symmetric_ijk = 4

   TYPE distribution_3d_type
      TYPE(distribution_2d_type), POINTER :: dist_2d_1 => NULL(), dist_2d_2 => NULL()
      INTEGER :: comm_3d, comm_2d_1, comm_2d_2
      LOGICAL :: owns_comm
   END TYPE distribution_3d_type

   TYPE neighbor_list_3c_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: ij_list => NULL(), jk_list => NULL()
      INTEGER :: sym
      TYPE(distribution_3d_type) :: dist_3d
      LOGICAL :: owns_dist
   END TYPE

   TYPE neighbor_list_3c_iterator_type
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iter_ij => NULL()
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iter_jk => NULL()
      INTEGER                                                    :: iter_level
      TYPE(neighbor_list_3c_type)                                :: ijk_nl
   END TYPE

END MODULE
