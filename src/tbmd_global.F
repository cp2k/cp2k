!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 2000  CP2K developers group                                  !
!------------------------------------------------------------------------------!
!
!!    Definition of some global variables for the TBMD code
    MODULE tbmd_global
!
!------------------------------------------------------------------------------!
      USE kinds, ONLY : dbl
      USE global_types, ONLY : global_environment_type
      USE tbmd_types, ONLY : tbatom_param_type, tb_hopping_type
!
      IMPLICIT NONE
!
      PRIVATE
      PUBLIC :: set_tbmd_global, tbmdpar, tbatom, tbhop
!
      TYPE (global_environment_type) :: tbmdpar
!
      TYPE (tbatom_param_type), DIMENSION (:), ALLOCATABLE :: tbatom
      TYPE (tb_hopping_type), DIMENSION (:,:), ALLOCATABLE :: tbhop
!
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
      SUBROUTINE set_tbmd_global(glopar)
        IMPLICIT NONE
        TYPE (global_environment_type) :: glopar

        tbmdpar = glopar
      END SUBROUTINE set_tbmd_global
!------------------------------------------------------------------------------!
    END MODULE tbmd_global
!------------------------------------------------------------------------------!

