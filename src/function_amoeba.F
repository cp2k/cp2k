!------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations          !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart               !
!------------------------------------------------------------------------------!
!*****************************************************************
!!
!!   NAME
!!     function_amoeba
!!
!*****************************************************************
!------------------------------------------------------------------------------!

    FUNCTION func(x,struc,inter,thermo,simpar,fit_param,box_ref, &
                  finfo,param,pstat,ewald_param, globenv)
      USE kinds, ONLY: dp
      USE structure_types, ONLY : structure_type, interaction_type
      USE pair_potential, ONLY : potentialparm_type
      USE particle_types, ONLY : particle_prop_type 
      USE md, ONLY : simulation_parameters_type, thermodynamic_type
      USE simulation_cell, ONLY : cell_type
      USE fit_parameters_types, ONLY : fit_parameter_type, linklist_parameters_type, &
                                       fit_info_type
      USE ewald_parameters_types, ONLY : ewald_parameters_type
      USE global_types, ONLY : global_environment_type
      USE fit_pair_potential, ONLY : local
      IMPLICIT NONE
      TYPE (structure_type), INTENT (INOUT) :: struc
      TYPE (interaction_type), INTENT (INOUT) :: inter
      TYPE (thermodynamic_type), INTENT (INOUT) :: thermo
      TYPE (simulation_parameters_type), INTENT (IN) :: simpar
      REAL(KIND=dp), INTENT (IN), DIMENSION (:) :: x
      TYPE (fit_parameter_type), INTENT(INOUT) :: fit_param
      TYPE (cell_type) , POINTER, DIMENSION (:) :: box_ref
      TYPE (fit_info_type), DIMENSION (:,:), POINTER :: finfo
      TYPE (linklist_parameters_type ), INTENT(INOUT), DIMENSION ( :, : ) :: param 
      TYPE (particle_prop_type ), DIMENSION ( : ), POINTER :: pstat 
      TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param  
      TYPE ( global_environment_type), INTENT ( IN ) :: globenv
      REAL(KIND=dp) :: func

      func = local(x,struc,inter,thermo,simpar,fit_param,box_ref, &
                  finfo,param,pstat,ewald_param, globenv)
      
   END FUNCTION func   


