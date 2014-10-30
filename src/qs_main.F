!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Perform a LCGTO-SCF-DF calculation using the Gaussian plane wave (GPW)
!>      method
!> \par Literature
!>      G. Lippert, J. Hutter, M. Parrinello, Mol. Phys. 92, 477 (1997)
!>      G. Lippert, J. Hutter, M. Parrinello, Theor. Chem. Acc. 103, 124 (1999)
!>      M. Krack, M. Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)
!> \par History
!>      none
!> \author Matthias Krack (03.04.2001)
! *****************************************************************************
MODULE qs_main
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_type
  USE f77_blas
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_p_type,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE qs_environment,                  ONLY: qs_init
  USE qs_environment_types,            ONLY: qs_env_create,&
                                             qs_env_release,&
                                             qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

    IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_main'

  PUBLIC :: quickstep_create_force_env

CONTAINS

! *****************************************************************************
!> \brief Main driver routine for a Quickstep calculation.
!>      QS_INIT:: Read the input and the database files and perform further
!>                initializations for the setup of the QUICKSTEP environment
!> \param force_env ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param subsys ...
!> \param cell ...
!> \param cell_ref ...
!> \param qmmm ...
!> \param qmmm_env_qm ...
!> \param qmmm_periodic ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \param error ...
!> \par History
!>     Creation (23.01.2002,MK)
!>     Modified for RESPA (07.2006)
! *****************************************************************************
  SUBROUTINE quickstep_create_force_env(force_env,root_section,para_env,globenv,&
       subsys,cell,cell_ref,qmmm,qmmm_env_qm,force_env_section,&
       subsys_section,use_motion_section,error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    TYPE(cell_type), OPTIONAL, POINTER       :: cell, cell_ref
    LOGICAL, INTENT(IN), OPTIONAL            :: qmmm
    TYPE(qmmm_env_qm_type), OPTIONAL, &
      POINTER                                :: qmmm_env_qm
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'quickstep_create_force_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(force_env_p_type), DIMENSION(:), &
      POINTER                                :: sub_force_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(qs_env, sub_force_env)
    !Reference qs_env
    CALL qs_env_create(qs_env,globenv,error=error)
    CALL qs_init(qs_env, para_env, globenv, root_section, cp_subsys=subsys, cell=cell, cell_ref=cell_ref,&
         qmmm=qmmm, qmmm_env_qm=qmmm_env_qm, force_env_section=force_env_section,&
         subsys_section=subsys_section, use_motion_section=use_motion_section, error=error)
    CALL force_env_create(force_env,root_section, qs_env=qs_env, para_env=para_env, globenv=globenv,&
         sub_force_env=sub_force_env, force_env_section=force_env_section,error=error)
    CALL qs_env_release(qs_env,error=error)

    CALL timestop(handle)
  END SUBROUTINE quickstep_create_force_env

END MODULE qs_main
