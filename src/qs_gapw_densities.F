!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_gapw_densities [1.0] *
!!
!!   NAME
!!     qs_gapw_densities
!!
!!   FUNCTION
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     
!!
!!
!!   SOURCE
!****************************************************************************
MODULE qs_gapw_densities
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_error_handling,               ONLY: cp_error_type
  USE kinds,                           ONLY: dp
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_rho0_methods,                 ONLY: calculate_rho0_atom
  USE qs_rho0_types,                   ONLY: rho0_mpole_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='qs_gapw_densities'

  PUBLIC :: prepare_gapw_den 


!******************************************************************************


CONTAINS

!******************************************************************************

 SUBROUTINE prepare_gapw_den(qs_env,error)

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(qs_environment_type), INTENT(INOUT) :: qs_env

    CHARACTER(LEN=*), PARAMETER         :: routine_name = "prepare_gapw_den"

    TYPE(atomic_kind_type), DIMENSION(:), &
                                 POINTER     :: atomic_kind_set
    TYPE(atomic_kind_type),      POINTER     :: atom_kind
    TYPE(grid_atom_type),      POINTER       :: grid_atom
    TYPE(rho0_mpole_type),      POINTER      :: rho0_mpole
    INTEGER                                  :: natom, nkind , ikind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: paw_atom

!******************************************************************************

    NULLIFY(atomic_kind_set)
    NULLIFY(atom_kind)
    NULLIFY(atom_list)
    NULLIFY(grid_atom)
    NULLIFY(rho0_mpole)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    rho0_mpole=rho0_mpole)

    nkind = SIZE(atomic_kind_set)

    DO ikind = 1,nkind
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
                           grid_atom=grid_atom,natom=natom,&
                           paw_atom=paw_atom)

!     Calculate rho1_h and rho1_s on the radial gris centered on the atomic position
      IF(paw_atom) &
        CALL calculate_rho_atom(qs_env,atom_kind,atom_list,grid_atom,natom)

!     Calculate rho0_h and rho0_s on the radial gris centered on the atomic position
      CALL calculate_rho0_atom(qs_env,rho0_mpole,atom_list,grid_atom,&
                               paw_atom,natom,ikind)
    ENDDO
!    stop 'cpc'
 END SUBROUTINE
!******************************************************************************

 END MODULE  qs_gapw_densities 
