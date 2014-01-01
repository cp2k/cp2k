!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a Kim-Gordon-like partitioning into molecular subunits 
!> \par History
!>       2012.06 created [Martin Haeufel]
!> \author Martin Haeufel and Florian Schiffmann
! *****************************************************************************
MODULE kg_correction
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE f77_blas
  USE kg_environment_types,            ONLY: kg_environment_type
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                             qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_release,&
                                             qs_rho_type
  USE qs_vxc,                          ONLY: qs_vxc_create
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_correction'

  PUBLIC kg_ekin_subset

CONTAINS

! *****************************************************************************
!> \brief Calculates the subsystem Hohenberg-Kohn kinetic energy and the forces
!> \par History
!>       2012.06 created [Martin Haeufel]
!> \author Martin Haeufel and Florian Schiffmann
! *****************************************************************************
  SUBROUTINE kg_ekin_subset(qs_env, ks_matrix, gapw, gapw_xc, v_rspace_new, &
       v_tau_rspace, ekin_mol, calc_force, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    LOGICAL, INTENT(in), OPTIONAL            :: gapw, gapw_xc
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, v_tau_rspace
    REAL(KIND=dp), INTENT(out)               :: ekin_mol
    LOGICAL                                  :: calc_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_ekin_subset', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, isub, natom, &
                                                nspins
    LOGICAL                                  :: my_gapw, my_gapw_xc
    REAL(KIND=dp)                            :: ekin_imol
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: density_matrix, matrix_h, &
                                                my_mat
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(kg_environment_type), POINTER       :: kg_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: vxc_rho, vxc_tau
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho_struct

    CALL timeset(routineN, handle)

    NULLIFY(vxc_rho, vxc_tau, rho_struct)
 
    kg_env => qs_env%kg_env

    my_gapw = .FALSE.
    IF (PRESENT(gapw)) my_gapw = gapw

    my_gapw_xc = .FALSE.
    IF (PRESENT(gapw_xc)) my_gapw_xc = gapw_xc

    ekin_mol = 0.0_dp

    ! get set of molecules, natom, dft_control, pw_env
    CALL get_qs_env(qs_env, matrix_h=matrix_h, natom=natom, dft_control=dft_control, pw_env=pw_env, error=error)

    nspins = dft_control%nspins

    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)

    ! get the density matrix
    density_matrix => qs_env%rho%rho_ao

    ! allocate and initialize the density 
    CALL qs_rho_rebuild(rho_struct, qs_env, gapw_xc=my_gapw_xc, rebuild_grids=.TRUE., error=error)

    ! loop over all subsets
    DO isub=1,kg_env%nsubsets

      my_mat => rho_struct%rho_ao

      ! set the density matrix to the blocked matrix  
      rho_struct%rho_ao =>  density_matrix  ! blocked_matrix

      ! for gapw the task_list_soft would be needed
      CALL cp_assert(.NOT.(my_gapw.OR.my_gapw_xc), cp_fatal_level, cp_assertion_failed, routineP, &
                     " KG is currently not compatible with all electron calculations.")  

      ! calculate the densities for the given blocked density matrix - pass the subset task_list
      CALL qs_rho_update_rho(rho_struct, qs_env, & 
                             !task_list_external=qs_env%task_list, &
                             task_list_external=kg_env%subset(isub)%task_list, &
                             error=error)

      ekin_imol=0.0_dp

      ! calc Hohenberg-Kohn kin. energy of the density corresp. to the remaining molecular block(s)
      CALL qs_vxc_create(qs_env=qs_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                        vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol, error=error)

      rho_struct%rho_ao => my_mat

      ekin_mol = ekin_mol + ekin_imol

      DO ispin=1,nspins

        vxc_rho(ispin)%pw%cr3d = -vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol*1.0_dp
 
        CALL integrate_v_rspace(v_rspace=vxc_rho(ispin), &
             p=density_matrix(ispin), & 
             h=ks_matrix(ispin), & 
             qs_env=qs_env, &
             calculate_forces=calc_force, &
             gapw=gapw, &
             task_list_external=kg_env%subset(isub)%task_list, &
             error=error)

        ! clean up vxc_rho
        CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw, error=error)

      END DO
      DEALLOCATE(vxc_rho)

      IF (ASSOCIATED(vxc_tau)) THEN
      ! this case has NOT been tested yet...
    
        CALL cp_assert(ASSOCIATED(vxc_tau), cp_fatal_level, cp_assertion_failed, routineP, &
                       " KG is currently not compatible with TPSS, Becke-Roussel, or other"//& 
                       " Meta-GGA functionals.")
    
        DO ispin=1,nspins

          ! clean up vxc_tau
          CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_tau(ispin)%pw, error=error)
       
        END DO

        DEALLOCATE(vxc_tau)

      END IF 

    END DO

    ! clean up rho_struct
    CALL qs_rho_release(rho_struct, error) 

    CALL timestop(handle)

  END SUBROUTINE kg_ekin_subset

END MODULE kg_correction
