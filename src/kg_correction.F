!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2012  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for a Kim-Gordon-like partitioning into molecular subunits 
!> \par History
!>       2012.06 created [Martin Haeufel]
!> \author Martin Haeufel and Florian Schiffmann
! *****************************************************************************
MODULE kg_correction
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_copy, cp_dbcsr_copy_into_existing, &
       cp_dbcsr_create, cp_dbcsr_finalize, cp_dbcsr_init, cp_dbcsr_release, &
       cp_dbcsr_reserve_blocks, cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_allocate_matrix_set,&
                                             cp_dbcsr_deallocate_matrix_set
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
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

  PUBLIC kg_ekin_subset, kg_ekin_mol

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

    INTEGER                                  :: handle, imol, ispin, isub, &
                                                jmol, natom, nspins
    LOGICAL                                  :: my_gapw, my_gapw_xc
    REAL(KIND=dp)                            :: ekin_imol
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: blocked_matrix, &
                                                density_matrix, &
                                                mat_ks_blocked, mat_tmp2, &
                                                matrix_h, my_mat
    TYPE(cp_dbcsr_type), POINTER             :: mat_tmp
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

      ! allocate dbcsrc matrix sets
      CALL cp_dbcsr_allocate_matrix_set(blocked_matrix, nspins, error=error)

      CALL cp_dbcsr_allocate_matrix_set(mat_ks_blocked, nspins, error=error)

      CALL cp_dbcsr_allocate_matrix_set(mat_tmp2, nspins, error=error)

      DO ispin=1,nspins

        ! build density matrix, all blocks not part of the current molecule set to zero        
        ALLOCATE(mat_tmp, mat_tmp2(ispin)%matrix, blocked_matrix(ispin)%matrix, mat_ks_blocked(ispin)%matrix)

        CALL cp_dbcsr_init(mat_tmp, error=error)

        CALL cp_dbcsr_init(blocked_matrix(ispin)%matrix, error=error)

        CALL cp_dbcsr_init(mat_ks_blocked(ispin)%matrix, error=error)

        CALL cp_dbcsr_init(mat_tmp2(ispin)%matrix, error=error)

        CALL cp_dbcsr_create(mat_tmp, 'MAT_TMP', template=density_matrix(ispin)%matrix, error=error)
  
        CALL cp_dbcsr_copy(blocked_matrix(ispin)%matrix, density_matrix(ispin)%matrix, 'BLOCKED MATRIX', error=error)

        CALL cp_dbcsr_copy(mat_tmp2(ispin)%matrix, density_matrix(ispin)%matrix, 'MAT_TMP2', error=error)
 
        ! zero all blocks
        CALL cp_dbcsr_set(blocked_matrix(ispin)%matrix, 0.0_dp, error=error)

        CALL cp_dbcsr_set(mat_tmp2(ispin)%matrix, 0.0_dp, error=error)

        ! loop over all molecules in the current subset
        DO jmol=1,SIZE(kg_env%subset(isub)%molecules)

          imol = kg_env%subset(isub)%molecules(jmol)

          CALL cp_dbcsr_reserve_blocks(mat_tmp, kg_env%molecular_blocks(imol)%rows, &
                                       kg_env%molecular_blocks(imol)%cols, error=error) 
  
        END DO

        CALL cp_dbcsr_finalize(mat_tmp, error=error)
  
        CALL cp_dbcsr_copy_into_existing(mat_tmp, density_matrix(ispin)%matrix, error)
  
        CALL cp_dbcsr_copy_into_existing(blocked_matrix(ispin)%matrix, mat_tmp, error)

        CALL cp_dbcsr_copy(mat_ks_blocked(ispin)%matrix, mat_tmp, 'MAT_KS_BLOCKED', error=error)

        CALL cp_dbcsr_release(mat_tmp, error)
  
        DEALLOCATE(mat_tmp)
 
      END DO

      my_mat => rho_struct%rho_ao

      ! set the density matrix to the blocked matrix  
      rho_struct%rho_ao => blocked_matrix

      ! calculate the densities for the given blocked density matrix
      CALL qs_rho_update_rho(rho_struct, qs_env, error=error)

      ekin_imol=0.0_dp

      ! calc Hohenberg-Kohn kin. energy of the density corresp. to the remaining molecular block(s)
      CALL qs_vxc_create(qs_env=qs_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                        vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol, error=error)

      rho_struct%rho_ao => my_mat

      ekin_mol = ekin_mol + ekin_imol

      DO ispin=1,nspins

        vxc_rho(ispin)%pw%cr3d = -vxc_rho(ispin)%pw%cr3d*vxc_rho(ispin)%pw%pw_grid%dvol*1.0_dp
 
        CALL integrate_v_rspace(v_rspace=vxc_rho(ispin), &
             p=blocked_matrix(ispin), &
             h=mat_tmp2(ispin), &
             qs_env=qs_env, &
             calculate_forces=calc_force, &
             gapw=gapw, &
             task_list_external=kg_env%subset(isub)%task_list, &
             error=error)

        CALL cp_dbcsr_copy_into_existing(mat_ks_blocked(ispin)%matrix, mat_tmp2(ispin)%matrix, error)

        CALL cp_dbcsr_add(ks_matrix(ispin)%matrix,mat_ks_blocked(ispin)%matrix,1.0_dp,1.0_dp,error=error)

        ! clean up vxc_rho
        CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw, error=error)

      END DO
      DEALLOCATE(vxc_rho)

      IF (ASSOCIATED(vxc_tau)) THEN
      ! this case has NOT been tested yet...
       
        DO ispin=1,nspins

          ! zero mat_tmp2
          CALL cp_dbcsr_set(mat_tmp2(ispin)%matrix, 0.0_dp, error=error)

          vxc_tau(ispin)%pw%cr3d = -vxc_tau(ispin)%pw%cr3d*vxc_tau(ispin)%pw%pw_grid%dvol*1.0_dp

          CALL integrate_v_rspace(v_rspace=vxc_tau(ispin),&
               p=blocked_matrix(ispin), h=mat_tmp2(ispin),&
               qs_env=qs_env, &
               calculate_forces=calc_force,&
               compute_tau=.TRUE.,&
               task_list_external=kg_env%subset(isub)%task_list, &
               gapw=gapw, error=error)
    
          CALL cp_dbcsr_copy_into_existing(mat_ks_blocked(ispin)%matrix, mat_tmp2(ispin)%matrix, error)

          CALL cp_dbcsr_add(ks_matrix(ispin)%matrix,mat_ks_blocked(ispin)%matrix,1.0_dp,1.0_dp,error=error)

          ! clean up vxc_tau
          CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_tau(ispin)%pw, error=error)
       
        END DO

        DEALLOCATE(vxc_tau)

      END IF 

      ! deallocate dbcsr matrix sets
      CALL cp_dbcsr_deallocate_matrix_set(blocked_matrix, error=error)    
 
      CALL cp_dbcsr_deallocate_matrix_set(mat_ks_blocked, error=error)  

      CALL cp_dbcsr_deallocate_matrix_set(mat_tmp2, error=error)    
    END DO

    ! clean up rho_struct
    CALL qs_rho_release(rho_struct, error) 

    CALL timestop(handle)

  END SUBROUTINE kg_ekin_subset


  SUBROUTINE kg_ekin_mol(qs_env, ks_matrix, gapw, gapw_xc, v_rspace_new, &
       v_tau_rspace, ekin_mol, calc_force, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    LOGICAL, INTENT(in), OPTIONAL            :: gapw, gapw_xc
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, v_tau_rspace
    REAL(KIND=dp), INTENT(out)               :: ekin_mol
    LOGICAL                                  :: calc_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'kg_ekin_mol', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, imol, ispin, natom, &
                                                nspins
    LOGICAL                                  :: my_gapw, my_gapw_xc
    REAL(KIND=dp)                            :: ekin_imol
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: blocked_matrix, &
                                                density_matrix, &
                                                mat_ks_blocked, mat_tmp2, &
                                                matrix_h, my_mat
    TYPE(cp_dbcsr_type), POINTER             :: mat_tmp
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

    ! loop over all molecules
    DO imol=1,SIZE(kg_env%molecule_set)

      ! allocate dbcsrc matrix sets
      CALL cp_dbcsr_allocate_matrix_set(blocked_matrix, nspins, error=error)

      CALL cp_dbcsr_allocate_matrix_set(mat_ks_blocked, nspins, error=error)

      CALL cp_dbcsr_allocate_matrix_set(mat_tmp2, nspins, error=error)

      DO ispin=1,nspins

        ! build density matrix, all blocks not part of the current molecule set to zero        
        ALLOCATE(mat_tmp, mat_tmp2(ispin)%matrix, blocked_matrix(ispin)%matrix, mat_ks_blocked(ispin)%matrix)

        CALL cp_dbcsr_init(mat_tmp, error=error)

        CALL cp_dbcsr_init(blocked_matrix(ispin)%matrix, error=error)

        CALL cp_dbcsr_init(mat_ks_blocked(ispin)%matrix, error=error)

        CALL cp_dbcsr_init(mat_tmp2(ispin)%matrix, error=error)

        CALL cp_dbcsr_create(mat_tmp, 'MAT_TMP', template=density_matrix(ispin)%matrix, error=error)
  
        CALL cp_dbcsr_copy(blocked_matrix(ispin)%matrix, density_matrix(ispin)%matrix, 'BLOCKED MATRIX', error=error)

        CALL cp_dbcsr_copy(mat_tmp2(ispin)%matrix, density_matrix(ispin)%matrix, 'MAT_TMP2', error=error)
 
        ! zero all blocks
        CALL cp_dbcsr_set(blocked_matrix(ispin)%matrix, 0.0_dp, error=error)

        CALL cp_dbcsr_set(mat_tmp2(ispin)%matrix, 0.0_dp, error=error)

        CALL cp_dbcsr_reserve_blocks(mat_tmp, kg_env%molecular_blocks(imol)%rows, &
                                     kg_env%molecular_blocks(imol)%cols, error=error) 

        CALL cp_dbcsr_finalize(mat_tmp, error=error)
  
        CALL cp_dbcsr_copy_into_existing(mat_tmp, density_matrix(ispin)%matrix, error)
  
        CALL cp_dbcsr_copy_into_existing(blocked_matrix(ispin)%matrix, mat_tmp, error)

        CALL cp_dbcsr_copy(mat_ks_blocked(ispin)%matrix, mat_tmp, 'MAT_KS_BLOCKED', error=error)

        CALL cp_dbcsr_release(mat_tmp, error)
  
        DEALLOCATE(mat_tmp)
 
      END DO

      my_mat => rho_struct%rho_ao

      ! set the density matrix to the blocked matrix  
      rho_struct%rho_ao => blocked_matrix

      ! calculate the densities for the given blocked density matrix
      CALL qs_rho_update_rho(rho_struct, qs_env, error=error)

      ekin_imol=0.0_dp

      ! calc Hohenberg-Kohn kin. energy of the density corresp. to the remaining molecular block(s)
      CALL qs_vxc_create(qs_env=qs_env, rho_struct=rho_struct, xc_section=kg_env%xc_section_kg, &
                        vxc_rho=vxc_rho, vxc_tau=vxc_tau, exc=ekin_imol, error=error)

      rho_struct%rho_ao => my_mat

      ekin_mol = ekin_mol + ekin_imol

      DO ispin=1,nspins

        ! clean up vxc_rho
        CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw, error=error)

      END DO
      DEALLOCATE(vxc_rho)

      IF (ASSOCIATED(vxc_tau)) THEN
      ! this case has NOT been tested yet...
       
        DO ispin=1,nspins

          ! clean up vxc_tau
          CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_tau(ispin)%pw, error=error)
       
        END DO

        DEALLOCATE(vxc_tau)

      END IF 

      ! deallocate dbcsr matrix sets
      CALL cp_dbcsr_deallocate_matrix_set(blocked_matrix, error=error)    
 
      CALL cp_dbcsr_deallocate_matrix_set(mat_ks_blocked, error=error)  

      CALL cp_dbcsr_deallocate_matrix_set(mat_tmp2, error=error)    
    END DO

    ! clean up rho_struct
    CALL qs_rho_release(rho_struct, error) 

    CALL timestop(handle)

  END SUBROUTINE kg_ekin_mol

END MODULE kg_correction
