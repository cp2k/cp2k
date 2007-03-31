!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****** cp2k/qs_linres_epr_ownutils *
!!
!!   NAME
!!     qs_linres_epr_ownutils
!!
!!   FUNCTION
!!
!!   NOTE
!!
!!   AUTHOR
!!      RD
!!
!!   MODIFICATION HISTORY
!!      created 06-2006 [RD]
!!
!!   SOURCE
!!****

MODULE qs_linres_epr_ownutils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_type
  USE pw_types,                        ONLY: pw_integral_ab,&
                                             REALDATA3D,&
                                             REALSPACE
  USE qs_collocate_density,            ONLY: calculate_nablavks
  USE qs_core_energies,                ONLY: calculate_ecore
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_linres_epr_atom_current,      ONLY: calculate_jrho_atom_ang
  USE qs_linres_epr_op,                ONLY: set_vecp
  USE qs_linres_types,                 ONLY: epr_env_type,&
                                             get_epr_env,&
                                             jrho_atom_type,&
                                             nablavks_atom_type
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_atom_methods,             ONLY: calculate_rho_atom_coeff,&
                                             calculate_rho_atom
  USE qs_rho_types,                    ONLY: qs_rho_type,&
                                             qs_rho_p_type
  USE qs_vxc_atom,                     ONLY: calc_rho_angular
  USE realspace_grid_cube,             ONLY: rs_pw_to_cube
  USE sparse_matrix_types,             ONLY: allocate_matrix_set,&
                                             cp_sm_scale_and_add,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             transfer_matrix
  USE util,                            ONLY: get_limit

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: epr_g_print, epr_g_zke, epr_g_so, epr_nablavks, epr_g_soo

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_epr_ownutils'

CONTAINS

! *****************************************************************************
!!****** cp2k/qs_linres_epr_utils/epr_g_print
!!
!!   NAME
!!     epr_g_print
!!
!!   FUNCTION
!!     Prints the g tensor
!!
!!   NOTES
!!
!!   AUTHOR
!!     RD
!!
!!   MODIFICATION HISTORY
!!     06.2006 created [RD]
!!
!****************************************************************************

  SUBROUTINE epr_g_print(epr_env,qs_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_g_print', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: lr_section
    INTEGER                                  :: output_unit

    NULLIFY(logger, lr_section)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    IF (output_unit>0) THEN
       ! g_zke
       WRITE (UNIT=output_unit,FMT="(T2,A)") "g_zke"
       WRITE (UNIT=output_unit,FMT="(T18,F20.10,1X,F20.10,1X,F20.10)")&
              epr_env%g_zke,       0.0_dp,              0.0_dp,&
              0.0_dp,              epr_env%g_zke,       0.0_dp,&
              0.0_dp,              0.0_dp,              epr_env%g_zke
       ! g_so
       WRITE (UNIT=output_unit,FMT="(T2,A)") "g_so"
       WRITE (UNIT=output_unit,FMT="(T18,F20.10,1X,F20.10,1X,F20.10)")&
              epr_env%g_so(1,1),   epr_env%g_so(1,2),   epr_env%g_so(1,3),&
              epr_env%g_so(2,1),   epr_env%g_so(2,2),   epr_env%g_so(2,3),&
              epr_env%g_so(3,1),   epr_env%g_so(3,2),   epr_env%g_so(3,3)
       ! g_soo
       WRITE (UNIT=output_unit,FMT="(T2,A)") "g_soo"
       WRITE (UNIT=output_unit,FMT="(T18,F20.10,1X,F20.10,1X,F20.10)")&
              epr_env%g_soo(1,1),  epr_env%g_soo(1,2),  epr_env%g_soo(1,3),&
              epr_env%g_soo(2,1),  epr_env%g_soo(2,2),  epr_env%g_soo(2,3),&
              epr_env%g_soo(3,1),  epr_env%g_soo(3,2),  epr_env%g_soo(3,3)
       ! g_total
       WRITE (UNIT=output_unit,FMT="(T2,A)") "g_total"
       WRITE (UNIT=output_unit,FMT="(T18,F20.10,1X,F20.10,1X,F20.10)")&
              epr_env%g_total(1,1),epr_env%g_total(1,2),epr_env%g_total(1,3),&
              epr_env%g_total(2,1),epr_env%g_total(2,2),epr_env%g_total(2,3),&
              epr_env%g_total(3,1),epr_env%g_total(3,2),epr_env%g_total(3,3)
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_g_print

! *****************************************************************************
!!****** cp2k/qs_linres_epr_utils/epr_g_zke
!!
!!   NAME
!!     epr_g_zke
!!
!!   FUNCTION
!!     Calculate zke part of the g tensor
!!
!!   NOTES
!!
!!   AUTHOR
!!     RD
!!
!!   MODIFICATION HISTORY
!!     06.2006 created [RD]
!!
!****************************************************************************

  SUBROUTINE epr_g_zke(epr_env,qs_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_g_zke', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: lr_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: kinetic
    INTEGER                                  :: output_unit, ispin, i1
    REAL(KIND=dp)                            :: epr_g_zke_temp(2)

    NULLIFY(dft_control, logger, lr_section, rho, kinetic, para_env)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    kinetic=kinetic,rho=rho,para_env=para_env,error=error)

    DO ispin=1,dft_control%nspins
       CALL calculate_ecore(h=kinetic(1)%matrix,&
                            p=rho%rho_ao(ispin)%matrix,&
                            ecore=epr_g_zke_temp(ispin),&
                            para_env=para_env)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T60,F20.10)")&
          "epr(debug)|ZKE:Kinetic Energy of spin a/b ",epr_g_zke_temp(ispin)
       END IF
    END DO

    epr_env%g_zke = epr_env%g_zke_factor * ( epr_g_zke_temp(1) - epr_g_zke_temp(2) )
    DO i1 = 1,3
       epr_env%g_total(i1,i1) = epr_env%g_total(i1,i1) + epr_env%g_zke
    END DO

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(T2,A,T60,F20.10)")&
       "epr(debug)|ZKE:g_zke",epr_env%g_zke
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_g_zke

! *****************************************************************************
!!****** cp2k/qs_linres_epr_utils/epr_g_so
!!
!!   NAME
!!     epr_g_so
!!
!!   FUNCTION
!!     Calculates g_so
!!
!!   NOTES
!!
!!   AUTHOR
!!     RD
!!
!!   MODIFICATION HISTORY
!!     06.2006 created [RD]
!!
!****************************************************************************

  SUBROUTINE epr_g_so(epr_env,qs_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_so', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_unit, ispin, iB, idir1, idir2, idir3,&
                                                nkind, ikind, nr, na, natom, iat, iatom,&
                                                num_pe, mepos, ir, ia, nspins
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    REAL(dp)                                 :: temp_so_soft, hard_radius
    REAL(dp), DIMENSION(3,3)                 :: temp_so_gapw
    REAL(dp), DIMENSION(:,:), POINTER        :: g_total, g_so
    LOGICAL                                  :: gapw, paw_atom
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(jrho_atom_type),DIMENSION(:),&
      POINTER                                :: jrho1_atom_set
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    TYPE(qs_rho_p_type),DIMENSION(:,:),&
      POINTER                                :: jrho1_set
    TYPE(qs_rho_p_type),DIMENSION(:,:),&
      POINTER                                :: nablavks_set 
    TYPE(nablavks_atom_type),DIMENSION(:),&
      POINTER                                :: nablavks_atom_set 
    TYPE(nablavks_atom_type), POINTER        :: nablavks_atom
    TYPE(section_vals_type), POINTER         :: lr_section

    NULLIFY(logger, lr_section, para_env, dft_control)
    NULLIFY(jrho1_atom_set, jrho1_set, nablavks_set, nablavks_atom_set, nablavks_atom, g_total, g_so)
    NULLIFY(atomic_kind_set, grid_atom, harmonics, atom_kind, atom_list)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    para_env=para_env,error=error)

    CALL get_epr_env(epr_env=epr_env,&
                     nablavks_set=nablavks_set,jrho1_set=jrho1_set,&
                     nablavks_atom_set=nablavks_atom_set,jrho1_atom_set=jrho1_atom_set,&
                     g_total=g_total, g_so=g_so)

    gapw = dft_control%qs_control%gapw
    nkind = SIZE(atomic_kind_set,1)
    nspins = dft_control%nspins

    DO iB=1,3
       DO idir1=1,3
          CALL set_vecp(idir1,idir2,idir3)

          ! 1) j_glob x nabla_vks_glob
          temp_so_soft = 0.0_dp
          DO ispin=1,nspins
            temp_so_soft = temp_so_soft + (-1.0_dp)**(1+ispin) * ( &
                                pw_integral_ab ( jrho1_set(idir2,iB)%rho%rho_r(ispin)%pw, &
                                nablavks_set(idir3,ispin)%rho%rho_r(1)%pw ) - &
                                pw_integral_ab ( jrho1_set(idir3,iB)%rho%rho_r(ispin)%pw, &
                                nablavks_set(idir2,ispin)%rho%rho_r(1)%pw ) &
                                )
          END DO
          temp_so_soft = epr_env%g_so_factor * temp_so_soft
          IF (output_unit>0) THEN
               WRITE (UNIT=output_unit,FMT="(T2,A,T21,I1,I1,T60,F20.10)")&
               "epr(debug)|SO:soft",iB,idir1,temp_so_soft
          END IF
          g_so(iB,idir1) = temp_so_soft
       END DO !idir1
    END DO !iB

    IF (gapw) THEN
       DO ikind = 1,nkind
          NULLIFY (atom_kind,atom_list,grid_atom,harmonics)
          atom_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                               natom=natom, paw_atom=paw_atom,&
                               harmonics=harmonics,grid_atom=grid_atom,&
                               hard_radius=hard_radius)
          nr = grid_atom%nr
          na = grid_atom%ng_sphere
          IF (.NOT.paw_atom) CYCLE

          ! Distribute the atoms of this kind
    
          num_pe = para_env%num_pe
          mepos  = para_env%mepos
          bo = get_limit( natom, num_pe, mepos )

          DO iat = bo(1),bo(2)
             iatom = atom_list(iat)
             NULLIFY (jrho1_atom, nablavks_atom)
             jrho1_atom => jrho1_atom_set(iatom)
             nablavks_atom => nablavks_atom_set(iatom)
             DO iB = 1,3
                DO idir1 = 1,3
                   DO ispin = 1,nspins
                      jrho1_atom%jrho_vec_rad_h(idir1,iB,ispin)%r_coef = 0.0_dp
                      jrho1_atom%jrho_vec_rad_s(idir1,iB,ispin)%r_coef = 0.0_dp
                      CALL calculate_jrho_atom_ang(jrho1_atom,&
                           jrho1_atom%jrho_vec_rad_h(idir1,iB,ispin)%r_coef,&
                           jrho1_atom%jrho_vec_rad_s(idir1,iB,ispin)%r_coef,&
                           grid_atom,harmonics,iB,idir1,ispin,error=error)
                   END DO !ispin
                END DO

                DO idir1 = 1,3
                   temp_so_gapw(iB,idir1) = 0.0_dp 
                   CALL set_vecp(idir1,idir2,idir3)
                   DO ispin = 1,nspins
                      DO ia = 1,na
                         DO ir = 1,nr
                            IF (grid_atom%rad(ir) <= hard_radius) THEN
                               ! 2) sum_i j_loc_s_i x ( nabla_vks_loc_h_i - nabla_vks_loc_s_i )
                               temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + (-1.0_dp)**(1+ispin) * ( &
                                                   ( &
                                                   jrho1_atom%jrho_vec_rad_s(idir2,iB,ispin)%r_coef(ir,ia)*&
                                                   (nablavks_atom%nablavks_vec_rad_h(idir3,ispin)%r_coef(ir,ia)-&
                                                   nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia)) - &
                                                   jrho1_atom%jrho_vec_rad_s(idir3,iB,ispin)%r_coef(ir,ia)*&
                                                   (nablavks_atom%nablavks_vec_rad_h(idir2,ispin)%r_coef(ir,ia)-&
                                                   nablavks_atom%nablavks_vec_rad_s(idir2,ispin)%r_coef(ir,ia)) &
                                                   ) *grid_atom%weight(ia,ir)&
                                                   )
                               ! 3) sum_i ( j_loc_h_i - j_loc_s_i ) nabla_vks_loc_s_i
                               temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + (-1.0_dp)**(1+ispin) * ( &
                                                   ( &
                                                   (jrho1_atom%jrho_vec_rad_h(idir2,iB,ispin)%r_coef(ir,ia)-&
                                                   jrho1_atom%jrho_vec_rad_s(idir2,iB,ispin)%r_coef(ir,ia))*&
                                                   nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia) - &
                                                   (jrho1_atom%jrho_vec_rad_h(idir3,iB,ispin)%r_coef(ir,ia)-&
                                                   jrho1_atom%jrho_vec_rad_s(idir3,iB,ispin)%r_coef(ir,ia))*&
                                                   nablavks_atom%nablavks_vec_rad_s(idir2,ispin)%r_coef(ir,ia) &
                                                   ) *grid_atom%weight(ia,ir)&
                                                   )
                               ! 4) sum_i ( j_loc_h_i - j_loc_s_i ) ( nabla_vks_loc_h_i - nabla_vks_loc_s_i )
                               temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + (-1.0_dp)**(1+ispin) * ( &
                                                   ( &
                                                   (jrho1_atom%jrho_vec_rad_h(idir2,iB,ispin)%r_coef(ir,ia)-&
                                                   jrho1_atom%jrho_vec_rad_s(idir2,iB,ispin)%r_coef(ir,ia))*&
                                                   (nablavks_atom%nablavks_vec_rad_h(idir3,ispin)%r_coef(ir,ia)-&
                                                   nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia)) - &
                                                   (jrho1_atom%jrho_vec_rad_h(idir3,iB,ispin)%r_coef(ir,ia)-&
                                                   jrho1_atom%jrho_vec_rad_s(idir3,iB,ispin)%r_coef(ir,ia))*&
                                                   (nablavks_atom%nablavks_vec_rad_h(idir2,ispin)%r_coef(ir,ia)-&
                                                   nablavks_atom%nablavks_vec_rad_s(idir2,ispin)%r_coef(ir,ia)) &
                                                   ) *grid_atom%weight(ia,ir)&
                                                   )
                            END IF
                         END DO !ir
                      END DO !ia
                   END DO !ispin
                END DO !idir1
             END DO !iB
          END DO !iat
       END DO !ikind
       CALL mp_sum(temp_so_gapw,para_env%group)
       temp_so_gapw(:,:) = epr_env%g_so_factor * temp_so_gapw(:,:)
       IF (output_unit>0) THEN
          DO iB = 1,3
             DO idir1 = 1,3
                WRITE (UNIT=output_unit,FMT="(T2,A,T21,I1,I1,T60,F20.10)")&
                "epr(debug)|SO:gapw",iB,idir1,temp_so_gapw(iB,idir1)
             END DO
          END DO
       END IF
       g_so(:,:) = g_so(:,:) + temp_so_gapw(:,:)
    END IF
    g_total(:,:) = g_total(:,:) + g_so(:,:)

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_g_so

! *****************************************************************************
!!****** cp2k/qs_linres_epr_utils/epr_nablavks
!!
!!   NAME
!!     epr_nablavks
!!
!!   FUNCTION
!!     Calculates Nabla V_ks and collocates it on the real space global (GAPW
!!     and GPW) and atomic grids (GAPW only)
!!
!!   NOTES
!!
!!   AUTHOR
!!     RD
!!
!!   MODIFICATION HISTORY
!!     06.2006 created [RD]
!!
!****************************************************************************

  SUBROUTINE epr_nablavks(epr_env,qs_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_nablavks', &
      routineP = moduleN//':'//routineN

    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: lr_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: kinetic, matrix_ks
    INTEGER                                  :: output_unit, ispin, nspins, idir,&
                                                ikind, nkind, natom
    REAL(dp)                                 :: total_nablavks
    LOGICAL                                  :: gapw
    TYPE(coeff_type), POINTER                :: wf_r
    INTEGER                                  :: unit_nr, istat
    CHARACTER(LEN=80)                        :: ext, filename
    TYPE(section_vals_type), POINTER         :: epr_section
    LOGICAL                                  :: ionode
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_subsystem_type), POINTER         :: subsys

    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: vks_atom_set
    TYPE(rho_atom_type), POINTER             :: vks_atom
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    REAL(dp), DIMENSION(:), POINTER          :: rho1_h_tot, rho1_s_tot
    LOGICAL                                  :: paw_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(grid_atom_type), POINTER            :: grid_atom
    INTEGER                                  :: nr, na, ir, ia, num_pe, mepos, iat,&
                                                iatom, group
    REAL(dp), DIMENSION(:, :), POINTER       :: rho_h, rho_s
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: drho_h, drho_s
    INTEGER, DIMENSION(2)                    :: bo
    TYPE(nablavks_atom_type),DIMENSION(:),&
      POINTER                                :: nablavks_atom_set
    TYPE(nablavks_atom_type), POINTER        :: nablavks_atom

    NULLIFY(dft_control, logger, lr_section, kinetic, matrix_ks, para_env)
    NULLIFY(auxbas_pw_pool, epr_section)
    NULLIFY(particle_set, subsys, particles)
    NULLIFY(atomic_kind_set, atom_kind)
    NULLIFY(rho1_h_tot,rho1_s_tot)
    NULLIFY(atom_list)
    NULLIFY(vks_atom_set,vks_atom)
    NULLIFY(nablavks_atom_set,nablavks_atom)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)
    ionode = logger%para_env%mepos==logger%para_env%source

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    kinetic=kinetic,matrix_ks=matrix_ks,para_env=para_env,error=error,&
                    pw_env=pw_env,atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,subsys=subsys)
    CALL cp_subsys_get(subsys,particles=particles,error=error)

    ! V_ks matrix elements in matrix_v(ispin) for each spin channel
    ! and calculate Nabla V_ks on the global grid
    nspins = dft_control%nspins
    CALL allocate_matrix_set(epr_env%matrix_v,nspins,error=error)
    DO ispin=1,nspins
       CALL replicate_matrix_structure(source=matrix_ks(ispin)%matrix,&
                                       TARGET=epr_env%matrix_v(ispin)%matrix,&
                                       target_name="KOHN SHAM POTENTIAL MATRIX",&
                                       error=error)
       CALL transfer_matrix(matrix_ks(ispin)%matrix,epr_env%matrix_v(ispin)%matrix,&
               error=error)
       CALL cp_sm_scale_and_add(epr_env%matrix_v(ispin)%matrix,1.0_dp,&
               kinetic(1)%matrix,-1.0_dp,error=error)
       DO idir=1,3
          CALL calculate_nablavks(matrix_p=epr_env%matrix_v(ispin)%matrix,&
               rho=epr_env%nablavks_set(idir,ispin)%rho%rho_r(1),&
               rho_gspace=epr_env%nablavks_set(idir,ispin)%rho%rho_g(1),&
               total_rho=total_nablavks,&
               idir=idir,&
               qs_env=qs_env,soft_valid=dft_control%qs_control%gapw,error=error)
          IF (output_unit>0) THEN
               WRITE (UNIT=output_unit,FMT="(T2,A,T60,F20.10)")&
               "epr(debug)|SO:total_nablavks(idir,ispin)",total_nablavks
          END IF
       END DO
    END DO

    ! Cubes

    epr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES%EPR",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,epr_section,&
                  "PRINT%NABLAVKS_CUBES",error=error),cp_p_file)) THEN
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
       ALLOCATE(wf_r,STAT=istat)
       CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
       DO idir = 1,3
          CALL coeff_zero(wf_r)
          CALL coeff_copy(epr_env%nablavks_set(idir,1)%rho%rho_r(1),wf_r)
          filename="nablavks"
          WRITE(ext,'(a2,I1,a5)')  "_d",idir,".cube"
          unit_nr=cp_print_key_unit_nr(logger,epr_section,"PRINT%NABLAVKS_CUBES",&
                  extension=TRIM(ext),middle_name=TRIM(filename),&
                  log_filename=.FALSE.,file_position="REWIND",error=error)
          CALL rs_pw_to_cube(wf_r%pw,unit_nr,ionode,"NABLA V_KS ",&
                  particles=particles,&
                  stride=section_get_ival(epr_section,"PRINT%NABLAVKS_CUBES%STRIDE",error=error),&
                  error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,epr_section,&
                  "PRINT%NABLAVKS_CUBES",error=error)
      END DO
      DEALLOCATE(wf_r,STAT=istat)
    END IF

    ! gapw part

    gapw = dft_control%qs_control%gapw
    IF(gapw) THEN

       group = para_env%group

       CALL get_epr_env(epr_env=epr_env,&
                        vks_atom_set=vks_atom_set)

       ! Nabla_V_KS on the local grids
       CALL calculate_rho_atom_coeff(qs_env, epr_env%matrix_v, vks_atom_set, error=error)
        
       nkind  = SIZE(atomic_kind_set,1)

       ALLOCATE(rho1_h_tot(1:nspins), rho1_s_tot(1:nspins))
       rho1_h_tot = 0.0_dp
       rho1_s_tot = 0.0_dp

       DO ikind = 1,nkind
          atom_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atom_kind, atom_list=atom_list,&
                               grid_atom=grid_atom,natom=natom,&
                               paw_atom=paw_atom, harmonics=harmonics)
    
    !     Calculate rho1_h and rho1_s on the radial gris centered on the atomic position
          IF(paw_atom) &
            CALL calculate_rho_atom(para_env,vks_atom_set,atom_kind,atom_list,grid_atom,natom,&
                                    nspins,rho1_h_tot,rho1_s_tot,error=error)
       ENDDO
   
       CALL mp_sum(rho1_h_tot,group)
       CALL mp_sum(rho1_s_tot,group)
       IF (output_unit>0) THEN
          DO ispin = 1,nspins
             WRITE (UNIT=output_unit,FMT="(T2,A,T60,F20.10)")&
             "epr(debug)|SO:total_vks_hard(ispin)",rho1_h_tot(ispin)
             WRITE (UNIT=output_unit,FMT="(T2,A,T60,F20.10)")&
             "epr(debug)|SO:total_vks_soft(ispin)",rho1_s_tot(ispin)
          END DO
       END IF

       DEALLOCATE(rho1_h_tot,rho1_s_tot)

       CALL get_epr_env(epr_env=epr_env,&
                        nablavks_atom_set=nablavks_atom_set)

       ! Nullify some pointers for work-arrays

       NULLIFY (rho_h,drho_h,rho_s,drho_s)

       ! Here starts the loop over all the atoms
       DO ikind = 1,nkind
 
          atom_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                               natom=natom,paw_atom=paw_atom,&
                               harmonics=harmonics,grid_atom=grid_atom)
    
          IF (.NOT.paw_atom) CYCLE
    
          nr = grid_atom%nr
          na = grid_atom%ng_sphere
           
          CALL reallocate(rho_h,1,na,1,nspins)
          CALL reallocate(rho_s,1,na,1,nspins)
          CALL reallocate(drho_h,1,4,1,na,1,nr,1,nspins)
          CALL reallocate(drho_s,1,4,1,na,1,nr,1,nspins)

          ! Distribute the atoms of this kind
    
          num_pe = para_env%num_pe
          mepos  = para_env%mepos
          bo = get_limit( natom, num_pe, mepos )
    
          DO iat = bo(1),bo(2) !1,natom
            iatom = atom_list(iat)
            vks_atom => vks_atom_set(iatom)
            nablavks_atom => nablavks_atom_set(iatom)

            DO ir = 1,nr
              CALL calc_rho_angular(grid_atom, harmonics, nspins, .TRUE.,&
                                    ir, vks_atom%rho_rad_h, vks_atom%rho_rad_s,&
                                    rho_h, rho_s, vks_atom%drho_rad_h, vks_atom%drho_rad_s, &
                                    vks_atom%rho_rad_h_d, vks_atom%rho_rad_s_d,&
                                    drho_h, drho_s, error=error)
            END DO
            DO idir = 1,3
               DO ir = 1,nr
                  DO ia = 1,na
                     DO ispin = 1,nspins
                        nablavks_atom%nablavks_vec_rad_h(idir,ispin)%r_coef(ir,ia) = &
                                      drho_h(idir,ia,ir,ispin)
                        nablavks_atom%nablavks_vec_rad_s(idir,ispin)%r_coef(ir,ia) = &
                                      drho_s(idir,ia,ir,ispin)
                     END DO
                  END DO
               END DO
            END DO
    
          END DO ! iat
    
       END DO ! ikind
    
       CALL mp_sync(group)

    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_nablavks

! *****************************************************************************
!!****** cp2k/qs_linres_epr_utils/epr_g_soo
!!
!!   NAME
!!     epr_g_soo
!!
!!   FUNCTION
!!     Calculates g_soo (soft part only for now)
!!
!!   NOTES
!!
!!   AUTHOR
!!     RD
!!
!!   MODIFICATION HISTORY
!!     06.2006 created [RD]
!!
!****************************************************************************

  SUBROUTINE epr_g_soo(epr_env,qs_env,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_soo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_unit, ispin, iB, idir1, &
                                                nspins
    REAL(dp), DIMENSION(:,:), POINTER        :: g_total, g_soo, bind0
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: lr_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_rho_p_type), DIMENSION(:,:), &
      POINTER                                :: bind_set

    NULLIFY(logger, lr_section, para_env, dft_control, rho,&
            bind_set, g_soo, g_total, bind0)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    rho=rho,para_env=para_env,&
                    error=error)

    CALL get_epr_env(epr_env=epr_env,bind_set=bind_set,&
                     g_soo=g_soo,g_total=g_total,bind0=bind0)

    DO iB=1,3
       DO idir1=1,3
          DO ispin=1,dft_control%nspins
             g_soo(iB,idir1) = g_soo(iB,idir1) + (-1.0_dp)**(1+ispin) * ( &
                               pw_integral_ab ( bind_set(idir1,iB)%rho%rho_r(1)%pw, &
                               rho%rho_r(ispin)%pw ) &
                               )
          END DO
          g_soo(iB,idir1) = epr_env%g_soo_factor * g_soo(iB,idir1)
          IF (output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(T2,A,T31,i1,i1,T60,F20.10)")&
             "epr(debug)|SOO:soft",iB,idir1,g_soo(idir1,iB)
          END IF
       END DO
    END DO

    DO iB=1,3
       DO idir1=1,3
          bind0(idir1,iB) = 2.0_dp / 3.0_dp * epr_env%g_soo_factor *&
                            bind0(idir1,iB) * ( 1.0_dp / 137.03602_dp )
          g_soo(iB,idir1) = g_soo(iB,idir1) + bind0(idir1,iB)
          g_total(iB,idir1) = g_total(iB,idir1) + g_soo(iB,idir1)
          IF (output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(T2,A,T34,i1,i1,T60,F20.10)")&
             "epr(debug)|SOO:soft_g0_analytic",iB,idir1,epr_env%bind0(idir1,iB)
          END IF
       END DO
    END DO

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_g_soo

END MODULE  qs_linres_epr_ownutils

