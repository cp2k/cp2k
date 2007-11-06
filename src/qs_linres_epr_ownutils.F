!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_integral_ab
  USE qs_core_energies,                ONLY: calculate_ecore
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_linres_epr_op,                ONLY: set_vecp
  USE qs_linres_types,                 ONLY: epr_env_type,&
                                             get_epr_env,&
                                             jrho_atom_type,&
                                             nablavks_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_p_type,&
                                             qs_rho_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: epr_g_print, epr_g_zke, epr_g_so, epr_g_soo

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

    INTEGER                                  :: output_unit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: lr_section

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

    INTEGER                                  :: i1, ispin, output_unit
    REAL(KIND=dp)                            :: epr_g_zke_temp(2)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: kinetic
    TYPE(section_vals_type), POINTER         :: lr_section

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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_g_so', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iB, idir1, idir2, idir3, &
                                                ispin, nkind, nspins, &
                                                output_unit
    LOGICAL                                  :: gapw
    REAL(dp)                                 :: temp_so_soft
    REAL(dp), DIMENSION(:, :), POINTER       :: g_so, g_total
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      POINTER                                :: jrho1_set, nablavks_set
    TYPE(section_vals_type), POINTER         :: lr_section

    NULLIFY(logger, lr_section, para_env, dft_control)
    NULLIFY(jrho1_atom_set, jrho1_set, nablavks_set, nablavks_atom_set, &
         g_total, g_so)
    NULLIFY(atomic_kind_set)

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
            temp_so_soft = temp_so_soft + (-1.0_dp)**(1.0_dp + REAL(ispin,KIND=dp)) * ( &
                                pw_integral_ab ( jrho1_set(idir2,iB)%rho%rho_r(ispin)%pw, &
                                nablavks_set(idir3,ispin)%rho%rho_r(1)%pw ) - &
                                pw_integral_ab ( jrho1_set(idir3,iB)%rho%rho_r(ispin)%pw, &
                                nablavks_set(idir2,ispin)%rho%rho_r(1)%pw ) &
                                )
          END DO
          temp_so_soft = epr_env%g_so_factor * temp_so_soft
          IF (output_unit>0) THEN
               WRITE (UNIT=output_unit,FMT="(T2,A,T25,I1,I1,T60,F20.10)")&
               "epr(debug)|SO:soft",iB,idir1,temp_so_soft
          END IF
          g_so(iB,idir1) = temp_so_soft
       END DO !idir1
    END DO !iB

!R    IF (gapw) THEN
!R       DO ikind = 1,nkind
!R          NULLIFY (atom_kind,atom_list,grid_atom,harmonics)
!R          atom_kind => atomic_kind_set(ikind)
!R          CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
!R                               natom=natom, paw_atom=paw_atom,&
!R                               harmonics=harmonics,grid_atom=grid_atom,&
!R                               hard_radius=hard_radius)
!R          nr = grid_atom%nr
!R          na = grid_atom%ng_sphere
!R          IF (.NOT.paw_atom) CYCLE
!R
!R          ! Distribute the atoms of this kind
!R    
!R          num_pe = para_env%num_pe
!R          mepos  = para_env%mepos
!R          bo = get_limit( natom, num_pe, mepos )
!R
!R          DO iat = bo(1),bo(2)
!R             iatom = atom_list(iat)
!R             NULLIFY (jrho1_atom, nablavks_atom)
!R             jrho1_atom => jrho1_atom_set(iatom)
!R             nablavks_atom => nablavks_atom_set(iatom)
!R             DO iB = 1,3
!R                DO idir1 = 1,3
!R                   DO ispin = 1,nspins
!R                      jrho1_atom%jrho_vec_rad_h(idir1,iB,ispin)%r_coef = 0.0_dp
!R                      jrho1_atom%jrho_vec_rad_s(idir1,iB,ispin)%r_coef = 0.0_dp
!R                      CALL calculate_jrho_atom_ang(jrho1_atom,&
!R                           jrho1_atom%jrho_vec_rad_h(idir1,iB,ispin)%r_coef,&
!R                           jrho1_atom%jrho_vec_rad_s(idir1,iB,ispin)%r_coef,&
!R                           grid_atom,harmonics,iB,idir1,ispin,error=error)
!R                   END DO !ispin
!R                END DO
!R
!R                DO idir1 = 1,3
!R                   temp_so_gapw(iB,idir1) = 0.0_dp 
!R                   CALL set_vecp(idir1,idir2,idir3)
!R                   DO ispin = 1,nspins
!R                      DO ia = 1,na
!R                         DO ir = 1,nr
!R                            IF (grid_atom%rad(ir) <= hard_radius) THEN
!R                               ! 2) sum_i j_loc_s_i x ( nabla_vks_loc_h_i - nabla_vks_loc_s_i )
!R                               temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + (-1.0_dp)**(1+ispin) * ( &
!R                                                   ( &
!R                                                   jrho1_atom%jrho_vec_rad_s(idir2,iB,ispin)%r_coef(ir,ia)*&
!R                                                   (nablavks_atom%nablavks_vec_rad_h(idir3,ispin)%r_coef(ir,ia)-&
!R                                                   nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia)) - &
!R                                                   jrho1_atom%jrho_vec_rad_s(idir3,iB,ispin)%r_coef(ir,ia)*&
!R                                                   (nablavks_atom%nablavks_vec_rad_h(idir2,ispin)%r_coef(ir,ia)-&
!R                                                   nablavks_atom%nablavks_vec_rad_s(idir2,ispin)%r_coef(ir,ia)) &
!R                                                   ) *grid_atom%weight(ia,ir)&
!R                                                   )
!R                               ! 3) sum_i ( j_loc_h_i - j_loc_s_i ) nabla_vks_loc_s_i
!R                               temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + (-1.0_dp)**(1+ispin) * ( &
!R                                                   ( &
!R                                                   (jrho1_atom%jrho_vec_rad_h(idir2,iB,ispin)%r_coef(ir,ia)-&
!R                                                   jrho1_atom%jrho_vec_rad_s(idir2,iB,ispin)%r_coef(ir,ia))*&
!R                                                   nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia) - &
!R                                                   (jrho1_atom%jrho_vec_rad_h(idir3,iB,ispin)%r_coef(ir,ia)-&
!R                                                   jrho1_atom%jrho_vec_rad_s(idir3,iB,ispin)%r_coef(ir,ia))*&
!R                                                   nablavks_atom%nablavks_vec_rad_s(idir2,ispin)%r_coef(ir,ia) &
!R                                                   ) *grid_atom%weight(ia,ir)&
!R                                                   )
!R                               ! 4) sum_i ( j_loc_h_i - j_loc_s_i ) ( nabla_vks_loc_h_i - nabla_vks_loc_s_i )
!R                               temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + (-1.0_dp)**(1+ispin) * ( &
!R                                                   ( &
!R                                                   (jrho1_atom%jrho_vec_rad_h(idir2,iB,ispin)%r_coef(ir,ia)-&
!R                                                   jrho1_atom%jrho_vec_rad_s(idir2,iB,ispin)%r_coef(ir,ia))*&
!R                                                   (nablavks_atom%nablavks_vec_rad_h(idir3,ispin)%r_coef(ir,ia)-&
!R                                                   nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia)) - &
!R                                                   (jrho1_atom%jrho_vec_rad_h(idir3,iB,ispin)%r_coef(ir,ia)-&
!R                                                   jrho1_atom%jrho_vec_rad_s(idir3,iB,ispin)%r_coef(ir,ia))*&
!R                                                   (nablavks_atom%nablavks_vec_rad_h(idir2,ispin)%r_coef(ir,ia)-&
!R                                                   nablavks_atom%nablavks_vec_rad_s(idir2,ispin)%r_coef(ir,ia)) &
!R                                                   ) *grid_atom%weight(ia,ir)&
!R                                                   )
!R                            END IF
!R                         END DO !ir
!R                      END DO !ia
!R                   END DO !ispin
!R                END DO !idir1
!R             END DO !iB
!R          END DO !iat
!R       END DO !ikind
!R       CALL mp_sum(temp_so_gapw,para_env%group)
!R       temp_so_gapw(:,:) = epr_env%g_so_factor * temp_so_gapw(:,:)
!R       IF (output_unit>0) THEN
!R          DO iB = 1,3
!R             DO idir1 = 1,3
!R                WRITE (UNIT=output_unit,FMT="(T2,A,T21,I1,I1,T60,F20.10)")&
!R                "epr(debug)|SO:gapw",iB,idir1,temp_so_gapw(iB,idir1)
!R             END DO
!R          END DO
!R       END IF
!R       g_so(:,:) = g_so(:,:) + temp_so_gapw(:,:)
!R    END IF
    g_total(:,:) = g_total(:,:) + g_so(:,:)

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_g_so

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

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_g_soo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iB, idir1, ispin, nspins, &
                                                output_unit
    REAL(dp), DIMENSION(:, :), POINTER       :: bind0, g_soo, g_total
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      POINTER                                :: bind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: lr_section

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
             g_soo(iB,idir1) = g_soo(iB,idir1) + (-1.0_dp)**(1.0_dp + REAL(ispin,KIND=dp)) * ( &
                               pw_integral_ab ( bind_set(idir1,iB)%rho%rho_r(1)%pw, &
                               rho%rho_r(ispin)%pw ) &
                               )
          END DO
          g_soo(iB,idir1) = epr_env%g_soo_factor * g_soo(iB,idir1)
          IF (output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(T2,A,T25,i1,i1,T60,F20.10)")&
             "epr(debug)|SOO:soft",iB,idir1,g_soo(idir1,iB)
          END IF
       END DO
    END DO

    DO iB=1,3
       DO idir1=1,3
       WRITE(6,*) iB, idir1, bind0(idir1,iB)
          g_soo(iB,idir1) = g_soo(iB,idir1) + epr_env%g_soo_chicorr_factor * bind0(idir1,iB)
          g_total(iB,idir1) = g_total(iB,idir1) + g_soo(iB,idir1)
          IF (output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(T2,A,T25,i1,i1,T60,F20.10)")&
             "epr(debug)|SOO:soft_g0",iB,idir1,epr_env%g_soo_chicorr_factor * bind0(idir1,iB)
          END IF
       END DO
    END DO

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_g_soo

END MODULE  qs_linres_epr_ownutils

