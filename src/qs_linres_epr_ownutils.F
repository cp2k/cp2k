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
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_val_get,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: twopi
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_pw,&
                                             pw_pool_create_pw,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc,&
                                             find_coeffs,&
                                             pw_spline_do_precond,&
                                             pw_spline_precond_create,&
                                             pw_spline_precond_release,&
                                             pw_spline_precond_set_kind,&
                                             pw_spline_precond_type,&
                                             spl3_pbc
  USE pw_types,                        ONLY: pw_integral_ab,&
                                             pw_p_type,&
                                             REALDATA3D,&
                                             REALSPACE
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
  USE qs_rho_types,                    ONLY: qs_rho_p_type,&
                                             qs_rho_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE util,                            ONLY: get_limit

#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: epr_g_print, epr_g_zke, epr_g_so, epr_g_soo, calculate_jrho_atom

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

    CHARACTER(LEN=default_string_length)     :: title
    INTEGER                                  :: idir1,idir2,output_unit,&
                                                unit_nr
    REAL(KIND=dp)                            :: eigenv_g(3),g_sym(3,3),&
                                                gsum
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: lr_section

    NULLIFY(logger, lr_section)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    gsum = 0.0_dp

    DO idir1 = 1,3
       DO idir2 = 1,3
          gsum = gsum + epr_env%g_total(idir1,idir2)
       END DO
    END DO

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(T2,A,T56,E24.16)")&
       "epr|TOT:checksum",gsum
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,lr_section,&
              "EPR%PRINT%G_TENSOR",error=error),cp_p_file)) THEN

       unit_nr=cp_print_key_unit_nr(logger,lr_section,"EPR%PRINT%G_TENSOR",&
            &                          extension=".data",middle_name="GTENSOR",&
            &                          log_filename=.FALSE.,error=error)

       IF(unit_nr > 0) THEN

          WRITE(title,"(A)") "G tensor "
          WRITE(unit_nr,"(T2,A)") title
   
          WRITE(unit_nr,"(T2,A)") "gmatrix_zke"
          WRITE(unit_nr,"(3(A,f15.10))")" XX=",epr_env%g_zke,&
               " XY=",0.0_dp," XZ=",0.0_dp
          WRITE(unit_nr,"(3(A,f15.10))")" YX=",0.0_dp,&
               " YY=",epr_env%g_zke," YZ=",0.0_dp
          WRITE(unit_nr,"(3(A,f15.10))")" ZX=",0.0_dp,&
               " ZY=",0.0_dp," ZZ=",epr_env%g_zke
   
          WRITE(unit_nr,"(T2,A)") "gmatrix_so"
          WRITE(unit_nr,"(3(A,f15.10))")" XX=",epr_env%g_so(1,1),&
               " XY=",epr_env%g_so(1,2)," XZ=",epr_env%g_so(1,3)
          WRITE(unit_nr,"(3(A,f15.10))")" YX=",epr_env%g_so(2,1),&
               " YY=",epr_env%g_so(2,2)," YZ=",epr_env%g_so(2,3)
          WRITE(unit_nr,"(3(A,f15.10))")" ZX=",epr_env%g_so(3,1),&
               " ZY=",epr_env%g_so(3,2)," ZZ=",epr_env%g_so(3,3)
   
          WRITE(unit_nr,"(T2,A)") "gmatrix_soo"
          WRITE(unit_nr,"(3(A,f15.10))")" XX=",epr_env%g_soo(1,1),&
               " XY=",epr_env%g_soo(1,2)," XZ=",epr_env%g_soo(1,3)
          WRITE(unit_nr,"(3(A,f15.10))")" YX=",epr_env%g_soo(2,1),&
               " YY=",epr_env%g_soo(2,2)," YZ=",epr_env%g_soo(2,3)
          WRITE(unit_nr,"(3(A,f15.10))")" ZX=",epr_env%g_soo(3,1),&
               " ZY=",epr_env%g_soo(3,2)," ZZ=",epr_env%g_soo(3,3)
   
          WRITE(unit_nr,"(T2,A)") "gmatrix_total"
          WRITE(unit_nr,"(3(A,f15.10))")" XX=",epr_env%g_total(1,1)+epr_env%g_free_factor,&
               " XY=",epr_env%g_total(1,2)," XZ=",epr_env%g_total(1,3)
          WRITE(unit_nr,"(3(A,f15.10))")" YX=",epr_env%g_total(2,1),&
               " YY=",epr_env%g_total(2,2)+epr_env%g_free_factor," YZ=",epr_env%g_total(2,3)
          WRITE(unit_nr,"(3(A,f15.10))")" ZX=",epr_env%g_total(3,1),&
               " ZY=",epr_env%g_total(3,2)," ZZ=",epr_env%g_total(3,3)+epr_env%g_free_factor
   
          DO idir1 = 1,3
             DO idir2 = 1,3
                g_sym(idir1,idir2) = (epr_env%g_total(idir1,idir2) + &
                                      epr_env%g_total(idir2,idir1))/2.0_dp
             END DO
          END DO
   
          WRITE(unit_nr,"(T2,A)") "gtensor_total"
          WRITE(unit_nr,"(3(A,f15.10))")" XX=",g_sym(1,1)+epr_env%g_free_factor,&
               " XY=",g_sym(1,2)," XZ=",g_sym(1,3)
          WRITE(unit_nr,"(3(A,f15.10))")" YX=",g_sym(2,1),&
               " YY=",g_sym(2,2)+epr_env%g_free_factor," YZ=",g_sym(2,3)
          WRITE(unit_nr,"(3(A,f15.10))")" ZX=",g_sym(3,1),&
               " ZY=",g_sym(3,2)," ZZ=",g_sym(3,3)+epr_env%g_free_factor
   
          CALL diamat_all(g_sym,eigenv_g)
          eigenv_g(:) = eigenv_g(:)*1.0e6_dp
   
          WRITE(unit_nr,"(T2,A)") "delta_g principal values in ppm"
          WRITE(unit_nr,"(f15.3,3(A,f15.10))") eigenv_g(1)," X=",g_sym(1,1),&
               " Y=",g_sym(2,1)," Z=",g_sym(3,1)
          WRITE(unit_nr,"(f15.3,3(A,f15.10))") eigenv_g(2)," X=",g_sym(1,2),&
               " Y=",g_sym(2,2)," Z=",g_sym(3,2)
          WRITE(unit_nr,"(f15.3,3(A,f15.10))") eigenv_g(3)," X=",g_sym(1,3),&
               " Y=",g_sym(2,3)," Z=",g_sym(3,3)

       END IF

       CALL cp_print_key_finished_output(unit_nr,logger,lr_section,&
            &                            "EPR%PRINT%G_TENSOR",error=error)

    ENDIF

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

    INTEGER                                  :: output_unit, ispin, i1
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
    END DO

    epr_env%g_zke = epr_env%g_zke_factor * ( epr_g_zke_temp(1) - epr_g_zke_temp(2) )
    DO i1 = 1,3
       epr_env%g_total(i1,i1) = epr_env%g_total(i1,i1) + epr_env%g_zke
    END DO

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(T2,A,T56,E24.16)")&
       "epr|ZKE:g_zke",epr_env%g_zke
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_g_zke

! *****************************************************************************
!!****** cp2k/qs_linres_epr_utils/calculate_jrho_atom
!!
!!   NAME
!!     calculate_jrho_atom
!!
!!   FUNCTION
!!     Calculates atomic currents
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

  SUBROUTINE calculate_jrho_atom(epr_env,qs_env,iB,idir,error)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    INTEGER, INTENT(IN)                      :: iB, idir

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_jrho_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: output_unit, ispin,&
                                                nkind, ikind, natom, iat, iatom,&
                                                nspins
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
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
    TYPE(section_vals_type), POINTER         :: lr_section

    NULLIFY(logger, lr_section, para_env, dft_control)
    NULLIFY(jrho1_atom_set, grid_atom, harmonics)
    NULLIFY(atomic_kind_set, atom_kind, atom_list)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    para_env=para_env,error=error) 

    CALL get_epr_env(epr_env=epr_env,&
                     jrho1_atom_set=jrho1_atom_set)

    gapw = dft_control%qs_control%gapw
    nkind = SIZE(atomic_kind_set,1)
    nspins = dft_control%nspins

    IF (gapw) THEN
       DO ikind = 1,nkind
          NULLIFY (atom_kind,atom_list,grid_atom,harmonics)
          atom_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                               harmonics=harmonics,grid_atom=grid_atom,&
                               natom=natom, paw_atom=paw_atom)
          IF (.NOT.paw_atom) CYCLE

          ! Distribute the atoms of this kind
    
          bo = get_limit( natom, para_env%num_pe, para_env%mepos )

          DO iat = bo(1),bo(2)
             iatom = atom_list(iat)
             NULLIFY (jrho1_atom)
             jrho1_atom => jrho1_atom_set(iatom)
             DO ispin = 1,nspins
                jrho1_atom%jrho_vec_rad_h(idir,iB,ispin)%r_coef = 0.0_dp
                jrho1_atom%jrho_vec_rad_s(idir,iB,ispin)%r_coef = 0.0_dp
                CALL calculate_jrho_atom_ang(jrho1_atom,&
                     jrho1_atom%jrho_vec_rad_h(idir,iB,ispin)%r_coef,&
                     jrho1_atom%jrho_vec_rad_s(idir,iB,ispin)%r_coef,&
                     grid_atom,harmonics,iB,idir,ispin,error=error)
             END DO !ispin
          END DO !iat
       END DO !ikind
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE calculate_jrho_atom 

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

    INTEGER                                  :: &
      aint_precond, ia, iat, iatom, iB, idir1, idir2, idir3,&
      ikind, ispin, ir, istat, max_iter, precond_kind, natom,&
      nkind, nspins, output_unit
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: gapw, paw_atom, failure,&
                                                success
    REAL(dp)                                 :: eps_r, eps_x,&
                                                temp_so_soft,&
                                                vks_ra_idir2,&
                                                vks_ra_idir3
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(dp), DIMENSION(3,3)                 :: temp_so_gapw
    REAL(dp), DIMENSION(:,:), POINTER        :: g_total, g_so
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(pw_p_type), DIMENSION(:,:),&
      POINTER                                :: vks_pw_spline
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(jrho_atom_type),DIMENSION(:),&
      POINTER                                :: jrho1_atom_set
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    TYPE(nablavks_atom_type),DIMENSION(:),&
      POINTER                                :: nablavks_atom_set 
    TYPE(nablavks_atom_type), POINTER        :: nablavks_atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(qs_rho_p_type),DIMENSION(:,:),&
      POINTER                                :: jrho1_set, nablavks_set
    TYPE(section_vals_type), POINTER         :: interp_section,&
                                                lr_section

    NULLIFY(atomic_kind_set, atom_kind, atom_list, dft_control,&
      grid_atom, g_so, g_total, harmonics, interp_section, jrho1_atom_set,&
      jrho1_set, logger, lr_section, nablavks_atom, nablavks_atom_set,&
      nablavks_set, para_env, particle_set)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    para_env=para_env,pw_env=pw_env,&
                    particle_set=particle_set,&
                    error=error)

    CALL get_epr_env(epr_env=epr_env,&
                     nablavks_set=nablavks_set,jrho1_set=jrho1_set,&
                     nablavks_atom_set=nablavks_atom_set,jrho1_atom_set=jrho1_atom_set,&
                     g_total=g_total,g_so=g_so)

    gapw = dft_control%qs_control%gapw
    nkind = SIZE(atomic_kind_set,1)
    nspins = dft_control%nspins

    DO iB=1,3
       DO idir1=1,3
          CALL set_vecp(idir1,idir2,idir3)
          ! j_pw x nabla_vks_pw
          temp_so_soft = 0.0_dp
          DO ispin=1,nspins
            temp_so_soft = temp_so_soft + (-1.0_dp)**(1 + ispin) * ( &
                                pw_integral_ab ( jrho1_set(idir2,iB)%rho%rho_r(ispin)%pw, &
                                nablavks_set(idir3,ispin)%rho%rho_r(1)%pw ) - &
                                pw_integral_ab ( jrho1_set(idir3,iB)%rho%rho_r(ispin)%pw, &
                                nablavks_set(idir2,ispin)%rho%rho_r(1)%pw ) &
                                )
          END DO
          temp_so_soft = -1.0_dp * epr_env%g_so_factor * temp_so_soft
          IF (output_unit>0) THEN
               WRITE (UNIT=output_unit,FMT="(T2,A,T18,I1,I1,T56,E24.16)")&
               "epr|SOX:soft",iB,idir1,temp_so_soft
          END IF
          g_so(iB,idir1) = temp_so_soft
       END DO !idir1
    END DO !iB

    IF (gapw) THEN

       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
       ALLOCATE(vks_pw_spline(3,nspins),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       interp_section =>  section_vals_get_subs_vals(lr_section,&
                           "EPR%INTERPOLATOR",error=error)
       CALL section_vals_val_get(interp_section,"aint_precond", &
                           i_val=aint_precond, error=error)
       CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind,error=error)
       CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter,error=error)
       CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r,error=error)
       CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x,error=error)

       DO ispin = 1,nspins
          DO idir1 = 1,3
             CALL pw_pool_create_pw(auxbas_pw_pool,&
                                    vks_pw_spline(idir1,ispin)%pw,&
                                    use_data=REALDATA3D,in_space=REALSPACE,&
                                    error=error)
             ! calculate spline coefficients
             CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
                     pool=auxbas_pw_pool,pbc=.TRUE.,transpose=.FALSE.,error=error)
             CALL pw_spline_do_precond(precond,&
                     nablavks_set(idir1,ispin)%rho%rho_r(1)%pw,&
                     vks_pw_spline(idir1,ispin)%pw,error=error)
             CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
             success=find_coeffs(values=nablavks_set(idir1,ispin)%rho%rho_r(1)%pw,&
                     coeffs=vks_pw_spline(idir1,ispin)%pw,linOp=spl3_pbc,&
                     preconditioner=precond,pool=auxbas_pw_pool,&
                     eps_r=eps_r,eps_x=eps_x,max_iter=max_iter,&
                     print_section=interp_section,error=error)
             CPPostconditionNoFail(success,cp_warning_level,routineP,error)
             CALL pw_spline_precond_release(precond,error=error)
          END DO ! idir1
       END DO ! ispin

       temp_so_gapw = 0.0_dp

       DO ikind = 1, nkind
          NULLIFY (atom_kind,atom_list,grid_atom,harmonics)
          atom_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                               grid_atom=grid_atom,harmonics=harmonics,&
                               natom=natom, paw_atom=paw_atom)

          IF (.NOT.paw_atom) CYCLE

          ! Distribute the atoms of this kind
    
          bo = get_limit( natom, para_env%num_pe, para_env%mepos )

          DO iat = 1,natom !bo(1),bo(2)! this partitioning blocks the interpolation 
             !                         ! routines (i.e. waiting for mp_sum)
             iatom = atom_list(iat)
             NULLIFY (jrho1_atom, nablavks_atom)
             jrho1_atom => jrho1_atom_set(iatom)
             nablavks_atom => nablavks_atom_set(iatom)
             DO iB = 1,3
                DO idir1 = 1,3
                   CALL set_vecp(idir1,idir2,idir3)
                   DO ispin = 1,nspins
                      DO ir = 1,grid_atom%nr

                         IF (grid_atom%rad(ir) >= atom_kind%hard_radius) CYCLE

                         DO ia = 1,grid_atom%ng_sphere

                            ra = particle_set(iatom)%r
                            ra(:) = ra(:) + grid_atom%rad(ir)*harmonics%a(:,ia)
                            vks_ra_idir2 = Eval_Interp_Spl3_pbc(ra,&
                                           vks_pw_spline(idir2,ispin)%pw,error)
                            vks_ra_idir3 = Eval_Interp_Spl3_pbc(ra,&
                                           vks_pw_spline(idir3,ispin)%pw,error)

                            IF(iat.LT.bo(1).OR.iat.GT.bo(2))CYCLE!quick and dirty: 
                            !                                    !here take care of the partition

                            ! + sum_A j_loc_h_A x nabla_vks_s_A
                            temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + &
                               (-1.0_dp)**(1 + ispin) * ( &
                               jrho1_atom%jrho_vec_rad_h(idir2,iB,ispin)%r_coef(ir,ia) * &
                               vks_ra_idir3 - &
                               jrho1_atom%jrho_vec_rad_h(idir3,iB,ispin)%r_coef(ir,ia) * &
                               vks_ra_idir2 &
                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)

                            ! - sum_A j_loc_s_A x nabla_vks_s_A
                            temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) - &
                               (-1.0_dp)**(1 + ispin) * ( &
                               jrho1_atom%jrho_vec_rad_s(idir2,iB,ispin)%r_coef(ir,ia) * &
                               vks_ra_idir3 - &
                               jrho1_atom%jrho_vec_rad_s(idir3,iB,ispin)%r_coef(ir,ia) * &
                               vks_ra_idir2 &
                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)

                            ! + sum_A j_loc_h_A x nabla_vks_loc_h_A
                            temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + &
                               (-1.0_dp)**(1 + ispin) * ( &
                               jrho1_atom%jrho_vec_rad_h(idir2,iB,ispin)%r_coef(ir,ia) * &
                               nablavks_atom%nablavks_vec_rad_h(idir3,ispin)%r_coef(ir,ia) - &
                               jrho1_atom%jrho_vec_rad_h(idir3,iB,ispin)%r_coef(ir,ia) * &
                               nablavks_atom%nablavks_vec_rad_h(idir2,ispin)%r_coef(ir,ia) &
                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)
 
                            ! - sum_A j_loc_h_A x nabla_vks_loc_s_A
                            temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) - &
                               (-1.0_dp)**(1 + ispin) * ( &
                               jrho1_atom%jrho_vec_rad_h(idir2,iB,ispin)%r_coef(ir,ia) * &
                               nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia) - &
                               jrho1_atom%jrho_vec_rad_h(idir3,iB,ispin)%r_coef(ir,ia) * &
                               nablavks_atom%nablavks_vec_rad_s(idir2,ispin)%r_coef(ir,ia) &
                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)

! ORIGINAL
!                            ! + sum_A j_loc_h_A x nabla_vks_loc_h_A
!                            temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + &
!                               (-1.0_dp)**(1.0_dp + REAL(ispin,KIND=dp)) * ( &
!                               jrho1_atom%jrho_vec_rad_h(idir2,iB,ispin)%r_coef(ir,ia) * &
!                               nablavks_atom%nablavks_vec_rad_h(idir3,ispin)%r_coef(ir,ia) - &
!                               jrho1_atom%jrho_vec_rad_h(idir3,iB,ispin)%r_coef(ir,ia) * &
!                               nablavks_atom%nablavks_vec_rad_h(idir2,ispin)%r_coef(ir,ia) &
!                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)
!                            ! - sum_A j_loc_s_A x nabla_vks_loc_s_A
!                            temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) - &
!                               (-1.0_dp)**(1.0_dp + REAL(ispin,KIND=dp)) * ( &
!                               jrho1_atom%jrho_vec_rad_s(idir2,iB,ispin)%r_coef(ir,ia) * &
!                               nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia) - &
!                               jrho1_atom%jrho_vec_rad_s(idir3,iB,ispin)%r_coef(ir,ia) * &
!                               nablavks_atom%nablavks_vec_rad_s(idir2,ispin)%r_coef(ir,ia) &
!                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)
                         END DO !ia
                      END DO !ir
                   END DO !ispin
                END DO !idir1
             END DO !iB
          END DO !iat
       END DO !ikind
   
       CALL mp_sum(temp_so_gapw,para_env%group)
       temp_so_gapw(:,:) = -1.0_dp * epr_env%g_so_factor_gapw * temp_so_gapw(:,:)

       IF (output_unit>0) THEN
          DO iB = 1,3
             DO idir1 = 1,3
                WRITE (UNIT=output_unit,FMT="(T2,A,T18,I1,I1,T56,E24.16)")&
                "epr|SOX:gapw",iB,idir1,temp_so_gapw(iB,idir1)
             END DO
          END DO
       END IF

       g_so(:,:) = g_so(:,:) + temp_so_gapw(:,:)

       DO ispin = 1,nspins
          DO idir1 = 1,3
             CALL pw_pool_give_back_pw(auxbas_pw_pool,vks_pw_spline(idir1,ispin)%pw,&
                                       error=error)
          END DO
       END DO

       DEALLOCATE(vks_pw_spline,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF ! gapw

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

    INTEGER                                  :: &
       aint_precond, ia, iat, iatom, iB, idir1, ikind, iso, ispin, ir, istat,&
       max_iter, precond_kind, natom, nkind, nspins, output_unit
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: gapw, paw_atom, failure,&
                                                soo_rho_hard, success
    REAL(dp)                                 :: eps_r, eps_x,&
                                                temp_soo_soft,&
                                                bind_ra_idir1,&
                                                rho_spin
    REAL(KIND=dp), DIMENSION(3)              :: ra
    REAL(dp), DIMENSION(3,3)                 :: temp_soo_gapw
    REAL(dp), DIMENSION(:,:), POINTER        :: bind0, g_soo, g_total
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(pw_p_type), DIMENSION(:,:),&
      POINTER                                :: bind_pw_spline
    TYPE(cp_logger_type), POINTER            :: logger

    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_rho_p_type), DIMENSION(:,:), &
      POINTER                                :: bind_set
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: rho_rad_h, rho_rad_s
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(section_vals_type), POINTER         :: g_section,&
                                                interp_section,&
                                                lr_section

    NULLIFY(atomic_kind_set, atom_kind, atom_list, bind0, bind_set,&
      dft_control, grid_atom, g_section, g_soo, g_total, harmonics,&
      interp_section, logger, lr_section, para_env, particle_set, rho,&
      rho_atom, rho_atom_set)

    logger => cp_error_get_logger(error)
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES",error=error)

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog",error=error)

    g_section => section_vals_get_subs_vals(lr_section,&
                 "EPR%PRINT%G_TENSOR",error=error)

    CALL section_vals_val_get(g_section,"soo_rho_hard",l_val=soo_rho_hard,error=error)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control,para_env=para_env,particle_set=particle_set,&
                    pw_env=pw_env,rho=rho,rho_atom_set=rho_atom_set,error=error)

    CALL get_epr_env(epr_env=epr_env,bind_set=bind_set,&
                     g_soo=g_soo,g_total=g_total,bind0=bind0)

    gapw = dft_control%qs_control%gapw
    nkind = SIZE(atomic_kind_set,1)
    nspins = dft_control%nspins

    DO iB=1,3
       DO idir1=1,3
          temp_soo_soft = 0.0_dp
          DO ispin=1,nspins
             temp_soo_soft = temp_soo_soft + (-1.0_dp)**(1 +ispin) * ( &
                               pw_integral_ab ( bind_set(idir1,iB)%rho%rho_r(1)%pw, &
                               rho%rho_r(ispin)%pw ) &
                               )
          END DO
          temp_soo_soft = -1.0_dp * epr_env%g_soo_factor * temp_soo_soft
          IF (output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(T2,A,T18,i1,i1,T56,E24.16)")&
             "epr|SOO:soft",iB,idir1,temp_soo_soft
          END IF
          g_soo(iB,idir1) = temp_soo_soft
       END DO
    END DO

    DO iB=1,3
       DO idir1=1,3
          temp_soo_soft = -1.0_dp * epr_env%g_soo_chicorr_factor * bind0(idir1,iB) * &
                          (REAL(dft_control%multiplicity,KIND=dp)-1.0_dp)
          IF (output_unit>0) THEN
             WRITE (UNIT=output_unit,FMT="(T2,A,T18,i1,i1,T56,E24.16)")&
             "epr|SOO:soft_g0",iB,idir1,temp_soo_soft
          END IF
          g_soo(iB,idir1) = g_soo(iB,idir1) + temp_soo_soft
       END DO
    END DO

    IF (gapw .AND. soo_rho_hard) THEN

       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
       ALLOCATE(bind_pw_spline(3,3),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       interp_section =>  section_vals_get_subs_vals(lr_section,&
                           "EPR%INTERPOLATOR",error=error)
       CALL section_vals_val_get(interp_section,"aint_precond", &
                           i_val=aint_precond, error=error)
       CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind,error=error)
       CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter,error=error)
       CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r,error=error)
       CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x,error=error)

       DO iB = 1,3
          DO idir1 = 1,3
             CALL pw_pool_create_pw(auxbas_pw_pool,&
                                    bind_pw_spline(idir1,iB)%pw,&
                                    use_data=REALDATA3D,in_space=REALSPACE,&
                                    error=error)
             ! calculate spline coefficients
             CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
                     pool=auxbas_pw_pool,pbc=.TRUE.,transpose=.FALSE.,error=error)
             CALL pw_spline_do_precond(precond,&
                     bind_set(idir1,iB)%rho%rho_r(1)%pw,&
                     bind_pw_spline(idir1,iB)%pw,error=error)
             CALL pw_spline_precond_set_kind(precond,precond_kind,error=error)
             success=find_coeffs(values=bind_set(idir1,iB)%rho%rho_r(1)%pw,&
                     coeffs=bind_pw_spline(idir1,iB)%pw,linOp=spl3_pbc,&
                     preconditioner=precond,pool=auxbas_pw_pool,&
                     eps_r=eps_r,eps_x=eps_x,max_iter=max_iter,&
                     print_section=interp_section,error=error)
             CPPostconditionNoFail(success,cp_warning_level,routineP,error)
             CALL pw_spline_precond_release(precond,error=error)
          END DO ! idir1
       END DO ! iB

       temp_soo_gapw = 0.0_dp

       DO ikind = 1,nkind
          NULLIFY (atom_kind,atom_list,grid_atom,harmonics)
          atom_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                               grid_atom=grid_atom,harmonics=harmonics,&
                               natom=natom, paw_atom=paw_atom)

          IF (.NOT.paw_atom) CYCLE

          ! Distribute the atoms of this kind
    
          bo = get_limit( natom, para_env%num_pe, para_env%mepos )

          DO iat = 1,natom !bo(1),bo(2)! this partitioning blocks the interpolation
             !                         ! routines (i.e. waiting for mp_sum)
             iatom = atom_list(iat)
             rho_atom => rho_atom_set(iatom)
             NULLIFY(rho_rad_h,rho_rad_s)
             CALL get_rho_atom(rho_atom=rho_atom, rho_rad_h=rho_rad_h,&
                               rho_rad_s=rho_rad_s)
             DO iB = 1,3
                DO idir1 = 1,3
                   DO ispin = 1,nspins
                      DO ir = 1,grid_atom%nr
   
                         IF (grid_atom%rad(ir) >= atom_kind%hard_radius) CYCLE
   
                         DO ia = 1,grid_atom%ng_sphere
   
                            ra = particle_set(iatom)%r
                            ra(:) = ra(:) + grid_atom%rad(ir)*harmonics%a(:,ia)
                            bind_ra_idir1 = Eval_Interp_Spl3_pbc(ra,&
                                            bind_pw_spline(idir1,iB)%pw,error)
   
                            IF(iat.LT.bo(1).OR.iat.GT.bo(2))CYCLE!quick and dirty: 
                            !                                    !here take care of the partition

                            rho_spin = 0.0_dp
   
                            DO iso = 1,harmonics%max_iso_not0
                               rho_spin = rho_spin + &
                                  (rho_rad_h(ispin)%r_coef(ir,iso) - &
                                  rho_rad_s(ispin)%r_coef(ir,iso))* &
                                  harmonics%slm(ia,iso)
                            END DO
   
                            temp_soo_gapw(iB,idir1) = temp_soo_gapw(iB,idir1) + &
                               (-1.0_dp)**(1 + ispin) * ( &
                               bind_ra_idir1 * rho_spin & 
                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)
   
                         END DO !ia
                      END DO !ir
                   END DO ! ispin
                END DO !idir1
             END DO !iB
          END DO !iat
       END DO !ikind

       CALL mp_sum(temp_soo_gapw,para_env%group)
       temp_soo_gapw(:,:) = -1.0_dp * epr_env%g_soo_factor * temp_soo_gapw(:,:)

       IF (output_unit>0) THEN
          DO iB = 1,3
             DO idir1 = 1,3
                WRITE (UNIT=output_unit,FMT="(T2,A,T18,I1,I1,T56,E24.16)")&
                "epr|SOO:gapw",iB,idir1,temp_soo_gapw(iB,idir1)
             END DO
          END DO
       END IF

       g_soo(:,:) = g_soo(:,:) + temp_soo_gapw(:,:)

       DO iB = 1,3
          DO idir1 = 1,3
             CALL pw_pool_give_back_pw(auxbas_pw_pool,bind_pw_spline(idir1,iB)%pw,&
                                       error=error)
          END DO
       END DO

       DEALLOCATE(bind_pw_spline,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF ! gapw

    g_total(:,:) = g_total(:,:) + g_soo(:,:)

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO",error=error)

  END SUBROUTINE epr_g_soo

END MODULE  qs_linres_epr_ownutils

