!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>       created 06-2006 [RD]
!> \author RD
! *****************************************************************************
MODULE qs_linres_epr_ownutils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathlib,                         ONLY: diamat_all
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_integral_ab,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc,&
                                             find_coeffs,&
                                             pw_spline_do_precond,&
                                             pw_spline_precond_create,&
                                             pw_spline_precond_release,&
                                             pw_spline_precond_set_kind,&
                                             pw_spline_precond_type,&
                                             spl3_pbc
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_core_energies,                ONLY: calculate_ptrace
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_linres_nmr_epr_common_utils,  ONLY: mult_G_ov_G2_grid
  USE qs_linres_op,                    ONLY: fac_vecp,&
                                             set_vecp,&
                                             set_vecp_rev
  USE qs_linres_types,                 ONLY: current_env_type,&
                                             epr_env_type,&
                                             get_current_env,&
                                             get_epr_env,&
                                             jrho_atom_type,&
                                             nablavks_atom_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_p_type,&
                                             qs_rho_type
  USE realspace_grid_types,            ONLY: realspace_grid_desc_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: epr_g_print, epr_g_zke, epr_g_so, epr_g_soo, epr_ind_magnetic_field

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_linres_epr_ownutils'

CONTAINS

! *****************************************************************************
!> \brief Prints the g tensor
!> \param epr_env ...
!> \param qs_env ...
!> \par History
!>      06.2006 created [RD]
!> \author RD
! *****************************************************************************
  SUBROUTINE epr_g_print(epr_env,qs_env)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_g_print', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: title
    INTEGER                                  :: idir1, idir2, output_unit, &
                                                unit_nr
    REAL(KIND=dp)                            :: eigenv_g(3), g_sym(3,3), gsum
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: lr_section

    NULLIFY(logger, lr_section)

    logger => cp_get_default_logger()
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES")

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog")

    gsum = 0.0_dp

    DO idir1 = 1,3
       DO idir2 = 1,3
          gsum = gsum + epr_env%g_total(idir1,idir2)
       END DO
    END DO

    IF (output_unit>0) THEN
       WRITE (UNIT=output_unit,FMT="(T2,A,T56,E14.6)")&
       "epr|TOT:checksum",gsum
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO")

    IF (BTEST(cp_print_key_should_output(logger%iter_info,lr_section,&
              "EPR%PRINT%G_TENSOR"),cp_p_file)) THEN

       unit_nr=cp_print_key_unit_nr(logger,lr_section,"EPR%PRINT%G_TENSOR",&
                                       extension=".data",middle_name="GTENSOR",&
                                       log_filename=.FALSE.)

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
            &                            "EPR%PRINT%G_TENSOR")

    ENDIF

  END SUBROUTINE epr_g_print

! *****************************************************************************
!> \brief Calculate zke part of the g tensor
!> \param epr_env ...
!> \param qs_env ...
!> \par History
!>      06.2006 created [RD]
!> \author RD
! *****************************************************************************
  SUBROUTINE epr_g_zke(epr_env,qs_env)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_g_zke', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, ispin, output_unit
    REAL(KIND=dp)                            :: epr_g_zke_temp(2)
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: kinetic, rho_ao
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: lr_section

    NULLIFY(dft_control, logger, lr_section, rho, kinetic, para_env, rho_ao)

    logger => cp_get_default_logger()
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES")

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog")

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    kinetic=kinetic,rho=rho,para_env=para_env)

    CALL qs_rho_get(rho, rho_ao=rho_ao)

    DO ispin=1,dft_control%nspins
       CALL calculate_ptrace(kinetic(1)%matrix,rho_ao(ispin)%matrix,&
                             ecore=epr_g_zke_temp(ispin))
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
         "PRINT%PROGRAM_RUN_INFO")

  END SUBROUTINE epr_g_zke

! *****************************************************************************
!> \brief Calculates g_so
!> \param epr_env ...
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
!> \par History
!>      06.2006 created [RD]
!> \author RD
! *****************************************************************************
  SUBROUTINE epr_g_so(epr_env,current_env,qs_env,iB)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_g_so', &
      routineP = moduleN//':'//routineN

    INTEGER :: aint_precond, ia, iat, iatom, idir1, idir2, idir3, ikind, ir, &
      ispin, istat, max_iter, natom, nkind, nspins, output_unit, precond_kind
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, gapw, paw_atom, &
                                                success
    REAL(dp)                                 :: eps_r, eps_x, hard_radius, &
                                                temp_so_soft, vks_ra_idir2, &
                                                vks_ra_idir3
    REAL(dp), DIMENSION(3, 3)                :: temp_so_gapw
    REAL(dp), DIMENSION(:, :), POINTER       :: g_so, g_total
    REAL(KIND=dp), DIMENSION(3)              :: ra
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(jrho_atom_type), DIMENSION(:), &
      POINTER                                :: jrho1_atom_set
    TYPE(jrho_atom_type), POINTER            :: jrho1_atom
    TYPE(nablavks_atom_type), DIMENSION(:), &
      POINTER                                :: nablavks_atom_set
    TYPE(nablavks_atom_type), POINTER        :: nablavks_atom
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: jrho2_r, jrho3_r, nrho1_r, &
                                                nrho2_r, nrho3_r
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: vks_pw_spline
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_p_type), DIMENSION(:), &
      POINTER                                :: jrho1_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      POINTER                                :: nablavks_set
    TYPE(section_vals_type), POINTER         :: interp_section, lr_section

    NULLIFY(atomic_kind_set, qs_kind_set, atom_list, dft_control,&
      grid_atom, g_so, g_total, harmonics, interp_section, jrho1_atom_set,&
      jrho1_set, logger, lr_section, nablavks_atom, nablavks_atom_set,&
      nablavks_set, para_env, particle_set, jrho2_r, jrho3_r, nrho1_r, nrho2_r, nrho3_r)

    logger => cp_get_default_logger()
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES")

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog")

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    para_env=para_env,pw_env=pw_env,&
                    particle_set=particle_set)

    CALL get_epr_env(epr_env=epr_env,&
                     nablavks_set=nablavks_set,&
                     nablavks_atom_set=nablavks_atom_set,&
                     g_total=g_total,g_so=g_so)

    CALL get_current_env(current_env=current_env,&
                         jrho1_set=jrho1_set,jrho1_atom_set=jrho1_atom_set)

    gapw = dft_control%qs_control%gapw
    nkind = SIZE(qs_kind_set,1)
    nspins = dft_control%nspins

    DO idir1=1,3
       CALL set_vecp(idir1,idir2,idir3)
       ! j_pw x nabla_vks_pw
       temp_so_soft = 0.0_dp
       DO ispin=1,nspins
          CALL qs_rho_get(jrho1_set(idir2)%rho, rho_r=jrho2_r)
          CALL qs_rho_get(jrho1_set(idir3)%rho, rho_r=jrho3_r)
          CALL qs_rho_get(nablavks_set(idir2,ispin)%rho, rho_r=nrho2_r)
          CALL qs_rho_get(nablavks_set(idir3,ispin)%rho, rho_r=nrho3_r)
          temp_so_soft = temp_so_soft + (-1.0_dp)**(1 + ispin) * ( &
                         pw_integral_ab(jrho2_r(ispin)%pw, nrho3_r(1)%pw) - &
                         pw_integral_ab(jrho3_r(ispin)%pw, nrho2_r(1)%pw)   )
       END DO
       temp_so_soft = -1.0_dp * epr_env%g_so_factor * temp_so_soft
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T18,I1,I1,T56,E24.16)")&
               "epr|SOX:soft",iB,idir1,temp_so_soft
       END IF
       g_so(iB,idir1) = temp_so_soft
    END DO !idir1

    IF (gapw) THEN

       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
       ALLOCATE(vks_pw_spline(3,nspins),STAT=istat)
       CPASSERT(istat==0)

       interp_section =>  section_vals_get_subs_vals(lr_section,&
                           "EPR%INTERPOLATOR")
       CALL section_vals_val_get(interp_section,"aint_precond", &
                           i_val=aint_precond)
       CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind)
       CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter)
       CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r)
       CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x)

       DO ispin = 1,nspins
          DO idir1 = 1,3
             CALL pw_pool_create_pw(auxbas_pw_pool,&
                                    vks_pw_spline(idir1,ispin)%pw,&
                                    use_data=REALDATA3D,in_space=REALSPACE)
             ! calculate spline coefficients
             CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
                     pool=auxbas_pw_pool,pbc=.TRUE.,transpose=.FALSE.)
             CALL qs_rho_get(nablavks_set(idir1,ispin)%rho, rho_r=nrho1_r)
             CALL pw_spline_do_precond(precond, nrho1_r(1)%pw,&
                     vks_pw_spline(idir1,ispin)%pw)
             CALL pw_spline_precond_set_kind(precond,precond_kind)
             success=find_coeffs(values=nrho1_r(1)%pw,&
                     coeffs=vks_pw_spline(idir1,ispin)%pw,linOp=spl3_pbc,&
                     preconditioner=precond,pool=auxbas_pw_pool,&
                     eps_r=eps_r,eps_x=eps_x,max_iter=max_iter)
             CPASSERT(success)
             CALL pw_spline_precond_release(precond)
          END DO ! idir1
       END DO ! ispin

       temp_so_gapw = 0.0_dp

       DO ikind = 1, nkind
          NULLIFY (atom_list,grid_atom,harmonics)
          CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
          CALL get_qs_kind(qs_kind_set(ikind),&
                           hard_radius=hard_radius,&
                           grid_atom=grid_atom,&
                           harmonics=harmonics,&
                           paw_atom=paw_atom)

          IF (.NOT.paw_atom) CYCLE

          ! Distribute the atoms of this kind

          bo = get_limit( natom, para_env%num_pe, para_env%mepos )

          DO iat = 1,natom !bo(1),bo(2)! this partitioning blocks the interpolation
             !                         ! routines (i.e. waiting for mp_sum)
             iatom = atom_list(iat)
             NULLIFY (jrho1_atom, nablavks_atom)
             jrho1_atom => jrho1_atom_set(iatom)
             nablavks_atom => nablavks_atom_set(iatom)
             DO idir1 = 1,3
                CALL set_vecp(idir1,idir2,idir3)
                DO ispin = 1,nspins
                   DO ir = 1,grid_atom%nr

                      IF (grid_atom%rad(ir) >= hard_radius) CYCLE

                      DO ia = 1,grid_atom%ng_sphere

                         ra = particle_set(iatom)%r
                         ra(:) = ra(:) + grid_atom%rad(ir)*harmonics%a(:,ia)
                         vks_ra_idir2 = Eval_Interp_Spl3_pbc(ra,&
                                        vks_pw_spline(idir2,ispin)%pw)
                         vks_ra_idir3 = Eval_Interp_Spl3_pbc(ra,&
                                        vks_pw_spline(idir3,ispin)%pw)

                         IF(iat.LT.bo(1).OR.iat.GT.bo(2))CYCLE!quick and dirty:
                         !                                    !here take care of the partition

                         ! + sum_A j_loc_h_A x nabla_vks_s_A
                         temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + &
                               (-1.0_dp)**(1 + ispin) * ( &
                               jrho1_atom%jrho_vec_rad_h(idir2,ispin)%r_coef(ir,ia) * &
                               vks_ra_idir3 - &
                               jrho1_atom%jrho_vec_rad_h(idir3,ispin)%r_coef(ir,ia) * &
                               vks_ra_idir2 &
                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)

                         ! - sum_A j_loc_s_A x nabla_vks_s_A
                         temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) - &
                               (-1.0_dp)**(1 + ispin) * ( &
                               jrho1_atom%jrho_vec_rad_s(idir2,ispin)%r_coef(ir,ia) * &
                               vks_ra_idir3 - &
                               jrho1_atom%jrho_vec_rad_s(idir3,ispin)%r_coef(ir,ia) * &
                               vks_ra_idir2 &
                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)

                         ! + sum_A j_loc_h_A x nabla_vks_loc_h_A
                         temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) + &
                               (-1.0_dp)**(1 + ispin) * ( &
                               jrho1_atom%jrho_vec_rad_h(idir2,ispin)%r_coef(ir,ia) * &
                               nablavks_atom%nablavks_vec_rad_h(idir3,ispin)%r_coef(ir,ia) - &
                               jrho1_atom%jrho_vec_rad_h(idir3,ispin)%r_coef(ir,ia) * &
                               nablavks_atom%nablavks_vec_rad_h(idir2,ispin)%r_coef(ir,ia) &
                               ) * grid_atom%wr(ir)*grid_atom%wa(ia)

                         ! - sum_A j_loc_h_A x nabla_vks_loc_s_A
                         temp_so_gapw(iB,idir1) = temp_so_gapw(iB,idir1) - &
                               (-1.0_dp)**(1 + ispin) * ( &
                               jrho1_atom%jrho_vec_rad_h(idir2,ispin)%r_coef(ir,ia) * &
                               nablavks_atom%nablavks_vec_rad_s(idir3,ispin)%r_coef(ir,ia) - &
                               jrho1_atom%jrho_vec_rad_h(idir3,ispin)%r_coef(ir,ia) * &
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
          END DO !iat
       END DO !ikind

       CALL mp_sum(temp_so_gapw,para_env%group)
       temp_so_gapw(:,:) = -1.0_dp * epr_env%g_so_factor_gapw * temp_so_gapw(:,:)

       IF (output_unit>0) THEN
          DO idir1 = 1,3
             WRITE (UNIT=output_unit,FMT="(T2,A,T18,I1,I1,T56,E24.16)")&
                  "epr|SOX:gapw",iB,idir1,temp_so_gapw(iB,idir1)
          END DO
       END IF

       g_so(iB,:) = g_so(iB,:) + temp_so_gapw(iB,:)

       DO ispin = 1,nspins
          DO idir1 = 1,3
             CALL pw_pool_give_back_pw(auxbas_pw_pool,vks_pw_spline(idir1,ispin)%pw)
          END DO
       END DO

       DEALLOCATE(vks_pw_spline,STAT=istat)
       CPASSERT(istat==0)

    END IF ! gapw

    g_total(iB,:) = g_total(iB,:) + g_so(iB,:)

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO")

  END SUBROUTINE epr_g_so

! *****************************************************************************
!> \brief Calculates g_soo (soft part only for now)
!> \param epr_env ...
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
!> \par History
!>      06.2006 created [RD]
!> \author RD
! *****************************************************************************
  SUBROUTINE epr_g_soo(epr_env,current_env,qs_env,iB)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_g_soo', &
      routineP = moduleN//':'//routineN

    INTEGER :: aint_precond, ia, iat, iatom, idir1, ikind, ir, iso, ispin, &
      istat, max_iter, natom, nkind, nspins, output_unit, precond_kind
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, gapw, paw_atom, &
                                                soo_rho_hard, success
    REAL(dp)                                 :: bind_ra_idir1, &
                                                chi_tensor(3,3,2), eps_r, &
                                                eps_x, hard_radius, rho_spin, &
                                                temp_soo_soft
    REAL(dp), DIMENSION(3, 3)                :: temp_soo_gapw
    REAL(dp), DIMENSION(:, :), POINTER       :: g_soo, g_total
    REAL(KIND=dp), DIMENSION(3)              :: ra
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: brho1_r, rho_r
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: bind_pw_spline
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_spline_precond_type), POINTER    :: precond
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_p_type), DIMENSION(:, :), &
      POINTER                                :: bind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: rho_rad_h, rho_rad_s
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(section_vals_type), POINTER         :: g_section, interp_section, &
                                                lr_section

    NULLIFY(atomic_kind_set, qs_kind_set, atom_list, bind_set, dft_control, &
         grid_atom, g_section, g_soo, g_total, harmonics, interp_section, &
         logger, lr_section, para_env, particle_set, rho, rho_atom, &
         rho_atom_set, rho_r, brho1_r)

    logger => cp_get_default_logger()
    lr_section => section_vals_get_subs_vals(qs_env%input,"PROPERTIES%LINRES")

    output_unit = cp_print_key_unit_nr(logger,lr_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".linresLog")

    g_section => section_vals_get_subs_vals(lr_section,&
                 "EPR%PRINT%G_TENSOR")

    CALL section_vals_val_get(g_section,"soo_rho_hard",l_val=soo_rho_hard)

    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
                    dft_control=dft_control,para_env=para_env,particle_set=particle_set,&
                    pw_env=pw_env,rho=rho,rho_atom_set=rho_atom_set)

    CALL get_epr_env(epr_env=epr_env,bind_set=bind_set,&
                     g_soo=g_soo,g_total=g_total)

    CALL get_current_env(current_env=current_env,&
                         chi_tensor=chi_tensor)
    CALL qs_rho_get(rho, rho_r=rho_r)

    gapw = dft_control%qs_control%gapw
    nkind = SIZE(qs_kind_set,1)
    nspins = dft_control%nspins

    DO idir1=1,3
       temp_soo_soft = 0.0_dp
       DO ispin=1,nspins
          CALL qs_rho_get(bind_set(idir1,iB)%rho, rho_r=brho1_r)
          temp_soo_soft = temp_soo_soft + (-1.0_dp)**(1 +ispin) * &
                  pw_integral_ab(brho1_r(1)%pw, rho_r(ispin)%pw)
       END DO
       temp_soo_soft = 1.0_dp * epr_env%g_soo_factor * temp_soo_soft
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T18,i1,i1,T56,E24.16)")&
               "epr|SOO:soft",iB,idir1,temp_soo_soft
       END IF
       g_soo(iB,idir1) = temp_soo_soft
    END DO

    DO idir1=1,3
       temp_soo_soft = 1.0_dp * epr_env%g_soo_chicorr_factor * chi_tensor(idir1,iB,2) * &
                       (REAL(dft_control%multiplicity,KIND=dp)-1.0_dp)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T18,i1,i1,T56,E24.16)")&
               "epr|SOO:soft_g0",iB,idir1,temp_soo_soft
       END IF
       g_soo(iB,idir1) = g_soo(iB,idir1) + temp_soo_soft
    END DO

    IF (gapw .AND. soo_rho_hard) THEN

       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
       ALLOCATE(bind_pw_spline(3,3),STAT=istat)
       CPASSERT(istat==0)

       interp_section =>  section_vals_get_subs_vals(lr_section,&
                           "EPR%INTERPOLATOR")
       CALL section_vals_val_get(interp_section,"aint_precond", &
                           i_val=aint_precond)
       CALL section_vals_val_get(interp_section,"precond",i_val=precond_kind)
       CALL section_vals_val_get(interp_section,"max_iter",i_val=max_iter)
       CALL section_vals_val_get(interp_section,"eps_r",r_val=eps_r)
       CALL section_vals_val_get(interp_section,"eps_x",r_val=eps_x)

       DO idir1 = 1,3
          CALL pw_pool_create_pw(auxbas_pw_pool,&
                                 bind_pw_spline(idir1,iB)%pw,&
                                 use_data=REALDATA3D,in_space=REALSPACE)
          ! calculate spline coefficients
          CALL pw_spline_precond_create(precond,precond_kind=aint_precond,&
                   pool=auxbas_pw_pool,pbc=.TRUE.,transpose=.FALSE.)
          CALL qs_rho_get(bind_set(idir1,iB)%rho, rho_r=brho1_r)
          CALL pw_spline_do_precond(precond, brho1_r(1)%pw,&
                                    bind_pw_spline(idir1,iB)%pw)
          CALL pw_spline_precond_set_kind(precond,precond_kind)
          success=find_coeffs(values=brho1_r(1)%pw,&
                  coeffs=bind_pw_spline(idir1,iB)%pw,linOp=spl3_pbc,&
                  preconditioner=precond,pool=auxbas_pw_pool,&
                  eps_r=eps_r,eps_x=eps_x,max_iter=max_iter)
          CPASSERT(success)
          CALL pw_spline_precond_release(precond)
       END DO ! idir1

       temp_soo_gapw = 0.0_dp

       DO ikind = 1,nkind
          NULLIFY (atom_list,grid_atom,harmonics)
          CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
          CALL get_qs_kind(qs_kind_set(ikind),&
                           hard_radius=hard_radius,&
                           grid_atom=grid_atom,&
                           harmonics=harmonics,&
                           paw_atom=paw_atom)

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
             DO idir1 = 1,3
                DO ispin = 1,nspins
                   DO ir = 1,grid_atom%nr

                      IF (grid_atom%rad(ir) >= hard_radius) CYCLE

                      DO ia = 1,grid_atom%ng_sphere

                         ra = particle_set(iatom)%r
                         ra(:) = ra(:) + grid_atom%rad(ir)*harmonics%a(:,ia)
                         bind_ra_idir1 = Eval_Interp_Spl3_pbc(ra,&
                                         bind_pw_spline(idir1,iB)%pw)

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
          END DO !iat
       END DO !ikind

       CALL mp_sum(temp_soo_gapw,para_env%group)
       temp_soo_gapw(:,:) = 1.0_dp * epr_env%g_soo_factor * temp_soo_gapw(:,:)

       IF (output_unit>0) THEN
          DO idir1 = 1,3
             WRITE (UNIT=output_unit,FMT="(T2,A,T18,I1,I1,T56,E24.16)")&
                  "epr|SOO:gapw",iB,idir1,temp_soo_gapw(iB,idir1)
          END DO
       END IF

       g_soo(iB,:) = g_soo(iB,:) + temp_soo_gapw(iB,:)

       DO idir1 = 1,3
          CALL pw_pool_give_back_pw(auxbas_pw_pool,bind_pw_spline(idir1,iB)%pw)
       END DO

       DEALLOCATE(bind_pw_spline,STAT=istat)
       CPASSERT(istat==0)

    END IF ! gapw

    g_total(iB,:) = g_total(iB,:) + g_soo(iB,:)

    CALL cp_print_key_finished_output(output_unit,logger,lr_section,&
         "PRINT%PROGRAM_RUN_INFO")

  END SUBROUTINE epr_g_soo

! *****************************************************************************
!> \brief ...
!> \param epr_env ...
!> \param current_env ...
!> \param qs_env ...
!> \param iB ...
! *****************************************************************************
  SUBROUTINE epr_ind_magnetic_field(epr_env,current_env,qs_env,iB)

    TYPE(epr_env_type)                       :: epr_env
    TYPE(current_env_type)                   :: current_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: iB

    CHARACTER(LEN=*), PARAMETER :: routineN = 'epr_ind_magnetic_field', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idir, idir2, idir3, &
                                                iiB, iiiB, ispin, istat, &
                                                natom, nspins
    LOGICAL                                  :: failure, gapw
    REAL(dp)                                 :: scale_fac
    TYPE(cell_type), POINTER                 :: cell
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: pw_gspace_work, &
                                                shift_pw_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: epr_rho_r, jrho1_g
    TYPE(pw_p_type), DIMENSION(:, :), &
      POINTER                                :: shift_pw_gspace
    TYPE(pw_p_type), POINTER                 :: rho_gspace
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(realspace_grid_desc_type), POINTER  :: auxbas_rs_desc

    failure = .FALSE.
    CALL timeset(routineN,handle)

    NULLIFY(cell, dft_control, pw_env, auxbas_rs_desc, auxbas_pw_pool, &
         rho_gspace, pw_pools, particle_set, jrho1_g, epr_rho_r)

    CALL get_qs_env(qs_env=qs_env,cell=cell,dft_control=dft_control,&
                    particle_set=particle_set)

    gapw = dft_control%qs_control%gapw
    natom = SIZE(particle_set,1)
    nspins = dft_control%nspins

    CALL get_epr_env(epr_env=epr_env)

    CALL get_current_env(current_env=current_env)


    CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
    CALL pw_env_get(pw_env,auxbas_rs_desc=auxbas_rs_desc,&
                    auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools)
    !
    ! Initialize
    ! Allocate grids for the calculation of jrho and the shift
    ALLOCATE(shift_pw_gspace(3,nspins),STAT=istat)
    CPASSERT(istat==0)
    DO ispin = 1,nspins
       DO idir = 1,3
          CALL pw_pool_create_pw(auxbas_pw_pool,shift_pw_gspace(idir,ispin)%pw,&
                                  use_data=COMPLEXDATA1D,&
                                  in_space=RECIPROCALSPACE)
          CALL pw_zero(shift_pw_gspace(idir,ispin)%pw)
       ENDDO
    ENDDO
    CALL pw_pool_create_pw(auxbas_pw_pool,shift_pw_rspace%pw,&
                           use_data=REALDATA3D,in_space=REALSPACE)
    CALL pw_zero(shift_pw_rspace%pw)
    CALL pw_pool_create_pw(auxbas_pw_pool,pw_gspace_work%pw,&
                           use_data=COMPLEXDATA1D,&
                           in_space=RECIPROCALSPACE)

    CALL pw_zero(pw_gspace_work%pw)
    !
    CALL set_vecp(iB,iiB,iiiB)
    !
    DO ispin = 1,nspins
       !
       DO idir3=1,3
          ! set to zero for the calculation of the shift
          CALL pw_zero(shift_pw_gspace(idir3,ispin)%pw)
       ENDDO
       DO idir = 1,3
          CALL qs_rho_get(current_env%jrho1_set(idir)%rho, rho_g=jrho1_g)
          rho_gspace => jrho1_g(ispin)
          ! Field gradient
          ! loop over the Gvec  components: x,y,z
          DO idir2 = 1,3
             IF(idir /= idir2) THEN
                ! in reciprocal space multiply (G_idir2(i)/G(i)^2)J_(idir)(G(i))
                CALL mult_G_ov_G2_grid(cell,auxbas_pw_pool,rho_gspace,&
                                       pw_gspace_work,idir2,0.0_dp)
                !
                ! scale and add to the correct component of the shift column
                CALL set_vecp_rev(idir,idir2,idir3)
                scale_fac=fac_vecp(idir3,idir2,idir)
                CALL pw_scale(pw_gspace_work%pw,scale_fac)
                CALL pw_axpy(pw_gspace_work%pw,shift_pw_gspace(idir3,ispin)%pw)
             ENDIF
          ENDDO
          !
       ENDDO ! idir
    ENDDO ! ispin

    ! Store the total soft induced magnetic field (corrected for sic)
    IF (dft_control%nspins == 2) THEN
       DO idir = 1,3
          CALL qs_rho_get(epr_env%bind_set(idir,iB)%rho, rho_r=epr_rho_r)
          CALL pw_transfer(shift_pw_gspace(idir,2)%pw,epr_rho_r(1)%pw)
       ENDDO
    END IF
    !
    ! Dellocate grids for the calculation of jrho and the shift
    CALL pw_pool_give_back_pw(auxbas_pw_pool,pw_gspace_work%pw)
    DO ispin = 1,dft_control%nspins
       DO idir = 1,3
          CALL pw_pool_give_back_pw(auxbas_pw_pool,shift_pw_gspace(idir,ispin)%pw)
       ENDDO
    ENDDO
    DEALLOCATE(shift_pw_gspace,STAT=istat)
    CPASSERT(istat==0)
    CALL pw_pool_give_back_pw(auxbas_pw_pool,shift_pw_rspace%pw)
    !
    ! Finalize
    CALL timestop(handle)
    !
  END SUBROUTINE epr_ind_magnetic_field

END MODULE qs_linres_epr_ownutils

