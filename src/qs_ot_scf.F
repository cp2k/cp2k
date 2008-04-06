!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief basic functionality for using ot in the scf routines.
!> \par History
!>      01.2003 : Joost VandeVondele : adapted for LSD
!> \author Joost VandeVondele (25.08.2002)
! *****************************************************************************
MODULE qs_ot_scf
  USE bibliography,                    ONLY: Weber2008,&
                                             cite_reference
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE f77_blas
  USE input_constants,                 ONLY: &
       ls_2pnt, ls_3pnt, ls_gold, ls_none, ot_algo_irac, &
       ot_algo_taylor_or_diag, ot_chol_irac, ot_lwdn_irac, ot_mini_cg, &
       ot_mini_diis, ot_mini_sd, ot_poly_irac, ot_precond_full_all, &
       ot_precond_full_kinetic, ot_precond_full_single, &
       ot_precond_full_single_inverse, ot_precond_none, ot_precond_s_inverse, &
       ot_precond_sparse_diag
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_restrict
  USE qs_ot,                           ONLY: &
       qs_ot_allocate, qs_ot_destroy, qs_ot_get_orbitals, &
       qs_ot_get_orbitals_ref, qs_ot_get_p, qs_ot_get_scp_coeffs, qs_ot_init, &
       qs_ot_settings_init, qs_ot_type
  USE qs_ot_minimizer,                 ONLY: ot_mini
  USE scp_coeff_types,                 ONLY: aux_coeff_set_type,&
                                             aux_coeff_type,&
                                             get_aux_coeff
  USE sparse_matrix_types,             ONLY: real_matrix_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ot_scf'
  ! *** Public subroutines ***

  PUBLIC :: ot_scf_init
  PUBLIC :: ot_scf_mini
  PUBLIC :: ot_scf_destroy
  PUBLIC :: ot_scf_read_input

CONTAINS

! *****************************************************************************
  SUBROUTINE ot_scf_read_input(qs_ot_env,para_env,scf_section,error)
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: scf_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=40)                        :: value
    INTEGER                                  :: handle, ispin, ls_method, &
                                                nspin, ot_algorithm, &
                                                ot_method, ot_ortho_irac, &
                                                output_unit
    LOGICAL                                  :: explicit
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: ot_section

    CALL timeset('ot_scf_read_input','I','',handle)

    logger => cp_error_get_logger(error)
    output_unit=cp_print_key_unit_nr(logger,scf_section,"PRINT%PROGRAM_RUN_INFO",&
         extension=".log",error=error)

    ! decide default settings
    CALL qs_ot_settings_init(qs_ot_env(1)%settings)

    ! use ot input new style
    ot_section=>section_vals_get_subs_vals(scf_section,"OT",error=error)
    CALL section_vals_get(ot_section,explicit=explicit,error=error)

    ! choose algorithm
    CALL section_vals_val_get(ot_section,"ALGORITHM",i_val=ot_algorithm,error=error)
    SELECT CASE(ot_algorithm)
    CASE (ot_algo_taylor_or_diag)
       qs_ot_env(1)%settings%ot_algorithm="TOD"
    CASE (ot_algo_irac)
       CALL cite_reference(Weber2008)
       qs_ot_env(1)%settings%ot_algorithm="REF"
    CASE DEFAULT
       CALL stop_program("READ OT ALGORITHM","VALUE UNKNOWN")
    END SELECT

    ! irac input
    CALL section_vals_val_get(ot_section,"IRAC_DEGREE",i_val=qs_ot_env(1)%settings%irac_degree,&
         &                    error=error)
    IF(   qs_ot_env(1)%settings%irac_degree.LT.2.OR.&
         &qs_ot_env(1)%settings%irac_degree.GT.4) THEN
       CALL stop_program("READ OT IRAC_DEGREE","VALUE UNKNOWN")
    ENDIF
    CALL section_vals_val_get(ot_section,"MAX_IRAC",i_val=qs_ot_env(1)%settings%max_irac,&
         &                    error=error)
    IF(qs_ot_env(1)%settings%max_irac.LT.1) THEN
       CALL stop_program("READ OT MAX_IRAC","VALUE MUST BE GREATER THAN ZERO")
    ENDIF
    CALL section_vals_val_get(ot_section,"EPS_IRAC",r_val=qs_ot_env(1)%settings%eps_irac,&
         &                    error=error)
    IF(qs_ot_env(1)%settings%eps_irac.LT.0.0_dp) THEN
       CALL stop_program("READ OT EPS_IRAC","VALUE MUST BE GREATER THAN ZERO")
    ENDIF
    CALL section_vals_val_get(ot_section,"EPS_IRAC_QUICK_EXIT",r_val=qs_ot_env(1)%settings%eps_irac_quick_exit,&
         &                    error=error)
    IF(qs_ot_env(1)%settings%eps_irac_quick_exit.LT.0.0_dp) THEN
       CALL stop_program("READ OT EPS_IRAC_QUICK_EXIT","VALUE MUST BE GREATER THAN ZERO")
    ENDIF

    CALL section_vals_val_get(ot_section,"EPS_IRAC_SWITCH",r_val=qs_ot_env(1)%settings%eps_irac_switch,&
         &                    error=error)
    IF(qs_ot_env(1)%settings%eps_irac_switch.LT.0.0_dp) THEN
       CALL stop_program("READ OT EPS_IRAC_SWITCH","VALUE MUST BE GREATER THAN ZERO")
    ENDIF

    CALL section_vals_val_get(ot_section,"ORTHO_IRAC",i_val=ot_ortho_irac,error=error)
    SELECT CASE(ot_ortho_irac)
    CASE(ot_chol_irac)
       qs_ot_env(1)%settings%ortho_irac="CHOL"
    CASE(ot_poly_irac)
       qs_ot_env(1)%settings%ortho_irac="POLY"
    CASE(ot_lwdn_irac)
       qs_ot_env(1)%settings%ortho_irac="LWDN"
    CASE DEFAULT
       CALL stop_program("READ OT ORTHO_IRAC","VALUE UNKNOWN")
    END SELECT

    CALL section_vals_val_get(ot_section,"ON_THE_FLY_LOC",l_val=qs_ot_env(1)%settings%on_the_fly_loc, &
         &                    error=error)

    CALL section_vals_val_get(ot_section,"MINIMIZER",i_val=ot_method,error=error)
    ! compatibility
    SELECT CASE(ot_method)
    CASE (ot_mini_sd)
       qs_ot_env(1)%settings%ot_method="SD"
    CASE (ot_mini_cg)
       qs_ot_env(1)%settings%ot_method="CG"
    CASE (ot_mini_diis)
       qs_ot_env(1)%settings%ot_method="DIIS"
       CALL section_vals_val_get(ot_section,"N_DIIS",i_val=qs_ot_env(1)%settings%diis_m,error=error)
    CASE DEFAULT
       CALL stop_program("READ OTSCF MINIMIZER","VALUE UNKNOWN")
    END SELECT
    CALL section_vals_val_get(ot_section,"SAFER_DIIS",l_val=qs_ot_env(1)%settings%safer_diis,error=error)
    CALL section_vals_val_get(ot_section,"LINESEARCH",i_val=ls_method,error=error)
    SELECT CASE(ls_method)
    CASE (ls_none)
       qs_ot_env(1)%settings%line_search_method="NONE"
    CASE (ls_2pnt)
       qs_ot_env(1)%settings%line_search_method="2PNT"
    CASE (ls_3pnt)
       qs_ot_env(1)%settings%line_search_method="3PNT"
    CASE (ls_gold)
       qs_ot_env(1)%settings%line_search_method="GOLD"
       CALL section_vals_val_get(ot_section,"GOLD_TARGET",r_val=qs_ot_env(1)%settings%gold_target,error=error)
    CASE DEFAULT
       CALL stop_program("READ OTSCF LS","VALUE UNKNOWN")
    END SELECT
    CALL section_vals_val_get(ot_section,"PRECONDITIONER",i_val=qs_ot_env(1)%settings%preconditioner_type,error=error)
    SELECT CASE(qs_ot_env(1)%settings%preconditioner_type)
    CASE(ot_precond_none)
       qs_ot_env(1)%settings%preconditioner_name="NONE"
    CASE(ot_precond_full_single)
       qs_ot_env(1)%settings%preconditioner_name="FULL_SINGLE"
    CASE(ot_precond_full_single_inverse)
       qs_ot_env(1)%settings%preconditioner_name="FULL_SINGLE_INVERSE"
    CASE(ot_precond_full_all)
       qs_ot_env(1)%settings%preconditioner_name="FULL_ALL"
    CASE(ot_precond_full_kinetic)
       qs_ot_env(1)%settings%preconditioner_name="FULL_KINETIC"
    CASE(ot_precond_s_inverse)
       qs_ot_env(1)%settings%preconditioner_name="FULL_S_INVERSE"
    CASE(ot_precond_sparse_diag)
       qs_ot_env(1)%settings%preconditioner_name="SPARSE_DIAG"
    CASE DEFAULT
       CALL stop_program("READ OTSCF PRECONDITIONER","VALUE UNKNOWN")
    END SELECT
    CALL section_vals_val_get(ot_section,"STEPSIZE",r_val=qs_ot_env(1)%settings%ds_min,error=error)
    CALL section_vals_val_get(ot_section,"ENERGY_GAP",r_val=qs_ot_env(1)%settings%energy_gap,error=error)
    CALL section_vals_val_get(ot_section,"EPS_TAYLOR",r_val=qs_ot_env(1)%settings%eps_taylor,error=error)
    CALL section_vals_val_get(ot_section,"MAX_TAYLOR",i_val=qs_ot_env(1)%settings%max_taylor,error=error)
    CALL section_vals_val_get(ot_section,"ROTATION",l_val=qs_ot_env(1)%settings%do_rotation,error=error)
    CALL section_vals_val_get(ot_section,"SCP",l_val=qs_ot_env(1)%settings%scp,error=error)

    ! write OT output

    IF (output_unit>0) THEN
       WRITE(output_unit,'(A)') " "
       WRITE(output_unit,'(A)') "  ----------------------------------- OT --------------------------------------"
       WRITE(output_unit,'(A)') " "
       WRITE(output_unit,'(A,L2)') "  Allowing for rotations: ",qs_ot_env(1)%settings%do_rotation
       SELECT CASE (qs_ot_env(1)%settings%OT_METHOD)
       CASE ("SD")
          WRITE(output_unit,'(A)') "  minimizer      : SD                  : steepest descent"
       CASE ("CG")
          WRITE(output_unit,'(A)') "  minimizer      : CG                  : conjugate gradient"
       CASE ("DIIS")
          WRITE(output_unit,'(A)') "  minimizer      : DIIS                : direct inversion "
          WRITE(output_unit,'(A)') "                                         in the iterative subspace"
          WRITE(output_unit,'(A,I3,A)')    "                            using      : - ",&
               qs_ot_env(1)%settings%diis_m," diis vectors"
          IF (qs_ot_env(1)%settings%safer_diis) THEN
             WRITE(output_unit,'(A,I3,A)') "                                         - safer DIIS on"
          ELSE
             WRITE(output_unit,'(A,I3,A)') "                                         - safer DIIS off"
          ENDIF
       CASE DEFAULT
          WRITE(output_unit,'(3A)') "  minimizer      :      ",qs_ot_env(1)%settings%OT_METHOD," : UNKNOWN"
       END SELECT
       SELECT CASE (qs_ot_env(1)%settings%preconditioner_name)
       CASE ("FULL_SINGLE")
          WRITE(output_unit,'(A)') "  preconditioner : FULL_SINGLE         : diagonalization based "
       CASE ("FULL_SINGLE_INVERSE")
          WRITE(output_unit,'(A)') "  preconditioner : FULL_SINGLE_INVERSE : cholesky inversion of H + eS "
       CASE ("FULL_ALL")
          WRITE(output_unit,'(A)') "  preconditioner : FULL_ALL            : diagonalization, state selective "
       CASE ("FULL_KINETIC")
          WRITE(output_unit,'(A)') "  preconditioner : FULL_KINETIC        : cholesky inversion of T + eS"
       CASE ("FULL_S_INVERSE")
          WRITE(output_unit,'(A)')            "  preconditioner : FULL_S_INVERSE      : cholesky inversion of S"
       CASE ("SPARSE_DIAG")
          WRITE(output_unit,'(A)')  &
               "  preconditioner : SPARSE_DIAG    : diagonal atomic block diagonalization"
       CASE ("NONE")
          WRITE(output_unit,'(A)') "  preconditioner : NONE  "
       CASE DEFAULT
          WRITE(output_unit,'(3A)') "  preconditioner : ",qs_ot_env(1)%settings%preconditioner_name," : UNKNOWN"
       END SELECT
       IF (qs_ot_env(1)%settings%OT_METHOD .EQ."SD".OR.qs_ot_env(1)%settings%OT_METHOD.EQ."CG") THEN
          SELECT CASE (qs_ot_env(1)%settings%line_search_method)
          CASE ("2PNT")
             WRITE(output_unit,'(A)') "  line search    : 2PNT                : 2 energies, one gradient "
          CASE ("3PNT")
             WRITE(output_unit,'(A)') "  line search    : 3PNT                : 3 energies "
          CASE ("GOLD")
             WRITE(output_unit,'(A)')       "  line search    : GOLD                : bracketing and golden section search "
             WRITE(output_unit,'(A,F14.8)') "                   target rel accuracy : ", &
                  qs_ot_env(1)%settings%gold_target
          CASE ("NONE")
             WRITE(output_unit,'(A)') "  line search    : NONE "
          CASE DEFAULT
             WRITE(output_unit,'(3A)') "  line search : ", &
                  qs_ot_env(1)%settings%line_search_method," : UNKNOWN"
          END SELECT
       ENDIF
       WRITE(output_unit,'(A,F14.8)') "  stepsize       :",qs_ot_env(1)%settings%ds_min
       WRITE(output_unit,'(A,F14.8)') "  energy_gap     :",qs_ot_env(1)%settings%energy_gap
       WRITE(output_unit,'(A,E14.5)') "  eps_taylor     :",qs_ot_env(1)%settings%eps_taylor
       WRITE(output_unit,'(A,I14)')   "  max_taylor     :",qs_ot_env(1)%settings%max_taylor
       WRITE(output_unit,'(A)') "  "
       IF(qs_ot_env(1)%settings%ot_algorithm.EQ.'REF') THEN
          WRITE(output_unit,'(A,1X,A)') "  ortho_irac         :",&
               & qs_ot_env(1)%settings%ortho_irac
          WRITE(output_unit,'(A,I14)')  "  irac_degree        :",&
               & qs_ot_env(1)%settings%irac_degree
          WRITE(output_unit,'(A,I14)')  "  max_irac           :",&
               & qs_ot_env(1)%settings%max_irac
          WRITE(output_unit,'(A,E14.5)')"  eps_irac           :",&
               & qs_ot_env(1)%settings%eps_irac
          WRITE(output_unit,'(A,E14.5)')"  eps_irac_switch    :",&
               & qs_ot_env(1)%settings%eps_irac_switch
          WRITE(output_unit,'(A,E14.5)')"  eps_irac_quick_exit:",&
               & qs_ot_env(1)%settings%eps_irac_quick_exit
          WRITE(output_unit,'(A,L2)')   "  on_the_fly_loc     :",&
               & qs_ot_env(1)%settings%on_the_fly_loc
          WRITE(output_unit,'(A)') "  "
       ENDIF
       WRITE(output_unit,'(A)') "  ----------------------------------- OT --------------------------------------"
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,scf_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

    ! copy the ot settings type so it is identical
    nspin=SIZE(qs_ot_env)
    DO ispin=2,nspin
       qs_ot_env(ispin)%settings=qs_ot_env(1)%settings
    ENDDO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE ot_scf_read_input
  !
  ! performs the actual minimisation, needs only limited info
  ! updated for restricted calculations
  ! matrix_dedc is the derivative of the energy with respect to the orbitals
  ! a null pointer for matrix_s implies that matrix_s is the unit matrix
  !
  !
! *****************************************************************************
  SUBROUTINE ot_scf_mini(mo_array,matrix_dedc,matrix_s,energy, &
       aux_coeff_set, energy_only, delta, qs_ot_env,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dedc
    TYPE(real_matrix_type), POINTER          :: matrix_s
    REAL(KIND=dp)                            :: energy
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    LOGICAL, INTENT(INOUT)                   :: energy_only
    REAL(KIND=dp)                            :: delta
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ot_scf_mini', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, k, n, nspin
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset('ot_scf_mini','I','',handle)

    ! we use the size of the ot_env to indicate the spin.
    ! this needs correction in case of restricted calculations
    nspin=SIZE(qs_ot_env)

    IF ( .NOT. energy_only) THEN
       IF (qs_ot_env(1)%settings%do_rotation) THEN
          DO ispin=1,nspin
             CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
             CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k,error=error)
             CALL cp_fm_gemm('T','N',k,k,n,1.0_dp,mo_coeff,matrix_dedc(ispin)%matrix, &
                  0.0_dp,qs_ot_env(ispin)%rot_mat_chcf,error=error)
          ENDDO
       ENDIF
    ENDIF

    ! notice we use (1) for driving all output and the minimization in case of LSD
    qs_ot_env(1)%etotal=energy
    CALL ot_mini(qs_ot_env,matrix_dedc,aux_coeff_set, error=error)
    delta           =qs_ot_env(1)%delta
    energy_only     =qs_ot_env(1)%energy_only

    ! generate the orbitals using the new matrix_x
    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)
       CALL cp_fm_get_info(mo_coeff,nrow_global=n,ncol_global=k,error=error)
       SELECT CASE(qs_ot_env(1)%settings%ot_algorithm)
       CASE("TOD")
          IF (ASSOCIATED(matrix_s)) THEN
             CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(ispin)%matrix_x, &
                  qs_ot_env(ispin)%matrix_sx, k,error=error)
          ELSE
             CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,error=error)
          ENDIF
          CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin),error=error)
          CALL qs_ot_get_orbitals(mo_coeff,qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin),error=error)
       CASE("REF")
          CALL qs_ot_get_orbitals_ref(mo_coeff,matrix_s,qs_ot_env(ispin)%matrix_x,&
               &           qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin)%matrix_gx_old,&
               &           qs_ot_env(ispin)%matrix_dx,qs_ot_env(ispin),qs_ot_env(1),&
               &           error=error)
       CASE DEFAULT
          CALL stop_program("ot_scf_mini","ALGORITHM NYI")
       END SELECT
    ENDDO
    IF (qs_ot_env(1)%restricted) THEN
       CALL mo_set_restrict(mo_array,error=error)
    ENDIF
    ! **** SCP
    ! generate the SCP aux coeffs using the new x
    IF (qs_ot_env(1)%settings%scp) CALL qs_ot_get_scp_coeffs ( qs_ot_env ( 1 ), aux_coeff_set, error )
    ! **** SCP
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE ot_scf_mini
  !
  ! initialises qs_ot_env so that mo_coeff is the current point
  ! and that the mimizization can be started.
  !
! *****************************************************************************
  SUBROUTINE ot_scf_init ( mo_array, matrix_s, aux_coeff_set, qs_ot_env, error )

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(real_matrix_type), POINTER          :: matrix_s
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    TYPE(qs_ot_type), DIMENSION(:), POINTER  :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: handle, i, icoef, icoef_atom, &
                                                ikind, ispin, k, n, n_els, &
                                                ncoef_atom, nkind, nspin
    REAL(dp), DIMENSION(:, :), POINTER       :: c, fc
    TYPE(aux_coeff_type), POINTER            :: local_coeffs
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    CALL timeset('ot_scf_init','I','',handle)

    ! **** SCP
    NULLIFY ( local_coeffs, c, fc )
    ! **** SCP
    ! adapted for work with the restricted keyword
    nspin=SIZE(qs_ot_env)

    DO ispin=1,nspin
       CALL get_mo_set(mo_set=mo_array(ispin)%mo_set,mo_coeff=mo_coeff)

       CALL cp_fm_get_info(mo_coeff, nrow_global=n, ncol_global=k,error=error)

       ! **** SCP
       ! allocate
       CALL qs_ot_allocate(qs_ot_env(ispin),mo_coeff%matrix_struct,  &
            aux_coeff_set = aux_coeff_set,    &
            error=error)
       ! **** SCP

       ! set c0,sc0
       CALL cp_fm_to_fm(mo_coeff,qs_ot_env(ispin)%matrix_c0,error=error)
       IF (ASSOCIATED(matrix_s)) THEN
          CALL cp_sm_fm_multiply(matrix_s,qs_ot_env(ispin)%matrix_c0, &
               qs_ot_env(ispin)%matrix_sc0,k,error=error)
       ELSE
          CALL cp_fm_to_fm ( qs_ot_env ( ispin ) % matrix_c0,  &
               qs_ot_env ( ispin ) % matrix_sc0, error = error )
       ENDIF

       ! init
       CALL qs_ot_init(qs_ot_env(ispin),error=error)

       ! set x
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_x,0.0_dp,error=error)
       CALL cp_fm_set_all(qs_ot_env(ispin)%matrix_sx,0.0_dp,error=error)

       IF (qs_ot_env(ispin)%settings%do_rotation) THEN
          CALL cp_fm_set_all(qs_ot_env(ispin)%rot_mat_x,0.0_dp,error=error)
       ENDIF
       IF (qs_ot_env(ispin)%settings%scp) THEN
          qs_ot_env ( ispin ) % x ( : ) = 0.0_dp
       ENDIF

       SELECT CASE(qs_ot_env(1)%settings%ot_algorithm)
       CASE("TOD")
          ! get c
          CALL qs_ot_get_p(qs_ot_env(ispin)%matrix_x,qs_ot_env(ispin)%matrix_sx,qs_ot_env(ispin),error=error)
       CASE("REF")
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_c0,qs_ot_env(ispin)%matrix_x,&
               &           error=error)
          CALL cp_fm_to_fm(qs_ot_env(ispin)%matrix_sc0,qs_ot_env(ispin)%matrix_sx,&
               &           error=error)
       CASE DEFAULT
          CALL stop_program("ot_scf_init","ALGORITHM NYI")
       END SELECT
       ! **** SCP
       ! set c = qs_ot_env % x and -fc = qs_ot_env % gx  for SCP
       IF (qs_ot_env(1)%settings%scp .AND. ispin == 1) THEN
          icoef = 0
          nkind = SIZE ( aux_coeff_set % coeffs_of_kind )
          DO ikind = 1, nkind
             local_coeffs => aux_coeff_set % coeffs_of_kind ( ikind ) % coeffs
             IF ( ASSOCIATED ( local_coeffs ) ) THEN
                CALL get_aux_coeff ( coeffs = local_coeffs, c = c, fc = fc,  & 
                     n_els = n_els, ncoef_atom = ncoef_atom, &
                     error = error  )
                DO i = 1, n_els
                   DO icoef_atom = 1, ncoef_atom
                      icoef = icoef + 1
                      qs_ot_env ( 1 ) % x ( icoef ) = c ( i, icoef_atom )
                      qs_ot_env ( 1 ) % gx ( icoef ) = -fc ( i, icoef_atom )
                   END DO
                END DO
             END IF
          END DO
       END IF
       ! **** SCP

    END DO
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE ot_scf_init

! *****************************************************************************
  SUBROUTINE ot_scf_destroy(qs_ot_env,error)

    TYPE(qs_ot_type)                         :: qs_ot_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CALL qs_ot_destroy(qs_ot_env,error=error)

  END SUBROUTINE ot_scf_destroy

END MODULE qs_ot_scf

