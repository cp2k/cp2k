!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types needed for MP2 calculations
!> \par History
!>       2011.05 created [Mauro Del Ben]
!> \author MDB
! *****************************************************************************
MODULE mp2_setup
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE mp2_types,                       ONLY: mp2_method_direct,&
                                             mp2_method_gpw,&
                                             mp2_method_none,&
                                             mp2_ri_optimize_basis,&
                                             mp2_type,&
                                             ri_mp2_laplace,&
                                             ri_mp2_method_gpw,&
                                             ri_rpa_method_gpw
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_setup'

  PUBLIC :: read_mp2_section

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param input ...
!> \param mp2_env ...
! *****************************************************************************
  SUBROUTINE read_mp2_section(input,mp2_env)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(mp2_type), INTENT(INOUT)            :: mp2_env

    CHARACTER(len=*), PARAMETER :: routineN = 'read_mp2_section', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ival, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: mp2_section

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    failure=.FALSE.

    mp2_section => section_vals_get_subs_vals(input,"DFT%XC%WF_CORRELATION")


    ! should come from input
    CALL section_vals_val_get(mp2_section,"METHOD",i_val=mp2_env%method)
    CALL section_vals_val_get(mp2_section,"MEMORY",r_val=mp2_env%mp2_memory)
    CALL section_vals_val_get(mp2_section,"SCALE_S",r_val=mp2_env%scale_S)
    CALL section_vals_val_get(mp2_section,"SCALE_T",r_val=mp2_env%scale_T)
    CALL section_vals_val_get(mp2_section,"GROUP_SIZE",i_val=mp2_env%mp2_num_proc)
    CALL section_vals_val_get(mp2_section,"DIRECT_CANONICAL%BIG_SEND",l_val=mp2_env%direct_canonical%big_send)
    CALL section_vals_val_get(mp2_section,"ROW_BLOCK",i_val=mp2_env%block_size_row)
    CALL section_vals_val_get(mp2_section,"COL_BLOCK",i_val=mp2_env%block_size_col)
    CALL section_vals_val_get(mp2_section,"CALC_COND_NUM",l_val=mp2_env%calc_PQ_cond_num)

    CALL section_vals_val_get(mp2_section,"INTERACTION_POTENTIAL%POTENTIAL_TYPE",&
                              i_val=mp2_env%potential_parameter%potential_type)
    CALL section_vals_val_get(mp2_section,"INTERACTION_POTENTIAL%TRUNCATION_RADIUS",&
                              r_val=mp2_env%potential_parameter%truncation_radius)
    CALL section_vals_val_get(mp2_section,"INTERACTION_POTENTIAL%TShPSC_DATA",&
                              c_val=mp2_env%potential_parameter%filename)

    CALL section_vals_val_get(mp2_section,"WFC_GPW%EPS_FILTER",&
                              r_val=mp2_env%mp2_gpw%eps_filter)
    CALL section_vals_val_get(mp2_section,"WFC_GPW%EPS_GRID",&
                              r_val=mp2_env%mp2_gpw%eps_grid)
    CALL section_vals_val_get(mp2_section,"WFC_GPW%CUTOFF",&
                              r_val=mp2_env%mp2_gpw%cutoff)
    CALL section_vals_val_get(mp2_section,"WFC_GPW%REL_CUTOFF",&
                              r_val=mp2_env%mp2_gpw%relative_cutoff)
    CALL section_vals_val_get(mp2_section,"WFC_GPW%MULTIPOLE_TWO_CENT_INT",&
                              l_val=mp2_env%mp2_gpw%do_mult_2c)
    CALL section_vals_val_get(mp2_section,"WFC_GPW%PRINT_LEVEL",&
                              i_val=mp2_env%mp2_gpw%print_level)

    CALL section_vals_val_get(mp2_section,"RI_RPA%RPA_NUM_QUAD_POINTS",&
                              i_val=mp2_env%ri_rpa%rpa_num_quad_points)
    CALL section_vals_val_get(mp2_section,"RI_RPA%SIZE_FREQ_INTEG_GROUP",&
                              i_val=mp2_env%ri_rpa%rpa_integ_group_size)
    CALL section_vals_val_get(mp2_section,"RI_RPA%MM_STYLE",&
                              i_val=mp2_env%ri_rpa%mm_style)
    CALL section_vals_val_get(mp2_section,"RI_RPA%MINIMAX_QUADRATURE",&
                              l_val=mp2_env%ri_rpa%minimax_quad)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0",&
                              l_val=mp2_env%ri_rpa%do_ri_g0w0)

    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%CORR_MOS_OCC",&
                              i_val=mp2_env%ri_g0w0%corr_mos_occ)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%CORR_MOS_VIRT",&
                              i_val=mp2_env%ri_g0w0%corr_mos_virt)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%SCALING",&
                              r_val=mp2_env%ri_g0w0%scaling)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%NUMB_POLES",&
                              i_val=mp2_env%ri_g0w0%num_poles)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%OMEGA_MAX_FIT",&
                              r_val=mp2_env%ri_g0w0%omega_max_fit)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%STOP_CRIT",&
                              r_val=mp2_env%ri_g0w0%stop_crit)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%MAX_ITER_FIT",&
                              i_val=mp2_env%ri_g0w0%max_iter_fit)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%CHECK_FIT",&
                              l_val=mp2_env%ri_g0w0%check_fit)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%CROSSING_SEARCH",&
                              i_val=mp2_env%ri_g0w0%crossing_search)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%FERMI_LEVEL_OFFSET",&
                              r_val=mp2_env%ri_g0w0%fermi_level_offset)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%CUTOFF_RADIUS",&
                              r_val=mp2_env%ri_g0w0%cutoff_rad_gw)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%TRUNCATION",&
                              l_val=mp2_env%ri_g0w0%do_truncation)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%EV_SC_ITER",&
                              i_val=mp2_env%ri_g0w0%iter_ev_sc)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%HF_LIKE_EV_START",&
                              l_val=mp2_env%ri_g0w0%hf_like_ev_start)
    CALL section_vals_val_get(mp2_section,"RI_RPA%RI_G0W0%PRINT_GW_DETAILS",&
                              l_val=mp2_env%ri_g0w0%print_gw_details)


    CALL section_vals_val_get(mp2_section,"RI_LAPLACE%QUADRATURE_POINTS",&
                              i_val=mp2_env%ri_laplace%n_quadrature)
    CALL section_vals_val_get(mp2_section,"RI_LAPLACE%SIZE_INTEG_GROUP",&
                              i_val=mp2_env%ri_laplace%integ_group_size)

    CALL section_vals_val_get(mp2_section,"RI_MP2%BLOCK_SIZE",&
                              i_val=mp2_env%ri_mp2%block_size)
    CALL section_vals_val_get(mp2_section,"RI_MP2%EPS_CANONICAL",&
                              r_val=mp2_env%ri_mp2%eps_canonical)
    CALL section_vals_val_get(mp2_section,"RI_MP2%FREE_HFX_BUFFER",&
                              l_val=mp2_env%ri_mp2%free_hfx_buffer)

    CALL section_vals_val_get(mp2_section,"CPHF%MAX_ITER",&
                              i_val=mp2_env%ri_grad%cphf_max_num_iter)
    CALL section_vals_val_get(mp2_section,"CPHF%EPS_CONV",&
                              r_val=mp2_env%ri_grad%cphf_eps_conv)

    CALL section_vals_val_get(mp2_section,"OPT_RI_BASIS%DELTA_I_REL",&
                              r_val=mp2_env%ri_opt_param%DI_rel)
    CALL section_vals_val_get(mp2_section,"OPT_RI_BASIS%DELTA_RI",&
                              r_val=mp2_env%ri_opt_param%DRI)
    CALL section_vals_val_get(mp2_section,"OPT_RI_BASIS%EPS_DERIV",&
                              r_val=mp2_env%ri_opt_param%eps_step)
    CALL section_vals_val_get(mp2_section,"OPT_RI_BASIS%MAX_ITER",&
                              i_val=mp2_env%ri_opt_param%max_num_iter)
    CALL section_vals_val_get(mp2_section,"OPT_RI_BASIS%BASIS_SIZE",&
                              i_val=mp2_env%ri_opt_param%basis_quality)
    NULLIFY(tmplist)
    CALL section_vals_val_get(mp2_section,"OPT_RI_BASIS%NUM_FUNC",&
                              i_vals=tmplist)
    IF(tmplist(1)>0) THEN
      ALLOCATE(mp2_env%ri_opt_param%RI_nset_per_l(0:SIZE(tmplist)-1))
      mp2_env%ri_opt_param%RI_nset_per_l=0
      DO ival=1, SIZE(tmplist)
        mp2_env%ri_opt_param%RI_nset_per_l(ival-1)=tmplist(ival)
      END DO
    END IF

    ! print some info about the MP2 parameters
    unit_nr = cp_print_key_unit_nr(logger,mp2_section,"MP2_INFO",&
         extension=".mp2Log")
    IF ((mp2_env%method .NE. mp2_method_none) .AND. unit_nr>0) THEN
       WRITE(unit_nr,'(T2,A)') ""
       SELECT CASE(mp2_env%method)
       CASE(mp2_method_direct)
         WRITE(unit_nr,'(T2,A)') "MP2| using direct canonical MP2"
       CASE(mp2_method_gpw)
         WRITE(unit_nr,'(T2,A)') "MP2| using MP2 GPW style"
       CASE(ri_mp2_method_gpw)
         WRITE(unit_nr,'(T2,A)') "MP2| using RI-MP2-GPW"
       CASE(ri_rpa_method_gpw)
         WRITE(unit_nr,'(T2,A)') "RI-RPA| using GPW style"
       CASE(ri_mp2_laplace)
         WRITE(unit_nr,'(T2,A)') "RI-SOS-Laplace-MP2| using GPW style"
       CASE(mp2_ri_optimize_basis)
         WRITE(unit_nr,'(T2,A)') "MP2| Optimize RI auxiliary basis"
       CASE DEFAULT
         CPASSERT(.FALSE.)
       END SELECT
       WRITE(unit_nr,'(T2,A)') ""
    ENDIF
    CALL cp_print_key_finished_output(unit_nr,logger,mp2_section,&
                                     "MP2_INFO")

    CALL timestop(handle)

  END SUBROUTINE read_mp2_section

END MODULE mp2_setup


