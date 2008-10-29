!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Subroutine input_torsions changed (DG) 05-Dec-2000
!>      Output formats changed (DG) 05-Dec-2000
!>      JGH (26-01-2002) : force field parameters stored in tables, not in
!>        matrices. Input changed to have parameters labeled by the position
!>        and not atom pairs (triples etc)
!>      Teo (11.2005) : Moved all information on force field  pair_potential to
!>                      a much lighter memory structure
!>      Teo 09.2006   : Splitted all routines force_field I/O in a separate file
!> \author CJM
! *****************************************************************************
MODULE force_fields_input
  USE bibliography,                    ONLY: Tersoff1988,&
                                             Tosi1964a,&
                                             Tosi1964b,&
                                             Yamada2000,&
                                             cite_reference
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line,&
                                             parser_get_object,&
                                             parser_search_string,&
                                             parser_test_next_token
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE force_field_types,               ONLY: charmm_info_type,&
                                             force_field_type,&
                                             gromos_info_type,&
                                             input_info_type
  USE force_fields_util,               ONLY: get_generic_info
  USE input_constants,                 ONLY: do_ff_charmm,&
                                             do_ff_g87,&
                                             do_ff_g96,&
                                             do_ff_undef
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE memory_utilities,                ONLY: reallocate
  USE pair_potential_types,            ONLY: &
       b4_type, bm_type, do_potential_single_allocation, ea_type, &
       eam_pot_type, ft_pot_type, ft_type, gp_type, gw_type, ip_type, &
       ipbv_pot_type, lj_charmm_type, no_potential_single_allocation, &
       pair_potential_p_type, pair_potential_reallocate, &
       potential_single_allocation, tersoff_type, wl_type
  USE shell_potential_types,           ONLY: shell_p_create,&
                                             shell_p_type
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields_input'

  PRIVATE
  PUBLIC :: read_force_field_section,&
            read_force_field_charmm,&
            read_force_field_gromos,&
            read_lj_section,&
            read_wl_section,&
            read_gd_section,&
            read_gp_section,&
            read_chrg_section

CONTAINS

! *****************************************************************************
!> \brief Reads the force_field input section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_force_field_section1 ( ff_section, mm_section, ff_type, para_env, error )
    TYPE(section_vals_type), POINTER         :: ff_section, mm_section
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_force_field_section1', &
      routineP = moduleN//':'//routineN

    INTEGER :: nb4, nbends, nbm, nbmhft, nbonds, nchg, neam, ngd, ngp, nipbv, &
      nlj, nshell, ntersoff, ntors, ntot, nubs, nwl
    LOGICAL                                  :: explicit, failure, &
                                                unique_spline
    REAL(KIND=dp)                            :: min_eps_spline_allowed
    TYPE(input_info_type), POINTER           :: inp_info
    TYPE(section_vals_type), POINTER         :: tmp_section, tmp_section2

    NULLIFY(tmp_section, tmp_section2)
    inp_info => ff_type%inp_info
    failure = .FALSE.
    IF (.NOT.failure) THEN
       CALL section_vals_val_get(ff_section,"PARMTYPE",i_val=ff_type%ff_type,error=error)
       CALL section_vals_val_get(ff_section,"EI_SCALE14",r_val=ff_type%ei_scale14,error=error)
       CALL section_vals_val_get(ff_section,"VDW_SCALE14",r_val=ff_type%vdw_scale14,error=error)
       CALL section_vals_val_get(ff_section,"SPLINE%RCUT_NB",r_val=ff_type%rcut_nb,error=error)
       CALL section_vals_val_get(ff_section,"SPLINE%R0_NB",r_val=ff_type%rlow_nb,error=error)
       CALL section_vals_val_get(ff_section,"SPLINE%EPS_SPLINE",r_val=ff_type%eps_spline,error=error)
       CALL section_vals_val_get(ff_section,"SPLINE%EMAX_SPLINE",r_val=ff_type%emax_spline,error=error)
       CALL section_vals_val_get(ff_section,"SPLINE%EMAX_ACCURACY",r_val=ff_type%max_energy,error=error)
       CPPostcondition(ff_type%max_energy<=ff_type%emax_spline,cp_failure_level,routineP,error,failure)
       ! Read the parameter file name only if the force field type requires it..
       SELECT CASE(ff_type%ff_type)
       CASE(do_ff_charmm,do_ff_g96,do_ff_g87)
          CALL section_vals_val_get(ff_section,"PARM_FILE_NAME",c_val=ff_type%ff_file_name,error=error)
          CALL cp_assert(TRIM(ff_type%ff_file_name)/="",cp_failure_level,cp_assertion_failed,routineP,&
               "Force Field Parameter's filename is empty! Please check your input file. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       CASE(do_ff_undef)
          ! Do Nothing
       CASE DEFAULT
          CALL stop_program ("force_fields_input","force field type not implemented")
       END SELECT 
       ! Numerical Accuracy:
       ! the factors here should depend on the energy and on the shape of each potential mapped
       ! with splines. this would make everything un-necessarily complicated. Let's just be safe 
       ! and assume that we cannot achieve an accuracy on the spline 2 orders of magnitude more 
       ! than the smallest representable number (taking into account also the max_energy for the 
       ! spline generation       
       min_eps_spline_allowed = 20.0_dp*MAX(ff_type%max_energy,10.0_dp)*EPSILON(0.0_dp)
       IF (ff_type%eps_spline < min_eps_spline_allowed) THEN
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Requested spline accuracy ("//TRIM(cp_to_string(ff_type%eps_spline))//" ) "//&
               "is smaller than the minimum value allowed ("//TRIM(cp_to_string(min_eps_spline_allowed))//&
               " ) with the present machine precision ("//TRIM(cp_to_string(EPSILON(0.0_dp)))//" ). "//&
               "New EPS_SPLINE value ("//TRIM(cp_to_string(min_eps_spline_allowed))//" ). "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ff_type%eps_spline = min_eps_spline_allowed
       END IF
       CALL section_vals_val_get(ff_section,"SCALE_CUTOFF",l_val=ff_type%scale_cutoff,error=error)
       CALL section_vals_val_get(ff_section,"SPLINE%UNIQUE_SPLINE",l_val=unique_spline,error=error)
       ! Single spline
       potential_single_allocation = no_potential_single_allocation
       IF (unique_spline) potential_single_allocation = do_potential_single_allocation

       CALL section_vals_val_get(ff_section,"MULTIPLE_POTENTIAL",l_val=ff_type%multiple_potential,error=error)
       CALL section_vals_val_get(ff_section,"NONBONDED",l_val=ff_type%nonbonded,error=error)
       tmp_section  => section_vals_get_subs_vals(ff_section,"NONBONDED",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,error=error)
       IF (explicit.AND.ff_type%nonbonded) THEN
          tmp_section2 => section_vals_get_subs_vals(tmp_section,"LENNARD-JONES",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nlj,error=error)
          ntot = 0
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nlj,lj_charmm=.TRUE.,error=error)
             CALL read_lj_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"WILLIAMS",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nwl,error=error)
          ntot=nlj
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nwl,williams=.TRUE.,error=error)
             CALL read_wl_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"EAM",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=neam,error=error)
          ntot=nlj+nwl
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+neam,eam=.TRUE.,error=error)
             CALL read_eam_section(inp_info%nonbonded,tmp_section2,ntot,para_env,mm_section,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"GOODWIN",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=ngd,error=error)
          ntot=nlj+nwl+neam
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+ngd,goodwin=.TRUE.,error=error)
             CALL read_gd_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"IPBV",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nipbv,error=error)
          ntot=nlj+nwl+neam+ngd
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nipbv,ipbv=.TRUE.,error=error)
             CALL read_ipbv_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"BMHFT",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nbmhft,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nbmhft,bmhft=.TRUE.,error=error)
             CALL read_bmhft_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"BUCK4RANGES",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nb4,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nb4,buck4r=.TRUE.,error=error)
             CALL read_b4_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"BUCKMORSE",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nbm,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nb4
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nbm,buckmo=.TRUE.,error=error)
             CALL read_bm_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"GENPOT",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=ngp,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nb4+nbm
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+ngp,gp=.TRUE.,error=error)
             CALL read_gp_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF
          tmp_section2 => section_vals_get_subs_vals(tmp_section,"TERSOFF",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=ntersoff,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nb4+nbm+ngp
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+ntersoff,tersoff=.TRUE.,error=error)
             CALL read_tersoff_section(inp_info%nonbonded,tmp_section2,ntot,tmp_section2,error)
          END IF
       END IF

       CALL section_vals_val_get(ff_section,"NONBONDED14",l_val=ff_type%nonbonded14,error=error)
       tmp_section => section_vals_get_subs_vals(ff_section,"NONBONDED14",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,error=error)
       IF (explicit.AND.ff_type%nonbonded14) THEN
          tmp_section2 => section_vals_get_subs_vals(tmp_section,"LENNARD-JONES",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nlj,error=error)
          ntot = 0
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded14,1,ntot+nlj,lj_charmm=.TRUE.,error=error)
             CALL read_lj_section(inp_info%nonbonded14,tmp_section2,ntot,error)
          END IF
          tmp_section2 => section_vals_get_subs_vals(tmp_section,"WILLIAMS",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nwl,error=error)
          ntot=nlj
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded14,1,ntot+nwl, williams=.TRUE.,error=error)
             CALL read_wl_section(inp_info%nonbonded14,tmp_section2,ntot,error)
          END IF
          tmp_section2 => section_vals_get_subs_vals(tmp_section,"GOODWIN",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=ngd,error=error)
          ntot=nlj+nwl
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded14,1,ntot+ngd, goodwin=.TRUE.,error=error)
             CALL read_gd_section(inp_info%nonbonded14,tmp_section2,ntot,error)
          END IF
          tmp_section2 => section_vals_get_subs_vals(tmp_section,"GENPOT",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=ngp,error=error)
          ntot=nlj+nwl+ngd
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded14,1,ntot+ngp,gp=.TRUE.,error=error)
             CALL read_gp_section(inp_info%nonbonded14,tmp_section2,ntot,error)
          END IF
       END IF

       tmp_section => section_vals_get_subs_vals(ff_section,"CHARGE",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nchg,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%charge_atm,1,nchg)
          CALL reallocate(inp_info%charge,1,nchg)
          CALL read_chrg_section(inp_info%charge_atm,inp_info%charge,tmp_section,ntot,error)
       END IF
       tmp_section => section_vals_get_subs_vals(ff_section,"SHELL",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nshell,error=error)
       IF (explicit) THEN
          ntot=0
          CALL shell_p_create(inp_info%shell_list,nshell,error)
          CALL read_shell_section(inp_info%shell_list,tmp_section,ntot,error)
       END IF

       tmp_section => section_vals_get_subs_vals(ff_section,"BOND",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nbonds,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%bond_kind,1,nbonds)
          CALL reallocate(inp_info%bond_a,1,nbonds)
          CALL reallocate(inp_info%bond_b,1,nbonds)
          CALL reallocate(inp_info%bond_k,1,3,1,nbonds)
          CALL reallocate(inp_info%bond_r0,1,nbonds)
          CALL reallocate(inp_info%bond_cs,1,nbonds)
          CALL read_bonds_section(inp_info%bond_kind, inp_info%bond_a, inp_info%bond_b, inp_info%bond_k,&
               inp_info%bond_r0, inp_info%bond_cs, tmp_section, ntot, error)
       END IF
       tmp_section => section_vals_get_subs_vals(ff_section,"BEND",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nbends,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%bend_kind,1,nbends)
          CALL reallocate(inp_info%bend_a,1,nbends)
          CALL reallocate(inp_info%bend_b,1,nbends)
          CALL reallocate(inp_info%bend_c,1,nbends)
          CALL reallocate(inp_info%bend_k,1,nbends)
          CALL reallocate(inp_info%bend_theta0,1,nbends)
          CALL reallocate(inp_info%bend_cb,1,nbends)
          CALL reallocate(inp_info%bend_r012,1,nbends)
          CALL reallocate(inp_info%bend_r032,1,nbends)
          CALL reallocate(inp_info%bend_kbs12,1,nbends)
          CALL reallocate(inp_info%bend_kbs32,1,nbends)
          CALL reallocate(inp_info%bend_kss,1,nbends)
          CALL read_bends_section(inp_info%bend_kind,inp_info%bend_a, inp_info%bend_b, inp_info%bend_c,&
               inp_info%bend_k, inp_info%bend_theta0, inp_info%bend_cb,&
               inp_info%bend_r012,inp_info%bend_r032,inp_info%bend_kbs12,inp_info%bend_kbs32,inp_info%bend_kss,&
               tmp_section, ntot, error)
       END IF
       tmp_section => section_vals_get_subs_vals(ff_section,"BEND",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nubs,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%ub_kind,1,nubs)
          CALL reallocate(inp_info%ub_a,1,nubs)
          CALL reallocate(inp_info%ub_b,1,nubs)
          CALL reallocate(inp_info%ub_c,1,nubs)
          CALL reallocate(inp_info%ub_k,1,3,1,nubs)
          CALL reallocate(inp_info%ub_r0,1,nubs)
          CALL read_ubs_section(inp_info%ub_kind,inp_info%ub_a, inp_info%ub_b, inp_info%ub_c,&
               inp_info%ub_k, inp_info%ub_r0, tmp_section, ntot, error)
       END IF
       tmp_section => section_vals_get_subs_vals(ff_section,"TORSION",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=ntors,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%torsion_kind,1,ntors)
          CALL reallocate(inp_info%torsion_a,1,ntors)
          CALL reallocate(inp_info%torsion_b,1,ntors)
          CALL reallocate(inp_info%torsion_c,1,ntors)
          CALL reallocate(inp_info%torsion_d,1,ntors)
          CALL reallocate(inp_info%torsion_k,1,ntors)
          CALL reallocate(inp_info%torsion_m,1,ntors)
          CALL reallocate(inp_info%torsion_cosphi0,1,ntors)
          CALL read_torsions_section(inp_info%torsion_kind,inp_info%torsion_a, inp_info%torsion_b, inp_info%torsion_c,&
               inp_info%torsion_d,inp_info%torsion_k,inp_info%torsion_cosphi0, inp_info%torsion_m,&
               tmp_section, ntot, error )
       END IF
    END IF

  END SUBROUTINE read_force_field_section1

! *****************************************************************************
!> \brief Set up of the IPBV force fields
!> \author teo
! *****************************************************************************
  SUBROUTINE set_IPBV_ff(at1, at2, ipbv )
    CHARACTER(LEN=*), INTENT(IN)             :: at1, at2
    TYPE(ipbv_pot_type), POINTER             :: ipbv

    IF  ( ( at1(1:1) == 'O' ) .AND. ( at2 (1:1) == 'O' ) ) THEN
       ipbv%rcore = 0.9_dp                  ! a.u.
       ipbv%m = -1.2226442563398141E+11_dp   ! Kelvin/a.u.
       ipbv%b =  1.1791292385486696E+11_dp  ! Hartree

       ! Hartree*a.u.^2
       ipbv% a ( 2 ) = 4.786380682394_dp
       ipbv% a ( 3 ) = -1543.407053545_dp
       ipbv% a ( 4 ) = 88783.31188529_dp
       ipbv% a ( 5 ) = -2361200.155376_dp
       ipbv% a ( 6 ) = 35940504.84679_dp
       ipbv% a ( 7 ) = -339762743.6358_dp
       ipbv% a ( 8 ) = 2043874926.466_dp
       ipbv% a ( 9 ) = -7654856796.383_dp
       ipbv% a ( 10 ) = 16195251405.65_dp
       ipbv% a ( 11 ) = -13140392992.18_dp
       ipbv% a ( 12 ) = -9285572894.245_dp
       ipbv% a ( 13 ) = 8756947519.029_dp
       ipbv% a ( 14 ) = 15793297761.67_dp
       ipbv% a ( 15 ) = 12917180227.21_dp
    ELSEIF(((at1(1:1) == 'O').AND.(at2 (1:1) == 'H')).OR.&
         ((at1(1:1) == 'H').AND.(at2 (1:1) == 'O')))  THEN
       ipbv%rcore = 2.95_dp       ! a.u.

       ipbv%m = -0.004025691139759147_dp  ! Hartree/a.u.
       ipbv%b =  -2.193731138097428_dp    ! Hartree
       ! Hartree*a.u.^2
       ipbv% a ( 2 ) = -195.7716013277_dp
       ipbv% a ( 3 ) = 15343.78613395_dp
       ipbv% a ( 4 ) = -530864.4586516_dp
       ipbv% a ( 5 ) = 10707934.39058_dp
       ipbv% a ( 6 ) = -140099704.7890_dp
       ipbv% a ( 7 ) = 1250943273.785_dp
       ipbv% a ( 8 ) = -7795458330.676_dp
       ipbv% a ( 9 ) = 33955897217.31_dp
       ipbv% a ( 10 ) = -101135640744.0_dp
       ipbv% a ( 11 ) = 193107995718.7_dp
       ipbv% a ( 12 ) = -193440560940.0_dp
       ipbv% a ( 13 ) =-4224406093.918E0_dp
       ipbv% a ( 14 ) = 217192386506.5E0_dp
       ipbv% a ( 15 ) = -157581228915.5_dp
    ELSEIF ( ( at1(1:1) == 'H' ) .AND. ( at2 (1:1) == 'H' ) ) THEN
       ipbv%rcore = 3.165_dp     ! a.u.
       ipbv%m = 0.002639704108787555_dp  ! Hartree/a.u.
       ipbv%b = -0.2735482611857583_dp ! Hartree
       ! Hartree*a.u.^2
       ipbv% a ( 2 ) = -26.29456010782_dp
       ipbv% a ( 3 ) =  2373.352548248_dp
       ipbv% a ( 4 ) = -93880.43551360_dp
       ipbv% a ( 5 ) =  2154624.884809_dp
       ipbv% a ( 6 ) = -31965151.34955_dp
       ipbv% a ( 7 ) = 322781785.3278_dp
       ipbv% a ( 8 ) = -2271097368.668_dp
       ipbv% a ( 9 ) = 11169163192.90_dp
       ipbv% a ( 10 ) =-37684457778.47_dp
       ipbv% a ( 11 ) = 82562104256.03_dp
       ipbv% a ( 12 ) = -100510435213.4_dp
       ipbv% a ( 13 ) = 24570342714.65E0_dp
       ipbv% a ( 14 ) = 88766181532.94E0_dp
       ipbv% a ( 15 ) = -79705131323.98_dp
    ELSE
       CALL stop_program ("input_nonbonded","IPBV only for WATER")
    ENDIF
  END SUBROUTINE set_IPBV_ff

! *****************************************************************************
!> \brief Set up of the BMHFT force fields
!> \author teo
! *****************************************************************************
  SUBROUTINE set_BMHFT_ff(at1, at2, ft, error)
    CHARACTER(LEN=*), INTENT(IN)             :: at1, at2
    TYPE(ft_pot_type), POINTER               :: ft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_BMHFT_ff', &
      routineP = moduleN//':'//routineN

    ft%b = cp_unit_to_cp2k(3.1545_dp,"angstrom^-1",error=error)
    IF  ( ( at1(1:2) == 'NA' ) .AND. ( at2 (1:2) == 'NA' ) ) THEN
       ft%a = cp_unit_to_cp2k(424.097_dp,"eV",error=error)
       ft%c = cp_unit_to_cp2k(1.05_dp,"eV*angstrom^6",error=error)
       ft%d = cp_unit_to_cp2k(0.499_dp,"eV*angstrom^8",error=error)
    ELSEIF(((at1(1:2) == 'NA').AND.(at2(1:2) == 'CL')).OR.&
         ((at1(1:2) == 'CL').AND.(at2(1:2) == 'NA'))) THEN
       ft%a = cp_unit_to_cp2k(1256.31_dp,"eV",error=error)
       ft%c = cp_unit_to_cp2k(7.00_dp,"eV*angstrom^6",error=error)
       ft%d = cp_unit_to_cp2k(8.676_dp,"eV*angstrom^8",error=error)
    ELSEIF ( ( at1(1:2) == 'CL' ) .AND. ( at2 (1:2) == 'CL' ) ) THEN
       ft%a = cp_unit_to_cp2k(3488.998_dp,"eV",error=error)
       ft%c = cp_unit_to_cp2k(72.50_dp,"eV*angstrom^6",error=error)
       ft%d = cp_unit_to_cp2k(145.427_dp,"eV*angstrom^8",error=error)
    ELSE
       CALL stop_program ("input_nonbonded","BMHFT only for NaCl")
    ENDIF

  END SUBROUTINE set_BMHFT_ff

! *****************************************************************************
!> \brief Reads the EAM section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_eam_section(nonbonded, section, start, para_env, mm_section, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_eam_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)

       nonbonded%pot(start+isec)%pot %type = ea_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       CALL section_vals_val_get(section,"PARM_FILE_NAME",i_rep_section=isec,&
            c_val=nonbonded%pot(start+isec)%pot%set(1)%eam%eam_file_name,error=error)
       CALL read_eam_data ( nonbonded%pot(start+isec)%pot%set(1)%eam, para_env,  mm_section,error)
       nonbonded%pot(start+isec)%pot%rcutsq  = nonbonded%pot(start+isec)%pot%set(1)%eam%acutal**2
    END DO
  END SUBROUTINE read_eam_section

! *****************************************************************************
!> \brief Reads the LJ section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_lj_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_lj_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items, n_rep
    REAL(KIND=dp)                            :: epsilon, rcut, sigma

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       CALL section_vals_val_get(section,"EPSILON",i_rep_section=isec,r_val=epsilon,error=error)
       CALL section_vals_val_get(section,"SIGMA",i_rep_section=isec,r_val=sigma,error=error)
       CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)

       nonbonded%pot(start+isec)%pot %type = lj_charmm_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       nonbonded%pot(start+isec)%pot %set(1)%lj%epsilon = epsilon
       nonbonded%pot(start+isec)%pot %set(1)%lj%sigma6  = sigma**6
       nonbonded%pot(start+isec)%pot %set(1)%lj%sigma12 = sigma**12
       nonbonded%pot(start+isec)%pot %rcutsq = rcut*rcut
       ! 
       CALL section_vals_val_get ( section, "RMIN", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,error=error)
       CALL section_vals_val_get ( section, "RMAX", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,error=error)
    END DO
  END SUBROUTINE read_lj_section

! *****************************************************************************
!> \brief Reads the WILLIAMS section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_wl_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_wl_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items, n_rep
    REAL(KIND=dp)                            :: a, b, c, rcut

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       CALL section_vals_val_get(section,"A",i_rep_section=isec,r_val=a,error=error)
       CALL section_vals_val_get(section,"B",i_rep_section=isec,r_val=b,error=error)
       CALL section_vals_val_get(section,"C",i_rep_section=isec,r_val=c,error=error)
       CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)

       nonbonded%pot(start+isec)%pot %type = wl_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       nonbonded%pot(start+isec)%pot %set(1)%willis%a = a
       nonbonded%pot(start+isec)%pot %set(1)%willis%b = b
       nonbonded%pot(start+isec)%pot %set(1)%willis%c = c
       nonbonded%pot(start+isec)%pot %rcutsq = rcut*rcut
       ! 
       CALL section_vals_val_get ( section, "RMIN", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,error=error)
       CALL section_vals_val_get ( section, "RMAX", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,error=error)
    END DO
  END SUBROUTINE read_wl_section

! *****************************************************************************
!> \brief Reads the GOODWIN section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_gd_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_gd_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, m, mc, n_items, n_rep
    REAL(KIND=dp)                            :: d, dc, rcut, vr0

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       CALL section_vals_val_get(section,"VR0",i_rep_section=isec,r_val=vr0,error=error)
       CALL section_vals_val_get(section,"D",i_rep_section=isec,r_val=d,error=error)
       CALL section_vals_val_get(section,"DC",i_rep_section=isec,r_val=dc,error=error)
       CALL section_vals_val_get(section,"M",i_rep_section=isec,i_val=m,error=error)
       CALL section_vals_val_get(section,"MC",i_rep_section=isec,i_val=mc,error=error)
       CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)

       nonbonded%pot(start+isec)%pot %type = gw_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       nonbonded%pot(start+isec)%pot %set(1)%goodwin%vr0 = vr0
       nonbonded%pot(start+isec)%pot %set(1)%goodwin%d = d
       nonbonded%pot(start+isec)%pot %set(1)%goodwin%dc = dc
       nonbonded%pot(start+isec)%pot %set(1)%goodwin%m = m
       nonbonded%pot(start+isec)%pot %set(1)%goodwin%mc = mc
       nonbonded%pot(start+isec)%pot %rcutsq = rcut*rcut
       ! 
       CALL section_vals_val_get ( section, "RMIN", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,error=error)
       CALL section_vals_val_get ( section, "RMAX", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,error=error)
    END DO
  END SUBROUTINE read_gd_section

! *****************************************************************************
!> \brief Reads the IPBV section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_ipbv_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_ipbv_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items, n_rep
    REAL(KIND=dp)                            :: rcut

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       nonbonded%pot(start+isec)%pot %type = ip_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       CALL set_IPBV_ff(nonbonded%pot(start+isec)%pot %at1, nonbonded%pot(start+isec)%pot %at2,&
            nonbonded%pot(start+isec)%pot%set(1)%ipbv )
       CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)
       nonbonded%pot(start+isec)%pot%rcutsq     = rcut**2
       ! 
       CALL section_vals_val_get ( section, "RMIN", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,error=error)
       CALL section_vals_val_get ( section, "RMAX", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,error=error)
    END DO
  END SUBROUTINE read_ipbv_section

! *****************************************************************************
!> \brief Reads the BMHFT section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_bmhft_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_bmhft_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(2)                           :: map_atoms
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: i, isec, n_items, n_rep
    REAL(KIND=dp)                            :: rcut

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL cite_reference(Tosi1964a)
       CALL cite_reference(Tosi1964b)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       nonbonded%pot(start+isec)%pot %type = ft_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       
       CALL section_vals_val_get(section,"A",i_rep_section=isec,n_rep_val=i,error=error)
       IF (i==1) THEN
          CALL section_vals_val_get(section,"A",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ft%a,error=error)
          CALL section_vals_val_get(section,"B",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ft%b,error=error)
          CALL section_vals_val_get(section,"C",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ft%c,error=error)
          CALL section_vals_val_get(section,"D",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ft%d,error=error)
       ELSE
          CALL section_vals_val_get(section,"MAP_ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
          map_atoms = atm_names
          CALL uppercase(map_atoms(1))
          CALL uppercase(map_atoms(2))
          CALL set_BMHFT_ff(map_atoms(1), map_atoms(2), nonbonded%pot(start+isec)%pot%set(1)%ft,&
                            error=error)
       END IF
       CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)
       nonbonded%pot(start+isec)%pot%rcutsq     = rcut**2
       ! 
       CALL section_vals_val_get ( section, "RMIN", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,error=error)
       CALL section_vals_val_get ( section, "RMAX", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,error=error)
    END DO
  END SUBROUTINE read_bmhft_section

! *****************************************************************************
!> \brief Reads the Buckingham 4 Ranges potential section
!> \author MI
! *****************************************************************************
  SUBROUTINE read_b4_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_b4_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: i, ir, isec, istat, n_items, &
                                                n_rep, np1, np2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, b, c, r1, r2, r3, rcut
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coeff1, coeff2, list

    failure = .FALSE.
    NULLIFY(coeff1,coeff2)
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       CALL section_vals_val_get(section,"A",i_rep_section=isec,r_val=a,error=error)
       CALL section_vals_val_get(section,"B",i_rep_section=isec,r_val=b,error=error)
       CALL section_vals_val_get(section,"C",i_rep_section=isec,r_val=c,error=error)
       CALL section_vals_val_get(section,"R1",i_rep_section=isec,r_val=r1,error=error)
       CALL section_vals_val_get(section,"R2",i_rep_section=isec,r_val=r2,error=error)
       CALL section_vals_val_get(section,"R3",i_rep_section=isec,r_val=r3,error=error)
       CALL section_vals_val_get(section,"POLY1",n_rep_val=n_rep,error=error)
       np1 = 0
       DO ir = 1,n_rep
         NULLIFY(list)
         CALL section_vals_val_get(section,"POLY1",i_rep_val=ir,r_vals=list,error=error)
         IF(ASSOCIATED(list)) THEN
           CALL reallocate(coeff1,0,np1+SIZE(list)-1)
           DO i = 1, SIZE(list)
             coeff1(i+np1-1) = list(i)
           END DO  ! i
           np1 = np1 + SIZE(list)
         END IF
       END DO  ! ir
       CALL section_vals_val_get(section,"POLY2",n_rep_val=n_rep,error=error)
       np2 = 0
       DO ir = 1,n_rep
         NULLIFY(list)
         CALL section_vals_val_get(section,"POLY2",i_rep_val=ir,r_vals=list,error=error)
         IF(ASSOCIATED(list)) THEN
           CALL reallocate(coeff2,0,np2+SIZE(list)-1)
           DO i = 1, SIZE(list)
             coeff2(i+np2-1) = list(i)
           END DO  ! i
           np2 = np2 + SIZE(list)
         END IF
       END DO  ! ir
       CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)

       nonbonded%pot(start+isec)%pot %type = b4_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%a = a
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%b = b
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%c = c
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%r1 = r1
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%r2 = r2
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%r3 = r3
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%npoly1 = np1-1
       IF( np1-1 > 10) THEN
         CALL stop_program ("input_nonbonded","BUCK4r: too high polynomial order (>0)")
       END IF
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%poly1(0:np1-1) = coeff1(0:np1-1)
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%npoly2 = np2-1
       IF( np2-1 > 10) THEN
         CALL stop_program ("input_nonbonded","BUCK4r: too high polynomial order (>0)")
       END IF
       nonbonded%pot(start+isec)%pot %set(1)%buck4r%poly2(0:np2-1) = coeff2(0:np2-1)
       nonbonded%pot(start+isec)%pot %rcutsq = rcut*rcut
       DEALLOCATE(coeff1,coeff2,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ! 
       CALL section_vals_val_get ( section, "RMIN", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,error=error)
       CALL section_vals_val_get ( section, "RMAX", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,error=error)
    END DO
  END SUBROUTINE read_b4_section

! *****************************************************************************
!> \brief Reads the GENPOT - generic potential section
!> \author Teodoro Laino - 10.2006
! *****************************************************************************
  SUBROUTINE read_gp_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_gp_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items, n_rep
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rcut
    REAL(KIND=dp), DIMENSION(:), POINTER     :: values

    failure = .FALSE.
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       NULLIFY(atm_names)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)
       nonbonded%pot(start+isec)%pot%type = gp_type
       nonbonded%pot(start+isec)%pot%at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot%at2  = atm_names(2)
       nonbonded%pot(start+isec)%pot%rcutsq = rcut*rcut
       CALL uppercase(nonbonded%pot(start+isec)%pot%at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot%at2)
       ! Parse the genpot info
       CALL get_generic_info(section,"FUNCTION",nonbonded%pot(start+isec)%pot%set(1)%gp%potential,&
            nonbonded%pot(start+isec)%pot%set(1)%gp%parameters,nonbonded%pot(start+isec)%pot%set(1)%gp%values,&
            size_variables=1,i_rep_sec=isec,error=error)
       nonbonded%pot(start+isec)%pot%set(1)%gp%variables = nonbonded%pot(start+isec)%pot%set(1)%gp%parameters(1)
       ! 
       CALL section_vals_val_get(section, "RMIN", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,error=error)
       CALL section_vals_val_get(section, "RMAX", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,error=error)
    END DO
  END SUBROUTINE read_gp_section

! *****************************************************************************
!> \brief Reads the tersoff section
!> \author ikuo
! *****************************************************************************
  SUBROUTINE read_tersoff_section(nonbonded, section, start, tersoff_section,error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(section_vals_type), POINTER         :: tersoff_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_tersoff_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items, n_rep
    REAL(KIND=dp)                            :: rcut, rcutsq

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL cite_reference(Tersoff1988)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)

       nonbonded%pot(start+isec)%pot %type = tersoff_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)

       CALL section_vals_val_get(tersoff_section,"A",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%A,error=error)
       CALL section_vals_val_get(tersoff_section,"B",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%B,error=error)
       CALL section_vals_val_get(tersoff_section,"lambda1",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%lambda1,error=error)
       CALL section_vals_val_get(tersoff_section,"lambda2",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%lambda2,error=error)
       CALL section_vals_val_get(tersoff_section,"alpha",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%alpha,error=error)
       CALL section_vals_val_get(tersoff_section,"beta",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%beta,error=error)
       CALL section_vals_val_get(tersoff_section,"n",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%n,error=error)
       CALL section_vals_val_get(tersoff_section,"c",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%c,error=error)
       CALL section_vals_val_get(tersoff_section,"d",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%d,error=error)
       CALL section_vals_val_get(tersoff_section,"h",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%h,error=error)
       CALL section_vals_val_get(tersoff_section,"lambda3",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%lambda3,error=error)
       CALL section_vals_val_get(tersoff_section,"bigR",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%bigR,error=error)
       CALL section_vals_val_get(tersoff_section,"bigD",r_val=nonbonded%pot(start+isec)%pot%set(1)%tersoff%bigD,error=error)

       rcutsq = (nonbonded%pot(start+isec)%pot%set(1)%tersoff%bigR+&
                 nonbonded%pot(start+isec)%pot%set(1)%tersoff%bigD)**2
       nonbonded%pot(start+isec)%pot%set(1)%tersoff%rcutsq =  rcutsq
       nonbonded%pot(start+isec)%pot%rcutsq = rcutsq

       ! In case it is defined override the standard specification of RCUT
       CALL section_vals_val_get ( tersoff_section, "RCUT", n_rep_val=n_rep, error=error)
       IF(n_rep==1) THEN 
          CALL section_vals_val_get(tersoff_section,"RCUT",r_val=rcut,error=error)
          nonbonded%pot(start+isec)%pot %rcutsq = rcut**2
       END IF
    END DO
  END SUBROUTINE read_tersoff_section

! *****************************************************************************
!> \brief Reads the Buckingham plus Morse potential section
!> \author MI
! *****************************************************************************
  SUBROUTINE read_bm_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_bm_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items, n_rep
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a1, a2, b1, b2, beta, c, d, &
                                                f0, r0, rcut

    failure = .FALSE.
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL cite_reference(Yamada2000)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       CALL section_vals_val_get(section,"F0",i_rep_section=isec,r_val=f0,error=error)
       CALL section_vals_val_get(section,"A1",i_rep_section=isec,r_val=a1,error=error)
       CALL section_vals_val_get(section,"A2",i_rep_section=isec,r_val=a2,error=error)
       CALL section_vals_val_get(section,"B1",i_rep_section=isec,r_val=b1,error=error)
       CALL section_vals_val_get(section,"B2",i_rep_section=isec,r_val=b2,error=error)
       CALL section_vals_val_get(section,"C",i_rep_section=isec,r_val=c,error=error)
       CALL section_vals_val_get(section,"D",i_rep_section=isec,r_val=d,error=error)
       CALL section_vals_val_get(section,"R0",i_rep_section=isec,r_val=r0,error=error)
       CALL section_vals_val_get(section,"Beta",i_rep_section=isec,r_val=beta,error=error)
       CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)

       nonbonded%pot(start+isec)%pot %type = bm_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%f0 = f0
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%a1 = a1
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%a2 = a2
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%b1 = b1
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%b2 = b2
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%c = c
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%d = d
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%r0 = r0
       nonbonded%pot(start+isec)%pot %set(1)%buckmo%beta = beta
       nonbonded%pot(start+isec)%pot %rcutsq = rcut*rcut
       ! 
       CALL section_vals_val_get ( section, "RMIN", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,error=error)
       CALL section_vals_val_get ( section, "RMAX", i_rep_section=isec, n_rep_val=n_rep, error=error)
       IF(n_rep==1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
            r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,error=error)
    END DO
  END SUBROUTINE read_bm_section

! *****************************************************************************
!> \brief Reads the CHARGE section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_chrg_section(charge_atm, charge, section, start, error)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: charge_atm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charge
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_chrg_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atm_name
    INTEGER                                  :: isec, n_items

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOM",i_rep_section=isec,c_val=atm_name,&
            error=error)
       charge_atm(start+isec) = atm_name
       CALL uppercase(charge_atm(start+isec))
       CALL section_vals_val_get(section,"CHARGE",i_rep_section=isec,r_val=charge(start+isec),&
            error=error)
    END DO
  END SUBROUTINE read_chrg_section

! *****************************************************************************
!> \brief Reads the SHELL section
!> \author MI
! *****************************************************************************
  SUBROUTINE read_shell_section(shell_list, section, start, error)
    TYPE(shell_p_type), DIMENSION(:), &
      POINTER                                :: shell_list
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_shell_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atm_name
    INTEGER                                  :: i_rep, n_rep
    REAL(dp)                                 :: ccharge, cutoff, maxdist, &
                                                mfrac, scharge, spring

    CALL section_vals_get(section,n_repetition=n_rep,error=error)
    DO i_rep = 1, n_rep
       CALL section_vals_val_get(section,"_SECTION_PARAMETERS_",&
            c_val=atm_name,i_rep_section=i_rep,error=error)
       CALL uppercase(atm_name)
       shell_list(start+i_rep)%atm_name = atm_name

       CALL section_vals_val_get(section,"CORE_CHARGE",i_rep_section=i_rep,r_val=ccharge,&
            error=error)
       shell_list(start+i_rep)%shell%charge_core=ccharge
       CALL section_vals_val_get(section,"SHELL_CHARGE",i_rep_section=i_rep,r_val=scharge,&
            error=error)
       shell_list(start+i_rep)%shell%charge_shell=scharge
       CALL section_vals_val_get(section,"MASS_FRACTION",i_rep_section=i_rep,r_val=mfrac,&
            error=error)
       shell_list(start+i_rep)%shell%massfrac=mfrac
       CALL section_vals_val_get(section,"SPRING",i_rep_section=i_rep,r_val=spring,&
            error=error)
       shell_list(start+i_rep)%shell%k_spring=spring
       CALL section_vals_val_get(section,"MAX_DISTANCE",i_rep_section=i_rep,r_val=maxdist,&
            error=error)
       shell_list(start+i_rep)%shell%max_dist=maxdist
       CALL section_vals_val_get(section,"SHELL_CUTOFF",i_rep_section=i_rep,r_val=cutoff,&
            error=error)
       shell_list(start+i_rep)%shell%shell_cutoff=cutoff

    END DO

  END SUBROUTINE read_shell_section

! *****************************************************************************
!> \brief Reads the BONDS section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_bonds_section(bond_kind,bond_a, bond_b, bond_k, bond_r0, bond_cs, section, start, error)
    INTEGER, DIMENSION(:), POINTER           :: bond_kind
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: bond_a, bond_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: bond_k
    REAL(KIND=dp), DIMENSION(:), POINTER     :: bond_r0, bond_cs
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_bonds_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, k, n_items
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Kvals

    failure = .FALSE.
    NULLIFY(Kvals, atm_names)
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"KIND",i_rep_section=isec,i_val=bond_kind(start+isec),error=error)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       bond_a(start+isec) = atm_names(1)
       bond_b(start+isec) = atm_names(2)
       CALL uppercase(bond_a(start+isec))
       CALL uppercase(bond_b(start+isec))
       CALL section_vals_val_get(section,"K",i_rep_section=isec,r_vals=Kvals,error=error)
       CPPostcondition(SIZE(Kvals) <= 3, cp_failure_level, routineP, error, failure)
       bond_k(:,start+isec) = 0.0_dp
       DO k=1,SIZE(Kvals)
          bond_k(k,start+isec) = Kvals(k)
       END DO
       CALL section_vals_val_get(section,"R0",i_rep_section=isec,r_val=bond_r0(start+isec),error=error)
       CALL section_vals_val_get(section,"CS",i_rep_section=isec,r_val=bond_cs(start+isec),error=error)
    END DO
  END SUBROUTINE read_bonds_section

! *****************************************************************************
!> \brief Reads the BENDS section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_bends_section(bend_kind,bend_a, bend_b, bend_c, bend_k, bend_theta0, bend_cb, &
                                bend_r012,bend_r032,bend_kbs12,bend_kbs32,bend_kss, section, start, error)
    INTEGER, DIMENSION(:), POINTER           :: bend_kind
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: bend_a, bend_b, bend_c
    REAL(KIND=dp), DIMENSION(:), POINTER     :: bend_k, bend_theta0, bend_cb, &
                                                bend_r012, bend_r032, &
                                                bend_kbs12, bend_kbs32, &
                                                bend_kss
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_bends_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Kvals

    failure = .FALSE.
    NULLIFY(Kvals, atm_names)
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"KIND",i_rep_section=isec,i_val=bend_kind(start+isec),error=error)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       bend_a(start+isec) = atm_names(1)
       bend_b(start+isec) = atm_names(2)
       bend_c(start+isec) = atm_names(3)
       CALL uppercase(bend_a(start+isec))
       CALL uppercase(bend_b(start+isec))
       CALL uppercase(bend_c(start+isec))
       CALL section_vals_val_get(section,"K",i_rep_section=isec,r_vals=Kvals,error=error)
       CPPostcondition(SIZE(Kvals) == 1, cp_failure_level, routineP, error, failure)
       bend_k(start+isec) = Kvals(1)
       CALL section_vals_val_get(section,"THETA0",i_rep_section=isec,r_val=bend_theta0(start+isec),error=error)
       CALL section_vals_val_get(section,"CB",i_rep_section=isec,r_val=bend_cb(start+isec),error=error)
       CALL section_vals_val_get(section,"R012",i_rep_section=isec,r_val=bend_r012(start+isec),error=error)
       CALL section_vals_val_get(section,"R032",i_rep_section=isec,r_val=bend_r032(start+isec),error=error)
       CALL section_vals_val_get(section,"KBS12",i_rep_section=isec,r_val=bend_kbs12(start+isec),error=error)
       CALL section_vals_val_get(section,"KBS32",i_rep_section=isec,r_val=bend_kbs32(start+isec),error=error)
       CALL section_vals_val_get(section,"KSS",i_rep_section=isec,r_val=bend_kss(start+isec),error=error)
    END DO
  END SUBROUTINE read_bends_section

! *****************************************************************************
  SUBROUTINE read_ubs_section(ub_kind,ub_a, ub_b, ub_c, ub_k, ub_r0, section, start, error)
    INTEGER, DIMENSION(:), POINTER           :: ub_kind
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: ub_a, ub_b, ub_c
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ub_k
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ub_r0
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_ubs_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, k, n_items
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: Kvals
    TYPE(section_vals_type), POINTER         :: subsection

    failure = .FALSE.
    NULLIFY(atm_names)
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       subsection => section_vals_get_subs_vals(section,"UB",i_rep_section=isec,error=error)
       CALL section_vals_get(subsection,explicit=explicit,error=error)
       IF (explicit) THEN
          CALL section_vals_val_get(subsection,"KIND",i_val=ub_kind(start+isec),error=error)
          CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
          ub_a(start+isec) = atm_names(1)
          ub_b(start+isec) = atm_names(2)
          ub_c(start+isec) = atm_names(3)
          CALL uppercase(ub_a(start+isec))
          CALL uppercase(ub_b(start+isec))
          CALL uppercase(ub_c(start+isec))
          CALL section_vals_val_get(subsection,"K",r_vals=Kvals,error=error)
          CPPostcondition(SIZE(Kvals) <= 3, cp_failure_level, routineP, error, failure)
          ub_k(:,start+isec) = 0.0_dp
          DO k=1,SIZE(Kvals)
             ub_k(k,start+isec) = Kvals(k)
          END DO
          CALL section_vals_val_get(subsection,"R0",r_val=ub_r0(start+isec),error=error)
       END IF
    END DO
  END SUBROUTINE read_ubs_section

! *****************************************************************************
!> \brief Reads the TORSIONS section
!> \author teo
! *****************************************************************************
  SUBROUTINE read_torsions_section(torsion_kind,torsion_a, torsion_b, torsion_c, torsion_d, torsion_k,&
       torsion_cosphi0, torsion_m, section, start, error )
    INTEGER, DIMENSION(:), POINTER           :: torsion_kind
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: torsion_a, torsion_b, &
                                                torsion_c, torsion_d
    REAL(KIND=dp), DIMENSION(:), POINTER     :: torsion_k, torsion_cosphi0
    INTEGER, DIMENSION(:), POINTER           :: torsion_m
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_torsions_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY( atm_names)
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"KIND",i_rep_section=isec,i_val=torsion_kind(start+isec),error=error)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       torsion_a(start+isec) = atm_names(1)
       torsion_b(start+isec) = atm_names(2)
       torsion_c(start+isec) = atm_names(3)
       torsion_d(start+isec) = atm_names(4)
       CALL uppercase(torsion_a(start+isec))
       CALL uppercase(torsion_b(start+isec))
       CALL uppercase(torsion_c(start+isec))
       CALL uppercase(torsion_d(start+isec))
       CALL section_vals_val_get(section,"K",i_rep_section=isec,r_val=torsion_k(start+isec),error=error)
       CALL section_vals_val_get(section,"PHI0",i_rep_section=isec,r_val=torsion_cosphi0(start+isec),error=error)
       CALL section_vals_val_get(section,"M",i_rep_section=isec,i_val=torsion_m(start+isec),error=error)
    END DO
  END SUBROUTINE read_torsions_section

! *****************************************************************************
!> \brief Reads the force_field input section
!> \par History
!>      JGH (30.11.2001) : moved determination of setup variables to
!>                         molecule_input
!> \author CJM
! *****************************************************************************
  SUBROUTINE read_force_field_section ( ff_type , para_env, mm_section, error )
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(section_vals_type), POINTER         :: ff_section

    NULLIFY(ff_section)
    ff_section => section_vals_get_subs_vals(mm_section,"FORCEFIELD",error=error)
    CALL read_force_field_section1(ff_section, mm_section, ff_type, para_env, error)
  END SUBROUTINE read_force_field_section

! *****************************************************************************
!> \brief Reads the GROMOS force_field
!> \author ikuo
! *****************************************************************************
  SUBROUTINE read_force_field_gromos ( ff_type , para_env,  mm_section, error )

    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_force_field_gromos', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: ekt = 2.5_dp

    CHARACTER(LEN=default_string_length)     :: label
    CHARACTER(LEN=default_string_length), &
      DIMENSION(21)                          :: avail_section
    CHARACTER(LEN=default_string_length), &
      POINTER                                :: namearray(:)
    INTEGER                                  :: handle, iatom, icon, itemp, &
                                                itype, iw, jatom, ncon, &
                                                ntype, offset, stat
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: cost2, csq, sdet
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser
    TYPE(gromos_info_type), POINTER          :: gro_info

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger,parser)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    avail_section( 1) = "TITLE"
    avail_section( 2) = "TOPPHYSCON"
    avail_section( 3) = "TOPVERSION"
    avail_section( 4) = "ATOMTYPENAME"
    avail_section( 5) = "RESNAME"
    avail_section( 6) = "SOLUTEATOM"
    avail_section( 7) = "BONDTYPE"
    avail_section( 8) = "BONDH"
    avail_section( 9) = "BOND"
    avail_section(10) = "BONDANGLETYPE"
    avail_section(11) = "BONDANGLEH"
    avail_section(12) = "BONDANGLE"
    avail_section(13) = "IMPDIHEDRALTYPE"
    avail_section(14) = "IMPDIHEDRALH"
    avail_section(15) = "IMPDIHEDRAL"
    avail_section(16) = "DIHEDRALTYPE"
    avail_section(17) = "DIHEDRALH"
    avail_section(18) = "DIHEDRAL"
    avail_section(19) = "LJPARAMETERS"
    avail_section(20) = "SOLVENTATOM"
    avail_section(21) = "SOLVENTCONSTR"

    gro_info => ff_type%gro_info
    gro_info%ff_type = ff_type%ff_type
    NULLIFY(namearray)
    ! ATOMTYPENAME SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GTOP_INFO| Parsing the ATOMTYPENAME section'
    CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(4))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(namearray,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,namearray(itype),lower_to_upper=.TRUE.,error=error)
          IF(iw>0) WRITE(iw,*) "GTOP_INFO|  ",TRIM(namearray(itype))
       END DO
    END IF
    CALL parser_release(parser,error=error)

    ! SOLVENTCONSTR SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the SOLVENTATOM section'
    CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)

    label = TRIM(avail_section(21))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ncon,error=error)
       CALL reallocate(gro_info%solvent_k,1,ncon)
       CALL reallocate(gro_info%solvent_r0,1,ncon)
       DO icon=1,ncon
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,itemp,error=error)
          CALL parser_get_object(parser,itemp,error=error)
          CALL parser_get_object(parser,gro_info%solvent_r0(icon),error=error)
          gro_info%solvent_k(icon)=0.0_dp
       END DO
    END IF
    CALL parser_release(parser,error=error)

    CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
    ! BONDTYPE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the BONDTYPE section'
    label = TRIM(avail_section(7))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(gro_info%bond_k,1,ntype)
       CALL reallocate(gro_info%bond_r0,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,gro_info%bond_k(itype),error=error)
          CALL parser_get_object(parser,gro_info%bond_r0(itype),error=error)
          IF(ff_type%ff_type==do_ff_g96) THEN
             gro_info%bond_k(itype) = cp_unit_to_cp2k(gro_info%bond_k(itype),"kjmol*nm^-4",error=error)
          ELSE ! Assume its G87
             gro_info%bond_k(itype) = (2.0_dp) * gro_info%bond_k(itype) * gro_info%bond_r0(itype)**2
             gro_info%bond_k(itype) = cp_unit_to_cp2k(gro_info%bond_k(itype),"kjmol*nm^-2",error=error)
          END IF
          gro_info%bond_r0(itype)= cp_unit_to_cp2k(gro_info%bond_r0(itype),"nm",error=error)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT BONDTYPE INFO HERE!!!!"
       END DO
    END IF

    ! BONDANGLETYPE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the BONDANGLETYPE section'
    label = TRIM(avail_section(10))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(gro_info%bend_k,1,ntype)
       CALL reallocate(gro_info%bend_theta0,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,gro_info%bend_k(itype),error=error)
          CALL parser_get_object(parser,gro_info%bend_theta0(itype),error=error)
          gro_info%bend_theta0(itype) = cp_unit_to_cp2k(gro_info%bend_theta0(itype),"deg",error=error)
          IF(ff_type%ff_type==do_ff_g96) THEN
             gro_info%bend_theta0(itype) = COS(gro_info%bend_theta0(itype))
          ELSE ! Assume its G87
             cost2 = COS(gro_info%bend_theta0(itype))*COS(gro_info%bend_theta0(itype))
             sdet = cost2*cost2 - (2.0_dp*cost2-1.0_dp)*(1.0_dp-ekt/gro_info%bend_k(itype))
             csq = (cost2-SQRT(sdet))/(2.0_dp*cost2-1.0_dp)
             gro_info%bend_k(itype) = ekt/ACOS(csq)**2
          END IF
          gro_info%bend_k(itype) = cp_unit_to_cp2k(gro_info%bend_k(itype),"kjmol",error=error)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT BONDANGLETYPE INFO HERE!!!!"
       END DO
    END IF

    ! IMPDIHEDRALTYPE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the IMPDIHEDRALTYPE section'
    label = TRIM(avail_section(13))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(gro_info%impr_k,1,ntype)
       CALL reallocate(gro_info%impr_phi0,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,gro_info%impr_k(itype),error=error)
          CALL parser_get_object(parser,gro_info%impr_phi0(itype),error=error)
          gro_info%impr_phi0(itype) = cp_unit_to_cp2k(gro_info%impr_phi0(itype),"deg",error=error)
          gro_info%impr_k(itype) = cp_unit_to_cp2k(gro_info%impr_k(itype),"kjmol*deg^-2",error=error)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT IMPDIHEDRALTYPE INFO HERE!!!!"
       END DO
    END IF

    ! DIHEDRALTYPE SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the DIHEDRALTYPE section'
    label = TRIM(avail_section(16))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       CALL reallocate(gro_info%torsion_k,1,ntype)
       CALL reallocate(gro_info%torsion_m,1,ntype)
       CALL reallocate(gro_info%torsion_cosphi0,1,ntype)
       DO itype=1,ntype
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,gro_info%torsion_k(itype),error=error)
          CALL parser_get_object(parser,gro_info%torsion_cosphi0(itype),error=error)
          CALL parser_get_object(parser,gro_info%torsion_m(itype),error=error)
          gro_info%torsion_k(itype) = cp_unit_to_cp2k(gro_info%torsion_k(itype),"kjmol",error=error) 
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT DIHEDRALTYPE INFO HERE!!!!"
       END DO
    END IF

    CALL parser_release(parser,error=error)

    ! LJPARAMETERS SECTION
    IF(iw>0) WRITE(iw,'(T2,A)') 'GROMOS_FF| Parsing the LJPARAMETERS section'    
    CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
    label = TRIM(avail_section(19))
    CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
    IF(found) THEN
       CALL parser_get_next_line(parser,1,error=error)
       CALL parser_get_object(parser,ntype,error=error)
       offset = 0
       IF(ASSOCIATED(gro_info%nonbond_a)) offset = SIZE(gro_info%nonbond_a)
       ntype=SIZE(namearray)
       CALL reallocate(gro_info%nonbond_a,1,ntype)
       CALL reallocate(gro_info%nonbond_a_14,1,ntype)
       CALL reallocate(gro_info%nonbond_c6,1,ntype,1,ntype)
       CALL reallocate(gro_info%nonbond_c12,1,ntype,1,ntype)
       CALL reallocate(gro_info%nonbond_c6_14,1,ntype,1,ntype)
       CALL reallocate(gro_info%nonbond_c12_14,1,ntype,1,ntype)

       gro_info%nonbond_c12 = 0._dp
       gro_info%nonbond_c6  = 0._dp
       gro_info%nonbond_c12_14 = 0._dp
       gro_info%nonbond_c6_14  = 0._dp

       DO itype=1,ntype*(ntype+1)/2
          CALL parser_get_next_line(parser,1,error=error)
          CALL parser_get_object(parser,iatom,error=error)
          CALL parser_get_object(parser,jatom,error=error)
          IF(iatom==jatom) THEN
             gro_info%nonbond_a(iatom) = namearray(iatom)
             gro_info%nonbond_a_14(iatom) = namearray(iatom)
          END IF
          CALL parser_get_object(parser,gro_info%nonbond_c12(iatom,jatom),error=error)
          CALL parser_get_object(parser,gro_info%nonbond_c6(iatom,jatom),error=error)
          CALL parser_get_object(parser,gro_info%nonbond_c12_14(iatom,jatom),error=error)
          CALL parser_get_object(parser,gro_info%nonbond_c6_14(iatom,jatom),error=error)
          gro_info%nonbond_c6(iatom,jatom) = cp_unit_to_cp2k(gro_info%nonbond_c6(iatom,jatom),"kjmol*nm^6",error=error)
          gro_info%nonbond_c12(iatom,jatom) = cp_unit_to_cp2k(gro_info%nonbond_c12(iatom,jatom),"kjmol*nm^12",error=error)
          gro_info%nonbond_c6_14(iatom,jatom) = cp_unit_to_cp2k(gro_info%nonbond_c6_14(iatom,jatom),"kjmol*nm^6",error=error)
          gro_info%nonbond_c12_14(iatom,jatom) = cp_unit_to_cp2k(gro_info%nonbond_c12_14(iatom,jatom),"kjmol*nm^12",error=error)
          
          gro_info%nonbond_c6_14(jatom,iatom)= gro_info%nonbond_c6_14(iatom,jatom)
          gro_info%nonbond_c12_14(jatom,iatom)= gro_info%nonbond_c12_14(iatom,jatom)
          gro_info%nonbond_c6(jatom,iatom)= gro_info%nonbond_c6(iatom,jatom)
          gro_info%nonbond_c12(jatom,iatom)= gro_info%nonbond_c12(iatom,jatom)
          IF(iw>0) WRITE(iw,*) "GROMOS_FF| PUT LJPARAMETERS INFO HERE!!!!"
       END DO
    END IF
    CALL parser_release(parser,error=error)

    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(handle)

    DEALLOCATE(namearray,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE read_force_field_gromos

! *****************************************************************************
!> \brief Reads the charmm force_field
!> \author ikuo
! *****************************************************************************
  SUBROUTINE read_force_field_charmm ( ff_type, para_env, mm_section, error )

    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_force_field_charmm', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: label, string, string2, &
                                                string3, string4
    CHARACTER(LEN=default_string_length), &
      DIMENSION(1)                           :: bond_section
    CHARACTER(LEN=default_string_length), &
      DIMENSION(19)                          :: avail_section
    CHARACTER(LEN=default_string_length), &
      DIMENSION(2)                           :: angl_section, impr_section, &
                                                nbon_section, thet_section
    INTEGER                                  :: dummy, handle, ilab, iw, &
                                                nbend, nbond, nimpr, &
                                                nnonbond, nonfo, ntorsion, nub
    LOGICAL                                  :: failure, found
    TYPE(charmm_info_type), POINTER          :: chm_info
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger,parser)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)

    avail_section(1) = "BOND"     ; bond_section(1) = avail_section(1)
    avail_section(11)= "BONDS"
    avail_section(2) = "ANGL"     ; angl_section(1) = avail_section(2)
    avail_section(3) = "THETA"    ; angl_section(2) = avail_section(3)
    avail_section(12)= "THETAS"
    avail_section(13)= "ANGLE"
    avail_section(14)= "ANGLES"
    avail_section(4) = "DIHEDRAL" ; thet_section(1) = avail_section(4)
    avail_section(15)= "DIHEDRALS"
    avail_section(5) = "PHI"      ; thet_section(2) = avail_section(5)
    avail_section(6) = "IMPROPER" ; impr_section(1) = avail_section(6)
    avail_section(7) = "IMPH"     ; impr_section(2) = avail_section(7)
    avail_section(16)= "IMPHI"
    avail_section(8) = "NONBONDED"; nbon_section(1) = avail_section(8)
    avail_section(9) = "NBOND"    ; nbon_section(2) = avail_section(9)
    avail_section(10)= "HBOND"
    avail_section(17)= "NBFIX"
    avail_section(18)= "CMAP"
    avail_section(19)= "END"

    chm_info => ff_type%chm_info
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 1. Read in all the Bonds info from the param file here
    !      Vbond = Kb(b-b0)^2
    !      UNITS for Kb: [(kcal/mol)/(A^2)] to [Eh/(AU^2)]
    !-----------------------------------------------------------------------------
    nbond = 0
    DO ilab = 1, SIZE(bond_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(bond_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                CALL parser_get_object(parser,string2,error=error)
                CALL uppercase ( string2 )
                nbond = nbond + 1
                CALL reallocate(chm_info%bond_a,1,nbond)
                CALL reallocate(chm_info%bond_b,1,nbond)
                CALL reallocate(chm_info%bond_k,1,nbond)
                CALL reallocate(chm_info%bond_r0,1,nbond)
                chm_info%bond_a(nbond) = string
                chm_info%bond_b(nbond) = string2
                CALL parser_get_object(parser,chm_info%bond_k(nbond),error=error)
                CALL parser_get_object(parser,chm_info%bond_r0(nbond),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM BOND ",nbond,&
                     TRIM(chm_info%bond_a(nbond))," ",&
                     TRIM(chm_info%bond_b(nbond))," ",&
                     chm_info%bond_k(nbond),&
                     chm_info%bond_r0(nbond)
                ! Do some units conversion into internal atomic units
                chm_info%bond_r0(nbond) = cp_unit_to_cp2k(chm_info%bond_r0(nbond),"angstrom",error=error)
                chm_info%bond_k(nbond) = cp_unit_to_cp2k(chm_info%bond_k(nbond),"kcalmol*angstrom^-2",error=error)
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2. Read in all the Bends and UB info from the param file here
    !      Vangle = Ktheta(theta-theta0)^2
    !      UNITS for Ktheta: [(kcal/mol)/(rad^2)] to [Eh/(rad^2)]
    !      FACTOR of "2" rolled into Ktheta
    !      Vub = Kub(S-S0)^2
    !      UNITS for Kub: [(kcal/mol)/(A^2)] to [Eh/(AU^2)]
    !-----------------------------------------------------------------------------
    nbend = 0
    nub = 0
    DO ilab = 1, SIZE(angl_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(angl_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                CALL parser_get_object(parser,string2,error=error)
                CALL parser_get_object(parser,string3,error=error)
                CALL uppercase ( string2 )
                CALL uppercase ( string3 )
                nbend = nbend + 1
                CALL reallocate(chm_info%bend_a,1,nbend)
                CALL reallocate(chm_info%bend_b,1,nbend)
                CALL reallocate(chm_info%bend_c,1,nbend)
                CALL reallocate(chm_info%bend_k,1,nbend)
                CALL reallocate(chm_info%bend_theta0,1,nbend)
                chm_info%bend_a(nbend) = string
                chm_info%bend_b(nbend) = string2
                chm_info%bend_c(nbend) = string3
                CALL parser_get_object(parser,chm_info%bend_k(nbend),error=error)
                CALL parser_get_object(parser,chm_info%bend_theta0(nbend),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM BEND ",nbend,&
                     TRIM(chm_info%bend_a(nbend))," ",&
                     TRIM(chm_info%bend_b(nbend))," ",&
                     TRIM(chm_info%bend_c(nbend))," ",&
                     chm_info%bend_k(nbend),&
                     chm_info%bend_theta0(nbend)
                ! Do some units conversion into internal atomic units
                chm_info%bend_theta0(nbend) = cp_unit_to_cp2k(chm_info%bend_theta0(nbend),"deg",error=error)
                chm_info%bend_k(nbend) = cp_unit_to_cp2k(chm_info%bend_k(nbend),"kcalmol*rad^-2",error=error)
                IF (parser_test_next_token(parser,error=error) == "FLT") THEN
                   nub = nub + 1
                   CALL reallocate(chm_info%ub_a,1,nub)
                   CALL reallocate(chm_info%ub_b,1,nub)
                   CALL reallocate(chm_info%ub_c,1,nub)
                   CALL reallocate(chm_info%ub_k,1,nub)
                   CALL reallocate(chm_info%ub_r0,1,nub)
                   chm_info%ub_a(nub) = string
                   chm_info%ub_b(nub) = string2
                   chm_info%ub_c(nub) = string3
                   CALL parser_get_object(parser,chm_info%ub_k(nub),error=error)
                   CALL parser_get_object(parser,chm_info%ub_r0(nub),error=error)
                   IF(iw>0) WRITE(iw,*) "    CHM UB ",nub,&
                        TRIM(chm_info%ub_a(nub))," ",&
                        TRIM(chm_info%ub_b(nub))," ",&
                        TRIM(chm_info%ub_c(nub))," ",&
                        chm_info%ub_k(nub),&
                        chm_info%ub_r0(nub)
                   ! Do some units conversion into internal atomic units
                   chm_info%ub_r0(nub) = cp_unit_to_cp2k(chm_info%ub_r0(nub),"angstrom",error=error)
                   chm_info%ub_k(nub) = cp_unit_to_cp2k(chm_info%ub_k(nub),"kcalmol*angstrom^-2",error=error)
                END IF
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3. Read in all the Dihedrals info from the param file here
    !      Vtorsion = Kphi(1+COS(n(phi)-delta))
    !      UNITS for Kphi: [(kcal/mol)] to [Eh]
    !-----------------------------------------------------------------------------
    ntorsion = 0
    DO ilab = 1, SIZE(thet_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(thet_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                CALL parser_get_object(parser,string2,error=error)
                CALL parser_get_object(parser,string3,error=error)
                CALL parser_get_object(parser,string4,error=error)
                CALL uppercase ( string2 )
                CALL uppercase ( string3 )
                CALL uppercase ( string4 )
                ntorsion = ntorsion + 1
                CALL reallocate(chm_info%torsion_a,1,ntorsion)
                CALL reallocate(chm_info%torsion_b,1,ntorsion)
                CALL reallocate(chm_info%torsion_c,1,ntorsion)
                CALL reallocate(chm_info%torsion_d,1,ntorsion)
                CALL reallocate(chm_info%torsion_k,1,ntorsion)
                CALL reallocate(chm_info%torsion_m,1,ntorsion)
                CALL reallocate(chm_info%torsion_cosphi0,1,ntorsion)
                chm_info%torsion_a(ntorsion) = string
                chm_info%torsion_b(ntorsion) = string2
                chm_info%torsion_c(ntorsion) = string3
                chm_info%torsion_d(ntorsion) = string4
                CALL parser_get_object(parser,chm_info%torsion_k(ntorsion),error=error)
                CALL parser_get_object(parser,chm_info%torsion_m(ntorsion),error=error)
                CALL parser_get_object(parser,chm_info%torsion_cosphi0(ntorsion),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM TORSION ",ntorsion,&
                     TRIM(chm_info%torsion_a(ntorsion))," ",&
                     TRIM(chm_info%torsion_b(ntorsion))," ",&
                     TRIM(chm_info%torsion_c(ntorsion))," ",&
                     TRIM(chm_info%torsion_d(ntorsion))," ",&
                     chm_info%torsion_k(ntorsion),&
                     chm_info%torsion_m(ntorsion),&
                     chm_info%torsion_cosphi0(ntorsion)
                ! Do some units conversion into internal atomic units
                chm_info%torsion_cosphi0(ntorsion) = cp_unit_to_cp2k(chm_info%torsion_cosphi0(ntorsion),&
                     "deg",error=error)
                chm_info%torsion_k(ntorsion) = cp_unit_to_cp2k(chm_info%torsion_k(ntorsion),"kcalmol",&
                     error=error)
                chm_info%torsion_cosphi0(ntorsion)=&
                     COS(chm_info%torsion_cosphi0(ntorsion))
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4. Read in all the Improper info from the param file here
    !      Vimpr = Kpsi(psi-psi0)^2
    !      UNITS for Kpsi: [(kcal/mol)/(rad^2)] to [Eh/(rad^2)]
    !-----------------------------------------------------------------------------
    nimpr = 0
    DO ilab = 1, SIZE(impr_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(impr_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                CALL parser_get_object(parser,string2,error=error)
                CALL parser_get_object(parser,string3,error=error)
                CALL parser_get_object(parser,string4,error=error)
                CALL uppercase ( string2 )
                CALL uppercase ( string3 )
                CALL uppercase ( string4 )
                nimpr = nimpr + 1
                CALL reallocate(chm_info%impr_a,1,nimpr)
                CALL reallocate(chm_info%impr_b,1,nimpr)
                CALL reallocate(chm_info%impr_c,1,nimpr)
                CALL reallocate(chm_info%impr_d,1,nimpr)
                CALL reallocate(chm_info%impr_k,1,nimpr)
                CALL reallocate(chm_info%impr_phi0,1,nimpr)
                chm_info%impr_a(nimpr) = string
                chm_info%impr_b(nimpr) = string2
                chm_info%impr_c(nimpr) = string3
                chm_info%impr_d(nimpr) = string4
                CALL parser_get_object(parser,chm_info%impr_k(nimpr),error=error)
                CALL parser_get_object(parser,dummy,error=error)
                CALL parser_get_object(parser,chm_info%impr_phi0(nimpr),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM IMPROPERS ",nimpr,&
                     TRIM(chm_info%impr_a(nimpr))," ",&
                     TRIM(chm_info%impr_b(nimpr))," ",&
                     TRIM(chm_info%impr_c(nimpr))," ",&
                     TRIM(chm_info%impr_d(nimpr))," ",&
                     chm_info%impr_k(nimpr),&
                     chm_info%impr_phi0(nimpr)
                ! Do some units conversion into internal atomic units
                chm_info%impr_phi0(nimpr) = cp_unit_to_cp2k(chm_info%impr_phi0(nimpr),"deg",error=error)
                chm_info%impr_k(nimpr) = cp_unit_to_cp2k(chm_info%impr_k(nimpr),"kcalmol",error=error)
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 5. Read in all the Nonbonded info from the param file here
    !-----------------------------------------------------------------------------
    nnonbond = 0
    nonfo = 0
    DO ilab = 1, SIZE(nbon_section)
       CALL parser_create(parser,ff_type%ff_file_name,para_env=para_env,error=error)
       label = TRIM(nbon_section(ilab))
       DO
          CALL parser_search_string(parser,label,.TRUE.,found,begin_line=.TRUE.,error=error)
          IF (found) THEN
             CALL parser_get_object(parser,string,error=error)
             IF (INDEX(string,TRIM(label)) /= 1) CYCLE
             CALL charmm_get_next_line(parser,1,error=error)
             DO
                CALL parser_get_object(parser,string,error=error)
                CALL uppercase ( string )
                IF(ANY(string == avail_section)) EXIT
                nnonbond = nnonbond + 1
                CALL reallocate(chm_info%nonbond_a,1,nnonbond)
                CALL reallocate(chm_info%nonbond_eps,1,nnonbond)
                CALL reallocate(chm_info%nonbond_rmin2,1,nnonbond)
                chm_info%nonbond_a(nnonbond) = string
                CALL parser_get_object(parser,chm_info%nonbond_eps(nnonbond),error=error)
                CALL parser_get_object(parser,chm_info%nonbond_eps(nnonbond),error=error)
                CALL parser_get_object(parser,chm_info%nonbond_rmin2(nnonbond),error=error)
                IF(iw>0) WRITE(iw,*) "    CHM NONBOND ",nnonbond,&
                     TRIM(chm_info%nonbond_a(nnonbond))," ",&
                     chm_info%nonbond_eps(nnonbond),&
                     chm_info%nonbond_rmin2(nnonbond)
                chm_info%nonbond_rmin2(nnonbond) = cp_unit_to_cp2k(chm_info%nonbond_rmin2(nnonbond),&
                     "angstrom", error=error)
                chm_info%nonbond_eps(nnonbond) = cp_unit_to_cp2k(chm_info%nonbond_eps(nnonbond),&
                     "kcalmol",error=error)
                IF (parser_test_next_token(parser,error=error) == "FLT") THEN
                   nonfo = nonfo + 1
                   CALL reallocate(chm_info%nonbond_a_14,1,nonfo)
                   CALL reallocate(chm_info%nonbond_eps_14,1,nonfo)
                   CALL reallocate(chm_info%nonbond_rmin2_14,1,nonfo)
                   chm_info%nonbond_a_14(nonfo) = chm_info%nonbond_a(nnonbond)
                   CALL parser_get_object(parser,chm_info%nonbond_eps_14(nonfo),error=error)
                   CALL parser_get_object(parser,chm_info%nonbond_eps_14(nonfo),error=error)
                   CALL parser_get_object(parser,chm_info%nonbond_rmin2_14(nonfo),error=error)
                   IF(iw>0) WRITE(iw,*) "    CHM ONFO ",nonfo,&
                        TRIM(chm_info%nonbond_a_14(nonfo))," ",&
                        chm_info%nonbond_eps_14(nonfo),&
                        chm_info%nonbond_rmin2_14(nonfo)
                   chm_info%nonbond_rmin2_14(nonfo) = cp_unit_to_cp2k(chm_info%nonbond_rmin2_14(nonfo),&
                        "angstrom",error=error)
                   chm_info%nonbond_eps_14(nonfo) = cp_unit_to_cp2k(chm_info%nonbond_eps_14(nonfo),&
                        "kcalmol",error=error)
                END IF
                CALL charmm_get_next_line(parser,1,error=error)
             END DO
          ELSE
             EXIT
          END IF
       END DO
       CALL parser_release(parser,error=error)
    END DO
    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%FF_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE read_force_field_charmm

! *****************************************************************************
!> \brief This function is simply a wrap to the parser_get_next_line..
!>      Comments: This routine would not be necessary if the continuation
!>                char for CHARMM would not be the "-".. How can you choose this
!>                character in a file of numbers as a continuation char????
!>                This sounds simply crazy....
!> \author Teodoro Laino - Zurich University - 06.2007
! *****************************************************************************
  SUBROUTINE charmm_get_next_line(parser, nline, error)
    TYPE(cp_parser_type), POINTER            :: parser
    INTEGER, INTENT(IN)                      :: nline
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'charmm_get_next_line', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=1), PARAMETER              :: continuation_char = "-"

    INTEGER                                  :: i, len_line

    DO i = 1, nline
       len_line = LEN_TRIM(parser%input_line)
       DO WHILE (parser%input_line(len_line:len_line)==continuation_char)
          CALL parser_get_next_line(parser,1,error=error)
          len_line = LEN_TRIM(parser%input_line)
       END DO
       CALL parser_get_next_line(parser,1,error=error)
    END DO

  END SUBROUTINE charmm_get_next_line

! *****************************************************************************
!> \brief reads EAM potential from library
! *****************************************************************************
  SUBROUTINE read_eam_data ( eam, para_env, mm_section, error)
    TYPE(eam_pot_type), POINTER              :: eam
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_eam_data', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iw
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_parser_type), POINTER            :: parser

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(parser, logger)
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%FF_INFO",&
         extension=".mmLog",error=error)
    IF (.NOT.failure) THEN
       IF (iw>0) WRITE ( iw, * ) "Reading EAM data from: ",TRIM(eam%eam_file_name)
       CALL parser_create(parser,TRIM(eam%eam_file_name),para_env=para_env,error=error)

       CALL parser_get_next_line(parser,1,error=error)
       IF (iw>0)  WRITE(iw,*) "Title: ",parser%input_line

       CALL parser_get_next_line(parser,2,error=error)
       READ(parser%input_line,*) eam%drar, eam%drhoar, eam%acutal, eam%npoints
       eam%drar = cp_unit_to_cp2k(eam%drar,"angstrom",error=error)
       eam%acutal = cp_unit_to_cp2k(eam%acutal,"angstrom",error=error)
       ! Relocating arrays with the right size
       CALL reallocate(eam%rho, 1, eam%npoints)
       CALL reallocate(eam%rhop, 1, eam%npoints)
       CALL reallocate(eam%rval, 1, eam%npoints)
       CALL reallocate(eam%rhoval, 1, eam%npoints)
       CALL reallocate(eam%phi, 1, eam%npoints)
       CALL reallocate(eam%phip, 1, eam%npoints)
       CALL reallocate(eam%frho, 1, eam%npoints)
       CALL reallocate(eam%frhop, 1, eam%npoints)
       ! Reading density and derivative of density (with respect to r)
       DO i = 1,eam%npoints
          CALL parser_get_next_line(parser,1,error=error)
          READ (parser%input_line,*) eam%rho(i), eam%rhop(i)
          eam%rhop(i) = cp_unit_to_cp2k(eam%rhop(i),"angstrom^-1",error=error)
          eam%rval(i)   = REAL(i-1,KIND=dp)*eam%drar
          eam%rhoval(i) = REAL(i-1,KIND=dp)*eam%drhoar
       END DO
       ! Reading pair potential PHI and its derivative (with respect to r)
       DO i = 1,eam%npoints
          CALL parser_get_next_line(parser,1,error=error)
          READ (parser%input_line,*) eam%phi(i), eam%phip(i)
          eam%phi(i)  = cp_unit_to_cp2k(eam%phi(i),"eV",error=error)
          eam%phip(i) = cp_unit_to_cp2k(eam%phip(i),"eV*angstrom^-1",error=error)
       END DO
       ! Reading embedded function and its derivative (with respect to density)
       DO i = 1,eam%npoints
          CALL parser_get_next_line(parser,1,error=error)
          READ (parser%input_line,*) eam%frho(i), eam%frhop(i)
          eam%frho(i) = cp_unit_to_cp2k(eam%frho(i),"eV",error=error)
          eam%frhop(i) = cp_unit_to_cp2k(eam%frhop(i),"eV",error=error)
       END DO

       IF (iw>0) WRITE(iw,*)"Finished EAM data"
       CALL parser_release(parser,error=error)
     END IF
    CALL cp_print_key_finished_output(iw,logger,mm_section,"PRINT%FF_INFO",error=error)
    CALL timestop(handle)

  END SUBROUTINE read_eam_data

END MODULE force_fields_input
