!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE bibliography,                    ONLY: Siepmann1995,&
                                             Tersoff1988,&
                                             Tosi1964a,&
                                             Tosi1964b,&
                                             Yamada2000,&
                                             cite_reference
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_parser_methods,               ONLY: parser_get_next_line
  USE cp_parser_types,                 ONLY: cp_parser_type,&
                                             parser_create,&
                                             parser_release
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE damping_dipole_types,            ONLY: damping_info_type
  USE f77_blas
  USE force_field_kind_types,          ONLY: do_ff_amber,&
                                             do_ff_charmm,&
                                             do_ff_g87,&
                                             do_ff_g96,&
                                             do_ff_undef
  USE force_field_types,               ONLY: force_field_type,&
                                             input_info_type
  USE force_fields_util,               ONLY: get_generic_info
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathlib,                         ONLY: invert_matrix
  USE memory_utilities,                ONLY: reallocate
  USE pair_potential_types,            ONLY: &
       b4_type, bm_type, do_potential_single_allocation, ea_type, &
       eam_pot_type, ft_pot_type, ft_type, ftd_pot_type, ftd_type, gp_type, &
       gw_type, ip_type, ipbv_pot_type, lj_charmm_type, &
       no_potential_single_allocation, pair_potential_p_type, &
       pair_potential_reallocate, potential_single_allocation, quip_type, &
       siepmann_type, tersoff_type, wl_type
  USE shell_potential_types,           ONLY: shell_p_create,&
                                             shell_p_type
  USE string_utilities,                ONLY: uppercase
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'force_fields_input'

  PRIVATE
  PUBLIC :: read_force_field_section,&
            read_lj_section,&
            read_wl_section,&
            read_gd_section,&
            read_gp_section,&
            read_chrg_section

CONTAINS

! *****************************************************************************
!> \brief Reads the force_field input section
!> \param ff_section ...
!> \param mm_section ...
!> \param ff_type ...
!> \param para_env ...
!> \param error ...
!> \author teo
! *****************************************************************************
  SUBROUTINE read_force_field_section1 ( ff_section, mm_section, ff_type, para_env, error )
    TYPE(section_vals_type), POINTER         :: ff_section, mm_section
    TYPE(force_field_type), INTENT(INOUT)    :: ff_type
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_force_field_section1', &
      routineP = moduleN//':'//routineN

    INTEGER :: nb4, nbends, nbm, nbmhft, nbmhftd, nbonds, nchg, neam, ngd, &
      ngp, nimpr, nipbv, nlj, nopbend, nquip, nshell, nsiepmann, ntersoff, &
      ntors, ntot, nubs, nwl
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
       CALL section_vals_val_get(ff_section,"SPLINE%NPOINTS",i_val=ff_type%npoints,error=error)
       CALL section_vals_val_get(ff_section,"IGNORE_MISSING_CRITICAL_PARAMS",l_val=ff_type%ignore_missing_critical,error=error)
       CPPostcondition(ff_type%max_energy<=ff_type%emax_spline,cp_failure_level,routineP,error,failure)
       ! Read the parameter file name only if the force field type requires it..
       SELECT CASE(ff_type%ff_type)
       CASE(do_ff_charmm,do_ff_amber,do_ff_g96,do_ff_g87)
          CALL section_vals_val_get(ff_section,"PARM_FILE_NAME",c_val=ff_type%ff_file_name,error=error)
          CALL cp_assert(TRIM(ff_type%ff_file_name)/="",cp_failure_level,cp_assertion_failed,routineP,&
               "Force Field Parameter's filename is empty! Please check your input file. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       CASE(do_ff_undef)
          ! Do Nothing
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Force field type not implemented",para_env)
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
       CALL section_vals_val_get(ff_section,"SHIFT_CUTOFF",l_val=ff_type%shift_cutoff,error=error)
       CALL section_vals_val_get(ff_section,"SPLINE%UNIQUE_SPLINE",l_val=unique_spline,error=error)
       ! Single spline
       potential_single_allocation = no_potential_single_allocation
       IF (unique_spline) potential_single_allocation = do_potential_single_allocation

       CALL section_vals_val_get(ff_section,"MULTIPLE_POTENTIAL",l_val=ff_type%multiple_potential,error=error)
       CALL section_vals_val_get(ff_section,"DO_NONBONDED",l_val=ff_type%do_nonbonded,error=error)
       tmp_section  => section_vals_get_subs_vals(ff_section,"NONBONDED",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,error=error)
       IF (explicit.AND.ff_type%do_nonbonded) THEN
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

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"BMHFTD",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nbmhftd,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nbmhftd,bmhftd=.TRUE.,error=error)
             CALL read_bmhftd_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"BUCK4RANGES",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nb4,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nbmhftd
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nb4,buck4r=.TRUE.,error=error)
             CALL read_b4_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"BUCKMORSE",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nbm,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nbmhftd+nb4
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nbm,buckmo=.TRUE.,error=error)
             CALL read_bm_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"GENPOT",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=ngp,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nbmhftd+nb4+nbm
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+ngp,gp=.TRUE.,error=error)
             CALL read_gp_section(inp_info%nonbonded,tmp_section2,ntot,error)
          END IF
          tmp_section2 => section_vals_get_subs_vals(tmp_section,"TERSOFF",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=ntersoff,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nbmhftd+nb4+nbm+ngp
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+ntersoff,tersoff=.TRUE.,error=error)
             CALL read_tersoff_section(inp_info%nonbonded,tmp_section2,ntot,tmp_section2,error)
          END IF
          tmp_section2 => section_vals_get_subs_vals(tmp_section,"SIEPMANN",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nsiepmann,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nbmhftd+nb4+nbm+ngp+ntersoff
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nsiepmann,siepmann=.TRUE.,error=error)
             CALL read_siepmann_section(inp_info%nonbonded,tmp_section2,ntot,tmp_section2,error)
          END IF

          tmp_section2 => section_vals_get_subs_vals(tmp_section,"quip",error=error)
          CALL section_vals_get(tmp_section2,explicit=explicit,n_repetition=nquip,error=error)
          ntot=nlj+nwl+neam+ngd+nipbv+nbmhft+nbmhftd+nb4+nbm+ngp+ntersoff+nsiepmann
          IF (explicit) THEN
             CALL pair_potential_reallocate(inp_info%nonbonded,1,ntot+nquip,quip=.TRUE.,error=error)
             CALL read_quip_section(inp_info%nonbonded,tmp_section2,ntot,para_env,mm_section,error)
          END IF
       END IF

       tmp_section => section_vals_get_subs_vals(ff_section,"NONBONDED14",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,error=error)
       IF (explicit.AND.ff_type%do_nonbonded) THEN
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
       tmp_section => section_vals_get_subs_vals(ff_section,"DIPOLE",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nchg,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%apol_atm,1,nchg)
          CALL reallocate(inp_info%apol,1,nchg)
          CALL read_apol_section(inp_info%apol_atm,inp_info%apol,inp_info%damping_list,&
                                 tmp_section,ntot,error)
       END IF
       tmp_section => section_vals_get_subs_vals(ff_section,"QUADRUPOLE",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nchg,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%cpol_atm,1,nchg)
          CALL reallocate(inp_info%cpol,1,nchg)
          CALL read_cpol_section(inp_info%cpol_atm,inp_info%cpol,tmp_section,ntot,error)
       END IF
       tmp_section => section_vals_get_subs_vals(ff_section,"SHELL",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nshell,error=error)
       IF (explicit) THEN
          ntot=0
          CALL shell_p_create(inp_info%shell_list,nshell,error)
          CALL read_shell_section(inp_info%shell_list,tmp_section,ntot,para_env,error)
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
          CALL reallocate(inp_info%torsion_phi0,1,ntors)
          CALL read_torsions_section(inp_info%torsion_kind,inp_info%torsion_a, inp_info%torsion_b,&
               inp_info%torsion_c,inp_info%torsion_d,inp_info%torsion_k,inp_info%torsion_phi0,&
               inp_info%torsion_m,tmp_section, ntot, error )
       END IF

       tmp_section => section_vals_get_subs_vals(ff_section,"IMPROPER",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nimpr,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%impr_kind,1,nimpr)
          CALL reallocate(inp_info%impr_a,1,nimpr)
          CALL reallocate(inp_info%impr_b,1,nimpr)
          CALL reallocate(inp_info%impr_c,1,nimpr)
          CALL reallocate(inp_info%impr_d,1,nimpr)
          CALL reallocate(inp_info%impr_k,1,nimpr)
          CALL reallocate(inp_info%impr_phi0,1,nimpr)
          CALL read_improper_section(inp_info%impr_kind,inp_info%impr_a, inp_info%impr_b,&
               inp_info%impr_c,inp_info%impr_d,inp_info%impr_k,inp_info%impr_phi0,&
               tmp_section, ntot, error )
       END IF

       tmp_section => section_vals_get_subs_vals(ff_section,"OPBEND",error=error)
       CALL section_vals_get(tmp_section,explicit=explicit,n_repetition=nopbend,error=error)
       IF (explicit) THEN
          ntot=0
          CALL reallocate(inp_info%opbend_kind,1,nopbend)
          CALL reallocate(inp_info%opbend_a,1,nopbend)
          CALL reallocate(inp_info%opbend_b,1,nopbend)
          CALL reallocate(inp_info%opbend_c,1,nopbend)
          CALL reallocate(inp_info%opbend_d,1,nopbend)
          CALL reallocate(inp_info%opbend_k,1,nopbend)
          CALL reallocate(inp_info%opbend_phi0,1,nopbend)
          CALL read_opbend_section(inp_info%opbend_kind,inp_info%opbend_a, inp_info%opbend_b,&
               inp_info%opbend_c,inp_info%opbend_d,inp_info%opbend_k,inp_info%opbend_phi0,&
               tmp_section, ntot, error )
       END IF
    END IF

  END SUBROUTINE read_force_field_section1

! *****************************************************************************
!> \brief Set up of the IPBV force fields
!> \param at1 ...
!> \param at2 ...
!> \param ipbv ...
!> \author teo
! *****************************************************************************
  SUBROUTINE set_IPBV_ff(at1, at2, ipbv )
    CHARACTER(LEN=*), INTENT(IN)             :: at1, at2
    TYPE(ipbv_pot_type), POINTER             :: ipbv

    CHARACTER(len=*), PARAMETER :: routineN = 'set_IPBV_ff', &
      routineP = moduleN//':'//routineN

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
       CALL stop_program (routineN,moduleN,__LINE__,"IPBV only for WATER")
    ENDIF
  END SUBROUTINE set_IPBV_ff

! *****************************************************************************
!> \brief Set up of the BMHFT force fields
!> \param at1 ...
!> \param at2 ...
!> \param ft ...
!> \param error ...
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
       CALL stop_program(routineN,moduleN,__LINE__,&
                         "BMHFT only for NaCl")
    ENDIF

  END SUBROUTINE set_BMHFT_ff

! *****************************************************************************
!> \brief Set up of the BMHFTD force fields
!> \param at1 ...
!> \param at2 ...
!> \param ftd ...
!> \param error ...
!> \author Mathieu Salanne 05.2010
! *****************************************************************************
  SUBROUTINE set_BMHFTD_ff(at1, at2, ftd, error)
    CHARACTER(LEN=*), INTENT(IN)             :: at1, at2
    TYPE(ftd_pot_type), POINTER              :: ftd
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_BMHFTD_ff', &
      routineP = moduleN//':'//routineN

    CALL stop_program(routineN,moduleN,__LINE__,&
                      "No default parameters present for BMHFTD")

  END SUBROUTINE set_BMHFTD_ff

! *****************************************************************************
!> \brief Reads the EAM section
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param para_env ...
!> \param mm_section ...
!> \param error ...
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
!> \brief Reads the QUIP section
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param para_env ...
!> \param mm_section ...
!> \param error ...
!> \author teo
! *****************************************************************************
  SUBROUTINE read_quip_section(nonbonded, section, start, para_env, mm_section, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_quip_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: args_str, atm_names
    INTEGER                                  :: is, isec, n_calc_args, n_items

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)

       nonbonded%pot(start+isec)%pot %type = quip_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)
       CALL section_vals_val_get(section,"PARM_FILE_NAME",i_rep_section=isec,&
            c_val=nonbonded%pot(start+isec)%pot%set(1)%quip%quip_file_name,error=error)
       CALL section_vals_val_get(section,"INIT_ARGS",i_rep_section=isec,&
            c_vals=args_str, error=error)
       nonbonded%pot(start+isec)%pot%set(1)%quip%init_args = ""
       DO is=1, SIZE(args_str)
          nonbonded%pot(start+isec)%pot%set(1)%quip%init_args = &
             TRIM(nonbonded%pot(start+isec)%pot%set(1)%quip%init_args) // &
             " "//TRIM(args_str(is))
       END DO ! is
       CALL section_vals_val_get(section,"CALC_ARGS",i_rep_section=isec,&
            n_rep_val=n_calc_args,error=error)
       IF (n_calc_args > 0) THEN
          CALL section_vals_val_get(section,"CALC_ARGS",i_rep_section=isec,&
               c_vals=args_str,error=error)
               DO is=1, SIZE(args_str)
                  nonbonded%pot(start+isec)%pot%set(1)%quip%calc_args =  &
                     TRIM(nonbonded%pot(start+isec)%pot%set(1)%quip%calc_args) // &
                     " "//TRIM(args_str(is))
               END DO ! is
       END IF
       nonbonded%pot(start+isec)%pot%rcutsq  = 0.0_dp
    END DO
  END SUBROUTINE read_quip_section

! *****************************************************************************
!> \brief Reads the LJ section
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \brief Reads the BMHFTD section
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
!> \author Mathieu Salanne 05.2010
! *****************************************************************************
  SUBROUTINE read_bmhftd_section(nonbonded, section, start, error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_bmhftd_section', &
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
       nonbonded%pot(start+isec)%pot %type = ftd_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)

       CALL section_vals_val_get(section,"A",i_rep_section=isec,n_rep_val=i,error=error)
       IF (i==1) THEN
          CALL section_vals_val_get(section,"A",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ftd%a,error=error)
          CALL section_vals_val_get(section,"B",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ftd%b,error=error)
          CALL section_vals_val_get(section,"C",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ftd%c,error=error)
          CALL section_vals_val_get(section,"D",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ftd%d,error=error)
          CALL section_vals_val_get(section,"BD",i_rep_section=isec,&
               r_val=nonbonded%pot(start+isec)%pot%set(1)%ftd%bd,error=error)

       ELSE
          CALL section_vals_val_get(section,"MAP_ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
          map_atoms = atm_names
          CALL uppercase(map_atoms(1))
          CALL uppercase(map_atoms(2))
          CALL set_BMHFTD_ff(map_atoms(1), map_atoms(2), nonbonded%pot(start+isec)%pot%set(1)%ftd,&
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
  END SUBROUTINE read_bmhftd_section

! *****************************************************************************
!> \brief Reads the Buckingham 4 Ranges potential section
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
!> \par History
!>      MK (11.11.2010): Automatic fit of the (default) polynomial coefficients
!> \author MI,MK
! *****************************************************************************
  SUBROUTINE read_b4_section(nonbonded,section,start,error)

    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_b4_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: i, ir, isec, istat, n_items, &
                                                n_rep, np1, np2
    LOGICAL                                  :: explicit_poly1, &
                                                explicit_poly2, failure
    REAL(KIND=dp)                            :: a, b, c, eval_error, r1, r2, &
                                                r3, rcut
    REAL(KIND=dp), DIMENSION(10)             :: v, x
    REAL(KIND=dp), DIMENSION(10, 10)         :: p, p_inv
    REAL(KIND=dp), DIMENSION(:), POINTER     :: coeff1, coeff2, list

    failure = .FALSE.

    NULLIFY (coeff1)
    NULLIFY (coeff2)

    CALL section_vals_get(section,n_repetition=n_items,error=error)

    DO isec=1,n_items
      CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
      CALL section_vals_val_get(section,"A",i_rep_section=isec,r_val=a,error=error)
      CALL section_vals_val_get(section,"B",i_rep_section=isec,r_val=b,error=error)
      CALL section_vals_val_get(section,"C",i_rep_section=isec,r_val=c,error=error)
      CALL section_vals_val_get(section,"R1",i_rep_section=isec,r_val=r1,error=error)
      CALL section_vals_val_get(section,"R2",i_rep_section=isec,r_val=r2,error=error)
      CALL section_vals_val_get(section,"R3",i_rep_section=isec,r_val=r3,error=error)
      CALL section_vals_val_get(section,"POLY1",explicit=explicit_poly1,n_rep_val=n_rep,error=error)
      ! Check if polynomial coefficients were specified for range 2 and 3 explicitly
      IF (explicit_poly1) THEN
        np1 = 0
        DO ir=1,n_rep
          NULLIFY (list)
          CALL section_vals_val_get(section,"POLY1",i_rep_val=ir,r_vals=list,error=error)
          IF (ASSOCIATED(list)) THEN
            CALL reallocate(coeff1,0,np1+SIZE(list)-1)
            DO i=1,SIZE(list)
              coeff1(i+np1-1) = list(i)
            END DO
            np1 = np1 + SIZE(list)
          END IF
        END DO
      END IF
      CALL section_vals_val_get(section,"POLY2",explicit=explicit_poly2,n_rep_val=n_rep,error=error)
      IF (explicit_poly2) THEN
        np2 = 0
        DO ir=1,n_rep
          NULLIFY (list)
          CALL section_vals_val_get(section,"POLY2",i_rep_val=ir,r_vals=list,error=error)
          IF (ASSOCIATED(list)) THEN
            CALL reallocate(coeff2,0,np2+SIZE(list)-1)
            DO i=1,SIZE(list)
              coeff2(i+np2-1) = list(i)
            END DO
            np2 = np2 + SIZE(list)
          END IF
        END DO
      END IF
      ! Default is a 5th/3rd-order polynomial fit
      IF ((.NOT.explicit_poly1).OR.(.NOT.explicit_poly2)) THEN
        ! Build matrix p and vector v to calculate the polynomial coefficients
        ! in the vector x from p*x = v
        p(:,:) = 0.0_dp
        ! Row 1: Match the 5th-order polynomial and the potential at r1
        p(1,1) = 1.0_dp
        DO i=2,6
          p(1,i) = p(1,i-1)*r1
        END DO
        ! Row 2: Match the first derivatives of the 5th-order polynomial and the potential at r1
        DO i=2,6
          p(2,i) = REAL(i-1,KIND=dp)*p(1,i-1)
        END DO
        ! Row 3: Match the second derivatives of the 5th-order polynomial and the potential at r1
        DO i=3,6
          p(3,i) = REAL(i-1,KIND=dp)*p(2,i-1)
        END DO
        ! Row 4: Match the 5th-order and the 3rd-order polynomials at r2
        p(4,1) = 1.0_dp
        DO i=2,6
          p(4,i) = p(4,i-1)*r2
        END DO
        p(4,7) = -1.0_dp
        DO i=8,10
          p(4,i) = p(4,i-1)*r2
        END DO
        ! Row 5: Match the first derivatives of the 5th-order and the 3rd-order polynomials at r2
        DO i=2,6
          p(5,i) = REAL(i-1,KIND=dp)*p(4,i-1)
        END DO
        DO i=8,10
          p(5,i) = REAL(i-7,KIND=dp)*p(4,i-1)
        END DO
        ! Row 6: Match the second derivatives of the 5th-order and the 3rd-order polynomials at r2
        DO i=3,6
          p(6,i) = REAL(i-1,KIND=dp)*p(5,i-1)
        END DO
        DO i=9,10
          p(6,i) = REAL(i-7,KIND=dp)*p(5,i-1)
        END DO
        ! Row 7: Minimum at r2, ie. the first derivative of the 3rd-order polynomial has to be zero at r2
        DO i=8,10
          p(7,i) = -p(5,i)
        END DO
        ! Row 8: Match the 3rd-order polynomial and the potential at r3
        p(8,7) = 1.0_dp
        DO i=8,10
          p(8,i) = p(8,i-1)*r3
        END DO
        ! Row 9: Match the first derivatives of the 3rd-order polynomial and the potential at r3
        DO i=8,10
          p(9,i) = REAL(i-7,KIND=dp)*p(8,i-1)
        END DO
        ! Row 10: Match the second derivatives of the 3rd-order polynomial and the potential at r3
        DO i=9,10
          p(10,i) = REAL(i-7,KIND=dp)*p(9,i-1)
        END DO
        ! Build the vector v
        v(1)   = a*EXP(-b*r1)
        v(2)   = -b*v(1)
        v(3)   = -b*v(2)
        v(4:7) = 0.0_dp
        v(8)   = -c/p(8,10)**2 ! = -c/r3**6
        v(9)   = -6.0_dp*v(8)/r3
        v(10)  = -7.0_dp*v(9)/r3
        ! Calculate p_inv the inverse of the matrix p
        p_inv(:,:) = 0.0_dp
        CALL invert_matrix(p,p_inv,eval_error,error=error)
        CALL cp_assert((eval_error < 1.0E-8_dp),cp_warning_level,cp_assertion_failed,routineP,&
                       "The polynomial fit for the BUCK4RANGES potential is only accurate to "//&
                       TRIM(cp_to_string(eval_error)),only_ionode=.TRUE.)
        ! Get the 6 coefficients of the 5th-order polynomial -> x(1:6)
        ! and the 4 coefficients of the 3rd-order polynomial -> x(7:10)
        x(:) = MATMUL(p_inv(:,:),v(:))
      END IF

      CALL section_vals_val_get(section,"RCUT",i_rep_section=isec,r_val=rcut,error=error)

      nonbonded%pot(start+isec)%pot%type = b4_type
      nonbonded%pot(start+isec)%pot%at1  = atm_names(1)
      nonbonded%pot(start+isec)%pot%at2  = atm_names(2)
      CALL uppercase(nonbonded%pot(start+isec)%pot%at1)
      CALL uppercase(nonbonded%pot(start+isec)%pot%at2)
      nonbonded%pot(start+isec)%pot%set(1)%buck4r%a = a
      nonbonded%pot(start+isec)%pot%set(1)%buck4r%b = b
      nonbonded%pot(start+isec)%pot%set(1)%buck4r%c = c
      nonbonded%pot(start+isec)%pot%set(1)%buck4r%r1 = r1
      nonbonded%pot(start+isec)%pot%set(1)%buck4r%r2 = r2
      nonbonded%pot(start+isec)%pot%set(1)%buck4r%r3 = r3
      IF ((.NOT.explicit_poly1).OR.(.NOT.explicit_poly2)) THEN
        nonbonded%pot(start+isec)%pot%set(1)%buck4r%npoly1 = 5
        nonbonded%pot(start+isec)%pot%set(1)%buck4r%poly1(0:5) = x(1:6)
        nonbonded%pot(start+isec)%pot%set(1)%buck4r%npoly2 = 3
        nonbonded%pot(start+isec)%pot%set(1)%buck4r%poly2(0:3) = x(7:10)
      ELSE
        nonbonded%pot(start+isec)%pot%set(1)%buck4r%npoly1 = np1-1
        CPPostcondition((np1 - 1 <= 10),cp_failure_level,routineP,error,failure)
        nonbonded%pot(start+isec)%pot%set(1)%buck4r%poly1(0:np1-1) = coeff1(0:np1-1)
        nonbonded%pot(start+isec)%pot%set(1)%buck4r%npoly2 = np2-1
        CPPostcondition((np2 - 1 <= 10),cp_failure_level,routineP,error,failure)
        nonbonded%pot(start+isec)%pot%set(1)%buck4r%poly2(0:np2-1) = coeff2(0:np2-1)
      END IF
      nonbonded%pot(start+isec)%pot%rcutsq = rcut*rcut

      IF (ASSOCIATED(coeff1)) THEN
        DEALLOCATE (coeff1,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF (ASSOCIATED(coeff2)) THEN
        DEALLOCATE (coeff2,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,n_rep_val=n_rep,error=error)
      IF (n_rep == 1) CALL section_vals_val_get(section,"RMIN",i_rep_section=isec,&
                                                r_val=nonbonded%pot(start+isec)%pot%set(1)%rmin,&
                                                error=error)
      CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,n_rep_val=n_rep,error=error)
      IF (n_rep == 1) CALL section_vals_val_get(section,"RMAX",i_rep_section=isec,&
                                                r_val=nonbonded%pot(start+isec)%pot%set(1)%rmax,&
                                                error=error)
    END DO

  END SUBROUTINE read_b4_section

! *****************************************************************************
!> \brief Reads the GENPOT - generic potential section
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param tersoff_section ...
!> \param error ...
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
!> \brief Reads the siepmann section
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param siepmann_section ...
!> \param error ...
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE read_siepmann_section(nonbonded, section, start, siepmann_section,error)
    TYPE(pair_potential_p_type), POINTER     :: nonbonded
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(section_vals_type), POINTER         :: siepmann_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_siepmann_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items, n_rep
    REAL(KIND=dp)                            :: rcut

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL cite_reference(Siepmann1995)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)

       nonbonded%pot(start+isec)%pot %type = siepmann_type
       nonbonded%pot(start+isec)%pot %at1  = atm_names(1)
       nonbonded%pot(start+isec)%pot %at2  = atm_names(2)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at1)
       CALL uppercase(nonbonded%pot(start+isec)%pot %at2)

       CALL section_vals_val_get(siepmann_section,"B",r_val=nonbonded%pot(start+isec)%pot%set(1)%siepmann%B,error=error)
       CALL section_vals_val_get(siepmann_section,"D",r_val=nonbonded%pot(start+isec)%pot%set(1)%siepmann%D,error=error)
       CALL section_vals_val_get(siepmann_section,"E",r_val=nonbonded%pot(start+isec)%pot%set(1)%siepmann%E,error=error)
       CALL section_vals_val_get(siepmann_section,"F",r_val=nonbonded%pot(start+isec)%pot%set(1)%siepmann%F,error=error)
       CALL section_vals_val_get(siepmann_section,"beta",r_val=nonbonded%pot(start+isec)%pot%set(1)%siepmann%beta,error=error)

     ! ! In case it is defined override the standard specification of RCUT
      CALL section_vals_val_get ( siepmann_section, "RCUT", n_rep_val=n_rep, error=error)
      IF(n_rep==1) THEN
         CALL section_vals_val_get(siepmann_section,"RCUT",r_val=rcut,error=error)
         nonbonded%pot(start+isec)%pot %rcutsq = rcut**2
         nonbonded%pot(start+isec)%pot%set(1)%siepmann%rcutsq =  rcut**2
      END IF
    END DO
  END SUBROUTINE read_siepmann_section

! *****************************************************************************
!> \brief Reads the Buckingham plus Morse potential section
!> \param nonbonded ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \param charge_atm ...
!> \param charge ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \brief Reads the POLARIZABILITY section
!> \param apol_atm ...
!> \param apol ...
!> \param damping_list ...
!> \param section ...
!> \param start ...
!> \param error ...
!> \author Marcel Baer
! *****************************************************************************
  SUBROUTINE read_apol_section(apol_atm, apol, damping_list, section,&
                               start, error)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: apol_atm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: apol
    TYPE(damping_info_type), DIMENSION(:), &
      POINTER                                :: damping_list
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_apol_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atm_name
    INTEGER                                  :: isec, isec_damp, n_damp, &
                                                n_items, start_damp, stat, &
                                                tmp_damp
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: tmp_section

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    NULLIFY(tmp_section)
    n_damp=0
! *** Counts number of DIPOLE%DAMPING sections ****
    DO isec = 1, n_items
      tmp_section => section_vals_get_subs_vals(section,"DAMPING",&
                                     i_rep_section=isec,error=error)
      CALL section_vals_get(tmp_section,n_repetition=tmp_damp,error=error)
      n_damp=n_damp+tmp_damp

#ifdef DEBUG_DAMPING
      WRITE(*,*) 'DAMPING DEBUG',isec,tmp_damp
#endif

    END DO

    IF (n_damp > 0) THEN
      failure=.FALSE.
      ALLOCATE (damping_list(1:n_damp),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ENDIF
#ifdef DEBUG_DAMPING
    WRITE(*,*) 'DAMPING DEBUG',n_damp
#endif

! *** Reads DIPOLE sections *****
    start_damp=0
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOM",i_rep_section=isec,c_val=atm_name,&
            error=error)
       apol_atm(start+isec) = atm_name
       CALL uppercase(apol_atm(start+isec))
       CALL section_vals_val_get(section,"APOL",i_rep_section=isec,r_val=apol(start+isec),&
            error=error)

       tmp_section => section_vals_get_subs_vals(section,"DAMPING",&
                                     i_rep_section=isec,error=error)
       CALL section_vals_get(tmp_section,n_repetition=tmp_damp,error=error)
       DO isec_damp=1,tmp_damp
          damping_list(start_damp+isec_damp)%atm_name1=apol_atm(start+isec)
          CALL section_vals_val_get(tmp_section,"ATOM",i_rep_section=isec_damp,&
               c_val=atm_name,&
               error=error)
          damping_list(start_damp+isec_damp)%atm_name2=atm_name
          CALL uppercase(damping_list(start_damp+isec_damp)%atm_name2)
          CALL section_vals_val_get(tmp_section,"TYPE",i_rep_section=isec_damp,&
               c_val=atm_name,ignore_required=.TRUE.,&
               error=error)
          damping_list(start_damp+isec_damp)%dtype=atm_name
          CALL uppercase(damping_list(start_damp+isec_damp)%dtype)

          CALL section_vals_val_get(tmp_section,"ORDER",i_rep_section=isec_damp,&
               i_val=damping_list(start_damp+isec_damp)%order,&
               ignore_required=.TRUE., error=error)
          CALL section_vals_val_get(tmp_section,"BIJ",i_rep_section=isec_damp,&
               r_val=damping_list(start_damp+isec_damp)%bij, error=error)
          CALL section_vals_val_get(tmp_section,"CIJ",i_rep_section=isec_damp,&
               r_val=damping_list(start_damp+isec_damp)%cij, error=error)
       END DO
       start_damp=start_damp+tmp_damp

    END DO

#ifdef DEBUG_DAMPING
    DO start_damp=1,n_damp
     WRITE(*,*) 'DAMPING DEBUG', TRIM(damping_list(start_damp)%atm_name1),&
                                 TRIM(damping_list(start_damp)%atm_name2),&
                                 TRIM(damping_list(start_damp)%dtype),&
                                 damping_list(start_damp)%order,&
                                 damping_list(start_damp)%bij,&
                                 damping_list(start_damp)%cij
    END DO
#endif

  END SUBROUTINE read_apol_section

! *****************************************************************************
!> \brief Reads the QUADRUPOLE POLARIZABILITY section
!> \param cpol_atm ...
!> \param cpol ...
!> \param section ...
!> \param start ...
!> \param error ...
!> \author Marcel Baer
! *****************************************************************************
  SUBROUTINE read_cpol_section(cpol_atm, cpol, section, start, error)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: cpol_atm
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cpol
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_cpol_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atm_name
    INTEGER                                  :: isec, n_items

    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"ATOM",i_rep_section=isec,c_val=atm_name,&
            error=error)
       cpol_atm(start+isec) = atm_name
       CALL uppercase(cpol_atm(start+isec))
       CALL section_vals_val_get(section,"CPOL",i_rep_section=isec,r_val=cpol(start+isec),&
            error=error)
    END DO
  END SUBROUTINE read_cpol_section


! *****************************************************************************
!> \brief Reads the SHELL section
!> \param shell_list ...
!> \param section ...
!> \param start ...
!> \param para_env ...
!> \param error ...
!> \author Marcella Iannuzzi
! *****************************************************************************
  SUBROUTINE read_shell_section(shell_list,section,start,para_env,error)

    TYPE(shell_p_type), DIMENSION(:), &
      POINTER                                :: shell_list
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_shell_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: atm_name
    INTEGER                                  :: i_rep, n_rep
    REAL(dp)                                 :: ccharge, cutoff, k, maxdist, &
                                                mfrac, scharge

    CALL section_vals_get(section,n_repetition=n_rep,error=error)

    DO i_rep = 1, n_rep
       CALL section_vals_val_get(section,"_SECTION_PARAMETERS_",&
            c_val=atm_name,i_rep_section=i_rep,error=error)
       CALL uppercase(atm_name)
       shell_list(start+i_rep)%atm_name = atm_name
       CALL section_vals_val_get(section,"CORE_CHARGE",i_rep_section=i_rep,r_val=ccharge,error=error)
       shell_list(start+i_rep)%shell%charge_core = ccharge
       CALL section_vals_val_get(section,"SHELL_CHARGE",i_rep_section=i_rep,r_val=scharge,error=error)
       shell_list(start+i_rep)%shell%charge_shell = scharge
       CALL section_vals_val_get(section,"MASS_FRACTION",i_rep_section=i_rep,r_val=mfrac,error=error)
       shell_list(start+i_rep)%shell%massfrac = mfrac
       CALL section_vals_val_get(section,"K2_SPRING",i_rep_section=i_rep,r_val=k,error=error)
       IF (k < 0.0_dp) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An invalid value was specified for the force constant k2 of the core-shell "//&
                            "spring potential",para_env)
       END IF
       shell_list(start+i_rep)%shell%k2_spring = k
       CALL section_vals_val_get(section,"K4_SPRING",i_rep_section=i_rep,r_val=k,error=error)
       IF (k < 0.0_dp) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "An invalid value was specified for the force constant k4 of the core-shell "//&
                            "spring potential",para_env)
       END IF
       shell_list(start+i_rep)%shell%k4_spring = k
       CALL section_vals_val_get(section,"MAX_DISTANCE",i_rep_section=i_rep,r_val=maxdist,error=error)
       shell_list(start+i_rep)%shell%max_dist = maxdist
       CALL section_vals_val_get(section,"SHELL_CUTOFF",i_rep_section=i_rep,r_val=cutoff,error=error)
       shell_list(start+i_rep)%shell%shell_cutoff=cutoff
    END DO

  END SUBROUTINE read_shell_section

! *****************************************************************************
!> \brief Reads the BONDS section
!> \param bond_kind ...
!> \param bond_a ...
!> \param bond_b ...
!> \param bond_k ...
!> \param bond_r0 ...
!> \param bond_cs ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \param bend_kind ...
!> \param bend_a ...
!> \param bend_b ...
!> \param bend_c ...
!> \param bend_k ...
!> \param bend_theta0 ...
!> \param bend_cb ...
!> \param bend_r012 ...
!> \param bend_r032 ...
!> \param bend_kbs12 ...
!> \param bend_kbs32 ...
!> \param bend_kss ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \brief ...
!> \param ub_kind ...
!> \param ub_a ...
!> \param ub_b ...
!> \param ub_c ...
!> \param ub_k ...
!> \param ub_r0 ...
!> \param section ...
!> \param start ...
!> \param error ...
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
!> \param torsion_kind ...
!> \param torsion_a ...
!> \param torsion_b ...
!> \param torsion_c ...
!> \param torsion_d ...
!> \param torsion_k ...
!> \param torsion_phi0 ...
!> \param torsion_m ...
!> \param section ...
!> \param start ...
!> \param error ...
!> \author teo
! *****************************************************************************
  SUBROUTINE read_torsions_section(torsion_kind,torsion_a, torsion_b, torsion_c, torsion_d, torsion_k,&
       torsion_phi0, torsion_m, section, start, error )
    INTEGER, DIMENSION(:), POINTER           :: torsion_kind
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: torsion_a, torsion_b, &
                                                torsion_c, torsion_d
    REAL(KIND=dp), DIMENSION(:), POINTER     :: torsion_k, torsion_phi0
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
       CALL section_vals_val_get(section,"PHI0",i_rep_section=isec,r_val=torsion_phi0(start+isec),error=error)
       CALL section_vals_val_get(section,"M",i_rep_section=isec,i_val=torsion_m(start+isec),error=error)
    END DO
  END SUBROUTINE read_torsions_section

! *****************************************************************************
!> \brief Reads the IMPROPER section
!> \param impr_kind ...
!> \param impr_a ...
!> \param impr_b ...
!> \param impr_c ...
!> \param impr_d ...
!> \param impr_k ...
!> \param impr_phi0 ...
!> \param section ...
!> \param start ...
!> \param error ...
!> \author louis vanduyfhuys
! *****************************************************************************
  SUBROUTINE read_improper_section(impr_kind,impr_a, impr_b, impr_c, impr_d, impr_k,&
       impr_phi0, section, start, error )
    INTEGER, DIMENSION(:), POINTER           :: impr_kind
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: impr_a, impr_b, impr_c, impr_d
    REAL(KIND=dp), DIMENSION(:), POINTER     :: impr_k, impr_phi0
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_improper_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY( atm_names)
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"KIND",i_rep_section=isec,i_val=impr_kind(start+isec),error=error)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       impr_a(start+isec) = atm_names(1)
       impr_b(start+isec) = atm_names(2)
       impr_c(start+isec) = atm_names(3)
       impr_d(start+isec) = atm_names(4)
       CALL uppercase(impr_a(start+isec))
       CALL uppercase(impr_b(start+isec))
       CALL uppercase(impr_c(start+isec))
       CALL uppercase(impr_d(start+isec))
       CALL section_vals_val_get(section,"K",i_rep_section=isec,r_val=impr_k(start+isec),error=error)
       CALL section_vals_val_get(section,"PHI0",i_rep_section=isec,r_val=impr_phi0(start+isec),error=error)
    END DO
  END SUBROUTINE read_improper_section

! *****************************************************************************
!> \brief Reads the OPBEND section
!> \param opbend_kind ...
!> \param opbend_a ...
!> \param opbend_b ...
!> \param opbend_c ...
!> \param opbend_d ...
!> \param opbend_k ...
!> \param opbend_phi0 ...
!> \param section ...
!> \param start ...
!> \param error ...
!> \author louis vanduyfhuys
! *****************************************************************************
  SUBROUTINE read_opbend_section(opbend_kind,opbend_a, opbend_b, opbend_c, opbend_d, opbend_k,&
       opbend_phi0, section, start, error )
    INTEGER, DIMENSION(:), POINTER           :: opbend_kind
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: opbend_a, opbend_b, opbend_c, &
                                                opbend_d
    REAL(KIND=dp), DIMENSION(:), POINTER     :: opbend_k, opbend_phi0
    TYPE(section_vals_type), POINTER         :: section
    INTEGER, INTENT(IN)                      :: start
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'read_opbend_section', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: atm_names
    INTEGER                                  :: isec, n_items
    LOGICAL                                  :: failure

    failure = .FALSE.
    NULLIFY( atm_names)
    CALL section_vals_get(section,n_repetition=n_items,error=error)
    DO isec = 1, n_items
       CALL section_vals_val_get(section,"KIND",i_rep_section=isec,i_val=opbend_kind(start+isec),error=error)
       CALL section_vals_val_get(section,"ATOMS",i_rep_section=isec,c_vals=atm_names,error=error)
       opbend_a(start+isec) = atm_names(1)
       opbend_b(start+isec) = atm_names(2)
       opbend_c(start+isec) = atm_names(3)
       opbend_d(start+isec) = atm_names(4)
       CALL uppercase(opbend_a(start+isec))
       CALL uppercase(opbend_b(start+isec))
       CALL uppercase(opbend_c(start+isec))
       CALL uppercase(opbend_d(start+isec))
       CALL section_vals_val_get(section,"K",i_rep_section=isec,r_val=opbend_k(start+isec),error=error)
       CALL section_vals_val_get(section,"PHI0",i_rep_section=isec,r_val=opbend_phi0(start+isec),error=error)
    END DO
  END SUBROUTINE read_opbend_section

! *****************************************************************************
!> \brief Reads the force_field input section
!> \param ff_type ...
!> \param para_env ...
!> \param mm_section ...
!> \param error ...
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
!> \brief reads EAM potential from library
!> \param eam ...
!> \param para_env ...
!> \param mm_section ...
!> \param error ...
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
