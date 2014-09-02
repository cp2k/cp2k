!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set disperson types for DFT calculations
!> \author JGH (04.2014)
! *****************************************************************************
MODULE qs_dispersion_utils

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_constants,                 ONLY: vdw_nl_DRSLL,&
                                             vdw_nl_LMKLL,&
                                             vdw_nl_RVV10,&
                                             vdw_pairpot_dftd2,&
                                             vdw_pairpot_dftd3,&
                                             xc_vdw_fun_nonloc,&
                                             xc_vdw_fun_pairpot
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: bohr,&
                                             kjmol
  USE qs_dispersion_pairpot,           ONLY: qs_scaling_dftd3,&
                                             qs_scaling_init
  USE qs_dispersion_types,             ONLY: qs_atom_dispersion_type,&
                                             qs_dispersion_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_utils'

  PUBLIC :: qs_dispersion_env_set, qs_write_dispersion

! *****************************************************************************
  CONTAINS
! *****************************************************************************
!> \brief ...
!> \param dispersion_env ...
!> \param xc_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_dispersion_env_set(dispersion_env,xc_section,error)
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_dispersion_env_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: exfun, explicit, failure
    REAL(dp), POINTER                        :: params(:), scal(:)
    TYPE(section_vals_type), POINTER         :: nl_section, pp_section, &
                                                vdw_section, xc_fun_section

    CPPrecondition(ASSOCIATED(dispersion_env),cp_failure_level,routineP,error,failure)
    failure = .FALSE.

    ! set general defaults
    dispersion_env%doabc=.FALSE.
    dispersion_env%c9cnst=.FALSE.
    dispersion_env%lrc=.FALSE.
    dispersion_env%verbose=.FALSE.
    NULLIFY(dispersion_env%c6ab,dispersion_env%maxci,dispersion_env%r0ab,dispersion_env%rcov,&
            dispersion_env%r2r4,dispersion_env%cn,dispersion_env%cnkind,dispersion_env%cnlist)
    NULLIFY(dispersion_env%q_mesh,dispersion_env%kernel,dispersion_env%d2phi_dk2,&
            dispersion_env%d2y_dx2)
    NULLIFY(dispersion_env%sab_vdw,dispersion_env%sab_cn)
    NULLIFY(dispersion_env%dftd_section)
    NULLIFY(vdw_section,xc_fun_section)
    vdw_section =>  section_vals_get_subs_vals(xc_section,"vdw_potential",error=error)
    xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",error=error)
    CALL section_vals_val_get(vdw_section, "POTENTIAL_TYPE", i_val=dispersion_env%type, error=error)
    IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN
       NULLIFY(pp_section)
       pp_section =>  section_vals_get_subs_vals(vdw_section,"PAIR_POTENTIAL",error=error)
       CALL section_vals_val_get(pp_section, "VERBOSE_OUTPUT", l_val=dispersion_env%verbose, error=error)
       CALL section_vals_val_get(pp_section, "TYPE", i_val=dispersion_env%pp_type, error=error)
       IF ( dispersion_env%pp_type == vdw_pairpot_dftd2 ) THEN
          ! functional parameters for Grimme D2 type
          CALL section_vals_val_get(pp_section, "EXP_PRE", r_val=dispersion_env%exp_pre, error=error)
          CALL section_vals_val_get(pp_section, "SCALING" ,explicit=explicit,error=error)
          IF ( .NOT. explicit ) THEN
            CALL section_vals_val_get(pp_section, "REFERENCE_FUNCTIONAL" ,explicit=exfun,error=error)
            CPPostcondition(exfun,cp_failure_level,routineP,error,failure)
            CALL qs_scaling_init(dispersion_env%scaling,vdw_section,error)
          ELSE
            CALL section_vals_val_get(pp_section, "SCALING", r_val=dispersion_env%scaling, error=error)
          END IF
       ELSE
          dispersion_env%exp_pre=0._dp
          dispersion_env%scaling=0._dp
       END IF
       IF ( dispersion_env%pp_type == vdw_pairpot_dftd3 ) THEN
          ! functional parameters for Grimme DFT-D3 type
          CALL section_vals_val_get(pp_section, "EPS_CN", r_val=dispersion_env%eps_cn, error=error)
          CALL section_vals_val_get(pp_section, "D3_SCALING" ,explicit=explicit,error=error)
          CALL section_vals_val_get(pp_section, "CALCULATE_C9_TERM" ,l_val=dispersion_env%doabc,error=error)
          CALL section_vals_val_get(pp_section, "REFERENCE_C9_TERM" ,l_val=dispersion_env%c9cnst,error=error)
          CALL section_vals_val_get(pp_section, "LONG_RANGE_CORRECTION" ,l_val=dispersion_env%lrc,error=error)
          IF ( .NOT. explicit ) THEN
             CALL section_vals_val_get(pp_section, "REFERENCE_FUNCTIONAL" ,explicit=exfun,error=error)
             CPPostcondition(exfun,cp_failure_level,routineP,error,failure)
             CALL qs_scaling_dftd3(dispersion_env%s6,dispersion_env%sr6,dispersion_env%s8,vdw_section,error)
          ELSE
             CALL section_vals_val_get(pp_section, "D3_SCALING", r_vals=scal, error=error)
             dispersion_env%s6  =scal(1)
             dispersion_env%sr6 =scal(2)
             dispersion_env%s8  =scal(3)
          END IF
       ELSE
          dispersion_env%s6=0._dp
          dispersion_env%sr6=0._dp
          dispersion_env%s8=0._dp
          dispersion_env%eps_cn=0._dp
       END IF
       CALL section_vals_val_get(pp_section, "R_CUTOFF", r_val=dispersion_env%rc_disp, error=error)
       CALL section_vals_val_get(pp_section,"PARAMETER_FILE_NAME",&
            c_val=dispersion_env%parameter_file_name,error=error)
       ! set DFTD section for output handling
       dispersion_env%dftd_section => pp_section
    ELSE IF ( dispersion_env%type == xc_vdw_fun_nonloc ) THEN
       NULLIFY(nl_section)
       nl_section =>  section_vals_get_subs_vals(vdw_section,"NON_LOCAL",error=error)
       CALL section_vals_val_get(nl_section, "VERBOSE_OUTPUT", l_val=dispersion_env%verbose, error=error)
       CALL section_vals_val_get(nl_section,"KERNEL_FILE_NAME",&
            c_val=dispersion_env%kernel_file_name,error=error)
       CALL section_vals_val_get(nl_section, "TYPE", i_val=dispersion_env%nl_type, error=error)
       CALL section_vals_val_get(nl_section, "CUTOFF", r_val=dispersion_env%pw_cutoff, error=error)
       dispersion_env%pw_cutoff = 0.5_dp * dispersion_env%pw_cutoff
       CALL section_vals_val_get(nl_section, "PARAMETERS", r_vals=params, error=error)
       dispersion_env%b_value = params(1)
       dispersion_env%c_value = params(2)
    END IF
  END SUBROUTINE qs_dispersion_env_set

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param ounit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_write_dispersion(qs_env,dispersion_env,ounit,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_dispersion_type), POINTER        :: dispersion_env
    INTEGER, INTENT(in), OPTIONAL            :: ounit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_write_dispersion', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: symbol
    INTEGER                                  :: i, ikind, nkind, output_unit
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_atom_dispersion_type), POINTER   :: disp
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: dft_section

    failure=.FALSE.

    IF(PRESENT(ounit)) THEN
       output_unit = ounit
    ELSE
       NULLIFY(logger)
       logger => cp_error_get_logger(error)

       dft_section => section_vals_get_subs_vals(qs_env%input,"DFT",error=error)
       output_unit = cp_print_key_unit_nr(logger,dft_section,&
            "PRINT%DFT_CONTROL_PARAMETERS",extension=".Log",error=error)
    END IF

    IF (output_unit>0) THEN
       ! vdW type specific output
       IF ( dispersion_env%type == xc_vdw_fun_pairpot ) THEN
          WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T67,'Pair Potential')")
          ! Pair potentials
          IF ( dispersion_env%pp_type == vdw_pairpot_dftd2 ) THEN
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Potential Form: S. Grimme, JCC 27: 1787 (2006)')")
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Cutoff Radius [Bohr]:',T73,F8.2)") dispersion_env%rc_disp
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Scaling Factor:',T73,F8.4)") dispersion_env%scaling
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T35,'Exp Prefactor for Damping:',T73,F8.1)") dispersion_env%exp_pre
             CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,error=error)
             nkind  = SIZE(atomic_kind_set)
             DO ikind = 1, nkind
                CALL get_atomic_kind(atomic_kind_set(ikind),element_symbol=symbol)
                CALL get_qs_kind(qs_kind_set(ikind),dispersion=disp)
                IF ( disp%defined ) THEN
                   WRITE (output_unit, fmt="(' vdW PARAMETER| ',T18,'Atom=',A2,"//&
                        "T28,'C6[J*nm^6*mol^-1]=',F8.4,T63,'r(vdW)[A]=',F8.4)") &
                        symbol, disp%c6/(1000._dp*bohr**6/kjmol), disp%vdw_radii/bohr
                ELSE
                   WRITE (output_unit, fmt="(' vdW PARAMETER| ',T20,'Atom=',A2,T70,'not defined')")
                END IF
             END DO
          ELSE IF ( dispersion_env%pp_type == vdw_pairpot_dftd3 ) THEN
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Potential Form: S. Grimme et al, JCP 132: 154104 (2010)')")
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Cutoff Radius [Bohr]:',T73,F8.2)") dispersion_env%rc_disp
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'s6 Scaling Factor:',T73,F8.4)") dispersion_env%s6
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'sr6 Scaling Factor:',T73,F8.4)") dispersion_env%sr6
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'s8 Scaling Factor:',T73,F8.4)") dispersion_env%s8
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T26,'Cutoff for CN calculation:',T69,E12.4)") dispersion_env%eps_cn
          END IF
       ELSE IF ( dispersion_env%type == xc_vdw_fun_nonloc ) THEN
          WRITE (output_unit, fmt="(' vdW POTENTIAL| ',T61,'Non-local Functional')")
          WRITE (output_unit, &
            fmt="(' vdW POTENTIAL| ','Implementation: G. Roman-Perez, J. Soler, PRL 103: 096102 (2009)')")
          WRITE (output_unit, &
            fmt="(' vdW POTENTIAL| ',T38,' T. Thonhauser et al, PRB 76: 125112 (2007)')")
          WRITE (output_unit, &
            fmt="(' vdW POTENTIAL| ',T22,' R. Sabatini et al, J.Phys:Condens Matter 24: 424209 (2012)')")
          WRITE (output_unit, &
            fmt="(' vdW POTENTIAL| ',T16,' Based on QE implementation by Brian Kolb, Timo Thonhauser (2009)')")
          SELECT CASE (dispersion_env%nl_type)
             CASE DEFAULT
                ! unknown functional
                CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
             CASE (vdw_nl_DRSLL)
                WRITE (output_unit, &
                 fmt="(' vdW POTENTIAL| ','DRSLL Functional:           M. Dion et al, PRL 92: 246401 (2004)')")
             CASE (vdw_nl_LMKLL)
                WRITE (output_unit, &
                 fmt="(' vdW POTENTIAL| ','LMKLL Functional:            K. Lee et al, PRB 82: 081101 (2010)')")
             CASE (vdw_nl_RVV10)
                WRITE (output_unit, &
                 fmt="(' vdW POTENTIAL| ','RVV10 Functional:    R. Sabatini et al, PRB 87: 041108(R) (2013)')")
          END SELECT
          IF ( dispersion_env%verbose ) THEN
             WRITE (output_unit,&
              fmt="(' vdW POTENTIAL| ','         Carrying out vdW-DF run using the following parameters:')")
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ','Nqs =',I8,'        Nr_points =',I8,'       r_max =',F10.3)")&
                dispersion_env%nqs,dispersion_env%nr_points,dispersion_env%r_max
             WRITE (output_unit, fmt="(' vdW POTENTIAL| ','q_mesh =')")
             WRITE (output_unit, fmt="(8X,4F18.8)" ) (dispersion_env%q_mesh(i), i=1, dispersion_env%nqs)
             WRITE (output_unit,&
              fmt="(' vdW POTENTIAL| ','Density cutoff for comvolution:',T68,F10.1,' Ry')") &
              2.0_dp*dispersion_env%pw_cutoff
          END IF
       END IF
    END IF
    IF(.NOT.PRESENT(ounit)) THEN
       CALL cp_print_key_finished_output(output_unit,logger,dft_section,&
            "PRINT%DFT_CONTROL_PARAMETERS",error=error)
    END IF

  END SUBROUTINE qs_write_dispersion

! *****************************************************************************

END MODULE qs_dispersion_utils

