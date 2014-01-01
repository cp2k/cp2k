!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Working with the  semi empirical parameter types.
!> \author JGH (14.08.2004)
! *****************************************************************************
MODULE semi_empirical_utils
  USE basis_set_types,                 ONLY: allocate_sto_basis_set,&
                                             create_gto_from_sto_basis,&
                                             gto_basis_set_type,&
                                             set_sto_basis_set
  USE cell_types,                      ONLY: cell_type,&
                                             plane_distance
  USE cp_control_types,                ONLY: semi_empirical_control_type
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_method_am1, do_method_mndo, do_method_mndod, do_method_pchg, &
       do_method_pdg, do_method_pm3, do_method_pm6, do_method_pnnl, &
       do_method_rm1
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE semi_empirical_mpole_methods,    ONLY: semi_empirical_mpole_p_setup
  USE semi_empirical_par_utils,        ONLY: get_se_basis,&
                                             setup_1c_2el_int
  USE semi_empirical_parameters,       ONLY: am1_default_parameter,&
                                             mndo_default_parameter,&
                                             pcharge_default_parameter,&
                                             pdg_default_parameter,&
                                             pm3_default_parameter,&
                                             pm6_default_parameter,&
                                             pnnl_default_parameter,&
                                             rm1_default_parameter
  USE semi_empirical_types,            ONLY: se_taper_type,&
                                             semi_empirical_type
  USE termination,                     ONLY: stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_utils'

  PUBLIC :: init_se_param, se_param_set_default, get_se_type,&
            initialize_se_taper, finalize_se_taper, se_cutoff_compatible

CONTAINS
! *****************************************************************************
!> \brief  Reset cutoffs trying to be somehow a bit smarter
!> \author Teodoro Laino [tlaino] - 03.2009
! *****************************************************************************
  SUBROUTINE se_cutoff_compatible(se_control, se_section, cell, output_unit, error)
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    TYPE(section_vals_type), POINTER         :: se_section
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_cutoff_compatible', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: explicit1, explicit2, failure
    REAL(KIND=dp)                            :: rc

    failure = .FALSE.
    ! Coulomb Cutoff Taper
    CALL section_vals_val_get(se_section,"COULOMB%CUTOFF",explicit=explicit1,error=error)
    CALL section_vals_val_get(se_section,"COULOMB%RC_TAPER",explicit=explicit2,error=error)
    IF ((.NOT.explicit1).AND.se_control%do_ewald_gks) THEN
       rc = MAX(0.5*plane_distance(1,0,0,cell),&
                0.5*plane_distance(0,1,0,cell),&
                0.5*plane_distance(0,0,1,cell))
       IF (rc /= se_control%cutoff_cou) THEN
          IF (output_unit>0) THEN
             WRITE(output_unit,*)
             WRITE(output_unit,'(A,T37,A)')     " SEMIEMPIRICAL|",&
                  " Coulomb Integral cutoff/taper was redefined"
             WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Old value [a.u.]",&
                  se_control%cutoff_cou
             WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| New value [a.u.]",rc
             WRITE(output_unit,*)
          END IF
       END IF
                           se_control%cutoff_cou = rc
       IF (.NOT.explicit2) se_control%taper_cou  = rc
    ELSE IF ((.NOT.explicit1).AND.(ALL(cell%perd==0))) THEN
       rc = MAX(plane_distance(1,0,0,cell),&
                plane_distance(0,1,0,cell),&
                plane_distance(0,0,1,cell))
       IF (rc /= se_control%cutoff_cou) THEN
          IF (output_unit>0) THEN
             WRITE(output_unit,*)
             WRITE(output_unit,'(A,T37,A)')     " SEMIEMPIRICAL|",&
                  " Coulomb Integral cutoff/taper was redefined"
             WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Old value [a.u.]",&
                  se_control%cutoff_cou
             WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| New value [a.u.]",rc
             WRITE(output_unit,*)
          END IF
       END IF
                           se_control%cutoff_cou = rc
       IF (.NOT.explicit2) se_control%taper_cou  = rc
    END IF
    IF (output_unit>0) THEN
       WRITE(output_unit,*)
       WRITE(output_unit,'(A,T44,A)')     " SEMIEMPIRICAL|",&
            " Coulomb Integral cutoff/taper values"
       WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Cutoff [a.u.]",&
            se_control%cutoff_cou
       WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Taper  [a.u.]",&
            se_control%taper_cou
       WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Range  [a.u.]",&
            se_control%range_cou
       WRITE(output_unit,*)
    END IF
    ! Exchange Cutoff Taper
    CALL section_vals_val_get(se_section,"EXCHANGE%CUTOFF",explicit=explicit1,error=error)
    CALL section_vals_val_get(se_section,"EXCHANGE%RC_TAPER",explicit=explicit2,error=error)
    rc = se_control%cutoff_exc
    IF (.NOT.explicit1) THEN
       rc = MIN(rc,MAX(0.25_dp*plane_distance(1,0,0,cell),&
                       0.25_dp*plane_distance(0,1,0,cell),&
                       0.25_dp*plane_distance(0,0,1,cell)))

       IF (rc /= se_control%cutoff_exc) THEN
          IF (output_unit>0) THEN
             WRITE(output_unit,*)
             WRITE(output_unit,'(A,T36,A)')     " SEMIEMPIRICAL|",&
                  " Exchange Integral cutoff/taper was redefined"
             WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Default value [a.u.]",&
                  se_control%cutoff_exc
             WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| New value [a.u.]",rc
             WRITE(output_unit,*)
          END IF
       END IF
    END IF
                        se_control%cutoff_exc = rc
    IF (.NOT.explicit2) se_control%taper_exc  = rc

    IF (output_unit>0) THEN
       WRITE(output_unit,*)
       WRITE(output_unit,'(A,T43,A)')     " SEMIEMPIRICAL|",&
            " Exchange Integral cutoff/taper values"
       WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Cutoff [a.u.]",&
            se_control%cutoff_exc
       WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Taper  [a.u.]",&
            se_control%taper_exc
       WRITE(output_unit,'(A,T71,F10.3)') " SEMIEMPIRICAL| Range  [a.u.]",&
            se_control%range_exc
       WRITE(output_unit,*)
    END IF

  END SUBROUTINE se_cutoff_compatible

! *****************************************************************************
!> \brief  Initializes the semi-empirical taper for a chunk calculation
!> \author Teodoro Laino [tlaino] - 03.2009
! *****************************************************************************
  SUBROUTINE initialize_se_taper(se_taper, coulomb, exchange, lr_corr, error)
    TYPE(se_taper_type), POINTER             :: se_taper
    LOGICAL, INTENT(IN), OPTIONAL            :: coulomb, exchange, lr_corr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'initialize_se_taper', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check, failure, l_coulomb, &
                                                l_exchange, l_lrc

    failure = .FALSE.
    check = .NOT.ASSOCIATED(se_taper%taper)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    l_coulomb  = .FALSE.
    l_exchange = .FALSE.
    l_lrc      = .FALSE.
    IF (PRESENT(coulomb))  l_coulomb  = coulomb
    IF (PRESENT(exchange)) l_exchange = exchange
    IF (PRESENT(lr_corr))  l_lrc      = lr_corr
    IF (l_coulomb)  THEN
       check = (.NOT.l_exchange).AND.(.NOT.l_lrc)
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       se_taper%taper => se_taper%taper_cou
    END IF
    IF (l_exchange) THEN
       check = (.NOT.l_coulomb).AND.(.NOT.l_lrc)
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       se_taper%taper => se_taper%taper_exc
    END IF
    IF (l_lrc)      THEN
       check = (.NOT.l_coulomb).AND.(.NOT.l_exchange)
       CPPrecondition(check,cp_failure_level,routineP,error,failure)
       se_taper%taper => se_taper%taper_lrc
    END IF
  END SUBROUTINE initialize_se_taper

! *****************************************************************************
!> \brief  Finalizes the semi-empirical taper for a chunk calculation
!> \author Teodoro Laino [tlaino] - 03.2009
! *****************************************************************************
  SUBROUTINE finalize_se_taper(se_taper, error)
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'finalize_se_taper', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: check, failure

    failure = .FALSE.
    check = ASSOCIATED(se_taper%taper)
    CPPrecondition(check,cp_failure_level,routineP,error,failure)
    NULLIFY(se_taper%taper)
  END SUBROUTINE finalize_se_taper

! *****************************************************************************
!> \brief Initialize semi_empirical type
! *****************************************************************************
  SUBROUTINE init_se_param(sep,orb_basis_set,ngauss,error)
    TYPE(semi_empirical_type), POINTER       :: sep
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    INTEGER, INTENT(IN)                      :: ngauss
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_se_param', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: symbol
    INTEGER                                  :: l, nshell, stat
    INTEGER, DIMENSION(:), POINTER           :: lq, nq
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: zet

    failure = .FALSE.
    IF (ASSOCIATED(sep)) THEN
       CALL allocate_sto_basis_set (sep%basis,error)
       nshell = 0
       IF (sep%natorb ==  1) nshell=1
       IF (sep%natorb ==  4) nshell=2
       IF (sep%natorb ==  9) nshell=3
       ALLOCATE (nq(0:3),lq(0:3),zet(0:3),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       ALLOCATE (symbol(0:3),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       symbol=""
       nq=0
       lq=0
       zet=0._dp
       DO l=0,nshell-1
          nq(l)=get_se_basis(sep,l,error)
          lq(l)=l
          zet(l)=sep%sto_exponents(l)
          IF(l==0) WRITE(symbol(0),'(I1,A1)') nq(l),"S"
          IF(l==1) WRITE(symbol(1),'(I1,A1)') nq(l),"P"
          IF(l==2) WRITE(symbol(2),'(I1,A1)') nq(l),"D"
       END DO

       IF (nshell > 0) THEN
          sep%ngauss = ngauss
          CALL set_sto_basis_set(sep%basis,name=sep%name,nshell=nshell,symbol=symbol,&
               nq=nq,lq=lq,zet=zet)
          CALL create_gto_from_sto_basis(sep%basis,orb_basis_set,sep%ngauss,error=error)
       END IF

       DEALLOCATE (nq,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (lq,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (zet,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE (symbol,STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer sep is not associated")
    END IF

  END SUBROUTINE init_se_param

! *****************************************************************************
!> \brief Initialize parameter for a semi_empirival type
! *****************************************************************************
  SUBROUTINE se_param_set_default(sep,z,method,error)

    TYPE(semi_empirical_type), POINTER       :: sep
    INTEGER, INTENT(IN)                      :: z, method
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_param_set_default', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(sep)) THEN
       IF (z<0) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Atomic number < 0")
       END IF
       SELECT CASE (method)
       CASE (do_method_am1)
          CALL am1_default_parameter(sep,z,error)
       CASE (do_method_rm1)
          CALL rm1_default_parameter(sep,z,error)
       CASE (do_method_pm3)
          CALL pm3_default_parameter(sep,z,error)
       CASE (do_method_pm6)
          CALL pm6_default_parameter(sep,z,error)
       CASE (do_method_pdg)
          CALL pdg_default_parameter(sep,z,error)
       CASE (do_method_mndo)
          CALL mndo_default_parameter(sep,z,do_method_mndo, error)
       CASE (do_method_mndod)
          CALL mndo_default_parameter(sep,z,do_method_mndod,error)
       CASE (do_method_pnnl)
          CALL pnnl_default_parameter(sep,z,error)
       CASE (do_method_pchg)
          CALL pcharge_default_parameter(sep,z,error)
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
               "Semiempirical method unknown")
       END SELECT
    ELSE
       CALL stop_program(routineN,moduleN,__LINE__,&
            "The pointer sep is not associated")
    END IF

    ! Check if the element has been defined..
    CALL cp_assert(sep%defined,cp_fatal_level,cp_assertion_failed,routineP,&
         "Semiempirical type ("//TRIM(sep%name)//") cannot be defined for "//&
         "the requested parameterization."//&
CPSourceFileRef,&
         only_ionode=.TRUE.)

    ! Fill 1 center - 2 electron integrals
    CALL setup_1c_2el_int(sep, error)

    ! Fill multipolar expansion of atomic orbitals charge distributions
    CALL semi_empirical_mpole_p_setup(sep%w_mpole, sep, method, error)

    ! Get the value of the size of CORE integral array
    sep%core_size = 0
    IF (sep%natorb>0) sep%core_size = 1
    IF (sep%natorb>1) sep%core_size = 4
    IF (sep%dorb)     sep%core_size =10

    ! Get size of the all possible combinations of atomic orbitals
    sep%atm_int_size = (sep%natorb+1)*sep%natorb/2

  END SUBROUTINE se_param_set_default

! *****************************************************************************
!> \brief Gives back the unique semi_empirical METHOD type
! *****************************************************************************
  FUNCTION get_se_type(se_method) RESULT(se_type)

    INTEGER, INTENT(IN)                      :: se_method
    INTEGER                                  :: se_type

    CHARACTER(len=*), PARAMETER :: routineN = 'get_se_type', &
      routineP = moduleN//':'//routineN

    SELECT CASE (se_method)
    CASE DEFAULT
       se_type=se_method
    CASE (do_method_am1,do_method_rm1)
       se_type=do_method_am1
    END SELECT

  END FUNCTION get_se_type

END MODULE semi_empirical_utils

