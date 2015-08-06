!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Handles the type to compute averages during an MD
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
MODULE averages_types
  USE cell_types,                      ONLY: cell_type
  USE colvar_utils,                    ONLY: get_clv_force,&
                                             number_of_colvar
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE force_env_types,                 ONLY: force_env_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_remove_values,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE md_ener_types,                   ONLY: md_ener_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_create,&
                                             virial_release,&
                                             virial_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************
  TYPE average_quantities_type
     INTEGER                                       :: id_nr, ref_count, itimes_start
     LOGICAL                                       :: do_averages
     TYPE(section_vals_type), POINTER              :: averages_section
     ! Real Scalar Quantities
     REAL(KIND=dp)                                 :: avetemp, avepot, avekin,&
                                                      avevol, aveca, avecb, avecc
     REAL(KIND=dp)                                 :: avetemp_baro, avehugoniot, avecpu
     REAL(KIND=dp)                                 :: aveal, avebe, avega, avepress,&
                                                      avekinc, avetempc, avepxx
     REAL(KIND=dp)                                 :: avetemp_qm, avekin_qm, econs
     ! Virial
     TYPE(virial_type), POINTER                    :: virial
     ! Colvar
     REAL(KIND=dp), POINTER, DIMENSION(:)          :: avecolvar
     REAL(KIND=dp), POINTER, DIMENSION(:)          :: aveMmatrix
  END TYPE average_quantities_type

! *****************************************************************************
  INTERFACE get_averages
     MODULE PROCEDURE get_averages_rs, get_averages_rv, get_averages_rm
  END INTERFACE get_averages

! *** Public subroutines and data types ***
  PUBLIC :: average_quantities_type, create_averages,  release_averages,&
            retain_averages, compute_averages

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'averages_types'
  INTEGER, SAVE, PRIVATE :: last_avg_env_id=0

CONTAINS

! *****************************************************************************
!> \brief Creates averages environment
!> \param averages ...
!> \param averages_section ...
!> \param virial_avg ...
!> \param force_env ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE create_averages(averages, averages_section, virial_avg, force_env, error)
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(section_vals_type), POINTER         :: averages_section
    LOGICAL, INTENT(IN), OPTIONAL            :: virial_avg
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'create_averages', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nint, stat
    LOGICAL                                  :: do_colvars, failure

    failure = .FALSE.
    ALLOCATE (averages, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    NULLIFY(averages%virial)
    NULLIFY(averages%avecolvar)
    NULLIFY(averages%aveMmatrix)
    ! Point to the averages section
    averages%averages_section => averages_section
    ! Initialize averages
    last_avg_env_id=last_avg_env_id+1
    averages%id_nr=last_avg_env_id
    averages%ref_count      =  1
    averages%itimes_start   = -1
    averages%avetemp        = 0.0_dp
    averages%avepot         = 0.0_dp
    averages%avekin         = 0.0_dp
    averages%avevol         = 0.0_dp
    averages%aveca          = 0.0_dp
    averages%avecb          = 0.0_dp
    averages%avecc          = 0.0_dp
    averages%avetemp_baro   = 0.0_dp
    averages%avehugoniot    = 0.0_dp
    averages%avecpu         = 0.0_dp
    averages%aveal          = 0.0_dp
    averages%avebe          = 0.0_dp
    averages%avega          = 0.0_dp
    averages%avepress       = 0.0_dp
    averages%avekinc        = 0.0_dp
    averages%avetempc       = 0.0_dp
    averages%avepxx         = 0.0_dp
    averages%avetemp_qm     = 0.0_dp
    averages%avekin_qm      = 0.0_dp
    averages%econs          = 0.0_dp
    CALL section_vals_val_get(averages_section,"_SECTION_PARAMETERS_",l_val=averages%do_averages,&
         error=error)
    IF (averages%do_averages) THEN
       ! Setup Virial if requested
       IF (PRESENT(virial_avg)) THEN
          IF (virial_avg) CALL virial_create(averages%virial, error)
       END IF
       CALL section_vals_val_get(averages_section,"AVERAGE_COLVAR",l_val=do_colvars,error=error)
       ! Total number of COLVARs
       nint = 0
       IF (do_colvars) nint = number_of_colvar(force_env, error=error)
       ALLOCATE(averages%avecolvar(nint), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(averages%aveMmatrix(nint*nint), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, nint
          averages%avecolvar(i)         = 0.0_dp
       END DO
       DO i = 1, nint*nint
          averages%aveMmatrix(i)        = 0.0_dp
       END DO
    END IF
  END SUBROUTINE create_averages

! *****************************************************************************
!> \brief retains the given averages env
!> \param averages ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE retain_averages(averages, error)
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'retain_averages', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(averages),cp_failure_level,routineP,error,failure)
    CPPrecondition(averages%ref_count>0,cp_failure_level,routineP,error,failure)
    averages%ref_count=averages%ref_count+1
  END SUBROUTINE retain_averages

! *****************************************************************************
!> \brief releases the given averages env
!> \param averages ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE release_averages(averages, error)
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'release_averages', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: work_section

    failure=.FALSE.
    IF (ASSOCIATED(averages)) THEN
       CPPrecondition(averages%ref_count>0,cp_failure_level,routineP,error,failure)
       averages%ref_count=averages%ref_count-1
       IF (averages%ref_count==0) THEN
          CALL virial_release(averages%virial, error)
          IF (ASSOCIATED(averages%avecolvar)) THEN
             DEALLOCATE(averages%avecolvar, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          IF (ASSOCIATED(averages%aveMmatrix)) THEN
             DEALLOCATE(averages%aveMmatrix, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          END IF
          ! Removes restart values from the corresponding restart section..
          work_section => section_vals_get_subs_vals(averages%averages_section,"RESTART_AVERAGES",error=error)
          CALL section_vals_remove_values(work_section, error)
          NULLIFY(averages%averages_section)
          !
          DEALLOCATE(averages,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

  END SUBROUTINE release_averages

! *****************************************************************************
!> \brief computes the averages
!> \param averages ...
!> \param force_env ...
!> \param md_ener ...
!> \param cell ...
!> \param virial ...
!> \param pv_scalar ...
!> \param pv_xx ...
!> \param used_time ...
!> \param hugoniot ...
!> \param abc ...
!> \param cell_angle ...
!> \param nat ...
!> \param itimes ...
!> \param time ...
!> \param my_pos ...
!> \param my_act ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE compute_averages (averages, force_env, md_ener, cell,  virial, &
       pv_scalar, pv_xx, used_time, hugoniot, abc, cell_angle, nat, itimes, &
       time, my_pos, my_act, error)
    TYPE(average_quantities_type), POINTER   :: averages
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(md_ener_type), POINTER              :: md_ener
    TYPE(cell_type), POINTER                 :: cell
    TYPE(virial_type), POINTER               :: virial
    REAL(KIND=dp), INTENT(IN)                :: pv_scalar, pv_xx
    REAL(KIND=dp), POINTER                   :: used_time
    REAL(KIND=dp), INTENT(IN)                :: hugoniot
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: abc, cell_angle
    INTEGER, INTENT(IN)                      :: nat, itimes
    REAL(KIND=dp), INTENT(IN)                :: time
    CHARACTER(LEN=default_string_length), &
      INTENT(IN)                             :: my_pos, my_act
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_averages', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: ctmp
    INTEGER                                  :: delta_t, handle, i, nint, &
                                                output_unit, stat
    LOGICAL                                  :: failure, restart_averages
    REAL(KIND=dp)                            :: start_time
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cvalues, Mmatrix, tmp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: restart_section

    failure = .FALSE.
    CALL timeset(routineN,handle)
    CALL section_vals_val_get(averages%averages_section,"ACQUISITION_START_TIME",&
         r_val=start_time, error=error)
    IF (averages%do_averages) THEN
       NULLIFY(cvalues, Mmatrix, logger)
       logger => cp_error_get_logger(error)
       ! Determine the nr. of internal colvar (if any/requested)
       nint = 0
       IF (ASSOCIATED(averages%avecolvar)) nint = SIZE(averages%avecolvar)

       ! Evaluate averages if we collected enough statistics (user defined)
       IF (time>=start_time) THEN

          ! Handling properly the restart
          IF (averages%itimes_start==-1) THEN
             restart_section => section_vals_get_subs_vals(averages%averages_section,"RESTART_AVERAGES",error=error)
             CALL section_vals_get(restart_section, explicit=restart_averages, error=error)
             IF (restart_averages) THEN
                CALL section_vals_val_get(restart_section,"ITIMES_START",i_val=averages%itimes_start,error=error)
                CALL section_vals_val_get(restart_section,"AVECPU",r_val=averages%avecpu,error=error)
                CALL section_vals_val_get(restart_section,"AVEHUGONIOT",r_val=averages%avehugoniot,error=error)
                CALL section_vals_val_get(restart_section,"AVETEMP_BARO",r_val=averages%avetemp_baro,error=error)
                CALL section_vals_val_get(restart_section,"AVEPOT",r_val=averages%avepot,error=error)
                CALL section_vals_val_get(restart_section,"AVEKIN",r_val=averages%avekin,error=error)
                CALL section_vals_val_get(restart_section,"AVETEMP",r_val=averages%avetemp,error=error)
                CALL section_vals_val_get(restart_section,"AVEKIN_QM",r_val=averages%avekin_qm,error=error)
                CALL section_vals_val_get(restart_section,"AVETEMP_QM",r_val=averages%avetemp_qm,error=error)
                CALL section_vals_val_get(restart_section,"AVEVOL",r_val=averages%avevol,error=error)
                CALL section_vals_val_get(restart_section,"AVECELL_A",r_val=averages%aveca,error=error)
                CALL section_vals_val_get(restart_section,"AVECELL_B",r_val=averages%avecb,error=error)
                CALL section_vals_val_get(restart_section,"AVECELL_C",r_val=averages%avecc,error=error)
                CALL section_vals_val_get(restart_section,"AVEALPHA",r_val=averages%aveal,error=error)
                CALL section_vals_val_get(restart_section,"AVEBETA",r_val=averages%avebe,error=error)
                CALL section_vals_val_get(restart_section,"AVEGAMMA",r_val=averages%avega,error=error)
                CALL section_vals_val_get(restart_section,"AVE_ECONS",r_val=averages%econs,error=error)
                ! Virial
                IF (virial%pv_availability) THEN
                   CALL section_vals_val_get(restart_section,"AVE_PRESS",r_val=averages%avepress,error=error)
                   CALL section_vals_val_get(restart_section,"AVE_PXX",r_val=averages%avepxx,error=error)
                   IF (ASSOCIATED(averages%virial)) THEN
                      CALL section_vals_val_get(restart_section,"AVE_PV_TOT",r_vals=tmp,error=error)
                      averages%virial%pv_total = RESHAPE(tmp,(/3,3/))
                      CALL section_vals_val_get(restart_section,"AVE_PV_VIR",r_vals=tmp,error=error)
                      averages%virial%pv_virial = RESHAPE(tmp,(/3,3/))
                      CALL section_vals_val_get(restart_section,"AVE_PV_KIN",r_vals=tmp,error=error)
                      averages%virial%pv_kinetic = RESHAPE(tmp,(/3,3/))
                      CALL section_vals_val_get(restart_section,"AVE_PV_CNSTR",r_vals=tmp,error=error)
                      averages%virial%pv_constraint = RESHAPE(tmp,(/3,3/))
                      CALL section_vals_val_get(restart_section,"AVE_PV_XC",r_vals=tmp,error=error)
                      averages%virial%pv_xc = RESHAPE(tmp,(/3,3/))
                      CALL section_vals_val_get(restart_section,"AVE_PV_FOCK_4C",r_vals=tmp,error=error)
                      averages%virial%pv_fock_4c = RESHAPE(tmp,(/3,3/))
                   END IF
                END IF
                ! Colvars
                IF (nint>0) THEN
                   CALL section_vals_val_get(restart_section,"AVE_COLVARS",r_vals=cvalues,error=error)
                   CALL section_vals_val_get(restart_section,"AVE_MMATRIX",r_vals=Mmatrix,error=error)
                   CPPostcondition(nint==SIZE(cvalues),cp_failure_level,routineP,error,failure)
                   CPPostcondition(nint*nint==SIZE(Mmatrix),cp_failure_level,routineP,error,failure)
                   averages%avecolvar  = cvalues
                   averages%aveMmatrix = Mmatrix
                END IF
             ELSE
                averages%itimes_start  = itimes
             END IF
          END IF
          delta_t = itimes - averages%itimes_start + 1

          ! Perform averages
          SELECT CASE (delta_t)
          CASE (1)
             averages%avecpu        = used_time
             averages%avehugoniot   = hugoniot
             averages%avetemp_baro  = md_ener%temp_baro
             averages%avepot        = md_ener%epot
             averages%avekin        = md_ener%ekin
             averages%avetemp       = md_ener%temp_part
             averages%avekin_qm     = md_ener%ekin_qm
             averages%avetemp_qm    = md_ener%temp_qm
             averages%avevol        = cell%deth
             averages%aveca         = abc(1)
             averages%avecb         = abc(2)
             averages%avecc         = abc(3)
             averages%aveal         = cell_angle(3)
             averages%avebe         = cell_angle(2)
             averages%avega         = cell_angle(1)
             averages%econs         = 0._dp
             ! Virial
             IF (virial%pv_availability) THEN
                averages%avepress = pv_scalar
                averages%avepxx   = pv_xx
                IF (ASSOCIATED(averages%virial)) THEN
                   averages%virial%pv_total      = virial%pv_total
                   averages%virial%pv_virial     = virial%pv_virial
                   averages%virial%pv_kinetic    = virial%pv_kinetic
                   averages%virial%pv_constraint = virial%pv_constraint
                   averages%virial%pv_xc         = virial%pv_xc
                   averages%virial%pv_fock_4c    = virial%pv_fock_4c
                END IF
             END IF
             ! Colvars
             IF (nint>0) THEN
                CALL get_clv_force(force_env, nsize_xyz=nat*3, nsize_int=nint, &
                     cvalues=averages%avecolvar,Mmatrix=averages%aveMmatrix, error=error)
             END IF
          CASE DEFAULT
             CALL get_averages(averages%avecpu,used_time,delta_t,error)
             CALL get_averages(averages%avehugoniot,hugoniot,delta_t,error)
             CALL get_averages(averages%avetemp_baro,md_ener%temp_baro,delta_t,error)
             CALL get_averages(averages%avepot,md_ener%epot,delta_t,error)
             CALL get_averages(averages%avekin,md_ener%ekin,delta_t,error)
             CALL get_averages(averages%avetemp,md_ener%temp_part,delta_t,error)
             CALL get_averages(averages%avekin_qm,md_ener%ekin_qm,delta_t,error)
             CALL get_averages(averages%avetemp_qm,md_ener%temp_qm,delta_t,error)
             CALL get_averages(averages%avevol,cell%deth,delta_t,error)
             CALL get_averages(averages%aveca,abc(1),delta_t,error)
             CALL get_averages(averages%avecb,abc(2),delta_t,error)
             CALL get_averages(averages%avecc,abc(3),delta_t,error)
             CALL get_averages(averages%aveal,cell_angle(3),delta_t,error)
             CALL get_averages(averages%avebe,cell_angle(2),delta_t,error)
             CALL get_averages(averages%avega,cell_angle(1),delta_t,error)
             CALL get_averages(averages%econs,md_ener%delta_cons,delta_t,error)
             ! Virial
             IF (virial%pv_availability) THEN
                CALL get_averages(averages%avepress,pv_scalar,delta_t,error)
                CALL get_averages(averages%avepxx,pv_xx,delta_t,error)
                IF (ASSOCIATED(averages%virial)) THEN
                   CALL get_averages(averages%virial%pv_total,virial%pv_total,delta_t,error)
                   CALL get_averages(averages%virial%pv_virial,virial%pv_virial,delta_t,error)
                   CALL get_averages(averages%virial%pv_kinetic,virial%pv_kinetic,delta_t,error)
                   CALL get_averages(averages%virial%pv_constraint,virial%pv_constraint,delta_t,error)
                   CALL get_averages(averages%virial%pv_xc,virial%pv_xc,delta_t,error)
                   CALL get_averages(averages%virial%pv_fock_4c,virial%pv_fock_4c,delta_t,error)
                END IF
             END IF
             ! Colvars
             IF (nint>0) THEN
                ALLOCATE(cvalues(nint),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                ALLOCATE(Mmatrix(nint*nint),stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                CALL get_clv_force(force_env, nsize_xyz=nat*3, nsize_int=nint, cvalues=cvalues,&
                     Mmatrix=Mmatrix, error=error)
                CALL get_averages(averages%avecolvar, cvalues, delta_t,error)
                CALL get_averages(averages%aveMmatrix, Mmatrix, delta_t,error)
                DEALLOCATE(cvalues,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                DEALLOCATE(Mmatrix,stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END SELECT
       END IF

       ! Possibly print averages
       output_unit = cp_print_key_unit_nr(logger,averages%averages_section,"PRINT_AVERAGES",&
            extension=".avg", file_position=my_pos, file_action=my_act, error=error)
       IF (output_unit>0) THEN
          WRITE(output_unit,FMT='(A15,1X,"=",1X,G15.9,"   NSTEP #",I15)')&
               "AVECPU",averages%avecpu,itimes,&
               "AVEHUGONIOT",averages%avehugoniot,itimes,&
               "AVETEMP_BARO",averages%avetemp_baro,itimes,&
               "AVEPOT",averages%avepot,itimes,&
               "AVEKIN",averages%avekin,itimes,&
               "AVETEMP",averages%avetemp,itimes,&
               "AVEKIN_QM",averages%avekin_qm,itimes,&
               "AVETEMP_QM",averages%avetemp_qm,itimes,&
               "AVEVOL",averages%avevol,itimes,&
               "AVECELL_A",averages%aveca,itimes,&
               "AVECELL_B",averages%avecb,itimes,&
               "AVECELL_C",averages%avecc,itimes,&
               "AVEALPHA",averages%aveal,itimes,&
               "AVEBETA",averages%avebe,itimes,&
               "AVEGAMMA",averages%avega,itimes,&
               "AVE_ECONS",averages%econs,itimes
          ! Print the virial
          IF (virial%pv_availability) THEN
             WRITE(output_unit,FMT='(A15,1X,"=",1X,G15.9,"   NSTEP #",I15)')&
                  "AVE_PRESS",averages%avepress,itimes,&
                  "AVE_PXX",averages%avepxx,itimes
             IF (ASSOCIATED(averages%virial)) THEN
                WRITE(output_unit,FMT='(A15,1X,"=",1X,G15.9,"   NSTEP #",I15)')&
                     "AVE_PV_TOT",averages%virial%pv_total,itimes,&
                     "AVE_PV_VIR",averages%virial%pv_virial,itimes,&
                     "AVE_PV_KIN",averages%virial%pv_kinetic,itimes,&
                     "AVE_PV_CNSTR",averages%virial%pv_constraint,itimes,&
                     "AVE_PV_XC",averages%virial%pv_xc,itimes,&
                     "AVE_PV_FOCK_4C",averages%virial%pv_fock_4c,itimes
             END IF
          END IF
          DO i = 1, nint
             ctmp = cp_to_string(i)
             WRITE(output_unit,FMT='(A15,1X,"=",1X,G15.9,"   NSTEP #",I15)')&
                  TRIM("AVE_CV-"//ADJUSTL(ctmp)),averages%avecolvar(i),itimes
          END DO
          WRITE(output_unit,FMT='(/)')
       END IF
       CALL cp_print_key_finished_output(output_unit,logger,averages%averages_section,&
            "PRINT_AVERAGES", error=error)
    END IF
    CALL timestop(handle)
  END SUBROUTINE compute_averages

! *****************************************************************************
!> \brief computes the averages - low level for REAL
!> \param avg ...
!> \param add ...
!> \param delta_t ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 03.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE get_averages_rs(avg, add, delta_t, error)
    REAL(KIND=dp), INTENT(INOUT)             :: avg
    REAL(KIND=dp), INTENT(IN)                :: add
    INTEGER, INTENT(IN)                      :: delta_t
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_averages_rs', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    avg = (avg*REAL(delta_t-1,dp) + add)/REAL(delta_t,dp)
  END SUBROUTINE get_averages_rs

! *****************************************************************************
!> \brief computes the averages - low level for REAL vector
!> \param avg ...
!> \param add ...
!> \param delta_t ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 10.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE get_averages_rv(avg, add, delta_t, error)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: avg
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: add
    INTEGER, INTENT(IN)                      :: delta_t
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_averages_rv', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    check = SIZE(avg)==SIZE(add)
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    DO i = 1, SIZE(avg)
       avg(i) = (avg(i)*REAL(delta_t-1,dp) + add(i))/REAL(delta_t,dp)
    END DO
  END SUBROUTINE get_averages_rv

! *****************************************************************************
!> \brief computes the averages - low level for REAL matrix
!> \param avg ...
!> \param add ...
!> \param delta_t ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 10.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE get_averages_rm(avg, add, delta_t, error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: avg
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: add
    INTEGER, INTENT(IN)                      :: delta_t
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_averages_rm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: check, failure

    failure = .FALSE.
    check = SIZE(avg,1)==SIZE(add,1)
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    check = SIZE(avg,2)==SIZE(add,2)
    CPPostcondition(check,cp_failure_level,routineP,error,failure)
    DO i = 1, SIZE(avg,2)
       DO j = 1, SIZE(avg,1)
          avg(j,i) = (avg(j,i)*REAL(delta_t-1,dp) + add(j,i))/REAL(delta_t,dp)
       END DO
    END DO
  END SUBROUTINE get_averages_rm

END MODULE averages_types
