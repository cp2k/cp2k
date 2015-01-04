!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief types used in the flexible partitioning scheme
!> \par History
!>      04.2006 [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE fp_types
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: fp_type
  PUBLIC :: fp_env_create,fp_env_release,fp_env_retain,fp_env_read,fp_env_write

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fp_types'

! *****************************************************************************
  TYPE fp_type

    ! input related objects
    LOGICAL                        :: use_fp
    INTEGER                        :: ref_count

    INTEGER                        :: central_atom
    INTEGER, DIMENSION(:), POINTER :: inner_atoms, outer_atoms
    REAL(KIND=dp)                  :: inner_radius, outer_radius
    REAL(KIND=dp)                  :: strength, smooth_width
    LOGICAL                        :: bias
    REAL(KIND=dp)                  :: temperature
    TYPE(section_vals_type), POINTER   :: print_section

    ! computed during runs
    INTEGER                        :: i1,i2,o1,o2
    REAL(KIND=dp)                  :: ri1,ri2,ro1,ro2
    REAL(KIND=dp)                  :: weight, comb_weight, bias_weight
    REAL(KIND=dp)                  :: energy, bias_energy, restraint_energy
  END TYPE fp_type

CONTAINS

! *****************************************************************************
!> \brief create retain release the flexible partitioning environment
!> \param fp_env ...
!> \param error ...
!> \par History
!>      04.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE fp_env_create(fp_env,error)
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fp_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
      failure=.FALSE.

      CPPrecondition(.NOT.ASSOCIATED(fp_env),cp_failure_level,routineP,error,failure)
      ALLOCATE(fp_env,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      fp_env%ref_count=1

      fp_env%use_fp=.FALSE.
      NULLIFY(fp_env%inner_atoms)
      NULLIFY(fp_env%outer_atoms)
      NULLIFY(fp_env%print_section)
    CALL timestop(handle)

  END SUBROUTINE fp_env_create

! *****************************************************************************
!> \brief ...
!> \param fp_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fp_env_release(fp_env,error)
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fp_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

     failure=.FALSE.
     IF (ASSOCIATED(fp_env)) THEN
        CPPrecondition(fp_env%ref_count>0,cp_failure_level,routineP,error,failure)
        fp_env%ref_count=fp_env%ref_count-1
        IF (fp_env%ref_count==0) THEN
           IF (ASSOCIATED(fp_env%inner_atoms)) DEALLOCATE(fp_env%inner_atoms)
           IF (ASSOCIATED(fp_env%outer_atoms)) DEALLOCATE(fp_env%outer_atoms)
           IF (ASSOCIATED(fp_env%print_section)) CALL section_vals_release(fp_env%print_section,error)
           fp_env%use_fp=.FALSE.
           DEALLOCATE(fp_env)
        ENDIF
     ENDIF
    CALL timestop(handle)

  END SUBROUTINE fp_env_release

! *****************************************************************************
!> \brief ...
!> \param fp_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fp_env_retain(fp_env,error)
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fp_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

      failure=.FALSE.
      CPPrecondition(ASSOCIATED(fp_env),cp_failure_level,routineP,error,failure)
      fp_env%ref_count=fp_env%ref_count+1

  END SUBROUTINE fp_env_retain

! *****************************************************************************
!> \brief reads the corresponding input section and stores it in the fp_env
!> \param fp_env ...
!> \param fp_section ...
!> \param error ...
!> \par History
!>      04.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE fp_env_read(fp_env,fp_section,error)
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(section_vals_type), POINTER         :: fp_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fp_env_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    INTEGER, DIMENSION(:), POINTER           :: tmplist
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(fp_env),cp_failure_level,routineP,error,failure)
    CALL section_vals_get(fp_section,explicit=fp_env%use_fp,error=error)
    IF (fp_env%use_fp) THEN
       CALL section_vals_val_get(fp_section,"CENTRAL_ATOM",i_val=fp_env%central_atom,error=error)

       CALL section_vals_val_get(fp_section,"INNER_ATOMS",i_vals=tmplist,error=error)
       ALLOCATE(fp_env%inner_atoms(SIZE(tmplist,1)))
       fp_env%inner_atoms=tmplist

       CALL section_vals_val_get(fp_section,"OUTER_ATOMS",i_vals=tmplist,error=error)
       ALLOCATE(fp_env%outer_atoms(SIZE(tmplist,1)))
       fp_env%outer_atoms=tmplist

       CALL section_vals_val_get(fp_section,"INNER_RADIUS",r_val=fp_env%inner_radius,error=error)
       CALL section_vals_val_get(fp_section,"OUTER_RADIUS",r_val=fp_env%outer_radius,error=error)
       CALL section_vals_val_get(fp_section,"STRENGTH",r_val=fp_env%strength,error=error)
       CALL section_vals_val_get(fp_section,"SMOOTH_WIDTH",r_val=fp_env%smooth_width,error=error)
       CALL section_vals_val_get(fp_section,"BIAS",l_val=fp_env%bias,error=error)
       CALL section_vals_val_get(fp_section,"TEMPERATURE",r_val=fp_env%temperature,error=error)

       fp_env%print_section=>section_vals_get_subs_vals(fp_section,"WEIGHTS",error=error)
       CALL section_vals_retain(fp_env%print_section,error=error)
    ENDIF
    CALL timestop(handle)

  END SUBROUTINE fp_env_read

! *****************************************************************************
!> \brief writes information concerning the fp_env to the output
!> \param fp_env ...
!> \param fp_section ...
!> \param error ...
!> \par History
!>      04.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE fp_env_write(fp_env,fp_section,error)
    TYPE(fp_type), POINTER                   :: fp_env
    TYPE(section_vals_type), POINTER         :: fp_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fp_env_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: PRINT

    CALL timeset(routineN,handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(fp_env),cp_failure_level,routineP,error,failure)

    IF (fp_env%use_fp) THEN
       PRINT=>section_vals_get_subs_vals(fp_section,"CONTROL",error=error)
       output_unit=cp_print_key_unit_nr(logger,PRINT,"",extension=".Log",error=error)
       IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(T2,A,T79,A)")&
             "FP| Flexible partitioning is ","ON"
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
             "FP| Central atom ", fp_env%central_atom
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
             "FP| number of inner atoms", SIZE(fp_env%inner_atoms,1)
          WRITE (UNIT=output_unit,FMT="(1(T2,8I8))") fp_env%inner_atoms
          WRITE (UNIT=output_unit,FMT="(T2,A,T71,I10)")&
             "FP| number of outer atoms", SIZE(fp_env%outer_atoms,1)
          WRITE (UNIT=output_unit,FMT="(1(T2,8I8))") fp_env%outer_atoms
          WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
             "FP| inner radius [a.u.] ", fp_env%inner_radius
          WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
             "FP| outer radius [a.u.] ", fp_env%outer_radius
          WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
             "FP| reflecting restraint strength ", fp_env%strength
          IF (fp_env%bias) THEN
             WRITE (UNIT=output_unit,FMT="(T2,A,T79,A)")&
                     "FP| Flexible partitioning bias is " ,"ON"
             WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
                     "FP| bias temperature [kT a.u.]", fp_env%temperature
             WRITE (UNIT=output_unit,FMT="(T2,A,T61,F20.10)")&
                     "FP| smooth width [a.u.] ", fp_env%smooth_width
          ELSE
             WRITE (UNIT=output_unit,FMT="(T2,A,T78,A)")&
                     "FP| Flexible partitioning bias is" ,"OFF"
          ENDIF
       ENDIF
       CALL cp_print_key_finished_output(output_unit,logger,PRINT,"",error=error)
    ENDIF
    CALL timestop(handle)

  END SUBROUTINE fp_env_write

END MODULE fp_types
