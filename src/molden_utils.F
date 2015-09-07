!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Functions handling the MOLDEN format. Split from mode_selective.
!> \author Teodoro Laino, 03.2009
! *****************************************************************************
MODULE molden_utils
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molden_utils'
  LOGICAL, PARAMETER                   :: debug_this_module=.FALSE.

  PUBLIC :: molden_out

CONTAINS
! *****************************************************************************
!> \brief writes the output for vibrational analysis in MOLDEN format
!> \param input ...
!> \param particles ...
!> \param freq ...
!> \param eigen_vec ...
!> \param intensities ...
!> \param calc_intens ...
!> \param dump_only_positive ...
!> \param logger ...
!> \param error ...
!> \author Florian Schiffmann 11.2007
! *****************************************************************************
  SUBROUTINE molden_out(input,particles,freq,eigen_vec,intensities,calc_intens,&
             dump_only_positive,logger,error)

    TYPE(section_vals_type), POINTER         :: input
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(:)              :: freq
    REAL(KIND=dp), DIMENSION(:, :)           :: eigen_vec
    REAL(KIND=dp), DIMENSION(:), POINTER     :: intensities
    LOGICAL                                  :: calc_intens, &
                                                dump_only_positive
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'molden_out', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: element_symbol
    INTEGER                                  :: i, iw, j, k, l
    LOGICAL                                  :: failure

    iw=cp_print_key_unit_nr(logger,input,"VIBRATIONAL_ANALYSIS%PRINT%MOLDEN_VIB",&
            extension=".mol",file_status='REPLACE',error=error)

    IF(iw.GT.0)THEN
       CPPostcondition(MOD(SIZE(eigen_vec,1),3)==0,cp_failure_level,routineP,error,failure)
       CPPostcondition(SIZE(particles)==SIZE(eigen_vec,1)/3,cp_failure_level,routineP,error,failure)
       CPPostcondition(SIZE(freq,1)==SIZE(eigen_vec,2),cp_failure_level,routineP,error,failure)
       WRITE(iw,'(T2,A)')"[Molden Format]"
       WRITE(iw,'(T2,A)')"[FREQ]"
       DO i=1,SIZE(freq,1)
          IF((.NOT.dump_only_positive).OR.(freq(i)>=0._dp))WRITE(iw,'(T5,F12.6)') freq(i)
       END DO
       WRITE(iw,'(T2,A)')"[FR-COORD]"
       DO i=1,SIZE(particles)
          CALL get_atomic_kind(atomic_kind=particles(i)%atomic_kind,&
               element_symbol=element_symbol)
               WRITE(iw,'(T2,A2,3X,3(F12.6,3X))')&
                    element_symbol, particles((i))%r(:)
       END DO
       WRITE(iw,'(T2,A)')"[FR-NORM-COORD]"
       l=0
       DO i=1,SIZE(eigen_vec,2)
          IF ((.NOT.dump_only_positive).OR.(freq(i)>=0._dp)) THEN
             l=l+1
             WRITE(iw,'(T2,A,1X,I6)')"vibration",l
             DO j=1,SIZE(eigen_vec,1)/3
                k=(j-1)*3
                WRITE(iw,'(T2,3(F12.6,3X))')eigen_vec(k+1,i),eigen_vec(k+2,i),eigen_vec(k+3,i)
             END DO
          END IF
       END DO
       IF(calc_intens)THEN
          WRITE(iw,'(T2,A)')"[INT]"
          DO i=1,SIZE(intensities)
             IF((.NOT.dump_only_positive).OR.(freq(i)>=0._dp))WRITE(iw,'(3X,F18.6)')intensities(i)
          END DO
       END IF
    END IF
    CALL cp_print_key_finished_output(iw,logger,input,"VIBRATIONAL_ANALYSIS%PRINT%MOLDEN_VIB",error=error)
  END SUBROUTINE molden_out

END MODULE molden_utils
