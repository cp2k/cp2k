!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Barostat utils
!> \author teo [tlaino] - University of Zurich - 02.2008
! *****************************************************************************
MODULE barostat_utils
  USE barostat_types,                  ONLY: barostat_type
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE extended_system_types,           ONLY: npt_info_type
  USE input_constants,                 ONLY: npe_f_ensemble,&
                                             npe_i_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             npt_f_ensemble,&
                                             npt_i_ensemble
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE physcon,                         ONLY: angstrom,&
                                             femtoseconds,&
                                             kelvin
  USE simpar_types,                    ONLY: simpar_type
#include "../../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  get_baro_energies, print_barostat_status

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'barostat_utils'

CONTAINS
! *****************************************************************************
!> \brief Calculates kinetic energy and potential of barostat
!> \param cell ...
!> \param simpar ...
!> \param npt ...
!> \param baro_kin ...
!> \param baro_pot ...
!> \param error ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE get_baro_energies ( cell, simpar, npt, baro_kin, baro_pot, error )

    TYPE(cell_type), POINTER                 :: cell
    TYPE(simpar_type), INTENT(IN)            :: simpar
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(IN)                             :: npt
    REAL(KIND=dp), INTENT(OUT)               :: baro_kin, baro_pot
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_baro_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    REAL(dp)                                 :: iv0, v0, v_shock

    IF ( simpar%ensemble==npt_i_ensemble .OR. simpar%ensemble==npe_i_ensemble) THEN
       baro_pot = simpar%p_ext * cell%deth
       baro_kin = 0.5_dp * npt ( 1, 1 )%v**2 * npt ( 1, 1 )%mass
    ELSE IF (simpar%ensemble==npt_f_ensemble .OR. simpar%ensemble==npe_f_ensemble) THEN
       baro_pot =  simpar%p_ext * cell%deth
       baro_kin = 0.0_dp
       DO i = 1, 3
          DO j = 1, 3
             baro_kin = baro_kin + 0.5_dp * npt(i,j)%v**2  * npt ( i, j )%mass
          END DO
       END DO
    ELSEIF ( simpar%ensemble == nph_uniaxial_ensemble .OR. simpar%ensemble == nph_uniaxial_damped_ensemble ) THEN
       v0 = simpar%v0
       iv0 = 1._dp / v0
       v_shock = simpar%v_shock

       ! Valid only for orthorhombic cell
       baro_pot = -0.5_dp * v_shock * v_shock * ( 1._dp - cell%deth * iv0 )**2 - simpar%p0 * ( v0 - cell%deth )
       ! Valid only for orthorhombic cell
       baro_kin = 0.5_dp * npt ( 1, 1 )%v * npt ( 1, 1 )%v * npt ( 1, 1 )%mass
    END IF

  END SUBROUTINE get_baro_energies

! *****************************************************************************
!> \brief Prints status of barostat during an MD run
!> \param barostat ...
!> \param simpar ...
!> \param my_pos ...
!> \param my_act ...
!> \param cell ...
!> \param itimes ...
!> \param time ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 02.2008 - University of Zurich
! *****************************************************************************
  SUBROUTINE print_barostat_status ( barostat, simpar, my_pos, my_act, cell, itimes, time, error )
    TYPE(barostat_type), POINTER             :: barostat
    TYPE(simpar_type), INTENT(IN)            :: simpar
    CHARACTER(LEN=default_string_length)     :: my_pos, my_act
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: itimes
    REAL(KIND=dp), INTENT(IN)                :: time
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'print_barostat_status', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: baro, nfree
    LOGICAL                                  :: failure, new_file
    REAL(KIND=dp)                            :: baro_kin, baro_pot, temp
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF (ASSOCIATED(barostat)) THEN
       baro = cp_print_key_unit_nr(logger,barostat%section,"PRINT%ENERGY",&
            extension=".bener",file_position=my_pos, file_action=my_act, is_new_file=new_file,&
            error=error)
       CALL  get_baro_energies ( cell, simpar, barostat%npt, baro_kin, baro_pot, error )
       nfree = SIZE(barostat%npt,1) * SIZE(barostat%npt, 2)
       temp = 2.0_dp * baro_kin / REAL(nfree,dp)*kelvin
       IF(baro > 0) THEN
          IF (new_file) THEN
             WRITE (baro,'("#",3X,A,10X,A,8X,3(5X,A,5X),3X,A)')"Step Nr.","Time[fs]","Kin.[a.u.]",&
                    "Temp.[K]","Pot.[a.u.]","Vol[Ang.^3]"
          END IF
          WRITE(UNIT=baro,FMT="(I10, F20.3,4F20.10)") itimes,time*femtoseconds,&
               baro_kin, temp, baro_pot, cell%deth*angstrom*angstrom*angstrom
          CALL m_flush(baro)
       END IF
       CALL cp_print_key_finished_output(baro,logger,barostat%section,"PRINT%ENERGY", error=error)
    END IF

  END SUBROUTINE print_barostat_status

END MODULE barostat_utils
