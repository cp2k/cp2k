!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/md_conserved_quantities [1.0] *
!!
!!   NAME
!!     conserved_quantities
!!
!!   FUNCTION
!!     computes the conserved quantities for a given md ensemble
!!     and also kinetic energies, thermo/barostat stuff
!!
!!   AUTHOR
!!     gtb, 05.02.2003
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE md_conserved_quantities

  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE dynamical_coeff_types,           ONLY: dyn_coeff_set_type,&
                                             dyn_coeff_type
  USE extended_system_types,           ONLY: npt_info_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: &
       dyn_coeff_mass_global, isokin_ensemble, langevin_ensemble, &
       nph_ensemble, nph_uniaxial_damped_ensemble, nph_uniaxial_ensemble, &
       npt_f_ensemble, npt_i_ensemble, nve_ensemble, nvt_ensemble, &
       reftraj_ensemble
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: zero
  USE md,                              ONLY: simulation_parameters_type
  USE md_environment_types,            ONLY: get_md_env,&
                                             md_environment_type,&
                                             set_md_env
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE shell_potential_types,           ONLY: shell_kind_type
  USE termination,                     ONLY: stop_program
  USE thermostat_types,                ONLY: thermostat_type
  USE thermostat_utils,                ONLY: get_thermostat_energies
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  TYPE md_ener_type
     REAL(KIND=dp) :: epot
     REAL(KIND=dp) :: ekin, ekin_qm
     REAL(KIND=dp) :: temp_part, temp_qm
     REAL(KIND=dp) :: temp_baro
     REAL(KIND=dp) :: ekin_coefs
     REAL(KIND=dp) :: temp_coefs
     REAL(KIND=dp) :: ekin_shell, temp_shell
     REAL(KIND=dp), DIMENSION(:), POINTER :: ekin_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: temp_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: ekin_shell_kind
     REAL(KIND=dp), DIMENSION(:), POINTER :: temp_shell_kind
  END TYPE md_ener_type


  PUBLIC :: compute_conserved_quantity
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'md_conserved_quantities'

  !******************************************************************************
CONTAINS

  !!*****************************************************************************
  !!******
  !!
  !!   NAME
  !!     compute_conserved_quantity
  !!
  !!   FUNCTION
  !!     calculates conserved quantity.
  !!
  !!   INPUT ARGUMENTS
  !!    md_env is the md_environment
  !!    epot is the total potential energy
  !!   OUTPUT ARGUMENTS
  !!    cons is the conserved quantity
  !!   OUTPUT OPTIONAL ARGUMENTS
  !!    cons_rel : relative cons. quantity (to the first md step)
  !!    ekin : kinetic energy of particles
  !!    temp : temperature
  !!    temp_qm : temperature of the QM system in a QM/MM calculation
  !!    ekin_coefs : ke of the coefficients (KG code)
  !!    temp_coefs : temperature of the coefficients (KG code)
  !!
  !!   AUTHOR
  !!    gloria
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** *************************************************************************
  SUBROUTINE compute_conserved_quantity (md_env, epot, cons, cons_rel, ekin,&
       temp_part, ekin_qm, temp_qm, temp_baro, ekin_coefs, temp_coefs,  &
       ekin_shell, temp_shell, temp_kind, temp_shell_kind, qmmm, &
       shell_adiabatic, shell_present, temperature_per_kind, error )
    TYPE(md_environment_type), POINTER       :: md_env
    REAL(KIND=dp), INTENT(IN)                :: epot
    REAL(KIND=dp), INTENT(OUT)               :: cons
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: cons_rel, ekin, temp_part, &
      ekin_qm, temp_qm, temp_baro, ekin_coefs, temp_coefs, ekin_shell, &
      temp_shell
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: temp_kind, temp_shell_kind
    LOGICAL, INTENT(OUT), OPTIONAL           :: qmmm
    LOGICAL, INTENT(INOUT), OPTIONAL         :: shell_adiabatic
    LOGICAL, INTENT(OUT), OPTIONAL           :: shell_present
    LOGICAL, INTENT(IN), OPTIONAL            :: temperature_per_kind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_conserved_quantity', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ikind, int_grp, istat, &
                                                nfree_qm, nfree_shell, nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: nfree_kind, nfree_shell_kind
    INTEGER, POINTER                         :: itimes
    LOGICAL                                  :: failure, my_shell_adiabatic, &
                                                my_shell_present, &
                                                my_temperature_per_kind, &
                                                myqmmm
    REAL(KIND=dp), POINTER                   :: constant
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(md_ener_type)                       :: md_ener
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

    failure = .FALSE.  

    md_ener%temp_part   = 0.0_dp
    md_ener%temp_baro   = 0.0_dp
    md_ener%ekin_coefs  = 0.0_dp
    md_ener%temp_coefs  = 0.0_dp
    md_ener%temp_qm     = 0.0_dp
    md_ener%ekin_shell  = 0.0_dp
    md_ener%temp_shell  = 0.0_dp
    md_ener%epot=epot
    myqmmm = .FALSE.
    my_shell_adiabatic =.FALSE.
    IF(PRESENT(shell_adiabatic)) my_shell_adiabatic = shell_adiabatic
    my_shell_present = .FALSE.
    my_temperature_per_kind = .FALSE.
    IF(PRESENT(temperature_per_kind)) my_temperature_per_kind = temperature_per_kind
    
    NULLIFY(md_ener%ekin_kind,md_ener%ekin_shell_kind)
    NULLIFY(md_ener%temp_kind,md_ener%temp_shell_kind)
    IF(my_temperature_per_kind) THEN
      CPPostcondition(ASSOCIATED(temp_kind),cp_failure_level,routineP,error,failure)
      nkind = SIZE(temp_kind,1)
      temp_kind = 0.0_dp
      md_ener%temp_kind => temp_kind
      ALLOCATE(md_ener%ekin_kind(nkind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      md_ener%ekin_kind = 0.0_dp
      ALLOCATE(nfree_kind(nkind),STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      nfree_kind = 0
      IF(my_shell_adiabatic) THEN 
        CPPostcondition(ASSOCIATED(temp_shell_kind),cp_failure_level,routineP,error,failure)
        temp_shell_kind = 0.0_dp
        md_ener%temp_shell_kind => temp_shell_kind
        ALLOCATE(md_ener%ekin_shell_kind(nkind),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        md_ener%ekin_shell_kind = 0.0_dp
        ALLOCATE(nfree_shell_kind(nkind),STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        nfree_shell_kind = 0
      ENDIF  
    END IF
    
    CALL get_md_env (md_env=md_env,&
                     constant=constant,&
                     itimes=itimes,&
                     simpar=simpar,&
                     para_env = para_env, &
                     error=error)

    int_grp =  para_env % group

    IF(my_temperature_per_kind) THEN
      CALL get_part_ke ( md_env, md_ener%ekin, md_ener%ekin_qm, md_ener%ekin_shell,&
           nfree_qm, nfree_shell, simpar%nfree, int_grp , my_shell_present,&
           temperature_per_kind=my_temperature_per_kind, ekin_kind=md_ener%ekin_kind,&
           ekin_shell_kind=md_ener%ekin_shell_kind, nfree_kind=nfree_kind, &
           nfree_shell_kind=nfree_shell_kind, error=error)
    ELSE        
      CALL get_part_ke ( md_env, md_ener%ekin, md_ener%ekin_qm, md_ener%ekin_shell,&
           nfree_qm, nfree_shell, simpar%nfree, int_grp , my_shell_present, &
           temperature_per_kind=my_temperature_per_kind, error=error)
    END IF
           
    IF (simpar%nfree /= 0) &
         md_ener % temp_part = 2.0_dp * md_ener % ekin / REAL ( simpar%nfree,KIND=dp)

    IF ( nfree_qm > 0) THEN
       myqmmm = .TRUE.
       md_ener % temp_qm   = 2.0_dp * md_ener % ekin_qm / REAL ( nfree_qm,KIND=dp)
    END IF

    IF ( nfree_shell > 0) THEN
       my_shell_adiabatic = .TRUE.
       md_ener % temp_shell   = 2.0_dp * md_ener % ekin_shell / REAL ( nfree_shell,KIND=dp)
    ELSEIF(my_shell_present) THEN
        md_ener % temp_shell  = 0.0_dp
    END IF

    IF(my_temperature_per_kind) THEN
       DO ikind = 1,nkind
         temp_kind(ikind) = 2.0_dp * &
                md_ener%ekin_kind(ikind)/REAL(nfree_kind(ikind),KIND=dp)
       END DO
       IF(my_shell_adiabatic) THEN
         DO ikind = 1,nkind
            temp_shell_kind(ikind) = 2.0_dp * &
                    md_ener%ekin_shell_kind(ikind)/REAL(nfree_shell_kind(ikind),KIND=dp)
         END DO
       END IF
    END IF

    SELECT CASE ( simpar % ensemble  )
    CASE DEFAULT
       ! hence, should be implemented
       CALL stop_program( RoutineN, 'Unknown Ensemble')
    CASE ( isokin_ensemble )
       cons = md_ener % ekin
    CASE ( reftraj_ensemble ) ! no constant of motion available
       cons = md_ener % epot
    CASE ( nve_ensemble)
       CALL get_econs_nve (md_env, md_ener, cons, int_grp,error=error)
    CASE ( nvt_ensemble)
       CALL get_econs_nvt (md_env, md_ener, cons, int_grp,error=error)
    CASE (  npt_i_ensemble,  npt_f_ensemble )
       CALL get_econs_npt (md_env, md_ener, cons, int_grp, error=error)
    CASE (  nph_uniaxial_ensemble )
       CALL get_econs_nph_uniaxial (md_env, md_ener, cons, int_grp, error=error)
    CASE ( nph_uniaxial_damped_ensemble )
       CALL get_econs_nph_uniaxial (md_env, md_ener, cons, int_grp, error=error)
    CASE ( langevin_ensemble )
       cons = md_ener%ekin + md_ener%epot
    END SELECT

    IF (itimes==1) THEN
       constant = cons
       CALL set_md_env (md_env=md_env, constant=constant,error=error)
    END IF

    IF ((PRESENT(cons_rel)).AND.(itimes /= 0 )) THEN
       cons_rel = ABS((cons-constant)/constant)
    END IF
    IF (PRESENT(ekin))            ekin       = md_ener%ekin
    IF (PRESENT(ekin_qm))         ekin_qm    = md_ener%ekin_qm
    IF (PRESENT(temp_part))       temp_part  = md_ener%temp_part
    IF (PRESENT(temp_qm))         temp_qm    = md_ener%temp_qm
    IF (PRESENT(temp_baro))       temp_baro  = md_ener%temp_baro
    IF (PRESENT(ekin_coefs))      ekin_coefs = md_ener%ekin_coefs
    IF (PRESENT(temp_coefs))      temp_coefs = md_ener%temp_coefs
    IF (PRESENT(qmmm))            qmmm       = myqmmm
    IF (PRESENT(ekin_shell))      ekin_shell = md_ener%ekin_shell
    IF (PRESENT(temp_shell))      temp_shell = md_ener%temp_shell
    IF (PRESENT(shell_adiabatic)) shell_adiabatic = my_shell_adiabatic
    IF (PRESENT(shell_present))   shell_present   = my_shell_present

    IF(my_temperature_per_kind) THEN
      NULLIFY(md_ener%temp_kind)      
      DEALLOCATE(md_ener%ekin_kind,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      DEALLOCATE(md_ener%ekin_shell_kind,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      IF(my_shell_adiabatic) THEN
        NULLIFY(md_ener%temp_shell_kind)      
        DEALLOCATE(nfree_kind,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        DEALLOCATE(nfree_shell_kind,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
    END IF
    
  END SUBROUTINE compute_conserved_quantity

  !******************************************************************************
  !!******
  !!
  !!   NAME
  !!     get_econs_nve
  !!
  !!   FUNCTION
  !!     calculates conserved quantity for nvt ensemble
  !!
  !!   AUTHOR
  !!    gloria
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** *************************************************************************
  SUBROUTINE get_econs_nve ( md_env, md_ener, cons, int_grp, error )
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)               :: cons
    INTEGER, INTENT(IN)                      :: int_grp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_econs_nve', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: thermostat_kin, thermostat_pot
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(thermostat_type), POINTER           :: thermostat

    NULLIFY(dyn_coeff_set, force_env, thermostat)

    CALL get_md_env (md_env, force_env=force_env, thermostat_coeff=thermostat, error=error)
    cons = md_ener % ekin + md_ener % epot + md_ener % ekin_shell

    ! Correct for KG coeff
    CALL force_env_get (force_env = force_env, dyn_coeff_set=dyn_coeff_set,error=error)

    IF (ASSOCIATED (dyn_coeff_set)) THEN
       CALL get_coef_ke (dyn_coeff_set, md_ener%ekin_coefs, md_ener%temp_coefs, int_grp, error)
       cons = cons + md_ener%ekin_coefs
       CALL get_thermostat_energies ( thermostat, thermostat_pot, thermostat_kin, int_grp, error )
       cons = cons + thermostat_kin + thermostat_pot
    END IF
  END SUBROUTINE get_econs_nve

  !******************************************************************************
  !!******
  !!
  !!   NAME
  !!     get_econs_nvt
  !!
  !!   FUNCTION
  !!     calculates conserved quantity for nvt ensemble
  !!
  !!   AUTHOR
  !!    gloria
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** *************************************************************************
  SUBROUTINE get_econs_nvt ( md_env, md_ener, cons, int_grp, error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)               :: cons
    INTEGER, INTENT(IN)                      :: int_grp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_econs_nvt', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: thermo_coef_kin, thermo_coef_pot, thermo_part_kin, &
      thermo_part_pot, thermo_shell_kin, thermo_shell_pot
    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(thermostat_type), POINTER           :: thermostat_coeff, &
                                                thermostat_part, &
                                                thermostat_shell

    NULLIFY(dyn_coeff_set, force_env, thermostat_part, thermostat_coeff, thermostat_shell)
    CALL get_md_env (md_env, force_env=force_env, thermostat_part=thermostat_part,&
         thermostat_coeff=thermostat_coeff, thermostat_shell=thermostat_shell, error=error )
    CALL get_thermostat_energies ( thermostat_part, thermo_part_pot, thermo_part_kin, int_grp, error )
    cons = md_ener% ekin + md_ener%epot +  md_ener % ekin_shell + thermo_part_kin + thermo_part_pot
    CALL force_env_get (force_env=force_env, dyn_coeff_set=dyn_coeff_set,error=error)
    IF (ASSOCIATED (dyn_coeff_set)) THEN
       CALL get_coef_ke (dyn_coeff_set, md_ener%ekin_coefs, md_ener%temp_coefs, int_grp, error)
       cons = cons + md_ener%ekin_coefs
       CALL get_thermostat_energies ( thermostat_coeff, thermo_coef_pot, thermo_coef_kin, int_grp, error )
       cons = cons + thermo_coef_kin + thermo_coef_pot
    END IF

    CALL get_thermostat_energies ( thermostat_shell, thermo_shell_pot, thermo_shell_kin, int_grp, error )
    cons = cons + thermo_shell_kin + thermo_shell_pot

  END SUBROUTINE get_econs_nvt

  !******************************************************************************
  !!******
  !!
  !!   NAME
  !!     get_econs_npt
  !!
  !!   FUNCTION
  !!     calculates conserved quantity for npt ensemble
  !!
  !!   AUTHOR
  !!     gloria
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*****************************************************************************
  SUBROUTINE get_econs_npt ( md_env, md_ener, cons, int_grp, error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)               :: cons
    INTEGER, INTENT(IN)                      :: int_grp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_econs_npt', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp) :: baro_kin, baro_pot, nfree, thermo_baro_kin, &
      thermo_baro_pot, thermo_part_kin, thermo_part_pot, thermo_shell_kin, &
      thermo_shell_pot
    TYPE(cell_type), POINTER                 :: box
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar
    TYPE(thermostat_type), POINTER           :: thermostat_baro, &
                                                thermostat_part, &
                                                thermostat_shell

    NULLIFY(thermostat_baro, thermostat_part, thermostat_shell, npt)
    CALL get_md_env ( md_env, thermostat_part=thermostat_part, thermostat_baro=thermostat_baro, &
         simpar = simpar, npt = npt, cell = box,  thermostat_shell=thermostat_shell, error=error)
    CALL get_thermostat_energies(thermostat_part, thermo_part_pot, thermo_part_kin, int_grp, error=error )
    CALL get_thermostat_energies(thermostat_baro, thermo_baro_pot, thermo_baro_kin, int_grp, error)
    CALL get_baro_energies ( box, simpar, npt, baro_kin, baro_pot, int_grp, error=error )
    nfree = SIZE ( npt, 1 ) * SIZE ( npt, 2 )
    md_ener % temp_baro = 2.0_dp * baro_kin / nfree

    cons = md_ener%ekin + md_ener%epot + md_ener%ekin_shell &
         + thermo_part_kin + thermo_part_pot  &
         + thermo_baro_kin + thermo_baro_pot  &
         + baro_kin + baro_pot

    CALL get_thermostat_energies ( thermostat_shell, thermo_shell_pot, thermo_shell_kin, int_grp, error )
    cons = cons + thermo_shell_kin + thermo_shell_pot

  END SUBROUTINE get_econs_npt

  !******************************************************************************
  !!******
  !!
  !!   NAME
  !!     get_econs_nph_uniaxial
  !!
  !!   FUNCTION
  !!     calculates conserved quantity for nph_uniaxial
  !!
  !!   AUTHOR
  !!     cjm
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*****************************************************************************
  SUBROUTINE get_econs_nph_uniaxial ( md_env, md_ener, cons, int_grp, error)
    TYPE(md_environment_type), POINTER       :: md_env
    TYPE(md_ener_type), INTENT(inout)        :: md_ener
    REAL(KIND=dp), INTENT(OUT)               :: cons
    INTEGER, INTENT(IN)                      :: int_grp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_econs_nph_uniaxial', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: baro_kin, baro_pot
    TYPE(cell_type), POINTER                 :: box
    TYPE(npt_info_type), POINTER             :: npt( :, : )
    TYPE(simulation_parameters_type), &
      POINTER                                :: simpar

    CALL get_md_env ( md_env, simpar = simpar, npt = npt, cell = box, error=error )

    CALL get_baro_energies ( box, simpar, npt, baro_kin, baro_pot, int_grp, error=error )
    md_ener % temp_baro = 2.0_dp * baro_kin
    cons = md_ener%ekin + md_ener%epot + baro_kin + baro_pot
  END SUBROUTINE get_econs_nph_uniaxial

  !******************************************************************************
  !!****** md_conserved_quantity/get_baro_energies [1.0] *
  !!
  !!   NAME
  !!     get_baro_energies
  !!
  !!   FUNCTION
  !!     Calculates kinetic energy and potential of barostat
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** *************************************************************************
  SUBROUTINE get_baro_energies ( box, simpar, npt, baro_kin, baro_pot, int_grp, error )

    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(npt_info_type), DIMENSION(:, :), &
      INTENT(IN)                             :: npt
    REAL(KIND=dp), INTENT(OUT)               :: baro_kin, baro_pot
    INTEGER, INTENT(IN)                      :: int_grp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_baro_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    REAL(dp)                                 :: iv0, v0, v_shock

    IF ( simpar%ensemble == npt_i_ensemble ) THEN
       baro_pot = simpar % p_ext * box % deth
       baro_kin = 0.5_dp * npt ( 1, 1 ) % v * npt ( 1, 1 ) % v * &
            npt ( 1, 1 ) % mass
    ELSE IF ( simpar%ensemble== npt_f_ensemble ) THEN
       baro_pot =  simpar % p_ext * box % deth
       baro_kin = 0.0_dp
       DO i = 1, 3
          DO j = 1, 3
             baro_kin = baro_kin + 0.5_dp * npt ( i, j ) % v *  &
                  npt ( i, j ) % v * npt ( i, j ) % mass
          END DO
       END DO
    ELSEIF ( simpar%ensemble == nph_ensemble .OR.&
         simpar%ensemble == nph_uniaxial_ensemble .OR.&
         simpar%ensemble == nph_uniaxial_damped_ensemble ) THEN
       v0 = simpar % v0
       iv0 = 1._dp / v0
       v_shock = simpar % v_shock

!!! Valid only for orthorhombic box !!!
       baro_pot = -0.5_dp * v_shock * v_shock * ( 1._dp - box % deth * iv0 )**2 - &
            simpar % p0 * ( v0 - box % deth )
!!! Valid only for orthorhombic box !!!
       baro_kin = 0.5_dp * npt ( 1, 1 ) % v * npt ( 1, 1 ) % v * &
            npt ( 1, 1 ) % mass
    END IF

  END SUBROUTINE get_baro_energies

  !******************************************************************************
  !!****** md_conserved_quantity/get_part_ke [1.0] *
  !!
  !!   NAME
  !!     get_part_ke
  !!
  !!   FUNCTION
  !!     Calculates kinetic energy of particles
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** *************************************************************************
  SUBROUTINE get_part_ke ( md_env, ekin, ekin_qm, ekin_shell, nfree_qm, &
             nfree_shell, nfree, int_grp, shell_present, temperature_per_kind, &
             nfree_kind, nfree_shell_kind, ekin_kind, ekin_shell_kind, error )
    TYPE(md_environment_type), POINTER       :: md_env
    REAL(KIND=dp), INTENT(OUT)               :: ekin, ekin_qm, ekin_shell
    INTEGER, INTENT(OUT)                     :: nfree_qm, nfree_shell
    INTEGER, INTENT(IN)                      :: nfree, int_grp
    LOGICAL, INTENT(OUT)                     :: shell_present
    LOGICAL, INTENT(IN)                      :: temperature_per_kind
    INTEGER, DIMENSION(:), INTENT(OUT), &
      OPTIONAL                               :: nfree_kind, nfree_shell_kind
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: ekin_kind, ekin_shell_kind
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_part_ke', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, iparticle, iparticle_kind, iparticle_local, nparticle_kind, &
      nparticle_local, shell_index
    LOGICAL                                  :: failure, is_shell, &
                                                shell_adiabatic
    REAL(KIND=dp)                            :: ekin_c, ekin_com, ekin_s, mass
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(shell_kind_type), POINTER           :: shell

    failure = .FALSE.
    
    NULLIFY(qmmm_env)
    CALL get_md_env ( md_env, force_env = force_env, error=error )
    CALL force_env_get (force_env, subsys = subsys, qmmm_env=qmmm_env, error=error)
    CALL cp_subsys_get(subsys=subsys,&
         atomic_kinds=atomic_kinds,&
         local_particles=local_particles,&
         particles=particles,shell_particles=shell_particles,&
         core_particles=core_particles,error=error)

    nparticle_kind = atomic_kinds%n_els
    atomic_kind_set => atomic_kinds%els

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                             shell_present=shell_present, &
                             shell_adiabatic=shell_adiabatic)

    nfree_qm = -HUGE(0)
    ekin     = zero
    ekin_qm  = zero
    ekin_shell  = zero
    ekin_s  = zero
    ekin_c  = zero
    ekin_com  = zero
    nfree_shell = -HUGE(0)
    IF(temperature_per_kind) THEN
      ekin_kind = zero
      nfree_kind = 0
      IF(shell_adiabatic) THEN
         ekin_shell_kind = zero
         nfree_shell_kind = 0
      END IF
    END IF
    
    particle_set => particles%els
    IF(shell_adiabatic) THEN
      shell_particle_set => shell_particles%els
      core_particle_set => core_particles%els
      nfree_shell = 3*SIZE(shell_particle_set)
      IF(temperature_per_kind) THEN
        DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass, &
                shell_active=is_shell, shell=shell)
           nparticle_local = local_particles%n_el(iparticle_kind)
           IF(is_shell) THEN
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                shell_index = particle_set(iparticle)%shell_index
                ekin_com =  0.5_dp * mass *  &
                     ( particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                     + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                     + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
                ekin = ekin +  ekin_com
                ekin_kind(iparticle_kind) = ekin_kind(iparticle_kind) + ekin_com
                nfree_kind(iparticle_kind) = nfree_kind(iparticle_kind) + 3
                ekin_c =  0.5_dp * shell%mass_core *  &
                   ( core_particle_set(shell_index)%v(1)* core_particle_set(shell_index)%v(1) &
                   + core_particle_set(shell_index)%v(2)* core_particle_set(shell_index)%v(2) &
                   + core_particle_set(shell_index)%v(3)* core_particle_set(shell_index)%v(3) )
                ekin_s =  0.5_dp * shell%mass_shell *  &
                   ( shell_particle_set(shell_index)%v(1)* shell_particle_set(shell_index)%v(1) &
                   + shell_particle_set(shell_index)%v(2)* shell_particle_set(shell_index)%v(2) &
                   + shell_particle_set(shell_index)%v(3)* shell_particle_set(shell_index)%v(3) )
                ekin_shell = ekin_shell + ekin_c + ekin_s - ekin_com
                ekin_shell_kind(iparticle_kind) = ekin_shell_kind(iparticle_kind) +&
                             ekin_c + ekin_s - ekin_com
                nfree_shell_kind(iparticle_kind) = nfree_shell_kind(iparticle_kind) + 3
             END DO  ! iparticle_local
           ELSE
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                ekin_com =  0.5_dp * mass *  &
                     ( particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                     + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                     + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
                ekin = ekin +  ekin_com
                ekin_kind(iparticle_kind) = ekin_kind(iparticle_kind) + ekin_com
                nfree_kind(iparticle_kind) = nfree_kind(iparticle_kind) + 3
              END DO  ! iparticle_local
           END IF
         END DO  ! iparticle_kind
         CALL mp_sum(ekin_kind, int_grp )
         CALL mp_sum(nfree_kind, int_grp )
         CALL mp_sum(ekin_shell_kind, int_grp )
         CALL mp_sum(nfree_shell_kind, int_grp )
      ELSE
        DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass, &
                shell_active=is_shell, shell=shell)
           nparticle_local = local_particles%n_el(iparticle_kind)
           IF(is_shell) THEN
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                shell_index =  particle_set(iparticle)%shell_index
                ekin_com =  0.5_dp * mass *  &
                     ( particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                     + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                     + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
                ekin = ekin +  ekin_com
                ekin_c =  0.5_dp * shell%mass_core *  &
                   ( core_particle_set(shell_index)%v(1)* core_particle_set(shell_index)%v(1) &
                   + core_particle_set(shell_index)%v(2)* core_particle_set(shell_index)%v(2) &
                   + core_particle_set(shell_index)%v(3)* core_particle_set(shell_index)%v(3) )
                ekin_s =  0.5_dp * shell%mass_shell *  &
                   ( shell_particle_set(shell_index)%v(1)* shell_particle_set(shell_index)%v(1) &
                   + shell_particle_set(shell_index)%v(2)* shell_particle_set(shell_index)%v(2) &
                   + shell_particle_set(shell_index)%v(3)* shell_particle_set(shell_index)%v(3) )
                ekin_shell = ekin_shell + ekin_c + ekin_s - ekin_com
             END DO  ! iparticle_local
           ELSE
             DO iparticle_local=1,nparticle_local
                iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                ekin_com =  0.5_dp * mass *  &
                     ( particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                     + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                     + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
                ekin = ekin +  ekin_com
              END DO  ! iparticle_local
           END IF
         END DO  ! iparticle_kind
      END IF   
      ! sum all contributions to energy over calculated parts on all processors
      CALL mp_sum ( ekin_shell, int_grp )
    ELSE
      IF(temperature_per_kind) THEN
        DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
              iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
              ekin_com =  0.5_dp * mass *  &
                   ( particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                   + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                   + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
              ekin = ekin +  ekin_com
              ekin_kind(iparticle_kind) = ekin_kind(iparticle_kind) + ekin_com
              nfree_kind(iparticle_kind) = nfree_kind(iparticle_kind) + 3
           END DO
        END DO  ! iparticle_kind
        CALL mp_sum(ekin_kind, int_grp )
        CALL mp_sum(nfree_kind, int_grp )
      ELSE
        DO iparticle_kind=1,nparticle_kind
           atomic_kind => atomic_kind_set(iparticle_kind)
           CALL get_atomic_kind(atomic_kind=atomic_kind,mass=mass)
           nparticle_local = local_particles%n_el(iparticle_kind)
           DO iparticle_local=1,nparticle_local
              iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
              ekin = ekin +  0.5_dp * mass *  &
                   ( particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                   + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                   + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
           END DO
        END DO
      END IF
    END IF

    ! sum all contributions to energy over calculated parts on all processors
    CALL mp_sum ( ekin, int_grp )
    !
    ! Compute the QM/MM kinetic energy
    IF (ASSOCIATED(qmmm_env)) THEN
       ! The degrees of freedom for the quantum part of the system
       ! are set to 3*Number of QM atoms and to simpar%nfree in case all the MM
       ! system is treated at QM level (not really QM/MM, just for consistency).
       ! The degree of freedom will not be correct if 1-3 atoms are treated only
       ! MM. In this case we should take care of rotations
       nfree_qm = 3*SIZE(qmmm_env%qm_atom_index)
       IF (nfree_qm == 3*SIZE(particle_set)) nfree_qm = nfree
       DO i = 1, SIZE(qmmm_env%qm_atom_index)
          iparticle = qmmm_env%qm_atom_index(i)
          mass = particle_set(iparticle)%atomic_kind%mass
          ekin_qm = ekin_qm +  0.5_dp * mass *  &
               ( particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
               + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
               + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3) )
       END DO
    END IF
  END SUBROUTINE get_part_ke

  !******************************************************************************
  !!****** md_conserved_quantity/get_coef_ke [1.0] *
  !!
  !!   NAME
  !!     get_coef_ke
  !!
  !!   FUNCTION
  !!     Calculates kinetic energy of coefs
  !!
  !!   AUTHOR
  !!     CJM
  !!
  !!   MODIFICATION HISTORY
  !!     gt 06-02-04 adapted to new dynamical coeffs type
  !!
  !!*****************************************************************************
  SUBROUTINE get_coef_ke ( dyn_coeff_set, ekinc, ctemp, int_group, error )

    TYPE(dyn_coeff_set_type), POINTER        :: dyn_coeff_set
    REAL(KIND=dp), INTENT(out)               :: ekinc, ctemp
    INTEGER                                  :: int_group
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_coef_ke', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: icoef, iel, ikind, ncoefs
    LOGICAL                                  :: gmass
    TYPE(dyn_coeff_type), POINTER            :: coeffs

    NULLIFY(coeffs)
    ekinc = zero
    ncoefs = 0
    gmass=(dyn_coeff_set%mass_storage==dyn_coeff_mass_global)

    DO ikind=1,SIZE(dyn_coeff_set%coeffs_of_kind)
       IF (ASSOCIATED(dyn_coeff_set%coeffs_of_kind(ikind)%coeffs)) THEN
          coeffs=>dyn_coeff_set%coeffs_of_kind(ikind)%coeffs
          DO iel = 1, coeffs%n_els
             DO icoef = 1, coeffs%ncoef_atom
                ncoefs=ncoefs+1
                IF(gmass) THEN
                   ekinc = ekinc + 0.5_dp*dyn_coeff_set%global_mass* &
                        coeffs%vel(iel,icoef)*coeffs%vel(iel,icoef)
                ELSE
                   ekinc = ekinc + 0.5_dp*coeffs%masses(iel,icoef)* &
                        coeffs%vel(iel,icoef)*coeffs%vel(iel,icoef)
                END IF
             END DO
          END DO
       END IF
    END DO

    ! sum all contributions to energy over calculated parts on all processors
    CALL mp_sum ( ekinc, int_group )
    CALL mp_sum ( ncoefs, int_group )

    ! get coefficient temperature

    ctemp=2.0_dp*ekinc/REAL(ncoefs,dp)

  END SUBROUTINE get_coef_ke
  !******************************************************************************

END MODULE md_conserved_quantities
!******************************************************************************
