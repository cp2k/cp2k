!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Set occupation of molecular orbitals
!> \par History
!>      - set_mo_occupation subroutines moved from qs_mo_types (11.12.2014 MI)
!> \author  MI
! *****************************************************************************

MODULE qs_mo_occupation

  USE fermi_utils,                     ONLY: FermiFixed,&
                                             FermiFixedDeriv
  USE input_constants,                 ONLY: smear_energy_window,&
                                             smear_fermi_dirac,&
                                             smear_list
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type,&
                                             set_mo_set
  USE scf_control_types,               ONLY: smear_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  USE xas_env_types,                   ONLY: get_xas_env,&
                                             xas_environment_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mo_occupation'

  PUBLIC :: set_mo_occupation


  INTERFACE set_mo_occupation
    MODULE PROCEDURE set_mo_occupation_1,set_mo_occupation_2
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief  Occupation for smeared spin polarized electronic structures
!>         with relaxed multiplicity
!>
!> \param mo_array ...
!> \param smear ...
!> \param error ...
!> \date    10.03.2011 (MI)
!> \author  MI
!> \version 1.0
! *****************************************************************************
   SUBROUTINE set_mo_occupation_3(mo_array,smear,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(smear_type), POINTER                :: smear
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_mo_occupation_3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: all_nmo, handle, homo_a, &
                                                homo_b, i, lfomo_a, lfomo_b, &
                                                nmo_a, nmo_b, stat, xas_estate
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: all_index
    LOGICAL                                  :: failure, is_large
    REAL(KIND=dp)                            :: all_nelec, kTS, mu, nelec_a, &
                                                nelec_b, occ_estate
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: all_eigval, all_occ
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigval_a, eigval_b, occ_a, &
                                                occ_b

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mo_array),cp_failure_level,routineP,error,failure)
    CPPrecondition((SIZE(mo_array) == 2),cp_failure_level,routineP,error,failure)
    CALL timeset(routineN,handle)

    NULLIFY(eigval_a,eigval_b,occ_a,occ_b)
    CALL get_mo_set(mo_set=mo_array(1)%mo_set,nmo=nmo_a,eigenvalues=eigval_a, &
         occupation_numbers=occ_a)
    CALL get_mo_set(mo_set=mo_array(2)%mo_set,nmo=nmo_b,eigenvalues=eigval_b, &
         occupation_numbers=occ_b)
    all_nmo = nmo_a+nmo_b
    ALLOCATE(all_eigval(all_nmo), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(all_occ(all_nmo), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(all_index(all_nmo), STAT=stat)

    all_eigval(1:nmo_a) = eigval_a(1:nmo_a)
    all_eigval(nmo_a+1:all_nmo) = eigval_b(1:nmo_b)

    CALL sort(all_eigval,all_nmo,all_index)

    xas_estate =  -1
    occ_estate = 0.0_dp

    nelec_a=0.0_dp
    nelec_b=0.0_dp
    all_nelec=0.0_dp
    nelec_a=accurate_sum(occ_a(:))
    nelec_b=accurate_sum(occ_b(:))
    all_nelec = nelec_a+nelec_b

    DO i = 1,all_nmo
       IF(all_index(i)<=nmo_a) THEN
         all_occ(i) = occ_a(all_index(i))
       ELSE
         all_occ(i) = occ_b(all_index(i)-nmo_a)
       END IF
    END DO

    CALL FermiFixed(all_occ,mu,kTS,all_eigval,all_nelec, &
                   smear%electronic_temperature,1._dp,xas_estate,occ_estate)

    is_large=ABS(MAXVAL(all_occ)-1.0_dp)> smear%eps_fermi_dirac
    ! this is not a real problem, but the temperature might be a bit large
        CALL cp_assert(.NOT.is_large,cp_warning_level,cp_assertion_failed,routineP,&
                       "Fermi-Dirac smearing includes the first MO",&
                       error,failure)

    is_large=ABS(MINVAL(all_occ))> smear%eps_fermi_dirac
    CALL cp_assert(.NOT.is_large,cp_warning_level,cp_assertion_failed,routineP,&
                   "Fermi-Dirac smearing includes the last MO => "//&
                    "Add more MOs for proper smearing.",error,failure)

    ! check that the total electron count is accurate
    is_large=(ABS(all_nelec - accurate_sum(all_occ(:))) > smear%eps_fermi_dirac*all_nelec)
    CALL cp_assert(.NOT.is_large,&
                   cp_warning_level,cp_assertion_failed,routineP,&
                   "Total number of electrons is not accurate",&
                    error,failure)

    DO i = 1,all_nmo
       IF(all_index(i)<=nmo_a) THEN
         occ_a(all_index(i)) = all_occ(i)
         eigval_a(all_index(i)) = all_eigval(i)
       ELSE
         occ_b(all_index(i)-nmo_a) = all_occ(i)
         eigval_b(all_index(i)-nmo_a) = all_eigval(i)
       END IF
    END DO

    nelec_a = accurate_sum(occ_a(:))
    nelec_b = accurate_sum(occ_b(:))

    DO i=1,nmo_a
       IF (occ_a(i) < 1.0_dp) THEN
          lfomo_a = i
          EXIT
       END IF
    END DO
    DO i=1,nmo_b
       IF (occ_b(i) < 1.0_dp) THEN
          lfomo_b = i
          EXIT
       END IF
    END DO
    homo_a = lfomo_a - 1
    DO i=nmo_a,lfomo_a,-1
       IF (occ_a(i) > smear%eps_fermi_dirac) THEN
          homo_a = i
          EXIT
       END IF
    END DO
    homo_b = lfomo_b - 1
    DO i=nmo_b,lfomo_b,-1
       IF (occ_b(i) > smear%eps_fermi_dirac) THEN
          homo_b = i
          EXIT
       END IF
    END DO

    CALL set_mo_set(mo_set=mo_array(1)%mo_set,kTS=kTS/2.0_dp,mu=mu,n_el_f=nelec_a,&
         lfomo=lfomo_a,homo=homo_a,uniform_occupation=.FALSE.,error=error)
    CALL set_mo_set(mo_set=mo_array(2)%mo_set,kTS=kTS/2.0_dp,mu=mu,n_el_f=nelec_b,&
         lfomo=lfomo_b,homo=homo_b,uniform_occupation=.FALSE.,error=error)

    CALL timestop(handle)

  END SUBROUTINE set_mo_occupation_3

! *****************************************************************************
!> \brief   Prepare an occupation of alpha and beta MOs following an Aufbau
!>          principle, i.e. allowing a change in multiplicity.
!> \param mo_array ...
!> \param smear ...
!> \param eval_deriv ...
!> \param error ...
!> \date    25.01.2010 (MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_mo_occupation_2(mo_array,smear,eval_deriv,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(smear_type), POINTER                :: smear
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eval_deriv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_mo_occupation_2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lumo_a, lumo_b, &
                                                multiplicity_new, &
                                                multiplicity_old, nelec
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: nelec_f, threshold
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigval_a, eigval_b
    TYPE(mo_set_type), POINTER               :: mo_set_a, mo_set_b

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mo_array),cp_failure_level,routineP,error,failure)
    mo_set_a => mo_array(1)%mo_set
    CPPrecondition(ASSOCIATED(mo_set_a),cp_failure_level,routineP,error,failure)
    ! Fall back for the case that we have only one MO set
    IF (SIZE(mo_array) == 1) THEN
       IF (PRESENT(eval_deriv)) THEN
          CALL set_mo_occupation_1(mo_set_a,smear=smear,eval_deriv=eval_deriv,error=error)
       ELSE
          CALL set_mo_occupation_1(mo_set_a,smear=smear,error=error)
       END IF
       CALL timestop(handle)
       RETURN
    END IF
    CPPrecondition((SIZE(mo_array) == 2),cp_failure_level,routineP,error,failure)
    mo_set_b => mo_array(2)%mo_set

    IF (smear%do_smear) THEN
       IF (smear%fixed_mag_mom < 0.0_dp) THEN
          IF (smear%fixed_mag_mom /= -1.0_dp) THEN
             CPPrecondition(.NOT.(PRESENT(eval_deriv)),cp_failure_level,routineP,error,failure)
             CALL set_mo_occupation_3(mo_array,smear=smear,error=error)
             CALL timestop(handle)
             RETURN
          END IF
       ELSE
          nelec_f = mo_set_a%n_el_f + mo_set_b%n_el_f
          IF (ABS((mo_set_a%n_el_f - mo_set_b%n_el_f) - smear%fixed_mag_mom) > smear%eps_fermi_dirac*nelec_f) THEN
             mo_set_a%n_el_f = nelec_f/2.0_dp + smear%fixed_mag_mom/2.0_dp
             mo_set_b%n_el_f = nelec_f/2.0_dp - smear%fixed_mag_mom/2.0_dp
          END IF
          CPPrecondition(.NOT.(PRESENT(eval_deriv)),cp_failure_level,routineP,error,failure)
          CALL set_mo_occupation_1(mo_set_a,smear=smear,error=error)
          CALL set_mo_occupation_1(mo_set_b,smear=smear,error=error)
       END IF
    END IF

    IF (.NOT.((mo_set_a%flexible_electron_count > 0.0_dp).AND.&
              (mo_set_b%flexible_electron_count > 0.0_dp))) THEN
       IF (PRESENT(eval_deriv)) THEN
          CALL set_mo_occupation_1(mo_set_a,smear=smear,eval_deriv=eval_deriv,error=error)
          CALL set_mo_occupation_1(mo_set_b,smear=smear,eval_deriv=eval_deriv,error=error)
       ELSE
          CALL set_mo_occupation_1(mo_set_a,smear=smear,error=error)
          CALL set_mo_occupation_1(mo_set_b,smear=smear,error=error)
       END IF
       CALL timestop(handle)
       RETURN
    END IF

    nelec = mo_set_a%nelectron + mo_set_b%nelectron

    multiplicity_old = mo_set_a%nelectron - mo_set_b%nelectron + 1

    CALL cp_assert((mo_set_a%nelectron < mo_set_a%nmo),&
                   cp_warning_level,cp_assertion_failed,routineP,&
                   "All alpha MOs are occupied. Add more alpha MOs to "//&
                   "allow for a higher multiplicity",only_ionode=.TRUE.,&
                   error=error)
    CALL cp_assert(((mo_set_b%nelectron < mo_set_b%nmo).OR.&
                    (mo_set_b%nelectron == mo_set_a%nelectron)),&
                   cp_warning_level,cp_assertion_failed,routineP,&
                   "All beta MOs are occupied. Add more beta MOs to "//&
                   "allow for a lower multiplicity",only_ionode=.TRUE.,&
                   error=error)

    eigval_a => mo_set_a%eigenvalues
    eigval_b => mo_set_b%eigenvalues

    lumo_a = 1
    lumo_b = 1

    ! Apply Aufbau principle
    DO i=1,nelec
       ! Threshold is needed to ensure a preference for alpha occupation in the case
       ! of degeneracy
       threshold = MAX(mo_set_a%flexible_electron_count,mo_set_b%flexible_electron_count)
       IF ((eigval_a(lumo_a) - threshold) < eigval_b(lumo_b)) THEN
          lumo_a = lumo_a + 1
        ELSE
          lumo_b = lumo_b + 1
       END IF
       IF (lumo_a > mo_set_a%nmo) THEN
          CALL cp_assert((i == nelec),&
                         cp_warning_level,cp_assertion_failed,routineP,&
                         "All alpha MOs are occupied. Add more alpha MOs to "//&
                         "allow for a higher multiplicity",only_ionode=.TRUE.,&
                         error=error)
          IF (i < nelec) THEN
             lumo_a = lumo_a - 1
             lumo_b = lumo_b + 1
          END IF
       END IF
       IF (lumo_b > mo_set_b%nmo) THEN
          CALL cp_assert((lumo_b >= lumo_a),&
                         cp_warning_level,cp_assertion_failed,routineP,&
                         "All beta MOs are occupied. Add more beta MOs to "//&
                         "allow for a lower multiplicity",only_ionode=.TRUE.,&
                         error=error)
          IF (i < nelec) THEN
             lumo_a = lumo_a + 1
             lumo_b = lumo_b - 1
          END IF
       END IF
    END DO

    mo_set_a%homo = lumo_a - 1
    mo_set_b%homo = lumo_b - 1

    IF (mo_set_b%homo > mo_set_a%homo) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "More beta ("//&
                      TRIM(ADJUSTL(cp_to_string(mo_set_b%homo)))//&
                      ") than alpha ("//&
                      TRIM(ADJUSTL(cp_to_string(mo_set_a%homo)))//&
                      ") MOs are occupied. Resorting to low spin state",&
                      only_ionode=.TRUE.,error=error)
       mo_set_a%homo = nelec/2 + MODULO(nelec,2)
       mo_set_b%homo = nelec/2
    END IF

    mo_set_a%nelectron = mo_set_a%homo
    mo_set_b%nelectron = mo_set_b%homo
    multiplicity_new = mo_set_a%nelectron - mo_set_b%nelectron + 1

    CALL cp_assert((multiplicity_new == multiplicity_old),&
                   cp_note_level,cp_assertion_failed,routineP,&
                   "Multiplicity changed from "//&
                   TRIM(ADJUSTL(cp_to_string(multiplicity_old)))//" to "//&
                   TRIM(ADJUSTL(cp_to_string(multiplicity_new))),&
                   only_ionode=.TRUE.,error=error,failure=failure)

    IF (PRESENT(eval_deriv)) THEN
       CALL set_mo_occupation_1(mo_set_a,smear=smear,eval_deriv=eval_deriv,error=error)
       CALL set_mo_occupation_1(mo_set_b,smear=smear,eval_deriv=eval_deriv,error=error)
    ELSE
       CALL set_mo_occupation_1(mo_set_a,smear=smear,error=error)
       CALL set_mo_occupation_1(mo_set_b,smear=smear,error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE set_mo_occupation_2

! *****************************************************************************
!> \brief   Smearing of the MO occupation with all kind of occupation numbers
!> \param   mo_set MO dataset structure
!> \param   smear optional smearing information
!> \param   eval_deriv on entry the derivative of the KS energy wrt to the occupation number
!>                     on exit  the derivative of the full free energy (i.e. KS and entropy) wrt to the eigenvalue
!> \param xas_env ...
!> \param error ...
!> \date    17.04.2002 (v1.0), 26.08.2008 (v1.1)
!> \author  Matthias Krack
!> \version 1.1
! *****************************************************************************
  SUBROUTINE set_mo_occupation_1(mo_set,smear,eval_deriv,xas_env, error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(smear_type), OPTIONAL, POINTER      :: smear
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eval_deriv
    TYPE(xas_environment_type), OPTIONAL, &
      POINTER                                :: xas_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_mo_occupation_1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i_first, imo, nmo, &
                                                nomo, stat, xas_estate
    LOGICAL                                  :: equal_size, failure, is_large
    REAL(KIND=dp)                            :: e1, e2, edelta, edist, &
                                                el_count, lengthscale, nelec, &
                                                occ_estate, xas_nelectron
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dfde

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(mo_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mo_set%eigenvalues),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(mo_set%occupation_numbers),cp_failure_level,routineP,error,failure)
    mo_set%occupation_numbers(:) = 0.0_dp

    ! Quick return, if no electrons are available
    IF (mo_set%nelectron == 0) THEN
      CALL timestop(handle)
      RETURN
    END IF

    xas_estate =  -1
    occ_estate = 0.0_dp
    IF(PRESENT(xas_env)) THEN
      CALL get_xas_env(xas_env=xas_env, xas_nelectron=xas_nelectron, occ_estate=occ_estate,xas_estate=xas_estate,error=error)
      nomo = INT(xas_nelectron + 1 - occ_estate)
      IF(MOD(xas_nelectron + 1 - occ_estate,1.0_dp)>EPSILON(0.0_dp)) nomo = nomo+1

      mo_set%occupation_numbers(1:nomo) = mo_set%maxocc
      IF(xas_estate>0) mo_set%occupation_numbers(xas_estate) = occ_estate
      el_count = SUM( mo_set%occupation_numbers(1:nomo))
      IF(el_count>xas_nelectron) mo_set%occupation_numbers(nomo) = mo_set%occupation_numbers(nomo) - (el_count-xas_nelectron)
      el_count = SUM( mo_set%occupation_numbers(1:nomo))
      CPPrecondition((el_count==xas_nelectron),cp_failure_level,routineP,error,failure)
    ELSE
     IF (MODULO(mo_set%nelectron,INT(mo_set%maxocc)) == 0) THEN
        nomo = NINT(mo_set%nelectron/mo_set%maxocc)
        ! Initialize MO occupations
        mo_set%occupation_numbers(1:nomo) = mo_set%maxocc
      ELSE
        nomo = INT(mo_set%nelectron/mo_set%maxocc) + 1
        ! Initialize MO occupations
        mo_set%occupation_numbers(1:nomo-1) = mo_set%maxocc
        mo_set%occupation_numbers(nomo) = mo_set%nelectron -(nomo-1)* mo_set%maxocc
      END IF
    END IF
    nmo  = SIZE(mo_set%eigenvalues)

    CPPrecondition((nmo >= nomo),cp_failure_level,routineP,error,failure)
    CPPrecondition((SIZE(mo_set%occupation_numbers) == nmo),cp_failure_level,routineP,error,failure)

    mo_set%homo = nomo
    mo_set%lfomo = nomo + 1
    mo_set%mu = mo_set%eigenvalues(nomo)

    ! Check consistency of the array lengths
    IF (PRESENT(eval_deriv)) THEN
      equal_size = (SIZE(mo_set%occupation_numbers,1) == SIZE(eval_deriv,1))
      CPPrecondition(equal_size,cp_failure_level,routineP,error,failure)
    END IF

    ! Quick return, if no smearing information is supplied (TO BE FIXED, smear should become non-optional...)
    IF (.NOT.PRESENT(smear)) THEN
      ! there is no dependence of the energy on the eigenvalues
      mo_set%uniform_occupation = .TRUE.
      IF (PRESENT(eval_deriv)) THEN
        eval_deriv = 0.0_dp
      END IF
      CALL timestop(handle)
      RETURN
    END IF

    ! Check if proper eigenvalues are already available
    IF (smear%method /= smear_list) THEN
      IF ((ABS(mo_set%eigenvalues(1)) < 1.0E-12_dp).AND.&
          (ABS(mo_set%eigenvalues(nmo)) < 1.0E-12_dp)) THEN
        CALL timestop(handle)
        RETURN
      END IF
    END IF

    ! Perform smearing
    IF (smear%do_smear) THEN
      IF(PRESENT(xas_env) ) THEN
         i_first = xas_estate + 1
         nelec = xas_nelectron
      ELSE
         i_first = 1
         IF (smear%fixed_mag_mom == -1.0_dp) THEN
            nelec = REAL(mo_set%nelectron,dp)
         ELSE
            nelec = mo_set%n_el_f
         END IF
      END IF
      SELECT CASE (smear%method)
      CASE (smear_fermi_dirac)
        IF (.NOT. PRESENT(eval_deriv)) THEN
          CALL FermiFixed(mo_set%occupation_numbers,mo_set%mu,mo_set%kTS,mo_set%eigenvalues,Nelec, &
                          smear%electronic_temperature,mo_set%maxocc,xas_estate,occ_estate)
        ELSE
          ! could be a relatively large matrix, but one could get rid of it by never storing it
          ! we only need dE/df * df/de, one could equally parallelize over entries, this could become expensive
          ALLOCATE(dfde(nmo,nmo),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! lengthscale could become a parameter, but this is pretty good
          lengthscale=10*smear%electronic_temperature

          CALL FermiFixedDeriv(dfde,mo_set%occupation_numbers,mo_set%mu,mo_set%kTS,mo_set%eigenvalues,Nelec, &
                               smear%electronic_temperature,mo_set%maxocc,xas_estate,occ_estate,lengthscale)

          ! deriv of E_{KS}-kT*S wrt to f_i
          eval_deriv=eval_deriv - mo_set%eigenvalues + mo_set%mu
          ! correspondingly the deriv of  E_{KS}-kT*S wrt to e_i
          eval_deriv=MATMUL(TRANSPOSE(dfde),eval_deriv)

          DEALLOCATE(dfde,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        END IF

        ! Find the lowest fractional occupied MO (LFOMO)
        DO imo=i_first,nmo
          IF (mo_set%occupation_numbers(imo) < mo_set%maxocc) THEN
            mo_set%lfomo = imo
            EXIT
          END IF
        END DO
        is_large=ABS(MAXVAL(mo_set%occupation_numbers)-mo_set%maxocc)> smear%eps_fermi_dirac
        ! this is not a real problem, but the temperature might be a bit large
        CALL cp_assert(.NOT.is_large,cp_warning_level,cp_assertion_failed,routineP,&
                       "Fermi-Dirac smearing includes the first MO",&
                       error,failure)

        ! Find the highest (fractional) occupied MO which will be now the HOMO
        DO imo=nmo,mo_set%lfomo,-1
          IF (mo_set%occupation_numbers(imo) > smear%eps_fermi_dirac) THEN
            mo_set%homo = imo
            EXIT
          END IF
        END DO
        is_large=ABS(MINVAL(mo_set%occupation_numbers))> smear%eps_fermi_dirac
        CALL cp_assert(.NOT.is_large,cp_warning_level,cp_assertion_failed,routineP,&
                       "Fermi-Dirac smearing includes the last MO => "//&
                       "Add more MOs for proper smearing.",error,failure)

        ! check that the total electron count is accurate
        is_large=(ABS(nelec - accurate_sum(mo_set%occupation_numbers(:))) > smear%eps_fermi_dirac*nelec)
        CALL cp_assert(.NOT.is_large,&
                       cp_warning_level,cp_assertion_failed,routineP,&
                       "Total number of electrons is not accurate",&
                       error,failure)

      CASE (smear_energy_window)
        ! not implemented
        CPPrecondition(.NOT.PRESENT(eval_deriv),cp_failure_level,routineP,error,failure)

        ! Define the energy window for the eigenvalues
        e1 = mo_set%eigenvalues(mo_set%homo) - 0.5_dp*smear%window_size
        CALL cp_assert((e1 > mo_set%eigenvalues(1)),cp_warning_level,cp_assertion_failed,routineP,&
                       "Energy window for smearing includes the first MO",&
                       error,failure)

        e2 = mo_set%eigenvalues(mo_set%homo) + 0.5_dp*smear%window_size
        CALL cp_assert((e2 < mo_set%eigenvalues(nmo)),cp_warning_level,cp_assertion_failed,routineP,&
                       "Energy window for smearing includes the last MO => "//&
                       "Add more MOs for proper smearing.",error,failure)

        ! Find the lowest fractional occupied MO (LFOMO)
        DO imo=i_first,nomo
          IF (mo_set%eigenvalues(imo) > e1) THEN
            mo_set%lfomo = imo
            EXIT
          END IF
        END DO

        ! Find the highest fractional occupied (non-zero) MO which will be the HOMO
        DO imo=nmo,nomo,-1
          IF (mo_set%eigenvalues(imo) < e2) THEN
            mo_set%homo = imo
            EXIT
          END IF
        END DO

        ! Get the number of electrons to be smeared
        edist = 0.0_dp
        nelec = 0.0_dp

        DO imo=mo_set%lfomo,mo_set%homo
          nelec = nelec + mo_set%occupation_numbers(imo)
          edist = edist + ABS(e2 - mo_set%eigenvalues(imo))
        END DO

        ! Smear electrons inside the energy window
        DO imo=mo_set%lfomo,mo_set%homo
          edelta = ABS(e2 - mo_set%eigenvalues(imo))
          mo_set%occupation_numbers(imo) = MIN(mo_set%maxocc,nelec*edelta/edist)
          nelec = nelec - mo_set%occupation_numbers(imo)
          edist = edist - edelta
        END DO

      CASE (smear_list)
        equal_size = SIZE(mo_set%occupation_numbers,1)==SIZE(smear%list,1)
        CPPrecondition(equal_size,cp_failure_level,routineP,error,failure)
        mo_set%occupation_numbers = smear%list
        ! there is no dependence of the energy on the eigenvalues
        IF (PRESENT(eval_deriv)) THEN
          eval_deriv = 0.0_dp
        END IF
        ! most general case
        mo_set%lfomo=1
        mo_set%homo =nmo
      END SELECT

      ! Check, if the smearing involves more than one MO
      IF (mo_set%lfomo == mo_set%homo) THEN
        mo_set%homo = nomo
        mo_set%lfomo = nomo + 1
      ELSE
        mo_set%uniform_occupation = .FALSE.
      END IF

    END IF ! do smear

    ! zeros don't count as uniform
    mo_set%uniform_occupation = ALL(mo_set%occupation_numbers==mo_set%maxocc)

    CALL timestop(handle)

  END SUBROUTINE set_mo_occupation_1

END MODULE qs_mo_occupation
