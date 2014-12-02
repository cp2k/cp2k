!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Definition and initialisation of the mo data type.
!> \par History
!>      - adapted to the new QS environment data structure (02.04.2002,MK)
!>      - set_mo_occupation added (17.04.02,MK)
!>      - correct_mo_eigenvalues added (18.04.02,MK)
!>      - calculate_density_matrix moved from qs_scf to here (22.04.02,MK)
!>      - mo_set_p_type added (23.04.02,MK)
!>      - PRIVATE attribute set for TYPE mo_set_type (23.04.02,MK)
!>      - started conversion to LSD (1.2003, Joost VandeVondele)
!> \author Matthias Krack (09.05.2001,MK)
! *****************************************************************************
MODULE qs_mo_types

  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_init_p,&
                                             cp_dbcsr_release_p,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_copy_columns_hack
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                             fm_pool_create_fm
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE fermi_utils,                     ONLY: FermiFixed,&
                                             FermiFixedDeriv
  USE input_constants,                 ONLY: smear_energy_window,&
                                             smear_fermi_dirac,&
                                             smear_list
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE scf_control_types,               ONLY: smear_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  USE xas_env_types,                   ONLY: get_xas_env,&
                                             xas_environment_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mo_types'

  TYPE mo_set_type
    ! the actual MO coefficients as a matrix
    TYPE(cp_fm_type), POINTER                 :: mo_coeff
    TYPE(cp_dbcsr_type), POINTER              :: mo_coeff_b
    ! we are using the dbcsr mo_coeff_b
    LOGICAL                                   :: use_mo_coeff_b
    ! number of molecular orbitals (# cols in mo_coeff)
    INTEGER                                   :: nmo
    ! number of atomic orbitals (# rows in mo_coeff)
    INTEGER                                   :: nao
    ! occupation - eigenvalues  of the nmo states (if eigenstates)
    REAL(KIND = dp), DIMENSION(:), POINTER    :: eigenvalues,occupation_numbers
    ! maximum allowed occupation number of an MO (1 or 2)
    REAL(KIND = dp)                           :: maxocc
    ! number of electrons (taking occupation into account)
    INTEGER                                   :: nelectron
    REAL(KIND=dp)                             :: n_el_f
    ! highest non-zero occupied orbital
    INTEGER                                   :: homo
    ! lowest non maxocc occupied orbital (e.g. fractional or zero)
    INTEGER                                   :: lfomo
    ! flag that indicates if the MOS have the same occupation number
    LOGICAL                                   :: uniform_occupation
    ! the entropic energy contribution
    REAL(KIND=dp)                             :: kTS
    ! Fermi energy level
    REAL(KIND=dp)                             :: mu
    ! Threshold value for multiplicity change
    REAL(KIND=dp)                             :: flexible_electron_count
  END TYPE mo_set_type

  TYPE mo_set_p_type
    TYPE(mo_set_type), POINTER :: mo_set
  END TYPE mo_set_p_type

  PUBLIC :: mo_set_p_type,&
            mo_set_type

  PUBLIC :: allocate_mo_set,&
            correct_mo_eigenvalues,&
            deallocate_mo_set,&
            get_mo_set,&
            init_mo_set,&
            set_mo_occupation,&
            set_mo_set,&
            mo_set_restrict,&
            duplicate_mo_set

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
!> \brief allocate a new mo_set, and copy the old data
!> \param mo_set_new ...
!> \param mo_set_old ...
!> \param error ...
!> \date 2009-7-19
!> \par History
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE duplicate_mo_set(mo_set_new,mo_set_old,error)
    TYPE(mo_set_type), POINTER               :: mo_set_new, mo_set_old
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'duplicate_mo_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, nmo
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE (mo_set_new,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    mo_set_new%maxocc = mo_set_old%maxocc
    mo_set_new%nelectron = mo_set_old%nelectron
    mo_set_new%n_el_f = mo_set_old%n_el_f
    mo_set_new%nao = mo_set_old%nao
    mo_set_new%nmo = mo_set_old%nmo
    mo_set_new%homo = mo_set_old%homo
    mo_set_new%lfomo = mo_set_old%lfomo
    mo_set_new%uniform_occupation = mo_set_old%uniform_occupation
    mo_set_new%kTS = mo_set_old%kTS
    mo_set_new%mu = mo_set_old%mu
    mo_set_new%flexible_electron_count = mo_set_old%flexible_electron_count

    nmo = mo_set_new%nmo

    NULLIFY(mo_set_new%mo_coeff)
    CALL cp_fm_create(mo_set_new%mo_coeff,mo_set_old%mo_coeff%matrix_struct,error=error)
    CALL cp_fm_to_fm(mo_set_old%mo_coeff,mo_set_new%mo_coeff,error=error)

    NULLIFY(mo_set_new%mo_coeff_b)
    IF(ASSOCIATED(mo_set_old%mo_coeff_b)) THEN
       CALL cp_dbcsr_init_p(mo_set_new%mo_coeff_b,error=error)
       CALL cp_dbcsr_copy(mo_set_new%mo_coeff_b,mo_set_old%mo_coeff_b,error=error)
    ENDIF
    mo_set_new%use_mo_coeff_b = mo_set_old%use_mo_coeff_b

    ALLOCATE (mo_set_new%eigenvalues(nmo),STAT=istat)
    CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
    mo_set_new%eigenvalues = mo_set_old%eigenvalues

    ALLOCATE (mo_set_new%occupation_numbers(nmo),STAT=istat)
    CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
    mo_set_new%occupation_numbers = mo_set_old%occupation_numbers

  END SUBROUTINE duplicate_mo_set

! *****************************************************************************
!> \brief Allocates a mo set and partially initializes it (nao,nmo,nelectron,
!>        and flexible_electron_count are vaild).
!>        For the full initialization you need to call init_mo_set
!> \param mo_set the mo_set to allocate
!> \param nao number of atom orbitals
!> \param nmo number of molecular orbitals
!> \param nelectron number of electrons
!> \param n_el_f ...
!> \param maxocc maximum occupation of an orbital (LDA: 2, LSD:1)
!> \param flexible_electron_count the number of electrons can be changed
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \date 15.05.2001
!> \par History
!>      11.2002 splitted initialization in two phases [fawzi]
!> \author Matthias Krack
! *****************************************************************************
  SUBROUTINE allocate_mo_set(mo_set,nao,nmo,nelectron,n_el_f,maxocc,&
                             flexible_electron_count,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    INTEGER, INTENT(IN)                      :: nao, nmo, nelectron
    REAL(KIND=dp), INTENT(IN)                :: n_el_f, maxocc, &
                                                flexible_electron_count
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_mo_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(mo_set)) CALL deallocate_mo_set(mo_set,error)

    ALLOCATE (mo_set,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    mo_set%maxocc = maxocc
    mo_set%nelectron = nelectron
    mo_set%n_el_f = n_el_f
    mo_set%nao = nao
    mo_set%nmo = nmo
    mo_set%homo = 0
    mo_set%lfomo = 0
    mo_set%uniform_occupation = .TRUE.
    mo_set%kTS = 0.0_dp
    mo_set%mu = 0.0_dp
    mo_set%flexible_electron_count = flexible_electron_count

    NULLIFY (mo_set%eigenvalues)
    NULLIFY (mo_set%occupation_numbers)
    NULLIFY (mo_set%mo_coeff)
    NULLIFY (mo_set%mo_coeff_b)
    mo_set%use_mo_coeff_b = .FALSE.

  END SUBROUTINE allocate_mo_set

! *****************************************************************************
!> \brief initializes an allocated mo_set.
!>      eigenvalues, mo_coeff, occupation_numbers are valid only
!>      after this call.
!> \param mo_set the mo_set to initialize
!> \param fm_pool a pool out which you initialize the mo_set
!> \param name ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 rewamped [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE init_mo_set(mo_set,fm_pool,name,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_pool_type), POINTER           :: fm_pool
    CHARACTER(LEN=*), INTENT(in)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_mo_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat, nao, nmo
    LOGICAL                                  :: failure

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(mo_set),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(fm_pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(mo_set%eigenvalues),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(mo_set%occupation_numbers),cp_failure_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(mo_set%mo_coeff),cp_failure_level,routineP,error,failure)

    CALL fm_pool_create_fm(fm_pool,mo_set%mo_coeff,name=name,error=error)
    CALL cp_fm_get_info(mo_set%mo_coeff,nrow_global=nao,ncol_global=nmo,error=error)
    CPPostcondition((nao >= mo_set%nao),cp_failure_level,routineP,error,failure)
    CPPostcondition((nmo >= mo_set%nmo),cp_failure_level,routineP,error,failure)

    ALLOCATE (mo_set%eigenvalues(nmo),STAT=istat)
    CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
    mo_set%eigenvalues(:) = 0.0_dp

    ALLOCATE (mo_set%occupation_numbers(nmo),STAT=istat)
    CPPostcondition((istat == 0),cp_failure_level,routineP,error,failure)
    CALL set_mo_occupation_1(mo_set=mo_set,error=error)

  END SUBROUTINE init_mo_set

! *****************************************************************************
!> \brief make the beta orbitals explicitly equal to the alpha orbitals
!>       effectively copying the orbital data
!> \param mo_array ...
!> \param convert_dbcsr ...
!> \param error ...
!> \par History
!>      10.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE mo_set_restrict(mo_array,convert_dbcsr,error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    LOGICAL, INTENT(in), OPTIONAL            :: convert_dbcsr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mo_set_restrict', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure, my_convert_dbcsr

    CALL timeset(routineN,handle)

    failure = .FALSE.

    my_convert_dbcsr = .FALSE.
    IF(PRESENT(convert_dbcsr)) my_convert_dbcsr = convert_dbcsr

    CPPrecondition(ASSOCIATED(mo_array),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(mo_array).EQ.2,cp_failure_level,routineP,error,failure)
    CPPrecondition(mo_array(1)%mo_set%nmo>=mo_array(2)%mo_set%nmo,cp_failure_level,routineP,error,failure)

    ! first nmo_beta orbitals are copied from alpha to beta
    IF (.NOT.failure) THEN
       IF(my_convert_dbcsr) THEN!fm->dbcsr
          CALL cp_dbcsr_copy_columns_hack(mo_array(2)%mo_set%mo_coeff_b,mo_array(1)%mo_set%mo_coeff_b,&!fm->dbcsr
               mo_array(2)%mo_set%nmo,1,1,&!fm->dbcsr
               para_env=mo_array(1)%mo_set%mo_coeff%matrix_struct%para_env,&!fm->dbcsr
               blacs_env=mo_array(1)%mo_set%mo_coeff%matrix_struct%context,error=error)!fm->dbcsr
       ELSE!fm->dbcsr
          CALL cp_fm_to_fm(mo_array(1)%mo_set%mo_coeff,mo_array(2)%mo_set%mo_coeff,mo_array(2)%mo_set%nmo)
       ENDIF
    END IF

    CALL timestop(handle)

  END SUBROUTINE mo_set_restrict

! *****************************************************************************
!> \brief   Correct MO eigenvalues after MO level shifting.
!> \param mo_set ...
!> \param level_shift ...
!> \date    19.04.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE correct_mo_eigenvalues(mo_set,level_shift)

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(IN)                :: level_shift

    CHARACTER(LEN=*), PARAMETER :: routineN = 'correct_mo_eigenvalues', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imo

    IF (level_shift /= 0.0_dp) THEN
       DO imo=mo_set%homo+1,mo_set%nmo
          mo_set%eigenvalues(imo) = mo_set%eigenvalues(imo) - level_shift
       END DO
    END IF

  END SUBROUTINE correct_mo_eigenvalues

! *****************************************************************************
!> \brief   Deallocate a wavefunction data structure.
!> \param mo_set ...
!> \param error ...
!> \date    15.05.2001
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE deallocate_mo_set(mo_set,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_mo_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(mo_set)) THEN
      IF (ASSOCIATED(mo_set%eigenvalues)) THEN
        DEALLOCATE (mo_set%eigenvalues,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      IF (ASSOCIATED(mo_set%occupation_numbers)) THEN
        DEALLOCATE (mo_set%occupation_numbers,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      CALL cp_fm_release(mo_set%mo_coeff,error=error)
      IF(ASSOCIATED(mo_set%mo_coeff_b)) CALL cp_dbcsr_release_p(mo_set%mo_coeff_b, error=error)
      DEALLOCATE (mo_set,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE deallocate_mo_set

! *****************************************************************************
!> \brief   Get the components of a MO set data structure.
!> \param mo_set ...
!> \param maxocc ...
!> \param homo ...
!> \param lfomo ...
!> \param nao ...
!> \param nelectron ...
!> \param n_el_f ...
!> \param nmo ...
!> \param eigenvalues ...
!> \param occupation_numbers ...
!> \param mo_coeff ...
!> \param mo_coeff_b ...
!> \param uniform_occupation ...
!> \param kTS ...
!> \param mu ...
!> \param flexible_electron_count ...
!> \date    22.04.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE get_mo_set(mo_set,maxocc,homo,lfomo,nao,nelectron,n_el_f,nmo,&
                        eigenvalues,occupation_numbers,mo_coeff,mo_coeff_b,&
                        uniform_occupation,kTS,mu,flexible_electron_count)

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: maxocc
    INTEGER, INTENT(OUT), OPTIONAL           :: homo, lfomo, nao, nelectron
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: n_el_f
    INTEGER, INTENT(OUT), OPTIONAL           :: nmo
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eigenvalues, &
                                                occupation_numbers
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: mo_coeff
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: mo_coeff_b
    LOGICAL, INTENT(OUT), OPTIONAL           :: uniform_occupation
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: kTS, mu, &
                                                flexible_electron_count

    IF (PRESENT(maxocc)) maxocc = mo_set%maxocc
    IF (PRESENT(homo)) homo = mo_set%homo
    IF (PRESENT(lfomo)) lfomo = mo_set%lfomo
    IF (PRESENT(nao)) nao = mo_set%nao
    IF (PRESENT(nelectron)) nelectron = mo_set%nelectron
    IF (PRESENT(n_el_f)) n_el_f = mo_set%n_el_f
    IF (PRESENT(nmo)) nmo = mo_set%nmo
    IF (PRESENT(eigenvalues)) eigenvalues => mo_set%eigenvalues
    IF (PRESENT(occupation_numbers)) THEN
      occupation_numbers => mo_set%occupation_numbers
    END IF
    IF (PRESENT(mo_coeff)) mo_coeff => mo_set%mo_coeff
    IF (PRESENT(mo_coeff_b)) mo_coeff_b => mo_set%mo_coeff_b
    IF (PRESENT(uniform_occupation)) uniform_occupation = mo_set%uniform_occupation
    IF (PRESENT(kTS)) kTS = mo_set%kTS
    IF (PRESENT(mu)) mu = mo_set%mu
    IF (PRESENT(flexible_electron_count)) flexible_electron_count = mo_set%flexible_electron_count

  END SUBROUTINE get_mo_set

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

! *****************************************************************************
!> \brief   Set the components of a MO set data structure.
!> \param mo_set ...
!> \param maxocc ...
!> \param homo ...
!> \param lfomo ...
!> \param nao ...
!> \param nelectron ...
!> \param n_el_f ...
!> \param nmo ...
!> \param eigenvalues ...
!> \param occupation_numbers ...
!> \param uniform_occupation ...
!> \param kTS ...
!> \param mu ...
!> \param flexible_electron_count ...
!> \param error ...
!> \date    22.04.2002
!> \author  MK
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_mo_set(mo_set,maxocc,homo,lfomo,nao,nelectron,n_el_f,nmo,&
                        eigenvalues,occupation_numbers,uniform_occupation,&
                        kTS,mu,flexible_electron_count,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: maxocc
    INTEGER, INTENT(IN), OPTIONAL            :: homo, lfomo, nao, nelectron
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: n_el_f
    INTEGER, INTENT(IN), OPTIONAL            :: nmo
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eigenvalues, &
                                                occupation_numbers
    LOGICAL, INTENT(IN), OPTIONAL            :: uniform_occupation
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: kTS, mu, &
                                                flexible_electron_count
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_mo_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    IF (PRESENT(maxocc)) mo_set%maxocc = maxocc
    IF (PRESENT(homo)) mo_set%homo = homo
    IF (PRESENT(lfomo)) mo_set%lfomo = lfomo
    IF (PRESENT(nao)) mo_set%nao = nao
    IF (PRESENT(nelectron)) mo_set%nelectron = nelectron
    IF (PRESENT(n_el_f)) mo_set%n_el_f = n_el_f
    IF (PRESENT(nmo)) mo_set%nmo = nmo
    IF (PRESENT(eigenvalues)) THEN
      IF (ASSOCIATED(mo_set%eigenvalues)) THEN
        DEALLOCATE(mo_set%eigenvalues,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      mo_set%eigenvalues => eigenvalues
    END IF
    IF (PRESENT(occupation_numbers)) THEN
      IF (ASSOCIATED(mo_set%occupation_numbers)) THEN
        DEALLOCATE(mo_set%occupation_numbers,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF
      mo_set%occupation_numbers => occupation_numbers
    END IF
    IF(PRESENT(uniform_occupation)) mo_set%uniform_occupation = uniform_occupation
    IF(PRESENT(kTS)) mo_set%kTS = kTS
    IF(PRESENT(mu)) mo_set%mu = mu
    IF(PRESENT(flexible_electron_count)) mo_set%flexible_electron_count = flexible_electron_count

  END SUBROUTINE set_mo_set

END MODULE qs_mo_types
