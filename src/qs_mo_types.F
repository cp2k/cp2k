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

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_init_p,&
                                             cp_dbcsr_release_p
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_copy_columns_hack
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                             fm_pool_create_fm
  USE cp_fm_types,                     ONLY: &
       cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, &
       cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, &
       cp_fm_type, cp_fm_write_unformatted
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_generate_filename,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE fermi_utils,                     ONLY: FermiFixed,&
                                             FermiFixedDeriv
  USE input_constants,                 ONLY: smear_energy_window,&
                                             smear_fermi_dirac,&
                                             smear_list
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast
  USE orbital_pointers,                ONLY: indco,&
                                             nco,&
                                             nso
  USE orbital_symbols,                 ONLY: cgf_symbol,&
                                             sgf_symbol
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE scf_control_types,               ONLY: smear_type
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: sort
  USE xas_env_types,                   ONLY: get_xas_env,&
                                             xas_environment_type
#include "cp_common_uses.h"

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
            wfn_restart_file_name,&
            read_mo_set,&
            read_mos_restart_low,&
            set_mo_occupation,&
            set_mo_set, &
            write_mo_set, &
            write_mo_set_low,&
            mo_set_restrict,&
            write_rt_mos_to_restart,&
            read_rt_mos_from_restart,&
            duplicate_mo_set

  INTERFACE read_mo_set
    MODULE PROCEDURE read_mo_set_from_restart
  END INTERFACE

  INTERFACE set_mo_occupation
    MODULE PROCEDURE set_mo_occupation_1,set_mo_occupation_2
  END INTERFACE

  INTERFACE write_mo_set
    MODULE PROCEDURE write_mo_set_to_output_unit,write_mo_set_to_restart
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief  Occupation for smeared spin polarized electronic structures
!>         with relaxed multiplicity
!>
!> \author  MI
!> \date    10.03.2011 (MI)
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
!> \author  Matthias Krack (MK)
!> \date    25.01.2010 (MK)
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
    REAL(KIND=dp)                            :: threshold
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigval_a, eigval_b
    TYPE(mo_set_type), POINTER               :: mo_set_a, mo_set_b

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
       RETURN
    END IF
    CPPrecondition((SIZE(mo_array) == 2),cp_failure_level,routineP,error,failure)
    mo_set_b => mo_array(2)%mo_set

    IF (smear%do_smear) THEN
       IF (smear%fixed_mag_mom < 0.0_dp) THEN
          IF (smear%fixed_mag_mom /= -1.0_dp) THEN
             CPPrecondition(.NOT.(PRESENT(eval_deriv)),cp_failure_level,routineP,error,failure)
             CALL set_mo_occupation_3(mo_array,smear=smear,error=error)
             RETURN
          END IF
       ELSE
          nelec = mo_set_a%n_el_f + mo_set_b%n_el_f
          IF (ABS((mo_set_a%n_el_f - mo_set_b%n_el_f) - smear%fixed_mag_mom) > smear%eps_fermi_dirac*nelec) THEN
             mo_set_a%n_el_f = nelec/2.0_dp + smear%fixed_mag_mom/2.0_dp
             mo_set_b%n_el_f = nelec/2.0_dp - smear%fixed_mag_mom/2.0_dp
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
       RETURN
    END IF

    CALL timeset(routineN,handle)

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
!> \par History
!> \author Joost VandeVondele
!> \date 2009-7-19
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
!> \param maxocc maximum occupation of an orbital (LDA: 2, LSD:1)
!> \param flexible_electron_count the number of electrons can be changed
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 splitted initialization in two phases [fawzi]
!> \author Matthias Krack
!> \date 15.05.2001
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
!> \param mo_se t
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
!> \author  MK
!> \date    19.04.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE correct_mo_eigenvalues(mo_set,level_shift)

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(IN)                :: level_shift

    CHARACTER(LEN=*), PARAMETER :: routineN = 'correct_mo_eigenvalues', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imo

    IF (level_shift == 0.0_dp) RETURN

    DO imo=mo_set%homo+1,mo_set%nmo
      mo_set%eigenvalues(imo) = mo_set%eigenvalues(imo) - level_shift
    END DO

  END SUBROUTINE correct_mo_eigenvalues

! *****************************************************************************
!> \brief   Deallocate a wavefunction data structure.
!> \author  MK
!> \date    15.05.2001
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
!> \author  MK
!> \date    22.04.2002
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
!> \author  Matthias Krack
!> \date    17.04.2002 (v1.0), 26.08.2008 (v1.1)
!> \version 1.1
!> \param   mo_set MO dataset structure
!> \param   smear optional smearing information
!> \param   eval_deriv on entry the derivative of the KS energy wrt to the occupation number
!>                     on exit  the derivative of the full free energy (i.e. KS and entropy) wrt to the eigenvalue
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

    INTEGER                                  :: handle, i, i_first, imo, nmo, &
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
    nmo  = SIZE(mo_set%eigenvalues)
    xas_estate =  -1
    occ_estate = 0.0_dp

    CPPrecondition((nmo >= nomo),cp_failure_level,routineP,error,failure)
    CPPrecondition((SIZE(mo_set%occupation_numbers) == nmo),cp_failure_level,routineP,error,failure)

    IF(PRESENT(xas_env) )THEN
      CALL get_xas_env(xas_env=xas_env,occ_estate=occ_estate,xas_estate=xas_estate,&
                  xas_nelectron=xas_nelectron,error=error)
      IF(xas_estate>-1) THEN
        mo_set%occupation_numbers(xas_estate) = occ_estate
        el_count = 0.0_dp
        DO i = 1,nomo
          el_count = el_count + mo_set%occupation_numbers(i)
          IF((el_count-xas_nelectron) >EPSILON(0.0_dp)) THEN
            nomo = i
            mo_set%occupation_numbers(nomo) = mo_set%occupation_numbers(nomo) - &
                 (el_count-xas_nelectron)
            EXIT
           END IF
        END DO
        IF((xas_nelectron - el_count) > EPSILON(0.0_dp) ) THEN
          nomo = nomo + 1
          mo_set%occupation_numbers(nomo) = xas_nelectron -  el_count
        END IF
      END IF
    ENDIF
    ! zeros don't count as uniform
    !MK mo_set%uniform_occupation = ALL(mo_set%occupation_numbers==mo_set%maxocc)

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

    CALL timestop(handle)

  END SUBROUTINE set_mo_occupation_1

! *****************************************************************************
!> \brief   Set the components of a MO set data structure.
!> \author  MK
!> \date    22.04.2002
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

! *****************************************************************************
  SUBROUTINE write_mo_set_to_restart(mo_array,particle_set,dft_section,&
                                     atomic_kind_set,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_to_restart', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=30), DIMENSION(2), PARAMETER :: &
      keys = (/"SCF%PRINT%RESTART_HISTORY","SCF%PRINT%RESTART        "/)

    INTEGER                                  :: handle, ikey, ires, ispin
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(1),error=error),cp_p_file) .OR.  &
         BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(2),error=error),cp_p_file) ) THEN

       IF(mo_array(1)%mo_set%use_mo_coeff_b) THEN
          ! we are using the dbcsr mo_coeff
          ! we copy it to the fm for anycase
          DO ispin=1,SIZE(mo_array)
             IF(.not.ASSOCIATED(mo_array(ispin)%mo_set%mo_coeff_b)) THEN
                CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
             ENDIF
             CALL copy_dbcsr_to_fm(mo_array(ispin)%mo_set%mo_coeff_b,&
                  mo_array(ispin)%mo_set%mo_coeff,error=error)!fm->dbcsr
          ENDDO
       ENDIF

       DO ikey=1,SIZE(keys)

          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               dft_section,keys(ikey),error=error),cp_p_file)) THEN

             ires = cp_print_key_unit_nr(logger,dft_section,keys(ikey),&
                  extension=".wfn", file_status="REPLACE", file_action="WRITE",&
                  do_backup=.TRUE., file_form="UNFORMATTED", error=error)

             CALL write_mo_set_low(mo_array, particle_set=particle_set, ires=ires, error=error)

             CALL cp_print_key_finished_output(ires,logger,dft_section,TRIM(keys(ikey)), error=error)
          END IF
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE write_mo_set_to_restart

! *****************************************************************************
  SUBROUTINE write_rt_mos_to_restart(mo_array,rt_mos,particle_set,dft_section,&
                                     atomic_kind_set,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rt_mos
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_rt_mos_to_restart', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=43), DIMENSION(2), PARAMETER :: keys = (/&
      "REAL_TIME_PROPAGATION%PRINT%RESTART_HISTORY",&
      "REAL_TIME_PROPAGATION%PRINT%RESTART        "/)

    INTEGER                                  :: handle, ikey, ires
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    logger => cp_error_get_logger(error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(1),error=error),cp_p_file) .OR.  &
         BTEST(cp_print_key_should_output(logger%iter_info,&
         dft_section,keys(2),error=error),cp_p_file) ) THEN

       DO ikey=1,SIZE(keys)

          IF (BTEST(cp_print_key_should_output(logger%iter_info,&
               dft_section,keys(ikey),error=error),cp_p_file)) THEN

             ires = cp_print_key_unit_nr(logger,dft_section,keys(ikey),&
                  extension=".rtpwfn", file_status="REPLACE", file_action="WRITE",&
                  do_backup=.TRUE., file_form="UNFORMATTED", error=error)

             CALL write_mo_set_low(mo_array, rt_mos=rt_mos, particle_set=particle_set, ires=ires, error=error)
             CALL cp_print_key_finished_output(ires,logger,dft_section,TRIM(keys(ikey)), error=error)
          END IF
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE write_rt_mos_to_restart

! *****************************************************************************
  SUBROUTINE write_mo_set_low(mo_array, particle_set, ires, rt_mos, error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: ires
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rt_mos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, imat, iset, ishell, ispin, istat, lmax, lshell, &
      max_block, nao, natom, nmo, nset, nset_max, nshell_max, nspin
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info
    LOGICAL                                  :: failure
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    CALL timeset(routineN,handle)
    nspin = SIZE(mo_array)
    nao = mo_array(1)%mo_set%nao

    IF (ires>0) THEN
       !     *** create some info about the basis set first ***
       natom = SIZE(particle_set,1)
       nset_max = 0
       nshell_max = 0

       DO iatom=1,natom
          NULLIFY(orb_basis_set,dftb_parameter)
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
               orb_basis_set=orb_basis_set,dftb_parameter=dftb_parameter)
          IF (ASSOCIATED(orb_basis_set)) THEN
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  nset=nset,&
                  nshell=nshell,&
                  l=l)
             nset_max = MAX(nset_max,nset)
             DO iset=1,nset
                nshell_max = MAX(nshell_max,nshell(iset))
             END DO
          ELSEIF (ASSOCIATED(dftb_parameter)) THEN
             CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
             nset_max = MAX(nset_max,1)
             nshell_max = MAX(nshell_max,lmax+1)
          ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "Unknown basis type. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          END IF
       END DO

       ALLOCATE (nso_info(nshell_max,nset_max,natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       nso_info(:,:,:) = 0

       ALLOCATE (nshell_info(nset_max,natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       nshell_info(:,:) = 0

       ALLOCATE (nset_info(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       nset_info(:) = 0

       DO iatom=1,natom
          NULLIFY(orb_basis_set,dftb_parameter)
          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
               orb_basis_set=orb_basis_set,dftb_parameter=dftb_parameter)
          IF (ASSOCIATED(orb_basis_set)) THEN
             CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                  nset=nset,&
                  nshell=nshell,&
                  l=l)
             nset_info(iatom) = nset
             DO iset=1,nset
                nshell_info(iset,iatom)=nshell(iset)
                DO ishell=1,nshell(iset)
                   lshell = l(ishell,iset)
                   nso_info(ishell,iset,iatom) = nso(lshell)
                END DO
             END DO
          ELSEIF (ASSOCIATED(dftb_parameter)) THEN
             CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
             nset_info(iatom) = 1
             nshell_info(1,iatom)=lmax+1
             DO ishell=1,lmax+1
                lshell = ishell-1
                nso_info(ishell,1,iatom) = nso(lshell)
             END DO
          ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
               "Unknown basis set type. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          END IF
       END DO

       WRITE (ires) natom,nspin,nao,nset_max,nshell_max
       WRITE (ires) nset_info
       WRITE (ires) nshell_info
       WRITE (ires) nso_info

       DEALLOCATE (nset_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DEALLOCATE (nshell_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

       DEALLOCATE (nso_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

    ! use the scalapack block size as a default for buffering columns
    CALL cp_fm_get_info(mo_array(1)%mo_set%mo_coeff,ncol_block=max_block,error=error)


    DO ispin=1,nspin
       nmo=mo_array(ispin)%mo_set%nmo
       IF ((ires>0).AND.(nmo > 0)) THEN
          WRITE (ires) nmo,&
               mo_array(ispin)%mo_set%homo,&
               mo_array(ispin)%mo_set%lfomo,&
               mo_array(ispin)%mo_set%nelectron
          WRITE (ires) mo_array(ispin)%mo_set%eigenvalues(1:nmo),&
               mo_array(ispin)%mo_set%occupation_numbers(1:nmo)
       END IF
       IF(PRESENT(rt_mos))THEN
          DO imat=2*ispin-1,2*ispin
             CALL cp_fm_write_unformatted(rt_mos(imat)%matrix,ires,error)
          END DO
       ELSE
          CALL cp_fm_write_unformatted(mo_array(ispin)%mo_set%mo_coeff,ires,error)
       END IF
    END DO

    CALL timestop(handle)

  END SUBROUTINE write_mo_set_low

! *****************************************************************************
  SUBROUTINE wfn_restart_file_name(filename,exist,section,logger,xas,rtp,error)
    CHARACTER(LEN=default_path_length), &
      INTENT(OUT)                            :: filename
    LOGICAL, INTENT(OUT)                     :: exist
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL, INTENT(IN), OPTIONAL            :: xas, rtp
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: n_rep_val
    LOGICAL                                  :: my_rtp, my_xas
    TYPE(section_vals_type), POINTER         :: print_key

    my_xas = .FALSE.
    my_rtp = .FALSE.
    IF(PRESENT(xas)) my_xas = xas
    IF(PRESENT(rtp)) my_rtp = rtp

    exist = .FALSE.
    CALL section_vals_val_get(section,"WFN_RESTART_FILE_NAME",n_rep_val=n_rep_val,error=error)
    IF (n_rep_val>0) THEN
      CALL section_vals_val_get(section,"WFN_RESTART_FILE_NAME",c_val=filename,error=error)
    ELSE
      IF(my_xas) THEN
       ! try to read from the filename that is generated automatically from the printkey
        print_key => section_vals_get_subs_vals(section,"PRINT%RESTART",error=error)
        filename = cp_print_key_generate_filename(logger,print_key, &
                    extension="",my_local=.FALSE., error=error)
      ELSE IF (my_rtp)THEN
       ! try to read from the filename that is generated automatically from the printkey
        print_key => section_vals_get_subs_vals(section,"REAL_TIME_PROPAGATION%PRINT%RESTART",error=error)
        filename = cp_print_key_generate_filename(logger,print_key, &
                    extension=".rtpwfn",my_local=.FALSE., error=error)
      ELSE
        ! try to read from the filename that is generated automatically from the printkey
        print_key => section_vals_get_subs_vals(section,"SCF%PRINT%RESTART",error=error)
        filename = cp_print_key_generate_filename(logger,print_key, &
                    extension=".wfn", my_local=.FALSE., error=error)
      END IF
    ENDIF
    IF(.NOT.my_xas) THEN
      INQUIRE(FILE=filename,exist=exist)
    END IF

  END SUBROUTINE wfn_restart_file_name

! *****************************************************************************
  SUBROUTINE read_mo_set_from_restart(mo_array,atomic_kind_set,particle_set,&
       para_env,id_nr,multiplicity,dft_section,natom_mismatch,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: id_nr, multiplicity
    TYPE(section_vals_type), POINTER         :: dft_section
    LOGICAL, INTENT(OUT), OPTIONAL           :: natom_mismatch
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mo_set_from_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: group, handle, ispin, natom, &
                                                nspin, restart_unit, source
    LOGICAL                                  :: exist, failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    failure = .FALSE.

    nspin = SIZE(mo_array)
    restart_unit = -1

    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN

      natom = SIZE(particle_set,1)
      CALL wfn_restart_file_name(file_name,exist,dft_section,logger,error=error)
      IF (id_nr/=0) THEN
         ! Is it one of the backup files?
         file_name = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=restart_unit)

    END IF

    CALL read_mos_restart_low (mo_array,para_env=para_env, particle_set=particle_set,&
         natom=natom,&
         rst_unit=restart_unit, multiplicity=multiplicity, natom_mismatch=natom_mismatch, &
         error=error)
    IF (PRESENT(natom_mismatch)) THEN
       ! read_mos_restart_low only the io_node returns natom_mismatch, must broadcast it
       CALL mp_bcast(natom_mismatch,source,group)
       IF (natom_mismatch) THEN
          IF (para_env%ionode) CALL close_file(unit_number=restart_unit)
          CALL timestop(handle)
          RETURN
       ENDIF
    ENDIF

    ! Close restart file
    IF (para_env%ionode) CALL close_file(unit_number=restart_unit)

    DO ispin = 1,nspin
      CALL write_mo_set(mo_array(ispin)%mo_set,atomic_kind_set,particle_set,4,&
                        dft_section,error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE read_mo_set_from_restart
! *****************************************************************************
  SUBROUTINE read_rt_mos_from_restart(mo_array,rt_mos,atomic_kind_set,particle_set,&
       para_env,id_nr,multiplicity,dft_section,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rt_mos
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: id_nr, multiplicity
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_rt_mos_from_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: group, handle, ispin, natom, &
                                                nspin, restart_unit, source
    LOGICAL                                  :: exist, failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    failure = .FALSE.

    nspin = SIZE(mo_array)
    restart_unit = -1

    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN

      natom = SIZE(particle_set,1)
      CALL wfn_restart_file_name(file_name,exist,dft_section,logger,rtp=.TRUE.,error=error)
      IF (id_nr/=0) THEN
         ! Is it one of the backup files?
         file_name = TRIM(file_name)//".bak-"//ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=restart_unit)

    END IF

    CALL read_mos_restart_low (mo_array, rt_mos=rt_mos,para_env= para_env,&
         particle_set=particle_set, natom=natom,&
         rst_unit=restart_unit, multiplicity=multiplicity, error=error)

    ! Close restart file
    IF (para_env%ionode) CALL close_file(unit_number=restart_unit)

    DO ispin = 1,nspin
      CALL write_mo_set(mo_array(ispin)%mo_set,atomic_kind_set,particle_set,4,&
                        dft_section,error=error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE read_rt_mos_from_restart
! *****************************************************************************
!> \brief Reading the mos from apreviously defined restart file
!> \par History
!>      12.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE read_mos_restart_low (mos,  para_env, particle_set, natom, rst_unit, &
        multiplicity, rt_mos, natom_mismatch, error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: natom, rst_unit
    INTEGER, INTENT(in), OPTIONAL            :: multiplicity
    TYPE(cp_fm_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rt_mos
    LOGICAL, INTENT(OUT), OPTIONAL           :: natom_mismatch
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mos_restart_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: group, homo, homo_read, i, iatom, imat, irow, iset, iset_read, &
      ishell, ishell_read, iso, ispin, istat, lfomo_read, lmax, lshell, &
      my_mult, nao, nao_read, natom_read, nelectron, nelectron_read, nmo, &
      nmo_read, nnshell, nset, nset_max, nshell_max, nspin, nspin_read, &
      offset_read, source
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info, offset_info
    LOGICAL                                  :: failure, minbas, natom_match, &
                                                use_this
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eig_read, occ_read
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer, vecbuffer_read
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    logger => cp_error_get_logger(error)

    nspin = SIZE(mos)
    nao = mos(1)%mo_set%nao
    my_mult = 0
    IF(PRESENT(multiplicity)) my_mult = multiplicity
    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN
       READ (rst_unit) natom_read,nspin_read,nao_read,nset_max,nshell_max
       IF(PRESENT(rt_mos)) THEN
          IF (nspin_read /= nspin) THEN
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                 "To change nspin is not possible. "//&
                 CPSourceFileRef,&
                 only_ionode=.TRUE.)
         END IF
      ELSE
         ! we should allow for restarting with different spin settings
         IF (nspin_read /= nspin) THEN
            WRITE(cp_logger_get_default_unit_nr(logger),*)  &
                 "READ RESTART : WARNING : nspin is not equal "
         END IF
         ! this case needs fixing of homo/lfomo/nelec/occupations ...
         IF (nspin_read > nspin) THEN
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                 "Reducing nspin is not possible. "//&
                 CPSourceFileRef,&
                 only_ionode=.TRUE.)
         ENDIF
      END IF

      natom_match = (natom_read == natom)

      IF (natom_match) THEN ! actually do the read read

         ! Let's make it possible to change the basis set
         ALLOCATE (nso_info(nshell_max,nset_max,natom_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE (nshell_info(nset_max,natom_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE (nset_info(natom_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE (offset_info(nshell_max,nset_max,natom_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

         IF (nao_read /= nao) THEN
            WRITE(cp_logger_get_default_unit_nr(logger),*) &
                 " READ RESTART : WARNING : DIFFERENT # AOs ",nao,nao_read
            IF(PRESENT(rt_mos))&
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                    "To change basis is not possible. "//&
                    CPSourceFileRef,&
                    only_ionode=.TRUE.)
         END IF

         READ (rst_unit) nset_info
         READ (rst_unit) nshell_info
         READ (rst_unit) nso_info

         i=1
         DO iatom=1,natom
            DO iset=1,nset_info(iatom)
               DO ishell=1,nshell_info(iset,iatom)
                  offset_info(ishell,iset,iatom) = i
                  i=i+nso_info(ishell,iset,iatom)
               END DO
            END DO
         END DO

         ALLOCATE(vecbuffer_read(1,nao_read),STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

      END IF ! natom_match
   END IF ! ionode

   ! make natom_match and natom_mismatch uniform across all nodes
   CALL mp_bcast(natom_match,source,group)
   IF (PRESENT(natom_mismatch)) natom_mismatch = .NOT. natom_match
   ! handle natom_match false
   IF (.NOT. natom_match) THEN
      IF (PRESENT(natom_mismatch)) THEN
         WRITE(cp_logger_get_default_unit_nr(logger),*) &
              " READ RESTART : WARNING : DIFFERENT natom, returning ",natom,natom_read
         RETURN
      ELSE
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
              "Incorrect number of atoms in restart file. "//&
              CPSourceFileRef,&
              only_ionode=.TRUE.)
      ENDIF
   ENDIF

   CALL mp_bcast(nspin_read,source,group)

   ALLOCATE (vecbuffer(1,nao),STAT=istat)
   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

   DO ispin=1,nspin

      nmo=mos(ispin)%mo_set%nmo
      homo=mos(ispin)%mo_set%homo
      mos(ispin)%mo_set%eigenvalues(:) = 0.0_dp
      mos(ispin)%mo_set%occupation_numbers(:) = 0.0_dp
      CALL cp_fm_set_all(mos(ispin)%mo_set%mo_coeff,0.0_dp,error=error)

      IF (para_env%ionode.AND.(nmo > 0)) THEN
         READ (rst_unit) nmo_read, homo_read, lfomo_read, nelectron_read
         ALLOCATE(eig_read(nmo_read), occ_read(nmo_read), STAT=istat)
         CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
         eig_read = 0.0_dp
         occ_read = 0.0_dp

         nmo = MIN(nmo,nmo_read)
         CALL cp_assert((nmo_read >= nmo),cp_warning_level,cp_assertion_failed,routineP,&
              "The number of MOs on the restart unit is smaller than the number of "//&
              "the allocated MOs. The MO set will be padded with zeros!"//&
              CPSourceFileRef,&
              only_ionode=.TRUE.)
        CALL cp_assert((nmo_read<=nmo),cp_warning_level,cp_assertion_failed,routineP,&
             "The number of MOs on the restart unit is greater than the number of "//&
             "the allocated MOs. The read MO set will be truncated!"//&
             CPSourceFileRef,&
             only_ionode=.TRUE.)

        READ (rst_unit) eig_read(1:nmo_read), occ_read(1:nmo_read)
        mos(ispin)%mo_set%eigenvalues(1:nmo) = eig_read(1:nmo)
        mos(ispin)%mo_set%occupation_numbers(1:nmo) = occ_read(1:nmo)
        DEALLOCATE(eig_read, occ_read, STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

        mos(ispin)%mo_set%homo = homo_read
        mos(ispin)%mo_set%lfomo = lfomo_read
        IF (homo_read > nmo) THEN
          IF(nelectron_read==mos(ispin)%mo_set%nelectron) THEN
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                 "The number of occupied MOs on the restart unit is larger than "//&
                 "the allocated MOs. The read MO set will be truncated and the occupation numbers recalculated!"//&
CPSourceFileRef,&
                 only_ionode=.TRUE.)
             CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,error=error)
          ELSE
              ! can not make this a warning i.e. homo must be smaller than nmo
              ! otherwise e.g. set_mo_occupation will go out of bounds
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                "Number of occupied MOs on restart unit larger than allocated MOs. "//&
CPSourceFileRef,&
                 only_ionode=.TRUE.)
          END IF
        END IF
      END IF

      CALL mp_bcast(nmo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%homo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%lfomo,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%nelectron,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%eigenvalues,source,group)
      CALL mp_bcast(mos(ispin)%mo_set%occupation_numbers,source,group)
      IF(PRESENT(rt_mos))THEN
         DO imat=2*ispin-1,2*ispin
            DO i=1,nmo
               IF (para_env%ionode) THEN
                  READ (rst_unit) vecbuffer
               ELSE
                  vecbuffer(1,:) = 0.0_dp
               END IF
               CALL mp_bcast(vecbuffer,source,group)
               CALL cp_fm_set_submatrix(rt_mos(imat)%matrix,&
                    vecbuffer,1,i,nao,1,transpose=.TRUE.,error=error)
            END DO
         END DO
      ELSE
         DO i=1,nmo
            IF (para_env%ionode) THEN
               READ (rst_unit) vecbuffer_read
               ! now, try to assign the read to the real vector
               ! in case the basis set changed this involves some guessing
               irow=1
               DO iatom=1,natom
                  NULLIFY(orb_basis_set,dftb_parameter,l,nshell)
                  CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                       orb_basis_set=orb_basis_set,dftb_parameter=dftb_parameter)
                  IF (ASSOCIATED(orb_basis_set)) THEN
                     CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                          nset=nset,&
                          nshell=nshell,&
                          l=l)
                     minbas=.FALSE.
                  ELSEIF (ASSOCIATED(dftb_parameter)) THEN
                     CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
                     nset = 1
                     minbas=.TRUE.
                  ELSE
                     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                          "Unknown basis set type. "//&
                          CPSourceFileRef,&
                          only_ionode=.TRUE.)
                  END IF

                  use_this = .TRUE.
                  iset_read = 1
                  DO iset=1,nset
                     ishell_read = 1
                     IF(minbas) THEN
                        nnshell = lmax+1
                     ELSE
                        nnshell = nshell(iset)
                     END IF
                     DO ishell=1,nnshell
                        IF(minbas) THEN
                           lshell = ishell-1
                        ELSE
                           lshell = l(ishell,iset)
                        END IF
                        IF (iset_read > nset_info(iatom)) use_this = .FALSE.
                        IF (use_this) THEN ! avoids out of bound access of the lower line if false
                           IF (nso(lshell) == nso_info(ishell_read,iset_read,iatom)) THEN
                              offset_read=offset_info(ishell_read,iset_read,iatom)
                              ishell_read=ishell_read+1
                              IF (ishell_read > nshell_info(iset,iatom)) THEN
                                 ishell_read = 1
                                 iset_read = iset_read+1
                              END IF
                           ELSE
                              use_this = .FALSE.
                           END IF
                        END IF
                        DO iso=1,nso(lshell)
                           IF (use_this) THEN
                              IF (offset_read-1+iso.LT.1 .OR. offset_read-1+iso.GT.nao_read) THEN
                                 vecbuffer(1,irow)=0.0_dp
                              ELSE
                                 vecbuffer(1,irow)=vecbuffer_read(1,offset_read-1+iso)
                              END IF
                           ELSE
                              vecbuffer(1,irow) = 0.0_dp
                           END IF
                           irow = irow + 1
                        END DO
                        use_this = .TRUE.
                     END DO
                  END DO
               END DO

            ELSE

               vecbuffer(1,:) = 0.0_dp

            END IF

            CALL mp_bcast(vecbuffer,source,group)
            CALL cp_fm_set_submatrix(mos(ispin)%mo_set%mo_coeff,&
                 vecbuffer,1,i,nao,1,transpose=.TRUE.,error=error)
         END DO
      END IF
      ! Skip extra MOs if there any
      IF (para_env%ionode) THEN
        !ignore nmo = 0
        IF(nmo>0) THEN
          DO i=nmo+1,nmo_read
            READ (rst_unit) vecbuffer_read
          END DO
        END IF
      END IF

      IF(.NOT.PRESENT(rt_mos)) THEN
        IF(ispin==1 .AND. nspin_read<nspin) THEN

          mos(ispin+1)%mo_set%homo = mos(ispin)%mo_set%homo
          mos(ispin+1)%mo_set%lfomo = mos(ispin)%mo_set%lfomo
          nelectron = mos(ispin)%mo_set%nelectron
          IF (my_mult.NE.1) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "Restarting an LSD calculation from an LDA wfn only works for multiplicity=1 (singlets)."//&
                  CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
          IF(mos(ispin+1)%mo_set%nelectron < 0)THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                  "LSD: too few electrons for this multiplisity. "//&
                  CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
          mos(ispin+1)%mo_set%eigenvalues = mos(ispin)%mo_set%eigenvalues
          mos(ispin)%mo_set%occupation_numbers = mos(ispin)%mo_set%occupation_numbers/2.0_dp
          mos(ispin+1)%mo_set%occupation_numbers = mos(ispin)%mo_set%occupation_numbers
          CALL cp_fm_to_fm(mos(ispin)%mo_set%mo_coeff,mos(ispin+1)%mo_set%mo_coeff,error=error)
          EXIT
       END IF
    END IF
 END DO   ! ispin

    DEALLOCATE(vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    IF (para_env%ionode) THEN
       DEALLOCATE(vecbuffer_read,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(offset_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(nso_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(nshell_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(nset_info,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE read_mos_restart_low

! *****************************************************************************
!> \brief   Write the MO eigenvalues, MO occupation numbers and
!>          MO mo_coeff.
!> \author  MK
!> \date    15.05.2001
!> \par History:
!>       - Optionally print Cartesian MOs (20.04.2005,MK)
!> \par Variables
!>       - after : Number of digits after point.
!>       - before: Number of digits before point.
!> \version 1.0
! *****************************************************************************
  SUBROUTINE write_mo_set_to_output_unit(mo_set,atomic_kind_set,particle_set,&
       before,dft_section,spin,last,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: before
    TYPE(section_vals_type), POINTER         :: dft_section
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: spin
    LOGICAL, INTENT(IN), OPTIONAL            :: last
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_to_output_unit', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: symbol
    CHARACTER(LEN=12), DIMENSION(:), POINTER :: bcgf_symbol
    CHARACTER(LEN=16)                        :: fmtstr5, fmtstr6
    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=2*default_string_length)   :: name
    CHARACTER(LEN=22)                        :: fmtstr2
    CHARACTER(LEN=25)                        :: fmtstr1, fmtstr7
    CHARACTER(LEN=27)                        :: fmtstr4
    CHARACTER(LEN=38)                        :: fmtstr3
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: bsgf_symbol
    INTEGER :: after, first_mo, from, iatom, icgf, ico, icol, imo, irow, &
      iset, isgf, ishell, iso, istat, iw, jcol, last_mo, left, lmax, lshell, &
      natom, ncgf, ncol, ncol_global, nrow_global, nset, nsgf, right, &
      scf_step, to, width
    INTEGER, DIMENSION(:), POINTER           :: mo_index_range, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l
    LOGICAL                                  :: failure, ionode, my_last, &
                                                p_cart, p_eval, p_evec, &
                                                p_occ, should_output
    REAL(KIND=dp)                            :: gap
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cmatrix, smatrix
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    NULLIFY (bcgf_symbol)
    NULLIFY (bsgf_symbol)
    NULLIFY (logger)
    NULLIFY (mo_index_range)
    NULLIFY (nshell)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    failure = .FALSE.
    CALL section_vals_val_get(dft_section,"PRINT%MO%EIGENVALUES",l_val=p_eval,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%EIGENVECTORS",l_val=p_evec,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%OCCUPATION_NUMBERS",l_val=p_occ,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%CARTESIAN",l_val=p_cart,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%MO_INDEX_RANGE",i_vals=mo_index_range,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%NDIGITS",i_val=after,error=error)
    after = MIN(MAX(after,1),16)
    should_output = BTEST(cp_print_key_should_output(logger%iter_info,dft_section,&
         "PRINT%MO",error=error),cp_p_file)

    IF ((.NOT.should_output).OR.(.NOT.(p_eval.OR.p_evec.OR.p_occ))) RETURN

    IF (PRESENT(last)) THEN
       my_last = last
    ELSE
       my_last = .FALSE.
    END IF

    scf_step = logger%iter_info%iteration(logger%iter_info%n_rlevel) - 1

    IF (p_evec) THEN
       CALL cp_fm_get_info(mo_set%mo_coeff,&
            nrow_global=nrow_global,&
            ncol_global=ncol_global,error=error)
       ALLOCATE(smatrix(nrow_global,ncol_global),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL cp_fm_get_submatrix(mo_set%mo_coeff,smatrix,error=error)
       IF (.NOT.ionode) THEN
          DEALLOCATE(smatrix,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

    iw = cp_print_key_unit_nr(logger,dft_section,"PRINT%MO",&
         ignore_should_output=should_output,&
         extension=".MOLog",error=error)

    IF (iw > 0) THEN

       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
            natom=natom,&
            ncgf=ncgf,&
            nsgf=nsgf)

       ! Definition of the variable formats

       fmtstr1 = "(/,T2,21X,  (  X,I5,  X))"
       fmtstr2 = "(T2,21X,  (1X,F  .  ))"
       fmtstr3 = "(T2,I5,1X,I5,1X,A,1X,A6,  (1X,F  .  ))"

       width = before + after + 3
       ncol = INT(56/width)

       right = MAX((after-2),1)
       left =  width - right - 5

       WRITE (UNIT=fmtstr1(11:12),FMT="(I2)") ncol
       WRITE (UNIT=fmtstr1(14:15),FMT="(I2)") left
       WRITE (UNIT=fmtstr1(21:22),FMT="(I2)") right

       WRITE (UNIT=fmtstr2(9:10),FMT="(I2)") ncol
       WRITE (UNIT=fmtstr2(16:17),FMT="(I2)") width - 1
       WRITE (UNIT=fmtstr2(19:20),FMT="(I2)") after

       WRITE (UNIT=fmtstr3(25:26),FMT="(I2)") ncol
       WRITE (UNIT=fmtstr3(32:33),FMT="(I2)") width - 1
       WRITE (UNIT=fmtstr3(35:36),FMT="(I2)") after

       IF (p_evec) THEN

          IF (p_cart) THEN

             IF (my_last) THEN
                name = "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "CARTESIAN MO EIGENVECTORS"
             ELSE
                WRITE (UNIT=name,FMT="(A,I6)")&
                     "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "CARTESIAN MO EIGENVECTORS AFTER SCF STEP",scf_step
             END IF

             ALLOCATE(cmatrix(ncgf,ncgf),STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

             cmatrix = 0.0_dp

             ! Transform spherical MOs to Cartesian MOs

             icgf = 1
             isgf = 1
             DO iatom=1,natom
                NULLIFY (orb_basis_set,dftb_parameter)
                CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                     orb_basis_set=orb_basis_set,&
                     dftb_parameter=dftb_parameter)
                IF (ASSOCIATED(orb_basis_set)) THEN
                   CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                        nset=nset,&
                        nshell=nshell,&
                        l=l)
                   DO iset=1,nset
                      DO ishell=1,nshell(iset)
                         lshell = l(ishell,iset)
                         CALL dgemm("T","N",nco(lshell),nsgf,nso(lshell),1.0_dp,&
                              orbtramat(lshell)%s2c,nso(lshell),&
                              smatrix(isgf,1),nsgf,0.0_dp,&
                              cmatrix(icgf,1),ncgf)
                         icgf = icgf + nco(lshell)
                         isgf = isgf + nso(lshell)
                      END DO
                   END DO
                ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                   CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
                   DO ishell=1,lmax+1
                      lshell = ishell-1
                      CALL dgemm("T","N",nco(lshell),nsgf,nso(lshell),1.0_dp,&
                           orbtramat(lshell)%s2c,nso(lshell),&
                           smatrix(isgf,1),nsgf,0.0_dp,&
                           cmatrix(icgf,1),ncgf)
                      icgf = icgf + nco(lshell)
                      isgf = isgf + nso(lshell)
                   END DO
                ELSE
                   CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                        "Unknown basis set type. "//&
                        CPSourceFileRef,&
                        only_ionode=.TRUE.)
                END IF
             END DO ! iatom

          ELSE

             IF (my_last) THEN
                name = "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "SPHERICAL MO EIGENVECTORS"
             ELSE
                WRITE (UNIT=name,FMT="(A,I6)")&
                     "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "SPHERICAL MO EIGENVECTORS AFTER SCF STEP",scf_step
             END IF

          END IF ! p_cart

       ELSE IF (p_occ.OR.p_eval) THEN

          IF (my_last) THEN
             name = "MO EIGENVALUES AND MO OCCUPATION NUMBERS"
          ELSE
             WRITE (UNIT=name,FMT="(A,I6)")&
                  "MO EIGENVALUES AND MO OCCUPATION NUMBERS AFTER "//&
                  "SCF STEP",scf_step
          END IF

       END IF ! p_evec

       CALL compress(name)

       ! Print headline

       IF (PRESENT(spin)) THEN
          WRITE (UNIT=iw,FMT="(/,/,T2,A)") spin//" "//TRIM(name)
       ELSE
          WRITE (UNIT=iw,FMT="(/,/,T2,A)") TRIM(name)
       END IF

       ! Check if only a subset of the MOs has to be printed
       IF ((mo_index_range(1) > 0).AND.&
            (mo_index_range(2) > 0).AND.&
            (mo_index_range(2) >= mo_index_range(1))) THEN
          first_mo = MAX(1,mo_index_range(1))
          last_mo = MIN(mo_set%nmo,mo_index_range(2))
       ELSE
          first_mo = 1
          last_mo = mo_set%nmo
       END IF

       IF (p_evec) THEN

          ! Print full MO information

          DO icol=first_mo,last_mo,ncol

             from = icol
             to = MIN((from+ncol-1),last_mo)

             WRITE (UNIT=iw,FMT=fmtstr1) (jcol,jcol=from,to)
             WRITE (UNIT=iw,FMT=fmtstr2) (mo_set%eigenvalues(jcol),jcol=from,to)
             WRITE (UNIT=iw,FMT="(A)") ""

             WRITE (UNIT=iw,FMT=fmtstr2) (mo_set%occupation_numbers(jcol),jcol=from,to)
             WRITE (UNIT=iw,FMT="(A)") ""

             irow = 1

             DO iatom=1,natom

                IF (iatom /= 1) WRITE (UNIT=iw,FMT="(A)") ""

                NULLIFY(orb_basis_set,dftb_parameter)
                CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                     element_symbol=element_symbol,&
                     orb_basis_set=orb_basis_set,&
                     dftb_parameter=dftb_parameter)

                IF (p_cart) THEN

                   IF (ASSOCIATED(orb_basis_set)) THEN
                      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                           nset=nset,&
                           nshell=nshell,&
                           l=l,&
                           cgf_symbol=bcgf_symbol)

                      icgf = 1
                      DO iset=1,nset
                         DO ishell=1,nshell(iset)
                            lshell = l(ishell,iset)
                            DO ico=1,nco(lshell)
                               WRITE (UNIT=iw,FMT=fmtstr3)&
                                    irow,iatom,ADJUSTR(element_symbol),bcgf_symbol(icgf),&
                                    (cmatrix(irow,jcol),jcol=from,to)
                               icgf = icgf + 1
                               irow = irow + 1
                            END DO
                         END DO
                      END DO
                   ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                      CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
                      icgf = 1
                      DO ishell=1,lmax+1
                         lshell = ishell-1
                         DO ico=1,nco(lshell)
                            symbol = cgf_symbol(1,indco(1:3,icgf))
                            symbol(1:2) = "  "
                            WRITE (UNIT=iw,FMT=fmtstr3)&
                                 irow,iatom,ADJUSTR(element_symbol),symbol,&
                                 (cmatrix(irow,jcol),jcol=from,to)
                            icgf = icgf + 1
                            irow = irow + 1
                         END DO
                      END DO
                   ELSE
                      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                           "Unknown basis set type. "//&
                           CPSourceFileRef,&
                           only_ionode=.TRUE.)
                   END IF

                ELSE

                   IF (ASSOCIATED(orb_basis_set)) THEN
                      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                           nset=nset,&
                           nshell=nshell,&
                           l=l,&
                           sgf_symbol=bsgf_symbol)
                      isgf = 1
                      DO iset=1,nset
                         DO ishell=1,nshell(iset)
                            lshell = l(ishell,iset)
                            DO iso=1,nso(lshell)
                               WRITE (UNIT=iw,FMT=fmtstr3)&
                                    irow,iatom,ADJUSTR(element_symbol),bsgf_symbol(isgf),&
                                    (smatrix(irow,jcol),jcol=from,to)
                               isgf = isgf + 1
                               irow = irow + 1
                            END DO
                         END DO
                      END DO
                   ELSE IF (ASSOCIATED(dftb_parameter)) THEN
                      CALL get_dftb_atom_param(dftb_parameter,lmax=lmax)
                      isgf = 1
                      DO ishell=1,lmax+1
                         lshell = ishell-1
                         DO iso=1,nso(lshell)
                            symbol = sgf_symbol(1,lshell,-lshell+iso-1)
                            symbol(1:2) = "  "
                            WRITE (UNIT=iw,FMT=fmtstr3)&
                                 irow,iatom,ADJUSTR(element_symbol),symbol,&
                                 (smatrix(irow,jcol),jcol=from,to)
                            isgf = isgf + 1
                            irow = irow + 1
                         END DO
                      END DO
                   ELSE
                      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                           "Unknown basis type. "//&
                           CPSourceFileRef,&
                           only_ionode=.TRUE.)
                   END IF

                END IF ! p_cart

             END DO ! iatom

          END DO ! icol

          WRITE (UNIT=iw,FMT="(/)")

          ! Release work storage

          DEALLOCATE (smatrix,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          IF (p_cart) THEN
             DEALLOCATE (cmatrix,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF

       ELSE IF (p_occ.OR.p_eval) THEN

          fmtstr4 = "(T2,I9,2X,F28.  ,1X,F24.  )"
          WRITE (UNIT=fmtstr4(15:16),FMT="(I2)") after
          WRITE (UNIT=fmtstr4(25:26),FMT="(I2)") after
          WRITE (UNIT=iw,FMT="(/,A)")&
               "# MO index          MO eigenvalue [a.u.]            MO occupation"
          DO imo=first_mo,last_mo
             WRITE (UNIT=iw,FMT=fmtstr4)&
                  imo,mo_set%eigenvalues(imo),mo_set%occupation_numbers(imo)
          END DO
          fmtstr5 = "(A,T42,F24.  ,/)"
          WRITE (UNIT=fmtstr5(12:13),FMT="(I2)") after
          WRITE (UNIT=iw,FMT=fmtstr5)&
               "# Sum",accurate_sum(mo_set%occupation_numbers(:))

       END IF ! p_evec

       fmtstr6 = "(A,T17,F24.  ,/)"
       WRITE (UNIT=fmtstr6(12:13),FMT="(I2)") after
       WRITE (UNIT=iw,FMT=fmtstr6) "  Fermi energy:",mo_set%mu
       IF ((mo_set%uniform_occupation).AND.(last_mo > mo_set%homo)) THEN
          gap = mo_set%eigenvalues(mo_set%homo+1) -&
                mo_set%eigenvalues(mo_set%homo)
          fmtstr7 = "(A,T17,F24.  ,A,F6.2,A,/)"
          WRITE (UNIT=fmtstr7(12:13),FMT="(I2)") after
          WRITE (UNIT=iw,FMT=fmtstr7)&
               "  HOMO-LUMO gap:",gap," = ",gap*evolt," eV"
       END IF

    END IF ! iw

    CALL cp_print_key_finished_output(iw,logger,dft_section,"PRINT%MO",&
         ignore_should_output=should_output,&
         error=error)

  END SUBROUTINE write_mo_set_to_output_unit

END MODULE qs_mo_types
