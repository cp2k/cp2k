!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
!>      - calculate_w_matrix added for the calculation of the Pulay forces
!>        (23.05.02,MK)
!>      - started conversion to LSD (1.2003, Joost VandeVondele)
!> \author Matthias Krack (09.05.2001,MK)
! *****************************************************************************
MODULE qs_mo_types

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                             fm_pool_create_fm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
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
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE scf_control_types,               ONLY: smear_type
  USE scp_environment_types,           ONLY: scp_environment_type
  USE scp_restarts,                    ONLY: read_aux_coeff_set,&
                                             write_scp_coeff_set
  USE string_utilities,                ONLY: compress
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_mo_types'

  TYPE mo_set_type
    ! the actual MO coefficients as a matrix
    TYPE(cp_fm_type), POINTER                 :: mo_coeff
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
    ! highest non-zero occupied orbital
    INTEGER                                   :: homo
    ! lowest non maxocc occupied orbital (e.g. fractional or zero)
    INTEGER                                   :: lfomo
    ! flag that indicates if the MOS have the same occupation number
    LOGICAL                                   :: uniform_occupation
    ! the entropic energy contribution
    REAL(KIND=dp)                             :: kTS
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
            set_mo_occupation_old,&
            set_mo_set, &
            write_mo_set, &
            write_mo_set_low,&
            mo_set_restrict

  INTERFACE read_mo_set
    MODULE PROCEDURE read_mo_set_from_restart
  END INTERFACE

  INTERFACE write_mo_set
    MODULE PROCEDURE write_mo_set_to_output_unit,write_mo_set_to_restart
  END INTERFACE

CONTAINS

! *****************************************************************************
!> \brief allocate a mo set and partially initializes it (nao,nmo,nelectron
!>      are vaild).
!>      For the full initialization you need to call init_mo_set
!> \param mo_set the mo_set to allocate
!> \param nao number of atom orbitals
!> \param nmo number of molecular orbitals
!> \param nelectrons number of electrons
!> \param maxocc maximum occupation of an orbital (LDA: 2, LSD:1)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2002 splitted initialization in two phases [fawzi]
!> \author Matthias Krack
!> \date 15.05.2001
! *****************************************************************************
  SUBROUTINE allocate_mo_set(mo_set,nao,nmo,nelectron,maxocc,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    INTEGER, INTENT(IN)                      :: nao, nmo, nelectron
    REAL(KIND=dp), INTENT(IN)                :: maxocc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_mo_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    failure = .FALSE.

    ALLOCATE (mo_set,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 

    mo_set%maxocc = maxocc
    mo_set%nelectron = nelectron
    mo_set%nao = nao
    mo_set%nmo = nmo
    mo_set%homo = 0
    mo_set%lfomo = 0
    mo_set%uniform_occupation = .TRUE.
    mo_set%kTS = 0.0_dp

    NULLIFY (mo_set%eigenvalues)
    NULLIFY (mo_set%occupation_numbers)
    NULLIFY (mo_set%mo_coeff)

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
    CPPrecondition(.NOT.ASSOCIATED(mo_set%mo_coeff),cp_failure_level,routineP,error,failure)

    ALLOCATE (mo_set%eigenvalues(mo_set%nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
    mo_set%eigenvalues(:) = 0.0_dp

    ALLOCATE (mo_set%occupation_numbers(mo_set%nmo),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
    CALL set_mo_occupation(mo_set=mo_set,error=error)
    CALL fm_pool_create_fm(fm_pool,mo_set%mo_coeff,name=name,error=error)
    CALL cp_fm_get_info(mo_set%mo_coeff,nrow_global=nao,ncol_global=nmo,error=error)
    CPPostcondition((nao >= mo_set%nao),cp_failure_level,routineP,error,failure)
    CPPostcondition((nmo >= mo_set%nmo),cp_failure_level,routineP,error,failure)

  END SUBROUTINE init_mo_set

! *****************************************************************************
!> \brief make the beta orbitals explicitly equal to the alpha orbitals
!>       effectively copying the orbital data
!> \param mo_se t
!> \par History
!>      10.2004 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE mo_set_restrict(mo_array,error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mo_set_restrict', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,"I","",handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(mo_array),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(mo_array).EQ.2,cp_failure_level,routineP,error,failure)
    CPPrecondition(mo_array(1)%mo_set%nmo>=mo_array(2)%mo_set%nmo,cp_failure_level,routineP,error,failure)

    ! first nmo_beta orbitals are copied from alpha to beta
    IF (.NOT.failure) THEN
      CALL cp_fm_to_fm(mo_array(1)%mo_set%mo_coeff,mo_array(2)%mo_set%mo_coeff,mo_array(2)%mo_set%nmo)
    END IF

    CALL timestop(0.0_dp,handle)

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
  SUBROUTINE get_mo_set(mo_set,maxocc,homo,lfomo,nao,nelectron,nmo,&
                        eigenvalues,occupation_numbers,mo_coeff,&
                        uniform_occupation,kTS)

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: maxocc
    INTEGER, INTENT(OUT), OPTIONAL           :: homo, lfomo, nao, nelectron, &
                                                nmo
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eigenvalues, &
                                                occupation_numbers
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: mo_coeff
    LOGICAL, INTENT(OUT), OPTIONAL           :: uniform_occupation
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: kTS

    IF (PRESENT(maxocc)) maxocc = mo_set%maxocc
    IF (PRESENT(homo)) homo = mo_set%homo
    IF (PRESENT(lfomo)) lfomo = mo_set%lfomo
    IF (PRESENT(nao)) nao = mo_set%nao
    IF (PRESENT(nelectron)) nelectron = mo_set%nelectron
    IF (PRESENT(nmo)) nmo = mo_set%nmo
    IF (PRESENT(eigenvalues)) eigenvalues => mo_set%eigenvalues
    IF (PRESENT(occupation_numbers)) THEN
      occupation_numbers => mo_set%occupation_numbers
    END IF
    IF (PRESENT(mo_coeff)) mo_coeff => mo_set%mo_coeff
    IF (PRESENT(uniform_occupation)) uniform_occupation = mo_set%uniform_occupation
    IF (PRESENT(kTS)) kTS = mo_set%kTS

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
  SUBROUTINE set_mo_occupation(mo_set,smear,eval_deriv,error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(smear_type), INTENT(IN), OPTIONAL   :: smear
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eval_deriv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_mo_occupation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imo, nmo, nomo, stat
    LOGICAL                                  :: equal_size, failure, is_large
    REAL(KIND=dp)                            :: e1, e2, edelta, edist, &
                                                lengthscale, mu, nelec
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: dfde

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(mo_set),cp_failure_level,routineP,error,failure)

    nomo = NINT(mo_set%nelectron/mo_set%maxocc)
    nmo  = SIZE(mo_set%eigenvalues)

    ! Check and initialize MO occupations

    CPPrecondition(ASSOCIATED(mo_set%occupation_numbers),cp_failure_level,routineP,error,failure)
    CPPrecondition((SIZE(mo_set%occupation_numbers) >= nmo),cp_failure_level,routineP,error,failure)
    mo_set%occupation_numbers(:) = 0.0_dp
    mo_set%occupation_numbers(1:nomo) =  mo_set%maxocc
    ! zeros don't count as uniform
    mo_set%uniform_occupation = ALL(mo_set%occupation_numbers==mo_set%maxocc)

    mo_set%homo = nomo
    mo_set%lfomo = nomo + 1

    ! check consistency of the array lengths
    IF (PRESENT(eval_deriv)) THEN
        equal_size=SIZE(mo_set%occupation_numbers,1)==SIZE(eval_deriv,1)
        CPPrecondition(equal_size,cp_failure_level,routineP,error,failure)
    ENDIF

    ! Quick return, if no smearing information is supplied (TO BE FIXED, smear should become non-optional...)
    IF (.NOT.PRESENT(smear)) THEN
      ! there is no dependence of the energy on the eigenvalues
      IF (PRESENT(eval_deriv)) THEN
         eval_deriv = 0.0_dp
      ENDIF
      RETURN
    ENDIF

    ! Perform smearing

    IF (smear%do_smear) THEN

      SELECT CASE (smear%method)

      CASE (smear_fermi_dirac)

        Nelec = mo_set%nelectron
        IF (.NOT. PRESENT(eval_deriv)) THEN
          CALL FermiFixed(mo_set%occupation_numbers,mu,mo_set%kTS,mo_set%eigenvalues,Nelec, &
                          smear%electronic_temperature,mo_set%maxocc)
        ELSE
          ! could be a relatively large matrix, but one could get rid of it by never storing it
          ! we only need dE/df * df/de, one could equally parallelize over entries, this could become expensive
          ALLOCATE(dfde(nmo,nmo),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ! lengthscale could become a parameter, but this is pretty good
          lengthscale=10*smear%electronic_temperature

          CALL FermiFixedDeriv(dfde,mo_set%occupation_numbers,mu,mo_set%kTS,mo_set%eigenvalues,Nelec, &
                               smear%electronic_temperature,mo_set%maxocc,lengthscale)

          ! deriv of E_{KS}-kT*S wrt to f_i
          eval_deriv=eval_deriv - mo_set%eigenvalues + mu
          ! correspondingly the deriv of  E_{KS}-kT*S wrt to e_i
          eval_deriv=MATMUL(TRANSPOSE(dfde),eval_deriv)

          DEALLOCATE(dfde,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        ENDIF

        ! Find the lowest fractional occupied MO (LFOMO)
        DO imo=1,nmo
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
        DO imo=nmo,1,-1
          IF (mo_set%occupation_numbers(imo) > 0.0_dp) THEN
            mo_set%homo = imo
            EXIT
          END IF
        END DO
        is_large=ABS(MINVAL(mo_set%occupation_numbers))> smear%eps_fermi_dirac
        CALL cp_assert(.NOT.is_large,cp_warning_level,cp_assertion_failed,routineP,&
                       "Fermi-Dirac smearing includes the last MO => "//&
                       "Add more MOs for proper smearing.",error,failure)

        ! check that the total electron count is accurate
        is_large=(ABS(mo_set%nelectron - accurate_sum(mo_set%occupation_numbers(:))) > smear%eps_fermi_dirac*mo_set%nelectron)
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

        DO imo=1,nomo
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

      CASE(smear_list)
        equal_size=SIZE(mo_set%occupation_numbers,1)==SIZE(smear%list,1)
        CPPrecondition(equal_size,cp_failure_level,routineP,error,failure)
        mo_set%occupation_numbers=smear%list
        ! there is no dependence of the energy on the eigenvalues
        IF (PRESENT(eval_deriv)) THEN
           eval_deriv = 0.0_dp
        ENDIF
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

  END SUBROUTINE set_mo_occupation

! *****************************************************************************
!> \brief   Set the components of a MO set data structure.
!> \author  MK
!> \date    22.04.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_mo_set(mo_set,maxocc,homo,lfomo,nao,nelectron,nmo,&
                        eigenvalues,occupation_numbers,uniform_occupation,&
                        kTS, error)

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: maxocc
    INTEGER, INTENT(IN), OPTIONAL            :: homo, lfomo, nao, nelectron, &
                                                nmo
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eigenvalues, &
                                                occupation_numbers
    LOGICAL, INTENT(IN), OPTIONAL            :: uniform_occupation
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: kTS
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'set_mo_set', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

    IF (PRESENT(maxocc)) mo_set%maxocc = maxocc
    IF (PRESENT(homo)) mo_set%homo = homo
    IF (PRESENT(lfomo)) mo_set%lfomo = lfomo
    IF (PRESENT(nao)) mo_set%nao = nao
    IF (PRESENT(nelectron)) mo_set%nelectron = nelectron
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

  END SUBROUTINE set_mo_set

! *****************************************************************************
  SUBROUTINE write_mo_set_to_restart(mo_array,particle_set,dft_section,scp,&
                                     scp_env,atomic_kind_set,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: dft_section
    LOGICAL, INTENT(IN)                      :: scp
    TYPE(scp_environment_type), OPTIONAL, &
      POINTER                                :: scp_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_to_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=30), DIMENSION(2) :: &
      keys = (/"SCF%PRINT%RESTART_HISTORY","SCF%PRINT%RESTART        "/)
    INTEGER                                  :: handle, ikey, ires
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,"I","",handle)
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
                  extension=".wfn", file_status="REPLACE", file_action="WRITE",&
                  do_backup=.TRUE., file_form="UNFORMATTED", error=error)

             CALL write_mo_set_low(mo_array, particle_set, ires, error)
             IF (scp) CALL write_scp_coeff_set(ires, scp_env, atomic_kind_set, particle_set, error)

             CALL cp_print_key_finished_output(ires,logger,dft_section,TRIM(keys(ikey)), error=error)
          END IF
       END DO
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE write_mo_set_to_restart

! *****************************************************************************
  SUBROUTINE write_mo_set_low(mo_array, particle_set, ires, error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER                                  :: ires
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, i_block, iatom, iset, ishell, ispin, istat, j, &
      lmax, lshell, max_block, nao, natom, nmo, nset, nset_max, nshell_max, &
      nspin
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    CALL timeset(routineN,"I","",handle)
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
    ALLOCATE(vecbuffer(nao,max_block),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 

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

       DO i=1,nmo,MAX(max_block,1)
          i_block=MIN(max_block,nmo-i+1)
          CALL cp_fm_get_submatrix(mo_array(ispin)%mo_set%mo_coeff,vecbuffer,1,i,nao,i_block,error=error)
          ! doing this in one write would increase efficiency, but breaks RESTART compatibility.
          ! to old ones, and in cases where max_block is different between runs, as might happen during
          ! restarts with a different number of CPUs
          DO j=1,i_block
             IF (ires>0) WRITE (ires) vecbuffer(1:nao,j)
          ENDDO
       END DO
    END DO

    DEALLOCATE (vecbuffer,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE write_mo_set_low

! *****************************************************************************
  SUBROUTINE wfn_restart_file_name(filename,exist,section,logger,xas,error)
    CHARACTER(LEN=default_path_length), &
      INTENT(OUT)                            :: filename
    LOGICAL, INTENT(OUT)                     :: exist
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL, INTENT(IN), OPTIONAL            :: xas
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: n_rep_val
    LOGICAL                                  :: my_xas
    TYPE(section_vals_type), POINTER         :: print_key

    my_xas = .FALSE.
    IF(PRESENT(xas)) my_xas = xas

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

  END SUBROUTINE

! *****************************************************************************
  SUBROUTINE read_mo_set_from_restart(mo_array,atomic_kind_set,particle_set,&
       para_env,id_nr,multiplicity,dft_section,scp,scp_env,error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER, INTENT(IN)                      :: id_nr, multiplicity
    TYPE(section_vals_type), POINTER         :: dft_section
    LOGICAL, INTENT(IN)                      :: scp
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mo_set_from_restart', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: file_name
    INTEGER                                  :: handle, ispin, natom, nspin, &
                                                restart_unit
    LOGICAL                                  :: exist
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,"I","",handle)
    logger => cp_error_get_logger(error)

    nspin = SIZE(mo_array)
    restart_unit = -1

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

    CALL read_mos_restart_low (mo_array, para_env, particle_set, natom,&
         restart_unit, multiplicity=multiplicity, error=error)
    IF (scp) CALL read_aux_coeff_set(restart_unit, scp_env%aux_coeff_set, para_env, error)

    ! Close restart file
    IF (para_env%ionode) CALL close_file(unit_number=restart_unit)

    DO ispin = 1,nspin
        CALL write_mo_set(mo_array(ispin)%mo_set,atomic_kind_set,particle_set,&
                            4,6,dft_section,error=error)
    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE read_mo_set_from_restart

! *****************************************************************************
!> \brief Reading the mos from apreviously defined restart file
!> \par History
!>      12.2007 created [MI]
!> \author MI
! *****************************************************************************
  SUBROUTINE read_mos_restart_low (mos, para_env, particle_set, natom, rst_unit, multiplicity, xas, error)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: natom, rst_unit
    INTEGER, INTENT(in), OPTIONAL            :: multiplicity
    LOGICAL, INTENT(IN), OPTIONAL            :: xas
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mos_restart_low', &
      routineP = moduleN//':'//routineN

    INTEGER :: group, homo, homo_read, i, iatom, irow, iset, iset_read, &
      ishell, ishell_read, iso, ispin, istat, lfomo_read, lmax, lshell, &
      my_mult, nao, nao_read, natom_read, nelectron, nelectron_read, nmo, &
      nmo_read, nnshell, nset, nset_max, nshell_max, nspin, nspin_read, &
      offset_read, source
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info, offset_info
    LOGICAL                                  :: failure, minbas, my_xas, &
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
    my_xas = .FALSE.
    IF(PRESENT(xas)) my_xas = xas
    group = para_env%group
    source = para_env%source

    IF (para_env%ionode) THEN
      READ (rst_unit) natom_read,nspin_read,nao_read,nset_max,nshell_max
      IF(my_xas) THEN
        IF (nspin_read /= nspin) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "To change nsin is not possible. "//&
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
          "To change nsin is not possible. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
        ENDIF
      END IF

      IF (natom_read /= natom) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "Incorrect number of atoms in restart file. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
      END IF

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
    END IF

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
        IF(my_xas) THEN
            CALL cp_assert((homo==mos(ispin)%mo_set%homo),cp_warning_level,&
                 cp_assertion_failed,routineP,&
                 "The homo has changed,  you are probably using a different "//&
                 " core-hole potential."//&
CPSourceFileRef,&
                only_ionode=.TRUE.)
        END IF

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

      ! Skip extra MOs if there any
      IF (para_env%ionode) THEN
        !ignore nmo = 0
        IF(nmo>0) THEN 
          DO i=nmo+1,nmo_read
            READ (rst_unit) vecbuffer_read
          END DO
        END IF
      END IF

      IF(.NOT. my_xas) THEN

        IF(ispin==1 .AND. nspin_read<nspin) THEN

          mos(ispin+1)%mo_set%homo = mos(ispin)%mo_set%homo
          mos(ispin+1)%mo_set%lfomo = mos(ispin)%mo_set%lfomo
          nelectron = mos(ispin)%mo_set%nelectron
          IF (MODULO(nelectron + my_mult - 1,2) /= 0) THEN
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "LSD: try to use different multiplicity. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          END IF
          mos(ispin)%mo_set%nelectron = (nelectron + my_mult - 1)/2
          mos(ispin+1)%mo_set%nelectron = (nelectron - my_mult + 1)/2
          IF(mos(ispin+1)%mo_set%nelectron < 0)THEN
            CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
            "LSD: too few electrons for this multiplisity. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          END IF
          mos(ispin+1)%mo_set%eigenvalues = mos(ispin)%mo_set%eigenvalues
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
                                         before,after,dft_section,spin,last,&
                                         error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(section_vals_type), POINTER         :: dft_section
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: spin
    LOGICAL, INTENT(IN), OPTIONAL            :: last
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'write_mo_set_to_output_unit', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=12)                        :: symbol
    CHARACTER(LEN=12), DIMENSION(:), POINTER :: bcgf_symbol
    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=2*default_string_length)   :: name
    CHARACTER(LEN=22)                        :: fmtstr2
    CHARACTER(LEN=25)                        :: fmtstr1
    CHARACTER(LEN=38)                        :: fmtstr3
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: bsgf_symbol
    INTEGER :: from, iatom, icgf, ico, icol, irow, iset, isgf, ishell, iso, &
      istat, iw, jcol, left, lmax, lshell, natom, ncgf, ncol, ncol_global, &
      nrow_global, nset, nsgf, right, scf_step, to, width
    INTEGER, DIMENSION(:), POINTER           :: nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l
    LOGICAL                                  :: failure, ionode, my_last, &
                                                p_cart, p_eval, p_evec, &
                                                p_occ, should_output
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cmatrix, smatrix
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_parameter

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    failure = .FALSE.
    CALL section_vals_val_get(dft_section,"PRINT%MO%EIGENVALUES",l_val=p_eval,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%EIGENVECTORS",l_val=p_evec,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%OCCUPATION_NUMBERS",l_val=p_occ,error=error)
    CALL section_vals_val_get(dft_section,"PRINT%MO%CARTESIAN",l_val=p_cart,error=error)
    should_output = BTEST(cp_print_key_should_output(logger%iter_info,dft_section,&
                          "PRINT%MO",error=error),cp_p_file)

    IF (PRESENT(last)) THEN
      my_last = last
    ELSE
      my_last = .FALSE.
    END IF

    should_output = (should_output.OR.my_last)

    IF ((.NOT.should_output).OR.(.NOT.(p_eval.OR.p_evec.OR.p_occ))) RETURN

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
            natom=natom, ncgf=ncgf, nsgf=nsgf)

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

!            ! Transform spherical MOs to Cartesian MOs

             icgf = 1
             isgf = 1
             DO iatom=1,natom
                NULLIFY(orb_basis_set,dftb_parameter)
                CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                     orb_basis_set=orb_basis_set,dftb_parameter=dftb_parameter)
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
                ELSEIF (ASSOCIATED(dftb_parameter)) THEN
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
             END DO

          ELSE

            IF (my_last) THEN
              name = "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                     "SPHERICAL MO EIGENVECTORS"
            ELSE
              WRITE (UNIT=name,FMT="(A,I6)")&
                "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND "//&
                "SPHERICAL MO EIGENVECTORS AFTER SCF STEP",scf_step
            END IF

          END IF

       ELSE IF (p_occ) THEN

         IF (my_last) THEN
           name = "MO EIGENVALUES AND MO OCCUPATION NUMBERS"
         ELSE
           WRITE (UNIT=name,FMT="(A,I6)")&
             "MO EIGENVALUES AND MO OCCUPATION NUMBERS AFTER "//&
             "SCF STEP",scf_step
         END IF

       ELSE

         IF (my_last) THEN
           name = "MO EIGENVALUES"
         ELSE
           WRITE (UNIT=name,FMT="(A,I6)")&
             "MO EIGENVALUES AFTER SCF STEP",scf_step
         END IF

       END IF

       CALL compress(name)

       ! Print headline

       IF (PRESENT(spin)) THEN
         WRITE (UNIT=iw,FMT="(/,/,T2,A)") spin//" "//TRIM(name)
       ELSE
         WRITE (UNIT=iw,FMT="(/,/,T2,A)") TRIM(name)
       END IF

       ! Print MO information
       DO icol=1,mo_set%nmo,ncol

          from = icol
          to = MIN((from+ncol-1),mo_set%nmo)

          WRITE (UNIT=iw,FMT=fmtstr1)&
               (jcol,jcol=from,to)
          WRITE (UNIT=iw,FMT=fmtstr2)&
               (mo_set%eigenvalues(jcol),jcol=from,to)
          WRITE (UNIT=iw,FMT="(A)") ""

          IF (p_occ) THEN
             WRITE (UNIT=iw,FMT=fmtstr2)&
                  (mo_set%occupation_numbers(jcol),jcol=from,to)
             WRITE (UNIT=iw,FMT="(A)") ""
          END IF

          IF (p_evec) THEN

             irow = 1

             DO iatom=1,natom

                IF (iatom /= 1) WRITE (UNIT=iw,FMT="(A)") ""

                NULLIFY(orb_basis_set,dftb_parameter)
                CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                     element_symbol=element_symbol,&
                     orb_basis_set=orb_basis_set,dftb_parameter=dftb_parameter)

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

                END IF

             END DO

          END IF

       END DO

       WRITE (UNIT=iw,FMT="(/)")

       ! Release work storage

       IF (p_evec) THEN
          DEALLOCATE (smatrix,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
          IF (p_cart) THEN
             DEALLOCATE (cmatrix,STAT=istat)
             CPPostcondition(istat==0,cp_failure_level,routineP,error,failure) 
          END IF
       END IF

    END IF

    CALL cp_print_key_finished_output(iw,logger,dft_section,"PRINT%MO",&
                                      ignore_should_output=should_output,&
                                      error=error)

  END SUBROUTINE write_mo_set_to_output_unit

  !MK Just kept for a while to allow a smooth transition to the new routine
  SUBROUTINE set_mo_occupation_old(occupation,eigenvalues,homo,lfomo,maxocc,nelectron,&
       uocc,smear,error)

    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation, eigenvalues
    INTEGER, INTENT(OUT)                     :: homo, lfomo
    REAL(KIND=dp), INTENT(IN)                :: maxocc
    INTEGER, INTENT(IN)                      :: nelectron
    LOGICAL, INTENT(INOUT)                   :: uocc
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: smear
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_mo_occupation_old', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imo, nmo, nomo, output_unit
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: e1, e2, edelta, edist, &
                                                my_smear, nelec
    TYPE(cp_logger_type), POINTER            :: logger

!   ---------------------------------------------------------------------------

    logger =>cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
       output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
       output_unit = -1
    END IF

    my_smear = 0.0_dp
    IF(PRESENT(smear)) my_smear = smear

    nomo = NINT(nelectron/maxocc)
    nmo  = SIZE(eigenvalues)

    occupation(:) = 0.0_dp
    occupation(1:nomo) =  maxocc

    homo = nomo
    lfomo = nomo + 1

    perform_smear: IF (my_smear /= 0.0_dp) THEN

       !     *** Define the energy window for the eigenvalues ***

       e1 = eigenvalues(nomo) - 0.5_dp*smear
       e2 = eigenvalues(nomo) + 0.5_dp*smear

       !     *** Find the lowest fractional occupied MO ***

       DO imo=1,nomo
          IF (eigenvalues(imo) > e1) THEN
             lfomo = imo
             EXIT
          END IF
       END DO

       !     *** Find the highest fractional occupied (non-zero) MO which will be the HOMO ***
       IF(e2 > eigenvalues(nmo)) THEN
          IF(ionode) THEN
            WRITE(output_unit,*) "WARNING might be too few states for the smearing !"
          END IF
       END IF
       DO imo=nmo,nomo,-1
          IF (eigenvalues(imo) < e2) THEN
             homo = imo
             EXIT
          END IF
       END DO

       !     *** Get the number of electrons to be smeared ***

       edist = 0.0_dp
       nelec = 0.0_dp

       DO imo=lfomo,homo
          nelec = nelec + occupation(imo)
          edist = edist + ABS(e2 - eigenvalues(imo))
       END DO

       !     *** Smear electrons inside the window ***

       DO imo=lfomo,homo
          edelta = ABS(e2 - eigenvalues(imo))
          occupation(imo) = MIN(maxocc,nelec*edelta/edist)
          nelec = nelec - occupation(imo)
          edist = edist - edelta
       END DO

       !     *** Check, if the smearing involves more than one MO ***

       IF (lfomo == homo) THEN
          homo = nomo
          lfomo = nomo + 1
       END IF
       uocc = .FALSE.
    END IF perform_smear

  END SUBROUTINE set_mo_occupation_old

END MODULE qs_mo_types
