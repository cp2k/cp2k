!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qs_mo_types [1.0] *
!!
!!   NAME
!!     qs_mo_types
!!
!!   FUNCTION
!!     Definition and initialisation of the mo data type.
!!
!!   AUTHOR
!!     Matthias Krack (09.05.2001,MK)
!!
!!   MODIFICATION HISTORY
!!     - adapted to the new QS environment data structure (02.04.2002,MK)
!!     - set_mo_occupation added (17.04.02,MK)
!!     - correct_mo_eigenvalues added (18.04.02,MK)
!!     - calculate_density_matrix moved from qs_scf to here (22.04.02,MK)
!!     - mo_set_p_type added (23.04.02,MK)
!!     - PRIVATE attribute set for TYPE mo_set_type (23.04.02,MK)
!!     - calculate_w_matrix added for the calculation of the Pulay forces
!!       (23.05.02,MK)
!!     - started conversion to LSD (1.2003, Joost VandeVondele)
!!
!!   SOURCE
!******************************************************************************

MODULE qs_mo_types
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                             fm_pool_create_fm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_set_submatrix,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: LOW,&
                                             global_environment_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast
  USE orbital_pointers,                ONLY: nso
  USE output_utilities,                ONLY: print_warning
  USE particle_types,                  ONLY: particle_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: module_name = "qs_mo_types"

  TYPE mo_set_type
    ! the actual MO coefficients as a matrix
    TYPE(cp_fm_type), POINTER :: mo_coeff
    ! number of molecular orbitals (# cols in mo_coeff)
    INTEGER                                   :: nmo       
    ! number of atomic orbitals (# rows in mo_coeff)
    INTEGER                                   :: nao       
    ! occupation - eigenvalues  of the nmo states (if eigenstates)
    REAL(KIND = dp), DIMENSION(:), POINTER    :: eigenvalues,occupation_numbers
    ! maximum allowed occupation number of an MO (1-2)
    REAL(KIND = dp)                           :: maxocc
    ! number of electrons (taking occupation into account)
    INTEGER                                   :: nelectron 
    ! highest non-zero occupied orbital
    INTEGER                                   :: homo 
    ! lowest non maxocc occupied orbital (e.g. fractional or zero)
    INTEGER                                   :: lfomo
  END TYPE mo_set_type

  TYPE mo_set_p_type
    TYPE(mo_set_type), POINTER :: mo_set
  END TYPE mo_set_p_type

! *** Public data types ***

  PUBLIC :: mo_set_p_type,&
            mo_set_type

! *** Public subroutines ***

  PUBLIC :: allocate_mo_set,&
            correct_mo_eigenvalues,&
            deallocate_mo_set,&
            get_mo_set,&
            init_mo_set,&
            read_mo_set,&
            set_mo_occupation,&
            set_mo_set, &
            write_mo_set, &
            mo_set_restrict

! *****************************************************************************

  INTERFACE read_mo_set
    MODULE PROCEDURE read_mo_set_from_restart
  END INTERFACE

  INTERFACE write_mo_set
    MODULE PROCEDURE write_mo_set_to_output_unit,write_mo_set_to_restart
  END INTERFACE

!!***
! *****************************************************************************

CONTAINS

!!****f* qs_mo_types/allocate_mo_set [1.0] *
!!
!!   NAME
!!     allocate_mo_set
!!
!!   SYNOPSIS
!!     Subroutine allocate_mo_set(mo_set, nao, nmo, nelectron, maxocc,&
!!         error)
!!       Type(mo_set_type), Pointer:: mo_set
!!       Real(KIND = dp), Intent (IN):: maxocc
!!       Integer, Intent (IN):: nao, nelectron, nmo
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine allocate_mo_set
!!
!!   FUNCTION
!!     allocate a mo set and partially initializes it (nao,nmo,nelectron
!!     are vaild).
!!     For the full initialization you need to call init_mo_set
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - mo_set: the mo_set to allocate
!!     - nao: number of atom orbitals
!!     - nmo: number of molecular orbitals
!!     - nelectrons: number of electrons
!!     - maxocc: maximum occupation of an orbital (LDA: 2, LSD:1)
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Matthias Krack
!!
!!   MODIFICATION HISTORY
!!     11.2002 splitted initialization in two phases [fawzi]
!!
!!*** *********************************************************************

  SUBROUTINE allocate_mo_set(mo_set,nao,nmo,nelectron,maxocc, error)

!   Purpose: Allocate a wavefunction data structure.

!   History: - Creation (15.05.2001, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set
    INTEGER, INTENT(IN)                      :: nao, nmo, nelectron
    REAL(KIND=dp), INTENT(IN)                :: maxocc
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE allocate_mo_set (MODULE qs_mo_types)"

    INTEGER                                  :: istat
    LOGICAL                                  :: failure

!   ---------------------------------------------------------------------------

    failure=.FALSE.

    ALLOCATE (mo_set,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"mo_set",0)

    mo_set%maxocc = maxocc
    mo_set%nelectron = nelectron
    mo_set%nao = nao
    mo_set%nmo = nmo
    mo_set%homo = 0
    mo_set%lfomo = 0

    NULLIFY (mo_set%eigenvalues)
    NULLIFY (mo_set%occupation_numbers)
    NULLIFY (mo_set%mo_coeff)

  END SUBROUTINE allocate_mo_set

!****************************************************************************

!!****f* qs_mo_types/init_mo_set [1.0] *
!!
!!   NAME
!!     init_mo_set
!!
!!   SYNOPSIS
!!     Subroutine init_mo_set(mo_set, fm_pool, name, error)
!!       Type(mo_set_type), Pointer:: mo_set
!!       Type(cp_fm_pool_type), Pointer:: fm_pool
!!       Character(Len=*), Intent (IN):: name
!!       Type(cp_error_type), Intent (INOUT), Optional:: error
!!     End Subroutine init_mo_set
!!
!!   FUNCTION
!!     initializes an allocated mo_set.
!!     eigenvalues, mo_coeff, occupation_numbers are valid only
!!     after this call.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - mo_set: the mo_set to initialize
!!     - fm_pool: a pool out which you initialize the mo_set
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 rewamped [fawzi]
!!
!!*** *********************************************************************

  SUBROUTINE init_mo_set(mo_set, fm_pool, name, error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_pool_type), POINTER           :: fm_pool
    CHARACTER(len=*), INTENT(in)             :: name
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_mo_set', &
      routineP = module_name//':'//routineN

    INTEGER                                  :: istat, nao, nmo
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(mo_set%mo_coeff),cp_failure_level,routineP,error,failure)

    ALLOCATE (mo_set%eigenvalues(mo_set%nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"mo_set%eigenvalues",mo_set%nmo*dp_size)
    END IF
    mo_set%eigenvalues(:) = 0.0_dp

    ALLOCATE (mo_set%occupation_numbers(mo_set%nmo),STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineP,"mo_set%occupation_numbers",mo_set%nmo*dp_size)
    END IF
    CALL set_mo_occupation(mo_set%occupation_numbers,&
                           mo_set%eigenvalues,mo_set%homo,mo_set%lfomo,&
                           mo_set%maxocc,mo_set%nelectron,smear=0.0_dp)

    CALL fm_pool_create_fm(fm_pool,mo_set%mo_coeff,&
         name=name,error=error)
    CALL cp_fm_get_info(mo_set%mo_coeff, nrow_global=nao, &
         ncol_global=nmo,error=error)
    CPPostcondition(nao>=mo_set%nao,cp_failure_level,routineP,error,failure)
    CPPostcondition(nmo>=mo_set%nmo,cp_failure_level,routineP,error,failure)

  END SUBROUTINE init_mo_set

!!****f* module/mo_set_restrict *
!!
!!   NAME
!!      mo_set_restrict
!!
!!   FUNCTION
!!      make the beta orbitals explicitly equal to the alpha orbitals
!!      effectively copying the orbital data
!!
!!   NOTES
!!
!!   INPUTS
!!    - mo_set
!!
!!   MODIFICATION HISTORY
!!     10.2004 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
  SUBROUTINE mo_set_restrict(mo_array,error)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    LOGICAL :: failure


    CHARACTER(LEN=*), PARAMETER :: routine_name="mo_set_restrict"
    INTEGER :: handle

    CALL timeset(routine_name,"I","",handle)
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(mo_array),cp_failure_level,routine_name,error,failure)
    CPPrecondition(SIZE(mo_array).EQ.2,cp_failure_level,routine_name,error,failure)
    CPPrecondition(mo_array(1)%mo_set%nmo>=mo_array(2)%mo_set%nmo,cp_failure_level,routine_name,error,failure)

    ! first nmo_beta orbitals are copied from alpha to beta
    IF (.NOT. failure) THEN
       CALL cp_fm_to_fm(mo_array(1)%mo_set%mo_coeff,mo_array(2)%mo_set%mo_coeff,mo_array(2)%mo_set%nmo)
    ENDIF

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE mo_set_restrict

! *****************************************************************************


  SUBROUTINE correct_mo_eigenvalues(mo_set,level_shift)

!   Purpose: Correct MO eigenvalues after MO level shifting.

!   History: - Creation (19.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(IN)                :: level_shift

    CHARACTER(LEN=*), PARAMETER :: routine_name = "correct_mo_eigenvalues", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: imo

!   ---------------------------------------------------------------------------

    IF (level_shift == 0.0_dp) RETURN

    DO imo=mo_set%homo+1,mo_set%nmo
      mo_set%eigenvalues(imo) = mo_set%eigenvalues(imo) - level_shift
    END DO

  END SUBROUTINE correct_mo_eigenvalues

! *****************************************************************************

  SUBROUTINE deallocate_mo_set(mo_set)

!   Purpose: Deallocate a wavefunction data structure.

!   History: - Creation (15.05.2001,MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_mo_set", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(mo_set)) THEN
       IF (ASSOCIATED(mo_set%eigenvalues)) THEN
          DEALLOCATE (mo_set%eigenvalues,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"mo_set%eigenvalues")
       END IF

       IF (ASSOCIATED(mo_set%occupation_numbers)) THEN
          DEALLOCATE (mo_set%occupation_numbers,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routine,"mo_set%occupation_numbers")
       END IF

       CALL cp_fm_release(mo_set%mo_coeff)

       DEALLOCATE (mo_set,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routine,"mo_set")
    END IF

  END SUBROUTINE deallocate_mo_set

! *****************************************************************************

  SUBROUTINE get_mo_set(mo_set,maxocc,homo,lfomo,nao,nelectron,nmo,&
                        eigenvalues,occupation_numbers,mo_coeff)

!   Purpose: Get the components of a MO set data structure.

!   History: - Creation (22.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: maxocc
    INTEGER, INTENT(OUT), OPTIONAL           :: homo, lfomo, nao, nelectron, &
                                                nmo
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eigenvalues, &
                                                occupation_numbers
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: mo_coeff

!   ---------------------------------------------------------------------------

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

  END SUBROUTINE get_mo_set

! *****************************************************************************

  SUBROUTINE set_mo_occupation(occupation,eigenvalues,homo,lfomo,maxocc,nelectron,smear)

!   Purpose: Smearing of the MO occupation due to the given smear value.

!   History: - Creation (17.04.2002, MK)

!   ***************************************************************************

    REAL(KIND=dp), DIMENSION(:), POINTER     :: occupation
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues
    INTEGER, INTENT(IN)                      :: nelectron 
    INTEGER, INTENT(OUT)                     :: homo,lfomo
    REAL(KIND=dp), INTENT(IN)                :: maxocc,smear

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_mo_occupation", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: imo, nmo, nomo
    REAL(KIND=dp)                            :: e1, e2, edelta, edist, nelec

!   ---------------------------------------------------------------------------

    nomo = NINT(nelectron/maxocc)
    nmo  = SIZE(eigenvalues)

    occupation(:) = 0.0_dp
    occupation(1:nomo) =  maxocc

    homo = nomo
    lfomo = nomo + 1

    perform_smear: IF (smear /= 0.0_dp) THEN

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
    END IF perform_smear

  END SUBROUTINE set_mo_occupation

! *****************************************************************************

  SUBROUTINE set_mo_set(mo_set,maxocc,homo,lfomo,nao,nelectron,nmo,&
                        eigenvalues,occupation_numbers,mo_coeff)

!   Purpose: Set the components of a MO set data structure.

!   History: - Creation (22.04.2002, MK)

!   ***************************************************************************

    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: maxocc
    INTEGER, INTENT(IN), OPTIONAL            :: homo, lfomo, nao, nelectron, &
                                                nmo
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: eigenvalues, &
                                                occupation_numbers
    TYPE(cp_fm_type), OPTIONAL, POINTER      :: mo_coeff

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_mo_set", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (PRESENT(maxocc)) mo_set%maxocc = maxocc
    IF (PRESENT(homo))  mo_set%homo = homo
    IF (PRESENT(lfomo))  mo_set%lfomo = lfomo
    IF (PRESENT(nao))  mo_set%nao = nao
    IF (PRESENT(nelectron))  mo_set%nelectron = nelectron
    IF (PRESENT(nmo)) mo_set%nmo = nmo
    IF (PRESENT(eigenvalues)) THEN
      IF(ASSOCIATED(mo_set%eigenvalues)) DEALLOCATE(mo_set%eigenvalues,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"mo_set%eigenvalues")
      mo_set%eigenvalues => eigenvalues
    END IF
    IF (PRESENT(occupation_numbers)) THEN
      IF(ASSOCIATED(mo_set%occupation_numbers)) &
                    DEALLOCATE(mo_set%occupation_numbers,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"mo_set%occupation_numbers")
      mo_set%occupation_numbers => occupation_numbers
    END IF

  END SUBROUTINE set_mo_set

! *****************************************************************************
  SUBROUTINE write_mo_set_to_restart(mo_array,atomic_kind_set,particle_set,&
                                     globenv,id_nr)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, INTENT(in)                      :: id_nr

    CHARACTER(LEN=*), PARAMETER :: routine_name = "write_mo_set_to_restart", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    CHARACTER(LEN=200)                       :: file_name
    INTEGER :: handle, i, iatom, iset, ishell, ispin, istat, lshell, nao, &
      natom, nmo, nset, nset_max, nsgf, nshell_max, nspin, restart_unit
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    INTEGER                                  :: max_block,i_block,j

!   *** Local paramters ***
!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    nspin = SIZE(mo_array)

    nao = mo_array(1)%mo_set%nao

    IF (globenv%ionode) THEN

      IF (id_nr==1) THEN
        file_name = globenv%restart_file_name
      ELSE
        file_name = TRIM(globenv%restart_file_name)//"-"//&
                    ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="WRITE",&
                     file_form="UNFORMATTED",&
                     file_status="REPLACE",&
                     unit_number=restart_unit)

!     *** create some info about the basis set first ***

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom=natom,&
                               nsgf=nsgf)

      nset_max = 0
      nshell_max = 0

      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             orb_basis_set=orb_basis_set)
        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                               nset=nset,&
                               nshell=nshell,&
                               l=l)
        nset_max = MAX(nset_max,nset)
        DO iset=1,nset
          nshell_max = MAX(nshell_max,nshell(iset))
        END DO
      END DO

      ALLOCATE (nso_info(nshell_max,nset_max,natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info",&
                                       nshell_max*nset_max*natom*int_size)
      nso_info(:,:,:) = 0

      ALLOCATE (nshell_info(nset_max,natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nshell_info",&
                                       nset_max*natom*int_size)
      nshell_info(:,:) = 0

      ALLOCATE (nset_info(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nset_info",natom*int_size)
      nset_info(:) = 0

      DO iatom=1,natom
        CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                             orb_basis_set=orb_basis_set)
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
      END DO

      WRITE (restart_unit) natom,nspin,nao,nset_max,nshell_max
      WRITE (restart_unit) nset_info
      WRITE (restart_unit) nshell_info
      WRITE (restart_unit) nso_info

      DEALLOCATE (nset_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nset_info")

      DEALLOCATE (nshell_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nshell_info")

      DEALLOCATE (nso_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info")

    END IF

    ! use the scalapack block size as a default for buffering columns
    CALL cp_fm_get_info(mo_array(1)%mo_set%mo_coeff,ncol_block=max_block)
    ALLOCATE(vecbuffer(nao,max_block),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"vecbuffer",nao*max_block*dp_size)

    DO ispin=1,nspin

      nmo=mo_array(ispin)%mo_set%nmo

      IF (globenv%ionode.AND.(nmo > 0)) THEN
        WRITE (restart_unit) nmo,&
                             mo_array(ispin)%mo_set%homo,&
                             mo_array(ispin)%mo_set%lfomo,&
                             mo_array(ispin)%mo_set%nelectron
        WRITE (restart_unit) mo_array(ispin)%mo_set%eigenvalues(1:nmo),&
                             mo_array(ispin)%mo_set%occupation_numbers(1:nmo)
      END IF
      DO i=1,nmo,max_block
        i_block=MIN(max_block,nmo-i+1)
        CALL cp_fm_get_submatrix(mo_array(ispin)%mo_set%mo_coeff,vecbuffer,1,i,nao,i_block)
        ! doing this in one write would increase efficiency, but breaks RESTART compatibility.
        ! to old ones, and in cases where max_block is different between runs, as might happen during
        ! restarts with a different number of CPUs
        DO j=1,i_block
           IF (globenv%ionode) WRITE (restart_unit) vecbuffer(1:nao,j)
        ENDDO
      END DO

    END DO

    DEALLOCATE (vecbuffer,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"vecbuffer")

    IF (globenv%ionode) CALL close_file(unit_number=restart_unit)

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE write_mo_set_to_restart

! *****************************************************************************

  SUBROUTINE read_mo_set_from_restart(mo_array,atomic_kind_set,particle_set,&
                                      globenv,id_nr)

    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mo_array
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    INTEGER, INTENT(in)                      :: id_nr

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_mo_set_from_restart", &
      routine = "SUBROUTINE "//routine_name//" (MODULE "//module_name//")"

    CHARACTER(LEN=200)                       :: file_name
    INTEGER :: group, handle, i, iatom, irow, iset, iset_read, ishell, &
      ishell_read, iso, ispin, istat, lshell, nao, nao_read, natom, &
      natom_read, nmo, nmo_read, nset, nset_max, nsgf, nshell_max, nspin, &
      nspin_read, offset_read, restart_unit, source
    INTEGER, DIMENSION(:), POINTER           :: nset_info, nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l, nshell_info
    INTEGER, DIMENSION(:, :, :), POINTER     :: nso_info, offset_info
    LOGICAL                                  :: did_change, use_this
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: vecbuffer, vecbuffer_read
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

!   *** Local paramters ***
!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I","",handle)

    group = globenv%group
    source = globenv%source

    nspin = SIZE(mo_array)
    nao = mo_array(1)%mo_set%nao

    IF (globenv%ionode) THEN

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               natom=natom,&
                               nsgf=nsgf)
      IF (id_nr==1) THEN
        file_name=globenv%restart_file_name
      ELSE
        file_name = TRIM(globenv%restart_file_name)//"-"//&
                         ADJUSTL(cp_to_string(id_nr))
      END IF

      CALL open_file(file_name=file_name,&
                     file_action="READ",&
                     file_form="UNFORMATTED",&
                     file_status="OLD",&
                     unit_number=restart_unit)

      READ (restart_unit) natom_read,nspin_read,nao_read,nset_max,nshell_max

      ! we should allow for restarting with different spin settings
      IF (nspin_read /= nspin) THEN
        CALL stop_program(routine,"nspin is not equal (program this !!)")
      END IF

      IF (natom_read /= natom) THEN
        CALL stop_program(routine,"incorrect natom")
      END IF

      ALLOCATE (nso_info(nshell_max,nset_max,natom_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info",&
                                       nshell_max*nset_max*natom_read*int_size)

      ALLOCATE (nshell_info(nset_max,natom_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nshell_info",&
                                       nset_max*natom_read*int_size)

      ALLOCATE (nset_info(natom_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nset_info",natom_read*int_size)


      ALLOCATE (offset_info(nshell_max,nset_max,natom_read),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info",&
                                       nshell_max*nset_max*natom_read*int_size)

      IF (nao_read /= nao) THEN
        WRITE(globenv%scr,*) " READ RESTART : WARNING : DIFFERENT # AOs ",nao,nao_read
        did_change=.TRUE.
      END IF

      READ (restart_unit) nset_info
      READ (restart_unit) nshell_info
      READ (restart_unit) nso_info

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
      IF (istat /= 0) CALL stop_memory(routine,"vecbuffer_read",&
                                       nao_read*dp_size)

    END IF

    ALLOCATE (vecbuffer(1,nao),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"vecbuffer",nao*dp_size)

    DO ispin=1,nspin

      nmo=mo_array(ispin)%mo_set%nmo

      mo_array(ispin)%mo_set%eigenvalues(:) = 0.0_dp
      mo_array(ispin)%mo_set%occupation_numbers(:) = 0.0_dp
      CALL cp_fm_set_all(mo_array(ispin)%mo_set%mo_coeff,0.0_dp)

      IF (globenv%ionode.AND.(nmo > 0)) THEN
        READ (restart_unit) nmo_read,&
                            mo_array(ispin)%mo_set%homo,&
                            mo_array(ispin)%mo_set%lfomo,&
                            mo_array(ispin)%mo_set%nelectron
        IF (nmo_read < nmo) THEN
          nmo = nmo_read
          CALL print_warning(routine,"The number of MOs on the restart "//&
                                     "unit is smaller than the number of "//&
                                     "the allocated MOs. The MO set will "//&
                                     "be padded with zeros!")
        ELSE IF (nmo_read > nmo) THEN
          CALL print_warning(routine,"The number of MOs on the restart "//&
                                     "unit is greater than the number of "//&
                                     "the allocated MOs. The read MO set "//&
                                     "will be truncated!")
        END IF
        READ (restart_unit) mo_array(ispin)%mo_set%eigenvalues(1:nmo),&
                            mo_array(ispin)%mo_set%occupation_numbers(1:nmo)
        IF (mo_array(ispin)%mo_set%homo > nmo) THEN
          ! can not make this a warning i.e. homo must be smaller than nmo
          ! otherwise e.g. set_mo_occupation will go out of bounds
          CALL stop_program(routine,"The number of the occupied MOs on the "//&
                                     "restart unit is greater than the "//&
                                     "number of the allocated MOs.")
        END IF
      END IF

      CALL mp_bcast(nmo,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%homo,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%lfomo,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%nelectron,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%eigenvalues,source,group)
      CALL mp_bcast(mo_array(ispin)%mo_set%occupation_numbers,source,group)

      DO i=1,nmo
        IF (globenv%ionode) THEN
          READ (restart_unit) vecbuffer_read
          ! now, try to assign the read to the real vector
          ! in case the basis set changed this involves some guessing
          irow=1
          DO iatom=1,natom
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                                 orb_basis_set=orb_basis_set)
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                   nset=nset,&
                                   nshell=nshell,&
                                   l=l)

            use_this = .TRUE.
            iset_read = 1
            DO iset=1,nset
              ishell_read = 1
              DO ishell=1,nshell(iset)
                lshell = l(ishell,iset)
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
                      ! this indicates a bug
                      ! IF (did_change) WRITE(globenv%scr,*) "-?-",irow,offset_read-1+iso
                      vecbuffer(1,irow)=0.0_dp
                    ELSE
                      vecbuffer(1,irow)=vecbuffer_read(1,offset_read-1+iso)
                      ! IF (did_change) WRITE(globenv%scr,*) "   ",irow,offset_read-1+iso
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
        CALL cp_fm_set_submatrix(mo_array(ispin)%mo_set%mo_coeff,vecbuffer,1,i,nao,1,transpose=.TRUE.)
        did_change=.FALSE.

      END DO

      IF (globenv%print%level > LOW) THEN
        CALL write_mo_set(mo_array(ispin)%mo_set,atomic_kind_set,particle_set,&
                          4,6,globenv)
      END IF

    END DO

    DEALLOCATE(vecbuffer,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"vecbuffer")

    IF (globenv%ionode) THEN
      DEALLOCATE(vecbuffer_read,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"vecbuffer_read")
      DEALLOCATE(offset_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"offset_info")
      DEALLOCATE(nso_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nso_info")
      DEALLOCATE(nshell_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nshell_info")
      DEALLOCATE(nset_info,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"nset_info")
      CALL close_file(unit_number=restart_unit)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE read_mo_set_from_restart

! *****************************************************************************

  SUBROUTINE write_mo_set_to_output_unit(mo_set,atomic_kind_set,particle_set,&
                                         before,after,globenv)

!   Purpose: Write the MO eigenvalues, MO occupation numbers and
!            MO mo_coeff.

!   History: - Creation (15.05.2001,MK)

!   ***************************************************************************

!   List of variables:

!   after : Number of digits after point.
!   before: Number of digits before point.

!   ---------------------------------------------------------------------------

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    INTEGER, INTENT(IN)                      :: before, after
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE write_mo_set_to_output_unit (MODULE qs_mo_types)"

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=22)                        :: fmtstr2
    CHARACTER(LEN=25)                        :: fmtstr1
    CHARACTER(LEN=35)                        :: fmtstr3
    CHARACTER(LEN=6), DIMENSION(:), POINTER  :: sgf_symbol
    CHARACTER(LEN=60)                        :: name
    INTEGER :: from, iatom, icol, irow, iset, isgf, ishell, iso, istat, jcol, &
      left, lshell, natom, ncol, ncol_global, nrow_global, nset, nsgf, &
      output_unit, right, to, width
    INTEGER, DIMENSION(:), POINTER           :: nshell
    INTEGER, DIMENSION(:, :), POINTER        :: l
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set

!   ---------------------------------------------------------------------------

    IF (.NOT.globenv%print%mo_eigenvalues) RETURN

    output_unit = globenv%scr

    IF (globenv%print%mo_eigenvectors) THEN
      CALL cp_fm_get_info(mo_set%mo_coeff,nrow_global=nrow_global,ncol_global=ncol_global)
      ALLOCATE(matrix(nrow_global,ncol_global))
      CALL cp_fm_get_submatrix(mo_set%mo_coeff,matrix)
      IF (.NOT.globenv%ionode) DEALLOCATE(matrix)
    END IF

    IF (.NOT.globenv%ionode) RETURN

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom=natom,&
                             nsgf=nsgf)

!   *** Definition of the variable formats ***

    fmtstr1 = "(/,T2,23X,  (  X,I5,  X))"
    fmtstr2 = "(T2,23X,  (1X,F  .  ))"
    fmtstr3 = "(T2,2I5,2X,A2,1X,A8,   (1X,F  .  ))"

!   *** Write headline ***

    IF (globenv%print%mo_eigenvectors) THEN
      name = "MO EIGENVALUES, MO OCCUPATION NUMBERS AND MO EIGENVECTORS"
    ELSE IF (globenv%print%mo_occupation_numbers) THEN
      name = "MO EIGENVALUES AND MO OCCUPATION NUMBERS"
    ELSE
      name = "MO EIGENVALUES"
    END IF

    WRITE (UNIT=output_unit,FMT="(/,/,T2,A)") TRIM(name)

!   *** Write the variable format strings ***

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

    WRITE (UNIT=fmtstr3(22:23),FMT="(I2)") ncol
    WRITE (UNIT=fmtstr3(29:30),FMT="(I2)") width - 1
    WRITE (UNIT=fmtstr3(32:33),FMT="(I2)") after

!   *** Write the matrix in the selected format ***

    DO icol=1,nsgf,ncol

      from = icol
      to = MIN((from+ncol-1),nsgf)

      WRITE (UNIT=output_unit,FMT=fmtstr1)&
        (jcol,jcol=from,to)
      WRITE (UNIT=output_unit,FMT=fmtstr2)&
        (mo_set%eigenvalues(jcol),jcol=from,to)
      WRITE (UNIT=output_unit,FMT="(A)") ""

      IF (globenv%print%mo_occupation_numbers) THEN
        WRITE (UNIT=output_unit,FMT=fmtstr2)&
          (mo_set%occupation_numbers(jcol),jcol=from,to)
        WRITE (UNIT=output_unit,FMT="(A)") ""
      END IF

      IF (globenv%print%mo_eigenvectors) THEN

        irow = 1

        DO iatom=1,natom

          IF (iatom /= 1) WRITE (UNIT=output_unit,FMT="(A)") ""

          CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind,&
                               element_symbol=element_symbol,&
                               orb_basis_set=orb_basis_set)
          CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                 nset=nset,&
                                 nshell=nshell,&
                                 l=l,&
                                 sgf_symbol=sgf_symbol)

          isgf = 1

          DO iset=1,nset
            DO ishell=1,nshell(iset)
              lshell = l(ishell,iset)
              DO iso=1,nso(lshell)
                WRITE (UNIT=output_unit,FMT=fmtstr3)&
                  irow,iatom,element_symbol,sgf_symbol(isgf),&
                  (matrix(irow,jcol),jcol=from,to)
                isgf = isgf + 1
                irow = irow + 1
              END DO
            END DO
          END DO

        END DO

      END IF

    END DO

    WRITE (UNIT=output_unit,FMT="(/)")

!   *** Release work storage ***

    IF (globenv%print%mo_eigenvectors) THEN
      DEALLOCATE (matrix,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"matrix")
    END IF

  END SUBROUTINE write_mo_set_to_output_unit

! *****************************************************************************

END MODULE qs_mo_types
