!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief types that represent a subsys, i.e. a part of the system
!> \par History
!>      07.2003 created [fawzi]
!>      09.2007 cleaned [tlaino] - University of Zurich
!>      22.11.2010 pack/unpack particle routines added (MK)
!> \author Fawzi Mohamed
! *****************************************************************************
MODULE cp_subsys_types
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_release,&
                                             atomic_kind_list_retain,&
                                             atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE atprop_types,                    ONLY: atprop_release,&
                                             atprop_type
  USE cell_types,                      ONLY: cell_release,&
                                             cell_retain,&
                                             cell_type,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE colvar_types,                    ONLY: colvar_p_release,&
                                             colvar_p_type
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: cp_result_release,&
                                             cp_result_retain,&
                                             cp_result_type
  USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                             distribution_1d_retain,&
                                             distribution_1d_type
  USE kinds,                           ONLY: dp
  USE mol_kind_new_list_types,         ONLY: mol_kind_new_list_release,&
                                             mol_kind_new_list_retain,&
                                             mol_kind_new_list_type
  USE mol_new_list_types,              ONLY: mol_new_list_release,&
                                             mol_new_list_retain,&
                                             mol_new_list_type
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: deallocate_global_constraint,&
                                             global_constraint_type,&
                                             molecule_type
  USE multipole_types,                 ONLY: multipole_type,&
                                             release_multipole_type,&
                                             retain_multipole_type
  USE particle_list_types,             ONLY: particle_list_release,&
                                             particle_list_retain,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE virial_types,                    ONLY: virial_release,&
                                             virial_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_subsys_types'
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

  PUBLIC :: cp_subsys_type,&
            cp_subsys_p_type

  PUBLIC :: cp_subsys_retain,&
            cp_subsys_release,&
            cp_subsys_get,&
            cp_subsys_set,&
            pack_subsys_particles,&
            unpack_subsys_particles

! *****************************************************************************
!> \brief represents a system: atoms, molecules, their pos,vel,...
!> \param atomic_kinds list with all the kinds in the actual subsys
!> \param particles list with the particles of the actual subsys
!> \param local_particles the particles that are local to the actual processor
!> \param molecule_kinds list with the molecule kinds
!> \param local_molecules the molecule structures of the actual subsys
!>        that are local to this processor
!> \param para_env the parallel environment of the actual subsys
!> \param shell_particles list with the shells of the actual subsys if shell-model is used
!> \param core_particles list with the shells of the actual subsys if shell-model is used
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_subsys_type
     INTEGER :: ref_count = 1
     TYPE (atomic_kind_list_type), POINTER       :: atomic_kinds => Null()
     TYPE (particle_list_type), POINTER          :: particles => Null()
     TYPE (particle_list_type), POINTER          :: shell_particles => Null()
     TYPE (particle_list_type), POINTER          :: core_particles => Null()
     TYPE (distribution_1d_type), POINTER        :: local_particles => Null()
     TYPE (cp_para_env_type), POINTER            :: para_env => Null()
     ! New molecules kinds
     TYPE (mol_new_list_type), POINTER           :: molecules_new => Null()
     TYPE (mol_kind_new_list_type), POINTER      :: molecule_kinds_new => Null()
     TYPE (distribution_1d_type), POINTER        :: local_molecules_new => Null()
     ! Definitions of the collective variables
     TYPE (colvar_p_type), DIMENSION(:), POINTER :: colvar_p => Null()
     ! Intermolecular constraints
     TYPE (global_constraint_type), POINTER      :: gci => Null()
     ! Multipoles
     TYPE (multipole_type), POINTER              :: multipoles => Null()
     TYPE (atprop_type), POINTER                 :: atprop => Null()
     TYPE (virial_type), POINTER                 :: virial => Null()
     TYPE (cp_result_type), POINTER              :: results => Null()
     TYPE(cell_type), POINTER                    :: cell => Null()
  END TYPE cp_subsys_type

! *****************************************************************************
!> \brief represent a pointer to a subsys, to be able to create arrays
!>      of pointers
!> \param subsys the pointer to the subsys
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE cp_subsys_p_type
     TYPE(cp_subsys_type), POINTER :: subsys
  END TYPE cp_subsys_p_type

CONTAINS

! *****************************************************************************
!> \brief retains a subsys (see doc/ReferenceCounting.html)
!> \param subsys the subsys to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_retain(subsys, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(subsys%ref_count>0,cp_failure_level,routineP,error)
       subsys%ref_count=subsys%ref_count+1
    END IF
  END SUBROUTINE cp_subsys_retain

! *****************************************************************************
!> \brief releases a subsys (see doc/ReferenceCounting.html)
!> \param subsys the subsys to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_release(subsys, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_release', &
      routineP = moduleN//':'//routineN

    IF (ASSOCIATED(subsys)) THEN
       CPPreconditionNoFail(subsys%ref_count>0,cp_failure_level,routineP,error)
       subsys%ref_count=subsys%ref_count-1
       IF (subsys%ref_count==0) THEN
          CALL atomic_kind_list_release(subsys%atomic_kinds,error=error)
          CALL particle_list_release(subsys%particles, error=error)
          CALL particle_list_release(subsys%shell_particles, error=error)
          CALL particle_list_release(subsys%core_particles, error=error)
          CALL distribution_1d_release(subsys%local_particles, error=error)
          CALL mol_kind_new_list_release(subsys%molecule_kinds_new, error=error)
          CALL mol_new_list_release(subsys%molecules_new, error=error)
          CALL distribution_1d_release(subsys%local_molecules_new,error=error)
          CALL cp_para_env_release(subsys%para_env, error=error)
          CALL release_multipole_type(subsys%multipoles, error)
          CALL colvar_p_release(subsys%colvar_p, error)
          CALL deallocate_global_constraint(subsys%gci, error)
          CALL atprop_release(subsys%atprop,error=error)
          CALL virial_release(subsys%virial,error=error)
          CALL cp_result_release(subsys%results,error=error)
          CALL cell_release(subsys%cell,error=error)
          DEALLOCATE(subsys)
       END IF
    END IF
  END SUBROUTINE cp_subsys_release

! *****************************************************************************
!> \brief sets various propreties of the subsys
!> \param subsys the subsys you want to modify
!> \param atomic_kinds ...
!> \param particles ...
!> \param local_particles ...
!> \param molecules_new ...
!> \param molecule_kinds_new ...
!> \param local_molecules_new ...
!> \param para_env ...
!> \param colvar_p ...
!> \param shell_particles ...
!> \param core_particles ...
!> \param gci ...
!> \param multipoles ...
!> \param results ...
!> \param cell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>      For the description of the other arguments see cp_subsys_type
!>      attributes.
!> \par History
!>      08.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_set(subsys, atomic_kinds, particles, local_particles,&
       molecules_new, molecule_kinds_new, local_molecules_new, para_env,&
       colvar_p, shell_particles, core_particles, gci, multipoles, results, cell, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(atomic_kind_list_type), OPTIONAL, &
      POINTER                                :: atomic_kinds
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: molecules_new
    TYPE(mol_kind_new_list_type), OPTIONAL, &
      POINTER                                :: molecule_kinds_new
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_molecules_new
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(colvar_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: colvar_p
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: shell_particles, &
                                                core_particles
    TYPE(global_constraint_type), OPTIONAL, &
      POINTER                                :: gci
    TYPE(multipole_type), OPTIONAL, POINTER  :: multipoles
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_subsys_set', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    CPPrecondition(subsys%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (PRESENT(multipoles)) THEN
          CALL retain_multipole_type(multipoles, error)
          CALL release_multipole_type(subsys%multipoles, error)
          subsys%multipoles => multipoles
       END IF
       IF (PRESENT(atomic_kinds)) THEN
          CALL atomic_kind_list_retain(atomic_kinds,error=error)
          CALL atomic_kind_list_release(subsys%atomic_kinds, error=error)
          subsys%atomic_kinds => atomic_kinds
       END IF
       IF (PRESENT(particles)) THEN
          CALL particle_list_retain(particles, error=error)
          CALL particle_list_release(subsys%particles, error=error)
          subsys%particles => particles
       END IF
       IF (PRESENT(local_particles)) THEN
          CALL distribution_1d_retain(local_particles,error=error)
          CALL distribution_1d_release(subsys%local_particles,error=error)
          subsys%local_particles => local_particles
       END IF
       IF (PRESENT(local_molecules_new)) THEN
          CALL distribution_1d_retain(local_molecules_new,error=error)
          CALL distribution_1d_release(subsys%local_molecules_new,error=error)
          subsys%local_molecules_new => local_molecules_new
       END IF
       IF (PRESENT(molecule_kinds_new)) THEN
          CALL mol_kind_new_list_retain(molecule_kinds_new, error=error)
          CALL mol_kind_new_list_release(subsys%molecule_kinds_new, error=error)
          subsys%molecule_kinds_new => molecule_kinds_new
       END IF
       IF (PRESENT(molecules_new)) THEN
          CALL mol_new_list_retain(molecules_new, error=error)
          CALL mol_new_list_release(subsys%molecules_new, error=error)
          subsys%molecules_new => molecules_new
       END IF
       IF (PRESENT(para_env)) THEN
          CALL cp_para_env_retain(para_env, error=error)
          CALL cp_para_env_release(subsys%para_env, error=error)
          subsys%para_env => para_env
       END IF
       IF (PRESENT(colvar_p)) THEN
          CPPrecondition(.NOT.ASSOCIATED(subsys%colvar_p),cp_failure_level,routineP,error,failure)
          subsys%colvar_p=>colvar_p
       ENDIF
       IF (PRESENT(shell_particles)) THEN
          IF(ASSOCIATED(shell_particles)) THEN
             CALL particle_list_retain(shell_particles, error=error)
             CALL particle_list_release(subsys%shell_particles, error=error)
             subsys%shell_particles => shell_particles
          END IF
       END IF
       IF (PRESENT(core_particles)) THEN
          IF(ASSOCIATED(core_particles)) THEN
             CALL particle_list_retain(core_particles, error=error)
             CALL particle_list_release(subsys%core_particles, error=error)
             subsys%core_particles => core_particles
          END IF
       END IF
       IF (PRESENT(gci)) THEN
          CPPrecondition(.NOT.ASSOCIATED(subsys%gci),cp_failure_level,routineP,error,failure)
          subsys%gci => gci
       ENDIF
       IF (PRESENT(results)) THEN
          IF(ASSOCIATED(results)) THEN
             CALL cp_result_retain(results, error=error)
             CALL cp_result_release(subsys%results, error=error)
             subsys%results => results
          END IF
       END IF
       IF (PRESENT(cell)) THEN
          IF(ASSOCIATED(cell)) THEN
             CALL cell_retain(cell, error=error)
             CALL cell_release(subsys%cell, error=error)
             subsys%cell => cell
          END IF
       END IF
    END IF
  END SUBROUTINE cp_subsys_set

! *****************************************************************************
!> \brief returns information about various attributes of the given subsys
!> \param subsys the subsys you want info about
!> \param ref_count ...
!> \param atomic_kinds ...
!> \param atomic_kind_set ...
!> \param particles ...
!> \param particle_set ...
!> \param local_particles ...
!> \param molecules_new ...
!> \param molecule_set ...
!> \param molecule_kinds_new ...
!> \param molecule_kind_set ...
!> \param local_molecules_new ...
!> \param para_env ...
!> \param colvar_p ...
!> \param shell_particles ...
!> \param core_particles ...
!> \param gci ...
!> \param multipoles ...
!> \param natom ...
!> \param nparticle ...
!> \param ncore ...
!> \param nshell ...
!> \param nkind ...
!> \param atprop ...
!> \param virial ...
!> \param results ...
!> \param cell ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>      For the description of the other arguments see cp_subsys_type
!>      attributes.
!> \par History
!>      08.2003 created [fawzi]
!>      22.11.2010 (MK)
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE cp_subsys_get(subsys, ref_count, atomic_kinds, atomic_kind_set, &
                           particles, particle_set,&
                           local_particles, molecules_new,molecule_set, molecule_kinds_new,&
                           molecule_kind_set, local_molecules_new, para_env, colvar_p,&
                           shell_particles, core_particles, gci, multipoles,&
                           natom, nparticle, ncore, nshell, nkind, atprop, virial, &
                           results, cell, error)
    TYPE(cp_subsys_type), POINTER            :: subsys
    INTEGER, INTENT(out), OPTIONAL           :: ref_count
    TYPE(atomic_kind_list_type), OPTIONAL, &
      POINTER                                :: atomic_kinds
    TYPE(atomic_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: atomic_kind_set
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    TYPE(particle_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: particle_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_particles
    TYPE(mol_new_list_type), OPTIONAL, &
      POINTER                                :: molecules_new
    TYPE(molecule_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_set
    TYPE(mol_kind_new_list_type), OPTIONAL, &
      POINTER                                :: molecule_kinds_new
    TYPE(molecule_kind_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: molecule_kind_set
    TYPE(distribution_1d_type), OPTIONAL, &
      POINTER                                :: local_molecules_new
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(colvar_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: colvar_p
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: shell_particles, &
                                                core_particles
    TYPE(global_constraint_type), OPTIONAL, &
      POINTER                                :: gci
    TYPE(multipole_type), OPTIONAL, POINTER  :: multipoles
    INTEGER, INTENT(out), OPTIONAL           :: natom, nparticle, ncore, &
                                                nshell, nkind
    TYPE(atprop_type), OPTIONAL, POINTER     :: atprop
    TYPE(virial_type), OPTIONAL, POINTER     :: virial
    TYPE(cp_result_type), OPTIONAL, POINTER  :: results
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_subsys_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_atom, n_core, n_shell
    LOGICAL                                  :: failure

    failure = .FALSE.
    n_atom = 0
    n_core = 0
    n_shell = 0

    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)
    CPPrecondition(subsys%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       IF (PRESENT(ref_count)) ref_count = subsys%ref_count
       IF (PRESENT(atomic_kinds)) atomic_kinds => subsys%atomic_kinds
       IF (PRESENT(atomic_kind_set)) atomic_kind_set => subsys%atomic_kinds%els
       IF (PRESENT(particles)) particles => subsys%particles
       IF (PRESENT(particle_set)) particle_set => subsys%particles%els
       IF (PRESENT(local_particles)) local_particles => subsys%local_particles
       IF (PRESENT(molecules_new)) molecules_new => subsys%molecules_new
       IF (PRESENT(molecule_set)) molecule_set => subsys%molecules_new%els
       IF (PRESENT(molecule_kinds_new)) molecule_kinds_new => subsys%molecule_kinds_new
       IF (PRESENT(molecule_kind_set)) molecule_kind_set => subsys%molecule_kinds_new%els
       IF (PRESENT(local_molecules_new)) local_molecules_new => subsys%local_molecules_new
       IF (PRESENT(para_env)) para_env => subsys%para_env
       IF (PRESENT(colvar_p)) colvar_p => subsys%colvar_p
       IF (PRESENT(shell_particles)) shell_particles => subsys%shell_particles
       IF (PRESENT(core_particles)) core_particles => subsys%core_particles
       IF (PRESENT(gci)) gci => subsys%gci
       IF (PRESENT(multipoles)) multipoles => subsys%multipoles
       IF (PRESENT(virial)) virial => subsys%virial
       IF (PRESENT(atprop)) atprop => subsys%atprop
       IF (PRESENT(results)) results => subsys%results
       IF (PRESENT(cell)) cell => subsys%cell
       IF (PRESENT(nkind)) nkind=SIZE(subsys%atomic_kinds%els)

       IF (PRESENT(natom).OR.PRESENT(nparticle).OR.PRESENT(nshell)) THEN
          ! An atomic particle set should be present in each subsystem at the moment
          CPPrecondition(ASSOCIATED(subsys%particles),cp_failure_level,routineP,error,failure)
          n_atom = subsys%particles%n_els
          ! Check if we have other kinds of particles in this subsystem
          IF (ASSOCIATED(subsys%shell_particles)) THEN
             n_shell = subsys%shell_particles%n_els
             CPPrecondition(ASSOCIATED(subsys%core_particles),cp_failure_level,routineP,error,failure)
             n_core = subsys%core_particles%n_els
             ! The same number of shell and core particles is assumed
             CPPrecondition((n_core == n_shell),cp_failure_level,routineP,error,failure)
          ELSE IF (ASSOCIATED(subsys%core_particles)) THEN
             ! This case should not occur at the moment
             CPPrecondition(ASSOCIATED(subsys%shell_particles),cp_failure_level,routineP,error,failure)
          ELSE
             n_core = 0
             n_shell = 0
          END IF
          IF (PRESENT(natom)) natom = n_atom
          IF (PRESENT(nparticle)) nparticle = n_atom + n_shell
          IF (PRESENT(ncore)) ncore = n_core
          IF (PRESENT(nshell)) nshell = n_shell
       END IF
    END IF

  END SUBROUTINE cp_subsys_get

! *****************************************************************************
!> \brief   Pack components of a subsystem particle sets into a single vector
!> \param subsys ...
!> \param f ...
!> \param r ...
!> \param s ...
!> \param v ...
!> \param fscale ...
!> \param cell ...
!> \param error ...
!> \date    19.11.10
!> \author  Matthias Krack (MK)
!> \version 1.0
!> \note    It is assumed that f, r, s, or v are properly allocated already
! *****************************************************************************
  SUBROUTINE pack_subsys_particles(subsys,f,r,s,v,fscale,cell,error)

    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(OUT), OPTIONAL                  :: f, r, s, v
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: fscale
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pack_subsys_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, j, k, natom, &
                                                nparticle, nsize, shell_index
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: rs
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)

    IF (PRESENT(s)) THEN
      CPPrecondition(PRESENT(cell),cp_failure_level,routineP,error,failure)
    END IF

    NULLIFY (core_particles)
    NULLIFY (particles)
    NULLIFY (shell_particles)

    CALL cp_subsys_get(subsys,&
                       core_particles=core_particles,&
                       natom=natom,&
                       nparticle=nparticle,&
                       particles=particles,&
                       shell_particles=shell_particles,&
                       error=error)

    nsize = 3*nparticle

    ! Pack forces

    IF (PRESENT(f)) THEN
      CPPrecondition((SIZE(f) >= nsize),cp_failure_level,routineP,error,failure)
      j = 0
      DO iatom=1,natom
        shell_index = particles%els(iatom)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            j = j + 1
            f(j) = particles%els(iatom)%f(i)
          END DO
        ELSE
          DO i=1,3
            j = j + 1
            f(j) = core_particles%els(shell_index)%f(i)
          END DO
          k = 3*(natom + shell_index - 1)
          DO i=1,3
            f(k+i) = shell_particles%els(shell_index)%f(i)
          END DO
        END IF
      END DO
      IF (PRESENT(fscale)) f(1:nsize) = fscale*f(1:nsize)
    END IF

    ! Pack coordinates

    IF (PRESENT(r)) THEN
      CPPrecondition((SIZE(r) >= nsize),cp_failure_level,routineP,error,failure)
      j = 0
      DO iatom=1,natom
        shell_index = particles%els(iatom)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            j = j + 1
            r(j) = particles%els(iatom)%r(i)
          END DO
        ELSE
          DO i=1,3
            j = j + 1
            r(j) = core_particles%els(shell_index)%r(i)
          END DO
          k = 3*(natom + shell_index - 1)
          DO i=1,3
            r(k+i) = shell_particles%els(shell_index)%r(i)
          END DO
        END IF
      END DO
    END IF

    ! Pack as scaled coordinates

    IF (PRESENT(s)) THEN
      CPPrecondition((SIZE(s) >= nsize),cp_failure_level,routineP,error,failure)
      CPPrecondition(PRESENT(cell),cp_failure_level,routineP,error,failure)
      j = 0
      DO iatom=1,natom
        shell_index = particles%els(iatom)%shell_index
        IF (shell_index == 0) THEN
          CALL real_to_scaled(rs,particles%els(iatom)%r,cell)
          DO i=1,3
            j = j + 1
            s(j) = rs(i)
          END DO
        ELSE
          CALL real_to_scaled(rs,core_particles%els(shell_index)%r,cell)
          DO i=1,3
            j = j + 1
            s(j) = rs(i)
          END DO
          CALL real_to_scaled(rs,shell_particles%els(shell_index)%r,cell)
          k = 3*(natom + shell_index - 1)
          DO i=1,3
            s(k+i) = rs(i)
          END DO
        END IF
      END DO
    END IF

    ! Pack velocities

    IF (PRESENT(v)) THEN
      CPPrecondition((SIZE(v) >= nsize),cp_failure_level,routineP,error,failure)
      j = 0
      DO iatom=1,natom
        shell_index = particles%els(iatom)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            j = j + 1
            v(j) = particles%els(iatom)%v(i)
          END DO
        ELSE
          DO i=1,3
            j = j + 1
            v(j) = core_particles%els(shell_index)%v(i)
          END DO
          k = 3*(natom + shell_index - 1)
          DO i=1,3
            v(k+i) = shell_particles%els(shell_index)%v(i)
          END DO
        END IF
      END DO
    END IF

  END SUBROUTINE pack_subsys_particles

! *****************************************************************************
!> \brief   Unpack components of a subsystem particle sets into a single vector
!> \param subsys ...
!> \param f ...
!> \param r ...
!> \param s ...
!> \param v ...
!> \param fscale ...
!> \param cell ...
!> \param error ...
!> \date    19.11.10
!> \author  Matthias Krack (MK)
!> \version 1.0
! *****************************************************************************
  SUBROUTINE unpack_subsys_particles(subsys,f,r,s,v,fscale,cell,error)

    TYPE(cp_subsys_type), POINTER            :: subsys
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: f, r, s, v
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: fscale
    TYPE(cell_type), OPTIONAL, POINTER       :: cell
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'unpack_subsys_particles', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, j, k, natom, &
                                                nparticle, nsize, shell_index
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fc, fs, mass, my_fscale
    REAL(KIND=dp), DIMENSION(3)              :: rs
    TYPE(particle_list_type), POINTER        :: core_particles, particles, &
                                                shell_particles

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(subsys),cp_failure_level,routineP,error,failure)

    NULLIFY (core_particles)
    NULLIFY (particles)
    NULLIFY (shell_particles)

    CALL cp_subsys_get(subsys,&
                       core_particles=core_particles,&
                       natom=natom,&
                       nparticle=nparticle,&
                       particles=particles,&
                       shell_particles=shell_particles,&
                       error=error)

    nsize = 3*nparticle

    ! Unpack forces

    IF (PRESENT(f)) THEN
      CPPrecondition((SIZE(f) >= nsize),cp_failure_level,routineP,error,failure)
      IF (PRESENT(fscale)) THEN
        my_fscale = fscale
      ELSE
        my_fscale = 1.0_dp
      END IF
      j = 0
      DO iatom=1,natom
        shell_index = particles%els(iatom)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            j = j + 1
            particles%els(iatom)%f(i) = my_fscale*f(j)
          END DO
        ELSE
          DO i=1,3
            j = j + 1
            core_particles%els(shell_index)%f(i) = my_fscale*f(j)
          END DO
          k = 3*(natom + shell_index - 1)
          DO i=1,3
            shell_particles%els(shell_index)%f(i) = my_fscale*f(k+i)
          END DO
        END IF
      END DO
    END IF

    ! Unpack coordinates

    IF (PRESENT(r)) THEN
      CPPrecondition((SIZE(r) >= nsize),cp_failure_level,routineP,error,failure)
      j = 0
      DO iatom=1,natom
        shell_index = particles%els(iatom)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            j = j + 1
            particles%els(iatom)%r(i) = r(j)
          END DO
        ELSE
          DO i=1,3
            j = j + 1
            core_particles%els(shell_index)%r(i) = r(j)
          END DO
          k = 3*(natom + shell_index - 1)
          DO i=1,3
            shell_particles%els(shell_index)%r(i) = r(k+i)
          END DO
          ! Update atomic position due to core and shell motion
          mass = particles%els(iatom)%atomic_kind%mass
          fc = core_particles%els(shell_index)%atomic_kind%shell%mass_core/mass
          fs = shell_particles%els(shell_index)%atomic_kind%shell%mass_shell/mass
          particles%els(iatom)%r(1:3) = fc*core_particles%els(shell_index)%r(1:3) +&
                                        fs*shell_particles%els(shell_index)%r(1:3)
        END IF
      END DO
    END IF

    ! Unpack scaled coordinates

    IF (PRESENT(s)) THEN
      CPPrecondition((SIZE(s) >= nsize),cp_failure_level,routineP,error,failure)
      CPPrecondition(PRESENT(cell),cp_failure_level,routineP,error,failure)
      j = 0
      DO iatom=1,natom
        shell_index = particles%els(iatom)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            j = j + 1
            rs(i) = s(j)
          END DO
          CALL scaled_to_real(particles%els(iatom)%r,rs,cell)
        ELSE
          DO i=1,3
            j = j + 1
            rs(i) = s(j)
          END DO
          CALL scaled_to_real(core_particles%els(shell_index)%r,rs,cell)
          k = 3*(natom + shell_index - 1)
          DO i=1,3
            rs(i) = s(k+i)
          END DO
          CALL scaled_to_real(shell_particles%els(shell_index)%r,rs,cell)
          ! Update atomic position due to core and shell motion
          mass = particles%els(iatom)%atomic_kind%mass
          fc = core_particles%els(shell_index)%atomic_kind%shell%mass_core/mass
          fs = shell_particles%els(shell_index)%atomic_kind%shell%mass_shell/mass
          particles%els(iatom)%r(1:3) = fc*core_particles%els(shell_index)%r(1:3) +&
                                        fs*shell_particles%els(shell_index)%r(1:3)
        END IF
      END DO
    END IF

    ! Unpack velocities

    IF (PRESENT(v)) THEN
      CPPrecondition((SIZE(v) >= nsize),cp_failure_level,routineP,error,failure)
      j = 0
      DO iatom=1,natom
        shell_index = particles%els(iatom)%shell_index
        IF (shell_index == 0) THEN
          DO i=1,3
            j = j + 1
            particles%els(iatom)%v(i) = v(j)
          END DO
        ELSE
          DO i=1,3
            j = j + 1
            core_particles%els(shell_index)%v(i) = v(j)
          END DO
          k = 3*(natom + shell_index - 1)
          DO i=1,3
            shell_particles%els(shell_index)%v(i) = v(k+i)
          END DO
          ! Update atomic velocity due to core and shell motion
          mass = particles%els(iatom)%atomic_kind%mass
          fc = core_particles%els(shell_index)%atomic_kind%shell%mass_core/mass
          fs = shell_particles%els(shell_index)%atomic_kind%shell%mass_shell/mass
          particles%els(iatom)%v(1:3) = fc*core_particles%els(shell_index)%v(1:3) +&
                                        fs*shell_particles%els(shell_index)%v(1:3)
        END IF
      END DO
    END IF

  END SUBROUTINE unpack_subsys_particles

END MODULE cp_subsys_types
