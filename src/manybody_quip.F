MODULE manybody_quip

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE bibliography,                    ONLY: QUIP_ref,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_set,&
                                             fist_nonbond_env_type,&
                                             quip_data_type
  USE kinds,                           ONLY: dp
  USE pair_potential_types,            ONLY: pair_potential_pp_type,&
                                             pair_potential_single_type,&
                                             quip_pot_type,&
                                             quip_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom,&
                                             evolt
  USE termination,                     ONLY: stop_program
#ifdef __QUIP
  USE quip_unified_wrapper_module,    ONLY : quip_unified_wrapper
#endif

#include "./common/cp_common_uses.f90"

IMPLICIT NONE

PRIVATE

PUBLIC quip_energy_store_force_virial, quip_add_force_virial

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'manybody_quip'

CONTAINS

SUBROUTINE quip_energy_store_force_virial(particle_set, cell, atomic_kind_set, potparm, fist_nonbond_env, &
                                         pot_quip, para_env, error)
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    REAL(kind=dp)                            :: pot_quip
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'quip_energy_store_force_virial', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2), ALLOCATABLE            :: elem_symbol(:)
    INTEGER                                  :: i, iat, iat_use, ikind, &
                                                jkind, n_atoms, n_atoms_use, &
                                                output_unit, stat
    LOGICAL                                  :: do_parallel, failure
    LOGICAL, ALLOCATABLE                     :: use_atom(:)
    REAL(kind=dp)                            :: lattice(3,3), virial(3,3)
    REAL(kind=dp), ALLOCATABLE               :: force(:,:), pos(:,:)
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pair_potential_single_type), &
      POINTER                                :: pot
    TYPE(quip_data_type), POINTER            :: quip_data
    TYPE(quip_pot_type), POINTER             :: quip

    failure = .FALSE.
#ifndef __QUIP
    CALL stop_program(routineN,moduleN,__LINE__,&
        "In order to use QUIP you need to download and install the libAtoms/QUIP library (check CP2K manual for details)")
#else

    n_atoms = SIZE(particle_set)
    ALLOCATE(use_atom(n_atoms))
    use_atom = .FALSE.

    NULLIFY(quip)

    DO ikind=1, SIZE(atomic_kind_set)
    DO jkind=1, SIZE(atomic_kind_set)
      pot => potparm%pot(ikind,jkind)%pot
      DO i=1, SIZE(pot%type)
        IF (pot%type(i) /= quip_type) CYCLE
        IF (.NOT. ASSOCIATED(quip)) quip => pot%set(i)%quip
        DO iat=1, n_atoms
          IF (particle_set(iat)%atomic_kind%kind_number == ikind .OR. &
              particle_set(iat)%atomic_kind%kind_number == jkind) use_atom(iat) = .TRUE.
        END DO ! iat
      END DO ! i
    END DO ! jkind
    END DO ! ikind
    n_atoms_use = COUNT(use_atom)
    ALLOCATE(pos(3,n_atoms_use), force(3,n_atoms_use), elem_symbol(n_atoms_use))

    iat_use = 0
    DO iat=1, n_atoms
      IF (.NOT. use_atom(iat)) CYCLE
      iat_use = iat_use + 1
      pos(1:3, iat_use) = particle_set(iat)%r * angstrom
      elem_symbol(iat_use) = particle_set(iat)%atomic_kind%element_symbol
    END DO
    IF (iat_use > 0) CALL cite_reference(QUIP_ref)
    logger => cp_get_default_logger()
    output_unit = cp_logger_get_default_io_unit(logger)
    lattice = cell%hmat * angstrom
    do_parallel = .FALSE.
    IF (PRESENT(para_env)) THEN
      do_parallel = para_env%num_pe > 1
    END IF
    IF (do_parallel) THEN
      CALL quip_unified_wrapper(&
               N=n_atoms_use, pos=pos, lattice=lattice, symbol=elem_symbol, &
               quip_param_file=TRIM(quip%quip_file_name), &
               quip_param_file_len=LEN_TRIM(quip%quip_file_name), &
               init_args_str=TRIM(quip%init_args), &
               init_args_str_len=LEN_TRIM(quip%init_args), &
               calc_args_str=TRIM(quip%calc_args), &
               calc_args_str_len=LEN_TRIM(quip%calc_args), &
               energy=pot_quip, force=force, virial=virial, &
               output_unit=output_unit, mpi_communicator=para_env%group)
    ELSE
      CALL quip_unified_wrapper(&
               N=n_atoms_use, pos=pos, lattice=lattice, symbol=elem_symbol, &
               quip_param_file=TRIM(quip%quip_file_name), &
               quip_param_file_len=LEN_TRIM(quip%quip_file_name), &
               init_args_str=TRIM(quip%init_args), &
               init_args_str_len=LEN_TRIM(quip%init_args), &
               calc_args_str=TRIM(quip%calc_args), &
               calc_args_str_len=LEN_TRIM(quip%calc_args), &
               energy=pot_quip, force=force, virial=virial, output_unit=output_unit)
    ENDIF
    ! convert units
    pot_quip = pot_quip / evolt
    force = force / (evolt/angstrom)
    virial = virial / evolt
    ! account for double counting from multiple MPI processes
    IF (PRESENT(para_env)) pot_quip = pot_quip / REAL(para_env%num_pe, dp)
    IF (PRESENT(para_env)) force = force / REAL(para_env%num_pe, dp)
    IF (PRESENT(para_env)) virial = virial / REAL(para_env%num_pe, dp)
    ! get quip_data to save force, virial info
    CALL fist_nonbond_env_get ( fist_nonbond_env, quip_data=quip_data, error=error)
    IF (.NOT. ASSOCIATED(quip_data)) THEN
       ALLOCATE(quip_data,STAT=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL fist_nonbond_env_set (fist_nonbond_env,  quip_data=quip_data, error=error)
       NULLIFY(quip_data%use_indices, quip_data%force)
    ENDIF
    IF (ASSOCIATED(quip_data%force)) THEN
      IF (SIZE(quip_data%force,2) /= n_atoms_use) THEN
         DEALLOCATE(quip_data%force, quip_data%use_indices)
         NULLIFY(quip_data%force, quip_data%use_indices)
      END IF
    ENDIF
    IF (.not.ASSOCIATED(quip_data%force)) THEN
       ALLOCATE(quip_data%force(3, n_atoms_use))
       ALLOCATE(quip_data%use_indices(n_atoms_use))
    END IF
    ! save force, virial info
    iat_use = 0
    DO iat=1, n_atoms
      IF (use_atom(iat)) THEN
         iat_use = iat_use + 1
         quip_data%use_indices(iat_use) = iat
      ENDIF
    END DO
    quip_data%force = force
    quip_data%virial = virial

    DEALLOCATE(use_atom, pos, force, elem_symbol)
#endif
END SUBROUTINE quip_energy_store_force_virial

SUBROUTINE quip_add_force_virial(fist_nonbond_env, force, virial, error)
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    REAL(KIND=dp)                            :: force(:,:), virial(3,3)
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'quip_add_force_virial', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iat, iat_use
    LOGICAL                                  :: failure
    TYPE(quip_data_type), POINTER            :: quip_data

    failure = .FALSE.
#ifndef __QUIP
    RETURN
#else

    CALL fist_nonbond_env_get ( fist_nonbond_env, quip_data=quip_data, error=error)
    IF (.NOT. ASSOCIATED(quip_data)) RETURN

    DO iat_use=1, SIZE(quip_data%use_indices)
      iat = quip_data%use_indices(iat_use)
      CPPostcondition(iat >= 1 .AND. iat <= SIZE(force,2),cp_failure_level,routineP,error,failure)
      force(1:3,iat) = force(1:3,iat) + quip_data%force(1:3,iat_use)
    END DO
    virial = virial + quip_data%virial
#endif
END SUBROUTINE quip_add_force_virial

END MODULE manybody_quip
