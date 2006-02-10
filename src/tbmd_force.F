!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_force [1.0] *
!!
!!   NAME
!!     tbmd_force
!!
!!   FUNCTION
!!     TBMD force calculation
!!
!!   NOTES
!!     -
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_force
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type, &
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_fm_pool_types,                ONLY: fm_pool_create_fm,&
                                             fm_pool_give_back_fm
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output, &
                                             cp_print_key_unit_nr
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathlib,                         ONLY: symmetrize_matrix
  USE message_passing,                 ONLY: mp_sum, &
                                             mp_sync
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: init_mo_set,&
                                             mo_set_type
  USE qs_scf,                          ONLY: eigensolver,&
                                             simple_eigensolver
  USE sparse_matrix_types,             ONLY: cp_sm_sm_trace, &
                                             replicate_matrix_structure, &
                                             real_matrix_type, &
                                             get_matrix_info, &
                                             copy_local_sm_to_replicated_fm
  USE tbmd_energy,                     ONLY: tbmd_print_energies, &
                                             tbmd_print_hamiltonian, &
                                             tbmd_print_hamiltonian_pattern, &
                                             tbmd_print_eigenvalues, &
                                             tbmd_print_eigenvectors, &
                                             tbmd_print_occupation
  USE tbmd_energy_types,               ONLY: tbmd_energy_clear
  USE tbmd_environment_types,          ONLY: tbmd_env_get,&
                                             tbmd_environment_type
  USE tbmd_force_types,                ONLY: tbmd_force_clear,&
                                             tbmd_force_type
  USE tbmd_hamiltonian,                ONLY: tbmd_build_hamiltonian
  USE tbmd_neighbor_list_env_types,    ONLY: tbmd_neighbor_list_env_type
  USE tbmd_neighbor_lists,             ONLY: tbmd_neighbor_list_control
  USE termination,                     ONLY: stop_program, &
                                             stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tbmd_force'

! *** Public subroutines ***

  PUBLIC :: tbmd_forces, &
            tbmd_print_forces, &
            tbmd_print_deriv_hamiltonian, &
            tbmd_print_deriv_hamiltonian_pattern

!!***
! *****************************************************************************

CONTAINS

!!****f* tbmd_force/tbmd_forces [1.0] *
!!
!!   NAME
!!     tbmd_forces
!!
!!   FUNCTION
!!     Main routine for calculating the TBMD forces yet by using the 
!!     subtle Hellman-Feynman theorem.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - globenv: The global environment
!!     - orthogonal_basis: Logical which indicates in the underlying
!!                         basis-set is assumed to be orthogonal or not.
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_forces(tbmd_env, globenv, orthogonal_basis, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: orthogonal_basis
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, natom, stat, &
                                                iw, &
                                                nparticle_kind, &
                                                iparticle_kind, &
                                                nparticle_local, &
                                                iparticle_local, &
                                                iparticle
    LOGICAL                                  :: failure, ortho_basis, &
                                                use_cholesky, use_jacobi
    REAL(KIND=dp), DIMENSION(:), POINTER     :: trace_vec_ptr
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_type), POINTER                :: fm_matrix_h, fm_matrix_ortho, &
                                                fm_matrix_work
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(mo_set_type), POINTER               :: mos
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_matrix_pools_type), POINTER      :: matrix_pools
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(tbmd_force_type), POINTER           :: tbmd_force
    TYPE(tbmd_neighbor_list_env_type), &
      POINTER                                :: tbmd_neighbor_list_env
    TYPE(cp_logger_type), POINTER            :: logger
    REAL(KIND=dp)                            :: ekin, &
                                                mass
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(tbmd_force, atomic_kind_set, particle_set, local_particles, &
            local_molecules, molecule_kind_set, molecule_set, &
            tbmd_neighbor_list_env, cell, force_env_section, mos, &
            fm_matrix_h, fm_matrix_ortho, fm_matrix_work, matrix_pools, &
            trace_vec_ptr, logger, subsys, atomic_kinds, atomic_kind)

    failure = .FALSE.
    logger => cp_error_get_logger(error)

    IF (PRESENT(orthogonal_basis)) THEN
      ortho_basis = orthogonal_basis
    ELSE
      ortho_basis = .TRUE.
    END IF

    use_cholesky = .TRUE.
    use_jacobi   = .FALSE.

    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN
      force_env_section => section_vals_get_subs_vals(globenv%input_file, "FORCE_EVAL", &
                                                      error=error)

      CALL tbmd_env_get(tbmd_env=tbmd_env, tbmd_force=tbmd_force, &
                        local_particles=local_particles, particle_set=particle_set, &
                        atomic_kind_set=atomic_kind_set, molecule_set=molecule_set, &
                        local_molecules=local_molecules, mos=mos, cell=cell, &
                        molecule_kind_set=molecule_kind_set, &
                        tbmd_neighbor_list_env=tbmd_neighbor_list_env, &
                        matrix_pools=matrix_pools, error=error)

      CALL tbmd_neighbor_list_control(atomic_kind_set=atomic_kind_set, &
                                      particle_set=particle_set, &
                                      local_particles=local_particles, cell=cell, &
                                      tbmd_neighbor_list_env=tbmd_neighbor_list_env, &
                                      globenv=globenv, &
                                      force_env_section=force_env_section, error=error)

      CALL tbmd_force_clear(tbmd_env%tbmd_force, natom=natom, error=error)
      CALL tbmd_energy_clear(tbmd_env%tbmd_energy, error=error)

      CALL tbmd_build_hamiltonian(tbmd_env=tbmd_env, globenv=globenv, derivative=.TRUE., &
                                  error=error)

      natom = SIZE(particle_set)
      IF (ortho_basis) THEN
         CALL init_mo_set(mo_set=mos, fm_pool=matrix_pools%ao_mo_fm_pools(1)%pool, &
                          name="TBMD MO set", error=error)

         CALL fm_pool_create_fm(pool=matrix_pools%ao_ao_fm_pools(1)%pool, &
                                element=fm_matrix_h, name="Full H-Matrix", error=error)

         CALL copy_sm_to_fm(real_matrix=tbmd_env%matrix_h, fm=fm_matrix_h, error=error)
         CALL simple_eigensolver(matrix_ks=fm_matrix_h, mo_set=mos, &
                                 do_level_shift=.FALSE., level_shift=0.0_dp, &
                                 work_syevx=1.0_dp, smear=0.0_dp, error=error)
      ELSE
         CALL init_mo_set(mo_set=mos, fm_pool=matrix_pools%ao_mo_fm_pools(1)%pool, &
                          name="TBMD MO set", error=error)

         CALL fm_pool_create_fm(pool=matrix_pools%ao_ao_fm_pools(1)%pool, &
                                element=fm_matrix_h, name="Full H-Matrix", error=error)

         CALL copy_sm_to_fm(real_matrix=tbmd_env%matrix_h, fm=fm_matrix_h, error=error)
         CALL eigensolver(matrix_ks=fm_matrix_h, mo_set=mos, ortho=fm_matrix_ortho, &
                          work=fm_matrix_work, do_level_shift=.FALSE., &
                          level_shift=0.0_dp, use_cholesky=use_cholesky, &
                          work_syevx=1.0_dp, use_jacobi=use_jacobi, &
                          jacobi_threshold=1.0e-7_dp, smear=0.0_dp, error=error)
      END IF

      CALL replicate_matrix_structure(source=tbmd_env%matrix_h, TARGET=tbmd_env%matrix_p, &
                                      target_name="Density-Matrix", &
                                      target_symmetry="symmetric")

      CALL calculate_density_matrix(mo_set=mos, density_matrix=tbmd_env%matrix_p, &
                                    error=error)

      ! E=Tr[p*H]
      CALL cp_sm_sm_trace(matrix1=tbmd_env%matrix_p, matrix2=tbmd_env%matrix_h, &
                          trace=tbmd_env%tbmd_energy%E_bandstructure, &
                          para_env=tbmd_env%para_env, error=error)

      ! F=Tr[p*(dH/dR)]
      DO i = 1,3
        trace_vec_ptr => tbmd_env%tbmd_force%f_bandstruct(:,1)
        CALL cp_sm_sm_trace(matrix1=tbmd_env%matrix_p, matrix2=tbmd_env%matrix_dh(i)%matrix, &
                            trace_vec=trace_vec_ptr, &
                            para_env=tbmd_env%para_env, error=error)

         tbmd_env%tbmd_force%f_bandstruct(:,i) = (-1)*trace_vec_ptr(:)

      END DO


      !CALL get_part_ke(md_env, tbmd_energy%E_kinetic, int_grp=globalenv%para_env%group)
      CALL tbmd_env_get(tbmd_env=tbmd_env, subsystem=subsys, error=error)
      CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds)

      nparticle_kind = atomic_kinds%n_els

      DO iparticle_kind = 1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local = 1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          ekin = ekin +  0.5_dp * mass * &
                     (particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3))
        END DO
      END DO

      ! sum all contributions to energy over calculated parts on all processors
      CALL mp_sum(ekin, globenv%para_env%group)
      tbmd_env%tbmd_energy%E_kinetic = ekin

      tbmd_env%tbmd_energy%E_tb = tbmd_env%tbmd_energy%E_kinetic &
                                + tbmd_env%tbmd_energy%E_bandstructure &
                                + tbmd_env%tbmd_energy%E_correction &
                                + tbmd_env%tbmd_energy%E_repulsion

      tbmd_env%tbmd_force%f_tb(:,:) = tbmd_env%tbmd_force%f_bandstruct(:,:) &
                                    + tbmd_env%tbmd_force%f_repulsion(:,:)

      DO iatom = 1,natom
        particle_set(iatom)%f(1:3) = tbmd_env%tbmd_force%f_tb(iatom, 1:3)
      END DO

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                force_env_section, "TBMD%PRINT%H-MATRIX", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%H-MATRIX", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_hamiltonian(tbmd_env, globenv=globenv, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%H-MATRIX", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, force_env_section, &
                "TBMD%PRINT%SPARSITY", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%SPARSITY", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_hamiltonian_pattern(tbmd_env, globenv=globenv, output_unit=iw, &
                                            error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%SPARSITY", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                force_env_section, "TBMD%PRINT%H-MATRIX_DERIV", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%H-MATRIX_DERIV", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_deriv_hamiltonian(tbmd_env, globenv=globenv, output_unit=iw, &
                                          error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%H-MATRIX_DERIV", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, force_env_section, &
                "TBMD%PRINT%SPARSITY_DERIV", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%SPARSITY_DERIV", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_deriv_hamiltonian_pattern(tbmd_env, globenv=globenv, output_unit=iw, &
                                            error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%SPARSITY_DERIV", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                force_env_section, "TBMD%PRINT%EIGENVALUES", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%EIGENVALUES", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_eigenvalues(tbmd_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%EIGENVALUES", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                force_env_section, "TBMD%PRINT%OCCUPATION", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%OCCUPATION", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_occupation(tbmd_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%OCCUPATION", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                force_env_section, "TBMD%PRINT%EIGENVECTORS", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%EIGENVECTORS", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_eigenvectors(tbmd_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%EIGENVECTORS", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                force_env_section, "TBMD%PRINT%ENERGIES", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%ENERGIES", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_energies(tbmd_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%ENERGIES", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                force_env_section, "TBMD%PRINT%FORCES", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%FORCES", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_forces(tbmd_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%FORCES", error=error)
      END IF


      CALL fm_pool_give_back_fm(pool=matrix_pools%ao_mo_fm_pools(1)%pool, &
                                element=mos%mo_coeff, error=error)

      DEALLOCATE(mos%eigenvalues, stat=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      DEALLOCATE(mos%occupation_numbers, stat=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      CALL fm_pool_give_back_fm(pool=matrix_pools%ao_ao_fm_pools(1)%pool, &
                                element=fm_matrix_h, error=error)

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_forces
!***************************************************************************

!!****f* tbmd_force/tbmd_print_forces [1.0] *
!!
!!   NAME
!!     tbmd_print_forces
!!
!!   FUNCTION
!!     Print routine for the forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     02.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_print_forces(tbmd_env, output_unit, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_print_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, natom, unit_nr
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(tbmd_force_type), POINTER           :: tbmd_force

!   ------------------------------------------------------------------------

    NULLIFY(particle_set, tbmd_force)

    unit_nr = output_unit

    IF (unit_nr > 0) THEN

      CALL tbmd_env_get(tbmd_env=tbmd_env, tbmd_force=tbmd_force, &
                        particle_set=particle_set, error=error)

      natom = SIZE(particle_set)

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The Tight-Binding forces!"
      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "Bandstructure forces [Hartree/Bohr]"
      DO iatom = 1, natom
        WRITE (unit_nr,*) tbmd_force%f_bandstruct(iatom, 1:3)
      END DO
      WRITE(unit_nr, *) "Pair repulsion forces [Hartree/Bohr]"
      DO iatom = 1, natom
        WRITE (unit_nr,*) tbmd_force%f_repulsion(iatom, 1:3)
      END DO
      WRITE(unit_nr, *) "Total TBMD forces [Hartree/Bohr]"
      DO iatom = 1, natom
        WRITE (unit_nr,*) tbmd_force%f_tb(iatom, 1:3)
      END DO

    END IF

  END SUBROUTINE tbmd_print_forces
!***************************************************************************

!!****f* tbmd_force/tbmd_print_deriv_hamiltonian [1.0] *
!!
!!   NAME
!!     tbmd_print_deriv_hamiltonian
!!
!!   FUNCTION
!!     Print routine for the first derivative of the TBMD Hamiltonian-Matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - globenv: The global environment
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     02.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_print_deriv_hamiltonian(tbmd_env, globenv, output_unit, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_print_deriv_hamiltonian', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, &
                                                group, &
                                                natom, &
                                                size_x, &
                                                size_y, &
                                                i, &
                                                istat
    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    CHARACTER(LEN=40)                        :: matrix_symmetry
    CHARACTER(LEN=80)                        :: matrix_name
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!   ------------------------------------------------------------------------

    IF (output_unit > 0) THEN

      DO k = 1,3
        sparse_matrix => tbmd_env%matrix_dh(k)%matrix
        group = globenv%group

        CALL get_matrix_info(matrix=sparse_matrix, &
                             matrix_name=matrix_name, &
                             matrix_symmetry=matrix_symmetry)

        NULLIFY (matrix)

        CALL copy_local_sm_to_replicated_fm(sparse_matrix, matrix)

        IF (matrix_symmetry == "symmetric") THEN
          CALL symmetrize_matrix(matrix, "upper_to_lower")
        ELSE IF (matrix_symmetry == "antisymmetric") THEN
          CALL symmetrize_matrix(matrix, "anti_upper_to_lower")
        ELSE IF (matrix_symmetry /= "none") THEN
          CALL stop_program(routineN, moduleN, __LINE__, &
                            "Invalid matrix structure", globenv)
        END IF

        CALL mp_sum(matrix, group)
        CALL m_flush(output_unit)

        CALL tbmd_env_get(tbmd_env, particle_set=particle_set, error=error)
        natom = SIZE(particle_set)

!     *** Write headline ***

        WRITE (output_unit,*) ""; WRITE (output_unit,*) ""
        !WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") TRIM(matrix_name)
        SELECT CASE(k)
          CASE(1)
            WRITE (output_unit,*) "Derivative of the H-Matrix w.r.t. the coordinates in x-direction"
          CASE(2)
            WRITE (output_unit,*) "Derivative of the H-Matrix w.r.t. the coordinates in y-direction"
          CASE(3)
            WRITE (output_unit,*) "Derivative of the H-Matrix w.r.t. the coordinates in z-direction"
          CASE DEFAULT
            CALL stop_program(routineN, moduleN, __LINE__, &
                 "Invalid matrix derivative dimension", globenv)
        END SELECT

!     *** Write the variable format strings ***
        WRITE (output_unit,*) ""
        size_x = SIZE(matrix(1,:)); WRITE (output_unit,*) "sizeX", size_x
        size_y = SIZE(matrix(:,1)); WRITE (output_unit,*) "sizeY", size_y
        WRITE (output_unit,*) ""

        DO i = 1,SIZE(matrix(1,:))
          WRITE (UNIT=output_unit, FMT='(32F8.3)') matrix(i,:)
        END DO

        CALL mp_sync(group)
        CALL m_flush(output_unit)

        IF (ASSOCIATED(matrix)) THEN
          DEALLOCATE (matrix, STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                           "matrix")
        END IF

      END DO

    END IF

  END SUBROUTINE tbmd_print_deriv_hamiltonian
!***************************************************************************

!!****f* tbmd_force/tbmd_print_deriv_hamiltonian_pattern [1.0] *
!!
!!   NAME
!!     tbmd_print_deriv_hamiltonian_pattern
!!
!!   FUNCTION
!!     Print routine for the sparsity pattern of the first derivative of 
!!     the TBMD Hamiltonian-Matrix
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - globenv: The global environment
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     02.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE tbmd_print_deriv_hamiltonian_pattern(tbmd_env, globenv, output_unit, &
                                                  error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(global_environment_type), POINTER   :: globenv
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_print_deriv_hamiltonian_pattern', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: k, &
                                                group, &
                                                natom, &
                                                size_x, &
                                                size_y, &
                                                i, &
                                                j, &
                                                istat
    TYPE(real_matrix_type), POINTER          :: sparse_matrix
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: matrix
    CHARACTER(LEN=40)                        :: matrix_symmetry
    CHARACTER(LEN=80)                        :: matrix_name
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    CHARACTER, DIMENSION(:), ALLOCATABLE     :: sparsity_pattern

!   ------------------------------------------------------------------------

    IF (output_unit > 0) THEN

      DO k = 1,3
        sparse_matrix => tbmd_env%matrix_dh(k)%matrix
        group = globenv%group

        CALL get_matrix_info(matrix=sparse_matrix, &
                             matrix_name=matrix_name, &
                             matrix_symmetry=matrix_symmetry)

        NULLIFY (matrix)

        CALL copy_local_sm_to_replicated_fm(sparse_matrix, matrix)

        IF (matrix_symmetry == "symmetric") THEN
          CALL symmetrize_matrix(matrix, "upper_to_lower")
        ELSE IF (matrix_symmetry == "antisymmetric") THEN
          CALL symmetrize_matrix(matrix, "anti_upper_to_lower")
        ELSE IF (matrix_symmetry /= "none") THEN
          CALL stop_program(routineN, moduleN, __LINE__, &
                            "Invalid matrix structure", globenv)
        END IF

        CALL mp_sum(matrix, group)
        CALL m_flush(output_unit)

        CALL tbmd_env_get(tbmd_env, particle_set=particle_set, error=error)
        natom = SIZE(particle_set)

        ALLOCATE(sparsity_pattern(SIZE(matrix(:,1))), STAT=istat)
        CPPostconditionNoFail(istat==0, cp_warning_level, routineP, error)

!     *** Write headline ***

        WRITE (output_unit,*) ""; WRITE (output_unit,*) ""
        !WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") TRIM(matrix_name)
        SELECT CASE(k)
          CASE(1)
            WRITE (output_unit,*) "Sparsity Pattern of the derivative of the H-Matrix"//&
                                  " w.r.t. the coordinates in x-direction"
          CASE(2)
            WRITE (output_unit,*) "Sparsity Pattern of the derivative of the H-Matrix"//&
                                  " w.r.t. the coordinates in y-direction"
          CASE(3)
            WRITE (output_unit,*) "Sparsity Pattern of the derivative of the H-Matrix"//&
                                  " w.r.t. the coordinates in z-direction"
          CASE DEFAULT
            CALL stop_program(routineN, moduleN, __LINE__, &
                 "Invalid matrix derivative dimension", globenv)
        END SELECT

!     *** Write the variable format strings ***
        WRITE (output_unit,*) ""
        size_x = SIZE(matrix(1,:)); WRITE (output_unit,*) "sizeX", size_x
        size_y = SIZE(matrix(:,1)); WRITE (output_unit,*) "sizeY", size_y
        WRITE (output_unit,*) ""

        DO i = 1,SIZE(matrix(1,:))
          DO j = 1,SIZE(matrix(:,1))
            IF (matrix(i,j) == 0.0_dp) THEN
              sparsity_pattern(j) = "0"
            ELSE
              sparsity_pattern(j) = "X"
            END IF
          END DO
          WRITE (UNIT=output_unit, FMT='(32A3)') sparsity_pattern(:)
        END DO

        CALL mp_sync(group)
        CALL m_flush(output_unit)

        IF (ASSOCIATED(matrix)) THEN
          DEALLOCATE (matrix, STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN, moduleN, __LINE__, &
                                           "matrix")
        END IF

        DEALLOCATE(sparsity_pattern, stat=istat)
        CPPostconditionNoFail(istat==0, cp_warning_level, routineP, error)

      END DO

    END IF

  END SUBROUTINE tbmd_print_deriv_hamiltonian_pattern
!***************************************************************************

END MODULE tbmd_force
