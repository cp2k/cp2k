!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_energy [1.0] *
!!
!!   NAME
!!     tbmd_energy
!!
!!   FUNCTION
!!     TBMD energy calculation
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
MODULE tbmd_energy
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_fm_pool_types,                ONLY: fm_pool_create_fm,&
                                             fm_pool_give_back_fm
  USE cp_fm_types,                     ONLY: cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output, &
                                             cp_print_key_unit_nr
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_types,                  ONLY: particle_type
  USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: init_mo_set,&
                                             mo_set_type
  USE qs_scf,                          ONLY: eigensolver,&
                                             simple_eigensolver
  USE sparse_matrix_types,             ONLY: cp_sm_sm_trace,&
                                             replicate_matrix_structure
  USE tbmd_energy_types,               ONLY: tbmd_energy_clear,&
                                             tbmd_energy_type
  USE tbmd_environment_types,          ONLY: tbmd_env_get,&
                                             tbmd_environment_type
  USE tbmd_hamiltonian,                ONLY: tbmd_build_hamiltonian
  USE tbmd_neighbor_list_env_types,    ONLY: tbmd_neighbor_list_env_type
  USE tbmd_neighbor_lists,             ONLY: tbmd_neighbor_list_control
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tbmd_energy'

! *** Public subroutines ***

  PUBLIC :: tbmd_energies, &
            tbmd_print_energies, &
            tbmd_print_hamiltonian, &
            tbmd_print_eigenvalues, &
            tbmd_print_eigenvectors, &
            tbmd_print_occupation

!!***
! *****************************************************************************

CONTAINS

!!****f* tbmd_energy/tbmd_energies [1.0] *
!!
!!   NAME
!!     tbmd_energies
!!
!!   FUNCTION
!!     Main routine for calculating the TBMD energies.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - globenv: The global environment
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
  SUBROUTINE tbmd_energies(tbmd_env, globenv, orthogonal_basis, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    TYPE(global_environment_type), POINTER   :: globenv
    LOGICAL, INTENT(IN), OPTIONAL            :: orthogonal_basis
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_energies', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat, iw
    LOGICAL                                  :: failure, ortho_basis, &
                                                use_cholesky, use_jacobi
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
    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(tbmd_neighbor_list_env_type), &
      POINTER                                :: tbmd_neighbor_list_env
    TYPE(cp_logger_type), POINTER            :: logger

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(tbmd_energy, atomic_kind_set, particle_set, local_particles, &
            local_molecules, molecule_kind_set, molecule_set, tbmd_neighbor_list_env, &
            cell, force_env_section, mos, fm_matrix_h, fm_matrix_ortho, &
            fm_matrix_work, matrix_pools, logger)

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

      CALL tbmd_env_get(tbmd_env=tbmd_env, tbmd_energy=tbmd_energy, &
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

      CALL tbmd_energy_clear(tbmd_env%tbmd_energy, error=error)

      CALL tbmd_build_hamiltonian(tbmd_env=tbmd_env, derivative=.FALSE., error=error)

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


      tbmd_env%tbmd_energy%E_tb = tbmd_env%tbmd_energy%E_kinetic &
                                + tbmd_env%tbmd_energy%E_bandstructure &
                                + tbmd_env%tbmd_energy%E_correction &
                                + tbmd_env%tbmd_energy%E_repulsion

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                force_env_section, "TBMD%PRINT%H-MATRIX", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger,force_env_section, "TBMD%PRINT%H-MATRIX", &
                                  extension=".mmLog", error=error)

        CALL tbmd_print_hamiltonian(tbmd_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, force_env_section, &
                                          "TBMD%PRINT%H-MATRIX", error=error)
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

  END SUBROUTINE tbmd_energies
!***************************************************************************

!!****f* tbmd_energy/tbmd_print_energies [1.0] *
!!
!!   NAME
!!     tbmd_print_energies
!!
!!   FUNCTION
!!     Print routine for the TBMD energies
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
  SUBROUTINE tbmd_print_energies(tbmd_env, output_unit, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_print_energies', &
      routineP = moduleN//':'//routineN

!   ------------------------------------------------------------------------

    IF (output_unit > 0) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
             "Kinetic energy [Hartree]:        ",tbmd_env%tbmd_energy%E_kinetic,&
             "Bandstructure energy [Hartree]:  ",tbmd_env%tbmd_energy%E_bandstructure,&
             "Pair repulsion energy [Hartree]: ",tbmd_env%tbmd_energy%E_repulsion,&
             "Total TBMD energy [Hartree]:     ",tbmd_env%tbmd_energy%E_tb
    END IF

  END SUBROUTINE tbmd_print_energies
!***************************************************************************

!!****f* tbmd_energy/tbmd_print_hamiltonian [1.0] *
!!
!!   NAME
!!     tbmd_print_hamiltonian
!!
!!   FUNCTION
!!     Print routine for the TBMD Hamiltonian-Matrix
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
  SUBROUTINE tbmd_print_hamiltonian(tbmd_env, output_unit, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_print_hamiltonian', &
      routineP = moduleN//':'//routineN

!   ------------------------------------------------------------------------

    IF (output_unit > 0) THEN

    END IF

  END SUBROUTINE tbmd_print_hamiltonian
!***************************************************************************

!!****f* tbmd_energy/tbmd_print_eigenvalues [1.0] *
!!
!!   NAME
!!     tbmd_print_eigenvalues
!!
!!   FUNCTION
!!     Print routine for the eigenvalues of the TBMD H-Matrix
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
  SUBROUTINE tbmd_print_eigenvalues(tbmd_env, output_unit, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_print_eigenvalues', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: unit_nr, &
                                                i

!   ------------------------------------------------------------------------

    IF (output_unit > 0) THEN

      unit_nr=output_unit

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The eigenvalues of the H-Matrix!"
      WRITE (unit_nr,*) "nmo=", tbmd_env%mos%nmo 
      WRITE (unit_nr,*) "nao=", tbmd_env%mos%nao
      WRITE (unit_nr,*) ""
      DO i = 1,SIZE(tbmd_env%mos%eigenvalues)
        WRITE (unit_nr,*) tbmd_env%mos%eigenvalues(i)
      END DO

    END IF

  END SUBROUTINE tbmd_print_eigenvalues
!***************************************************************************

!!****f* tbmd_energy/tbmd_print_occupation [1.0] *
!!
!!   NAME
!!     tbmd_print_occupation
!!
!!   FUNCTION
!!     Print routine for the occupation numbers
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
  SUBROUTINE tbmd_print_occupation(tbmd_env, output_unit, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_print_occupation', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: unit_nr, &
                                                i

!   ------------------------------------------------------------------------

    IF (output_unit > 0) THEN

      unit_nr=output_unit

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The occupation numbers!"
      WRITE (unit_nr,*) ""
      DO i = 1,SIZE(tbmd_env%mos%occupation_numbers)
        WRITE (unit_nr,*) tbmd_env%mos%occupation_numbers(i)
      END DO

    END IF

  END SUBROUTINE tbmd_print_occupation
!***************************************************************************

!!****f* tbmd_energy/tbmd_print_eigenvectors [1.0] *
!!
!!   NAME
!!     tbmd_print_eigenvectors
!!
!!   FUNCTION
!!     Print routine for the eigenvectors of the TBMD H-Matrix
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
  SUBROUTINE tbmd_print_eigenvectors(tbmd_env, output_unit, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_print_eigenvectors', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr

!   ------------------------------------------------------------------------

    IF (output_unit > 0) THEN

      unit_nr=output_unit

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The eigenvectors of the H-Matrix!"
      WRITE (unit_nr,*) ""
      
    END IF

  END SUBROUTINE tbmd_print_eigenvectors
!***************************************************************************

END MODULE tbmd_energy
