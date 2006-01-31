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
  USE atomic_kind_types,                      ONLY: atomic_kind_type
  USE cell_types,                             ONLY: cell_type
  USE cp_fm_pool_types,                       ONLY: fm_pool_create_fm
  USE cp_sm_fm_interactions,                  ONLY: copy_sm_to_fm
  USE cp_fm_types,                            ONLY: cp_fm_type
  USE distribution_1d_types,                  ONLY: distribution_1d_type
  USE global_types,                           ONLY: global_environment_type
  USE input_section_types,                    ONLY: section_vals_type, &
                                                    section_vals_get_subs_vals
  USE kinds,                                  ONLY: dp
  USE molecule_kind_types,                    ONLY: molecule_kind_type
  USE molecule_types_new,                     ONLY: molecule_type
  USE particle_types,                         ONLY: particle_type
  USE qs_matrix_pools,                        ONLY: qs_matrix_pools_type
  USE qs_mo_methods,                          ONLY: calculate_density_matrix
  USE qs_mo_types,                            ONLY: write_mo_set, &
                                                    mo_set_type, &
                                                    init_mo_set
  USE qs_scf,                                 ONLY: eigensolver, &
                                                    simple_eigensolver
  USE sparse_matrix_types,                    ONLY: real_matrix_type
  USE tbmd_energy_types,                      ONLY: tbmd_energy_type
  USE tbmd_environment_types,                 ONLY: tbmd_environment_type, &
                                                    tbmd_env_get
  USE tbmd_hamiltonian,                       ONLY: tbmd_build_hamiltonian
  USE tbmd_neighbor_list_env_types,           ONLY: tbmd_neighbor_list_env_type
  USE tbmd_neighbor_lists,                    ONLY: tbmd_neighbor_list_control
  USE timings,                                ONLY: timeset, &
                                                    timestop

#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: moduleN = "tbmd_energy"

! *** Public subroutines ***

  PUBLIC :: tbmd_energies

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

    LOGICAL                                  :: failure, &
                                                ortho_basis, &
                                                use_cholesky, &
                                                use_jacobi
    INTEGER                                  :: handle
    TYPE(tbmd_energy_type), POINTER          :: tbmd_energy
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles, &
                                                local_molecules
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(tbmd_neighbor_list_env_type), &
      POINTER                                :: tbmd_neighbor_list_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(mo_set_type), POINTER               :: mos
    TYPE(cp_fm_type), POINTER                :: fm_matrix_h, &
                                                fm_matrix_ortho, &
                                                fm_matrix_work
    TYPE(qs_matrix_pools_type), POINTER      :: matrix_pools

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(tbmd_energy, atomic_kind_set, particle_set, local_particles, &
            local_molecules, molecule_kind_set, molecule_set, tbmd_neighbor_list_env, &
            cell, force_env_section, mos, fm_matrix_h, fm_matrix_ortho, &
            fm_matrix_work, matrix_pools)

    failure = .FALSE.

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

      CALL calculate_density_matrix(mo_set=mos, density_matrix=tbmd_env%matrix_p, &
                                    error=error)

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_energies
!***************************************************************************

END MODULE tbmd_energy
