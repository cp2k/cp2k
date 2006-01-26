!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/tbmd_hamiltonian [1.0] *
!!
!!   NAME
!!     tbmd_hamiltonian
!!
!!   FUNCTION
!!     Tight-Binding Hamiltonian built-up
!!
!!   NOTES
!!     Builds the sparse orthogonal tight-binding Hamiltonian matrix. For the 
!!     underlying basis set an orthogonal and minimal sp-basis is assumed.
!!     The linear scaling FIST neighbor-lists of "Professore Krack" are 
!!     heavily used :-)
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     01.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE tbmd_hamiltonian
  USE atomic_kind_types,                      ONLY: atomic_kind_type
  USE cell_types,                             ONLY: cell_type, &
                                                    scaled_to_real
  USE cp_subsystem_types,                     ONLY: cp_subsystem_type
  USE distribution_2d_types,                  ONLY: distribution_2d_type
  USE fist_neighbor_list_types,               ONLY: neighbor_list_set_type, &
                                                    neighbor_list_type, &
                                                    neighbor_node_type, &
                                                    first_list, &
                                                    first_node, &
                                                    next, &
                                                    get_neighbor_list, &
                                                    get_neighbor_list_set, &
                                                    get_neighbor_node
  USE global_types,                           ONLY: global_environment_type
  USE kinds,                                  ONLY: dp
  USE molecule_kind_types,                    ONLY: molecule_kind_type
  USE particle_types,                         ONLY: particle_type
  USE sparse_matrix_types,                    ONLY: real_matrix_type, &
                                                    allocate_matrix, &
                                                    replicate_matrix_structure, &
                                                    add_block_node
  USE tbmd_environment_types,                 ONLY: tbmd_environment_type, &
                                                    tbmd_env_get
  USE tbmd_neighbor_list_env_types,           ONLY: tbmd_neighbor_list_env_type, &
                                                    tbmd_neighbor_list_env_get
  USE timings,                                ONLY: timeset, &
                                                    timestop

#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: moduleN = 'tbmd_hamiltonian'

  ! *** Public subroutines ***
  PUBLIC :: tbmd_build_hamiltonian

!***
!****************************************************************************

CONTAINS

!!****f* tbmd_hamiltonian/tbmd_build_hamiltonian [1.0] *
!!
!!   NAME
!!     tbmd_build_hamiltonian
!!
!!   FUNCTION
!!     Builds up the orthogonal tight-binding Hamiltonian.
!!
!!   REFERENCES
!!     J.C. Slater & G.F. Koster: Simplified LCAO Method for Periodic Potential 
!!                                Problem; Phys. Rev. 94, 1498-1524 (1954)
!!     Mike Finnis: Interatomic Forces in Condensed Matter; Oxford University
!!                  Press (2004)
!!
!!     The orginal reference of Slater & Koster is brackish and not recommendable!
!!     Any textbook should be prefered instead, e.g. the indicated book of Mike Finnis.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - tbmd_env: The tbmd environment to retain
!!     - derivative: Logical which indicates if the dH/dR has to be build-up on-the 
!!                   fly for the force calculation.
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
  SUBROUTINE tbmd_build_hamiltonian(tbmd_env, derivative, error)
    TYPE(tbmd_environment_type), POINTER     :: tbmd_env
    LOGICAL, INTENT(IN)                      :: derivative
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    ! *** Local types *** !
    TYPE block_p_type
      REAL(KIND = dp), DIMENSION(:,:), POINTER :: block
    END TYPE block_p_type

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_build_hamiltonian', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, &
                                                new_atom_b
    INTEGER                                  :: handle, &
                                                ikind, &
                                                jkind, &
                                                nkinds, &
                                                natoms, &
                                                ab, &
                                                ilist, &
                                                nlist, &
                                                inode, &
                                                nnode, &
                                                atom_a, &
                                                atom_b, &
                                                last_atom_b, &
                                                stat, &
                                                irow, &
                                                icol, &
                                                i
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(tbmd_neighbor_list_env_type), &
      POINTER                                :: tbmd_neighbor_list_env
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    REAL(KIND=dp)                            :: rab2_max, &
                                                rab2
    REAL(KIND=dp), DIMENSION(3)              :: dra, &
                                                drb, &
                                                rab
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: r_last_update
    TYPE(cell_type), POINTER                 :: cell
    TYPE(real_matrix_type), POINTER          :: matrix_h, &
                                                matrix_dh
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    INTEGER, DIMENSION(:), ALLOCATABLE       :: first_index_vec, &
                                                last_index_vec
    TYPE(block_p_type), DIMENSION(:), &
      ALLOCATABLE                            :: sint

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(atomic_kind_set, particle_set, subsystem, tbmd_neighbor_list_env, &
            neighbor_list_set, neighbor_list, neighbor_node, r_last_update, &
            cell, matrix_h, matrix_dh, distribution_2d)

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN

      CALL tbmd_env_get(tbmd_env=tbmd_env, subsystem=subsystem, &
                        tbmd_neighbor_list_env=tbmd_neighbor_list_env, &
                        cell=cell, matrix_h=matrix_h, matrix_dh=matrix_dh, &
                        error=error)

      particle_set => subsystem%particles%els
      atomic_kind_set => subsystem%atomic_kinds%els
      nkinds = SIZE(atomic_kind_set)
      natoms = SIZE(particle_set)

      CALL tbmd_neighbor_list_env_get(tbmd_neighbor_list_env=tbmd_neighbor_list_env, &
                                      r_last_update=r_last_update, error=error)

      ALLOCATE(first_index_vec(natoms), STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      
      ALLOCATE(last_index_vec(natoms), STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      CALL allocate_matrix(matrix=tbmd_env%matrix_h, nrow=4*natoms, ncol=4*natoms, &
                           nblock_row=natoms, nblock_col=natoms, &
                           first_row=first_index_vec, last_row=last_index_vec, &
                           first_col=first_index_vec, last_col=last_index_vec, &
                           matrix_name="Tight-Binding Hamiltonian Matrix", &
                           matrix_symmetry="symmetric", &
                           sparsity_id=tbmd_neighbor_list_env%id_nr)
                           !distribution_2d=distribution_2d)

      DEALLOCATE(first_index_vec, STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      DEALLOCATE(last_index_vec, STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      IF (derivative) THEN
        CALL replicate_matrix_structure(source=tbmd_env%matrix_h, TARGET=tbmd_env%matrix_dh, &
                                        target_name="Derivative of the Tight-Binding " //&
                                        "Hamiltonian Matrix WITH respect to the nuclear "//&
                                        "coordinates", target_symmetry="symmetric") ! antisymmetric?
        
        ALLOCATE(sint(1), STAT=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

        NULLIFY(sint(1)%block)
      ELSE
        ALLOCATE(sint(2), STAT=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

        DO i = 1,2
          NULLIFY(sint(i)%block)
        END DO
      END IF

      ! Loop over the neighbor-list
      DO ikind = 1,nkinds
        DO jkind = 1,nkinds
          ab = ikind + nkinds*(jkind - 1)

          neighbor_list_set => tbmd_neighbor_list_env%tbmd_neighbor_list(ab)%neighbor_list_set

          IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE

          rab2_max = tbmd_neighbor_list_env%r_cutoff(ikind,jkind) &
                   * tbmd_neighbor_list_env%r_cutoff(ikind,jkind)

          CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set, &
                                     nlist=nlist)

          neighbor_list => first_list(neighbor_list_set)

          DO ilist=1,nlist

            CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                   atom=atom_a,&
                                   nnode=nnode)

            dra(:) = particle_set(atom_a)%r(:) - &
                     scaled_to_real(r_last_update(:,atom_a), cell)

            neighbor_node => first_node(neighbor_list)

            DO inode=1,nnode

              CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                     neighbor=atom_b,&
                                     r=rab)

              drb(:) = particle_set(atom_b)%r(:) -&
                       scaled_to_real(r_last_update(:,atom_b),cell)

              rab(:) = scaled_to_real(rab(:),cell) - dra(:) + drb(:)

              rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

              IF (rab2 <= rab2_max) THEN

                ! Tight-binding Hamiltonian build-up
                IF (atom_b /= last_atom_b) THEN
                  new_atom_b = .TRUE.
                  last_atom_b = atom_b
                ELSE
                  new_atom_b = .FALSE.
                END IF

                IF (new_atom_b) THEN
                  IF (atom_a <= atom_b) THEN
                    irow = atom_a
                    icol = atom_b
                  ELSE
                    irow = atom_b
                    icol = atom_a
                  END IF
                END IF

                CALL add_block_node(matrix=tbmd_env%matrix_h, block_row=irow, &
                                    block_col=icol, BLOCK=sint(1)%block, error=error)

                IF (derivative) THEN
                  CALL add_block_node(matrix=tbmd_env%matrix_dh, block_row=irow, &
                                      block_col=icol, BLOCK=sint(2)%block, error=error)
                END IF

              END IF

              neighbor_node => next(neighbor_node)

            END DO

            neighbor_list => next(neighbor_list)

          END DO

        END DO

      END DO

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE tbmd_build_hamiltonian
!***************************************************************************

END MODULE tbmd_hamiltonian
