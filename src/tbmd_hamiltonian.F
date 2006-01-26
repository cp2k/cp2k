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
                                                    real_matrix_p_type, &
                                                    allocate_matrix, &
                                                    allocate_matrix_set, &
                                                    replicate_matrix_structure, &
                                                    add_block_node
  USE tbmd_environment_types,                 ONLY: tbmd_environment_type, &
                                                    tbmd_env_get
  USE tbmd_model_types,                       ONLY: tbmd_model_type, &
                                                    tbmd_model_get
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
                                                i, &
                                                i_dim
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
                                                rab2, &
                                                Es, &
                                                dEpsilon, &
                                                ssS, &
                                                spS, &
                                                ppS, &
                                                ppP, &
                                                l, &
                                                m, &
                                                n, &
                                                distance, &
                                                temp_sm_st, &
                                                temp_dsm_st, &
                                                cos_temp, &
                                                dl, &
                                                dm, &
                                                dn
    REAL(KIND=dp), DIMENSION(3)              :: dra, &
                                                drb, &
                                                rab
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: r_last_update
    TYPE(cell_type), POINTER                 :: cell
    TYPE(real_matrix_type), POINTER          :: matrix_h
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_dh
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    INTEGER, DIMENSION(:), ALLOCATABLE       :: first_index_vec, &
                                                last_index_vec
    TYPE(block_p_type), DIMENSION(:), &
      ALLOCATABLE                            :: sint
    TYPE(tbmd_model_type), POINTER           :: tbmd_model

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(atomic_kind_set, particle_set, subsystem, tbmd_neighbor_list_env, &
            neighbor_list_set, neighbor_list, neighbor_node, r_last_update, &
            cell, matrix_h, matrix_dh, distribution_2d, tbmd_model)

    failure = .FALSE.
  
    CPPrecondition(ASSOCIATED(tbmd_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN

      CALL tbmd_env_get(tbmd_env=tbmd_env, subsystem=subsystem, &
                        tbmd_neighbor_list_env=tbmd_neighbor_list_env, &
                        cell=cell, matrix_h=matrix_h, matrix_dh=matrix_dh, &
                        tbmd_model=tbmd_model, error=error)

      particle_set => subsystem%particles%els
      atomic_kind_set => subsystem%atomic_kinds%els
      nkinds = SIZE(atomic_kind_set)
      natoms = SIZE(particle_set)

      CALL tbmd_neighbor_list_env_get(tbmd_neighbor_list_env=tbmd_neighbor_list_env, &
                                      r_last_update=r_last_update, error=error)
      
      CALL tbmd_model_get(tbmd_model=tbmd_model, E_s=Es, dEpsilon=dEpsilon, h_ssS=ssS, &
                          h_spS=spS, h_ppS=ppS, h_ppP=ppP, error=error)


      ALLOCATE(first_index_vec(natoms), STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      
      ALLOCATE(last_index_vec(natoms), STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      CALL tbmd_init_index_vec(first_index_vec=first_index_vec, last_index_vec=last_index_vec, &
                               natoms=natoms, error=error)

      CALL allocate_matrix(matrix=tbmd_env%matrix_h, nrow=4*natoms, ncol=4*natoms, &
                           nblock_row=natoms, nblock_col=natoms, &
                           first_row=first_index_vec, last_row=last_index_vec, &
                           first_col=first_index_vec, last_col=last_index_vec, &
                           matrix_name="Tight-Binding Hamiltonian Matrix", &
                           matrix_symmetry="symmetric", &
                           sparsity_id=tbmd_neighbor_list_env%id_nr, &
                           distribution_2d=distribution_2d)

      DEALLOCATE(first_index_vec, STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      DEALLOCATE(last_index_vec, STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      IF (derivative) THEN
        CALL allocate_matrix_set(tbmd_env%matrix_dh,3)

        DO i_dim = 1,3
          CALL replicate_matrix_structure(source=tbmd_env%matrix_h, &
                                          TARGET=tbmd_env%matrix_dh(i_dim)%matrix, &
                                          target_name="Derivative of the Tight-Binding " //&
                                          "Hamiltonian Matrix WITH respect to the nuclear "//&
                                          "coordinates", target_symmetry="symmetric") ! antisymmetric?
        END DO
        
        ALLOCATE(sint(4), STAT=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

        DO i = 1,4        
          NULLIFY(sint(4)%block)
        END DO
      ELSE
        ALLOCATE(sint(1), STAT=stat)
        CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

        NULLIFY(sint(1)%block)
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

                NULLIFY(sint(1)%block)
                CALL add_block_node(matrix=tbmd_env%matrix_h, block_row=irow, &
                                    block_col=icol, BLOCK=sint(1)%block, error=error)

                ! Cosines
                distance = SQRT(rab2)
                l = rab(1)/distance
                m = rab(2)/distance
                n = rab(3)/distance

                ! Scaling function
                temp_sm_st = 1.0_dp

                ! The 4x4 tight-binding hamiltonian matrix-blocks
                IF (atom_a == atom_b) THEN   ! diagonal-blocks
                  sint(1)%block(:,:) = 0.0_dp

                  sint(1)%block(1,1) = Es
                  sint(1)%block(2,2) = Es + dEpsilon
                  sint(1)%block(3,3) = Es + dEpsilon
                  sint(1)%block(4,4) = Es + dEpsilon
                ELSE                         ! off-diagonal-blocks
                  ! *** Diagonal elements of an off-diagonal-block *** !
                  ! H_ss = V_ssS(r_ij)
                  sint(1)%block(1,1) = temp_sm_st * ssS

                  ! H_xx = (l_ij)^2 * V_ppS(r_ij) + (1 - (l_ij)^2) * V_ppP(r_ij)
                  sint(1)%block(2,2) = temp_sm_st * (ppS * l*l + (1.0_dp - l*l) * ppP)
                  sint(1)%block(3,3) = temp_sm_st * (ppS * m*m + (1.0_dp - m*m) * ppP)
                  sint(1)%block(4,4) = temp_sm_st * (ppS * n*n + (1.0_dp - n*n) * ppP)

                  ! *** s-p orbtial interactions *** !
                  ! H_sx = l_ij * V_spS(r_ij)
                  sint(1)%block(1,2) = temp_sm_st * spS * l
                  sint(1)%block(1,3) = temp_sm_st * spS * m
                  sint(1)%block(1,4) = temp_sm_st * spS * n

                  ! H_xs = -H_sx
                  sint(1)%block(2,1) = -sint(1)%block(1,2)
                  sint(1)%block(3,1) = -sint(1)%block(1,3)
                  sint(1)%block(4,1) = -sint(1)%block(1,4)

                  ! *** p-p orbital interactions *** !
                  ! H_xy = l_ij*m_ij * V_ppS(r_ij) - l_ij*m_ij * V_ppP(r_ij)
                  !      = l_ij*m_ij * (V_ppS(r_ij) - V_ppP(r_ij))
                  sint(1)%block(2,3) = temp_sm_st * l*m * (ppS - ppP)
                  sint(1)%block(2,4) = temp_sm_st * l*n * (ppS - ppP)
                  sint(1)%block(3,4) = temp_sm_st * m*n * (ppS - ppP)

                  ! H_yx = H_xy
                  sint(1)%block(3,2) = sint(1)%block(2,3)
                  sint(1)%block(4,2) = sint(1)%block(2,4)
                  sint(1)%block(4,3) = sint(1)%block(3,4)
                END IF

                ! dH/dR Matrix built-up
                IF (derivative) THEN
                  DO i_dim = 1,3
                    NULLIFY(sint(i_dim+1)%block)
                    CALL add_block_node(matrix=tbmd_env%matrix_dh(i_dim)%matrix, block_row=irow, &
                                        block_col=icol, BLOCK=sint(i_dim+1)%block, error=error)


                    ! Derivatives of the cosines with respect to x_{i_dim}
                    SELECT CASE(i_dim)
                      CASE(1)
                        cos_temp = l
                        dl = 1.0_dp/distance * (1.0_dp - cos_temp*cos_temp)
                        dm = -cos_temp*m/distance
                        dn = -cos_temp*n/distance
                      CASE(2)
                        cos_temp = m
                        dl = -cos_temp*l/distance
                        dm = 1.0_dp/distance * (1.0_dp - cos_temp*cos_temp)
                        dn = -cos_temp*n/distance
                      CASE(3)
                        cos_temp = n
                        dl = -cos_temp*l/distance
                        dm = -cos_temp*m/distance
                        dn = 1.0_dp/distance * (1.0_dp - cos_temp*cos_temp)
                      CASE DEFAULT

                    END SELECT

                    ! Scaling function
                    temp_dsm_st = 1.0_dp

                    ! The 4x4 tight-binding dH/dR_i matrix-blocks
                    IF (atom_a == atom_b) THEN   ! diagonal-blocks
                      sint(i_dim+1)%block(:,:) = 0.0_dp
                    ELSE
                      ! *** Diagonal elements of an off-diagonal-block *** !
                      ! d/dR_i <H_ss> = d/dR_i <V_ssS(r_ij)>
                      sint(i_dim+1)%block(1,1) = temp_dsm_st * cos_temp * ssS

                      ! d/dR_i <H_xx> = d/dR_i <(l_ij)^2 * V_ppS(r_ij) + (1 - (l_ij)^2) * V_ppP(r_ij)>
                      sint(i_dim+1)%block(2,2) = temp_dsm_st * cos_temp &
                                                 * (ppS * l*l + (1.0_dp - l*l) * ppP) &
                                                 + temp_sm_st * 2.0_dp * l * dl * (ppS - ppP)
                      sint(i_dim+1)%block(3,3) = temp_dsm_st * cos_temp &
                                                 * (ppS * m*m + (1.0_dp - m*m) * ppP) &
                                                 + temp_sm_st * 2.0_dp * m * dm * (ppS - ppP)
                      sint(i_dim+1)%block(4,4) = temp_dsm_st * cos_temp &
                                                 * (ppS * n*n + (1.0_dp - n*n) * ppP) &
                                                 + temp_sm_st * 2.0_dp * n * dn * (ppS - ppP)

                      ! *** s-p orbtial interactions *** !
                      ! d/dR_i <H_sx> = d/dR_i <l_ij * V_spS(r_ij)>
                      sint(i_dim+1)%block(1,2) = temp_dsm_st * cos_temp * spS * l &
                                                 + temp_sm_st * spS * dl
                      sint(i_dim+1)%block(1,3) = temp_dsm_st * cos_temp * spS * m &
                                                 + temp_sm_st * spS * dm
                      sint(i_dim+1)%block(1,4) = temp_dsm_st * cos_temp * spS * n &
                                                 + temp_sm_st * spS * dn

                      ! d/dR_i <H_xs> = -d/dR_i <H_sx>
                      sint(i_dim+1)%block(2,1) = -sint(i_dim+1)%block(1,2)
                      sint(i_dim+1)%block(3,1) = -sint(i_dim+1)%block(1,3)
                      sint(i_dim+1)%block(4,1) = -sint(i_dim+1)%block(1,4)

                      ! *** p-p orbital interactions *** !
                      ! d/dR_i <H_xy> = d/dR_i <l_ij*m_ij * V_ppS(r_ij) - l_ij*m_ij * V_ppP(r_ij)>
                      !               = d/dR_i <l_ij*m_ij * (V_ppS(r_ij) - V_ppP(r_ij))>
                      sint(i_dim+1)%block(2,3) = temp_dsm_st * cos_temp * l * m * (ppS - ppP) &
                                                 + temp_sm_st * dl * m * (ppS - ppP) &
                                                 + temp_sm_st * l * dm * (ppS - ppP)
                      sint(i_dim+1)%block(2,4) = temp_dsm_st * cos_temp * l * n * (ppS - ppP) &
                                                 + temp_sm_st * dl * n * (ppS - ppP) &
                                                 + temp_sm_st * l * dn * (ppS - ppP)
                      sint(i_dim+1)%block(3,4) = temp_dsm_st * cos_temp * m * n * (ppS - ppP) &
                                                 + temp_sm_st * dm * n * (ppS - ppP) &
                                                 + temp_sm_st * m * dn * (ppS - ppP)

                      ! d/dR_i <H_yx> = d/dR_i <H_xy>
                      sint(i_dim+1)%block(3,2) = sint(i_dim+1)%block(2,3)
                      sint(i_dim+1)%block(4,2) = sint(i_dim+1)%block(2,4)
                      sint(i_dim+1)%block(4,3) = sint(i_dim+1)%block(3,4)
                    END IF
                  END DO
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

!!****f* tbmd_hamiltonian/tbmd_init_index_vec [1.0] *
!!
!!   NAME
!!     tbmd_init_index_vec
!!
!!   FUNCTION
!!     Correctly initializes the first and last index vectors, nedded for 
!!     the sparse H & dH/dR matrix allocation.
!!
!!   NOTES
!!     In this first orthogonal tight-binding implementation I'm only constructing a 
!!     s-p orbitals Hamiltonian matrix, which I consider as sufficient for my purposes 
!!     in mind. The resulting matrix blocks are therefore of dimensionality 4x4.
!!
!!   INPUTS
!!     - first_index_vec: Vector of the first entry indices of the individual
!!                        blocks of the sparse tight-binding Hamiltonian matrix.
!!     - last_index_vec: Vector consisting of the last indices of the individual
!!                       blocks of the sparse tight-binding Hamiltonian matrix.
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
  SUBROUTINE tbmd_init_index_vec(first_index_vec, last_index_vec, natoms, error)
    INTEGER, DIMENSION(:), INTENT(INOUT)     :: first_index_vec, &
                                                last_index_vec
    INTEGER, INTENT(IN)                      :: natoms
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'tbmd_init_index_vec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i

!   ------------------------------------------------------------------------

   DO i = 1,natoms
     first_index_vec(i) = 4*(i-1) + 1
     last_index_vec(i)  = 4*i
   END DO

  END SUBROUTINE tbmd_init_index_vec
!***************************************************************************

END MODULE tbmd_hamiltonian
