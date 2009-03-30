!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for data exchange between MPI processes 
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon 
! *****************************************************************************
MODULE hfx_communication 
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE hfx_types,                       ONLY: hfx_basis_type,&
                                             hfx_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE mathconstants
  USE message_passing,                 ONLY: mp_isendrecv,&
                                             mp_sync,&
                                             mp_waitall
  USE orbital_pointers
  USE particle_types,                  ONLY: particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE sparse_matrix_types,             ONLY: first_block_node,&
                                             get_block_node,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
   
  PUBLIC get_full_density, distribute_ks_matrix, scale_and_add_fock_to_ks_matrix
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_communication'

!***
  
  CONTAINS

! *****************************************************************************
!> \brief - Collects full density matrix from all CPUs
!> \param para_en v
!> \param full_density The full Density matrix
!> \param rho Distributed density
!> \param number_of_p_entries Maximal buffer size
!> \param is_assoc_atomic_block tells us which node owns which atomic pairs
!> \param particle_se t
!> \param nato m
!> \param last_sgf_global helper array for mapping local to global indices
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - Communication with left/right node only
!>        added a mp_sync before and after the ring of isendrecv. This *speed up* the
!>        communication, and might protect against idle neighbors flooding a busy node
!>        with messages [Joost]
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE get_full_density(para_env, full_density, rho, number_of_p_entries, is_assoc_atomic_block,&
                              natom, last_sgf_global, kind_of, basis_parameter, N_pmat, error)

    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(dp), DIMENSION(:), POINTER          :: full_density
    TYPE(real_matrix_type), POINTER          :: rho
    INTEGER, INTENT(IN)                      :: number_of_p_entries
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block
    INTEGER, INTENT(IN)                      :: natom
    INTEGER                                  :: last_sgf_global(0:natom), &
                                                kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    INTEGER(int_8), INTENT(IN)               :: N_pmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_full_density', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_from, dest, i, iatom, &
                                                icpu, ikind, istat, jatom, &
                                                jkind, mepos, ncpu, pa, pa1, &
                                                pb, pb1, req(2), source
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: recbuffer, sendbuffer, &
                                                swapbuffer
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

    failure = .FALSE.

    full_density = 0.0_dp
    ALLOCATE(sendbuffer(number_of_p_entries),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(recbuffer(number_of_p_entries),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    i=1
    DO iatom=1,rho%nblock_row
      block_node => first_block_node(rho,iatom)
      DO WHILE (ASSOCIATED(block_node))
        CALL get_block_node(block_node=block_node,&
                            block_col=jatom,&
                            BLOCK=sparse_block)
        DO pb=1,SIZE(sparse_block,2)
          DO pa=1,SIZE(sparse_block,1)
            sendbuffer(i)= sparse_block(pa,pb)
            i=i+1
          END DO
        END DO
        block_node => next_block_node(block_node)
      END DO
    END DO

    ! sync before/after a ring of isendrecv
    CALL mp_sync(para_env%group)
    ncpu = para_env%num_pe
    mepos = para_env%mepos
    dest  =MODULO(mepos+1,ncpu)
    source=MODULO(mepos-1,ncpu)
    DO icpu = 0, ncpu-1
      IF (icpu .NE. ncpu-1) THEN
        CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                          para_env%group,req(1),req(2),13)
      ENDIF
      i = 1
      data_from = MODULO(mepos - icpu,ncpu)
      DO iatom=1,natom
        ikind = kind_of(iatom)
        pa1 = last_sgf_global(iatom) - basis_parameter(ikind)%nsgf_total
        DO jatom=iatom,natom
          jkind =kind_of(jatom)
          pb1 = last_sgf_global(jatom) - basis_parameter(jkind)%nsgf_total
          IF(is_assoc_atomic_block(iatom,jatom)-1 == MODULO(data_from,ncpu)) THEN
            DO pb = pb1+1, pb1+ basis_parameter(jkind)%nsgf_total 
              DO pa = pa1+1,pa1+ basis_parameter(ikind)%nsgf_total
                full_density(get_1d_idx(pa,pb,N_pmat)) = sendbuffer(i)
                i = i + 1
              END DO
            END DO
          END IF
        END DO
      END DO
      IF (icpu .NE. ncpu -1) THEN
        CALL mp_waitall(req)
      ENDIF
      swapbuffer=>sendbuffer
      sendbuffer=>recbuffer
      recbuffer=>swapbuffer
    ENDDO
    DEALLOCATE(sendbuffer, recbuffer, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ! sync before/after a ring of isendrecv
    CALL mp_sync(para_env%group)

  END SUBROUTINE get_full_density

! *****************************************************************************
!> \brief - Distributes the local full Kohn-Sham matrix to all CPUS
!> \param para_en v
!> \param full_ks The full Kohn-Sham matrix
!> \param ks_matrix Distributed Kohn-Sham matrix
!> \param number_of_p_entries Maximal buffer size
!> \param is_assoc_atomic_block tells us which node owns which atomic pairs 
!> \param particle_se t
!> \param nato m
!> \param last_sgf_global helper array for mapping local to global indices
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - Communication with left/right node only
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE distribute_ks_matrix(para_env, full_ks, ks_matrix, number_of_p_entries, is_assoc_atomic_block,&
                                  natom, last_sgf_global, kind_of, basis_parameter, N_pmat, error)

    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(dp), DIMENSION(:)                   :: full_ks
    TYPE(real_matrix_type), POINTER          :: ks_matrix
    INTEGER, INTENT(IN)                      :: number_of_p_entries
    INTEGER, DIMENSION(:, :)                 :: is_assoc_atomic_block
    INTEGER, INTENT(IN)                      :: natom
    INTEGER                                  :: last_sgf_global(0:natom), &
                                                kind_of(*)
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    INTEGER(int_8), INTENT(IN)               :: N_pmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'distribute_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: data_to, dest, i, iatom, &
                                                icpu, ikind, istat, jatom, &
                                                jkind, mepos, ncpu, pa, pa1, &
                                                pb, pb1, req(2), source
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: recbuffer, sendbuffer, &
                                                swapbuffer
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(real_block_node_type), POINTER      :: block_node

    ALLOCATE(sendbuffer(number_of_p_entries), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    sendbuffer =0.0_dp
    ALLOCATE(recbuffer(number_of_p_entries), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    recbuffer =0.0_dp

    ncpu = para_env%num_pe
    mepos = para_env%mepos
    dest  =MODULO(mepos+1,ncpu)
    source=MODULO(mepos-1,ncpu)

    ! sync before/after a ring of isendrecv
    CALL mp_sync(para_env%group)
    DO icpu = 1,ncpu
      i = 1
      data_to = mepos-icpu
      DO iatom=1,natom
        ikind = kind_of(iatom)
        pa1 = last_sgf_global(iatom) - basis_parameter(ikind)%nsgf_total
        DO jatom=iatom,natom
          jkind = kind_of(jatom)
          pb1 = last_sgf_global(jatom) - basis_parameter(jkind)%nsgf_total
          IF(is_assoc_atomic_block(iatom,jatom)-1 == MODULO(data_to,ncpu)) THEN
            DO pb=pb1+1,pb1+ basis_parameter(jkind)%nsgf_total
              DO pa=pa1+1,pa1+ basis_parameter(ikind)%nsgf_total
                sendbuffer(i) = sendbuffer(i) + full_ks(get_1d_idx(pa,pb,N_pmat))
                i = i+1
              END DO
            END DO  
          END IF
        END DO
      END DO

      IF ( icpu .EQ. ncpu ) EXIT
      CALL mp_isendrecv(sendbuffer,dest,recbuffer,source, &
                        para_env%group,req(1),req(2),13)
      CALL mp_waitall(req)
      swapbuffer=>sendbuffer
      sendbuffer=>recbuffer
      recbuffer=>swapbuffer
    ENDDO
    ! sync before/after a ring of isendrecv
    CALL mp_sync(para_env%group)

    i = 1
    DO iatom=1,ks_matrix%nblock_row
      ikind = kind_of(iatom)
      block_node => first_block_node(ks_matrix,iatom)
      IF (.NOT. ASSOCIATED(block_node)) CYCLE
      pa1 = last_sgf_global(iatom) - basis_parameter(ikind)%nsgf_total 
      block_node => first_block_node(ks_matrix,iatom)
      DO WHILE (ASSOCIATED(block_node))
        CALL get_block_node(block_node=block_node,&
                            block_col=jatom,&
                            BLOCK=sparse_block)
        block_node => next_block_node(block_node)
        jkind = kind_of(jatom)
        pb1 = last_sgf_global(jatom) - basis_parameter(jkind)%nsgf_total
        DO pb = 1,SIZE(sparse_block,2)
          DO pa = 1,SIZE(sparse_block,1)
            sparse_block(pa,pb) = sendbuffer(i) + sparse_block(pa,pb)
            i = i + 1
          END DO
        END DO 
      END DO
    END DO
    
    DEALLOCATE(sendbuffer, recbuffer, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
            
  END SUBROUTINE distribute_ks_matrix

! *****************************************************************************
!> \brief - Distributes the local full Kohn-Sham matrix to all CPUS. Is called in
!>        case of adiabatic rescaling. This is just a refactored version of
!>        distribute_ks_matrix
!> \param para_en v
!> \param full_ks The full Kohn-Sham matrix
!> \param ks_matrix Distributed Kohn-Sham matrix
!> \param number_of_p_entries Maximal buffer size
!> \param is_assoc_atomic_block tells us which node owns which atomic pairs 
!> \param particle_se t
!> \param nato m
!> \param last_sgf_global helper array for mapping local to global indices
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      - Communication with left/right node only
!> \par History
!>      11.2007 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE scale_and_add_fock_to_ks_matrix(para_env, qs_env, ks_matrix, irep ,&
                                             scaling_factor, N_pmat, error)

    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_matrix
    INTEGER, INTENT(IN)                      :: irep
    REAL(dp), INTENT(IN)                     :: scaling_factor
    INTEGER(int_8), INTENT(IN)               :: N_pmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'scale_and_add_fock_to_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, ikind, istat, j, &
                                                natom, nspins, &
                                                number_of_p_entries, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of, last_sgf_global
    INTEGER, DIMENSION(:, :), POINTER        :: is_assoc_atomic_block
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: full_ks
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(hfx_basis_type), DIMENSION(:), &
      POINTER                                :: basis_parameter
    TYPE(hfx_type), POINTER                  :: actual_x_data
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!! All shared data is saved in i_thread = 1!

    actual_x_data => qs_env%x_data(irep, 1)
    number_of_p_entries = actual_x_data%number_of_p_entries
    is_assoc_atomic_block => actual_x_data%is_assoc_atomic_block
    basis_parameter => actual_x_data%basis_parameter
   
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,&
                    error=error)

    nspins = qs_env%dft_control%nspins    
    natom = SIZE(particle_set,1)
    ALLOCATE(kind_of(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             kind_of=kind_of)
    
    ALLOCATE(last_sgf_global(0:natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    last_sgf_global(0)=0
    DO iatom=1,natom
      ikind = kind_of(iatom)
      last_sgf_global(iatom) = last_sgf_global(iatom-1)+ basis_parameter(ikind)%nsgf_total
    END DO
    full_ks => actual_x_data%full_ks_alpha
    IF(scaling_factor /= 1.0_dp ) THEN
      DO j=1,SIZE(full_ks,1)
        full_ks(j) = full_ks(j) * scaling_factor
      END DO
    END IF
    CALL distribute_ks_matrix(para_env, full_ks, ks_matrix(1)%matrix, actual_x_data%number_of_p_entries, &
                              actual_x_data%is_assoc_atomic_block,&
                              natom, last_sgf_global, kind_of, basis_parameter, N_pmat, error)
    DEALLOCATE(actual_x_data%full_ks_alpha, STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF(nspins == 2) THEN
      full_ks => actual_x_data%full_ks_beta
      IF(scaling_factor /= 1.0_dp ) THEN
        DO j=1,SIZE(full_ks,1)
            full_ks(j) = full_ks(j) * scaling_factor
        END DO
      END IF

      CALL distribute_ks_matrix(para_env, full_ks, ks_matrix(2)%matrix, actual_x_data%number_of_p_entries, &
                                actual_x_data%is_assoc_atomic_block,&
                                natom, last_sgf_global, kind_of, basis_parameter, N_pmat, error)
      DEALLOCATE(actual_x_data%full_ks_beta, STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    DEALLOCATE(kind_of, last_sgf_global, STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE scale_and_add_fock_to_ks_matrix 

! *****************************************************************************
!> \brief Given a 2d index pair, this function returns a 1d index pair for
!>        a symmetric upper triangle NxN matrix
!>        The compiler should inline this function, therefore it appears in
!>        several modules
!> \param i,j 2d index
!> \param N matrix size
!> \par History
!>      03.2009 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  PURE FUNCTION get_1D_idx(i,j,N)
    INTEGER, INTENT(IN)                      :: i, j
    INTEGER(int_8), INTENT(IN)               :: N
    INTEGER(int_8)                           :: get_1D_idx

    INTEGER(int_8)                           :: min_ij

    min_ij = MIN(i,j)
    get_1D_idx = min_ij*N + MAX(i,j) - (min_ij-1)*min_ij/2 - N

  END FUNCTION get_1D_idx


END MODULE hfx_communication
