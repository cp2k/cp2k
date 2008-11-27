!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the SCP specific terms for semi-empirical methods
!> \author CJM
! *****************************************************************************
MODULE scp_nddo_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             get_matrix_diagonal,&
                                             get_matrix_info,&
                                             real_matrix_type,&
                                             set_matrix_diagonal, &
                                             mpsum_dbm
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_nddo_utils'

  PUBLIC :: evaluate_scp_kernel, enforce_scp_nddo_constraint


CONTAINS
! *****************************************************************************
!> \brief   Allocation of scp coefficient arrays for scp-nddo
!> \author  CJM
!> \par     History
!> \date    11.03.2002
!> \version 1.0
! *****************************************************************************
  SUBROUTINE enforce_scp_nddo_constraint ( pscp, ascp_diag, error )
    TYPE(real_matrix_type), POINTER          :: pscp
    REAL(dp), POINTER                        :: ascp_diag( : )
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'enforce_scp_nddo_constraint', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, n_el, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ratio, sum_adiag, sum_pdiag
    REAL(KIND=dp), ALLOCATABLE               :: pdiag( : )

    failure = .FALSE.
    CALL get_matrix_info ( pscp, nrow=n_el )
    ALLOCATE (pdiag(n_el),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_matrix_diagonal ( pscp, pdiag )
    sum_pdiag = SUM ( pdiag )
    sum_adiag = SUM ( ascp_diag )
    ratio = sum_pdiag/sum_adiag
    DO ii = 1, n_el
       pdiag ( ii ) = pdiag ( ii ) - ascp_diag ( ii ) * ratio 
    END DO
    CALL set_matrix_diagonal ( pscp, pdiag )
    DEALLOCATE (pdiag,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE enforce_scp_nddo_constraint

! *****************************************************************************
  SUBROUTINE evaluate_scp_kernel( atomic_kind_set, pscp, &
       ks_scp, para_env, energy_only, escp, error )

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(real_matrix_type), POINTER          :: pscp, ks_scp
    TYPE(cp_para_env_type), POINTER          :: para_env
    LOGICAL, INTENT(IN)                      :: energy_only
    REAL(dp), INTENT(out)                    :: escp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_scp_kernel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ikind, &
                                                j, natorb, &
                                                nkind, iparticle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: app, asp, ass
    REAL(KIND=dp), DIMENSION(9, 9)           :: a_scp
    REAL(KIND=dp), POINTER                   :: ks_block( :, : ), &
                                                p_block( :, : )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(semi_empirical_type), POINTER       :: se_kind
    INTEGER, POINTER                         :: atom_list ( : )

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    nkind = SIZE(atomic_kind_set)
    escp = 0.0_dp

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,se_parameter=se_kind, atom_list=atom_list)
       CALL get_se_param(se_kind, ass=ass, app=app, asp=asp, natorb=natorb)
       IF (natorb>4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="SCP-NDDO not implemened for d-orbitals!", &
               error=error, error_level=cp_failure_level)
       END IF
       a_scp( 1  ,  1) = ass
       a_scp( 1  ,2:4) = asp
       a_scp( 2:4,  1) = asp
       a_scp( 2:4,2:4) = app

       DO iparticle=1, SIZE ( atom_list )
          iatom = atom_list ( iparticle )


          CALL get_block_node(matrix=pscp,&
                              block_row=iatom,&
                              block_col=iatom,&
                              BLOCK=p_block)
          DO i=1,SIZE(p_block,1)
             DO j=1,SIZE(p_block,2)
                escp = escp + p_block(i,j)*p_block (i,j)/a_scp(i,j)
             END DO
          END DO
          IF(.NOT.energy_only) THEN
             CALL get_block_node(matrix=ks_scp,&
                                 block_row=iatom,&
                                 block_col=iatom,&
                                 BLOCK=ks_block)
             DO i=1,SIZE(ks_block,1) 
                DO j=1,SIZE(ks_block,2) 
                   ks_block(i,j) = a_scp (i,j) * ks_block(i,j)+ p_block(i,j)
                END DO
             END DO
          END IF
       END DO ! iparticle_local
    END DO ! ikind 
    ! collect escp from all processors
    escp=escp*0.5_dp
  END SUBROUTINE evaluate_scp_kernel

! *****************************************************************************

END MODULE scp_nddo_utils

