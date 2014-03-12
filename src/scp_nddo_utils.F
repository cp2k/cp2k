!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculation of the SCP specific terms for semi-empirical methods
!> \author CJM
! *****************************************************************************
MODULE scp_nddo_utils
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_get_diag,&
                                             cp_dbcsr_get_info,&
                                             cp_dbcsr_set_diag
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE kinds,                           ONLY: dp
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_nddo_utils'

  PUBLIC :: evaluate_scp_force


CONTAINS

! *****************************************************************************
  SUBROUTINE evaluate_scp_force( atomic_kind_set, pscp, &
       ks_scp, ascp_diag, energy_only, escp, error )

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_type), POINTER             :: pscp, ks_scp
    REAL(dp), POINTER                        :: ascp_diag( : )
    LOGICAL, INTENT(IN)                      :: energy_only
    REAL(dp), INTENT(out)                    :: escp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'evaluate_scp_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ii, ikind, &
                                                iparticle, j, n_el, natorb, &
                                                nkind, stat
    INTEGER, POINTER                         :: atom_list( : )
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: app, asp, ass, ratio, &
                                                sum_adiag, sum_ksdiag
    REAL(KIND=dp), ALLOCATABLE               :: ksdiag( : )
    REAL(KIND=dp), DIMENSION(9, 9)           :: a_scp
    REAL(KIND=dp), POINTER                   :: ks_block( :, : ), &
                                                p_block( :, : )
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(semi_empirical_type), POINTER       :: se_kind

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    nkind = SIZE(atomic_kind_set)
    escp = 0.0_dp

!    WRITE ( *, * ) 'XXXXXXXXX'
!    CALL cp_dbcsr_print ( ks_scp, matlab_format=.TRUE., error=error )
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,scp_se_parameter=se_kind, atom_list=atom_list)
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

          CALL cp_dbcsr_get_block_p(matrix=pscp,&
               row=iatom,col=iatom,BLOCK=p_block,found=found)
          CPPrecondition(ASSOCIATED(p_block),cp_failure_level,routineP,error,failure)
          DO i=1,SIZE(p_block,1)
             DO j=1,SIZE(p_block,2)
                escp = escp + p_block(i,j)*p_block (i,j)/a_scp(i,j)
             END DO
          END DO
          IF(.NOT.energy_only) THEN
             CALL cp_dbcsr_get_block_p(matrix=ks_scp,&
                  row=iatom,col=iatom,BLOCK=ks_block,found=found)
             CPPrecondition(ASSOCIATED(ks_block),cp_failure_level,routineP,error,failure)
             DO i=1,SIZE(ks_block,1)
                DO j=1,SIZE(ks_block,2)
                   ks_block(i,j) = a_scp (i,j) * ks_block(i,j)+ p_block(i,j)
                END DO
             END DO
          END IF
       END DO ! iparticle_local
    END DO ! ikind
    escp=escp*0.5_dp
! compute constraint force
    IF(.NOT.energy_only) THEN
       CALL cp_dbcsr_get_info ( ks_scp, nfullrows_total=n_el )
       ALLOCATE (ksdiag(n_el),STAT=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL cp_dbcsr_get_diag ( ks_scp, ksdiag, error=error )
      sum_ksdiag = SUM ( ksdiag )
      sum_adiag = SUM ( ascp_diag )
      ratio = sum_ksdiag/sum_adiag
      DO ii = 1, n_el
         ksdiag ( ii ) = ksdiag ( ii ) - ascp_diag ( ii ) * ratio
      END DO
!      CALL cp_dbcsr_print ( ks_scp, matlab_format=.TRUE., error=error )
      CALL cp_dbcsr_set_diag ( ks_scp, ksdiag, error=error )
      DEALLOCATE (ksdiag,STAT=stat)
      CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

  END SUBROUTINE evaluate_scp_force

! *****************************************************************************

END MODULE scp_nddo_utils

