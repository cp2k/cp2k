!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  Methods to apply GLE to PI runs.
!> \author michelec
!> \par History
!>      06.2010 created [michelec]
!> \note   trying to keep duplication at a minimum....
! *****************************************************************************

MODULE pint_gle
  USE gle_system_dynamics,             ONLY: gle_cholesky_stab
  USE gle_system_types,                ONLY: gle_type
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: next_random_number,&
                                             rng_stream_type
  USE pint_types,                      ONLY: pint_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pint_gle_step, pint_gle_init, pint_calc_gle_energy

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pint_gle'

  CONTAINS
! *****************************************************************************
!> \brief ...
!> \param pint_env ...
! *****************************************************************************
  SUBROUTINE pint_calc_gle_energy( pint_env)
    TYPE(pint_env_type), POINTER             :: pint_env

    INTEGER                                  :: i

    pint_env%e_gle=0._dp
    IF (ASSOCIATED(pint_env%gle)) THEN
      DO i = 1, pint_env%gle%loc_num_gle
        pint_env%e_gle=pint_env%e_gle+pint_env%gle%nvt(i)%thermostat_energy
      END DO
    END IF
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param pint_env ...
! *****************************************************************************
  SUBROUTINE pint_gle_init( pint_env)
    TYPE(pint_env_type), POINTER             :: pint_env

    INTEGER                                  :: i, ib, idim, imap, j
    LOGICAL                                  :: failure
    REAL(dp) :: mf, rr(pint_env%gle%ndim), &
      cc(pint_env%gle%ndim,pint_env%gle%ndim)
    TYPE(rng_stream_type), POINTER           :: rng_stream

    failure = .FALSE.

    CALL gle_cholesky_stab(pint_env%gle%c_mat, cc, pint_env%gle%ndim)
    DO i = 1, pint_env%gle%loc_num_gle
      imap = pint_env%gle%map_info%index(i)
      ib=1+(imap-1)/pint_env%ndim
      idim=1+MOD(imap-1,pint_env%ndim)
      mf=1.0_dp/SQRT(pint_env%mass_fict(ib,idim))
      rng_stream => pint_env%gle%nvt(i)%gaussian_rng_stream
      DO j = 1, pint_env%gle%ndim
        rr(j) =  next_random_number( rng_stream)*mf
      END DO
      pint_env%gle%nvt(i)%s = MATMUL(cc,rr)
    END DO

  END SUBROUTINE pint_gle_init

! *****************************************************************************
!> \brief ...
!> \param pint_env ...
! *****************************************************************************
  SUBROUTINE pint_gle_step( pint_env)
    TYPE(pint_env_type), POINTER             :: pint_env

    CHARACTER(len=*), PARAMETER :: routineN = 'pint_gle_step', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iadd, ib, ideg, idim, &
                                                imap, ndim, num
    REAL(dp)                                 :: alpha, beta, mf, rr
    REAL(dp), DIMENSION(:, :), POINTER       :: a_mat, e_tmp, h_tmp, s_tmp
    TYPE(gle_type), POINTER                  :: gle
    TYPE(rng_stream_type), POINTER           :: rng_stream

    CALL timeset(routineN,handle)

    gle => pint_env%gle
    NULLIFY( rng_stream)
    ndim = gle%ndim

    ALLOCATE (s_tmp(ndim,gle%loc_num_gle))
    s_tmp = 0.0_dp
    ALLOCATE (e_tmp(ndim,gle%loc_num_gle))
    ALLOCATE (h_tmp(ndim,gle%loc_num_gle))

    DO ideg = 1,gle%loc_num_gle
       imap = gle%map_info%index(ideg)
       ib=1+(imap-1)/pint_env%ndim
       idim=1+MOD(imap-1,pint_env%ndim)

       gle%nvt(ideg)%s(1)= pint_env%uv_t(ib,idim)
       gle%nvt(ideg)%thermostat_energy=gle%nvt(ideg)%thermostat_energy &
               +0.5_dp*pint_env%mass_fict(ib,idim)*gle%nvt(ideg)%s(1)**2
       s_tmp(1,imap) = gle%nvt(ideg)%s(1)
       rng_stream => gle%nvt(ideg)%gaussian_rng_stream
       rr =  next_random_number( rng_stream)
       mf=1.0_dp/SQRT(pint_env%mass_fict(ib,idim))
       e_tmp(1,imap)= rr*mf
       DO iadd = 2,ndim
          s_tmp(iadd,imap) = gle%nvt(ideg)%s(iadd)
          rr =  next_random_number( rng_stream)
          e_tmp(iadd,imap) = rr*mf
       END DO
    END DO
    num =  gle%loc_num_gle
    a_mat => gle%gle_s
    alpha = 1.0_dp
    beta = 0.0_dp


    CALL DGEMM('N','N',ndim,num,ndim,alpha,a_mat(1,1),ndim,e_tmp(1,1),ndim,beta,h_tmp(1,1),ndim)

    a_mat => gle%gle_t
    beta = 1.0_dp
    CALL DGEMM("N","N",ndim,num,ndim,alpha,a_mat(1,1),ndim,s_tmp(1,1),ndim,beta,h_tmp(1,1),ndim )


    DO ideg = 1,gle%loc_num_gle
       imap = gle%map_info%index(ideg)

       DO iadd = 1,ndim
          gle%nvt(ideg)%s(iadd) = h_tmp(iadd,imap)
       END DO
       ib=1+(imap-1)/pint_env%ndim
       idim=1+MOD(imap-1,pint_env%ndim)
       pint_env%uv_t(ib,idim)=gle%nvt(ideg)%s(1)
       gle%nvt(ideg)%thermostat_energy=gle%nvt(ideg)%thermostat_energy &
               -0.5_dp*pint_env%mass_fict(ib,idim)*gle%nvt(ideg)%s(1)**2
    END DO
       pint_env%e_kin_t=0.0_dp
    DEALLOCATE(e_tmp,s_tmp,h_tmp)
    CALL timestop(handle)
  END SUBROUTINE
END MODULE
