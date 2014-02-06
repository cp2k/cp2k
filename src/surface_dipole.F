!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
! *****************************************************************************
MODULE surface_dipole

  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE message_passing,                 ONLY: mp_sum
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'surface_dipole'

  PUBLIC :: calc_dipsurf_potential


CONTAINS

! *****************************************************************************
!> \brief compute the surface dipole and the correction to the hartree potential
!> \param qs_env the qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      01.2014 created [MI]
!> \author MI
! *****************************************************************************


  SUBROUTINE calc_dipsurf_potential(qs_env,ks_env,energy,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_dipsurf_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, idir_surfdip, &
                                                ilayer_min, ilow, irho, &
                                                ispin, isurf, iup, jsurf, stat
    INTEGER, DIMENSION(3)                    :: ngrid
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dh(3,3), dip_fac, dip_hh, &
                                                dsurf, height_min, hh, &
                                                rhoav_min, surfarea, vdip, &
                                                vdip_fac
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: rhoavsurf
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: vdip_r, wf_r
    TYPE(pw_p_type), POINTER                 :: rho0_s_gs, rho_core
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(cell, dft_control, rho, pw_env, rho_core, rho0_s_gs, auxbas_pw_pool, &
             pw_pools,subsys,particles)

      CALL get_qs_env(qs_env, dft_control=dft_control, rho=rho, rho_core=rho_core, &
               rho0_s_gs=rho0_s_gs, cell=cell, pw_env=pw_env,  subsys=subsys, error=error)
      CALL cp_subsys_get(subsys,particles=particles,error=error)

      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools,error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool,vdip_r%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)

      IF (dft_control%qs_control%gapw) THEN
          CALL pw_transfer(rho0_s_gs%pw,wf_r%pw,error=error)
          IF(dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
              CALL pw_axpy(rho_core%pw,wf_r%pw, error=error)
          END IF
      ELSE
          CALL pw_transfer(rho_core%pw,wf_r%pw,error=error)
      END IF
      DO ispin=1,dft_control%nspins
          CALL pw_axpy(rho%rho_r(ispin)%pw,wf_r%pw, error=error)
      END DO

      ngrid(1:3) = wf_r%pw%pw_grid%npts(1:3)
      idir_surfdip = dft_control%dir_surf_dip

      
      DO i = 1,3
       IF(i /= idir_surfdip) THEN
           IF(ABS(wf_r%pw%pw_grid%dh(idir_surfdip,i)) > 1.E-7_dp) THEN
             ! stop surface dipole defined only for slab perpendigular to one of the Cartesian axis
              CALL cp_assert(.FALSE., cp_failure_level,cp_assertion_failed,&
                   routineP, " Dipole correction only for surface perpendicular to "//&
                   " one Cartesian axis",error, only_ionode=.TRUE.)
!  not properly general, we assume that vectors A, B, and C are along x y and z respectively, 
!  in the ortorhombic cell, but in principle it does not need to be this way, importan
!  is that the cell angles are 90 degrees.
           END IF
       END IF
      ENDDO


      ilow = wf_r%pw%pw_grid%bounds(1,idir_surfdip)
      iup = wf_r%pw%pw_grid%bounds(2,idir_surfdip)

      ALLOCATE (rhoavsurf(ilow:iup),stat=stat)
      CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
      rhoavsurf = 0.0_dp

      bo => wf_r%pw%pw_grid%bounds_local
      dh = wf_r%pw%pw_grid%dh

      IF(idir_surfdip==3) THEN
        isurf = 1
        jsurf = 2

        DO i = bo(1,3),bo(2,3)
           rhoavsurf(i) = accurate_sum(wf_r%pw%cr3d(bo(1,1):bo(2,1),bo(1,2):bo(2,2),i))
        END DO

      ELSEIF(idir_surfdip==2) THEN
        isurf = 3
        jsurf = 1

        DO i =  bo(1,2),bo(2,2)
           rhoavsurf(i) = accurate_sum(wf_r%pw%cr3d(bo(1,1):bo(2,1),i, bo(1,3):bo(2,3)))
        END DO
      ELSE
        isurf = 2
        jsurf = 3

        DO i = bo(1,1),bo(2,1)
           rhoavsurf(i) = accurate_sum(wf_r%pw%cr3d(i,bo(1,2):bo(2,2),bo(1,3):bo(2,3)))
        END DO
      ENDIF

      surfarea = cell%hmat(isurf,isurf)*cell%hmat(jsurf,jsurf) -&
                  cell%hmat(isurf,jsurf)*cell%hmat(jsurf,isurf)
      dsurf = surfarea/REAL(ngrid(isurf)*ngrid(jsurf),dp)


      IF (wf_r%pw%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
         CALL mp_sum(rhoavsurf,wf_r%pw%pw_grid%para%group)
      END IF
      rhoavsurf(ilow:iup) = dsurf*rhoavsurf(ilow:iup)


      ! locate where the vacuum is, and set the reference point for the calculation of the dipole
      rhoavsurf(ilow:iup) = rhoavsurf(ilow:iup)/surfarea
      ilayer_min = ilow -1 + MINLOC(ABS(rhoavsurf(ilow:iup)),1)  
      rhoav_min = ABS(rhoavsurf(ilayer_min))
      CALL cp_assert(rhoav_min<1.E-5_dp, cp_failure_level,cp_assertion_failed,&
                routineP, " Dipole correction needs more vacuum space above the surface ",&
                error, only_ionode=.TRUE.)
      height_min =  REAL((ilayer_min - ilow),dp)*dh(idir_surfdip,idir_surfdip)


!   surface dipole form average rhoavsurf
!   \sum_i NjdjNkdkdi rhoav_i (i-imin)di
     dip_hh = 0.0_dp
     dip_fac = wf_r%pw%pw_grid%vol * dh(idir_surfdip,idir_surfdip)/REAL(ngrid(idir_surfdip),dp)

     DO i = ilayer_min+1,ilayer_min+ngrid(idir_surfdip)
        hh = REAL((i - ilayer_min),dp)
        IF(i>iup) THEN
          irho =  i - ngrid(idir_surfdip)
        ELSE
          irho = i
        END IF
        dip_hh = dip_hh + rhoavsurf(irho)*hh*dip_fac
     END DO 

     DEALLOCATE(rhoavsurf,stat=stat)
     CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

!    Calculation of the dipole potential as a function of the perpendicular coordinate
     CALL pw_zero(vdip_r%pw, error=error)
     vdip_fac = dip_hh*4.0_dp*pi

     DO i = ilayer_min+1,ilayer_min+ngrid(idir_surfdip)
       hh = REAL((i - ilayer_min),dp)*dh(idir_surfdip,idir_surfdip)
       vdip = vdip_fac*(-0.5_dp+(hh/ cell%hmat(idir_surfdip,idir_surfdip)))*&
                    ks_env%v_hartree_rspace%pw%pw_grid%dvol/surfarea
        IF(i>iup) THEN
          irho =  i - ngrid(idir_surfdip)
        ELSE
          irho = i
        END IF
        IF(idir_surfdip==3) THEN
          vdip_r%pw%cr3d(bo(1,1):bo(2,1),bo(1,2):bo(2,2),irho) =&
                  vdip_r%pw%cr3d(bo(1,1):bo(2,1),bo(1,2):bo(2,2),irho) + vdip
        ELSEIF(idir_surfdip==2) THEN
          IF(irho>=bo(1,2) .AND. irho<=bo(2,2)) THEN
            vdip_r%pw%cr3d(bo(1,1):bo(2,1),irho,bo(1,3):bo(2,3)) =&
                  vdip_r%pw%cr3d(bo(1,1):bo(2,1),irho,bo(1,3):bo(2,3)) + vdip
          END IF
        ELSE
          IF(irho>=bo(1,1) .AND. irho<=bo(2,1)) THEN
            vdip_r%pw%cr3d(irho,bo(1,2):bo(2,2),bo(1,3):bo(2,3)) =&
                  vdip_r%pw%cr3d(irho,bo(1,2):bo(2,2),bo(1,3):bo(2,3)) + vdip
          END IF
        END IF

     END DO

!    Dipole correction contribution to the energy
     energy%surf_dipole = 0.5_dp * accurate_sum(vdip_r%pw%cr3d* wf_r%pw%cr3d)
     IF (wf_r%pw%pw_grid%para%mode /= PW_MODE_LOCAL) THEN
         CALL mp_sum(energy%surf_dipole,wf_r%pw%pw_grid%para%group)
     END IF

!    Add the dipole potential to the hartree potential on the realspace grid
     ks_env%v_hartree_rspace%pw%cr3d =  ks_env%v_hartree_rspace%pw%cr3d + vdip_r%pw%cr3d

     CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
     CALL pw_pool_give_back_pw(auxbas_pw_pool,vdip_r%pw, error=error)

     CALL timestop(handle)

  END SUBROUTINE  calc_dipsurf_potential



END MODULE surface_dipole
