!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains information regarding the decoupling/recoupling method of Bloechl
!> \author Teodoro Laino
! *****************************************************************************
MODULE cp_ddapc_methods
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erfc
  USE f77_blas
  USE input_constants,                 ONLY: weight_type_mass,&
                                             weight_type_unit
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi,&
                                             oorootpi,&
                                             pi,&
                                             twopi
  USE mathlib,                         ONLY: diamat_all,&
                                             invert_matrix
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc
  USE pw_types,                        ONLY: pw_type
  USE spherical_harmonics,             ONLY: legendre
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc_methods'  
  PUBLIC :: ddapc_eval_gfunc,&
            build_b_vector,&
            build_der_b_vector,&
            build_A_matrix,&
            build_der_A_matrix,&
            ddapc_eval_AmI,&
            ewald_ddapc_pot,&
            solvation_ddapc_pot

CONTAINS

! *****************************************************************************
  SUBROUTINE ddapc_eval_gfunc(gfunc, w, gcut, rho_tot_g, radii, error)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: gfunc
    REAL(kind=dp), DIMENSION(:), POINTER     :: w
    REAL(KIND=dp), INTENT(IN)                :: gcut
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(kind=dp), DIMENSION(:), POINTER     :: radii
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_eval_gfunc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: e_dim, handle, ig, igauss, &
                                                s_dim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: g2, gcut2, rc, rc2

    CALL timeset(routineN,handle)
    gcut2 = gcut * gcut
    failure = .FALSE.
    IF (.NOT.failure) THEN
       !
       s_dim = rho_tot_g % pw_grid % first_gne0
       e_dim = rho_tot_g % pw_grid % ngpts_cut_local
       ALLOCATE(gfunc(s_dim:e_dim,SIZE(radii)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(w(s_dim:e_dim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       gfunc = 0.0_dp
       w     = 0.0_dp
       DO igauss = 1, SIZE(radii)
          rc  = radii(igauss)
          rc2 = rc*rc
          DO ig = s_dim, e_dim
             g2 = rho_tot_g % pw_grid % gsq ( ig )
             IF (g2 > gcut2) EXIT
             gfunc(ig,igauss) = EXP(-g2*rc2/4.0_dp)
          END DO
       END DO
       DO ig = s_dim, e_dim
          g2 = rho_tot_g % pw_grid % gsq ( ig )
          IF (g2 > gcut2) EXIT
          w(ig) =  fourpi*(g2 -gcut2)**2/(g2*gcut2)
       END DO
    END IF
    CALL timestop(handle)
  END SUBROUTINE ddapc_eval_gfunc

! *****************************************************************************
!> \brief Computes the B vector for the solution of the linear system
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE build_b_vector(bv, gfunc, w, particle_set, radii, rho_tot_g, gcut, error)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: bv
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: gfunc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: w
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_b_vector', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: phase
    INTEGER                                  :: e_dim, handle, idim, ig, &
                                                igauss, igmax, iparticle, &
                                                s_dim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, g2, gcut2, gvec(3), &
                                                rvec(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_bv, my_bvw

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (.NOT.failure) THEN
       NULLIFY(my_bv, my_bvw)
       gcut2 = gcut * gcut
       s_dim = rho_tot_g % pw_grid % first_gne0
       e_dim = rho_tot_g % pw_grid % ngpts_cut_local
       igmax = 0
       DO ig = s_dim, e_dim
             g2 = rho_tot_g % pw_grid % gsq ( ig )
             IF (g2 > gcut2) EXIT
             igmax = ig
       ENDDO
       ALLOCATE(my_bv(s_dim:igmax),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(my_bvw(s_dim:igmax),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !   
       DO iparticle = 1, SIZE(particle_set)
          rvec =  particle_set(iparticle)%r
          my_bv = 0.0_dp
          DO ig = s_dim, igmax
             gvec  = rho_tot_g % pw_grid % g (:,ig)
             arg   = DOT_PRODUCT(gvec,rvec)
             phase = CMPLX(COS(arg),-SIN(arg),KIND=dp)
             my_bv (ig) =  w(ig) * REAL(CONJG(rho_tot_g%cc(ig))*phase,KIND=dp)
          END DO
          DO igauss = 1, SIZE(radii)
             idim = (iparticle-1)*SIZE(radii)+igauss
             DO ig = s_dim, igmax
                my_bvw(ig) = my_bv(ig)* gfunc(ig,igauss)
             END DO
             bv (idim) = accurate_sum(my_bvw)
          END DO
       END DO
       DEALLOCATE(my_bvw,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(my_bv,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE build_b_vector
    
! *****************************************************************************
!> \brief Computes the A matrix for the solution of the linear system
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE build_A_matrix(Am, gfunc, w, particle_set, radii, rho_tot_g, gcut, error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: Am
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: gfunc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: w
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_A_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: e_dim, handle, idim1, idim2, ig, igauss1, igauss2, igmax, &
      iparticle1, iparticle2, istart_g, s_dim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, g2, gcut2, gvec(3), &
                                                rvec1(3), rvec2(3), &
                                                rvec_diff(3), tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: my_Am, my_Amw

    failure = .FALSE.

    CALL timeset(routineN,handle)
    IF (.NOT. failure) THEN
       gcut2 = gcut * gcut
       s_dim = rho_tot_g % pw_grid % first_gne0
       e_dim = rho_tot_g % pw_grid % ngpts_cut_local
       igmax=0
       DO ig = s_dim, e_dim
          g2 = rho_tot_g % pw_grid % gsq ( ig )
          IF (g2 > gcut2) EXIT
          igmax  = ig
       ENDDO
       ALLOCATE(my_Am(s_dim:igmax),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(my_Amw(s_dim:igmax),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iparticle1 = 1, SIZE(particle_set)
          rvec1=particle_set(iparticle1)%r
          DO iparticle2 = iparticle1, SIZE(particle_set)
             rvec2=particle_set(iparticle2)%r
             rvec_diff=rvec1-rvec2
             DO ig = s_dim, igmax
                gvec   = rho_tot_g % pw_grid % g (:,ig)
                arg    = DOT_PRODUCT(gvec,rvec_diff)
                my_Am (ig) = w(ig) * COS(arg)
             END DO
             DO igauss1 = 1, SIZE(radii)
                idim1 = (iparticle1-1)*SIZE(radii)+igauss1
                istart_g = 1
                IF (iparticle2==iparticle1) istart_g = igauss1
                DO igauss2 = istart_g, SIZE(radii)
                   idim2 = (iparticle2-1)*SIZE(radii)+igauss2
                   DO ig = s_dim, igmax
                      my_Amw(ig) = my_Am(ig)*gfunc(ig,igauss1)*gfunc(ig,igauss2)
                   END DO
                   tmp = accurate_sum(my_Amw)
                   Am (idim2,idim1) = tmp
                   Am (idim1,idim2) = tmp
                END DO
             END DO
          END DO
       END DO
       DEALLOCATE(my_Amw,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(my_Am,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE build_A_matrix

! *****************************************************************************
!> \brief Computes the derivative of B vector for the evaluation of the Pulay forces
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE build_der_b_vector(dbv, gfunc, w, particle_set, radii, rho_tot_g, gcut, iparticle0, error)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: dbv
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: gfunc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: w
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(IN)                      :: iparticle0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_der_b_vector', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: dphase
    INTEGER                                  :: e_dim, handle, idim, ig, &
                                                igauss, igmax, iparticle, &
                                                s_dim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, g2, gcut2
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: my_dbvw
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_dbv
    REAL(KIND=dp), DIMENSION(3)              :: gvec, rvec

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (.NOT.failure) THEN
       gcut2 = gcut * gcut
       s_dim = rho_tot_g % pw_grid % first_gne0
       e_dim = rho_tot_g % pw_grid % ngpts_cut_local
       igmax  = 0
       DO ig = s_dim, e_dim
             g2 = rho_tot_g % pw_grid % gsq ( ig )
             IF (g2 > gcut2) EXIT
             igmax  = ig
       ENDDO
       ALLOCATE(my_dbv(3,s_dim:igmax),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(my_dbvw(s_dim:igmax),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iparticle = 1, SIZE(particle_set)
          IF (iparticle /= iparticle0) CYCLE
          rvec =  particle_set(iparticle)%r
          DO ig = s_dim, igmax
             gvec   = rho_tot_g % pw_grid % g (:,ig)
             arg    = DOT_PRODUCT(gvec,rvec)
             dphase = - CMPLX(SIN(arg),COS(arg),KIND=dp) 
             my_dbv (:,ig) = w(ig) * REAL(CONJG(rho_tot_g%cc(ig))*dphase,KIND=dp) * gvec(:)
          END DO
          DO igauss = 1, SIZE(radii)
             idim = (iparticle-1)*SIZE(radii)+igauss
             DO ig = s_dim, igmax
                 my_dbvw(ig) = my_dbv(1,ig)*gfunc(ig,igauss)
             END DO
             dbv (idim,1) = accurate_sum(my_dbvw)
             DO ig = s_dim, igmax
                 my_dbvw(ig) = my_dbv(2,ig)*gfunc(ig,igauss)
             END DO
             dbv (idim,2) = accurate_sum(my_dbvw)
             DO ig = s_dim, igmax
                 my_dbvw(ig) = my_dbv(3,ig)*gfunc(ig,igauss)
             END DO
             dbv (idim,3) = accurate_sum(my_dbvw)
          END DO
       END DO
       DEALLOCATE(my_dbvw,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(my_dbv,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE build_der_b_vector
    
! *****************************************************************************
!> \brief Computes the derivative of the A matrix for the evaluation of the 
!>      Pulay forces
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE build_der_A_matrix(dAm, gfunc, w, particle_set, radii,&
       rho_tot_g, gcut, iparticle0, error)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(INOUT)                          :: dAm
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: gfunc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: w
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), INTENT(IN)                :: gcut
    INTEGER, INTENT(IN)                      :: iparticle0
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_der_A_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER :: e_dim, handle, idim1, idim2, ig, igauss1, igauss2, igmax, &
      iparticle1, iparticle2, istart_g, s_dim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg1, factor, g2, gcut2, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: my_Amw
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: my_Am
    REAL(KIND=dp), DIMENSION(3)              :: gvec, rvec1, rvec2, rvec_diff

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (.NOT. failure) THEN
       gcut2 = gcut * gcut
       s_dim = rho_tot_g % pw_grid % first_gne0
       e_dim = rho_tot_g % pw_grid % ngpts_cut_local
       igmax = 0
       DO ig = s_dim, e_dim
          g2 = rho_tot_g % pw_grid % gsq ( ig )
          IF (g2 > gcut2) EXIT
          igmax  = ig
       ENDDO
       ALLOCATE(my_Am(3,s_dim:igmax),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(my_Amw(s_dim:igmax),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       DO iparticle1 = 1, SIZE(particle_set)
          rvec1 =  particle_set(iparticle1)%r
          DO iparticle2 = iparticle1, SIZE(particle_set)
             rvec2 =  particle_set(iparticle2)%r
             IF ((iparticle1 /= iparticle0).AND.(iparticle2 /= iparticle0)) CYCLE
             IF    ((iparticle1 == iparticle0).AND.(iparticle2==iparticle0))  THEN
                   factor=2.0_dp
             ELSEIF (iparticle1 == iparticle0) THEN
                   factor=1.0_dp
             ELSEIF (iparticle2 == iparticle0) THEN
                   factor=-1.0_dp
             END IF
             rvec_diff=rvec2-rvec1

             DO ig = s_dim, igmax
                gvec   = rho_tot_g % pw_grid % g (:,ig)
                arg1 = DOT_PRODUCT(gvec,rvec_diff)
                my_Am(1:3,ig) = factor * w(ig) * SIN(arg1) * gvec
             END DO
             DO igauss1 = 1, SIZE(radii)
                idim1 = (iparticle1-1)*SIZE(radii)+igauss1
                istart_g = 1
                IF (iparticle2==iparticle1) istart_g = igauss1
                DO igauss2 = istart_g, SIZE(radii)
                   idim2 = (iparticle2-1)*SIZE(radii)+igauss2
                   DO ig = s_dim, igmax
                      my_Amw(ig) =  my_Am(1,ig)*gfunc(ig,igauss1)*gfunc(ig,igauss2)
                   END DO
                   tmp = accurate_sum(my_Amw)
                   dAm (idim1,idim2,1) = tmp
                   dAm (idim2,idim1,1) = tmp
                   DO ig = s_dim, igmax
                      my_Amw(ig) =  my_Am(2,ig)*gfunc(ig,igauss1)*gfunc(ig,igauss2)
                   END DO
                   tmp = accurate_sum(my_Amw)
                   dAm (idim1,idim2,2) = tmp
                   dAm (idim2,idim1,2) = tmp
                   DO ig = s_dim, igmax
                      my_Amw(ig) =  my_Am(3,ig)*gfunc(ig,igauss1)*gfunc(ig,igauss2)
                   END DO
                   tmp = accurate_sum(my_Amw)
                   dAm (idim1,idim2,3) = tmp
                   dAm (idim2,idim1,3) = tmp
                END DO
             END DO
          END DO
       END DO
       DEALLOCATE(my_Amw,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(my_Am,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE build_der_A_matrix

! *****************************************************************************
!> \brief Computes the inverse AmI of the Am matrix 
!> \par History
!>      12.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE  ddapc_eval_AmI(GAmI, c0,  gfunc, w, particle_set, gcut,&
       rho_tot_g, radii, iw, Vol, error)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: GAmI
    REAL(KIND=dp), INTENT(OUT)               :: c0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: gfunc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: w
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), INTENT(IN)                :: gcut
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    INTEGER, INTENT(IN)                      :: iw
    REAL(KIND=dp), INTENT(IN)                :: Vol
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_eval_AmI', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ndim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: condition_number, inv_error
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: AmE, cv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Am, AmI, Amw

    failure = .FALSE.
    CALL timeset(routineN,handle)
    IF (.NOT. failure) THEN
       ndim = SIZE(particle_set)*SIZE(radii)
       ALLOCATE(Am(ndim, ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(AmI(ndim, ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(GAmI(ndim, ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(cv(ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       Am  = 0.0_dp
       AmI = 0.0_dp
       cv  = 1.0_dp/Vol
       CALL build_A_matrix(Am, gfunc, w, particle_set, radii, rho_tot_g, gcut, error)
       Am = Am / (Vol*Vol)
       CALL mp_sum(Am,rho_tot_g%pw_grid%para%group)
       IF (iw>0) THEN
          ! Checking conditions numbers and eigenvalues
          ALLOCATE(Amw(ndim, ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(AmE(ndim), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          Amw = Am
          CALL diamat_all(Amw, AmE,error=error)
          condition_number = MAXVAL(ABS(AmE))/MINVAL(ABS(AmE))
          WRITE(iw,'(T3,A,3E15.9)')" Eigenvalues of Matrix A:",AmE
          WRITE(iw,'(T3,A,3E15.9)')" Condition number:",condition_number
          IF (condition_number > 1.0E12_dp) THEN
             WRITE (iw,FMT="(/,T2,A)")&
                  "WARNING: high condition number => possibly ill-conditioned matrix"
          END IF
          DEALLOCATE(Amw, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DEALLOCATE(AmE, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       CALL invert_matrix(Am, AmI, inv_error, "N",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(T3,A,F15.9)')" Error inverting the A matrix: ", inv_error
       END IF
       c0 = DOT_PRODUCT(cv,MATMUL(AmI,cv))       
       DEALLOCATE(Am, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(cv, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       GAmI = AmI
       DEALLOCATE(AmI, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(handle)
  END SUBROUTINE ddapc_eval_AmI

! *****************************************************************************
!> \brief Evaluates the Ewald term E2 and E3 energy term for the decoupling/coupling
!>      of periodic images
!> \par History
!>      08.2005 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  RECURSIVE SUBROUTINE ewald_ddapc_pot(coeff, factor, cell, multipole_section,&
       particle_set, M, radii, error)
    TYPE(pw_type), POINTER                   :: coeff
    REAL(KIND=dp), INTENT(IN)                :: factor
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: multipole_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: M
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_ddapc_pot', &
      routineP = moduleN//':'//routineN

    INTEGER :: ewmdim, handle, idim, idim1, idim2, idimo, igauss1, igauss2, &
      ip1, ip2, iparticle1, iparticle2, istart_g, k1, k2, k3, n_rep, ndim, &
      nmax1, nmax2, nmax3, r1, r2, r3, rmax1, rmax2, rmax3, sfact, stat
    LOGICAL                                  :: analyt, failure
    REAL(KIND=dp) :: alpha, eps, ew_neut, fac, fac3, fs, g_ewald, galpha, &
      gsq, gsqi, ij_fac, my_val, r, r2tmp, r_ewald, rc1, rc12, rc2, rc22, &
      rcut, rcut2, t1, tol, tol1
    REAL(KIND=dp), DIMENSION(3)              :: fvec, gvec, ra, rvec
    REAL(KIND=dp), DIMENSION(:), POINTER     :: EwM

    failure = .FALSE.
    NULLIFY(EwM)
    CALL timeset(routineN,handle)
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(M),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(radii),cp_failure_level,routineP,error,failure)
       CPPostcondition(cell%orthorhombic,cp_failure_level,routineP,error,failure)
       sfact = factor
       rcut = MIN(cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3))/2.0_dp 
       CALL section_vals_val_get(multipole_section,"RCUT",n_rep_val=n_rep,error=error)
       IF (n_rep==1) CALL section_vals_val_get(multipole_section,"RCUT",r_val=rcut,error=error)
       CALL section_vals_val_get(multipole_section,"EWALD_PRECISION",r_val=eps,error=error)
       CALL section_vals_val_get(multipole_section,"ANALYTICAL_GTERM",l_val=analyt,error=error)
       rcut2 = rcut**2
       !
       ! Setting-up parameters for Ewald summation
       !
       eps    = MIN(ABS(eps),0.5_dp)
       tol    = SQRT(ABS(LOG(eps*rcut)))
       alpha  = SQRT(ABS(LOG(eps*rcut*tol)))/rcut
       galpha = 1.0_dp/(4.0_dp*alpha*alpha)
       tol1   = SQRT(-LOG(eps*rcut*(2.0_dp*tol*alpha)**2))
       nmax1  = NINT(0.25_dp + cell%hmat(1,1)*alpha*tol1/pi)
       nmax2  = NINT(0.25_dp + cell%hmat(2,2)*alpha*tol1/pi)
       nmax3  = NINT(0.25_dp + cell%hmat(3,3)*alpha*tol1/pi)
       
       rmax1  = CEILING(rcut/cell%hmat(1,1))
       rmax2  = CEILING(rcut/cell%hmat(2,2))
       rmax3  = CEILING(rcut/cell%hmat(3,3))
       fac    = 1.e0_dp/cell%deth
       fac3   = fac*pi
       fvec   = twopi / (/cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3)/)
       ew_neut = - fac * pi / alpha ** 2
       !
       ewmdim = SIZE(particle_set) * (SIZE(particle_set)+1) / 2
       ndim   = SIZE(particle_set) * SIZE(radii)
       ALLOCATE(EwM(ewmdim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(M(ndim, ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
       M = 0.0_dp
       !
       idim = 0
       DO iparticle1 = 1, SIZE(particle_set)
          ip1 = (iparticle1-1)*SIZE(radii)          
          DO iparticle2 = 1, iparticle1
             ij_fac = 1.0_dp
             IF (iparticle1==iparticle2) ij_fac=0.5_dp

             ip2 = (iparticle2-1)*SIZE(radii)
             idim = idim + 1
             !
             ! Real-Space Contribution
             !
             my_val = 0.0_dp
             rvec   = particle_set(iparticle1)%r - particle_set(iparticle2)%r
             r_ewald = 0.0_dp
             IF (iparticle1 /= iparticle2) THEN
                ra = rvec 
                r2tmp  = DOT_PRODUCT(ra,ra)
                IF (r2tmp<=rcut2) THEN 
                   r  = SQRT(r2tmp)
                   t1 = erfc(alpha*r) / r
                   r_ewald = t1
                END IF
             END IF
             DO r1 = -rmax1, rmax1
                DO r2 = -rmax2, rmax2
                   DO r3 = -rmax3, rmax3
                      IF ((r1==0).AND.(r2==0).AND.(r3==0)) CYCLE
                      ra(1) = rvec(1) + cell%hmat(1,1)*r1 
                      ra(2) = rvec(2) + cell%hmat(2,2)*r2 
                      ra(3) = rvec(3) + cell%hmat(3,3)*r3 
                      r2tmp  = DOT_PRODUCT(ra,ra)
                      IF (r2tmp<=rcut2) THEN
                         r  = SQRT(r2tmp)
                         t1 = erfc(alpha*r) / r
                         r_ewald = r_ewald + t1*ij_fac           
                      END IF
                   END DO
                END DO
             END DO
             !
             ! G-space Contribution
             !
             IF (analyt) THEN
                g_ewald = 0.0_dp
                DO k1 = 0, nmax1
                   DO k2 = -nmax2, nmax2
                      DO k3 = -nmax3, nmax3
                         IF (k1 == 0.AND.k2 == 0.AND.k3 == 0) CYCLE
                         fs = 2.0_dp; IF (k1==0) fs = 1.0_dp
                         gvec = fvec * (/REAL(k1,KIND=dp),REAL(k2,KIND=dp),REAL(k3,KIND=dp)/)
                         gsq  = DOT_PRODUCT(gvec,gvec)
                         gsqi = fs/gsq
                         t1   = fac * gsqi * EXP(-galpha*gsq)
                         g_ewald = g_ewald + t1 * COS(DOT_PRODUCT(gvec,rvec))
                      END DO
                   END DO
                END DO
             ELSE
                g_ewald = Eval_Interp_Spl3_pbc(rvec, coeff, error)
             END IF
             !
             ! G-EWALD, R-EWALD
             !
             g_ewald = r_ewald + fourpi*g_ewald
             !
             ! Self Contribution
             !
             IF (iparticle1 == iparticle2) THEN
                g_ewald = g_ewald - 2.0_dp * alpha * oorootpi
             END IF
             !
             IF (iparticle1/=iparticle2) THEN 
                ra = rvec 
                r  = SQRT(DOT_PRODUCT(ra,ra))
                my_val = sfact / r
             END IF
             EwM(idim) = my_val - sfact * g_ewald
          END DO
       END DO
       idim = 0
       DO iparticle2 = 1, SIZE(particle_set)
          ip2 = (iparticle2-1)*SIZE(radii)
          idimo = (iparticle2-1)
          idimo = idimo *(idimo+1) / 2
          DO igauss2 = 1, SIZE(radii)
             idim2 = ip2 + igauss2
             rc2   = radii(igauss2)
             rc22  = rc2*rc2
             DO iparticle1 = 1, iparticle2
                ip1 = (iparticle1-1)*SIZE(radii)
                idim = idimo + iparticle1
                istart_g = 1                
                IF (iparticle1==iparticle2) istart_g = igauss2
                DO igauss1 = istart_g, SIZE(radii)
                   idim1  = ip1 + igauss1
                   rc1    = radii(igauss1)
                   rc12   = rc1*rc1
                   M(idim1,idim2) = EwM(idim)  - sfact * ew_neut - sfact * fac3 * (rc12+rc22)
                   M(idim2,idim1) = M(idim1,idim2)
                END DO
             END DO
          END DO
       END DO
       DEALLOCATE(EwM, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)       
    END IF
    CALL timestop(handle)
  END SUBROUTINE ewald_ddapc_pot

! *****************************************************************************
!> \brief Evaluates the electrostatic potential due to a simple solvation model
!>      Spherical cavity in a dieletric medium
!> \par History
!>      08.2006 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE solvation_ddapc_pot( solvation_section, particle_set, M,  radii, error)
    TYPE(section_vals_type), POINTER         :: solvation_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: M
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'solvation_ddapc_pot', &
      routineP = moduleN//':'//routineN

    INTEGER :: i, idim, idim1, idim2, igauss1, igauss2, ip1, ip2, iparticle1, &
      iparticle2, istart_g, j, l, lmax, n_rep1, n_rep2, ndim, output_unit, &
      stat, weight
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure, fixed_center, ionode
    REAL(KIND=dp)                            :: center(3), eps_in, eps_out, &
                                                factor, mass, mycos, r1, r2, &
                                                Rs, rvec(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: pos, R0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cost, LocP
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env

    failure = .FALSE.
    fixed_center = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    output_unit = -1
    IF (ionode) THEN
      output_unit= cp_logger_get_default_unit_nr(logger)
    END IF
    IF (.NOT.failure) THEN
       ndim   = SIZE(particle_set) * SIZE(radii)
       ALLOCATE(M(ndim, ndim), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
       M      = 0.0_dp
       eps_in = 1.0_dp
       CALL section_vals_val_get(solvation_section,"EPS_OUT",r_val=eps_out,error=error)
       CALL section_vals_val_get(solvation_section,"LMAX",i_val=lmax,error=error)
       CALL section_vals_val_get(solvation_section,"SPHERE%RADIUS",r_val=Rs,error=error)
       CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%XYZ",n_rep_val=n_rep1,&
            error=error)
       IF (n_rep1/=0) THEN
          CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%XYZ",r_vals=R0,&
               error=error)
          center = R0
       ELSE
          CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%ATOM_LIST",&
               n_rep_val=n_rep2,error=error)
          IF (n_rep2/=0) THEN
             CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%ATOM_LIST",i_vals=list,&
                  error=error)
             CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%WEIGHT_TYPE",i_val=weight,&
                  error=error)
             ALLOCATE(R0(SIZE(list)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             SELECT CASE(weight)
             CASE (weight_type_unit)
                R0 = 0.0_dp
                DO i = 1, SIZE(list)
                   R0 = R0 + particle_set(list(i))%r
                END DO
                R0 = R0 / REAL(SIZE(list),KIND=dp)
             CASE (weight_type_mass)
                R0 = 0.0_dp
                mass = 0.0_dp
                DO i = 1, SIZE(list)
                   R0 = R0 + particle_set(list(i))%r * particle_set(list(i))%atomic_kind%mass
                   mass = mass + particle_set(list(i))%atomic_kind%mass
                END DO
                R0 = R0 / mass
             END SELECT
             center = R0
             CALL section_vals_val_get(solvation_section,"SPHERE%CENTER%FIXED",l_val=fixed_center,&
                  error=error)
             IF (fixed_center) THEN             
                CALL section_vals_val_set(solvation_section,"SPHERE%CENTER%XYZ",&
                     r_vals_ptr=R0,error=error)
             ELSE
                DEALLOCATE(R0, stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure) 
             END IF
          END IF
       END IF
       CPPostcondition(n_rep1/=0.OR.n_rep2/=0,cp_failure_level,routineP,error,failure)
       ! Potential calculation
       ALLOCATE(LocP(0:lmax,SIZE(particle_set)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(pos(SIZE(particle_set)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(cost(SIZE(particle_set),SIZE(particle_set)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Determining the single atomic contribution to the dielectric dipole
       DO i = 1, SIZE(particle_set)
          rvec = particle_set(i)%r-center
          r2   = DOT_PRODUCT(rvec,rvec)
          r1   = SQRT(r2)
          IF (r1>=Rs) THEN
             IF (output_unit>0) THEN
                WRITE(output_unit,'(A,I6,A)')"Atom number :: ",i," is out of the solvation sphere"
                WRITE(output_unit,'(2(A,F12.6))')"Distance from the center::",r1," Radius of the sphere::",rs
             END IF
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
          LocP(:,i) = 0.0_dp
          IF (r1/=0.0_dp) THEN
             DO l = 0, lmax
                LocP(l,i) = ( r1**l * REAL(l+1,KIND=dp) * (eps_in-eps_out))/&
                     (Rs**(2*l+1)*eps_in*(REAL(l,KIND=dp)*eps_in+REAL(l+1,KIND=dp)*eps_out))
             END DO
          ELSE
             ! limit for r->0
             LocP(0,i) = (eps_in-eps_out)/(Rs*eps_in*eps_out)
          END IF
          pos(i) =  r1
       END DO
       ! Particle-Particle potential energy matrix
       cost = 0.0_dp
       DO i = 1, SIZE(particle_set)
          DO j = 1, i
             factor    = 0.0_dp
             IF (pos(i)*pos(j)/=0.0_dp) THEN
                mycos     = DOT_PRODUCT(particle_set(i)%r-center,particle_set(j)%r-center)/(pos(i)*pos(j))
                IF (ABS(mycos) > 1.0_dp) mycos = SIGN(1.0_dp, mycos)
                DO l = 0, lmax
                   factor = factor + LocP(l,i) * pos(j)**l * legendre(mycos,l,0)
                END DO
             ELSE
                factor = LocP(0,i)
             END IF
             cost(i,j) = factor
             cost(j,i) = factor
          END DO
       END DO
       ! Computes the full potential energy matrix
       idim = 0
       DO iparticle2 = 1, SIZE(particle_set)
          ip2  = (iparticle2-1)*SIZE(radii)
          DO igauss2 = 1, SIZE(radii)
             idim2 = ip2 + igauss2
             DO iparticle1 = 1, iparticle2
                ip1 = (iparticle1-1)*SIZE(radii)
                istart_g = 1
                IF (iparticle1==iparticle2) istart_g = igauss2
                DO igauss1 = istart_g, SIZE(radii)
                   idim1  = ip1 + igauss1
                   M(idim1,idim2) = cost(iparticle1,iparticle2)
                   M(idim2,idim1) = M(idim1,idim2)
                END DO
             END DO
          END DO
       END DO
       DEALLOCATE(cost,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(pos,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(LocP,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE solvation_ddapc_pot

END MODULE cp_ddapc_methods
