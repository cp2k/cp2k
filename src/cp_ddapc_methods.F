!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/cp_ddapc_methods [1.0] *
!!
!!   NAME
!!     cp_ddapc_methods
!!
!!   FUNCTION
!!     contains information regarding the decoupling/recoupling method of Bloechl
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!
!!
!!   SOURCE
!*****************************************************************************
MODULE cp_ddapc_methods
  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE erf_fn,                          ONLY: erfc
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE mathlib,                         ONLY: invert_matrix
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_spline_utils,                 ONLY: Eval_Interp_Spl3_pbc
  USE pw_types,                        ONLY: pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_ddapc_methods'  
  PUBLIC :: ddapc_eval_gfunc,&
            build_b_vector,&
            build_der_b_vector,&
            build_A_matrix,&
            build_der_A_matrix,&
            ddapc_eval_AmI,&
            ewald_ddapc_pot

CONTAINS

!!****f* cp_ddapc_methods/ddapc_eval_gfunc [1.0] *
!!
!!   NAME
!!     ddapc_eval_gfunc
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     
!!
!!*** **********************************************************************
  SUBROUTINE ddapc_eval_gfunc(gfunc, w, gcut, rho_tot_g, radii, error)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: gfunc
    REAL(kind=dp), DIMENSION(:), POINTER     :: w
    REAL(KIND=dp), INTENT(IN)                :: gcut
    TYPE(pw_type), POINTER                   :: rho_tot_g
    REAL(kind=dp), DIMENSION(:), POINTER     :: radii
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_eval_gfunc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ig, igauss, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: g2, gcut2, rc, rc2

    CALL timeset(routineN,"I","",handle)
    gcut2 = gcut * gcut
    failure = .FALSE.
    IF (.NOT.failure) THEN
       !   
       ALLOCATE(gfunc(rho_tot_g % pw_grid % first_gne0:rho_tot_g % pw_grid % ngpts_cut_local,&
            SIZE(radii)), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(w(rho_tot_g % pw_grid % first_gne0:rho_tot_g % pw_grid % ngpts_cut_local), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       gfunc = 0.0_dp
       w     = 0.0_dp
       DO igauss = 1, SIZE(radii)
          rc  = radii(igauss)
          rc2 = rc*rc
          DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
             g2 = rho_tot_g % pw_grid % gsq ( ig )
             IF (g2 > gcut2) EXIT
             gfunc(ig,igauss) = EXP(-g2*rc2/4.0_dp)
          END DO
       END DO
       DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
          g2 = rho_tot_g % pw_grid % gsq ( ig )
          IF (g2 > gcut2) EXIT
          w(ig) = 4.0_dp * pi * (g2 -gcut2)**2/(g2*gcut2)
       END DO
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE ddapc_eval_gfunc

!!****f* cp_ddapc/build_b_vector [1.0] *
!!
!!   NAME
!!     build_b_vector
!!
!!   FUNCTION
!!     Computes the B vector for the solution of the linear system
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_b_vector', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: phase
    INTEGER                                  :: handle, idim, ig, igauss, &
                                                igmax, iparticle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, g2, gcut2, gvec(3), &
                                                rvec(3), tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_bv

    failure = .FALSE.
    NULLIFY(my_bv)
    CALL timeset(routineN,"I","",handle)
    gcut2 = gcut * gcut
    IF (.NOT.failure) THEN
       ALLOCATE(my_bv(rho_tot_g % pw_grid % first_gne0:rho_tot_g % pw_grid % ngpts_cut_local),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !   
       idim  = 0
       DO iparticle = 1, SIZE(particle_set)
          rvec =  particle_set(iparticle)%r
          my_bv = 0.0_dp
          igmax = 0
          DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
             g2 = rho_tot_g % pw_grid % gsq ( ig )
             IF (g2 > gcut2) EXIT
             igmax = ig
             gvec  = rho_tot_g % pw_grid % g (:,ig)
             arg   = DOT_PRODUCT(gvec,rvec)
             phase = CMPLX(COS(arg),-SIN(arg),KIND=dp)
             my_bv (ig) =  w(ig) * REAL(CONJG(rho_tot_g%cc(ig))*phase,KIND=dp)
          END DO
          DO igauss = 1, SIZE(radii)
             idim = (iparticle-1)*SIZE(radii)+igauss
             tmp = 0.0_dp
             DO ig = rho_tot_g % pw_grid % first_gne0, igmax
                tmp = tmp + my_bv(ig)* gfunc(ig,igauss)
             END DO
             bv (idim) = tmp
          END DO
       END DO
       DEALLOCATE(my_bv,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE build_b_vector
    
!!****f* cp_ddapc/build_A_matrix [1.0] *
!!
!!   NAME
!!     build_A_matrix
!!
!!   FUNCTION
!!     Computes the A matrix for the solution of the linear system
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_A_matrix', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: phase1, phase2
    INTEGER                                  :: handle, idim1, idim2, ig, &
                                                igauss1, igauss2, igmax, &
                                                iparticle1, iparticle2, &
                                                istart_g, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, g2, gcut2, gvec(3), &
                                                rvec1(3), rvec2(3), tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: my_Am

    failure = .FALSE.
    NULLIFY(my_Am)
    CALL timeset(routineN,"I","",handle)
    gcut2 = gcut * gcut
    IF (.NOT. failure) THEN
       ALLOCATE(my_Am(rho_tot_g % pw_grid % first_gne0:rho_tot_g % pw_grid % ngpts_cut_local),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iparticle1 = 1, SIZE(particle_set)
          rvec1=particle_set(iparticle1)%r
          DO iparticle2 = iparticle1, SIZE(particle_set)
             rvec2=particle_set(iparticle2)%r
             my_Am = 0.0_dp
             igmax = 0
             DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
                g2 = rho_tot_g % pw_grid % gsq ( ig )
                IF (g2 > gcut2) EXIT
                igmax  = ig
                gvec   = rho_tot_g % pw_grid % g (:,ig)
                arg    = DOT_PRODUCT(gvec,rvec1)
                phase1 = CMPLX(COS(arg),-SIN(arg),KIND=dp)                     
                arg    = DOT_PRODUCT(gvec,rvec2)
                phase2 = CMPLX(COS(arg),-SIN(arg),KIND=dp)
                my_Am (ig) = w(ig) * REAL(CONJG(phase1)*phase2,KIND=dp)
             END DO
             DO igauss1 = 1, SIZE(radii)
                idim1 = (iparticle1-1)*SIZE(radii)+igauss1
                istart_g = 1
                IF (iparticle2==iparticle1) istart_g = igauss1
                DO igauss2 = istart_g, SIZE(radii)
                   idim2 = (iparticle2-1)*SIZE(radii)+igauss2
                   tmp = 0.0_dp
                   DO ig = rho_tot_g % pw_grid % first_gne0, igmax
                      tmp = tmp + my_Am(ig)*gfunc(ig,igauss1)*gfunc(ig,igauss2)
                   END DO
                   Am (idim2,idim1) = tmp
                   Am (idim1,idim2) = tmp
                END DO
             END DO
          END DO
       END DO
       DEALLOCATE(my_Am,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE build_A_matrix

!!****f* cp_ddapc/build_der_b_vector [1.0] *
!!
!!   NAME
!!     build_der_b_vector
!!
!!   FUNCTION
!!     Computes the derivative of B vector for the evaluation of the Pulay forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_der_b_vector', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: dphase
    INTEGER                                  :: handle, idim, ig, igauss, &
                                                igmax, iparticle, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg, g2, gcut2
    REAL(KIND=dp), DIMENSION(3)              :: gvec, rvec, tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_dbv

    failure = .FALSE.
    NULLIFY(my_dbv)
    CALL timeset(routineN,"I","",handle)
    gcut2 = gcut * gcut
    IF (.NOT.failure) THEN
       ALLOCATE(my_dbv(rho_tot_g % pw_grid % first_gne0:rho_tot_g % pw_grid % ngpts_cut_local,3),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       idim  = 0
       DO iparticle = 1, SIZE(particle_set)
          rvec =  particle_set(iparticle)%r
          IF (iparticle /= iparticle0) CYCLE
          my_dbv = 0.0_dp
          igmax  = 0
          DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
             g2 = rho_tot_g % pw_grid % gsq ( ig )
             IF (g2 > gcut2) EXIT
             igmax  = ig
             gvec   = rho_tot_g % pw_grid % g (:,ig)
             arg    = DOT_PRODUCT(gvec,rvec)
             dphase = - CMPLX(SIN(arg),COS(arg),KIND=dp) 
             my_dbv (ig,:) = w(ig) * REAL(CONJG(rho_tot_g%cc(ig))*dphase,KIND=dp) * gvec(:)
          END DO
          DO igauss = 1, SIZE(radii)
             idim = (iparticle-1)*SIZE(radii)+igauss
             tmp = 0.0_dp
             DO ig = rho_tot_g % pw_grid % first_gne0, igmax
                tmp = tmp + my_dbv(ig,1:3)*gfunc(ig,igauss)
             END DO
             dbv (idim,1:3) = tmp
          END DO
       END DO
       DEALLOCATE(my_dbv,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE build_der_b_vector
    
!!****f* cp_ddapc/build_der_A_matrix [1.0] *
!!
!!   NAME
!!     build_der_A_matrix
!!
!!   FUNCTION
!!     Computes the derivative of the A matrix for the evaluation of the 
!!     Pulay forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_der_A_matrix', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: dphase1, dphase2, phase1, &
                                                phase2
    INTEGER                                  :: handle, idim1, idim2, ig, &
                                                igauss1, igauss2, igmax, &
                                                iparticle1, iparticle2, &
                                                istart_g, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: arg1, arg2, g2, gcut2
    REAL(KIND=dp), DIMENSION(3)              :: gvec, rvec1, rvec2, tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: my_Am

    NULLIFY(my_Am)
    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
    gcut2 = gcut * gcut
    IF (.NOT. failure) THEN
       ALLOCATE(my_Am(rho_tot_g % pw_grid % first_gne0:rho_tot_g % pw_grid % ngpts_cut_local,3),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO iparticle1 = 1, SIZE(particle_set)
          rvec1 =  particle_set(iparticle1)%r
          DO iparticle2 = iparticle1, SIZE(particle_set)
             rvec2 =  particle_set(iparticle2)%r
             IF ((iparticle1 /= iparticle0).AND.(iparticle2 /= iparticle0)) CYCLE
             my_Am = 0.0_dp
             igmax = 0
             DO ig = rho_tot_g % pw_grid % first_gne0, rho_tot_g % pw_grid % ngpts_cut_local
                g2 = rho_tot_g % pw_grid % gsq ( ig )
                IF (g2 > gcut2) EXIT
                igmax  = ig
                gvec   = rho_tot_g % pw_grid % g (:,ig)
                arg1   = DOT_PRODUCT(gvec,rvec1)
                arg2   = DOT_PRODUCT(gvec,rvec2)
                IF    ((iparticle1 == iparticle0).AND.(iparticle2==iparticle0))  THEN
                   dphase1 = - CMPLX(SIN(arg1), COS(arg1),KIND=dp)                         
                   phase2  =   CMPLX(COS(arg2),-SIN(arg2),KIND=dp)
                   my_Am (ig,1:3) = 2.0_dp * REAL(CONJG(dphase1)*phase2,KIND=dp) * gvec
                ELSEIF (iparticle1 == iparticle0) THEN
                   dphase1 = - CMPLX(SIN(arg1), COS(arg1),KIND=dp)
                   phase2  =   CMPLX(COS(arg2),-SIN(arg2),KIND=dp)
                   my_Am (ig,1:3) = REAL(CONJG(dphase1)*phase2,KIND=dp) * gvec
                ELSEIF (iparticle2 == iparticle0) THEN
                   phase1  =   CMPLX(COS(arg1),-SIN(arg1),KIND=dp)
                   dphase2 = - CMPLX(SIN(arg2), COS(arg2),KIND=dp)
                   my_Am (ig,1:3) = REAL(CONJG(phase1)*dphase2,KIND=dp) * gvec
                END IF
             END DO
             DO igauss1 = 1, SIZE(radii)
                idim1 = (iparticle1-1)*SIZE(radii)+igauss1
                istart_g = 1
                IF (iparticle2==iparticle1) istart_g = igauss1
                DO igauss2 = istart_g, SIZE(radii)
                   idim2 = (iparticle2-1)*SIZE(radii)+igauss2
                   tmp = 0.0_dp
                   DO ig = rho_tot_g % pw_grid % first_gne0, igmax
                      tmp = tmp + my_Am(ig,1:3)*gfunc(ig,igauss1)*gfunc(ig,igauss2)*w(ig)
                   END DO
                   dAm (idim1,idim2,1:3) = tmp
                   dAm (idim2,idim1,1:3) = tmp
                END DO
             END DO
          END DO
       END DO
       DEALLOCATE(my_Am,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE build_der_A_matrix

!!****f* cp_ddapc/ddapc_eval_AmI [1.0] *
!!
!!   NAME
!!     ddapc_eval_AmI
!!
!!   FUNCTION
!!     Computes the inverse AmI of the Am matrix 
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [tlaino]
!!
!!*** **********************************************************************
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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ddapc_eval_AmI', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ndim, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: inv_error
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cv
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Am, AmI

    failure = .FALSE.
    CALL timeset(routineN,"I","",handle)
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
       CALL build_A_matrix(Am, gfunc, w, particle_set, radii, rho_tot_g, gcut, error);  Am = Am /(Vol*Vol)
       CALL mp_sum(Am,rho_tot_g%pw_grid%para%group)
       CALL invert_matrix(Am, AmI, inv_error, "N")
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
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE ddapc_eval_AmI

!!****f* cp_ddapc_methods/ewald_ddapc_pot [1.0] *
!!
!!   NAME
!!     ewald_ddapc_pot
!!
!!   FUNCTION
!!     Evaluates the Ewald term E2 and E3 energy term for the decoupling/coupling
!!     of periodic images
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2005 created [tlaino]
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE ewald_ddapc_pot(coeff, qmmm_periodic, cell, multipole_section,&
       particle_set, M, radii, error)
    TYPE(pw_type), POINTER                   :: coeff
    LOGICAL, INTENT(IN)                      :: qmmm_periodic
    TYPE(cell_type), POINTER                 :: cell
    TYPE(section_vals_type), POINTER         :: multipole_section
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: M
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radii
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ewald_ddapc_pot', &
      routineP = moduleN//':'//routineN

    INTEGER :: ewmdim, handle, idim, idim1, idim2, idimo, igauss1, igauss2, &
      ip1, ip2, iparticle1, iparticle2, istart_g, k1, k2, k3, n_rep, ndim, &
      nmax1, nmax2, nmax3, r1, r2, r3, rmax1, rmax2, rmax3, sfact, stat
    LOGICAL                                  :: analyt, failure
    REAL(KIND=dp) :: alpha, dfact, eps, eps0, ew_neut, fac, fac2, fac3, fs, &
      g_ewald, galpha, gsq, gsqi, my_val, r, r_ewald, rc1, rc12, rc2, rc22, &
      rcut, t1, tol, tol1
    REAL(KIND=dp), DIMENSION(3)              :: fvec, gvec, ra, rb, rvec
    REAL(KIND=dp), DIMENSION(:), POINTER     :: EwM
    TYPE(unit_convert_type), POINTER         :: units

    failure = .FALSE.
    NULLIFY( EwM, units)
    CALL timeset(routineN,"I","",handle)
    IF (.NOT.failure) THEN
       CPPostcondition(.NOT.ASSOCIATED(M),cp_failure_level,routineP,error,failure)
       CPPostcondition(ASSOCIATED(radii),cp_failure_level,routineP,error,failure)
       CALL get_cp2k_units ( units )
       sfact =  1.0_dp
       IF (qmmm_periodic) THEN 
          sfact = -1.0_dp
       END IF
       rcut = MIN(cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3))/2.0_dp
       CALL section_vals_val_get(multipole_section,"RCUT",n_rep_val=n_rep,error=error)
       IF (n_rep==1) CALL section_vals_val_get(multipole_section,"RCUT",r_val=rcut,error=error)
       CALL section_vals_val_get(multipole_section,"EWALD_PRECISION",r_val=eps,error=error)
       CALL section_vals_val_get(multipole_section,"ANALYTICAL_GTERM",l_val=analyt,error=error)
       !
       ! Setting-up parameters for Ewald summation
       !
       eps0   = units%eps0
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
       fac    = 1.d0/cell%deth
       fac2   = 1.0_dp / ( 4.0_dp * pi * eps0 )
       fac3   = fac*pi
       fvec   = 2.0_dp * pi / (/cell%hmat(1,1),cell%hmat(2,2),cell%hmat(3,3)/)
       ew_neut = - fac / ( 4.0_dp * eps0 * alpha ** 2 )
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
             ip2 = (iparticle2-1)*SIZE(radii)
             idim = idim + 1
             !
             ! Real-Space Contribution
             !
             my_val = 0.0_dp
             rvec   = particle_set(iparticle1)%r - particle_set(iparticle2)%r
             r_ewald = 0.0_dp
             dfact = 1.0_dp
             IF (iparticle1 == iparticle2) dfact = 0.5_dp
             IF (iparticle1 /= iparticle2) THEN
                ra = rvec 
                r  = SQRT(DOT_PRODUCT(ra,ra))
                IF (r<=rcut) THEN 
                   t1 = erfc(alpha*r) / r
                   r_ewald = t1
                END IF
             END IF
             DO r1 = -rmax1, rmax1
                DO r2 = -rmax2, rmax2
                   DO r3 = -rmax3, rmax3
                      IF ((r1==0).AND.(r2==0).AND.(r3==0)) CYCLE
                      rb = (/REAL(r1,KIND=dp),REAL(r2,KIND=dp),REAL(r3,KIND=dp)/)
                      ra = rvec + MATMUL(cell%hmat,rb) 
                      r  = SQRT(DOT_PRODUCT(ra,ra))
                      IF (r<=rcut) THEN 
                         t1 = erfc(alpha*r) / r
                         r_ewald = r_ewald + t1                           
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
             g_ewald = r_ewald*fac2 + g_ewald/eps0
             !
             ! Self Contribution
             !
             IF (iparticle1 == iparticle2) THEN
                g_ewald = g_ewald - alpha/(2.0_dp*eps0*pi*SQRT(pi))
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
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE ewald_ddapc_pot

END MODULE cp_ddapc_methods
