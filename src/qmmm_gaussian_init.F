!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qmmm_gaussian_init [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_init
!!
!!   FUNCTION
!!     Initialize the use of the gaussians to treat the QMMM 
!!     coupling potential
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     6.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_gaussian_init
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             dbl
  USE memory_utilities,                ONLY: reallocate
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type,&
                                             LongRange_type
  USE qmmm_gaussian_input,             ONLY: read_mm_potential
  USE global_types,                    ONLY: global_environment_type
  USE simulation_cell,                 ONLY: cell_type
  USE qs_util,                         ONLY: exp_radius
  USE mathconstants,                   ONLY: pi
  USE pw_pool_types, ONLY: pw_pool_type, pw_pool_p_type
  USE pw_env_types, ONLY: pw_env_type, pw_env_release, pw_env_get
  USE qmmm_types, only: qmmm_env_type
  USE gaussian_gridlevels, ONLY: gaussian_gridlevel, gridlevel_info_type
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_gaussian_init'
  PUBLIC :: qmmm_gaussian_initialize

CONTAINS

!!****f* qmmm_gaussian_init/qmmm_gaussian_initialize [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_initialize
!!
!!   FUNCTION
!!     Initialize the Gaussian QMMM Environment
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     06.2004 created [tlaino]
!!
!!*** **********************************************************************  
  SUBROUTINE qmmm_gaussian_initialize(qmmm_env,globenv,pw_env,&
                                      qm_cell_small,mm_cell,error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_env_type), pointer :: qmmm_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(pw_env_type), pointer :: pw_env
    TYPE(cell_type), POINTER                 :: qm_cell_small, mm_cell
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_gaussian_initialize', &
      routineP = moduleN//':'//routineN
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER                       :: qmmm_gaussian_fns    
    REAL(KIND=dbl), DIMENSION(:), POINTER            :: mm_el_pot_radius 
    CHARACTER(len=default_string_length) :: qmmm_coupl_type  
    REAL(KIND=dp)                 :: eps_mm_rspace

    REAL(KIND=dbl), DIMENSION(:), POINTER :: radius, maxradius
    INTEGER :: stat, i, j, Ndim ,ilevel
    LOGICAL :: Failure, Found
    TYPE(qmmm_gaussian_type), POINTER :: mypgf
    REAL(KIND=dbl) :: alpha, mymaxradius, Prefactor
    TYPE(pw_pool_p_type), DIMENSION(:), pointer :: pools
    TYPE(gridlevel_info_type), POINTER :: gridlevel_info
    ! Statements
    NULLIFY(radius,mypgf, maxradius, gridlevel_info)
    mm_el_pot_radius => qmmm_env%mm_el_pot_radius 
    qmmm_coupl_type  = qmmm_env%qmmm_coupl_type  
    eps_mm_rspace = qmmm_env%eps_mm_rspace
    SELECT CASE(qmmm_coupl_type)
    CASE("GAUSS")
       !
       ! Preprocessing...
       !
       ALLOCATE(radius(1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       Loop_on_all_values: DO I = 1, SIZE(mm_el_pot_radius)
          Found=.FALSE.
          Loop_on_found_values: DO J = 1, SIZE(radius)-1
             IF (mm_el_pot_radius(i).EQ.radius(j)) THEN
                Found=.TRUE.
                EXIT Loop_on_found_values
             END IF
          END DO Loop_on_found_values
          IF (.NOT.Found) THEN
             Ndim = SIZE(radius)
             radius(Ndim) = mm_el_pot_radius(i)
             Ndim = Ndim + 1
             CALL REALLOCATE(radius,1,Ndim)
          END IF
       END DO Loop_on_all_values
       CALL REALLOCATE(radius,1,Ndim-1)
       !
       ALLOCATE(qmmm_gaussian_fns(SIZE(radius)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)       
       DO I = 1, SIZE(radius)
          NULLIFY(qmmm_gaussian_fns(I)%pgf)
          ALLOCATE(qmmm_gaussian_fns(I)%pgf, stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          NULLIFY(qmmm_gaussian_fns(I)%pgf%mm_atom_index)
          NULLIFY(qmmm_gaussian_fns(I)%pgf%Ak)
          NULLIFY(qmmm_gaussian_fns(I)%pgf%Gk)
          NULLIFY(qmmm_gaussian_fns(I)%pgf%grid_level)
          ! 
          ! Default Values
          !
          qmmm_gaussian_fns(I)%pgf%Elp_Radius          = radius(I)

          ALLOCATE(qmmm_gaussian_fns(I)%pgf%mm_atom_index(1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          ! Build mm_atom_index List
          DO J= 1, SIZE(mm_el_pot_radius)
             IF (radius(I).EQ.mm_el_pot_radius(J)) THEN
                Ndim = SIZE(qmmm_gaussian_fns(I)%pgf%mm_atom_index) 
                qmmm_gaussian_fns(I)%pgf%mm_atom_index(Ndim) = J
                CALL reallocate(qmmm_gaussian_fns(I)%pgf%mm_atom_index,1, Ndim+1)
             ENDIF
          END DO
          CALL reallocate(qmmm_gaussian_fns(I)%pgf%mm_atom_index, 1, Ndim)          
       END DO
       DEALLOCATE(radius,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)      
       !
       CALL read_mm_potential(globenv, qmmm_gaussian_fns, qm_cell_small,&
                              error)
       !
       CALL pw_env_get(pw_env,pw_pools=pools, gridlevel_info=gridlevel_info,&
            error=error)
       ALLOCATE(maxradius(size(pools)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       maxradius = 0.0_dbl
       DO J = 1, SIZE(qmmm_gaussian_fns)
          mypgf => qmmm_gaussian_fns(J)%pgf
          ALLOCATE(mypgf%grid_level(SIZE(mypgf%Ak)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          mypgf%grid_level=0
          mymaxradius = 0.0_dbl
          DO I = 1, mypgf%number_of_gaussians
             IF ( mypgf%Gk(I) /= 0.0_dbl) THEN
                alpha = 1.0_dbl / mypgf%Gk(I)
                alpha = alpha*alpha
                ilevel=gaussian_gridlevel(gridlevel_info,alpha)
                Prefactor = mypgf%Ak(I)
                mymaxradius=exp_radius( 0, alpha, eps_mm_rspace, Prefactor)
                maxradius(ilevel)=MAX(maxradius(ilevel),mymaxradius)
                mypgf%grid_level(i)=ilevel
             END IF
          END DO
       END DO       
    CASE DEFAULT
       WRITE(*,'(A)')" QMMM Gaussian Data Not Initialized!"
    END SELECT

    CALL qmmm_gaussian_eval_LR(pgfs=qmmm_gaussian_fns,&
                               mm_cell=mm_cell,&
                               error=error)

    qmmm_env%pgfs => qmmm_gaussian_fns
    qmmm_env%maxradius => maxradius

  END SUBROUTINE qmmm_gaussian_initialize

!!****f* qmmm_gaussian_init/qmmm_gaussian_eval_LR [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_eval_LR
!!
!!   FUNCTION
!!     Initialize the vectors that store the LR information
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     08.2004 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE qmmm_gaussian_eval_LR(pgfs, mm_cell, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER                       :: pgfs
    TYPE(cell_type), POINTER                         :: mm_cell   
     TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_gaussian_eval_LR', &
      routineP = moduleN//':'//routineN
    LOGICAL :: failure
    TYPE(qmmm_gaussian_type), POINTER                :: pgf
    TYPE(LongRange_Type), POINTER :: myLR
    REAL(KIND=dbl) :: Rmax, r, r2
    REAL(KIND=dbl) :: r0, r0Nexp, r0Nexp1, Gterm, Gdpot
    REAL(KIND=dbl) :: alpha, height, W, dpot, rNexp, rNexp1, rNexp0
    INTEGER        :: I, J, Np, stat, IGauss
    REAL(KIND=dbl), PARAMETER :: dx=0.00001_dbl


    Rmax = SQRT(mm_cell%hmat(1,1)**2+&
                mm_cell%hmat(2,2)**2+&
                mm_cell%hmat(3,3)**2 ) 
    NULLIFY(myLR, pgf)
    DO I = 1, SIZE(pgfs)
       pgf => pgfs(I)%pgf
       myLR  => pgfs(I)%pgf%LR
       Np    =  NINT(Rmax / dx)
       myLR%dx   = dx
       myLR%Rmin = dx
       myLR%Rmax = Rmax
       myLR%Npts = Np
       !
       NULLIFY(myLR%Pot, myLR%DPot)
       ALLOCATE(myLR%Pot(Np),  stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       ALLOCATE(myLR%DPot(Np), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       r0      = pgf%Elp_Radius
       r0Nexp  = r0**2
       r0Nexp  = r0Nexp*r0Nexp
       r0Nexp1 = r0*r0Nexp       
       DO J = 1, Np
          r      = REAL(J,dbl)*dx
          r2     = r * r
          rNexp  = r**2
          rNexp0 = rNexp * r 
          rNexp  = rNexp*rNexp
          rNexp1 = rNexp*r
          Gterm  = 0.0_dbl
          GdPot  = 0.0_dbl
          DO IGauss = 1, pgf%Number_of_Gaussians
             alpha  = 1.0_dbl / pgf%Gk(IGauss)
             alpha  = alpha * alpha
             height = pgf%Ak(IGauss)
             W      = height
             Gterm  = Gterm + W * EXP(-alpha*r2)
             GdPot  = GdPot - 2.0_dbl * alpha * W * r * EXP(-alpha*r2)
          END DO
          Gterm = Gterm + pgf%A0 
          Gterm = Gterm - 1.0_dbl * (rNexp - r0Nexp)/(rNexp1 - r0Nexp1) 
          dpot = 4.0_dbl * (rNexp1-r0Nexp1)* rNexp0 - &
                 5.0_dbl * (rNexp -r0Nexp) * rNexp
          dpot  = GdPot - dpot / ((rNexp1 - r0Nexp1)**2 ) ! The Real derivative
          dpot  = dpot / r
          !
          myLR%Pot(J)  = Gterm
          myLR%DPot(J) = Dpot
       END DO

    END DO
  END SUBROUTINE qmmm_gaussian_eval_LR

END MODULE qmmm_gaussian_init
