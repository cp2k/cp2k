!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

MODULE qmmm_elpot
!!****h* cp2k/qmmm_elpot *
!!
!!   NAME
!!    qmmm_elpot
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     09.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************  
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
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_types,                      ONLY: qmmm_env_type,&
                                             qmmm_Pot_Type
  USE simulation_cell,                 ONLY: cell_type  
  USE memory_utilities,                ONLY: reallocate
  USE global_types,                    ONLY: global_environment_type
  USE erf_fn,                          ONLY: erf
  USE mathconstants,                   ONLY: rootpi
  USE pw_env_types, ONLY: pw_env_type, pw_env_release, pw_env_get
  USE input_cp2k, ONLY: do_qmmm_grid, do_qmmm_gauss, do_qmmm_none

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_elpot'
  PUBLIC :: qmmm_potential_init
  
CONTAINS
!!****f* qmmm_elpot/qmmm_potential_init [1.0] *
!!
!!   NAME
!!     qmmm_potential_init
!!
!!   FUNCTION
!!     Initialize the QMMM potential stored on vector, 
!!     according the qmmm_coupl_type
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
!!     09.2004 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE qmmm_potential_init(qmmm_env, qm_cell_small, mm_cell,&
                                 globenv, pw_env, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_env_type), POINTER                     :: qmmm_env
    TYPE(cell_type), POINTER                         :: mm_cell, qm_cell_small   
     TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(pw_env_type), POINTER :: pw_env    
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_potential_init', &
      routineP = moduleN//':'//routineN
    LOGICAL :: failure, found
    TYPE(qmmm_gaussian_type), POINTER                :: pgf
    REAL(KIND=dp) :: Rmax, Rmin
    REAL(KIND=dp) :: rc, x, t, A,G, t1, t2
    INTEGER        :: I, J, K, Np, stat, ndim, ig
    INTEGER,        DIMENSION(:), POINTER :: mm_atom_index
    REAL(KIND=dp), PARAMETER :: dx = 0.05_dp, d2x=dx*dx
    REAL(KIND=dp), DIMENSION(:,:),POINTER :: pot0_2
    REAL(KIND=dp), DIMENSION(:), POINTER :: radius
    INTEGER :: qmmm_coupl_type

    Rmin = 0.0_dp
    Rmax = SQRT(mm_cell%hmat(1,1)**2+&
                mm_cell%hmat(2,2)**2+&
                mm_cell%hmat(3,3)**2 ) 
    np   = CEILING(rmax/dx)+1
    qmmm_coupl_type  = qmmm_env%qmmm_coupl_type
    !
    ! Preprocessing
    !
    ALLOCATE(radius(1),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    radius(1) = qmmm_env%mm_el_pot_radius(1)
    Loop_on_all_values: DO I = 2, SIZE(qmmm_env%mm_el_pot_radius)
       Found=.FALSE.
       Loop_on_found_values: DO J = 1, SIZE(radius)
          IF (qmmm_env%mm_el_pot_radius(i).EQ.radius(j)) THEN
             Found=.TRUE.
             EXIT Loop_on_found_values
          END IF
       END DO Loop_on_found_values
       IF (.NOT.Found) THEN
          Ndim = SIZE(radius)
          Ndim = Ndim + 1
          CALL REALLOCATE(radius,1,Ndim)
          radius(Ndim) = qmmm_env%mm_el_pot_radius(i)          
       END IF
    END DO Loop_on_all_values
    !
    CPPrecondition(.NOT.ASSOCIATED(qmmm_env%potentials),cp_failure_level,routineP,error,failure)
    ALLOCATE(qmmm_env%potentials(SIZE(radius)),stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

    Potential_Type: DO K = 1, SIZE(radius)
       
       rc = radius(K)
       ALLOCATE(qmmm_env%Potentials(K)%Pot, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_grid)
          ALLOCATE(pot0_2(3,np), stat=stat)
       CASE(do_qmmm_gauss)
          ALLOCATE(pot0_2(2,np), stat=stat)
       END SELECT
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_grid)
          pot0_2(1,1)= 2.0_dp/(rootpi*rc)
          pot0_2(2,1)= 0.0_dp
          pot0_2(3,1)=-4.0_dp/(3.0_dp*rootpi*rc**3)
          x=0.0_dp
          DO i=2,np
             x=x+dx
             pot0_2(1,i)=erf(x/rc)/x
             t=2._dp/(rootpi*x*rc)*EXP(-(x/rc)**2)
             pot0_2(2,i)=(t-pot0_2(1,i)/x)*dx
             pot0_2(3,i)=(-t*2._dp*(x/rc**2+1._dp/x)+2._dp*pot0_2(1,i)/x**2)*d2x
          END DO
       CASE(do_qmmm_gauss)
          pot0_2(1,1)= 2.0_dp/(rootpi*rc)
          pot0_2(2,1)= 0.0_dp
          x=0.0_dp
          DO i=2,np
             x=x+dx
             pot0_2(1,i)=erf(x/rc)/x
             t=2._dp/(rootpi*x*rc)*EXP(-(x/rc)**2)
             pot0_2(2,i)=(t-pot0_2(1,i)/x)*dx
          END DO
          pgf => qmmm_env%pgfs(K)%pgf
          CPPostcondition(pgf%Elp_Radius==rc,cp_failure_level,routineP,error,failure)
          DO Ig = 1, pgf%number_of_gaussians
             A = pgf%Ak(Ig) 
             G = pgf%Gk(Ig)
             pot0_2(1,1) = pot0_2(1,1) - A
             x = 0.0_dp
             DO i=2,np
                x=x+dx
                t=EXP(-(x/G)**2)*A
                t1=1/G**2
                t2=t1*t
                pot0_2(1,i) = pot0_2(1,i)-t 
                pot0_2(2,i) = pot0_2(2,i)+2.0_dp*x*t2*dx
             END DO
          END DO
       CASE DEFAULT
          DEALLOCATE(qmmm_env%Potentials(K)%Pot)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          WRITE(*,'(A)')" QMMM Potential - Spline Interpolation - not Initialized!"
          EXIT Potential_Type
       END SELECT

       NULLIFY(mm_atom_index)
       ALLOCATE(mm_atom_index(1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       ! Build mm_atom_index List
       DO J= 1, SIZE(qmmm_env%mm_el_pot_radius)
          IF (rc.EQ.qmmm_env%mm_el_pot_radius(J)) THEN
             Ndim = SIZE(mm_atom_index) 
             mm_atom_index(Ndim) = J
             CALL reallocate(mm_atom_index,1, Ndim+1)
          ENDIF
       END DO
       CALL reallocate(mm_atom_index, 1, Ndim)
      
       NULLIFY(qmmm_env%Potentials(K)%Pot%pot0_2)
       CALL qmmm_pot_type_create(qmmm_env%Potentials(K)%Pot, pot0_2=pot0_2,&
                                 Rmax=Rmax, Rmin=Rmin, dx=dx, Rc=rc, npts=np,&
                                 mm_atom_index=mm_atom_index)

    END DO Potential_Type

    DEALLOCATE(radius, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
  END SUBROUTINE qmmm_potential_init

!!****f* qmmm_elpot/qmmm_pot_type_create [1.0] *
!!
!!   NAME
!!     qmmm_pot_type_create
!!
!!   FUNCTION
!!     Creates the qmmm_pot_type structure
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
!!     09.2004 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE qmmm_pot_type_create(Pot, pot0_2, Rmax, Rmin, dx, npts, rc,&
       mm_atom_index)
    IMPLICIT NONE
    TYPE(qmmm_Pot_Type), POINTER :: Pot
    REAL(KIND=dp), DIMENSION(:,:), POINTER :: pot0_2
    REAL(KIND=dp), INTENT(IN) :: Rc, Rmax, Rmin, dx
    INTEGER, INTENT(IN) :: npts
    INTEGER, DIMENSION(:), POINTER :: mm_atom_index
    
    Pot%pot0_2        => pot0_2
    Pot%mm_atom_index => mm_atom_index
    Pot%Rmax          =  Rmax
    Pot%Rmin          =  Rmin
    Pot%Rc            =  rc
    Pot%dx            =  dx
    Pot%npts          =  npts

  END SUBROUTINE qmmm_pot_type_create

END MODULE qmmm_elpot
