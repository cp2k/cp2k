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
                                             cp_warning_level,&
                                             cp_logger_type
  USE cp_output_handling, ONLY: cp_should_output
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type,&
                                             qmmm_Pot_Type,&
                                             qmmm_pot_p_type
  USE simulation_cell,                 ONLY: cell_type  
  USE memory_utilities,                ONLY: reallocate
  USE global_types,                    ONLY: global_environment_type
  USE erf_fn,                          ONLY: erf
  USE mathconstants,                   ONLY: rootpi
  USE pw_env_types,                    ONLY: pw_env_type, pw_env_release, pw_env_get
  USE input_constants,                 ONLY: do_qmmm_grid, do_qmmm_gauss, do_qmmm_none
  USE particle_types,                  ONLY: particle_type
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get, section_vals_val_get  

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_elpot'
  PUBLIC :: qmmm_potential_init,&
            qmmm_validate_this_run
  
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
  SUBROUTINE qmmm_potential_init(qmmm_coupl_type, mm_el_pot_radius, potentials,&
    pgfs, qm_cell_small, mm_cell, globenv, pw_env, compatibility, print_section,&
    error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_pot_p_type), DIMENSION(:), POINTER            :: potentials
    TYPE(qmmm_gaussian_p_type), &
         DIMENSION(:), POINTER                       :: pgfs    
    REAL(KIND=dp), DIMENSION(:), POINTER             :: mm_el_pot_radius 
    INTEGER, INTENT(IN)                              :: qmmm_coupl_type      
    TYPE(cell_type), POINTER                         :: mm_cell, qm_cell_small  
    TYPE(section_vals_type), POINTER         :: print_section
     TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(pw_env_type), POINTER :: pw_env    
    LOGICAL, INTENT(IN) :: compatibility
    ! Local Variables
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_potential_init', &
      routineP = moduleN//':'//routineN
    LOGICAL :: failure, found
    TYPE(qmmm_gaussian_type), POINTER                :: pgf
    REAL(KIND=dp) :: Rmax, Rmin, Fac
    REAL(KIND=dp) :: rc, x, t, A,G, t1, t2
    REAL(KIND=dp) :: rc2,rc3,rc4,rc5,rc6,rc7, x2,x3,x4,x5,x6,x7
    INTEGER        :: I, J, K, Np, stat, ndim, ig, ig_start, myind
    INTEGER,        DIMENSION(:), POINTER :: mm_atom_index
    REAL(KIND=dp), PARAMETER :: dx = 0.05_dp, d2x=dx*dx
    REAL(KIND=dp), DIMENSION(:,:),POINTER :: pot0_2
    REAL(KIND=dp), DIMENSION(:), POINTER :: radius
    TYPE(section_vals_type), POINTER         :: print_pot
    TYPE(cp_logger_type), POINTER :: logger
    CHARACTER(LEN=80)                        :: FileName
    CHARACTER(LEN=default_string_length)     :: myFormat

    logger => cp_error_get_logger(error)
    Rmin = 0.0_dp
    Rmax = SQRT(mm_cell%hmat(1,1)**2+&
                mm_cell%hmat(2,2)**2+&
                mm_cell%hmat(3,3)**2 ) 
    np   = CEILING(rmax/dx)+1
    print_pot => section_vals_get_subs_vals(print_section,"MM_POTENTIAL",error=error)
    !
    ! Preprocessing
    !
    IF (SIZE(mm_el_pot_radius) /= 0) THEN
       ALLOCATE(radius(1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       radius(1) = mm_el_pot_radius(1)
    ELSE
       ALLOCATE(radius(0),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
    END IF
    Loop_on_all_values: DO I = 2, SIZE(mm_el_pot_radius)
       Found=.FALSE.
       Loop_on_found_values: DO J = 1, SIZE(radius)
          IF (mm_el_pot_radius(i).EQ.radius(j)) THEN
             Found=.TRUE.
             EXIT Loop_on_found_values
          END IF
       END DO Loop_on_found_values
       IF (.NOT.Found) THEN
          Ndim = SIZE(radius)
          Ndim = Ndim + 1
          CALL REALLOCATE(radius,1,Ndim)
          radius(Ndim) = mm_el_pot_radius(i)          
       END IF
    END DO Loop_on_all_values
    !
    CPPrecondition(.NOT.ASSOCIATED(potentials),cp_failure_level,routineP,error,failure)
    ALLOCATE(potentials(SIZE(radius)),stat=stat)
    CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

    Potential_Type: DO K = 1, SIZE(radius)
       
       rc = radius(K)
       ALLOCATE(potentials(K)%Pot, stat=stat)
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
          IF (compatibility) THEN
             rc2=rc*rc; rc3=rc2*rc; rc4=rc3*rc; rc5=rc4*rc; rc6=rc5*rc; rc7=rc6*rc
             pot0_2(1,1)= 1.0_dp/rc
             pot0_2(2,1)= 0.0_dp
             pot0_2(3,1)= 0.0_dp
             x=0.0_dp
             DO i=2,np
                x=x+dx ; x2=x*x; x3=x2*x; x4=x3*x; x5=x4*x; x6=x5*x; x7=x6*x
                pot0_2(1,i)= (rc4-x4)/(rc5-x5)
                pot0_2(2,i)= x3*(-x5 + 5.0_dp*x*rc4 - 4.0_dp*rc5)/(x5-rc5)**2
                pot0_2(3,i)= 2.0_dp*x2*&
                     (-6.0_dp*rc7-8.0_dp*rc6*x-6.0_dp*rc5*x2+10.0_dp*rc3*x4+6.0_dp*rc2*x5+3.0_dp*rc*x6+x7)/&
                     (rc4+rc3*x+rc2*x2+rc*x3+x4)**3
             END DO
          ELSE
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
          END IF
          ! Print info on the unidimensional MM electrostatic potential
          IF (cp_should_output(logger%iter_info,print_pot,error=error).AND.&
               logger%para_env%mepos==logger%para_env%source) THEN
             CALL section_vals_val_get(print_pot,"FILENAME",c_val=FileName,error=error)
             CALL section_vals_val_get(print_pot,"FILENAME",c_val=FileName,error=error)
             myind = INDEX(FileName," ")
             CPPostcondition(myind /= 1,cp_failure_level,routineP,error,failure)
             WRITE(Filename(myind:),'(A1,F4.2)')"_",rc
             IF (globenv%ionode) THEN 
                OPEN(UNIT=21,FILE=FileName)
                WRITE(21,'(A)')"# MM ELECTROSTATIC POTENTIAL - UNIDIMENSIONAL - ATOMIC UNITS"
                WRITE(21,'(A,T10,A,T30,A)')"#","Xval","ElecPot"
                x=0.0_dp
                DO i=1,np
                   WRITE(21,'(T10,F15.9,T30,F15.9)')x,pot0_2(1,i)
                   x=x+dx
                END DO
                CLOSE(21)
             END IF
          END IF
          
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
          pgf => pgfs(K)%pgf
          CPPostcondition(pgf%Elp_Radius==rc,cp_failure_level,routineP,error,failure)
          ig_start = 1
          IF (compatibility) ig_start = 2
          DO Ig = ig_start, pgf%number_of_gaussians
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
          ! Print info on the unidimensional MM electrostatic potential
          IF (cp_should_output(logger%iter_info,print_pot,error=error).AND.&
               logger%para_env%mepos==logger%para_env%source) THEN
             CALL section_vals_val_get(print_pot,"FILENAME",c_val=FileName,error=error)
             myind = INDEX(FileName," ")
             CPPostcondition(myind /= 1,cp_failure_level,routineP,error,failure)
             WRITE(Filename(myind:),'(A1,F4.2)')"_",rc
             IF (globenv%ionode) THEN 
                OPEN(UNIT=21,FILE=FileName)
                WRITE(21,'(A)')"# MM ELECTROSTATIC POTENTIAL - UNIDIMENSIONAL - ATOMIC UNITS"
                WRITE(21,'(A,T10,A,T30,A,T130,A)')"#","Xval","Gaussians","LongRange"
                myFormat = "T10,F15.9,T30,"
                DO Ig=1,pgf%number_of_gaussians
                   myind = INDEX(myFormat," ")
                   WRITE(myFormat(myind:),'(A6)')"F12.6,"
                END DO
                myind = INDEX(myFormat," ")-1
                myFormat = myFormat(1:myind)//"T130,F15.9"
                myind = INDEX(myFormat," ")-1
                x=0.0_dp
                DO i=1,np
                   WRITE(21,'('//myFormat(1:myind)//')')&
                        x,(EXP(-(x/pgf%Gk(Ig))**2)*pgf%Ak(Ig),Ig=1,pgf%number_of_gaussians),pot0_2(1,i)
                   x=x+dx
                END DO
                CLOSE(21)
             END IF
          END IF

       CASE DEFAULT
          DEALLOCATE(potentials(K)%Pot)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          NULLIFY(potentials(K)%Pot)
          WRITE(*,'(A)')" QMMM Potential - Spline Interpolation - not Initialized!"
          CYCLE Potential_Type
       END SELECT
       NULLIFY(mm_atom_index)
       ALLOCATE(mm_atom_index(1),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
       ! Build mm_atom_index List
       DO J= 1, SIZE(mm_el_pot_radius)
          IF (rc.EQ.mm_el_pot_radius(J)) THEN
             Ndim = SIZE(mm_atom_index) 
             mm_atom_index(Ndim) = J
             CALL reallocate(mm_atom_index,1, Ndim+1)
          ENDIF
       END DO
       CALL reallocate(mm_atom_index, 1, Ndim)
      
       NULLIFY(potentials(K)%Pot%pot0_2)
       CALL qmmm_pot_type_create(potentials(K)%Pot, pot0_2=pot0_2,&
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

!!****f* qmmm_elpot/qmmm_validate_this_run [1.0] *
!!
!!   NAME
!!     qmmm_validate_this_run
!!
!!   FUNCTION
!!     Validates the qmmm run checking if  all setup for electrostatic
!!     potential are still fulfilled
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
!!     01.2005 created [tlaino]
!!
!!*** ********************************************************************** 
  SUBROUTINE qmmm_validate_this_run(particles_qm, particles_mm, mm_cell, qm_cell, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(particle_type), DIMENSION(:), POINTER :: particles_qm, particles_mm
    TYPE(cell_type), POINTER                   :: mm_cell, qm_cell
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error
    ! Local Variables
    REAL(KIND=dp)                 :: Rmax
    REAL(KIND=dp), DIMENSION(3,8) :: xb

    Rmax = SQRT(mm_cell%hmat(1,1)**2+&
                mm_cell%hmat(2,2)**2+&
                mm_cell%hmat(3,3)**2 ) 
    

  END SUBROUTINE qmmm_validate_this_run

END MODULE qmmm_elpot
