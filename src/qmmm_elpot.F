!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_elpot
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE erf_fn,                          ONLY: erf
  USE input_constants,                 ONLY: do_qmmm_coulomb,&
                                             do_qmmm_gauss,&
                                             do_qmmm_pcharge,&
                                             do_qmmm_swave
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: rootpi
  USE memory_utilities,                ONLY: reallocate
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qmmm_types,                      ONLY: qmmm_Pot_Type,&
                                             qmmm_pot_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_elpot'
  PUBLIC :: qmmm_potential_init

CONTAINS

! *****************************************************************************
!> \brief Initialize the QMMM potential stored on vector,
!>      according the qmmm_coupl_type
!> \param qmmm_coupl_type ...
!> \param mm_el_pot_radius ...
!> \param potentials ...
!> \param pgfs ...
!> \param mm_cell ...
!> \param compatibility ...
!> \param print_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_potential_init(qmmm_coupl_type, mm_el_pot_radius, potentials,&
    pgfs, mm_cell, compatibility, print_section,&
    error)
    INTEGER, INTENT(IN)                      :: qmmm_coupl_type
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius
    TYPE(qmmm_pot_p_type), DIMENSION(:), &
      POINTER                                :: potentials
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: pgfs
    TYPE(cell_type), POINTER                 :: mm_cell
    LOGICAL, INTENT(IN)                      :: compatibility
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_potential_init', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: dx = 0.05_dp, d2x = dx*dx

    CHARACTER(LEN=default_path_length)       :: myFormat
    CHARACTER(LEN=default_string_length)     :: rc_s
    INTEGER                                  :: I, ig, ig_start, J, K, myind, &
                                                ndim, Np, output_unit, stat, &
                                                unit_nr
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index
    LOGICAL                                  :: failure, found
    REAL(KIND=dp)                            :: A, G, rc, Rmax, Rmin, t, t1, &
                                                t2, x
    REAL(KIND=dp), DIMENSION(:), POINTER     :: radius
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pot0_2
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qmmm_gaussian_type), POINTER        :: pgf

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    Rmin = 0.0_dp
    Rmax = SQRT(mm_cell%hmat(1,1)**2+&
                mm_cell%hmat(2,2)**2+&
                mm_cell%hmat(3,3)**2 )
    np   = CEILING(rmax/dx)+1
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
       CASE(do_qmmm_coulomb)
          NULLIFY(pot0_2)
       CASE(do_qmmm_pcharge)
          NULLIFY(pot0_2)
       CASE(do_qmmm_gauss,do_qmmm_swave)
          ALLOCATE(pot0_2(2,np), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END SELECT

       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_coulomb,do_qmmm_pcharge)
          ! Do Nothing
       CASE(do_qmmm_gauss,do_qmmm_swave)
          IF     (qmmm_coupl_type==do_qmmm_gauss) THEN
             ! Smooth Coulomb Potential ::  Erf(x/rc)/x
             pot0_2(1,1)= 2.0_dp/(rootpi*rc)
             pot0_2(2,1)= 0.0_dp
             x=0.0_dp
             DO i=2,np
                x=x+dx
                pot0_2(1,i)=erf(x/rc)/x
                t=2._dp/(rootpi*x*rc)*EXP(-(x/rc)**2)
                pot0_2(2,i)=(t-pot0_2(1,i)/x)*dx
             END DO
          ELSEIF (qmmm_coupl_type==do_qmmm_swave) THEN
             ! S-wave expansion :: 1/x - exp(-2*x/rc) * ( 1/x - 1/rc )
             pot0_2(1,1)= 1.0_dp/rc
             pot0_2(2,1)= 0.0_dp
             x=0.0_dp
             DO i=2,np
                x=x+dx
                t=EXP(-2.0_dp*x/rc)/rc
                pot0_2(1,i)=(1.0_dp-t*(rc+x))/x
                pot0_2(2,i)=((t*(rc**2+2.0_dp*rc*x+2.0_dp*x**2)/rc-1.0_dp)/x**2)*dx
             END DO
          END IF
          pgf => pgfs(K)%pgf
          CPPostcondition(pgf%Elp_Radius==rc,cp_failure_level,routineP,error,failure)
          ig_start = 1
          IF (compatibility.AND.(qmmm_coupl_type==do_qmmm_gauss)) ig_start = 2
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
          IF (BTEST(cp_print_key_should_output(logger%iter_info,print_section,"MM_POTENTIAL",&
                    error=error),cp_p_file)) THEN
             WRITE(rc_s,'(F6.3)')rc
             unit_nr = cp_print_key_unit_nr(logger,print_section,"MM_POTENTIAL",&
                  extension="_rc="//TRIM(ADJUSTL(rc_s))//".data",error=error)
             IF (unit_nr>0) THEN
                WRITE(unit_nr,'(A)')"# MM ELECTROSTATIC POTENTIAL - UNIDIMENSIONAL - ATOMIC UNITS"
                WRITE(unit_nr,'(A,I5)')"# MM ELECTROSTATIC POTENTIAL - Nr. of Gaussians:",pgf%number_of_gaussians
                WRITE(unit_nr,'(A,T10,A,T30,A,T300,A)')"#","Xval","Gaussians","LongRange"
                myFormat = "T10,F15.9,T30,"
                DO Ig=1,pgf%number_of_gaussians
                   myind = INDEX(myFormat," ")
                   WRITE(myFormat(myind:),'(A6)')"F12.9,"
                END DO
                myind = INDEX(myFormat," ")-1
                myFormat = myFormat(1:myind)//"T300,F15.9"
                myind = INDEX(myFormat," ")-1
                x=0.0_dp
                DO i=1,np
                   WRITE(unit_nr,'('//myFormat(1:myind)//')')&
                        x,(EXP(-(x/pgf%Gk(Ig))**2)*pgf%Ak(Ig),Ig=1,pgf%number_of_gaussians),pot0_2(1,i)
                   x=x+dx
                END DO
             END IF
             CALL cp_print_key_finished_output(unit_nr,logger,print_section,&
                  "MM_POTENTIAL", error=error)
          END IF
       CASE DEFAULT
          DEALLOCATE(potentials(K)%Pot)
          CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
          NULLIFY(potentials(K)%Pot)
          IF (output_unit>0) WRITE(output_unit,'(A)')" QMMM Potential - Spline Interpolation - not Initialized!"
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

! *****************************************************************************
!> \brief Creates the qmmm_pot_type structure
!> \param Pot ...
!> \param pot0_2 ...
!> \param Rmax ...
!> \param Rmin ...
!> \param dx ...
!> \param npts ...
!> \param rc ...
!> \param mm_atom_index ...
!> \par History
!>      09.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_pot_type_create(Pot, pot0_2, Rmax, Rmin, dx, npts, rc,&
       mm_atom_index)
    TYPE(qmmm_Pot_Type), POINTER             :: Pot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pot0_2
    REAL(KIND=dp), INTENT(IN)                :: Rmax, Rmin, dx
    INTEGER, INTENT(IN)                      :: npts
    REAL(KIND=dp), INTENT(IN)                :: Rc
    INTEGER, DIMENSION(:), POINTER           :: mm_atom_index

    Pot%pot0_2        => pot0_2
    Pot%mm_atom_index => mm_atom_index
    Pot%Rmax          =  Rmax
    Pot%Rmin          =  Rmin
    Pot%Rc            =  rc
    Pot%dx            =  dx
    Pot%npts          =  npts

  END SUBROUTINE qmmm_pot_type_create

END MODULE qmmm_elpot
