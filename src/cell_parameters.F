!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/cell_parameters [1.0] *
!!
!!   NAME
!!     cell_parameters
!!
!!   FUNCTION
!!     Definition of the periodic cell parameters.
!!
!!   AUTHOR
!!     Matthias Krack (12.05.2000)
!!
!!   MODIFICATION HISTORY
!!     1) Modified to call init_cell (Joost VandeVondele, 22.02.02)
!!     2) Inserted the multiple grid code (Joost VandeVondele 28.02.02)
!!
!!   SOURCE
!******************************************************************************

MODULE cell_parameters

! *****************************************************************************

! Index:

! SUBROUTINE read_cell_parameters(start_section,end_section)
! SUBROUTINE write_cell_parameters(lunit)

! FUNCTION pbc(r) RESULT(r_pbc)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE simulation_cell, ONLY: cell_type,init_cell,get_hinv

  IMPLICIT NONE

  PRIVATE

  TYPE(cell_type)  :: cell
  CHARACTER(LEN=8) :: unit_of_length_name = "BOHR"
  REAL(wp)         :: omega = 0.0_wp,&
                      unit_of_length = 1.0_wp
  LOGICAL          :: scaled_coordinates = .FALSE.


  REAL(wp), DIMENSION(:), pointer :: cutoff => NULL()
  REAL(wp)                        :: rel_cutoff,progression_factor
  INTEGER          :: ngrid_levels

  REAL(wp), DIMENSION(3) :: abc = (/0.0_wp,0.0_wp,0.0_wp/)

! *** Public variables ***

  PUBLIC :: abc,&
            cell,&
            cutoff,&
            omega,&
            scaled_coordinates,&
            unit_of_length,&
            unit_of_length_name,&
            ngrid_levels, &
            rel_cutoff

! *** Public subroutines ***

  PUBLIC :: pbc,&
            read_cell_parameters,&
            write_cell_parameters

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_cell_parameters(start_section,end_section,globenv)

!   Purpose: Read the periodic cell parameters.

!   History: - Creation (19.05.2000, Matthias Krack)
!   modified to call the init_cell routine (22.02.2002, Joost VandeVondele)

!   ***************************************************************************

    USE global_types,          ONLY: global_environment_type
    USE input_utilities,       ONLY: finish_input_session,&
                                     read_object,&
                                     search,&
                                     start_input_session
    USE method_specifications, ONLY: gpw
    USE physcon,               ONLY: bohr
    USE string_utilities,      ONLY: uppercase
    USE termination,           ONLY: stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_cell_parameters (MODULE cell_parameters)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword,string
    LOGICAL           :: found,user_ngrid
    real(wp)          :: hmat(3,3)
    real(wp)          :: tmp
    integer           :: periodic(3),i
    CHARACTER ( LEN = 80 ) :: poisson_solver

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    unit_of_length_name = "BOHR"
    unit_of_length = 1.0_wp
    scaled_coordinates = .FALSE.
    abc = (/0.0_wp,0.0_wp,0.0_wp/)
    periodic = (/1,1,1/)
    poisson_solver = "ANALYTIC"
    ngrid_levels = 1
    user_ngrid=.false.
    progression_factor=2.0_wp

    allocate(cutoff(1))
    cutoff(1)=160 ! 320 Ry should be a nice, expensive default
! rel_cutoff : the cutoff that would be used for a gaussian with exponent 1
! gaussians with exponent alfa will use alfa*rel_cutoff
! unless this is larger than the maximum cutoff allowed in the system
! a fully consistent calcultion would use rel_cutoff=cutoff/max_exponent
! where max_exponent is twice the largest exponent in the basis set
! rel_cutoff of about 40 Ry seems to give very good results : see thesis G.R. Lippert
! pages 93 and 134
! however, this remains a very important parameter of the calculation
    rel_cutoff=20 ! that is 40 Ry

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the input section with the cell parameters ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        SELECT CASE (TRIM(keyword))
        CASE ("ISOLATED")
          periodic=(/0,0,0/)
        CASE ("ABC")
          CALL read_object(abc(1))
          CALL read_object(abc(2))
          CALL read_object(abc(3))
        CASE ("NGRID")
          CALL read_object(ngrid_levels)
          user_ngrid=.true.
        CASE ("PROGRESSION_FACTOR")
          CALL read_object(progression_factor)
        CASE ("REL_CUTOFF", "RELATIVE_CUTOFF")
          CALL read_object(rel_cutoff)
          rel_cutoff=rel_cutoff*0.5
        CASE ("CUTOFF")
             deallocate(cutoff)
             allocate(cutoff(ngrid_levels))
             DO I=1,ngrid_levels
                CALL read_object(cutoff(i))
                cutoff(i) = 0.5_wp*cutoff(i)
             ENDDO
        CASE ("UNIT")
          CALL read_object(string)
          SELECT CASE (TRIM(string))
          CASE ("ANGSTROM")
            unit_of_length = bohr
            unit_of_length_name = "ANGSTROM"
            scaled_coordinates = .FALSE.
          CASE ("BOHR")
            unit_of_length = 1.0_wp
            unit_of_length_name = "BOHR"
            scaled_coordinates = .FALSE.
          CASE ("SCALED_ANGSTROM")
            unit_of_length = bohr
            unit_of_length_name = "ANGSTROM"
            scaled_coordinates = .TRUE.
          CASE ("SCALED_BOHR")
            unit_of_length = 1.0_wp
            unit_of_length_name = "BOHR"
            scaled_coordinates = .TRUE.
          CASE DEFAULT
            CALL stop_program(routine,&
                              "Invalid name for the unit of length <"//&
                              TRIM(string)//&
                              "> found in the input section <"//&
                              TRIM(start_section)//">")
          END SELECT
        CASE DEFAULT
          IF (keyword == end_section) THEN
            EXIT
          ELSE
            CALL stop_program(routine,&
                              "Invalid keyword <"//TRIM(keyword)//&
                              "> found in the input section <"//&
                              TRIM(start_section)//">")
          END IF
        END SELECT

      END DO

      abc(:) = abc(:)*unit_of_length
      omega = PRODUCT(abc(:))

! atomatically generate the gridlevels using a geometric progression for the cutoff
      if (.not. user_ngrid) then
         ngrid_levels=MAX(1,1+CEILING(log(cutoff(1)/(0.05*rel_cutoff))/log(progression_factor)))
         if (ngrid_levels.gt.1) then
            tmp=cutoff(1)
            deallocate(cutoff)
            allocate(cutoff(ngrid_levels)) 
            cutoff(1)=tmp
            do i=2,ngrid_levels
               cutoff(i)=cutoff(i-1)/progression_factor 
            enddo
         endif
      endif

    ELSE

!     *** Stop program, if the end of file is reached ***

      IF (gpw) THEN
        CALL stop_program(routine,&
                          "No input section <"//TRIM(start_section)//&
                          "> found in the input file <"//&
                          TRIM(globenv%input_file_name)//">")
      END IF

    END IF

    CALL finish_input_session

!   *** Check the cell parameters ***

    IF (gpw.OR.scaled_coordinates) THEN
      IF ((abc(1) <= 0.0_wp).OR.&
          (abc(2) <= 0.0_wp).OR.&
          (abc(3) <= 0.0_wp)) THEN
        CALL stop_program(routine,"Invalid cell dimensions were specified")
      END IF
    END IF

    hmat(:,:) = 0.0_wp
    hmat(1,1) = abc(1)
    hmat(2,2) = abc(2)
    hmat(3,3) = abc(3)

    CALL init_cell(cell,hmat,periodic,poisson_solver)

  END SUBROUTINE read_cell_parameters

! *****************************************************************************

  SUBROUTINE write_cell_parameters(lunit)

!   Purpose: Write the periodic cell parameters to the logical unit number
!            "lunit".

!   History: - Creation (02.06.2000, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: lunit
    INTEGER :: i

!   ---------------------------------------------------------------------------

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/)") "CELL PARAMETERS"

    WRITE (UNIT=lunit,FMT="(T3,A,T26,3F10.3,/,/,T3,A,T31,F15.3,A)")&
      "Cell size in "//TRIM(unit_of_length_name)//":",&
      abc(:)/unit_of_length,&
      "Cell volume in "//TRIM(unit_of_length_name)//"**3:",&
      omega/unit_of_length**3
    WRITE (UNIT=lunit,FMT='(T3,A)') " "
    WRITE (UNIT=lunit,FMT='(T3,A,I6)')         "Number of grids     : ",ngrid_levels
    WRITE (UNIT=lunit,FMT='(T3,A,F15.3)')      "Relative Cutoff a.u.: ",rel_cutoff
    WRITE (UNIT=lunit,FMT='(T3,A,F15.3)')      "Maximum  Cutoff a.u.: ",cutoff(1)
    WRITE (UNIT=lunit,FMT='(T3,A,F15.3)')      "Progression Factor  : ",progression_factor
    do i=1,ngrid_levels
    WRITE (UNIT=lunit,FMT='(T3,A,I3,A,F15.3)') "Cutoff for grid ",i," : ",cutoff(i)
    enddo

  END SUBROUTINE write_cell_parameters

! *****************************************************************************

  FUNCTION pbc(r) RESULT(r_pbc)

!   Purpose: Apply the periodic boundary condition to the position vector r.

!   History: - Creation (22.05.2000, Matthias Krack)

!   ***************************************************************************

    USE method_specifications, ONLY: periodic

    REAL(wp), DIMENSION(3) :: r

    REAL(wp), DIMENSION(3) :: r_pbc

!   ---------------------------------------------------------------------------

    IF (periodic) THEN
      r_pbc(:) = r(:) - abc(:)*ANINT(r(:)/abc(:))
    ELSE
      r_pbc(:) = r(:)
    END IF

  END FUNCTION pbc

! *****************************************************************************

END MODULE cell_parameters
