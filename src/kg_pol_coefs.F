!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!***** cp2k/kg_pol_coefs [1.0] *
!!
!!   NAME
!!     kg_pol_coefs
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!!
!!   SOURCE
!******************************************************************************

MODULE kg_pol_coefs
  
  USE ao_types,                 ONLY : ao_allocate, &
                                       ao_type, & 
                                       ao_zero
  USE coefficient_types,        ONLY : coeff_type, &
                                       PLANEWAVES, &
                                       ATOMICORBITALS, &
                                       PW_RECIPROCALSPACE
  USE global_types,             ONLY : global_environment_type
  USE kinds,                    ONLY : dbl
  USE md,                       ONLY : simulation_parameters_type
  USE pw_grid_types,            ONLY : pw_grid_type
  USE pw_types,                 ONLY : COMPLEXDATA1D, &
                                       pw_allocate, &
                                       pw_zero
  USE termination,              ONLY : stop_memory
  USE util,                     ONLY : gasdev

  IMPLICIT NONE

  PUBLIC :: allocate_coefs, initialize_coefs, dump_coefs

  PRIVATE


  CONTAINS

!******************************************************************************
 SUBROUTINE allocate_coefs(n, c_pos, c_vel, c_force, pw_grid )

    IMPLICIT NONE
    INTEGER, INTENT (in) :: n
    TYPE ( coeff_type ), INTENT(OUT) :: c_pos  
    TYPE ( coeff_type ), INTENT(OUT) :: c_vel  
    TYPE ( coeff_type ), INTENT(OUT) :: c_force  
    TYPE ( pw_grid_type ), INTENT(IN), OPTIONAL :: pw_grid 

! locals

   INTEGER :: inuse
   LOGICAL :: mass_flag, grid_flag, eigr_flag
     
   IF ( PRESENT ( pw_grid ) ) THEN         ! basis_type = PW
     c_pos % in_use = PLANEWAVES
     inuse = COMPLEXDATA1D
     mass_flag = .FALSE.
     CALL pw_allocate ( c_pos % pw, pw_grid, inuse, mass_flag )  
     CALL pw_zero ( c_pos % pw )
     c_pos % pw % in_space = PW_RECIPROCALSPACE
     
!   allocate cvel 
    
     c_vel % in_use = PLANEWAVES
     inuse = COMPLEXDATA1D
     mass_flag = .TRUE.
     CALL pw_allocate ( c_vel % pw, pw_grid, inuse, mass_flag )  
     CALL pw_zero ( c_vel % pw )
     c_vel % pw % in_space = PW_RECIPROCALSPACE
     
!   allocate cforce
    
     c_force % in_use = PLANEWAVES
     inuse = COMPLEXDATA1D
     mass_flag = .FALSE.
     CALL pw_allocate ( c_force % pw, pw_grid, inuse, mass_flag )  
     CALL pw_zero ( c_force % pw )
     c_force % pw % in_space = PW_RECIPROCALSPACE

   ELSE                                    ! basis_type = AO

!   allocate cpos 
    
     c_pos % in_use = ATOMICORBITALS
     mass_flag = .FALSE.
     grid_flag = .FALSE.
     CALL ao_allocate ( c_pos % ao, n, mass_flag, grid_flag )  
     CALL ao_zero ( c_pos % ao)
     
!   allocate cvel 
    
     c_vel % in_use = ATOMICORBITALS
     mass_flag = .TRUE.
     grid_flag = .FALSE.
     CALL ao_allocate ( c_vel % ao, n, mass_flag, grid_flag )  
     CALL ao_zero ( c_vel % ao)
     
!   allocate cforce
    
     c_force % in_use = ATOMICORBITALS
     mass_flag = .FALSE.
     grid_flag = .FALSE.
     CALL ao_allocate ( c_force % ao, n, mass_flag, grid_flag )  
     CALL ao_zero ( c_force % ao)
   ENDIF
     

  END SUBROUTINE allocate_coefs

!***************************************************************************
!
  SUBROUTINE dump_coefs( c_pos, c_vel, c_force )

    IMPLICIT NONE
    TYPE (coeff_type), INTENT(IN) :: c_pos
    TYPE (coeff_type), INTENT(IN) :: c_vel
    TYPE (coeff_type), INTENT(IN) :: c_force

! locals
    INTEGER :: icoef, ncoef

! write to the restart file
    OPEN ( unit = 666, file = 'coefficent.rst', status = 'unknown' )
    SELECT CASE ( c_pos % in_use )
    CASE ( ATOMICORBITALS )
      ncoef = SIZE ( c_pos % ao % cr  )
      DO icoef = 1, ncoef
         WRITE ( 666, * ) c_pos % ao % cr ( icoef ), c_vel % ao % cr ( icoef ), &
                      c_force % ao % cr ( icoef )
       END DO
      CLOSE (666)
    CASE ( PLANEWAVES )
      ncoef = SIZE ( c_pos % pw % cc  )
      DO icoef = 1, ncoef
         WRITE ( 666, * ) c_pos % pw % cc ( icoef ), c_vel % pw % cc ( icoef ), &
                      c_force % pw % cc ( icoef )
       END DO
      CLOSE (666)
    END SELECT

  END SUBROUTINE dump_coefs
!
!***************************************************************************

 SUBROUTINE initialize_coefs ( simpar, c_pos, c_vel, c_force, globenv )

    IMPLICIT NONE
    TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
    TYPE ( coeff_type ), INTENT(INOUT) :: c_pos  
    TYPE ( coeff_type ), INTENT(INOUT) :: c_vel  
    TYPE ( coeff_type ), INTENT(INOUT) :: c_force  
    TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

!   locals 
    INTEGER :: ncoef, icoef
    REAL ( dbl ) :: dum1_real, dum2_real
    COMPLEX ( dbl ) :: dum1_cplx, dum2_cplx


    SELECT CASE (simpar%read_type_coef)   

    CASE ("INIT")

      CALL initialize_coef_velocities ( simpar, c_vel, globenv )
    
      SELECT CASE ( c_pos % in_use )
      CASE ( ATOMICORBITALS )
        ncoef = size (c_pos%ao%cr)
        c_pos % ao % cr ( : ) =  1.0_dbl
      CASE ( PLANEWAVES )
        ncoef = size (c_pos%pw%cc)
        c_pos % pw % cc ( : ) =  CMPLX ( 0.0_dbl, 0.0_dbl, dbl )
        c_pos % pw % cc ( 1 ) =  CMPLX ( 0.0_dbl, 0.0_dbl, dbl )
      END SELECT 
      
    CASE ("POS")

      CALL initialize_coef_velocities ( simpar, c_vel, globenv )

      SELECT CASE ( c_pos % in_use )
      CASE ( ATOMICORBITALS )
      ncoef = size (c_pos%ao%cr)
      OPEN (unit=666,file='coefficent.rst',status='unknown')
        DO icoef = 1, ncoef
          READ (666,*)  c_pos % ao % cr ( icoef ), dum1_real, dum2_real
        END DO
      CLOSE (666)                                         
      CASE ( PLANEWAVES )
      ncoef = size (c_pos%pw%cc)
      OPEN (unit=666,file='coefficent.rst',status='unknown')
        DO icoef = 1, ncoef
          READ ( 666, * ) c_pos % pw % cc ( icoef ), dum1_cplx, dum2_cplx
        END DO
      CLOSE (666)                                         
      END SELECT 

    CASE ("POSVEL")

      OPEN (unit=666,file='coefficent.rst',status='unknown')
      SELECT CASE ( c_pos % in_use )
      CASE ( ATOMICORBITALS )
      ncoef = size (c_pos%ao%cr)
        DO icoef = 1, ncoef
          READ (666,*)  c_pos % ao % cr ( icoef ),  &
                        c_vel % ao % cr ( icoef ), &
                        dum1_real
          c_vel % ao % mass ( icoef ) = simpar % mass_coef
        END DO
      CASE ( PLANEWAVES )
      ncoef = size (c_pos%pw%cc)
        DO icoef = 1, ncoef
          READ ( 666, * ) c_pos % pw % cc ( icoef ), &
                          c_vel % pw % cc ( icoef ), & 
                          dum1_cplx
          c_vel % pw % mass_cc ( icoef ) = CMPLX ( simpar % mass_coef,  &
                                                   simpar % mass_coef, dbl )
        END DO
      END SELECT 
      CLOSE (666)                                                          

    END SELECT

 END SUBROUTINE initialize_coefs 

!------------------------------------------------------------------------------!
 SUBROUTINE initialize_coef_velocities ( simpar, c_vel, globenv )

 IMPLICIT NONE
!arguments
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( coeff_type ), INTENT ( OUT ) :: c_vel  
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
! locals
  REAL (dbl) :: akin, cmass, ctemp, sc
  INTEGER :: ncoef, i
! initializing velocities

  SELECT CASE ( c_vel % in_use )
  CASE ( ATOMICORBITALS )
    ncoef = SIZE ( c_vel % ao % cr  )
    cmass = simpar % mass_coef         

    DO i = 1, ncoef
       c_vel % ao % cr ( i ) = gasdev ( globenv % idum )
       c_vel % ao % mass ( i ) = cmass                    
    END DO
    akin = 0.0_dbl
    DO i = 1, ncoef
       akin = akin + 0.5_dbl *  &
        ( cmass * c_vel % ao % cr ( i ) * c_vel % ao % cr ( i ) )
    END DO

    ctemp = 2.0_dbl * akin / REAL ( ncoef, dbl )

! scale velocities to get the correct initial temperature

    sc = SQRT ( simpar % temp_coef / ctemp )
    DO i = 1, ncoef 
       c_vel % ao % cr ( i ) = sc *  c_vel % ao % cr ( i )
    END DO
  CASE ( PLANEWAVES )
    ncoef = SIZE ( c_vel % pw % cc  )
    cmass = simpar % mass_coef         

    DO i = 1, ncoef
        c_vel % pw % cc ( i ) = CMPLX ( gasdev ( globenv % idum ), &
                                        gasdev ( globenv % idum ), dbl )
        c_vel % pw % mass_cc ( i ) = CMPLX ( cmass, cmass, dbl )
    END DO
    akin = 0.0_dbl
    DO i = 1, ncoef
       akin = akin + 0.5_dbl *  &
        ( cmass * REAL ( c_vel % pw % cc ( i ) ) * &
                  REAL ( c_vel % pw % cc ( i ) ) + &
          cmass * AIMAG ( c_vel % pw % cc ( i ) ) * &
                  AIMAG ( c_vel % pw % cc ( i ) ) )
    END DO

    ctemp = 2.0_dbl * akin / REAL ( ncoef, dbl ) / 2._dbl

! scale velocities to get the correct initial temperature

    sc = SQRT ( simpar % temp_coef / ctemp )
    DO i = 1, ncoef 
        c_vel % pw % cc ( i ) =  CMPLX ( sc * REAL ( c_vel % pw % cc ( i ) ), &
                                         sc * AIMAG ( c_vel % pw % cc ( i ) ), dbl )
    END DO
! zeroing the g=0
    c_vel % pw % cc ( 1 ) =  CMPLX ( 0._dbl, 0._dbl, dbl )
  END SELECT 
 END SUBROUTINE initialize_coef_velocities 

!------------------------------------------------------------------------------!
   END MODULE kg_pol_coefs
!------------------------------------------------------------------------------!
