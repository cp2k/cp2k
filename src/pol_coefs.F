!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!***** cp2k/pol_coefs [1.0] *
!!
!!   NAME
!!     pol_coefs
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

MODULE pol_coefs
  
  USE ao_types,                        ONLY: ao_allocate,&
                                             ao_type,&
                                             ao_zero
  USE atomic_kinds,                    ONLY: kind_info_type
  USE coefficient_types,               ONLY: ATOMICORBITALS,&
                                             PLANEWAVES,&
                                             PW_RECIPROCALSPACE,&
                                             coeff_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE md,                              ONLY: simulation_parameters_type
  USE particle_types,                  ONLY: particle_type
  USE pol_projectors_scale,            ONLY: get_proj_scale_factor
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             pw_allocate,&
                                             pw_zero
  USE util,                            ONLY: gasdev

  IMPLICIT NONE

  PUBLIC :: allocate_pol_coefs, initialize_pol_coefs, dump_pol_coefs

  PRIVATE


  CONTAINS

!******************************************************************************
 SUBROUTINE allocate_pol_coefs(n, c_pos, c_vel, c_force, projector, pw_grid )

    INTEGER, INTENT(in)                      :: n
    TYPE(coeff_type), INTENT(OUT)            :: c_pos, c_vel, c_force
    TYPE(coeff_type), INTENT(OUT), OPTIONAL  :: projector
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: pw_grid

    INTEGER                                  :: inuse
    LOGICAL                                  :: eigr_flag, grid_flag, &
                                                mass_flag

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

!   allocate projector

     projector % in_use = ATOMICORBITALS
     mass_flag = .FALSE.
     eigr_flag = .TRUE.
     grid_flag = .TRUE.
     CALL ao_allocate ( projector % ao, n, mass_flag, grid_flag, eigr_flag )  
     CALL ao_zero ( projector % ao )

   ELSE                                    ! basis_type = AO

!   allocate cpos 
    
     c_pos % in_use = ATOMICORBITALS
     mass_flag = .FALSE.
     grid_flag = .TRUE.
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
     

  END SUBROUTINE allocate_pol_coefs

!***************************************************************************
!
  SUBROUTINE dump_pol_coefs( c_pos, c_vel, c_force )

    TYPE(coeff_type), INTENT(IN)             :: c_pos, c_vel, c_force

    INTEGER                                  :: icoef, ncoef

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

  END SUBROUTINE dump_pol_coefs
!
!***************************************************************************

 SUBROUTINE initialize_pol_coefs ( simpar , ki, part, c_pos, c_vel, c_force, &
                                   globenv, proj )

    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: ki
    TYPE(particle_type), DIMENSION(:), &
      intent(in)                             :: part
    TYPE(coeff_type), INTENT(OUT)            :: c_pos, c_vel, c_force
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv
    TYPE(coeff_type), INTENT(OUT), OPTIONAL  :: proj

    COMPLEX(dbl)                             :: dum1_cplx, dum2_cplx
    INTEGER                                  :: icoef, ncoef
    LOGICAL                                  :: kind_flag
    REAL(dbl)                                :: dum1_real, dum2_real

    SELECT CASE (simpar%pol_read_type)   

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
          c_vel % ao % mass ( icoef ) = simpar % pol_mass
        END DO
      CASE ( PLANEWAVES )
      ncoef = size (c_pos%pw%cc)
        DO icoef = 1, ncoef
          READ ( 666, * ) c_pos % pw % cc ( icoef ), &
                          c_vel % pw % cc ( icoef ), & 
                          dum1_cplx
          c_vel % pw % mass_cc ( icoef ) = CMPLX ( simpar % pol_mass,  &
                                                   simpar % pol_mass, dbl )
        END DO
      END SELECT 
      CLOSE (666)                                                          

    END SELECT

! fills the remaining information...

   IF ( PRESENT ( proj ) ) THEN
     kind_flag = .TRUE.
     CALL get_kind_and_part_index ( proj % ao, ki, part, kind_flag)
     CALL get_proj_scale_factor ( proj % ao, ki, part)
   ELSE
     kind_flag = .TRUE.
     CALL get_kind_and_part_index ( c_pos % ao, ki, part, kind_flag)
   
     kind_flag = .FALSE.
     CALL get_kind_and_part_index ( c_vel % ao, ki, part, kind_flag)
   
     kind_flag = .FALSE.
     CALL get_kind_and_part_index ( c_force % ao, ki, part, kind_flag)
   ENDIF

 END SUBROUTINE initialize_pol_coefs 

!------------------------------------------------------------------------------!
 SUBROUTINE initialize_coef_velocities ( simpar, c_vel, globenv )

    TYPE(simulation_parameters_type), &
      INTENT(IN)                             :: simpar
    TYPE(coeff_type), INTENT(OUT)            :: c_vel
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    INTEGER                                  :: i, ncoef
    REAL(dbl)                                :: akin, cmass, ctemp, sc

! initializing velocities

  SELECT CASE ( c_vel % in_use )
  CASE ( ATOMICORBITALS )
    ncoef = SIZE ( c_vel % ao % cr  )
    cmass = simpar % pol_mass         

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

    sc = SQRT ( simpar % pol_temp / ctemp )
    DO i = 1, ncoef 
       c_vel % ao % cr ( i ) = sc *  c_vel % ao % cr ( i )
    END DO
  CASE ( PLANEWAVES )
    ncoef = SIZE ( c_vel % pw % cc  )
    cmass = simpar % pol_mass         

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

    sc = SQRT ( simpar % pol_temp / ctemp )
    DO i = 1, ncoef 
        c_vel % pw % cc ( i ) =  CMPLX ( sc * REAL ( c_vel % pw % cc ( i ) ), &
                                         sc * AIMAG ( c_vel % pw % cc ( i ) ), dbl )
    END DO
! zeroing the g=0
    c_vel % pw % cc ( 1 ) =  CMPLX ( 0._dbl, 0._dbl, dbl )
  END SELECT 
 END SUBROUTINE initialize_coef_velocities 

!------------------------------------------------------------------------------!
   
 SUBROUTINE get_kind_and_part_index (ao, ki, part, flag)
  
    TYPE(ao_type), intent(inout)             :: ao
    TYPE(kind_info_type), DIMENSION(:), &
      POINTER                                :: ki
    TYPE(particle_type), DIMENSION(:), &
      intent(in)                             :: part
    LOGICAL, intent(in)                      :: flag

    INTEGER                                  :: first_cgf, iat, icgf, icoef, &
                                                ii, ikind, ios, ipart, iset, &
                                                ishell, last_cgf, natoms, &
                                                ncoef, nset

   ncoef = size (ao % cr)

   IF (flag) THEN
      NULLIFY (ao % kind_info )
      IF (.NOT.ASSOCIATED(ao % kind_info)) THEN         
         ALLOCATE (ao % kind_info (ncoef), STAT = ios)
      END IF
   END IF

   DO ikind = 1 , size (ki)
   
      natoms = ki (ikind ) % natom

      DO ii = 1, natoms
 
        ipart = ki(ikind) % atom_list(ii)        
        nset = ki(ikind) % orb_basis_set % nset 

        DO iset = 1, nset

           DO ishell=1,ki(ikind) % orb_basis_set % nshell(iset)

              first_cgf = ki(ikind) % orb_basis_set %first_cgf(ishell,iset)
              last_cgf = ki(ikind) % orb_basis_set %last_cgf(ishell,iset)
              
              DO icgf = first_cgf, last_cgf
                 icoef  = part ( ipart ) % coef_list( icgf )
                 ao % coef_to_basis ( icoef ) = ikind 
                 ao % coef_to_set ( icoef ) = iset 
                 ao % coef_to_part ( icoef ) = ipart 

                 IF (flag) THEN
                   ao % kind_info (icoef) % orb_basis_set => ki (ikind) % orb_basis_set
                   ao % kind_info (icoef) % orb_basis_set_name = ki (ikind) % orb_basis_set_name
                   ao % kind_info (icoef) % number_of_grid_points = &
                        ki (ikind) % number_of_grid_points
                   ao % kind_info (icoef) % element_symbol = ki (ikind) % element_symbol
                   ao % kind_info (icoef) % natom = ki (ikind) % natom
                   ao % kind_info (icoef) % z = ki (ikind) % z
                   NULLIFY (ao % kind_info ( icoef ) % atom_list )
                   IF (.NOT.ASSOCIATED(ao % kind_info (icoef) % atom_list)) THEN         
                     ALLOCATE (ao % kind_info (icoef) % atom_list (natoms) , STAT = ios)
                   END IF
                   DO iat = 1, natoms
                      ao % kind_info (icoef) % atom_list (iat) = ki (ikind) % atom_list (iat)
                   END DO
                 END IF
               END DO                                
            END DO
         END DO
      END DO
   END DO

 END SUBROUTINE get_kind_and_part_index

!------------------------------------------------------------------------------!

   END MODULE pol_coefs
!------------------------------------------------------------------------------!

