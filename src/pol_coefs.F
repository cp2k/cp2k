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
  
  USE atomic_kinds, ONLY :  kind_info_type
  USE ao_types, ONLY : ao_type, ao_allocate, ao_zero
  USE coefficient_types, ONLY : coeff_type
  USE kinds, ONLY : dbl
  USE particle_types, ONLY : particle_type
  USE molecule_types , ONLY : molecule_type, read_basis_type 
  USE termination, ONLY : stop_memory

  IMPLICIT NONE

  PUBLIC :: allocate_pol_coefs, initialize_pol_coefs

  PRIVATE


  CONTAINS

!******************************************************************************
 SUBROUTINE allocate_pol_coefs(ncoef, c_pos, c_vel, c_force)

    IMPLICIT NONE
    INTEGER, INTENT (in) :: ncoef
    TYPE ( coeff_type ), INTENT(OUT) :: c_pos  
    TYPE ( coeff_type ), INTENT(OUT) :: c_vel  
    TYPE ( coeff_type ), INTENT(OUT) :: c_force  

! locals

   INTEGER, PARAMETER :: ATOMICORBITALS = 402
   LOGICAL :: mass_flag, grid_flag
     
!   allocate cpos 
    
   c_pos % in_use = ATOMICORBITALS
   mass_flag = .TRUE.
   grid_flag = .TRUE.
   CALL ao_allocate ( c_pos % ao, ncoef, mass_flag, grid_flag )  
   CALL ao_zero ( c_pos % ao)
     
!   allocate cvel 
    
   c_vel % in_use = ATOMICORBITALS
   mass_flag = .FALSE.
   grid_flag = .FALSE.
   CALL ao_allocate ( c_vel % ao, ncoef, mass_flag, grid_flag )  
   CALL ao_zero ( c_vel % ao)
     
!   allocate cforce
    
   c_force % in_use = ATOMICORBITALS
   mass_flag = .FALSE.
   grid_flag = .FALSE.
   CALL ao_allocate ( c_force % ao, ncoef, mass_flag, grid_flag )  
   CALL ao_zero ( c_force % ao)

  END SUBROUTINE allocate_pol_coefs

!***************************************************************************

 SUBROUTINE initialize_pol_coefs ( rtype , ki, part, c_pos, c_vel, c_force)

    IMPLICIT NONE
    CHARACTER (LEN = 12),  INTENT (in) :: rtype
    TYPE (kind_info_type), POINTER, DIMENSION(:) :: ki
    TYPE (particle_type), intent(in), DIMENSION(:) :: part
    TYPE ( coeff_type ), INTENT(OUT) :: c_pos  
    TYPE ( coeff_type ), INTENT(OUT) :: c_vel  
    TYPE ( coeff_type ), INTENT(OUT) :: c_force  

!   locals 
    INTEGER :: ncoef, icoef
    LOGICAL :: kind_flag

    ncoef = size (c_pos%ao%cr(:))

    SELECT CASE (rtype)   

    CASE ("INIT")

!!    CALL initialize_coeff_velocities  (c_vel, stuff)!!    
    
      c_pos % ao % cr ( : ) =  1.0_dbl
      
    CASE ("POS")

!!    CALL initialize_coeff_velocities  (c_vel, stuff)!!    

      OPEN (unit=666,file='coefficent.rst',status='unknown')
      DO icoef = 1, ncoef
        READ (666,*) c_pos % ao % cr (icoef)
      END DO
      CLOSE (666)                                         

    CASE ("POSVEL")

      OPEN (unit=666,file='coefficent.rst',status='unknown')
      DO icoef = 1, ncoef
        READ (666,*)  c_pos % ao % cr (icoef),  c_vel % ao % cr (icoef)
      END DO
      CLOSE (666)                                                          

    END SELECT	

! fills the remaining information...

   kind_flag = .TRUE.
   CALL get_kind_and_part_index ( c_pos % ao, ki, part, kind_flag)
   
   kind_flag = .FALSE.
   CALL get_kind_and_part_index ( c_vel % ao, ki, part, kind_flag)
   
   kind_flag = .FALSE.
   CALL get_kind_and_part_index ( c_force % ao, ki, part, kind_flag)

 END SUBROUTINE initialize_pol_coefs 

!------------------------------------------------------------------------------!
   
 SUBROUTINE get_kind_and_part_index (ao, ki, part, flag)
  
   IMPLICIT NONE
   TYPE (ao_type), INTENT (INOUT)  :: ao
   TYPE (kind_info_type), POINTER, DIMENSION(:) :: ki
   TYPE (particle_type), intent(in), DIMENSION(:) :: part
   LOGICAL, intent (in) :: flag

! locals
  INTEGER ::  ncoef, ikind, ipart, ii, icoef, icgf, ncgf, ios
  TYPE (kind_info_type), POINTER  :: kinfo

   ncoef = size (ao % cr)
   DO ikind = 1 , size (ki)
   
      DO ii = 1, ki(ikind)%natom
 
        ipart = ki(ikind)%atom_list(ii)        
	ncgf = ki(ikind) % orb_basis_set % ncgf 

          DO icgf = 1, ncgf
             icoef  = part ( ipart ) % coef_list( icgf )
	     ao % coef_to_basis_set ( icoef ) = ikind 
	     ao % coef_to_part ( icoef ) = ipart 
	     IF (flag) THEN
	       ALLOCATE (ao % kind_info (ncoef), STAT = ios)
	       kinfo =>   ao % kind_info ( icoef )
	       kinfo => ki (ikind)
	     END IF
         ENDDO                                
	 
      END DO	 
      
   END DO	 

 END SUBROUTINE get_kind_and_part_index
!------------------------------------------------------------------------------!

   END MODULE pol_coefs
!------------------------------------------------------------------------------!

