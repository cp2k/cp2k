!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_grid_atom [1.0] *
!!
!!   NAME
!!     qs_grid_atom
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE qs_grid_atom

  USE kinds,                           ONLY: default_string_length,&
                                             dp,&
                                             int_size
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate

  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_grid_atom"

  TYPE grid_atom_type

     INTEGER                         :: nr, ng_sphere
     REAL(dp), DIMENSION(:), POINTER :: rad, rad2,&
                                        wr, wa,&
                                        cos_theta, sin_theta,&
                                        usin_theta, cos_phi, sin_phi,&
                                        theta, cotan_theta, phi
     REAL(dp), DIMENSION(:,:), &
                             POINTER :: rad2l,oorad2l

  END TYPE grid_atom_type


! *** Public subroutines ***

  PUBLIC :: allocate_grid_atom, &
            create_grid_atom, &
            deallocate_grid_atom


! *** Public data types ***

  PUBLIC :: grid_atom_type


!!***
! *****************************************************************************

 CONTAINS

! *****************************************************************************


  SUBROUTINE allocate_grid_atom(grid_atom)


!   ***************************************************************************

    TYPE(grid_atom_type), POINTER     :: grid_atom

    CHARACTER(LEN=*), PARAMETER :: routine_name = "allocate_grid_atom"

    INTEGER                           :: iat, istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(grid_atom)) THEN
      CALL deallocate_grid_atom(grid_atom)
    END IF

    ALLOCATE (grid_atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                            "allocate_grid_atom",0)

    NULLIFY (grid_atom%rad)
    NULLIFY (grid_atom%rad2)
    NULLIFY (grid_atom%wr)
    NULLIFY (grid_atom%wa)
    NULLIFY (grid_atom%cos_theta)
    NULLIFY (grid_atom%sin_theta)
    NULLIFY (grid_atom%usin_theta)
    NULLIFY (grid_atom%cos_phi)
    NULLIFY (grid_atom%sin_phi)
    NULLIFY (grid_atom%cotan_theta)
    NULLIFY (grid_atom%theta)
    NULLIFY (grid_atom%phi)
    NULLIFY (grid_atom%rad2l)
    NULLIFY (grid_atom%oorad2l)

  END SUBROUTINE allocate_grid_atom

! *****************************************************************************

  SUBROUTINE deallocate_grid_atom(grid_atom)

!   Purpose: Deallocate a Gaussian-type orbital (GTO) basis set data set.

!   History: - Creation (03.11.2000,MK)

!   ***************************************************************************

    TYPE(grid_atom_type), POINTER        :: grid_atom

    CHARACTER(LEN=*), PARAMETER :: routine_name = "deallocate_grid_atom"

    INTEGER                                  :: istat

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(grid_atom)) THEN

      IF (ASSOCIATED(grid_atom%rad)) THEN
         DEALLOCATE (grid_atom%rad,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%rad")
      END IF

      IF (ASSOCIATED(grid_atom%rad2)) THEN
         DEALLOCATE (grid_atom%rad2,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%rad2")
      END IF

      IF (ASSOCIATED(grid_atom%wr)) THEN
         DEALLOCATE (grid_atom%wr,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%wr")
      END IF

      IF (ASSOCIATED(grid_atom%wa)) THEN
         DEALLOCATE (grid_atom%wa,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%wa")
      END IF

      IF (ASSOCIATED(grid_atom%cos_theta)) THEN
         DEALLOCATE (grid_atom%cos_theta,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%cos_theta")
      END IF

      IF (ASSOCIATED(grid_atom%sin_theta)) THEN
         DEALLOCATE (grid_atom%sin_theta,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%sin_theta")
      END IF

      IF (ASSOCIATED(grid_atom%usin_theta)) THEN
         DEALLOCATE (grid_atom%usin_theta,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%usin_theta")
      END IF

      IF (ASSOCIATED(grid_atom%cos_phi)) THEN
         DEALLOCATE (grid_atom%cos_phi,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%cos_phi")
      END IF

      IF (ASSOCIATED(grid_atom%sin_phi)) THEN
         DEALLOCATE (grid_atom%sin_phi,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%sin_phi")
      END IF

      IF (ASSOCIATED(grid_atom%cotan_theta)) THEN
         DEALLOCATE (grid_atom%cotan_theta,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%cotan_theta")
      END IF

      IF (ASSOCIATED(grid_atom%theta)) THEN
         DEALLOCATE (grid_atom%theta,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%theta")
      END IF

      IF (ASSOCIATED(grid_atom%phi)) THEN
         DEALLOCATE (grid_atom%phi,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%phi")
      END IF

      IF (ASSOCIATED(grid_atom%rad2l)) THEN
         DEALLOCATE (grid_atom%rad2l,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%rad2l")
      END IF

      IF (ASSOCIATED(grid_atom%oorad2l)) THEN
         DEALLOCATE (grid_atom%oorad2l,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom%oorad2l")
      END IF

      DEALLOCATE (grid_atom,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "grid_atom")
    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "The pointer grid_atom is not associated and "//&
                        "cannot be deallocated")
    END IF
  END SUBROUTINE deallocate_grid_atom

! *****************************************************************************

  SUBROUTINE create_grid_atom(grid_atom,nr,na,llmax,quadrature)

    TYPE(grid_atom_type), POINTER   :: grid_atom
    INTEGER, INTENT(IN)             :: llmax,na,nr, quadrature

    CHARACTER(LEN=*), PARAMETER         :: routine_name ="create_grid_atom"

    REAL(dp), DIMENSION(:), POINTER :: rad,rad2,wr 
    INTEGER                         :: l

! *****************************************************************************

    NULLIFY (rad,rad2,wr)
 
    IF(ASSOCIATED(grid_atom)) THEN
  
!     *** allocate the radial grid arrays *** 
      CALL reallocate(grid_atom%rad,1,nr)
      CALL reallocate(grid_atom%rad2,1,nr)
      CALL reallocate(grid_atom%wr,1,nr)
      CALL reallocate(grid_atom%wa,1,na)
      CALL reallocate(grid_atom%cos_theta,1,na)
      CALL reallocate(grid_atom%sin_theta,1,na)
      CALL reallocate(grid_atom%usin_theta,1,na)
      CALL reallocate(grid_atom%cos_phi,1,na)
      CALL reallocate(grid_atom%sin_phi,1,na)
      CALL reallocate(grid_atom%cotan_theta,1,na)
      CALL reallocate(grid_atom%theta,1,na)
      CALL reallocate(grid_atom%phi,1,na)
      CALL reallocate(grid_atom%rad2l,1,nr,0,llmax+1)
      CALL reallocate(grid_atom%oorad2l,1,nr,0,llmax+1)

!     *** calculate the radial grid for this kind ***
      rad => grid_atom%rad
      rad2 => grid_atom%rad2
      wr => grid_atom%wr

      CALL radial_grid(nr,rad,rad2,wr,quadrature) 

      grid_atom%rad2l(:,0) = 1._dp
      grid_atom%oorad2l(:,0) = 1._dp
      DO l = 1,llmax+1
        grid_atom%rad2l(:,l)=grid_atom%rad2l(:,l-1)*rad(:)
        grid_atom%oorad2l(:,l)=grid_atom%oorad2l(:,l-1)/rad(:)
      ENDDO

    ELSE
      CALL stop_program(routine_name,module_name,__LINE__,&
               "The pointer grid_atom is not associated")
    END IF
  END SUBROUTINE create_grid_atom

! *****************************************************************************

  SUBROUTINE radial_grid(n,r,r2,wr,radial_quadrature)

!   Purpose: Generate a radial grid with n points by a quadrature rule.

!   Literature: - A. D. Becke, J. Chem. Phys. 88, 2547 (1988)
!               - J. M. Perez-Jorda, A. D. Becke and E. San-Fabian,
!                 J. Chem. Phys. 100, 6520 (1994)
!               - M. Krack and A. M. Koester, J. Chem. Phys. 108, 3226 (1998)

!   History: - Creation (20.09.1999, Matthias Krack)

!   ***************************************************************************

    REAL(dp), DIMENSION(:), POINTER :: r,r2,wr
    INTEGER, INTENT(IN)             :: n,radial_quadrature

!   *** Local variables ***

    REAL(dp) :: cost,f,t,sint,sint2,w,x
    INTEGER  :: i

!   ---------------------------------------------------------------------------

    f = pi/REAL(n + 1,dp)

    SELECT CASE (radial_quadrature)

    CASE (1)

!     *** Gauss-Chebyshev quadrature formula of the second kind ***
!     *** u [-1,+1] -> r [0,infinity] =>  r = (1 + u)/(1 - u)   ***

      DO i=1,n
        t = REAL(i,dp)*f
        x = COS(t)
        w = f*SIN(t)**2
        r(i) = (1.0_dp + x)/(1.0_dp - x)
        r2(i) = r(i)**2
        wr(i) = w/SQRT(1.0_dp - x**2)
        wr(i) = 2.0_dp*wr(i)*r2(i)/(1.0_dp - x)**2
      END DO

    CASE (2)

!     *** Transformed Gauss-Chebyshev quadrature formula of the second kind ***
!     *** u [-1,+1] -> r [0,infinity] => r = (1 + u)/(1 - u)                ***

      DO i=1,n
        t = REAL(i,dp)*f
        cost = COS(t)
        sint = SIN(t)
        sint2 = sint**2
        x = REAL(2*i - n - 1,dp)/REAL(n + 1,dp) -&
            2.0_dp*(1.0_dp + 2.0_dp*sint2/3.0_dp)*cost*sint/pi
        w = 16.0_dp*sint2**2/REAL(3*(n + 1),dp)
        r(n+1-i) = (1.0_dp + x)/(1.0_dp - x)
        r2(n+1-i) = r(n+1-i)**2
        wr(n+1-i) = 2.0_dp*w*r2(n+1-i)/(1.0_dp - x)**2
      END DO

    CASE (3)

!     *** Transformed Gauss-Chebyshev quadrature formula of the second kind ***
!     *** u [-1,+1] -> r [0,infinity] => r = ln(2/(1 - u))/ln(2)            ***

      DO i=1,n
        t = REAL(i,dp)*f
        cost = COS(t)
        sint = SIN(t)
        sint2 = sint**2
        x = REAL(2*i - n - 1,dp)/REAL(n + 1,dp) -&
            2.0_dp*(1.0_dp + 2.0_dp*sint2/3.0_dp)*cost*sint/pi
        w = 16.0_dp*sint2**2/REAL(3*(n + 1),dp)
        r(n+1-i) = LOG(2.0_dp/(1.0_dp - x))/LOG(2.0_dp)
        r2(n+1-i) = r(n+1-i)**2
        wr(n+1-i) = w*r2(n+1-i)/(LOG(2.0_dp)*(1.0_dp - x))
      END DO

    CASE DEFAULT

      CALL stop_program("SUBROUTINE radial_grid (MODULE vxc_acsgrid)",&
                        "Invalid radial quadrature type specified")

    END SELECT

  END SUBROUTINE radial_grid

! *****************************************************************************

END MODULE qs_grid_atom

