!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/spherical_gaussians [1.0] *
!!
!!   NAME
!!     spherical_gaussians
!!
!!   FUNCTION
!!
!!   Normalized spherical Gaussian functions
!!
!!               X(lm,p) = N(lp) r^l exp(-p*r^2) Y_lm(theta,phi)
!!
!!               N(lp) = 2^(l+2) [(2l+1)!!]^(-1/2) (2Pi)^(-1/4) p^(l/2+3/4)
!!
!!   AUTHOR
!!     JGH : 16-06-2002
!!
!!   MODIFICATION HISTORY
!!     NONE
!!
!!*****
!******************************************************************************

MODULE spherical_gaussians

  USE ai_spherical_gaussian,           ONLY: sg_kinetic,&
                                             sg_nuclear,&
                                             sg_overlap
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  TYPE exponent_set_type
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: exponents
  END TYPE exponent_set_type

  TYPE sg_basis_type
     INTEGER :: lmax
     INTEGER, DIMENSION ( : ), POINTER :: nl
     TYPE(exponent_set_type), DIMENSION ( : ), POINTER :: expset
     TYPE(all_potential_type), POINTER :: allpot
     TYPE(gth_potential_type), POINTER :: gthpot
  END TYPE sg_basis_type

  TYPE sg_opmatrix_type
     REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: opmatrix
  END TYPE sg_opmatrix_type

  TYPE sg_operator_type
     CHARACTER ( LEN=20 ) :: nametag
     TYPE(sg_basis_type), POINTER :: basis
     LOGICAL :: spherical
     TYPE(sg_opmatrix_type), DIMENSION ( :, : ), POINTER :: opmat
  END TYPE sg_operator_type

  PRIVATE

!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE sg_construct_operator ( op, basis, spherical )
    TYPE(sg_operator_type), INTENT(out)      :: op
    TYPE(sg_basis_type), INTENT(in), TARGET  :: basis
    LOGICAL, INTENT(in)                      :: spherical

    INTEGER                                  :: ierr, k, l, lm, lmax, m, n

   op%nametag = ""
   op%basis => basis
   op%spherical = spherical
   lmax = basis%lmax
   IF ( spherical ) THEN
      lm = 0
   ELSE
      lm = lmax
   END IF
   ALLOCATE ( op%opmat(0:lmax,0:lm), STAT=ierr )
   IF ( ierr /= 0 ) CALL stop_memory ( "sg_construct_operator", &
       "op%opmat", (lmax+1)*(lm+1) )
   DO l = 0, lmax
      n = basis%nl(l)
      DO k = 0, lm
         m = basis%nl(k)
         IF ( spherical ) m = n
         ALLOCATE ( op%opmat(l,k)%opmatrix(n,m), STAT=ierr )
         IF ( ierr /= 0 ) CALL stop_memory ( "sg_construct_operator", &
            "op%opmat%opmatrix", n*m )
      END DO
   END DO

END SUBROUTINE sg_construct_operator

SUBROUTINE sg_deconstruct_operator ( op )
    TYPE(sg_operator_type), INTENT(inout)    :: op

    INTEGER                                  :: ierr, k, l, m, n

   n = SIZE ( op%opmat, 1 )
   m = SIZE ( op%opmat, 2 )

   DO l = 0, n-1
      DO k = 0, m-1
         DEALLOCATE ( op%opmat(l,k)%opmatrix, STAT=ierr )
         IF ( ierr /= 0 ) CALL stop_memory ( "sg_deconstruct_operator", &
            "op%opmat%opmatrix" )
      END DO
   END DO
   DEALLOCATE ( op%opmat, STAT=ierr )
   IF ( ierr /= 0 ) CALL stop_memory ( "sg_deconstruct_operator", "op%opmat" )

   NULLIFY ( op%basis )

END SUBROUTINE sg_deconstruct_operator

!******************************************************************************

SUBROUTINE sg_calculate_operator ( optype, op )
    CHARACTER(LEN=*), INTENT(in)             :: optype
    TYPE(sg_operator_type), INTENT(inout)    :: op

    INTEGER                                  :: ierr, k, l, lm, lmax, n, zval
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)  :: ops, opx
    TYPE(sg_basis_type), POINTER             :: basis

   basis => op%basis
   lmax = basis%lmax

   DO l = 0, lmax

     IF ( op%spherical ) THEN
       lm = 0
     ELSE
       lm = l
     END IF

     SELECT CASE ( optype )
     CASE DEFAULT
     CASE ( "OVERLAP" )
       CALL sg_overlap ( op%opmat(l,lm)%opmatrix, l, basis%expset(l)%exponents )
     CASE ( "KINETIC" )
       CALL sg_kinetic ( op%opmat(l,lm)%opmatrix, l, basis%expset(l)%exponents )
     CASE ( "EXTERNAL" )
       IF ( ASSOCIATED ( basis%allpot ) ) THEN
         CALL get_potential ( basis%allpot, z=zval )
         CALL sg_nuclear ( op%opmat(l,lm)%opmatrix, l, basis%expset(l)%exponents )
         op%opmat(l,lm)%opmatrix = -REAL(zval,dp) * op%opmat(l,lm)%opmatrix
       ELSE IF ( ASSOCIATED ( basis%gthpot ) ) THEN
         stop
       ELSE
         CALL stop_program ( "sg_calculate_operator", "external potential not specified" )
       END IF
     CASE ( "CORE" )
       n = SIZE ( op%opmat(l,lm)%opmatrix, 1 )
       ALLOCATE ( ops ( n, n ), STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "sg_calculate_operator", "ops", n*n )
       ALLOCATE ( opx ( n, n ), STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "sg_calculate_operator", "opx", n*n )
       CALL sg_overlap ( ops, l, basis%expset(l)%exponents )
       CALL sg_kinetic ( opx, l, basis%expset(l)%exponents, smat=ops )
       op%opmat(l,lm)%opmatrix = opx
       IF ( ASSOCIATED ( basis%allpot ) ) THEN
         CALL get_potential ( basis%allpot, z=zval )
         CALL sg_nuclear ( opx, l, basis%expset(l)%exponents, smat=ops )
         op%opmat(l,lm)%opmatrix = op%opmat(l,lm)%opmatrix - REAL(zval,dp) * opx
       ELSE IF ( ASSOCIATED ( basis%gthpot ) ) THEN
         stop
       ELSE
         CALL stop_program ( "sg_calculate_operator", "external potential not specified" )
       END IF
       DEALLOCATE ( ops, STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "sg_calculate_operator", "ops" )
       DEALLOCATE ( opx, STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "sg_calculate_operator", "opx" )
     CASE ( "XC" )
     CASE ( "COULOMB" )
     END SELECT
     
   END DO

   IF ( .NOT. op%spherical ) THEN
     DO l = 0, lmax
       DO k = l+1, lmax
         SELECT CASE ( optype )
         CASE DEFAULT
         CASE ( "XC" )
         CASE ( "COULOMB" )
         END SELECT
       END DO
     END DO
   END IF

END SUBROUTINE sg_calculate_operator

!******************************************************************************

END MODULE spherical_gaussians

!******************************************************************************
