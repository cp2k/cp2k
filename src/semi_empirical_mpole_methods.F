!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Setup and Methods for semi-empirical multipole types 
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
MODULE semi_empirical_mpole_methods
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_arrays,       ONLY: alm,&
                                             indexa,&
                                             indexb,&
                                             se_map_alm
  USE semi_empirical_mpole_types,      ONLY: nddo_mpole_create,&
                                             nddo_mpole_release,&
                                             nddo_mpole_type,&
                                             semi_empirical_mpole_p_create,&
                                             semi_empirical_mpole_p_type,&
                                             semi_empirical_mpole_type
  USE semi_empirical_par_utils,        ONLY: amn_l
  USE semi_empirical_types,            ONLY: semi_empirical_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_mpole_methods'

  PUBLIC :: semi_empirical_mpole_p_setup,&
            nddo_mpole_setup,&
            quadrupole_sph_to_cart

CONTAINS

! *****************************************************************************
!> \brief Setup semi-empirical mpole type
!>        This function setup for each semi-empirical type a structure containing
!>        the multipolar expansion for all possible combination on-site of atomic
!>        orbitals ( \mu \nu |
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   09.2008
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_p_setup(mpoles, se_parameter, error)
    TYPE(semi_empirical_mpole_p_type), &
      DIMENSION(:), POINTER                  :: mpoles
    TYPE(semi_empirical_type), POINTER       :: se_parameter
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_mpole_p_setup', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: Nb = 4, Nm = Nb+Nb*(Nb-1)/2
    INTEGER, DIMENSION(Nb, Nb), PARAMETER :: &
      M = RESHAPE( (/1,2,3,4,2,5,8,10,3,8,6,9,4,10,9,7/) , (/4,4/) )
    REAL(KIND=dp), DIMENSION(9), PARAMETER :: &
      loc_index = (/1,2,2,2,3,3,3,3,3/)

    INTEGER                                  :: a, b, i, ind1, ind2, j, k, &
                                                k1, k2, mu, natorb, ndim, nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dlm, tmp, wp, ws, zb, ZP, ZS, &
                                                zt
    REAL(KIND=dp), DIMENSION(6, 0:2)         :: amn
    REAL(KIND=dp), DIMENSION(Nm)             :: M0
    REAL(KIND=dp), DIMENSION(Nm, 3)          :: M1
    REAL(KIND=dp), DIMENSION(Nm, 3, 3)       :: M2
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(mpoles),cp_failure_level,routineP,error,failure)
    ! If there are atomic orbitals proceed with the expansion in multipoles
    natorb = se_parameter%natorb
    IF (natorb/=0) THEN
       ndim   = natorb*(natorb+1)/2
       CALL semi_empirical_mpole_p_create(mpoles, ndim, error)
       
       M0(:)=     0.0_dp
       M1(:,:)=   0.0_dp
       M2(:,:,:)= 0.0_dp

       DO mu=1,Nb
         M0(M(mu,mu))= 1.0_dp
       ENDDO

       ZS=se_parameter%sto_exponents(0)
       ZP=se_parameter%sto_exponents(1)
       nr=se_parameter%nr

       IF (ZP > 0._dp) THEN
         zt=SQRT(ZS*ZP)
         zb=0.5_dp*(ZS+ZP)
         DO k=1,3
           M1(M(1,1+k),k)=(zt/zb)**(2*nr+1)*REAL(2*nr+1,dp)/(2.0*zb*SQRT(3.0_dp))
         ENDDO
  
         ws=REAL((2*nr+2)*(2*nr+1),dp)/(24.0_dp*ZS**2)
         DO k=1,3
           M2(M(1,1),k,k)= ws
         ENDDO
         wp=REAL((2*nr+2)*(2*nr+1),dp)/(40.0_dp*ZP**2)
         DO k1=1,3
           DO k2=1,3
             IF(k1 == k2)THEN
               M2(M(1+k1,1+k1),k2,k2)= 3.0_dp*wp
             ELSE
               M2(M(1+k1,1+k1),k2,k2)= wp
             ENDIF
           ENDDO
         ENDDO
         M2(M(1+1,1+2),1,2)= wp
         M2(M(1+1,1+2),2,1)= wp
         M2(M(1+2,1+3),2,3)= wp
         M2(M(1+2,1+3),3,2)= wp
         M2(M(1+3,1+1),3,1)= wp
         M2(M(1+3,1+1),1,3)= wp
       ENDIF

       ! Fill in information on multipole expansion due to atomic orbitals charge 
       ! distribution
       NULLIFY(mpole)
       CALL amn_l(se_parameter, amn, error)
       DO i = 1, natorb
          DO j = 1, i
             ind1 = indexa(se_map_alm(i),se_map_alm(j))
             ind2 = indexb(i,j)
             ! the order in the mpoles structure is like the standard one for the
             ! integrals: s px py pz dx2-y2 dzx dz2 dzy dxy (lower triangular)
             ! which differs from the order of the Hamiltonian in CP2K. But I 
             ! preferred to keep this order for consistency with the integrals
             mpole => mpoles(ind2)%mpole
             mpole%indi = i
             mpole%indj = j
             a = loc_index(i)
             b = loc_index(j)
             mpole%c  = HUGE(0.0_dp)
             mpole%d  = HUGE(0.0_dp)
             mpole%qs = HUGE(0.0_dp)
             mpole%qc = HUGE(0.0_dp)

             ! Charge
             IF (alm(ind1,0,0)/=0.0_dp) THEN
                dlm = 1.0_dp/SQRT(REAL((2*0+1),KIND=dp))
                tmp = - dlm * amn(indexb(a,b),0)
                mpole%c      = tmp * alm(ind1,0,0)
                mpole%task(1)= .TRUE.
             END IF
             
             ! Dipole
             IF (ANY(alm(ind1,1,-1:1)/=0.0_dp)) THEN
                dlm = 1.0_dp/SQRT(REAL((2*1+1),KIND=dp))
                tmp = - dlm * amn(indexb(a,b),1)
                mpole%d(1)   = tmp * alm(ind1,1, 1)
                mpole%d(2)   = tmp * alm(ind1,1,-1)
                mpole%d(3)   = tmp * alm(ind1,1, 0)
                mpole%task(2)= .TRUE.
             END IF
             
             ! Quadrupole
             IF (ANY(alm(ind1,2,-2:2)/=0.0_dp)) THEN
                dlm = 1.0_dp/SQRT(REAL((2*2+1),KIND=dp))
                tmp = - dlm * amn(indexb(a,b),2)
                
                ! Spherical components
                mpole%qs(1)   = tmp * alm(ind1,2, 0) ! d3z2-r2
                mpole%qs(2)   = tmp * alm(ind1,2, 1) ! dzx 
                mpole%qs(3)   = tmp * alm(ind1,2,-1) ! dzy 
                mpole%qs(4)   = tmp * alm(ind1,2, 2) ! dx2-y2 
                mpole%qs(5)   = tmp * alm(ind1,2,-2) ! dxy

                ! Convert into cartesian components
                CALL quadrupole_sph_to_cart(mpole%qc, mpole%qs, error)
                mpole%task(3)= .TRUE.
             END IF

          END DO
       END DO
    END IF

    ! no d-function for Schenter type integrals
    IF (natorb>0.AND.natorb<=4) THEN
       DO i = 1, natorb
          DO j = 1, i
             ind1 = indexa(se_map_alm(i),se_map_alm(j))
             ind2 = indexb(i,j)
             mpole => mpoles(ind2)%mpole
             !Schenter definition
             mpole%cs = HUGE(0.0_dp)
             mpole%ds = HUGE(0.0_dp)
             mpole%qq = HUGE(0.0_dp)
             mpole%cs      = -M0(M(i,j))
             mpole%ds(:)   = -M1(M(i,j),:)
             mpole%qq(:,:) = -3._dp*M2(M(i,j),:,:)
          END DO
       END DO
    END IF

  END SUBROUTINE semi_empirical_mpole_p_setup

! *****************************************************************************
!> \brief  Transforms the quadrupole components from sphericals to cartesians
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   09.2008
! *****************************************************************************
  SUBROUTINE quadrupole_sph_to_cart(qcart, qsph, error)
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(OUT)                            :: qcart
    REAL(KIND=dp), DIMENSION(5), INTENT(IN)  :: qsph
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'quadrupole_sph_to_cart', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Cartesian components
    qcart(1,1) =    qsph(4) - qsph(1) * 0.5_dp
    qcart(2,1) =    qsph(5)
    qcart(3,1) =    qsph(2)
    qcart(2,2) = - (qsph(4) + qsph(1) * 0.5_dp)
    qcart(3,2) =    qsph(3)
    qcart(3,3) =    qsph(1)
    ! Symmetrize tensor
    qcart(1,2) = qcart(2,1)
    qcart(1,3) = qcart(3,1)
    qcart(2,3) = qcart(3,2)

  END SUBROUTINE quadrupole_sph_to_cart

! *****************************************************************************
!> \brief Setup NDDO multipole type
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date   09.2008
! *****************************************************************************
  SUBROUTINE nddo_mpole_setup(nddo_mpole, natom, error)
    TYPE(nddo_mpole_type), POINTER           :: nddo_mpole
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nddo_mpole_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(nddo_mpole)) THEN
       CALL nddo_mpole_release(nddo_mpole, error=error)
    END IF
    CALL nddo_mpole_create(nddo_mpole, error=error)
    ! Allocate Global Arrays
    ALLOCATE(nddo_mpole%charge(natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nddo_mpole%dipole(3,natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nddo_mpole%quadrupole(3,3,natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(nddo_mpole%efield0(  natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nddo_mpole%efield1(3,natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(nddo_mpole%efield2(9,natom),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END SUBROUTINE nddo_mpole_setup
  
END MODULE semi_empirical_mpole_methods
