!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Setup and Methods for semi-empirical multipole types
!> \author Teodoro Laino [tlaino] - 08.2008 Zurich University
! *****************************************************************************
MODULE semi_empirical_mpole_methods
  
  USE input_constants,                 ONLY: do_method_pnnl
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
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
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
!> \param mpoles ...
!> \param se_parameter ...
!> \param method ...
!> \param error ...
!> \date   09.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE semi_empirical_mpole_p_setup(mpoles, se_parameter, method, error)
    TYPE(semi_empirical_mpole_p_type), &
      DIMENSION(:), POINTER                  :: mpoles
    TYPE(semi_empirical_type), POINTER       :: se_parameter
    INTEGER, INTENT(IN)                      :: method
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_mpole_p_setup', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=3), DIMENSION(9), PARAMETER :: label_print_orb = (/"  s",&
      " px"," py"," pz","dx2","dzx","dz2","dzy","dxy"/)
    INTEGER, DIMENSION(9), PARAMETER :: loc_index = (/1,2,2,2,3,3,3,3,3/)

    INTEGER                                  :: a, b, i, ind1, ind2, j, k, &
                                                k1, k2, mu, natorb, ndim, nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dlm, tmp, wp, ws, zb, ZP, ZS, &
                                                zt
    REAL(KIND=dp), DIMENSION(3, 3, 45)       :: M2
    REAL(KIND=dp), DIMENSION(3, 45)          :: M1
    REAL(KIND=dp), DIMENSION(45)             :: M0
    REAL(KIND=dp), DIMENSION(6, 0:2)         :: amn
    TYPE(semi_empirical_mpole_type), POINTER :: mpole

    failure = .FALSE.
    CPPrecondition(.NOT.ASSOCIATED(mpoles),cp_failure_level,routineP,error,failure)
    ! If there are atomic orbitals proceed with the expansion in multipoles
    natorb = se_parameter%natorb
    IF (natorb/=0) THEN
       ndim   = natorb*(natorb+1)/2
       CALL semi_empirical_mpole_p_create(mpoles, ndim, error)

       ! Select method for multipolar expansion

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

             IF (debug_this_module) THEN
                WRITE(*,'(A,2I6,A)')"Orbitals ",i,j,&
                     " ("//label_print_orb(i)//","//label_print_orb(j)//")"
                IF(mpole%task(1)) WRITE(*,'(9F12.6)')mpole%c
                IF(mpole%task(2)) WRITE(*,'(9F12.6)')mpole%d
                IF(mpole%task(3)) WRITE(*,'(9F12.6)')mpole%qc
                WRITE(*,*)
             END IF
          END DO
       END DO

       IF (method == do_method_pnnl) THEN
          ! No d-function for Schenter type integrals
          CPPostcondition(natorb<=4,cp_failure_level,routineP,error,failure)

          M0 = 0.0_dp
          M1 = 0.0_dp
          M2 = 0.0_dp

          DO mu=1,se_parameter%natorb
             M0(indexb(mu,mu))= 1.0_dp
          ENDDO

          ZS=se_parameter%sto_exponents(0)
          ZP=se_parameter%sto_exponents(1)
          nr=se_parameter%nr

          ws=REAL((2*nr+2)*(2*nr+1),dp)/(24.0_dp*ZS**2)
          DO k=1,3
             M2(k,k,indexb(1,1))= ws
          ENDDO

          IF (ZP > 0._dp) THEN
             zt=SQRT(ZS*ZP)
             zb=0.5_dp*(ZS+ZP)
             DO k=1,3
                M1(k, indexb(1,1+k))=(zt/zb)**(2*nr+1)*REAL(2*nr+1,dp)/(2.0*zb*SQRT(3.0_dp))
             ENDDO

             wp=REAL((2*nr+2)*(2*nr+1),dp)/(40.0_dp*ZP**2)
             DO k1=1,3
                DO k2=1,3
                   IF(k1 == k2)THEN
                      M2(k2,k2,indexb(1+k1,1+k1))= 3.0_dp*wp
                   ELSE
                      M2(k2,k2,indexb(1+k1,1+k1))= wp
                   ENDIF
                ENDDO
             ENDDO
             M2(1,2,indexb(1+1,1+2))= wp
             M2(2,1,indexb(1+1,1+2))= wp
             M2(2,3,indexb(1+2,1+3))= wp
             M2(3,2,indexb(1+2,1+3))= wp
             M2(3,1,indexb(1+3,1+1))= wp
             M2(1,3,indexb(1+3,1+1))= wp
          ENDIF

          DO i = 1, natorb
             DO j = 1, i
                ind1 = indexa(se_map_alm(i),se_map_alm(j))
                ind2 = indexb(i,j)
                mpole => mpoles(ind2)%mpole
                mpole%indi = i
                mpole%indj = j
                ! Charge
                mpole%cs      = -M0(indexb(i,j))
                ! Dipole
                mpole%ds      = -M1(1:3, indexb(i,j))
                ! Quadrupole
                mpole%qq      = -3._dp*M2(1:3,1:3, indexb(i,j))
                IF (debug_this_module) THEN
                   WRITE(*,'(A,2I6,A)')"Orbitals ",i,j,&
                        " ("//label_print_orb(i)//","//label_print_orb(j)//")"
                   WRITE(*,'(9F12.6)')mpole%cs
                   WRITE(*,'(9F12.6)')mpole%ds
                   WRITE(*,'(9F12.6)')mpole%qq
                   WRITE(*,*)
                END IF
             END DO
          END DO
       ELSE
          mpole%cs = mpole%c
          mpole%ds = mpole%d
          mpole%qq = mpole%qc
       END IF
    END IF

  END SUBROUTINE semi_empirical_mpole_p_setup

! *****************************************************************************
!> \brief  Transforms the quadrupole components from sphericals to cartesians
!> \param qcart ...
!> \param qsph ...
!> \param error ...
!> \date   09.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
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
    ! Notation
    !          qs(1) - d3z2-r2
    !          qs(2) - dzx
    !          qs(3) - dzy
    !          qs(4) - dx2-y2
    !          qs(5) - dxy
    ! Cartesian components
    qcart(1,1) =   (qsph(4) - qsph(1) / SQRT(3.0_dp))*SQRT(3.0_dp)/2.0_dp
    qcart(2,1) =    qsph(5)*SQRT(3.0_dp)/2.0_dp
    qcart(3,1) =    qsph(2)*SQRT(3.0_dp)/2.0_dp
    qcart(2,2) = - (qsph(4) + qsph(1) / SQRT(3.0_dp))*SQRT(3.0_dp)/2.0_dp
    qcart(3,2) =    qsph(3)*SQRT(3.0_dp)/2.0_dp
    qcart(3,3) =    qsph(1)
    ! Symmetrize tensor
    qcart(1,2) = qcart(2,1)
    qcart(1,3) = qcart(3,1)
    qcart(2,3) = qcart(3,2)

  END SUBROUTINE quadrupole_sph_to_cart

! *****************************************************************************
!> \brief Setup NDDO multipole type
!> \param nddo_mpole ...
!> \param natom ...
!> \param error ...
!> \date   09.2008
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE nddo_mpole_setup(nddo_mpole, natom, error)
    TYPE(nddo_mpole_type), POINTER           :: nddo_mpole
    INTEGER, INTENT(IN)                      :: natom
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nddo_mpole_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

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

    CALL timestop(handle)

  END SUBROUTINE nddo_mpole_setup

END MODULE semi_empirical_mpole_methods
