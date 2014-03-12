!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Defines functions to perform rmsd in 3D
!> \author Teodoro Laino 09.2006
! *****************************************************************************
MODULE rmsd
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all,&
                                             matvec_3x3
  USE particle_types,                  ONLY: particle_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rmsd'
  REAL(KIND=dp), PARAMETER, PRIVATE    :: Epsi = EPSILON(0.0_dp)*1.0E4_dp

  PUBLIC :: rmsd3

CONTAINS

! *****************************************************************************
!> \brief Computes the RMSD in 3D. Provides also derivatives.
!> \note
!>      Optional arguments:
!>          my_val -> gives back the value of the RMSD
!>          transl -> provides the translational vector
!>          rotate -> if .true. gives back in r the coordinates rotated
!>                    in order to minimize the RMSD3
!>          rot    -> provides the rotational matrix
!>          drmsd3 -> derivatives of RMSD3 w.r.t. atomic positions
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE rmsd3( particle_set, r, r0, output_unit, weights, my_val,&
       rotate, transl, rot, drmsd3, error)
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: r, r0
    INTEGER, INTENT(IN)                      :: output_unit
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: weights
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: my_val
    LOGICAL, INTENT(IN), OPTIONAL            :: rotate
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(OUT), OPTIONAL                  :: transl
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: rot, drmsd3
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rmsd3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ix, j, k, natom, stat
    LOGICAL                                  :: failure, my_rotate
    REAL(KIND=dp) :: dm_r(4,4,3), lambda(4), loc_tr(3), M(4,4), mtot, &
      my_rot(3,3), Q(0:3), rr(3), rr0(3), rrsq, s, xx, yy, Z(4,4), zz
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: w
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r0p, rp

    failure   = .FALSE.
    IF (.NOT.failure) THEN
       CPPostcondition(SIZE(r)==SIZE(r0),cp_failure_level,routineP,error,failure)
       natom     = SIZE(particle_set)
       my_rotate = .FALSE.
       IF (PRESENT(rotate)) my_rotate=rotate
       ALLOCATE(w(natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (PRESENT(weights)) THEN
          w = weights
       ELSE
          ! All atoms have a weight proportional to their mass in the RMSD unless
          ! differently requested
          DO i=1,natom
             w(i) = particle_set(i)%atomic_kind%mass
          END DO
          mtot = MINVAL(w)
          IF (mtot /= 0.0_dp) w = w/mtot
       END IF
       ALLOCATE(rp(3,natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(r0p(3,natom),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Molecule given by coordinates R
       ! Find COM and center molecule in COM
       xx=0.0_dp
       yy=0.0_dp
       zz=0.0_dp
       mtot=0.0_dp
       DO i=1,natom
          mtot = mtot +particle_set(i)%atomic_kind%mass
          xx=xx+r((i-1)*3+1)*particle_set(i)%atomic_kind%mass
          yy=yy+r((i-1)*3+2)*particle_set(i)%atomic_kind%mass
          zz=zz+r((i-1)*3+3)*particle_set(i)%atomic_kind%mass
       END DO
       xx=xx/mtot
       yy=yy/mtot
       zz=zz/mtot
       DO i=1,natom
          rp(1,i)=r((i-1)*3+1)-xx
          rp(2,i)=r((i-1)*3+2)-yy
          rp(3,i)=r((i-1)*3+3)-zz
       END DO
       IF (PRESENT(transl)) THEN
          transl(1)=xx
          transl(2)=yy
          transl(3)=zz
       END IF
       ! Molecule given by coordinates R0
       ! Find COM and center molecule in COM
       xx=0.0_dp
       yy=0.0_dp
       zz=0.0_dp
       DO i=1,natom
          xx=xx+r0((i-1)*3+1)*particle_set(i)%atomic_kind%mass
          yy=yy+r0((i-1)*3+2)*particle_set(i)%atomic_kind%mass
          zz=zz+r0((i-1)*3+3)*particle_set(i)%atomic_kind%mass
       END DO
       xx=xx/mtot
       yy=yy/mtot
       zz=zz/mtot
       DO i=1,natom
          r0p(1,i)=r0((i-1)*3+1)-xx
          r0p(2,i)=r0((i-1)*3+2)-yy
          r0p(3,i)=r0((i-1)*3+3)-zz
       END DO
       loc_tr(1) = xx
       loc_tr(2) = yy
       loc_tr(3) = zz
       ! Give back the translational vector
       IF (PRESENT(transl)) THEN
          transl(1)=transl(1)-xx
          transl(2)=transl(2)-yy
          transl(3)=transl(3)-zz
       END IF
       M = 0.0_dp
       !
       DO i = 1, natom
          IF (w(i)==0.0_dp) CYCLE
          rr (1)=rp (1,I)
          rr (2)=rp (2,I)
          rr (3)=rp (3,I)
          rr0(1)=r0p(1,I)
          rr0(2)=r0p(2,I)
          rr0(3)=r0p(3,I)
          rrsq  =w(I)*(rr0(1)**2+rr0(2)**2+rr0(3)**2+rr(1)**2+rr(2)**2+rr(3)**2)
          rr0(1)=w(I)*rr0(1)
          rr0(2)=w(I)*rr0(2)
          rr0(3)=w(I)*rr0(3)
          M(1,1)=M(1,1) + rrsq + 2.0_dp*(-rr0(1)*rr(1)-rr0(2)*rr(2)-rr0(3)*rr(3))
          M(2,2)=M(2,2) + rrsq + 2.0_dp*(-rr0(1)*rr(1)+rr0(2)*rr(2)+rr0(3)*rr(3))
          M(3,3)=M(3,3) + rrsq + 2.0_dp*( rr0(1)*rr(1)-rr0(2)*rr(2)+rr0(3)*rr(3))
          M(4,4)=M(4,4) + rrsq + 2.0_dp*( rr0(1)*rr(1)+rr0(2)*rr(2)-rr0(3)*rr(3))
          M(1,2)=M(1,2) + 2.0_dp*(-rr0(2)*rr(3)+rr0(3)*rr(2))
          M(1,3)=M(1,3) + 2.0_dp*( rr0(1)*rr(3)-rr0(3)*rr(1))
          M(1,4)=M(1,4) + 2.0_dp*(-rr0(1)*rr(2)+rr0(2)*rr(1))
          M(2,3)=M(2,3) - 2.0_dp*( rr0(1)*rr(2)+rr0(2)*rr(1))
          M(2,4)=M(2,4) - 2.0_dp*( rr0(1)*rr(3)+rr0(3)*rr(1))
          M(3,4)=M(3,4) - 2.0_dp*( rr0(2)*rr(3)+rr0(3)*rr(2))
       END DO
       ! Symmetrize
       M(2,1) = M(1,2)
       M(3,1) = M(1,3)
       M(3,2) = M(2,3)
       M(4,1) = M(1,4)
       M(4,2) = M(2,4)
       M(4,3) = M(3,4)
       ! Solve the eigenvalue problem for M
       Z = M
       CALL diamat_all(Z,lambda,error=error)
       ! Pick the correct eigenvectors
       S =  1.0_dp
       IF (Z(1,1) .LT. 0.0_dp) S = -1.0_dp
       Q(0) = S*Z(1,1)
       Q(1) = S*Z(2,1)
       Q(2) = S*Z(3,1)
       Q(3) = S*Z(4,1)
       IF (PRESENT(my_val)) THEN
          IF (ABS(lambda(1))<epsi) THEN
             my_val = 0.0_dp
          ELSE
             my_val = lambda(1)/REAL(natom,KIND=dp)
          ENDIF
       END IF
       IF (ABS(lambda(1)-lambda(2))<epsi) THEN
          IF (output_unit>0) WRITE(output_unit,FMT='(T2,"RMSD3|",A)')&
               'NORMAL EXECUTION, NON-UNIQUE SOLUTION'
       END IF
       ! Computes derivatives w.r.t. the positions if requested
       IF (PRESENT(drmsd3)) THEN
          DO I = 1, natom
             IF (W(I)==0.0_dp) CYCLE
             rr(1)  = W(I)*2.0_dp*rp (1,I)
             rr(2)  = W(I)*2.0_dp*rp (2,I)
             rr(3)  = W(I)*2.0_dp*rp (3,I)
             rr0(1) = W(I)*2.0_dp*r0p(1,I)
             rr0(2) = W(I)*2.0_dp*r0p(2,I)
             rr0(3) = W(I)*2.0_dp*r0p(3,I)

             dm_r (1,1,1)=(rr(1)-rr0(1))
             dm_r (1,1,2)=(rr(2)-rr0(2))
             dm_r (1,1,3)=(rr(3)-rr0(3))

             dm_r (1,2,1)= 0.0_dp
             dm_r (1,2,2)= rr0(3)
             dm_r (1,2,3)=-rr0(2)

             dm_r (1,3,1)=-rr0(3)
             dm_r (1,3,2)= 0.0_dp
             dm_r (1,3,3)= rr0(1)

             dm_r (1,4,1)= rr0(2)
             dm_r (1,4,2)=-rr0(1)
             dm_r (1,4,3)= 0.0_dp

             dm_r (2,2,1)=(rr(1)-rr0(1))
             dm_r (2,2,2)=(rr(2)+rr0(2))
             dm_r (2,2,3)=(rr(3)+rr0(3))

             dm_r (2,3,1)=-rr0(2)
             dm_r (2,3,2)=-rr0(1)
             dm_r (2,3,3)= 0.0_dp

             dm_r (2,4,1)=-rr0(3)
             dm_r (2,4,2)= 0.0_dp
             dm_r (2,4,3)=-rr0(1)

             dm_r (3,3,1)=(rr(1)+rr0(1))
             dm_r (3,3,2)=(rr(2)-rr0(2))
             dm_r (3,3,3)=(rr(3)+rr0(3))

             dm_r (3,4,1)= 0.0_dp
             dm_r (3,4,2)=-rr0(3)
             dm_r (3,4,3)=-rr0(2)

             dm_r (4,4,1)=(rr(1)+rr0(1))
             dm_r (4,4,2)=(rr(2)+rr0(2))
             dm_r (4,4,3)=(rr(3)-rr0(3))

             DO ix=1,3
                dm_r(2,1,ix) = dm_r(1,2,ix)
                dm_r(3,1,ix) = dm_r(1,3,ix)
                dm_r(4,1,ix) = dm_r(1,4,ix)
                dm_r(3,2,ix) = dm_r(2,3,ix)
                dm_r(4,2,ix) = dm_r(2,4,ix)
                dm_r(4,3,ix) = dm_r(3,4,ix)
             END DO
             !
             DO ix=1,3
                drmsd3 (ix,I)=0.0_dp
                DO k=1,4
                   DO j=1,4
                      drmsd3(ix,i) = drmsd3(ix,i)+Q(K-1)*Q(j-1)*dm_r(j,k,ix)
                   END DO
                END DO
                drmsd3 (ix,I)=drmsd3(ix,I)/REAL(natom,KIND=dp)
             END DO
          END DO
       END IF
       ! Computes the rotation matrix in terms of quaternions
       my_rot(1,1)=-2.0_dp*Q(2)**2 - 2.0_dp*Q(3)**2 + 1.0_dp
       my_rot(1,2)= 2.0_dp*(-Q(0)*Q(3) + Q(1)*Q(2))
       my_rot(1,3)= 2.0_dp*( Q(0)*Q(2) + Q(1)*Q(3))
       my_rot(2,1)= 2.0_dp*( Q(0)*Q(3) + Q(1)*Q(2))
       my_rot(2,2)=-2.0_dp*Q(1)**2 - 2.0_dp*Q(3)**2 + 1.0_dp
       my_rot(2,3)= 2.0_dp*(-Q(0)*Q(1) + Q(2)*Q(3))
       my_rot(3,1)= 2.0_dp*(-Q(0)*Q(2) + Q(1)*Q(3))
       my_rot(3,2)= 2.0_dp*( Q(0)*Q(1) + Q(2)*Q(3))
       my_rot(3,3)=-2.0_dp*Q(1)**2 - 2.0_dp*Q(2)**2 + 1.0_dp
       IF (PRESENT(rot)) rot = my_rot
       ! Give back coordinates rotated in order to minimize the RMSD
       IF (my_rotate) THEN
          DO i = 1,natom
             CALL matvec_3x3(r((i-1)*3+1:(i-1)*3+3),TRANSPOSE(my_rot),rp(:,i))
             r((i-1)*3+1:(i-1)*3+3) = r((i-1)*3+1:(i-1)*3+3) + loc_tr
          END DO
       END IF
       DEALLOCATE(w,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(rp,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DEALLOCATE(r0p,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE rmsd3

END MODULE rmsd
