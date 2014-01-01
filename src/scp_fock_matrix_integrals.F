!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides the low level routines to build  Coulomb Fock matrices
!>        for SCP
!> \author CJM [tlaino] (12.2010)
!> \par History
! *****************************************************************************
MODULE scp_fock_matrix_integrals
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_arrays,       ONLY: se_orbital_pointer
  USE semi_empirical_integrals,        ONLY: drotint,&
                                             drotnuc,&
                                             rotint,&
                                             rotnuc
  USE semi_empirical_store_int_types,  ONLY: semi_empirical_si_type
  USE semi_empirical_types,            ONLY: se_int_control_type,&
                                             se_taper_type,&
                                             semi_empirical_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_fock_matrix_integrals'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.

  PUBLIC :: scp_1el, dscp_1el, scp_1el_ew, scp2C_ew, scp2C, dscp2C

CONTAINS

! *****************************************************************************
!> \brief  Construction of 2-center 1-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \param  DIMENSION(sepi%natorb, sepi%natorb) :: ksi_block
!> \param  DIMENSION(sepj%natorb, sepj%natorb) :: ksj_block
!> \date   04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE scp_1el (sepi, sepj, rij, ksi_block, ksj_block, pi_block, pj_block,&
       ecore, itype, anag, se_int_control, se_taper, store_int_env, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ksi_block, ksj_block
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), INTENT(INOUT)             :: ecore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(45)             :: e1b, e2a

    failure = .FALSE.
    ! Compute integrals
!    WRITE ( *, * ) 'SCP INTEGRAL CONTROL', (se_int_control%integral_screening==do_se_IS_slater)
    CALL rotnuc (sepi, sepj, rij, e1b=e1b, e2a=e2a, itype=itype, anag=anag,&
         se_int_control=se_int_control, se_taper=se_taper, store_int_env=store_int_env,&
         error=error)
    !
    ! Add the electron-nuclear attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          ksi_block(i1,j1) = ksi_block(i1,j1) + e1b(i2)
          ksi_block(j1,i1) = ksi_block(i1,j1)
          ecore = ecore + 2.0_dp * e1b(i2) * pi_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       ksi_block(i1,j1) = ksi_block(i1,j1) + e1b(i2)
       ecore = ecore + e1b(i2) * pi_block(i1,j1)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          ksj_block(i1,j1) = ksj_block(i1,j1) + e2a(i2)
          ksj_block(j1,i1) = ksj_block(i1,j1)
          ecore = ecore + 2.0_dp * e2a(i2) * pj_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       ksj_block(i1,j1) = ksj_block(i1,j1) + e2a(i2)
       ecore = ecore + e2a(i2) * pj_block(i1,j1)
    END DO

  END SUBROUTINE scp_1el

! *****************************************************************************
!> \brief Derivatives of 2-center 1-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dscp_1el (sepi, sepj, rij, pi_block, pj_block, itype, anag,&
       se_int_control, se_taper, force, delta, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dscp_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), DIMENSION(3, 45)          :: de1b, de2a

    failure = .FALSE.
    ! Compute integrals
    CALL drotnuc (sepi, sepj, rij, de1b=de1b, de2a=de2a, itype=itype, anag=anag,&
         se_int_control=se_int_control, se_taper=se_taper, delta=delta, error=error)
    !
    ! Add the electron-nuclear attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          tmp= 2.0_dp * pi_block(i1,j1)
          force(1) = force(1) + de1b(1,i2) * tmp
          force(2) = force(2) + de1b(2,i2) * tmp
          force(3) = force(3) + de1b(3,i2) * tmp
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       force(1) = force(1) + de1b(1,i2) * pi_block(i1,j1)
       force(2) = force(2) + de1b(2,i2) * pi_block(i1,j1)
       force(3) = force(3) + de1b(3,i2) * pi_block(i1,j1)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          tmp= 2.0_dp * pj_block(i1,j1)
          force(1) = force(1) + de2a(1,i2) * tmp
          force(2) = force(2) + de2a(2,i2) * tmp
          force(3) = force(3) + de2a(3,i2) * tmp
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       force(1) = force(1) + de2a(1,i2) * pj_block(i1,j1)
       force(2) = force(2) + de2a(2,i2) * pj_block(i1,j1)
       force(3) = force(3) + de2a(3,i2) * pj_block(i1,j1)
    END DO

  END SUBROUTINE dscp_1el

! *****************************************************************************
!> \brief Construction of 2-center 1-electron Fock Matrix (Ewald self term)
!> \author jgh - University of Zurich
!> \param  DIMENSION(sep%natorb, sep%natorb) :: ks_block
!> \date 04.2009 [jgh]
! *****************************************************************************
  SUBROUTINE scp_1el_ew (sep, rij, ks_block, p_block, ecore, itype, anag, &
                           se_int_control, se_taper, store_int_env, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ks_block
    REAL(KIND=dp), &
      DIMENSION(sep%natorb, sep%natorb), &
      INTENT(IN)                             :: p_block
    REAL(KIND=dp), INTENT(INOUT)             :: ecore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_1el_ew', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L, n
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(45)             :: e1b, e2a

    failure = .FALSE.
    ! Compute integrals
    CALL rotnuc (sep, sep, rij, e1b=e1b, e2a=e2a, itype=itype, anag=anag,&
         se_int_control=se_int_control, se_taper=se_taper, store_int_env=store_int_env,&
         error=error)
    !
    ! Add the electron-nuclear attraction term for atom sep
    ! e1b == e2a
    !
    n = (sep%natorb*(sep%natorb+1))/2
    i2 = 0
    DO i1L = 1, sep%natorb
       i1 = se_orbital_pointer(i1L)
       DO j1L = 1, i1L-1
          j1 = se_orbital_pointer(j1L)
          i2 = i2 + 1
          ks_block(i1,j1) = ks_block(i1,j1) + e1b(i2)
          ks_block(j1,i1) = ks_block(i1,j1)
          ecore = ecore + 2._dp * e1b(i2) * p_block(i1,j1)
       END DO
       ! i1L == j1L
       i2 = i2 + 1
       ks_block(i1,i1) = ks_block(i1,i1) + e1b(i2)
       ecore = ecore + e1b(i2) * p_block(i1,i1)
    END DO

  END SUBROUTINE scp_1el_ew

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Self Terms (Ewald)
!> \author jgh - University of Zurich
!> \param  DIMENSION(sep%natorb, sep%natorb) :: f_mat
!> \date 04.2009 [jgh]
! *****************************************************************************
  SUBROUTINE scp2C_ew(sep, rij, p_tot, f_mat, factor, anag, se_int_control, &
                       se_taper, store_int_env, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: p_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: f_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp2C_ew', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, aa, bb
    REAL(KIND=dp), DIMENSION(2025)           :: w

    failure = .FALSE.
    ! Evaluate integrals
    CALL rotint (sep,sep,rij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env, error=error)
    kr = 0
    natorb = sep%natorb
    DO iL = 1, natorb
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a = 0.5_dp*w(kr)*factor
                ! Coulomb
                f_mat(i,j) = f_mat(i,j) + bb * a * p_tot(k,l)
                f_mat(k,l) = f_mat(k,l) + aa * a * p_tot(i,j)
                f_mat(j,i) = f_mat(i,j)
                f_mat(l,k) = f_mat(k,l)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE scp2C_ew

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Terms
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \param  DIMENSION(sepi%natorb, sepi%natorb) :: fi_mat
!> \param  DIMENSION(sepj%natorb, sepj%natorb) :: fj_mat
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE scp2C(sepi, sepj, rij, switch, pi_tot, fi_mat, pj_tot, fj_mat, &
       factor, anag, se_int_control, se_taper, store_int_env, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fi_mat
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pj_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fj_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp2C', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, aa, bb, irij(3)
    REAL(KIND=dp), DIMENSION(2025)           :: w

    failure = .FALSE.
    ! Evaluate integrals
    IF (.NOT.switch) THEN
       CALL rotint (sepi,sepj, rij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env, error=error)
    ELSE
       irij = -rij
       CALL rotint (sepj,sepi,irij,w,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,store_int_env=store_int_env, error=error)
    END IF
    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a = w(kr)*factor
                ! Coulomb
                IF (.NOT.switch) THEN
                   fi_mat(i,j) = fi_mat(i,j) + bb * a * pj_tot(k,l)
                   fj_mat(k,l) = fj_mat(k,l) + aa * a * pi_tot(i,j)
                   fi_mat(j,i) = fi_mat(i,j)
                   fj_mat(l,k) = fj_mat(k,l)
                ELSE
                   fj_mat(i,j) = fj_mat(i,j) + bb * a * pi_tot(k,l)
                   fi_mat(k,l) = fi_mat(k,l) + aa * a * pj_tot(i,j)
                   fj_mat(j,i) = fj_mat(i,j)
                   fi_mat(l,k) = fi_mat(k,l)
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE scp2C

! *****************************************************************************
!> \brief Derivatives of 2-center Fock Matrix - Coulomb Terms
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dscp2C(sepi, sepj, rij, switch, pi_tot, pj_tot, factor, anag,&
       se_int_control, se_taper, force, delta, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot, pj_tot
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dscp2C', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, bb, tmp
    REAL(KIND=dp), DIMENSION(3)              :: a, irij
    REAL(KIND=dp), DIMENSION(3, 2025)        :: dw

    failure = .FALSE.
    ! Evaluate integrals' derivatives
    IF (.NOT.switch) THEN
       CALL drotint (sepi,sepj, rij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,error=error)
    ELSE
       irij = -rij
       CALL drotint (sepj,sepi,irij,dw,delta,anag=anag,se_int_control=se_int_control,&
            se_taper=se_taper,error=error)
    END IF

    kr = 0
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
    END IF
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       aa = 2.0_dp
       DO jL = 1, iL
          j = se_orbital_pointer(jL)
          IF (i == j) THEN
             aa = 1.0_dp
          END IF
          DO kL = 1, natorb(2)
             k = se_orbital_pointer(kL)
             bb = 2.0_dp
             DO lL = 1, kL
                l = se_orbital_pointer(lL)
                IF (k == l) THEN
                   bb = 1.0_dp
                END IF
                kr = kr + 1
                a(1) = dw(1,kr)*factor
                a(2) = dw(2,kr)*factor
                a(3) = dw(3,kr)*factor
                ! Coulomb
                IF (.NOT.switch) THEN
                   tmp = bb * aa * pj_tot(k,l) * pi_tot(i,j)
                ELSE
                   tmp = bb * aa * pi_tot(k,l) * pj_tot(i,j)
                END IF
                force(1) = force(1) + a(1) * tmp
                force(2) = force(2) + a(2) * tmp
                force(3) = force(3) + a(3) * tmp
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE dscp2C

END MODULE scp_fock_matrix_integrals

