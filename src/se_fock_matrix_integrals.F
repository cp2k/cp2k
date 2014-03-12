!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Provides the low level routines to build both the exchange and the
!>        Coulomb Fock matrices.. This routines support d-orbitals and should
!>        be changed only if one knows exactly what he is doing..
!> \author Teodoro Laino [tlaino] (05.2009) - Split and module reorganization
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : d-orbitals
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Speed-up
!>      Teodoro Laino (09.2008) [tlaino] - University of Zurich : Periodic SE
! *****************************************************************************
MODULE se_fock_matrix_integrals
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
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'se_fock_matrix_integrals'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module       = .FALSE.

  PUBLIC :: fock2_1el, dfock2_1el, fock1_2el, fock2_1el_ew, fock2C_ew,&
            fock2C, dfock2C, fock2E, dfock2E, fock2_1el_r3, dfock2_1el_r3,&
            fock2C_r3, dfock2C_r3, se_coulomb_ij_interaction

CONTAINS

! *****************************************************************************
!> \brief  Construction of 2-center 1-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \param  DIMENSION(sepi%natorb, sepi%natorb) :: ksi_block
!> \param  DIMENSION(sepj%natorb, sepj%natorb) :: ksj_block
!> \date   04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2_1el (sepi, sepj, rij, ksi_block, ksj_block, pi_block, pj_block,&
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
    REAL(KIND=dp), DIMENSION(2), &
      INTENT(INOUT)                          :: ecore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2_1el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2, j1, j1L
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(45)             :: e1b, e2a

    failure = .FALSE.
    ! Compute integrals
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
          ecore(1) = ecore(1) + 2.0_dp * e1b(i2) * pi_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       ksi_block(i1,j1) = ksi_block(i1,j1) + e1b(i2)
       ecore(1) = ecore(1) + e1b(i2) * pi_block(i1,j1)
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
          ecore(2) = ecore(2) + 2.0_dp * e2a(i2) * pj_block(i1,j1)
       END DO
       j1 = se_orbital_pointer(j1L)
       i2 = i2 + 1
       ksj_block(i1,j1) = ksj_block(i1,j1) + e2a(i2)
       ecore(2) = ecore(2) + e2a(i2) * pj_block(i1,j1)
    END DO

  END SUBROUTINE fock2_1el

! *****************************************************************************
!> \brief Derivatives of 2-center 1-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2_1el (sepi, sepj, rij, pi_block, pj_block, itype, anag,&
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

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2_1el', &
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

  END SUBROUTINE dfock2_1el

! *****************************************************************************
!> \brief Construction of 1-center 2-electron Fock Matrix
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \param  DIMENSION(sep%natorb, sep%natorb) :: f_mat
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock1_2el(sep, p_tot, p_mat, f_mat, factor, error)
    TYPE(semi_empirical_type), POINTER       :: sep
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: p_tot
    REAL(KIND=dp), &
      DIMENSION(sep%natorb, sep%natorb), &
      INTENT(IN)                             :: p_mat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: f_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock1_2el', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ijw, ikw, iL, im, j, jL, &
                                                jlw, jm, k, kL, klw, l, lL
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: sum

    failure = .FALSE.
    !   One-center coulomb and exchange terms for semiempirical_type sep
    !
    !  F(i,j)=F(i,j)+sum(k,l)((PA(k,l)+PB(k,l))*<i,j|k,l>
    !                        -(PA(k,l)        )*<i,k|j,l>), k,l on type sep.
    !
    DO iL = 1, sep%natorb
       i = se_orbital_pointer(iL)
       DO jL = 1, iL
          j = se_orbital_pointer(jL)

          !    `J' Address IJ in W
          ijw = (iL*(iL-1))/2 + jL
          sum = 0.0_dp
          DO kL = 1, sep%natorb
             k = se_orbital_pointer(kL)
             DO lL =  1, sep%natorb
                l = se_orbital_pointer(lL)

                !    `J' Address KL in W
                im = MAX(kL,lL)
                jm = MIN(kL,lL)
                klw = (im*(im-1))/2 + jm

                !    `K' Address IK in W
                im = MAX(kL,jL)
                jm = MIN(kL,jL)
                ikw = (im*(im-1))/2 + jm

                !    `K' Address JL in W
                im = MAX(lL,iL)
                jm = MIN(lL,iL)
                jlw = (im*(im-1))/2 + jm

                sum = sum + p_tot(k,l) * sep%w(ijw, klw) - p_mat(k,l) * sep%w(ikw, jlw)
             END DO
          END DO
          f_mat(i,j) = f_mat(i,j) + factor*sum
          f_mat(j,i) = f_mat(i,j)
       END DO
    END DO
  END SUBROUTINE fock1_2el

! *****************************************************************************
!> \brief Construction of 2-center 1-electron Fock Matrix (Ewald self term)
!> \author jgh - University of Zurich
!> \param  DIMENSION(sep%natorb, sep%natorb) :: ks_block
!> \date 04.2009 [jgh]
! *****************************************************************************
  SUBROUTINE fock2_1el_ew (sep, rij, ks_block, p_block, ecore, itype, anag, &
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

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2_1el_ew', &
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

  END SUBROUTINE fock2_1el_ew

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Self Terms (Ewald)
!> \author jgh - University of Zurich
!> \param  DIMENSION(sep%natorb, sep%natorb) :: f_mat
!> \date 04.2009 [jgh]
! *****************************************************************************
  SUBROUTINE fock2C_ew(sep, rij, p_tot, f_mat, factor, anag, se_int_control, &
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

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2C_ew', &
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

  END SUBROUTINE fock2C_ew

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Terms
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \param  DIMENSION(sepi%natorb, sepi%natorb) :: fi_mat
!> \param  DIMENSION(sepj%natorb, sepj%natorb) :: fj_mat
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2C(sepi, sepj, rij, switch, pi_tot, fi_mat, pj_tot, fj_mat, &
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

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2C', &
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

  END SUBROUTINE fock2C

! *****************************************************************************
!> \brief Derivatives of 2-center Fock Matrix - Coulomb Terms
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2C(sepi, sepj, rij, switch, pi_tot, pj_tot, factor, anag,&
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

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2C', &
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
  END SUBROUTINE dfock2C

! *****************************************************************************
!> \brief Construction of 2-center Fock Matrix - General Driver
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \param  DIMENSION(isize(1), isize(2)) :: fi_mat
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2E(sepi, sepj, rij, switch, isize, pi_tot, pi_mat, fi_mat, factor,&
       anag, se_int_control, se_taper, store_int_env, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    INTEGER, DIMENSION(2), INTENT(IN)        :: isize
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), &
      DIMENSION(isize(1), isize(2)), &
      INTENT(IN)                             :: pi_mat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fi_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(semi_empirical_si_type), POINTER    :: store_int_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2E', &
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
                ! Exchange
                a = a * aa * bb * 0.25_dp
                fi_mat(i,k) = fi_mat(i,k) - a * pi_mat(j,l)
                fi_mat(i,l) = fi_mat(i,l) - a * pi_mat(j,k)
                fi_mat(j,k) = fi_mat(j,k) - a * pi_mat(i,l)
                fi_mat(j,l) = fi_mat(j,l) - a * pi_mat(i,k)
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE fock2E

! *****************************************************************************
!> \brief Derivatives of 2-center Fock Matrix - General Driver
!> \author Teodoro Laino [tlaino] - University of Zurich
!> \date 04.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2E(sepi, sepj, rij, switch, isize, pi_tot, pi_mat, factor, anag,&
       se_int_control, se_taper, force, delta, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rij
    LOGICAL, INTENT(IN)                      :: switch
    INTEGER, DIMENSION(2), INTENT(IN)        :: isize
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), &
      DIMENSION(isize(1), isize(2)), &
      INTENT(IN)                             :: pi_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    LOGICAL, INTENT(IN)                      :: anag
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), INTENT(IN)                :: delta
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2E', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, j, jL, k, kL, kr, l, &
                                                lL, natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aa, bb, tmp, tmp1, tmp2, &
                                                tmp3, tmp4
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
                kr   = kr + 1
                tmp  = factor * aa * bb * 0.25_dp
                a(1) = dw(1,kr)*tmp
                a(2) = dw(2,kr)*tmp
                a(3) = dw(3,kr)*tmp
                ! Exchange
                tmp1 = pi_mat(j,l) * pi_mat(i,k)
                tmp2 = pi_mat(j,k) * pi_mat(i,l)
                tmp3 = pi_mat(i,l) * pi_mat(j,k)
                tmp4 = pi_mat(i,k) * pi_mat(j,l)

                force(1) = force(1) - a(1) * tmp1
                force(1) = force(1) - a(1) * tmp2
                force(1) = force(1) - a(1) * tmp3
                force(1) = force(1) - a(1) * tmp4

                force(2) = force(2) - a(2) * tmp1
                force(2) = force(2) - a(2) * tmp2
                force(2) = force(2) - a(2) * tmp3
                force(2) = force(2) - a(2) * tmp4

                force(3) = force(3) - a(3) * tmp1
                force(3) = force(3) - a(3) * tmp2
                force(3) = force(3) - a(3) * tmp3
                force(3) = force(3) - a(3) * tmp4
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE dfock2E

! *****************************************************************************
!> \brief  Construction of 2-center 1-electron Fock Matrix for the residual
!>         (1/R^3) integral part
!> \param  DIMENSION(sepi%natorb, sepi%natorb) :: ksi_block
!> \param  DIMENSION(sepj%natorb, sepj%natorb) :: ksj_block
!> \author Teodoro Laino [tlaino]
!> \date   12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2_1el_r3 (sepi, sepj, ksi_block, ksj_block, pi_block, pj_block,&
       e1b, e2a, ecore, rp, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ksi_block, ksj_block
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e1b, e2a
    REAL(KIND=dp), DIMENSION(2), &
      INTENT(INOUT)                          :: ecore
    REAL(KIND=dp), INTENT(IN)                :: rp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2_1el_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2
    LOGICAL                                  :: failure

    failure = .FALSE.
    !
    ! Add the electron-nuclear residual attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i2 = i2 + 1
       i1 = se_orbital_pointer(i1L)
       ksi_block(i1,i1) = ksi_block(i1,i1) + e1b(i2) * rp
       ecore(1) = ecore(1) + e1b(i2) * rp * pi_block(i1,i1)
    END DO
    !
    ! Add the electron-nuclear residual attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i2 = i2 + 1
       i1 = se_orbital_pointer(i1L)
       ksj_block(i1,i1) = ksj_block(i1,i1) + e2a(i2) * rp
       ecore(2) = ecore(2) + e2a(i2) * rp * pj_block(i1,i1)
    END DO

  END SUBROUTINE fock2_1el_r3

! *****************************************************************************
!> \brief  Derivatives of 2-center 1-electron Fock Matrix residual (1/R^3)
!>         integral part
!> \author Teodoro Laino [tlaino]
!> \date   12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2_1el_r3 (sepi, sepj, drp, pi_block, pj_block, force, e1b, e2a,&
       error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: drp
    REAL(KIND=dp), &
      DIMENSION(sepi%natorb, sepi%natorb), &
      INTENT(IN)                             :: pi_block
    REAL(KIND=dp), &
      DIMENSION(sepj%natorb, sepj%natorb), &
      INTENT(IN)                             :: pj_block
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: e1b, e2a
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2_1el_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i1L, i2
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp

    failure = .FALSE.
    !
    ! Add the electron-nuclear residual attraction term for atom sepi
    !
    i2 = 0
    DO i1L = 1, sepi%natorb
       i1 = se_orbital_pointer(i1L)
       i2 = i2 + 1
       tmp = e1b(i2) * pi_block(i1,i1)
       force(1) = force(1) + tmp * drp(1)
       force(2) = force(2) + tmp * drp(2)
       force(3) = force(3) + tmp * drp(3)
    END DO
    !
    ! Add the electron-nuclear attraction term for atom sepj
    !
    i2 = 0
    DO i1L = 1, sepj%natorb
       i1 = se_orbital_pointer(i1L)
       i2 = i2 + 1
       tmp = e2a(i2) * pj_block(i1,i1)
       force(1) = force(1) + tmp * drp(1)
       force(2) = force(2) + tmp * drp(2)
       force(3) = force(3) + tmp * drp(3)
    END DO

  END SUBROUTINE dfock2_1el_r3

! *****************************************************************************
!> \brief  Construction of 2-center Fock Matrix - Coulomb Terms for the residual
!>         (1/R^3) integral part
!> \author Teodoro Laino [tlaino]
!> \param  DIMENSION(sepi%natorb, sepi%natorb) :: fi_mat
!> \param  DIMENSION(sepj%natorb, sepj%natorb) :: fj_mat
!> \date   12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE fock2C_r3(sepi, sepj, switch, pi_tot, fi_mat, pj_tot, fj_mat, &
       factor, w, rp, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fi_mat
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pj_tot
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: fj_mat
    REAL(KIND=dp), INTENT(IN)                :: factor
    REAL(KIND=dp), DIMENSION(81), INTENT(IN) :: w
    REAL(KIND=dp), INTENT(IN)                :: rp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fock2C_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, ind, j, k, kL, kr, &
                                                natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, w_l(81)

    failure = .FALSE.
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
       ! Reshuffle the integral array (natural storage order is sepi/sepj)
       kr = 0
       DO i = 1, sepj%natorb
          DO j = 1, sepi%natorb
             kr  = kr + 1
             ind = (j-1)*sepj%natorb+i
             w_l(kr) = w(ind)
          END DO
       END DO
    ELSE
       w_l = w
    END IF

    ! Modify the Fock Matrix
    kr = 0
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       DO kL = 1, natorb(2)
          k = se_orbital_pointer(kL)
          kr = kr + 1
          a = w_l(kr) * factor * rp
          ! Coulomb
          IF (.NOT.switch) THEN
             fi_mat(i,i) = fi_mat(i,i) + a * pj_tot(k,k)
             fj_mat(k,k) = fj_mat(k,k) + a * pi_tot(i,i)
          ELSE
             fj_mat(i,i) = fj_mat(i,i) + a * pi_tot(k,k)
             fi_mat(k,k) = fi_mat(k,k) + a * pj_tot(i,i)
          END IF
       END DO
    END DO

  END SUBROUTINE fock2C_r3

! *****************************************************************************
!> \brief  Derivatives of 2-center Fock Matrix - Coulomb Terms for the residual
!>         (1/R^3) integral part
!> \author Teodoro Laino [tlaino]
!> \date   12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE dfock2C_r3(sepi, sepj, switch, pi_tot, pj_tot, factor,  w, drp,&
       force, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    LOGICAL, INTENT(IN)                      :: switch
    REAL(KIND=dp), DIMENSION(45, 45), &
      INTENT(IN)                             :: pi_tot, pj_tot
    REAL(KIND=dp), INTENT(IN)                :: factor
    REAL(KIND=dp), DIMENSION(81), INTENT(IN) :: w
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: drp
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dfock2C_r3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iL, ind, j, k, kL, kr, &
                                                natorb(2)
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a(3), tmp, w_l(81)

    failure = .FALSE.
    natorb(1) = sepi%natorb
    natorb(2) = sepj%natorb
    IF (switch) THEN
       natorb(1) = sepj%natorb
       natorb(2) = sepi%natorb
       ! Reshuffle the integral array (natural storage order is sepi/sepj)
       kr = 0
       DO i = 1, sepj%natorb
          DO j = 1, sepi%natorb
             kr  = kr + 1
             ind = (j-1)*sepj%natorb+i
             w_l(kr) = w(ind)
          END DO
       END DO
    ELSE
       w_l = w
    END IF

    ! Modify the Fock Matrix
    kr = 0
    DO iL = 1, natorb(1)
       i = se_orbital_pointer(iL)
       DO kL = 1, natorb(2)
          k = se_orbital_pointer(kL)
          kr = kr + 1
          tmp  = w_l(kr) * factor
          a(1) = tmp * drp(1)
          a(2) = tmp * drp(2)
          a(3) = tmp * drp(3)
          ! Coulomb
          IF (.NOT.switch) THEN
             tmp = pj_tot(k,k) * pi_tot(i,i)
          ELSE
             tmp = pi_tot(k,k) * pj_tot(i,i)
          END IF
          force(1) = force(1) + a(1) * tmp
          force(2) = force(2) + a(2) * tmp
          force(3) = force(3) + a(3) * tmp
       END DO
    END DO

  END SUBROUTINE dfock2C_r3

! *****************************************************************************
!> \brief  Coulomb interaction multipolar correction
!> \author Teodoro Laino [tlaino]
!> \date   05.2009 [tlaino]
! *****************************************************************************
  SUBROUTINE se_coulomb_ij_interaction (atom_a, atom_b, my_task, do_forces, do_efield,&
       do_stress, charges, dipoles, quadrupoles, force_ab, efield0, efield1, efield2, &
       rab2, rab, integral_value, ptens11, ptens12, ptens13, ptens21, ptens22, ptens23,&
       ptens31, ptens32, ptens33, error)
    INTEGER, INTENT(IN)                      :: atom_a, atom_b
    LOGICAL, DIMENSION(3)                    :: my_task
    LOGICAL, INTENT(IN)                      :: do_forces, do_efield, &
                                                do_stress
    REAL(KIND=dp), DIMENSION(:), POINTER     :: charges
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dipoles
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: quadrupoles
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: force_ab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: efield0
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: efield1, efield2
    REAL(KIND=dp), INTENT(IN)                :: rab2
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: integral_value
    REAL(KIND=dp), INTENT(INOUT)             :: ptens11, ptens12, ptens13, &
                                                ptens21, ptens22, ptens23, &
                                                ptens31, ptens32, ptens33
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'se_coulomb_ij_interaction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: a, b, c, d, e, i, j, k
    LOGICAL                                  :: do_efield0, do_efield1, &
                                                do_efield2, failure, &
                                                force_eval
    LOGICAL, DIMENSION(3)                    :: do_task
    LOGICAL, DIMENSION(3, 3)                 :: task
    REAL(KIND=dp) :: ch_i, ch_j, ef0_i, ef0_j, eloc, energy, fac, fac_ij, ir, &
      irab2, r, tij, tmp, tmp1, tmp11, tmp12, tmp13, tmp2, tmp21, tmp22, &
      tmp23, tmp31, tmp32, tmp33, tmp_ij, tmp_ji
    REAL(KIND=dp), DIMENSION(0:5)            :: f
    REAL(KIND=dp), DIMENSION(3)              :: dp_i, dp_j, ef1_i, ef1_j, fr, &
                                                tij_a
    REAL(KIND=dp), DIMENSION(3, 3)           :: ef2_i, ef2_j, qp_i, qp_j, &
                                                tij_ab
    REAL(KIND=dp), DIMENSION(3, 3, 3)        :: tij_abc
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3)     :: tij_abcd
    REAL(KIND=dp), DIMENSION(3, 3, 3, 3, 3)  :: tij_abcde

    failure = .FALSE.
    do_task  = my_task
    energy   = 0.0_dp
    DO i = 1, 3
       IF (do_task(i)) THEN
          SELECT CASE(i)
          CASE(1)
             do_task(1) = (charges(atom_a)/=0.0_dp).OR.(charges(atom_b)/=0.0_dp)
          CASE(2)
             do_task(2) = (ANY(dipoles(:,atom_a)/=0.0_dp)).OR.(ANY(dipoles(:,atom_b)/=0.0_dp))
          CASE(3)
             do_task(3) = (ANY(quadrupoles(:,:,atom_a)/=0.0_dp)).OR.(ANY(quadrupoles(:,:,atom_b)/=0.0_dp))
          END SELECT
       END IF
    END DO
    DO i = 1,3
       DO j =i,3
          task(j,i) = do_task(i).AND.do_task(j)
          task(i,j) = task(j,i)
       END DO
    END DO
    do_efield0 = do_efield.AND.ASSOCIATED(efield0)
    do_efield1 = do_efield.AND.ASSOCIATED(efield1)
    do_efield2 = do_efield.AND.ASSOCIATED(efield2)

    fac_ij = 1.0_dp
    IF (atom_a==atom_b) fac_ij = 0.5_dp
#define PURE_COULOMB
#include "ewalds_multipole_sr.f90"
#undef PURE_COULOMB

    IF (PRESENT(integral_value)) THEN
       integral_value = eloc
    END IF
    IF (do_forces) THEN
       force_ab = fr
    END IF
  END SUBROUTINE se_coulomb_ij_interaction

END MODULE se_fock_matrix_integrals

