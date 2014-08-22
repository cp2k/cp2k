!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Integrals for semi-empiric methods
!> \par History
!>      JGH (27.10.2004) : separate routine for nuclear attraction integrals
!>      JGH (13.03.2006) : tapering function
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!> \author JGH (11.10.2004)
! *****************************************************************************
MODULE semi_empirical_int_num

  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_pchg,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_undef,&
                                             do_multipole_none,&
                                             do_se_IS_kdso_d
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             evolt
  USE semi_empirical_int_arrays,       ONLY: &
       ijkl_ind, ijkl_sym, inddd, inddp, indexa, indexb, indpp, int2c_type, &
       l_index, rij_threshold
  USE semi_empirical_int_utils,        ONLY: ijkl_d,&
                                             ijkl_sp,&
                                             rot_2el_2c_first,&
                                             rotmat,&
                                             store_2el_2c_diag
  USE semi_empirical_types,            ONLY: rotmat_create,&
                                             rotmat_release,&
                                             rotmat_type,&
                                             se_int_control_type,&
                                             se_int_screen_type,&
                                             se_taper_type,&
                                             semi_empirical_type,&
                                             setup_se_int_control_type
  USE taper_types,                     ONLY: taper_eval
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
#include "semi_empirical_int_args.h"

  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_num'
  PUBLIC :: rotint_num,      rotnuc_num,     corecore_num,&
           drotint_num,     drotnuc_num,    dcorecore_num,&
       ssss_nucint_num, core_nucint_num,        terep_num,&
         nucint_sp_num,    terep_sp_num,      terep_d_num,&
          nucint_d_num, corecore_el_num, dcorecore_el_num

CONTAINS

! *****************************************************************************
!> \brief Computes the two particle interactions in the lab frame
!> \param sepi Atomic parameters of first atom
!> \param sepj Atomic parameters of second atom
!> \param rijv Coordinate vector i -> j
!> \param w Array of two-electron repulsion integrals.
!> \param se_int_control input parameters that control the calculation of SE
!>                           integrals (shortrange, R3 residual, screening type)
!> \param se_taper ...
!> \param error ...
!> \note routine adapted from mopac7 (rotate)
!>       written by Ernest R. Davidson, Indiana University.
!>       Teodoro Laino [tlaino] - University of Zurich 04.2008 : major rewriting
! *****************************************************************************
  SUBROUTINE rotint_num (sepi,sepj,rijv,w,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
    REAL(dp), DIMENSION(2025), INTENT(OUT)   :: w
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotint_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i1, ii, ij, ij1, iminus, &
                                                istep, iw_loc, j, j1, jj, k, &
                                                kk, kl, l, limij, limkl, mm
    LOGICAL, DIMENSION(45, 45)               :: logv
    REAL(dp)                                 :: rij
    REAL(KIND=dp)                            :: cc, wrepp
    REAL(KIND=dp), DIMENSION(2025)           :: ww
    REAL(KIND=dp), DIMENSION(45, 45)         :: v
    REAL(KIND=dp), DIMENSION(491)            :: rep
    TYPE(rotmat_type), POINTER               :: ij_matrix

    NULLIFY(ij_matrix)
    rij  =  DOT_PRODUCT(rijv,rijv)
    IF (rij > rij_threshold) THEN
       ! The repulsion integrals over molecular frame (w) are stored in the
       ! order in which they will later be used.  ie.  (i,j/k,l) where
       ! j.le.i  and  l.le.k     and l varies most rapidly and i least
       ! rapidly.  (anti-normal computer storage)
       rij = SQRT(rij)

       ! Create the rotation matrix
       CALL rotmat_create(ij_matrix, error)
       CALL rotmat (sepi, sepj, rijv, rij, ij_matrix, do_derivatives=.FALSE., error=error)

       ! Compute Integrals in Diatomic Frame
       CALL terep_num(sepi,sepj,rij,rep,se_taper=se_taper,se_int_control=se_int_control,error=error)

       ! Rotate Integrals
       ii = sepi%natorb
       kk = sepj%natorb
       IF (ii*kk > 0) THEN
          limij = sepi%atm_int_size
          limkl = sepj%atm_int_size
          istep = limkl * limij
          DO i1 = 1, istep
             ww(i1) = 0.0_dp
          END DO

          ! First step in rotation of integrals
          CALL rot_2el_2c_first(sepi, sepj, rijv, se_int_control, se_taper, .FALSE., ii, kk, rep, &
               logv, ij_matrix, v, lgrad=.FALSE.,error=error)

          ! Second step in rotation of integrals
          DO i1 = 1, ii
             DO j1 = 1, i1
                ij = indexa(i1, j1)
                jj = indexb(i1, j1)
                mm = int2c_type(jj)
                DO k = 1, kk
                   DO l = 1, k
                      kl = indexb(k, l)
                      IF (logv(ij, kl)) THEN
                         wrepp = v(ij, kl)
                         SELECT CASE (mm)
                         CASE(1)
                            ! (SS/)
                            i = 1
                            j = 1
                            iw_loc = (indexb(i,j)-1) * limkl + kl
                            ww(iw_loc) = wrepp
                         CASE(2)
                            ! (SP/)
                            j = 1
                            DO i = 1, 3
                               iw_loc = (indexb(i+1,j)-1) * limkl + kl
                               ww(iw_loc) = ww(iw_loc) + ij_matrix%sp(i1-1, i) * wrepp
                            END DO
                         CASE(3)
                            ! (PP/)
                            DO i = 1, 3
                               cc = ij_matrix%pp(i, i1-1, j1-1)
                               iw_loc = (indexb(i+1,i+1)-1) * limkl + kl
                               ww(iw_loc) = ww(iw_loc) + cc * wrepp
                               iminus = i - 1
                               IF (iminus /= 0) THEN
                                  DO j = 1, iminus
                                     cc = ij_matrix%pp(1+i+j, i1-1, j1-1)
                                     iw_loc = (indexb(i+1,j+1)-1) * limkl + kl
                                     ww(iw_loc) = ww(iw_loc) + cc * wrepp
                                  END DO
                               END IF
                            END DO
                         CASE(4)
                            ! (SD/)
                            j = 1
                            DO i = 1, 5
                               iw_loc = (indexb(i+4,j)-1) * limkl + kl
                               ww(iw_loc) = ww(iw_loc) + ij_matrix%sd(i1-4, i) * wrepp
                            END DO
                         CASE(5)
                            ! (DP/)
                            DO i = 1, 5
                               DO j = 1, 3
                                  iw_loc = (indexb(i+4,j+1)-1) * limkl + kl
                                  ij1 = 3 * (i-1) + j
                                  ww(iw_loc) = ww(iw_loc) + ij_matrix%pd(ij1, i1-4, j1-1) * wrepp
                               END DO
                            END DO
                         CASE(6)
                            ! (DD/)
                            DO i = 1, 5
                               cc = ij_matrix%dd(i, i1-4, j1-4)
                               iw_loc = (indexb(i+4,i+4)-1) * limkl + kl
                               ww(iw_loc) = ww(iw_loc) + cc * wrepp
                               iminus = i - 1
                               IF (iminus /= 0) THEN
                                  DO j = 1, iminus
                                     ij1 = inddd(i, j)
                                     cc = ij_matrix%dd(ij1, i1-4, j1-4)
                                     iw_loc = (indexb(i+4,j+4)-1) * limkl + kl
                                     ww(iw_loc) = ww(iw_loc) + cc * wrepp
                                  END DO
                               END IF
                            END DO
                         END SELECT
                      END IF
                   END DO
                END DO
             END DO
          END DO
       END IF
       CALL rotmat_release(ij_matrix, error)

       ! Store two electron integrals in the triangular format
       CALL store_2el_2c_diag(limij, limkl, ww, w, error=error)
    ENDIF
  END SUBROUTINE rotint_num

! *****************************************************************************
!> \brief Calculates the derivative pf two-electron repulsion integrals and the
!>      nuclear attraction integrals w.r.t. |r|
!> \param sepi paramters of atom i
!> \param sepj paramters of atom j
!> \param rij interatomic distance
!> \param rep array of two-electron repulsion integrals
!> \param se_taper ...
!> \param se_int_control input parameters that control the calculation of SE
!>                         integrals (shortrange, R3 residual, screening type)
!> \param error ...
!> \par History
!>      03.2008 created [tlaino]
!> \author Teodoro Laino [tlaino] - Zurich University
! *****************************************************************************
  SUBROUTINE terep_num(sepi, sepj, rij, rep, se_taper, se_int_control, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(491), INTENT(OUT)    :: rep
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'terep_num', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ft
    TYPE(se_int_screen_type)                 :: se_int_screen

    failure           = .FALSE.
    ft                = taper_eval(se_taper%taper, rij, error)
    ! In case of dumped integrals compute an additional taper term
    IF (se_int_control%integral_screening==do_se_IS_kdso_d) THEN
       se_int_screen%ft  = taper_eval(se_taper%taper_add, rij, error)
    END IF

    ! Contribution from sp shells
    CALL terep_sp_num(sepi, sepj, rij, rep, se_int_control, se_int_screen, ft,&
         error)

    IF (sepi%dorb.OR.sepj%dorb) THEN
       ! Compute the contribution from d shells
       CALL terep_d_num (sepi, sepj, rij, rep, se_int_control, se_int_screen,&
            ft, error)
    END IF

  END SUBROUTINE terep_num

! *****************************************************************************
!> \brief Calculates the  two-electron repulsion integrals - sp shell only
!> \param sepi paramters of atom i
!> \param sepj paramters of atom j
!> \param rij ...
!> \param rep array of two-electron repulsion integrals
!> \param se_int_control input parameters that control the calculation of SE
!>                         integrals (shortrange, R3 residual, screening type)
!> \param se_int_screen contains information for computing the screened
!>                         integrals KDSO-D
!> \param ft ...
!> \param error ...
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
! *****************************************************************************
  SUBROUTINE terep_sp_num (sepi, sepj, rij, rep, se_int_control, se_int_screen,&
       ft, error )
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(491), INTENT(OUT)    :: rep
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    REAL(dp), INTENT(IN)                     :: ft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'terep_sp_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij, j, k, kl, l, lasti, &
                                                lastj, li, lj, lk, ll, nold, &
                                                numb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp

    failure = .FALSE.
    lasti = sepi%natorb
    lastj = sepj%natorb
    DO i = 1, MIN(lasti,4)
       li = l_index(i)
       DO j = 1, i
          lj = l_index(j)
          ij = indexa(i, j)
          DO k = 1, MIN(lastj,4)
             lk = l_index(k)
             DO l = 1, k
                ll = l_index(l)
                kl = indexa(k, l)

                numb = ijkl_ind(ij, kl)
                IF (numb > 0) THEN
                   nold = ijkl_sym(numb)
                   IF (nold > 0) THEN
                      rep(numb) =  rep( nold)
                   ELSE IF (nold < 0) THEN
                      rep(numb) = -rep(-nold)
                   ELSE IF (nold == 0) THEN
                      tmp = ijkl_sp (sepi, sepj, ij, kl, li, lj, lk, ll, 0, rij,&
                                     se_int_control, se_int_screen, do_method_undef,&
                                     error) * ft
                      rep(numb) = tmp
                   END IF
                END IF
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE terep_sp_num

! *****************************************************************************
!> \brief Calculates the  two-electron repulsion integrals - d shell only
!> \param sepi ...
!> \param sepj ...
!> \param rij interatomic distance
!> \param rep array of two-electron repulsion integrals
!> \param se_int_control input parameters that control the calculation of SE
!>                         integrals (shortrange, R3 residual, screening type)
!> \param se_int_screen contains information for computing the screened
!>                         integrals KDSO-D
!> \param ft ...
!> \param error ...
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
! *****************************************************************************
  SUBROUTINE terep_d_num (sepi, sepj, rij, rep, se_int_control, se_int_screen,&
       ft, error )
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(491), INTENT(INOUT)  :: rep
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    REAL(dp), INTENT(IN)                     :: ft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'terep_d_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij, j, k, kl, l, lasti, &
                                                lastj, li, lj, lk, ll, nold, &
                                                numb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp

    failure = .FALSE.
    lasti = sepi%natorb
    lastj = sepj%natorb
    DO i = 1, lasti
       li = l_index(i)
       DO j = 1, i
          lj = l_index(j)
          ij = indexa(i, j)
          DO k = 1, lastj
             lk = l_index(k)
             DO l = 1, k
                ll = l_index(l)
                kl = indexa(k, l)

                numb = ijkl_ind(ij, kl)
                ! From 1 to 34 we store integrals involving sp shells
                IF (numb > 34) THEN
                   nold = ijkl_sym(numb)
                   IF (nold > 34) THEN
                      rep(numb) =  rep( nold)
                   ELSE IF (nold <-34) THEN
                      rep(numb) = -rep(-nold)
                   ELSE IF (nold == 0) THEN
                      tmp = ijkl_d (sepi, sepj, ij, kl, li, lj, lk, ll, 0, rij,&
                                    se_int_control, se_int_screen, do_method_undef,&
                                    error) * ft
                      rep(numb) = tmp
                   END IF
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE terep_d_num

! *****************************************************************************
!> \brief Computes the two-particle interactions.
!> \param sepi Atomic parameters of first atom
!> \param sepj Atomic parameters of second atom
!> \param rijv Coordinate vector i -> j
!> \param e1b Array of electron-nuclear attraction integrals,
!>                       e1b = Electron on atom ni attracting nucleus of nj.
!> \param e2a Array of electron-nuclear attraction integrals,
!>                       e2a = Electron on atom nj attracting nucleus of ni.
!> \param itype ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
!> \note routine adapted from mopac7 (rotate)
!>       written by Ernest R. Davidson, Indiana University.
!>       Teodoro Laino [tlaino] - University of Zurich 04.2008 : major rewriting
!>       Teodoro Laino [tlaino] - University of Zurich 04.2008 : removed the core-core part
! *****************************************************************************
  SUBROUTINE rotnuc_num(sepi,sepj,rijv,e1b,e2a,itype,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
    REAL(dp), DIMENSION(45), INTENT(OUT), &
      OPTIONAL                               :: e1b, e2a
    INTEGER, INTENT(IN)                      :: itype
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotnuc_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idd, idp, ind1, ind2, ipp, &
                                                j, last_orbital(2), m, n
    LOGICAL                                  :: l_e1b, l_e2a, task(2)
    REAL(dp)                                 :: rij
    REAL(dp), DIMENSION(10, 2)               :: core
    REAL(dp), DIMENSION(45)                  :: tmp
    TYPE(rotmat_type), POINTER               :: ij_matrix

    NULLIFY(ij_matrix)
    l_e1b  = PRESENT(e1b)
    l_e2a  = PRESENT(e2a)
    rij    =  DOT_PRODUCT(rijv,rijv)

    IF (rij > rij_threshold) THEN
       rij = SQRT(rij)
       ! Create the rotation matrix
       CALL rotmat_create(ij_matrix, error)
       CALL rotmat (sepi, sepj, rijv, rij, ij_matrix, do_derivatives=.FALSE., error=error)

       ! Compute Integrals in Diatomic Frame
       CALL core_nucint_num(sepi,sepj,rij,core=core,itype=itype,se_taper=se_taper,&
            se_int_control=se_int_control,error=error)

       ! Copy parameters over to arrays for do loop.
       last_orbital(1)  = sepi%natorb
       last_orbital(2)  = sepj%natorb
       task (1)         = l_e1b
       task (2)         = l_e2a
       DO n = 1, 2
          IF (.NOT.task(n)) CYCLE
          DO i = 1, last_orbital(n)
             ind1 = i - 1
             DO j = 1, i
                ind2 = j - 1
                m = (i*(i-1))/2 + j
                ! Perform Rotations ...
                IF (ind2 == 0) THEN
                   IF (ind1 == 0) THEN
                      ! Type of Integral (SS/)
                      tmp(m) = core(1, n)
                   ELSE IF (ind1 < 4) THEN
                      ! Type of Integral (SP/)
                      tmp(m) = ij_matrix%sp(1,   ind1) * core(2, n)
                   ELSE
                      ! Type of Integral (SD/)
                      tmp(m) = ij_matrix%sd(1, ind1-3) * core(5, n)
                   END IF
                ELSE IF (ind2 < 4) THEN
                   IF (ind1 < 4) THEN
                      ! Type of Integral (PP/)
                      ipp = indpp(ind1, ind2)
                      tmp(m) = core(3, n) *  ij_matrix%pp(ipp, 1, 1) + &
                               core(4, n) * (ij_matrix%pp(ipp, 2, 2)+ij_matrix%pp(ipp, 3, 3))
                   ELSE
                      ! Type of Integral (PD/)
                      idp = inddp(ind1-3, ind2)
                      tmp(m) = core(6, n) *  ij_matrix%pd(idp, 1, 1) + &
                               core(8, n) * (ij_matrix%pd(idp, 2, 2)+ij_matrix%pd(idp, 3, 3))
                   END IF
                ELSE
                   ! Type of Integral (DD/)
                   idd = inddd(ind1-3, ind2-3)
                   tmp(m) = core( 7, n) *  ij_matrix%dd(idd, 1, 1) + &
                            core( 9, n) * (ij_matrix%dd(idd, 2, 2)+ij_matrix%dd(idd, 3, 3)) + &
                            core(10, n) * (ij_matrix%dd(idd, 4, 4)+ij_matrix%dd(idd, 5, 5))
                END IF
             END DO
          END DO
          IF (n==1) THEN
             DO i = 1, sepi%atm_int_size
                e1b(i) = -tmp(i)
             END DO
          END IF
          IF (n==2) THEN
             DO i = 1, sepj%atm_int_size
                e2a(i) = -tmp(i)
             END DO
          END IF
       END DO
       CALL rotmat_release(ij_matrix, error)
    END IF
  END SUBROUTINE rotnuc_num

! *****************************************************************************
!> \brief Computes the core-core interactions.
!> \param sepi Atomic parameters of first atom
!> \param sepj Atomic parameters of second atom
!> \param rijv Coordinate vector i -> j
!> \param enuc nuclear-nuclear repulsion term.
!> \param itype ...
!> \param se_int_control input parameters that control the calculation of SE
!>                           integrals (shortrange, R3 residual, screening type)
!> \param se_taper ...
!> \param error ...
!> \note routine adapted from mopac7 (rotate)
!>       written by Ernest R. Davidson, Indiana University.
!>       Teodoro Laino [tlaino] - University of Zurich 04.2008 : major rewriting
!>       Teodoro Laino [tlaino] - University of Zurich 04.2008 : splitted from rotnuc
! *****************************************************************************
  SUBROUTINE corecore_num(sepi,sepj,rijv,enuc,itype,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
    REAL(dp), INTENT(OUT)                    :: enuc
    INTEGER, INTENT(IN)                      :: itype
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'corecore_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ig, nt
    REAL(dp) :: aab, alpi, alpj, apdg, ax, dai, daj, dbi, dbj, pai, paj, pbi, &
      pbj, qcorr, rij, rija, scale, ssss, ssss_sr, tmp, xab, zaf, zbf, zz
    REAL(dp), DIMENSION(4)                   :: fni1, fni2, fni3, fnj1, fnj2, &
                                                fnj3
    TYPE(se_int_control_type)                :: se_int_control_off

    rij    =  DOT_PRODUCT(rijv,rijv)
    enuc   =  0.0_dp
    IF (rij > rij_threshold) THEN
       rij = SQRT(rij)
       CALL setup_se_int_control_type(se_int_control_off, shortrange=.FALSE., do_ewald_r3=.FALSE.,&
            do_ewald_gks=.FALSE., integral_screening=se_int_control%integral_screening,&
            max_multipole=do_multipole_none, pc_coulomb_int=.FALSE.)
       CALL ssss_nucint_num(sepi,sepj,rij,ssss=ssss,itype=itype,se_taper=se_taper,&
            se_int_control=se_int_control_off,error=error)
       ! In case let's compute the short-range part of the (ss|ss) integral
       IF (se_int_control%shortrange) THEN
          CALL ssss_nucint_num(sepi,sepj,rij,ssss=ssss_sr,itype=itype,se_taper=se_taper,&
               se_int_control=se_int_control,error=error)
       ELSE
          ssss_sr = ssss
       END IF
       zz = sepi%zeff*sepj%zeff
       ! Nuclear-Nuclear electrostatic contribution
       enuc = zz*ssss_sr
       ! Method dependent correction terms
       IF(itype /= do_method_pm6) THEN
          alpi = sepi%alp
          alpj = sepj%alp
          scale = EXP(-alpi*rij)+EXP(-alpj*rij)

          nt=sepi%z+sepj%z
          IF(nt == 8 .OR. nt == 9) THEN
             IF(sepi%z == 7.OR.sepi%z == 8) scale=scale+(angstrom*rij-1._dp)*EXP(-alpi*rij)
             IF(sepj%z == 7.OR.sepj%z == 8) scale=scale+(angstrom*rij-1._dp)*EXP(-alpj*rij)
          ENDIF
          scale= ABS(scale*zz*ssss)
          zz   = zz/rij
          IF(itype == do_method_am1 .OR. itype == do_method_pm3 .OR. itype == do_method_pdg) THEN
             IF(itype == do_method_am1 .AND. sepi%z == 5) THEN
                !special case AM1 Boron
                SELECT CASE (sepj%z)
                CASE DEFAULT
                   nt=1
                CASE (1)
                   nt=2
                CASE (6)
                   nt=3
                CASE (9,17,35,53)
                   nt=4
                END SELECT
                fni1(:) = sepi%bfn1(:,nt)
                fni2(:) = sepi%bfn2(:,nt)
                fni3(:) = sepi%bfn3(:,nt)
             ELSE
                fni1(:) = sepi%fn1(:)
                fni2(:) = sepi%fn2(:)
                fni3(:) = sepi%fn3(:)
             END IF
             IF(itype == do_method_am1 .AND. sepj%z == 5) THEN
                !special case AM1 Boron
                SELECT CASE (sepi%z)
                CASE DEFAULT
                   nt=1
                CASE (1)
                   nt=2
                CASE (6)
                   nt=3
                CASE (9,17,35,53)
                   nt=4
                END SELECT
                fnj1(:) = sepj%bfn1(:,nt)
                fnj2(:) = sepj%bfn2(:,nt)
                fnj3(:) = sepj%bfn3(:,nt)
             ELSE
                fnj1(:) = sepj%fn1(:)
                fnj2(:) = sepj%fn2(:)
                fnj3(:) = sepj%fn3(:)
             END IF
             ! AM1/PM3/PDG correction to nuclear repulsion
             DO ig=1,SIZE(fni1)
                IF(ABS(fni1(ig)) > 0._dp) THEN
                   ax = fni2(ig)*(rij-fni3(ig))**2
                   IF(ax <= 25._dp) THEN
                      scale=scale +zz*fni1(ig)*EXP(-ax)
                   ENDIF
                ENDIF
                IF(ABS(fnj1(ig)) > 0._dp) THEN
                   ax = fnj2(ig)*(rij-fnj3(ig))**2
                   IF(ax <= 25._dp) THEN
                      scale=scale +zz*fnj1(ig)*EXP(-ax)
                   ENDIF
                ENDIF
             END DO
          ENDIF
          IF(itype == do_method_pdg) THEN
             ! PDDG function
             zaf=sepi%zeff/nt
             zbf=sepj%zeff/nt
             pai=sepi%pre(1)
             pbi=sepi%pre(2)
             paj=sepj%pre(1)
             pbj=sepj%pre(2)
             dai=sepi%d(1)
             dbi=sepi%d(2)
             daj=sepj%d(1)
             dbj=sepj%d(2)
             apdg=10._dp*angstrom**2
             qcorr = &
                  (zaf*pai+zbf*paj)*EXP(-apdg*(rij-dai-daj)**2)+ &
                  (zaf*pai+zbf*pbj)*EXP(-apdg*(rij-dai-dbj)**2)+ &
                  (zaf*pbi+zbf*paj)*EXP(-apdg*(rij-dbi-daj)**2)+ &
                  (zaf*pbi+zbf*pbj)*EXP(-apdg*(rij-dbi-dbj)**2)
          ELSEIF (itype==do_method_pchg) THEN
             qcorr = 0.0_dp
             scale = 0.0_dp
          ELSE
             qcorr = 0.0_dp
          END IF
       ELSE
          rija  = rij*angstrom
          scale = zz*ssss
          ! PM6 core-core terms
          xab = sepi%xab(sepj%z)
          aab = sepi%aab(sepj%z)
          IF ( (sepi%z == 1 .AND. (sepj%z == 6 .OR. sepj%z == 7 .OR. sepj%z == 8)) .OR. &
               (sepj%z == 1 .AND. (sepi%z == 6 .OR. sepi%z == 7 .OR. sepi%z == 8)) ) THEN
             ! Special Case O-H or N-H or C-H
             scale = scale*(2._dp*xab*EXP(-aab*rija*rija))
          ELSEIF ( sepi%z == 6 .AND. sepj%z == 6 ) THEN
             ! Special Case C-C
             scale = scale*(2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))+9.28_dp*EXP(-5.98_dp*rija))
          ELSEIF ( (sepi%z == 8 .AND. sepj%z == 14) .OR. &
               (sepj%z == 8 .AND. sepi%z == 14) ) THEN
             ! Special Case Si-O
             scale = scale*(2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))-0.0007_dp*EXP(-(rija-2.9_dp)**2))
          ELSE
             ! General Case
             ! Factor of 2 found by experiment
             scale = scale*(2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6)))
          END IF
          ! General correction term a*exp(-b*(rij-c)^2)
          qcorr = (sepi%a*EXP(-sepi%b*(rij-sepi%c)**2))*sepi%zeff*sepj%zeff/rij +&
                  (sepj%a*EXP(-sepj%b*(rij-sepj%c)**2))*sepi%zeff*sepj%zeff/rij
          ! Hard core repulsion
          tmp   = (REAL(sepi%z,dp)**(1._dp/3._dp) + REAL(sepj%z,dp)**(1._dp/3._dp))
          qcorr = qcorr + 1.e-8_dp/evolt * ( tmp/rija )**12
       END IF
       enuc=enuc+scale+qcorr
    END IF
  END SUBROUTINE corecore_num

! *****************************************************************************
!> \brief Computes the electrostatic core-core interactions only.
!> \param sepi Atomic parameters of first atom
!> \param sepj Atomic parameters of second atom
!> \param rijv Coordinate vector i -> j
!> \param enuc nuclear-nuclear repulsion term.
!> \param itype ...
!> \param se_int_control input parameters that control the calculation of SE
!>                           integrals (shortrange, R3 residual, screening type)
!> \param se_taper ...
!> \param error ...
!> \author Teodoro Laino [tlaino] - 05.2009
! *****************************************************************************
  SUBROUTINE corecore_el_num(sepi,sepj,rijv,enuc,itype,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
    REAL(dp), INTENT(OUT)                    :: enuc
    INTEGER, INTENT(IN)                      :: itype
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'corecore_el_num', &
      routineP = moduleN//':'//routineN

    REAL(dp)                                 :: rij, ssss, ssss_sr, zz
    TYPE(se_int_control_type)                :: se_int_control_off

    rij    =  DOT_PRODUCT(rijv,rijv)
    enuc   =  0.0_dp
    IF (rij > rij_threshold) THEN
       rij = SQRT(rij)
       CALL setup_se_int_control_type(se_int_control_off, shortrange=.FALSE., do_ewald_r3=.FALSE.,&
            do_ewald_gks=.FALSE., integral_screening=se_int_control%integral_screening,&
            max_multipole=do_multipole_none, pc_coulomb_int=.FALSE.)
       CALL ssss_nucint_num(sepi,sepj,rij,ssss=ssss,itype=itype,se_taper=se_taper,&
            se_int_control=se_int_control_off,error=error)
       ! In case let's compute the short-range part of the (ss|ss) integral
       IF (se_int_control%shortrange.OR.se_int_control%pc_coulomb_int) THEN
          CALL ssss_nucint_num(sepi,sepj,rij,ssss=ssss_sr,itype=itype,se_taper=se_taper,&
               se_int_control=se_int_control,error=error)
       ELSE
          ssss_sr = ssss
       END IF
       zz = sepi%zeff*sepj%zeff
       ! Nuclear-Nuclear electrostatic contribution
       enuc = zz*ssss_sr
    END IF
  END SUBROUTINE corecore_el_num

! *****************************************************************************
!> \brief Calculates the SSSS integrals (main driver)
!> \param sepi paramters of atom i
!> \param sepj paramters of atom j
!> \param rij interatomic distance
!> \param ssss derivative of (ssss) integral
!>                          derivatives are intended w.r.t. rij
!> \param itype type of semi_empirical model
!>                          extension to the original routine to compute qm/mm integrals
!> \param se_taper ...
!> \param se_int_control input parameters that control the calculation of SE
!>                          integrals (shortrange, R3 residual, screening type)
!> \param error ...
!> \par History
!>      03.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE ssss_nucint_num(sepi, sepj, rij, ssss, itype, se_taper, se_int_control, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(OUT)                    :: ssss
    INTEGER, INTENT(IN)                      :: itype
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ssss_nucint_num', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ft
    TYPE(se_int_screen_type)                 :: se_int_screen

    failure   = .FALSE.

    ! Computing Tapering function
    ft               = 1.0_dp
    IF (itype /=do_method_pchg) THEN
       ft = taper_eval(se_taper%taper, rij, error)
    END IF

    ! In case of dumped integrals compute an additional taper term
    IF (se_int_control%integral_screening==do_se_IS_kdso_d) THEN
       se_int_screen%ft = 1.0_dp
       IF (itype /=do_method_pchg) THEN
          se_int_screen%ft  = taper_eval(se_taper%taper_add, rij, error)
       END IF
    END IF

    ! Contribution from the sp shells
    CALL nucint_sp_num(sepi, sepj, rij, ssss=ssss, itype=itype,&
         se_int_control=se_int_control,se_int_screen=se_int_screen,&
         error=error)

    ! Tapering the integrals
    ssss = ft*ssss

  END SUBROUTINE ssss_nucint_num

! *****************************************************************************
!> \brief Calculates the nuclear attraction integrals CORE (main driver)
!> \param sepi paramters of atom i
!> \param sepj paramters of atom j
!> \param rij interatomic distance
!> \param core derivative of 4 X 2 array of electron-core attraction integrals
!>                          derivatives are intended w.r.t. rij
!> \param itype type of semi_empirical model
!>                          extension to the original routine to compute qm/mm integrals
!> \param se_taper ...
!> \param se_int_control input parameters that control the calculation of SE
!>                          integrals (shortrange, R3 residual, screening type)
!> \param error ...
!> \par History
!>      03.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE core_nucint_num(sepi, sepj, rij, core, itype, se_taper, se_int_control, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(10, 2), INTENT(OUT)  :: core
    INTEGER, INTENT(IN)                      :: itype
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'core_nucint_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ft
    TYPE(se_int_screen_type)                 :: se_int_screen

    failure   = .FALSE.
    ! Computing the Tapering function
    ft               = 1.0_dp
    IF (itype /=do_method_pchg) THEN
       ft = taper_eval(se_taper%taper, rij, error)
    END IF

    ! In case of dumped integrals compute an additional taper term
    IF (se_int_control%integral_screening==do_se_IS_kdso_d) THEN
       se_int_screen%ft = 1.0_dp
       IF (itype /=do_method_pchg) THEN
          se_int_screen%ft  = taper_eval(se_taper%taper_add, rij, error)
       END IF
    END IF

    ! Contribution from the sp shells
    CALL nucint_sp_num(sepi, sepj, rij, core=core, itype=itype,&
         se_int_control=se_int_control,se_int_screen=se_int_screen, error=error)

    IF (sepi%dorb.OR.sepj%dorb) THEN
       ! Compute the contribution from d shells
       CALL nucint_d_num(sepi, sepj, rij, core, itype, se_int_control,&
            se_int_screen, error=error)
    END IF

    ! Tapering the integrals
    DO i = 1, sepi%core_size
       core(i,1)  = ft*core(i,1)
    END DO
    DO i = 1, sepj%core_size
       core(i,2)  = ft*core(i,2)
    END DO

  END SUBROUTINE core_nucint_num

! *****************************************************************************
!> \brief ...
!> \param sepi ...
!> \param sepj ...
!> \param rij ...
!> \param ssss ...
!> \param core ...
!> \param itype ...
!> \param se_int_control ...
!> \param se_int_screen ...
!> \param error ...
!> \par History
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : new driver
!>                    for computing integrals
!>      Teodoro Laino (04.2008) [tlaino] - Totally rewritten: nothing to do with
!>                    the old version
! *****************************************************************************

  SUBROUTINE nucint_sp_num( sepi, sepj, rij, ssss, core, itype, se_int_control,&
       se_int_screen, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(INOUT), OPTIONAL        :: ssss
    REAL(dp), DIMENSION(10, 2), &
      INTENT(INOUT), OPTIONAL                :: core
    INTEGER, INTENT(IN)                      :: itype
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nucint_sp_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ij, kl
    LOGICAL                                  :: failure, l_core, l_ssss, si, &
                                                sj

    failure = .FALSE.
    l_core  = PRESENT(core)
    l_ssss  = PRESENT(ssss)
    IF (.NOT.(l_core.OR.l_ssss)) RETURN
    si  = (sepi%natorb > 1)
    sj  = (sepj%natorb > 1)

    ij = indexa(1, 1)
    IF (l_ssss) THEN
       ! Store the value for <S  S  | S  S  > (Used for computing the core-core interactions)
       ssss = ijkl_sp (sepi, sepj, ij, ij, 0, 0, 0, 0, -1, rij, CPPint_args)
    END IF

    IF (l_core) THEN
       !     <S  S  | S  S  >
       kl = indexa(1, 1)
       core(1,1)    = ijkl_sp (sepi, sepj, kl, ij, 0, 0, 0, 0, 2, rij, CPPint_args) * sepj%zeff
       IF (si) THEN
          !  <S  P  | S  S  >
          kl = indexa(2, 1)
          core(2,1) = ijkl_sp (sepi, sepj, kl, ij, 0, 1, 0, 0, 2, rij, CPPint_args) * sepj%zeff
          !  <P  P  | S  S  >
          kl = indexa(2, 2)
          core(3,1) = ijkl_sp (sepi, sepj, kl, ij, 1, 1, 0, 0, 2, rij, CPPint_args) * sepj%zeff
          !  <P+ P+ | S  S  >
          kl = indexa(3, 3)
          core(4,1) = ijkl_sp (sepi, sepj, kl, ij, 1, 1, 0, 0, 2, rij, CPPint_args) * sepj%zeff
       END IF

       !     <S  S  | S  S  >
       kl = indexa(1, 1)
       core(1,2)    = ijkl_sp (sepi, sepj, ij, kl, 0, 0, 0, 0, 1, rij, CPPint_args) * sepi%zeff
       IF (sj) THEN
          !  <S  S  | S  P  >
          kl = indexa(2, 1)
          core(2,2) = ijkl_sp (sepi, sepj, ij, kl, 0, 0, 0, 1, 1, rij, CPPint_args) * sepi%zeff
          !  <S  S  | P  P  >
          kl = indexa(2, 2)
          core(3,2) = ijkl_sp (sepi, sepj, ij, kl, 0, 0, 1, 1, 1, rij, CPPint_args) * sepi%zeff
          !  <S  S  | P+ P+ >
          kl = indexa(3, 3)
          core(4,2) = ijkl_sp (sepi, sepj, ij, kl, 0, 0, 1, 1, 1, rij, CPPint_args) * sepi%zeff
       END IF
    END IF

  END SUBROUTINE nucint_sp_num

! *****************************************************************************
!> \brief Calculates the nuclear attraction integrals involving d orbitals
!> \param sepi paramters of atom i
!> \param sepj paramters of atom j
!> \param rij interatomic distance
!> \param core 4 X 2 array of electron-core attraction integrals
!>
!>         The storage of the nuclear attraction integrals  core(kl/ij) iS
!>         (SS/)=1,   (SP/)=2,   (PP/)=3,   (P+P+/)=4,   (SD/)=5,
!>         (DP/)=6,   (DD/)=7,   (D+P+)=8,  (D+D+/)=9,   (D#D#)=10
!>
!>         where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
!> \param itype type of semi_empirical model
!>                         extension to the original routine to compute qm/mm integrals
!> \param se_int_control input parameters that control the calculation of SE
!>                         integrals (shortrange, R3 residual, screening type)
!> \param se_int_screen contains information for computing the screened
!>                         integrals KDSO-D
!> \param error ...
!> \author
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich
! *****************************************************************************
  SUBROUTINE nucint_d_num( sepi, sepj, rij, core, itype, se_int_control, &
       se_int_screen, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(10, 2), &
      INTENT(INOUT)                          :: core
    INTEGER, INTENT(IN)                      :: itype
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_int_screen_type), INTENT(IN)     :: se_int_screen
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'nucint_d_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ij, kl
    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Check if d-orbitals are present
    IF (sepi%dorb.OR.sepj%dorb) THEN
       ij = indexa(1, 1)
       IF (sepj%dorb) THEN
          !  <S S | D S>
          kl = indexa(5, 1)
          core(5, 2) = ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 0, 1, rij, CPPint_args) * sepi%zeff
          !  <S S | D P >
          kl = indexa(5, 2)
          core(6, 2) = ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 1, 1, rij, CPPint_args) * sepi%zeff
          !  <S S | D D >
          kl = indexa(5, 5)
          core(7, 2) = ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, rij, CPPint_args) * sepi%zeff
          !  <S S | D+P+>
          kl = indexa(6, 3)
          core(8, 2) = ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 1, 1, rij, CPPint_args) * sepi%zeff
          !  <S S | D+D+>
          kl = indexa(6, 6)
          core(9, 2) = ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, rij, CPPint_args) * sepi%zeff
          !  <S S | D#D#>
          kl = indexa(8, 8)
          core(10, 2)= ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, rij, CPPint_args) * sepi%zeff
       END IF
       IF (sepi%dorb) THEN
          !  <D S | S S>
          kl = indexa(5, 1)
          core(5, 1) = ijkl_d (sepi, sepj, kl, ij, 2, 0, 0, 0, 2, rij, CPPint_args) * sepj%zeff
          !  <D P | S S >
          kl = indexa(5, 2)
          core(6, 1) = ijkl_d (sepi, sepj, kl, ij, 2, 1, 0, 0, 2, rij, CPPint_args) * sepj%zeff
          !  <D D | S S >
          kl = indexa(5, 5)
          core(7, 1) = ijkl_d (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, rij, CPPint_args) * sepj%zeff
          !  <D+P+| S S >
          kl = indexa(6, 3)
          core(8, 1) = ijkl_d (sepi, sepj, kl, ij, 2, 1, 0, 0, 2, rij, CPPint_args) * sepj%zeff
          !  <D+D+| S S >
          kl = indexa(6, 6)
          core(9, 1) = ijkl_d (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, rij, CPPint_args) * sepj%zeff
          !  <D#D#| S S >
          kl = indexa(8, 8)
          core(10, 1)= ijkl_d (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, rij, CPPint_args) * sepj%zeff
       END IF

    END IF
  END SUBROUTINE nucint_d_num

! *****************************************************************************
!> \brief Numerical Derivatives for rotint
!> \param sepi ...
!> \param sepj ...
!> \param r ...
!> \param dw ...
!> \param delta ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE drotint_num(sepi,sepj,r,dw,delta,se_int_control,se_taper, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: r
    REAL(dp), DIMENSION(3, 2025), &
      INTENT(OUT)                            :: dw
    REAL(dp), INTENT(IN)                     :: delta
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotint_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, nsize
    REAL(dp)                                 :: od
    REAL(dp), DIMENSION(2025)                :: wm, wp
    REAL(dp), DIMENSION(3)                   :: rr

    od    = 0.5_dp/delta
    nsize = sepi%atm_int_size*sepj%atm_int_size
    DO i=1,3
       rr=r
       rr(i)=rr(i)+delta
       CALL rotint_num(sepi,sepj,rr,wp,se_int_control,se_taper=se_taper,error=error)
       rr(i)=rr(i)-2._dp*delta
       CALL rotint_num(sepi,sepj,rr,wm,se_int_control,se_taper=se_taper,error=error)
       DO j = 1, nsize
          dw(i,j) = od * (wp(j) - wm(j))
       END DO
    END DO

  END SUBROUTINE drotint_num

! *****************************************************************************
!> \brief Numerical Derivatives for rotnuc
!> \param sepi ...
!> \param sepj ...
!> \param r ...
!> \param de1b ...
!> \param de2a ...
!> \param itype ...
!> \param delta ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE drotnuc_num(sepi,sepj,r,de1b,de2a,itype,delta,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: r
    REAL(dp), DIMENSION(3, 45), &
      INTENT(OUT), OPTIONAL                  :: de1b, de2a
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'drotnuc_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j
    LOGICAL                                  :: l_de1b, l_de2a
    REAL(dp)                                 :: od
    REAL(dp), DIMENSION(3)                   :: rr
    REAL(dp), DIMENSION(45)                  :: e1m, e1p, e2m, e2p

    l_de1b  = PRESENT(de1b)
    l_de2a  = PRESENT(de2a)
    IF (.NOT.(l_de1b.OR.l_de2a)) RETURN
    od = 0.5_dp/delta
    DO i=1,3
       rr=r
       rr(i)=rr(i)+delta
       CALL rotnuc_num(sepi,sepj,rr,e1p,e2p,itype,se_int_control,se_taper=se_taper,error=error)
       rr(i)=rr(i)-2._dp*delta
       CALL rotnuc_num(sepi,sepj,rr,e1m,e2m,itype,se_int_control,se_taper=se_taper,error=error)
       IF (l_de1b) THEN
          DO j = 1, sepi%atm_int_size
             de1b(i,j) = od * (e1p(j) - e1m(j))
          END DO
       END IF
       IF (l_de2a) THEN
          DO j = 1, sepj%atm_int_size
             de2a(i,j) = od * (e2p(j) - e2m(j))
          END DO
       END IF
    END DO
  END SUBROUTINE drotnuc_num

! *****************************************************************************
!> \brief Numerical Derivatives for corecore
!> \param sepi ...
!> \param sepj ...
!> \param r ...
!> \param denuc ...
!> \param itype ...
!> \param delta ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dcorecore_num(sepi,sepj,r,denuc,itype,delta,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: r
    REAL(dp), DIMENSION(3), INTENT(OUT)      :: denuc
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dcorecore_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: enucm, enucp, od
    REAL(dp), DIMENSION(3)                   :: rr

    od = 0.5_dp/delta
    DO i=1,3
       rr=r
       rr(i)=rr(i)+delta
       CALL corecore_num(sepi,sepj,rr,enucp,itype,se_int_control,se_taper=se_taper,error=error)
       rr(i)=rr(i)-2._dp*delta
       CALL corecore_num(sepi,sepj,rr,enucm,itype,se_int_control,se_taper=se_taper,error=error)
       denuc(i)  = od * ( enucp - enucm )
    END DO
  END SUBROUTINE dcorecore_num

! *****************************************************************************
!> \brief Numerical Derivatives for corecore
!> \param sepi ...
!> \param sepj ...
!> \param r ...
!> \param denuc ...
!> \param itype ...
!> \param delta ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dcorecore_el_num(sepi,sepj,r,denuc,itype,delta,se_int_control,se_taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: r
    REAL(dp), DIMENSION(3), INTENT(OUT)      :: denuc
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(IN)                     :: delta
    TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
    TYPE(se_taper_type), POINTER             :: se_taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dcorecore_el_num', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    REAL(dp)                                 :: enucm, enucp, od
    REAL(dp), DIMENSION(3)                   :: rr

    od = 0.5_dp/delta
    DO i=1,3
       rr=r
       rr(i)=rr(i)+delta
       CALL corecore_el_num(sepi,sepj,rr,enucp,itype,se_int_control,se_taper=se_taper,error=error)
       rr(i)=rr(i)-2._dp*delta
       CALL corecore_el_num(sepi,sepj,rr,enucm,itype,se_int_control,se_taper=se_taper,error=error)
       denuc(i)  = od * ( enucp - enucm )
    END DO
  END SUBROUTINE dcorecore_el_num

END MODULE semi_empirical_int_num
