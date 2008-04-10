!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Analytical derivatives of Integrals for semi-empirical methods
!> \author Teodoro Laino - Zurich University 04.2007 [tlaino]
!> \par History
!>      23.11.2007 jhu   short range version of integrals
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
! *****************************************************************************
MODULE semi_empirical_int_ana

  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_pchg,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             evolt
  USE semi_empirical_int_arrays,       ONLY: &
       inddd, inddp, indexa, indexb, indpp, int2c_type, map_22_to_491, &
       map_full_sp, pxy, rij_threshold
  USE semi_empirical_int_utils,        ONLY: dtaper_eval,&
                                             rot_2el_2c_first,&
                                             rotmat,&
                                             store_2el_2c_diag,&
                                             taper_eval
  USE semi_empirical_types,            ONLY: rotmat_create,&
                                             rotmat_release,&
                                             rotmat_type,&
                                             semi_empirical_type,&
                                             taper_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
#include "semi_empirical_int_debug.h"

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_int_ana'
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module=.FALSE.
  PUBLIC :: rotnuc_ana, rotint_ana, corecore_ana

CONTAINS

! *****************************************************************************
!> \brief Computes analytical gradients for semiempirical integrals
!> \param on input  sepi   = Atomic parameters of first atom
!>                  sepj   = Atomic parameters of second atom
!>                  rijv   = Coordinate vector i -> j 
!> 
!>    on output e1b,e2a= Array of electron-nuclear attraction integrals,
!>                       e1b = Electron on atom ni attracting nucleus of nj.
!>                       e2a = Electron on atom nj attracting nucleus of ni.
!>                       enuc= nuclear-nuclear repulsion term.
!> 
!>                      de1b = derivative of e1b term
!>                      de2a = derivative of e2a term
!>                      denuc= derivative of nuclear-nuclear repulsion term. 
!> \note
!>      Analytical version of the MOPAC rotnuc routine
!> \par History
!>      04.2007 created [tlaino]
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!>       Teodoro Laino [tlaino] - University of Zurich 04.2008 : removed the core-core part
!> \author Teodoro Laino [tlaino] - Zurich University
! *****************************************************************************
  RECURSIVE SUBROUTINE rotnuc_ana (sepi,sepj,rijv,itype,e1b,e2a,de1b,de2a,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), DIMENSION(45), INTENT(OUT), &
      OPTIONAL                               :: e1b, e2a
    REAL(dp), DIMENSION(45, 3), &
      INTENT(OUT), OPTIONAL                  :: de1b, de2a
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotnuc_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idd, idp, ind1, ind2, ipp, &
                                                j, last_orbital(2), m, n
    LOGICAL                                  :: invert, l_de1b, l_de2a, &
                                                l_e1b, l_e2a, lgrad, task(2)
    REAL(KIND=dp)                            :: rij, xtmp
    REAL(KIND=dp), DIMENSION(10, 2)          :: core, dcore
    REAL(KIND=dp), DIMENSION(3)              :: drij
    REAL(KIND=dp), DIMENSION(45)             :: tmp
    REAL(KIND=dp), DIMENSION(45, 3)          :: tmp_d
    TYPE(rotmat_type), POINTER               :: ij_matrix

    NULLIFY(ij_matrix)
    rij=DOT_PRODUCT(rijv,rijv)
    ! Initialization
    l_e1b  = PRESENT(e1b)
    l_e2a  = PRESENT(e2a)
    l_de1b = PRESENT(de1b)
    l_de2a = PRESENT(de2a)
    lgrad  = l_de1b.OR.l_de2a
    ! Zeros all arrays
    IF (l_e1b) THEN
       e1b =   0._dp
    END IF
    IF (l_e2a) THEN
       e2a =   0._dp
    END IF
    IF (l_de1b) THEN
       de1b =  0._dp
    END IF
    IF (l_de2a) THEN
       de2a =  0._dp
    END IF
    IF (rij > rij_threshold) THEN
       ! Compute Integrals in diatomic frame opportunely inverted
       rij  = SQRT(rij)
       ! Create the rotation matrix 
       CALL rotmat_create(ij_matrix, error)
       CALL rotmat (sepi, sepj, -rijv, rij, ij_matrix, do_derivatives=lgrad, do_invert=invert,&
            error=error)

       IF (lgrad) THEN
          drij(:) = rijv/rij
          ! Possibly Invert Frame
          IF (invert) THEN
             xtmp=drij(3)
             drij(3)=drij(2)
             drij(2)=xtmp
          END IF
       END IF
       
       CALL dnucint_ana(sepi,sepj,rij,core=core,dcore=dcore,itype=itype,taper=taper,&
            shortrange=shortrange,error=error)
      
       ! Copy parameters over to arrays for do loop.
       last_orbital(1)  = sepi%natorb
       last_orbital(2)  = sepj%natorb
       task (1)         = l_e1b
       task (2)         = l_e2a
       DO n = 1, 2
          IF (.NOT.task(n)) CYCLE
          tmp = 0.0_dp
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
          IF (n==1) e1b = -tmp
          IF (n==2) e2a = -tmp
       END DO
       IF (invert.AND.l_e1b)           CALL invert_integral  ( e1b,3,error)
       IF (invert.AND.l_e2a)           CALL invert_integral  ( e2a,3,error)

       ! Possibly compute derivatives
       task (1)         = l_de1b
       task (2)         = l_de2a
       DO n = 1, 2
          IF (.NOT.task(n)) CYCLE
          tmp_d = 0.0_dp
          DO i = 1, last_orbital(n)
             ind1 = i - 1
             DO j = 1, i
                ind2 = j - 1
                m = (i*(i-1))/2 + j
                ! Perform Rotations ...
                IF (ind2 == 0) THEN
                   IF (ind1 == 0) THEN
                      ! Type of Integral (SS/)
                      tmp_d(m,:) = dcore(1, n)*drij
                   ELSE IF (ind1 < 4) THEN
                      ! Type of Integral (SP/)
                      tmp_d(m,:) = ij_matrix%sp_d(1, ind1, :) *  core(2, n)           + &
                                   ij_matrix%sp  (1, ind1   ) * dcore(2, n) * drij
                   ELSE
                      ! Type of Integral (SD/)
                      tmp_d(m,:) = ij_matrix%sd_d(1, ind1-3, :) *  core(5, n)         + &
                                   ij_matrix%sd  (1, ind1-3   ) * dcore(5, n) * drij
                   END IF
                ELSE IF (ind2 < 4) THEN
                   IF (ind1 < 4) THEN
                      ! Type of Integral (PP/)
                      ipp = indpp(ind1, ind2)
                      tmp_d(m,:) = dcore(3, n) * drij *  ij_matrix%pp  (ipp, 1, 1)    + &
                                    core(3, n)        *  ij_matrix%pp_d(ipp, 1, 1, :) + &
                                   dcore(4, n) * drij * (ij_matrix%pp  (ipp, 2, 2   )+ij_matrix%pp  (ipp, 3, 3   )) + &
                                    core(4, n) *        (ij_matrix%pp_d(ipp, 2, 2, :)+ij_matrix%pp_d(ipp, 3, 3, :))
                   ELSE
                      ! Type of Integral (PD/)
                      idp = inddp(ind1-3, ind2)
                      tmp_d(m,:) = dcore(6, n) * drij *  ij_matrix%pd  (idp, 1, 1   ) + &
                                    core(6, n) *         ij_matrix%pd_d(idp, 1, 1, :) + &
                                   dcore(8, n) * drij * (ij_matrix%pd  (idp, 2, 2   )+ij_matrix%pd  (idp, 3, 3   )) + &
                                    core(8, n) *        (ij_matrix%pd_d(idp, 2, 2, :)+ij_matrix%pd_d(idp, 3, 3, :)) 
                   END IF
                ELSE
                   ! Type of Integral (DD/)
                   idd = inddd(ind1-3, ind2-3)
                   tmp_d(m,:) = dcore( 7, n) * drij *  ij_matrix%dd  (idd, 1, 1   ) + &
                                 core( 7, n) *         ij_matrix%dd_d(idd, 1, 1, :) + &
                                dcore( 9, n) * drij * (ij_matrix%dd  (idd, 2, 2   )+ij_matrix%dd  (idd, 3, 3   )) + &
                                 core( 9, n) *        (ij_matrix%dd_d(idd, 2, 2, :)+ij_matrix%dd_d(idd, 3, 3, :)) + &
                                dcore(10, n) * drij * (ij_matrix%dd  (idd, 4, 4   )+ij_matrix%dd  (idd, 5, 5   )) + &
                                 core(10, n) *        (ij_matrix%dd_d(idd, 4, 4, :)+ij_matrix%dd_d(idd, 5, 5, :))
                END IF
             END DO
          END DO
          IF (n==1) de1b = -tmp_d
          IF (n==2) de2a = -tmp_d
       END DO
       IF (invert.AND.l_de1b) CALL invert_derivative(de1b,3,error)
       IF (invert.AND.l_de2a) CALL invert_derivative(de2a,3,error)
       CALL rotmat_release(ij_matrix, error)

       ! Possibly debug the analytical values versus the numerical ones
       IF (debug_this_module) THEN
          CALL check_drotnuc_ana(sepi, sepj, rijv, itype,shortrange, taper, e1b, e2a, de1b, de2a, error)
       END IF
    END IF
  END SUBROUTINE rotnuc_ana

! *****************************************************************************
!> \brief Computes analytical gradients for semiempirical integrals
!> \param on input  sepi   = Atomic parameters of first atom
!>                  sepj   = Atomic parameters of second atom
!>                  rijv   = Coordinate vector i -> j 
!> 
!>    on output e1b,e2a= Array of electron-nuclear attraction integrals,
!>                       e1b = Electron on atom ni attracting nucleus of nj.
!>                       e2a = Electron on atom nj attracting nucleus of ni.
!>                       enuc= nuclear-nuclear repulsion term.
!> 
!>                      de1b = derivative of e1b term
!>                      de2a = derivative of e2a term
!>                      denuc= derivative of nuclear-nuclear repulsion term. 
!> \note
!>      Analytical version of the MOPAC rotnuc routine
!> \par History
!>      04.2007 created [tlaino]
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!>       Teodoro Laino [tlaino] - University of Zurich 04.2008 : removed the core-core part
!> \author Teodoro Laino [tlaino] - Zurich University
! *****************************************************************************
  RECURSIVE SUBROUTINE corecore_ana (sepi,sepj,rijv,itype,enuc,denuc,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
    INTEGER, INTENT(IN)                      :: itype
    REAL(dp), INTENT(OUT), OPTIONAL          :: enuc
    REAL(dp), DIMENSION(3), INTENT(OUT), &
      OPTIONAL                               :: denuc
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'corecore_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ig, nt
    LOGICAL                                  :: failure, l_denuc, l_enuc
    REAL(dp) :: a, aab, alpi, alpj, apdg, ax, dai, daj, dax, dbi, dbj, &
      denuc_loc, dqcorr, drija, dscale, dssss, dzz, enuc_loc, pai, paj, pbi, &
      pbj, poi, poj, qcorr, rij, rija, scale, ssss, xab, xtmp, zaf, zbf, zz
    REAL(dp), DIMENSION(3)                   :: drij, z
    REAL(dp), DIMENSION(4)                   :: fni1, fni2, fni3, fnj1, fnj2, &
                                                fnj3

    failure = .FALSE.
    rij=DOT_PRODUCT(rijv,rijv)
    ! Initialization
    l_enuc = PRESENT(enuc)
    l_denuc= PRESENT(denuc)
    IF (l_enuc) THEN
       enuc =  0._dp
    END IF
    IF (l_denuc) THEN
       denuc = 0._dp
    END IF
    IF (rij > rij_threshold) THEN
       ! Compute Integrals in diatomic frame opportunely inverted
       rij  = SQRT(rij)

       IF (l_enuc.OR.l_denuc) THEN
          CALL dnucint_ana(sepi,sepj,rij,ssss=ssss,dssss=dssss,itype=itype,taper=taper,&
               shortrange=shortrange,error=error)

          ! Zeroing local method dependent core-core corrections
          qcorr = 0.0_dp
          scale = 0.0_dp
          dscale= 0.0_dp
          dqcorr= 0.0_dp
          ! Method dependent code
          IF(itype /= do_method_pm6) THEN
             alpi = sepi%alp
             alpj = sepj%alp
             scale = EXP(-alpi*rij)+EXP(-alpj*rij)
             IF (l_denuc) THEN
                dscale= -alpi* EXP(-alpi*rij)-alpj*EXP(-alpj*rij)
             END IF
             nt=sepi%z+sepj%z
             IF(nt == 8 .OR. nt == 9) THEN
                IF(sepi%z == 7.OR.sepi%z == 8) THEN
                   scale=scale+(0.5291772083_dp*rij-1._dp)*EXP(-alpi*rij)
                   IF (l_denuc) THEN
                      dscale=dscale+0.5291772083_dp*EXP(-alpi*rij)-&
                           (0.5291772083_dp*rij-1._dp)*alpi*EXP(-alpi*rij)
                   END IF
                END IF
                IF(sepj%z == 7.OR.sepj%z == 8) THEN
                   scale=scale+(0.5291772083_dp*rij-1._dp)*EXP(-alpj*rij)
                   IF (l_denuc) THEN
                      dscale=dscale+0.5291772083_dp*EXP(-alpj*rij)-&
                           (0.5291772083_dp*rij-1._dp)*alpj*EXP(-alpj*rij)
                   END IF
                END IF
             ENDIF
             zz = sepi%zeff*sepj%zeff
             enuc_loc = zz*ssss
             IF (l_denuc) THEN
                denuc_loc= zz*dssss
                dscale=SIGN(1.0_dp,scale*enuc_loc)*(dscale*enuc_loc+scale*denuc_loc)
                dzz=-zz/rij**2
             END IF
             scale=ABS(scale*enuc_loc)
             zz=zz/rij
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
                         IF (l_denuc) THEN
                            dax = fni2(ig)*2.0_dp*(rij-fni3(ig))
                            dscale = dscale + dzz*fni1(ig)*EXP(-ax)-dax*zz*fni1(ig)*EXP(-ax)
                         END IF
                      ENDIF
                   ENDIF
                   IF(ABS(fnj1(ig)) > 0._dp) THEN
                      ax = fnj2(ig)*(rij-fnj3(ig))**2
                      IF(ax <= 25._dp) THEN
                         scale=scale +zz*fnj1(ig)*EXP(-ax)
                         IF (l_denuc) THEN
                            dax = fnj2(ig)*2.0_dp*(rij-fnj3(ig))
                            dscale = dscale + dzz*fnj1(ig)*EXP(-ax)-dax*zz*fnj1(ig)*EXP(-ax)
                         END IF
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
                apdg=10._dp*0.5291772083_dp*0.5291772083_dp
                qcorr = (zaf*pai+zbf*paj)*EXP(-apdg*(rij-dai-daj)**2)+ &
                        (zaf*pai+zbf*pbj)*EXP(-apdg*(rij-dai-dbj)**2)+ &
                        (zaf*pbi+zbf*paj)*EXP(-apdg*(rij-dbi-daj)**2)+ &
                        (zaf*pbi+zbf*pbj)*EXP(-apdg*(rij-dbi-dbj)**2)
                IF (l_denuc) THEN
                   dqcorr = (zaf*pai+zbf*paj)*EXP(-apdg*(rij-dai-daj)**2)*(-2.0_dp*apdg*(rij-dai-daj))+ &
                            (zaf*pai+zbf*pbj)*EXP(-apdg*(rij-dai-dbj)**2)*(-2.0_dp*apdg*(rij-dai-dbj))+ &
                            (zaf*pbi+zbf*paj)*EXP(-apdg*(rij-dbi-daj)**2)*(-2.0_dp*apdg*(rij-dbi-daj))+ &
                            (zaf*pbi+zbf*pbj)*EXP(-apdg*(rij-dbi-dbj)**2)*(-2.0_dp*apdg*(rij-dbi-dbj))
                END IF
             ELSEIF (itype==do_method_pchg) THEN
                qcorr = 0.0_dp
                scale = 0.0_dp
                dscale= 0.0_dp
                dqcorr= 0.0_dp
             ELSE
                qcorr = 0.0_dp
                dqcorr= 0.0_dp
             END IF
          ELSE
             ! PM6 core-core terms
             IF ( ABS(sepi%rho) > EPSILON(0.0_dp) ) THEN
               poi = sepi%rho
             ELSE
               poi = sepi%ko(1)
             END IF
             IF ( ABS(sepj%rho) > EPSILON(0.0_dp) ) THEN
               poj = sepj%rho
             ELSE
               poj = sepj%ko(1)
             END IF
             scale = sepi%zeff*sepj%zeff/SQRT(rij*rij+(poi+poj)**2)
             IF(shortrange) scale = scale - sepi%zeff*sepj%zeff/rij
             IF (l_denuc) THEN
               dscale = -rij*sepi%zeff*sepj%zeff/(rij*rij+(poi+poj)**2)**(1.5_dp)
               IF(shortrange) dscale = dscale + sepi%zeff*sepj%zeff/rij/rij
             END IF
             drija  = angstrom
             rija   = rij*drija
             xab    = sepi%xab(sepj%z)
             aab    = sepi%aab(sepj%z)
             IF ( (sepi%z == 1 .AND. (sepj%z == 6 .OR. sepj%z == 7 .OR. sepj%z == 8)) .OR. &
                  (sepj%z == 1 .AND. (sepi%z == 6 .OR. sepi%z == 7 .OR. sepi%z == 8)) ) THEN
                ! Special Case O-H or N-H or C-H
                enuc_loc=scale*(1._dp+2._dp*xab*EXP(-aab*rija*rija))
                IF (l_denuc) THEN
                   denuc_loc = dscale*(1._dp+2._dp*xab*EXP(-aab*rija*rija))-&
                        scale*2._dp*xab*EXP(-aab*rija*rija)*(2.0_dp*aab*rija)*drija
                END IF
             ELSEIF ( sepi%z == 6 .AND. sepj%z == 6 ) THEN
                ! Special Case C-C
                enuc_loc=scale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))+9.28_dp*EXP(-5.98*rija))
                IF (l_denuc) THEN
                   denuc_loc = dscale * (1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))+9.28_dp*EXP(-5.98*rija))&
                        -scale*2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))*(aab+6.0_dp*0.0003_dp*rija**5)*drija &
                        -scale*9.28_dp*EXP(-5.98*rija)*5.98_dp*drija
                END IF
             ELSEIF ( (sepi%z == 8 .AND. sepj%z == 14) .OR. &
                  (sepj%z == 8 .AND. sepi%z == 14) ) THEN
                ! Special Case Si-O
                enuc_loc=scale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))-0.0007_dp*EXP(-(rija-2.9_dp)**2))
                IF (l_denuc) THEN
                   denuc_loc = dscale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))-0.0007_dp*EXP(-(rija-2.9_dp)**2))&
                        -scale*2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))*(aab+6.0_dp*0.0003_dp*rija**5)*drija+&
                        scale*0.0007_dp*EXP(-(rija-2.9_dp)**2)*(2.0_dp*(rija-2.9_dp)*drija)
                END IF
             ELSE
                ! General Case
                ! Factor of 2 found by experiment 
                enuc_loc=scale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6)))
                IF (l_denuc) THEN
                   denuc_loc = dscale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6)))&
                        -scale*2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))*(aab+6.0_dp*0.0003_dp*rija**5)*drija
                END IF
             END IF
             ! General correction term a*exp(-b*(rij-c)^2)
             enuc_loc=enuc_loc+(sepi%a*EXP(-sepi%b*(rij-sepi%c)**2))*sepi%zeff*sepj%zeff/rij
             IF (l_denuc) THEN
                denuc_loc=denuc_loc+(sepi%a*EXP(-sepi%b*(rij-sepi%c)**2)*(-2.0_dp*sepi%b*(rij-sepi%c)))&
                     *sepi%zeff*sepj%zeff/rij-(sepi%a*EXP(-sepi%b*(rij-sepi%c)**2))*sepi%zeff*sepj%zeff/rij**2
             END IF
             enuc_loc=enuc_loc+(sepj%a*EXP(-sepj%b*(rij-sepj%c)**2))*sepi%zeff*sepj%zeff/rij
             IF (l_denuc) THEN
                denuc_loc=denuc_loc+(sepj%a*EXP(-sepj%b*(rij-sepj%c)**2)*(-2.0_dp*sepj%b*(rij-sepj%c)))&
                     *sepi%zeff*sepj%zeff/rij-(sepj%a*EXP(-sepj%b*(rij-sepj%c)**2))*sepi%zeff*sepj%zeff/rij**2
             END IF
             ! Hard core repulsion
             xtmp = 1.e-8_dp/evolt * ( (REAL(sepi%z,dp)**(1._dp/3._dp) + REAL(sepj%z,dp)**(1._dp/3._dp))/rija )**12
             enuc_loc=enuc_loc+xtmp
             IF (l_denuc) THEN
                denuc_loc=denuc_loc-12.0_dp*xtmp/rija*drija
             END IF
             scale  = 0.0_dp
             dscale = 0.0_dp
          ENDIF

          ! Only at the very end let's sum-up the several contributions energy/derivatives
          ! This assignment should be method indipendent
          IF (l_enuc) THEN
             enuc=enuc_loc+scale+qcorr
          END IF
          IF (l_denuc) THEN
             drij(:) = rijv/rij
             denuc = (denuc_loc + dscale + dqcorr)*drij
          END IF
          IF (debug_this_module) THEN
             CALL check_dcorecore_ana(sepi, sepj, rijv, itype,shortrange, taper, enuc, denuc, error)
          END IF
       ENDIF
    END IF
  END SUBROUTINE corecore_ana

! *****************************************************************************
!> \brief Exploits inversion symmetry to avoid divergence
!> \par History
!>      04.2007 created [tlaino]
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE invert_integral(array, np, error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: array
    INTEGER, INTENT(IN)                      :: np
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'invert_integral', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i1, i2
    INTEGER, DIMENSION(2, 42)                :: map
    REAL(KIND=dp)                            :: tmp

    map = RESHAPE(map_full_sp,(/2,42/))
    DO i = 1, np
       i1= map(1,i)
       i2= map(2,i)

       tmp=array(i1)
       array(i1)=array(i2)
       array(i2)=tmp
    END DO

  END SUBROUTINE invert_integral

! *****************************************************************************
!> \brief Exploits inversion symmetry to avoid divergence
!> \par History
!>      04.2007 created [tlaino]
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE invert_derivative(array, np, error)
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: array
    INTEGER, INTENT(IN)                      :: np
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'invert_derivative', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i1, i2, j, m
    INTEGER, DIMENSION(2, 42)                :: map
    REAL(KIND=dp)                            :: tmp

    map = RESHAPE(map_full_sp,(/2,42/))
    DO j=1,3
       DO i =1, np
          i1= map(1,i)
          i2= map(2,i)

          tmp=array(i1,j)
          array(i1,j)=array(i2,j)
          array(i2,j)=tmp
       END DO
    END DO

    DO m=1,SIZE(array,1)
       tmp=array(m,2)
       array(m,2)=array(m,3)
       array(m,3)=tmp
    ENDDO
  END SUBROUTINE invert_derivative

! *****************************************************************************
!> \brief Calculates the nuclear attraction integrals and analytical integrals (main driver)
!> \param on input rij     = interatomic distance
!>                sepi    = paramters of atom i
!>                sepj    = paramters of atom j
!>                itype   = type of semi_empirical model
!>                          extension to the original routine to compute qm/mm integrals
!>      on output core    = derivative of 4 X 2 array of electron-core attraction integrals
!>                ssss    = derivative of (ssss) integral
!>                          derivatives are intended w.r.t. rij
!>         The storage of the nuclear attraction integrals  core(kl/ij) iS
!>         (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
!>         where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
!> \note
!>      Analytical version - Analytical evaluation of gradients
!>      Teodoro Laino - Zurich University 04.2007
!>      
!> \par History
!>      03.2007 created [tlaino]
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  RECURSIVE SUBROUTINE dnucint_ana ( sepi, sepj, rij, ssss, core, dssss, dcore, itype,&
       taper, shortrange, do_longrange, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(OUT), OPTIONAL          :: ssss
    REAL(dp), DIMENSION(10, 2), &
      INTENT(OUT), OPTIONAL                  :: core
    REAL(dp), INTENT(OUT), OPTIONAL          :: dssss
    REAL(dp), DIMENSION(10, 2), &
      INTENT(OUT), OPTIONAL                  :: dcore
    INTEGER, INTENT(IN)                      :: itype
    TYPE(taper_type), POINTER                :: taper
    LOGICAL, INTENT(IN)                      :: shortrange
    LOGICAL, INTENT(IN), OPTIONAL            :: do_longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dnucint_ana', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, l_core, l_ssss, &
                                                longrange
    REAL(KIND=dp)                            :: dft, dssss_lr, ft, ssss_lr
    REAL(KIND=dp), DIMENSION(10, 2)          :: core_lr, dcore_lr

    failure = .FALSE.


    longrange = .FALSE.
    IF (PRESENT(do_longrange))  longrange = do_longrange
    CPPrecondition(PRESENT(core).EQV.PRESENT(dcore),cp_failure_level,routineP,error,failure)
    CPPrecondition(PRESENT(ssss).EQV.PRESENT(dssss),cp_failure_level,routineP,error,failure)
    l_core  = PRESENT(core)
    l_ssss  = PRESENT(ssss)
    IF (l_core) THEN
       core  = 0.0_dp
       dcore = 0.0_dp
    END IF
    IF (l_ssss) THEN
       ssss  = 0.0_dp
       dssss = 0.0_dp
    END IF

    CALL dnucint_ana_spcore(sepi, sepj, rij, ssss, core, dssss, dcore, itype, taper,&
         longrange,error=error)

    IF (sepi%dorb.OR.sepj%dorb) THEN
       ! compute the contribution from d-orbitals
       STOP "D-ORBITALS: dnucint_ana"
    END IF

    IF (.NOT.longrange) THEN
       ! Tapering function
       ft = 1.0_dp
       dft= 0.0_dp
       IF (itype/=do_method_pchg) THEN
          ft = taper_eval (taper, rij, error)
          dft= dtaper_eval (taper, rij, error)
       END IF

       IF (l_core) THEN
          dcore(:,:) = ft*dcore(:,:)+dft*core(:,:)
          core(:,:)  = ft*core(:,:)
       END IF
       IF (l_ssss) THEN
          dssss      = ft*dssss+dft*ssss
          ssss       = ft*ssss
       END IF
    END IF

    ! Modify integrals in case of short-range request
    IF ( shortrange ) THEN
       IF (l_ssss) THEN
          CALL dnucint_ana( sepi, sepj, rij, ssss=ssss_lr, dssss=dssss_lr, itype=itype,&
               taper=taper, shortrange=.FALSE., do_longrange=.TRUE., error=error )
           ssss =  ssss -  ssss_lr
          dssss = dssss - dssss_lr
       END IF
       IF (l_core) THEN
          CALL dnucint_ana( sepi, sepj, rij, core=core_lr, dcore=dcore_lr, itype=itype,&
               taper=taper, shortrange=.FALSE., do_longrange=.TRUE., error=error )
           core =  core -  core_lr
          dcore = dcore - dcore_lr
       END IF
    END IF


    IF (l_core) THEN
       ! If no orbitals are present on the atom let's zero the integral contribution
       IF (sepi%natorb == -1) dcore(:,1) = 0.0_dp
       IF (sepj%natorb == -1) dcore(:,2) = 0.0_dp
       
       ! If no orbitals are present on the atom let's zero the integral contribution
       IF (sepi%natorb == -1) core(:,1) = 0.0_dp
       IF (sepj%natorb == -1) core(:,2) = 0.0_dp
    END IF

    ! Debug Procedure.. Check valifity of analytical gradients of nucint
    IF (debug_this_module) THEN
       CALL check_dnucint_ana(sepi,sepj,rij,dssss,dcore,itype,shortrange,taper=taper,&
            error=error)
    END IF

  END SUBROUTINE dnucint_ana

! *****************************************************************************
!> \brief Calculates the nuclear attraction integrals for sp basis
!> \param on input rij     = interatomic distance
!>                sepi    = paramters of atom i
!>                sepj    = paramters of atom j
!>                itype   = type of semi_empirical model
!>                          extension to the original routine to compute qm/mm integrals
!>      on output core    = derivative of 4 X 2 array of electron-core attraction integrals
!>                ssss    = derivative of (ssss) integral
!>                          derivatives are intended w.r.t. rij
!>         The storage of the nuclear attraction integrals  core(kl/ij) iS
!>         (SS/)=1,   (SO/)=2,   (OO/)=3,   (PP/)=4
!>         where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
!> \note
!>      Analytical version - Analytical evaluation of gradients
!>      Teodoro Laino - Zurich University 04.2007
!>      routine adapted from mopac7 (repp)
!>      vector version written by Ernest R. Davidson, Indiana University
!> \par History
!>      04.2007 created [tlaino]
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE dnucint_ana_spcore ( sepi, sepj, rij, ssss, core, dssss, dcore, itype,&
       taper, longrange, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(INOUT), OPTIONAL        :: ssss
    REAL(dp), DIMENSION(10, 2), &
      INTENT(INOUT), OPTIONAL                :: core
    REAL(dp), INTENT(INOUT), OPTIONAL        :: dssss
    REAL(dp), DIMENSION(10, 2), &
      INTENT(INOUT), OPTIONAL                :: dcore
    INTEGER, INTENT(IN)                      :: itype
    TYPE(taper_type), POINTER                :: taper
    LOGICAL, INTENT(IN)                      :: longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dnucint_ana_spcore', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, l_core, l_ssss, si, &
                                                sj
    REAL(dp)                                 :: ade, adi, adj, aed, aee, aeq, &
                                                ami, amj, aqe, aqi, aqj, da, &
                                                db, ddi, ddj, qa, qb, qqi, &
                                                qqj, r, rsq, xxx, zi, zj
    REAL(dp), DIMENSION(7)                   :: arg, darg, dri, dsqr, ri, sqr

    failure   = .FALSE.
    CPPrecondition(PRESENT(core).EQV.PRESENT(dcore),cp_failure_level,routineP,error,failure)
    CPPrecondition(PRESENT(ssss).EQV.PRESENT(dssss),cp_failure_level,routineP,error,failure)
    l_core  = PRESENT(core)
    l_ssss  = PRESENT(ssss)

    IF (.NOT.(l_core.OR.l_ssss)) RETURN
    ri  = 0.0_dp
    dri = 0.0_dp
    r   = rij
    si  = (sepi%natorb > 1)
    sj  = (sepj%natorb > 1)
    zi  = sepi%zeff
    zj  = sepj%zeff
    ddi = sepi%cs(2)
    ddj = sepj%cs(2)
    qqi = sepi%cs(3)
    qqj = sepj%cs(3)

    !  hydrogen - hydrogen  (SS/SS)
    aee = 0.0_dp
    IF (itype/=do_method_pchg) THEN
       ami = sepi%ko(1)
       amj = sepj%ko(1)
       aee = ami + amj
       aee = aee * aee
    END IF
    IF (longrange)  aee = 0.0_dp
    rsq     = r**2
    arg(1)  = rsq + aee
    sqr(1)  = pxy(1)/SQRT(arg(1))
    darg(1) = 2.0_dp *r
    dsqr(1) = -(0.5_dp/arg(1))*darg(1)*sqr(1)
    ri(1)      = sqr(1)
    dri(1)     = dsqr(1)
    ! Store the value for (SS|SS)
    IF (l_ssss) THEN
       ssss   =  ri(1)
       dssss  = dri(1)
    END IF

    IF (l_core) THEN
       core(1,1)  = zj*ri(1)
       core(1,2)  = zi*ri(1)
       dcore(1,1) = zj*dri(1)
       dcore(1,2) = zi*dri(1)
       IF (si .OR.sj) THEN
          ! In case (i) is an heavy atom
          IF (si) THEN
             da  = ddi
             qa  = qqi
             ade = 0.0_dp
             aqe = 0.0_dp
             IF (itype /=do_method_pchg) THEN
                adi = sepi%ko(2)
                aqi = sepi%ko(3)
                ade = adi + amj
                aqe = aqi + amj
                ade = ade * ade
                aqe = aqe * aqe
             END IF
             IF (longrange) THEN
                ade = 0.0_dp
                aqe = 0.0_dp
             END IF
             xxx = r+da
             arg(2) = xxx*xxx + ade
             darg(2)= 2.0_dp *xxx
             xxx = r-da
             arg(3) = xxx*xxx + ade
             darg(3)= 2.0_dp *xxx
             xxx = r+qa
             arg(4) = xxx*xxx + aqe
             darg(4)= 2.0_dp *xxx
             xxx = r-qa
             arg(5) = xxx*xxx + aqe
             darg(5)= 2.0_dp *xxx
             arg(6) = rsq + aqe
             darg(6)= 2.0_dp *r
             arg(7) = arg(6) + qa*qa
             darg(7)= darg(6)
             sqr(2:7) = pxy(2:7)/SQRT(arg(2:7))
             dsqr     = -(0.5_dp/arg)*darg*sqr
             ri(2)      = sqr(2) + sqr(3)
             ri(3)      = sqr(1) + sqr(4) + sqr(5) + sqr(6)
             ri(4)      = sqr(1) + sqr(6) + sqr(7)
             core(2,1)  = zj*ri(2)
             core(3,1)  = zj*ri(3)
             core(4,1)  = zj*ri(4)
             dri(2)     = dsqr(2) + dsqr(3)
             dri(3)     = dsqr(1) + dsqr(4) + dsqr(5) + dsqr(6)
             dri(4)     = dsqr(1) + dsqr(6) + dsqr(7)
             dcore(2,1) = zj*dri(2)
             dcore(3,1) = zj*dri(3)
             dcore(4,1) = zj*dri(4)
          END IF
          
          ! In case (j) is an heavy atom
          IF (sj) THEN
             db  = ddj
             qb  = qqj
             aed = 0.0_dp
             aeq = 0.0_dp
             IF (itype /=do_method_pchg) THEN
                adj = sepj%ko(2)
                aqj = sepj%ko(3)
                aed = ami + adj
                aeq = ami + aqj
                aed = aed * aed
                aeq = aeq * aeq
             END IF
             IF (longrange) THEN
                ade = 0.0_dp
                aqe = 0.0_dp
             END IF
             xxx = r-db
             arg(2) = xxx*xxx + aed
             darg(2)= 2.0_dp *xxx
             xxx = r+db
             arg(3) = xxx*xxx + aed
             darg(3)= 2.0_dp *xxx
             xxx = r-qb
             arg(4) = xxx*xxx + aeq
             darg(4)= 2.0_dp *xxx
             xxx = r+qb
             arg(5) = xxx*xxx + aeq
             darg(5)= 2.0_dp *xxx
             arg(6) = rsq + aeq
             darg(6)= 2.0_dp *r
             arg(7) = arg(6) + qb*qb
             darg(7)= darg(6)
             sqr(2:7) = pxy(2:7)/SQRT(arg(2:7))
             dsqr     = -(0.5_dp/arg)*darg*sqr
             ri(2)      = sqr(2) + sqr(3)
             ri(3)      = sqr(1) + sqr(4) + sqr(5) + sqr(6)
             ri(4)      = sqr(1) + sqr(6) + sqr(7)
             core(2,2)  = zi*ri(2)
             core(3,2)  = zi*ri(3)
             core(4,2)  = zi*ri(4)
             dri(2)     = dsqr(2) + dsqr(3)
             dri(3)     = dsqr(1) + dsqr(4) + dsqr(5) + dsqr(6)
             dri(4)     = dsqr(1) + dsqr(6) + dsqr(7)
             dcore(2,2) = zi*dri(2)
             dcore(3,2) = zi*dri(3)
             dcore(4,2) = zi*dri(4)
          END IF
       END IF
    END IF

  END SUBROUTINE dnucint_ana_spcore

! *****************************************************************************
!> \brief calculates the derivative of the two-particle interactions
!> \param on input  sepi   = Atomic parameters of first atom
!>               sepj   = Atomic parameters of second atom
!>               rijv   = Coordinate vector i -> j 
!> 
!>     on output w      = Array of two-electron repulsion integrals.
!> \note
!>      Analytical version - Analytical evaluation of gradients
!>      Teodoro Laino - Zurich University 04.2007
!>      routine adapted from mopac7 (repp)
!>      vector version written by Ernest R. Davidson, Indiana University
!> \par History
!>      04.2007 created [tlaino]
!>      Teodoro Laino (04.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  RECURSIVE SUBROUTINE rotint_ana (sepi,sepj,rijv,w,dw,shortrange,taper,error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
    REAL(dp), DIMENSION(2025), INTENT(OUT), &
      OPTIONAL                               :: w
    REAL(dp), DIMENSION(2025, 3), &
      INTENT(OUT), OPTIONAL                  :: dw
    LOGICAL, INTENT(IN)                      :: shortrange
    TYPE(taper_type), POINTER                :: taper
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rotint_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i1, ii, ij, ij1, iminus, &
                                                istep, iw_loc, j, j1, jj, k, &
                                                kk, kl, l, limij, limkl, mm
    LOGICAL                                  :: invert, l_w, lgrad
    LOGICAL, DIMENSION(45, 45)               :: logv, logv_d
    REAL(dp)                                 :: rij, xtmp
    REAL(dp), DIMENSION(3)                   :: drij
    REAL(KIND=dp)                            :: cc, cc_d(3), wrepp, wrepp_d(3)
    REAL(KIND=dp), DIMENSION(2025)           :: ww
    REAL(KIND=dp), DIMENSION(2025, 3)        :: ww_d
    REAL(KIND=dp), DIMENSION(45, 45)         :: v
    REAL(KIND=dp), DIMENSION(45, 45, 3)      :: v_d
    REAL(KIND=dp), DIMENSION(491)            :: rep, rep_d
    TYPE(rotmat_type), POINTER               :: ij_matrix

    NULLIFY(ij_matrix)
    l_w   = PRESENT(w)
    lgrad = PRESENT(dw)
    IF (l_w)    w = 0.0_dp
    IF (lgrad) dw = 0.0_dp
    IF (.NOT.(l_w.OR.lgrad)) RETURN
    rij=DOT_PRODUCT(rijv, rijv)
    IF (rij > rij_threshold) THEN
       ! The repulsion integrals over molecular frame (w) are stored in the
       ! order in which they will later be used.  ie.  (i,j/k,l) where
       ! j.le.i  and  l.le.k     and l varies most rapidly and i least
       ! rapidly.  (anti-normal computer storage)
       rij = SQRT(rij)

       ! Create the rotation matrix 
       CALL rotmat_create(ij_matrix, error)
       CALL rotmat (sepi, sepj, -rijv, rij, ij_matrix, do_derivatives=lgrad, do_invert=invert, error=error)          

       ! Compute integrals in diatomic frame as well their derivatives (if requested)
       CALL dterep_ana(sepi,sepj,rij,rep,rep_d,taper,shortrange,error=error)

       IF (lgrad) THEN
          drij(:) = rijv/rij
          ! Possibly Invert Frame
          IF (invert) THEN
             xtmp=drij(3)
             drij(3)=drij(2)
             drij(2)=xtmp
          END IF
       END IF

       ii = sepi%natorb
       kk = sepj%natorb
       ! First step in rotation of integrals
       CALL rot_2el_2c_first(sepi, sepj, rijv, shortrange, taper, invert, ii, kk, rep, logv, ij_matrix,&
            v, lgrad, rep_d, v_d, logv_d, drij, error)

       ! Integrals if requested
       IF (l_w) THEN
          ! Rotate Integrals
          IF (ii*kk > 0) THEN
             limij = indexb(ii, ii)
             limkl = indexb(kk, kk)
             istep = limkl * limij
             ww(1:istep) = 0.0_dp
             
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
             ! Store two electron integrals in the triangular format
             limij = (ii*(ii+1))/2
             limkl = (kk*(kk+1))/2
             CALL store_2el_2c_diag(limij, limkl, ww(1:istep), w, error=error)
             IF (invert) CALL invert_integral(w,42,error)
          END IF

          IF (debug_this_module) THEN
             ! Check value of integrals
             CALL check_rotint_ana(sepi,sepj,rijv,w,shortrange=shortrange,taper=taper,error=error)
          END IF
       END IF

       ! Gradients if requested
       IF (lgrad) THEN
          dw = 0.0_dp
          ! Rotate Integrals derivatives
          IF (ii*kk > 0) THEN
             limij = indexb(ii, ii)
             limkl = indexb(kk, kk)
             istep = limkl * limij
             ww_d(1:istep,:) = 0.0_dp
             
             ! Second step in rotation of integrals
             DO i1 = 1, ii
                DO j1 = 1, i1
                   ij = indexa(i1, j1)
                   jj = indexb(i1, j1)
                   mm = int2c_type(jj)
                   DO k = 1, kk
                      DO l = 1, k
                         kl = indexb(k, l)
                         IF (logv_d(ij, kl)) THEN
                            wrepp_d = v_d(ij, kl, :)
                            wrepp   =   v(ij, kl)
                            SELECT CASE (mm)
                            CASE(1)
                               ! (SS/)
                               i = 1
                               j = 1
                               iw_loc = (indexb(i,j)-1) * limkl + kl
                               ww_d(iw_loc, :) = wrepp_d
                            CASE(2)
                               ! (SP/)
                               j = 1
                               DO i = 1, 3
                                  iw_loc = (indexb(i+1,j)-1) * limkl + kl
                                  ww_d(iw_loc, :) = ww_d(iw_loc, :) + ij_matrix%sp_d(i1-1, i,:) * wrepp +&
                                                                      ij_matrix%sp(i1-1, i) * wrepp_d
                               END DO
                            CASE(3)
                               ! (PP/)
                               DO i = 1, 3
                                  cc   = ij_matrix%pp(i, i1-1, j1-1)
                                  cc_d = ij_matrix%pp_d(i, i1-1, j1-1, :)
                                  iw_loc = (indexb(i+1,i+1)-1) * limkl + kl
                                  ww_d(iw_loc, :) = ww_d(iw_loc, :) + cc_d * wrepp + cc *wrepp_d
                                  iminus = i - 1
                                  IF (iminus /= 0) THEN
                                     DO j = 1, iminus
                                        cc   = ij_matrix%pp(1+i+j, i1-1, j1-1)
                                        cc_d = ij_matrix%pp_d(1+i+j, i1-1, j1-1, :)
                                        iw_loc = (indexb(i+1,j+1)-1) * limkl + kl
                                        ww_d(iw_loc, :) = ww_d(iw_loc, :) + cc_d * wrepp + cc * wrepp_d
                                     END DO
                                  END IF
                               END DO
                            CASE(4)
                               ! (SD/)
                               j = 1
                               DO i = 1, 5
                                  iw_loc = (indexb(i+4,j)-1) * limkl + kl
                                  ww_d(iw_loc, :) = ww_d(iw_loc, :) + ij_matrix%sd_d(i1-4, i, :) * wrepp +&
                                                                      ij_matrix%sd(i1-4, i)      * wrepp_d
                               END DO
                            CASE(5)
                               ! (DP/)
                               DO i = 1, 5
                                  DO j = 1, 3
                                     iw_loc = (indexb(i+4,j+1)-1) * limkl + kl
                                     ij1 = 3 * (i-1) + j
                                     ww_d(iw_loc, :) = ww_d(iw_loc, :) + ij_matrix%pd_d(ij1, i1-4, j1-1, :) * wrepp +&
                                                                         ij_matrix%pd(ij1, i1-4, j1-1)      * wrepp_d
                                  END DO
                               END DO
                            CASE(6)
                               ! (DD/)
                               DO i = 1, 5
                                  cc   = ij_matrix%dd(i, i1-4, j1-4)
                                  cc_d = ij_matrix%dd_d(i, i1-4, j1-4, :)
                                  iw_loc = (indexb(i+4,i+4)-1) * limkl + kl
                                  ww_d(iw_loc, :) = ww_d(iw_loc, :) + cc_d * wrepp + cc * wrepp_d
                                  iminus = i - 1
                                  IF (iminus /= 0) THEN
                                     DO j = 1, iminus
                                        ij1 = inddd(i, j)
                                        cc   = ij_matrix%dd(ij1, i1-4, j1-4)
                                        cc_d = ij_matrix%dd_d(ij1, i1-4, j1-4, :)
                                        iw_loc = (indexb(i+4,j+4)-1) * limkl + kl
                                        ww_d(iw_loc, :) = ww_d(iw_loc, :) + cc_d * wrepp + cc * wrepp_d
                                     END DO
                                  END IF
                               END DO
                            END SELECT
                         END IF
                      END DO
                   END DO
                END DO
             END DO
             ! Store two electron integrals in the triangular format
             limij = (ii*(ii+1))/2
             limkl = (kk*(kk+1))/2
             CALL store_2el_2c_diag(limij, limkl, ww_dx=ww_d(1:istep,1), ww_dy=ww_d(1:istep,2), ww_dz=ww_d(1:istep,3),&
                  dw=dw, error=error)
             IF (invert) CALL invert_derivative(dw,42,error)
          END IF

          IF (debug_this_module) THEN
             ! Check derivatives
             CALL check_rotint_ana(sepi,sepj,rijv,dw=dw,shortrange=shortrange,taper=taper,error=error)
          END IF
       END IF
       CALL rotmat_release(ij_matrix, error)
    ENDIF

  END SUBROUTINE rotint_ana

! *****************************************************************************
!> \brief Calculates the derivative and the value of two-electron repulsion 
!>      integrals and the nuclear attraction integrals w.r.t. |r| 
!> \param on input rij   = interatomic distance
!>               sepi    = paramters of atom i
!>               sepj    = paramters of atom j
!>        short_range    = requests only the evaluation of the short range part
!>
!>     on output  ri     = array of two-electron repulsion integrals
!>     on output dri     = array of two-electron repulsion integrals derivatives
!>
!> \par History
!>      03.2008 created [tlaino]
!> \author Teodoro Laino [tlaino] - Zurich University
! *****************************************************************************
  RECURSIVE SUBROUTINE dterep_ana (sepi, sepj, rij, rep, rep_d, taper, shortrange,&
       do_longrange, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), INTENT(IN)                :: rij
    REAL(KIND=dp), DIMENSION(491), &
      INTENT(OUT)                            :: rep, rep_d
    TYPE(taper_type), POINTER                :: taper
    LOGICAL, INTENT(IN)                      :: shortrange
    LOGICAL, INTENT(IN), OPTIONAL            :: do_longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dterep_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ind
    LOGICAL                                  :: failure, longrange
    REAL(dp), DIMENSION(491)                 :: rep_d_lr, rep_lr
    REAL(KIND=dp)                            :: dft, ft
    REAL(KIND=dp), DIMENSION(22)             :: dri, ri

    failure = .FALSE.
    longrange = .FALSE.
    IF (PRESENT(do_longrange))  longrange = do_longrange
    rep   = 0.0_dp
    rep_d = 0.0_dp
    ri    = 0.0_dp
    dri   = 0.0_dp

    CALL dterep_sp_ana ( sepi, sepj, rij, ri, dri, taper, longrange, error=error)

    ! Store the value of the integrals in the proper order
    DO i = 1, 34
       ind      = map_22_to_491(i)
       rep(i)   =  ri(ind)
       rep_d(i) = dri(ind)
    END DO

    IF (sepi%dorb.OR.sepi%dorb) THEN
       ! compute the contribution from d-orbitals
       STOP "D-ORBITALS: dterep_ana"
    END IF

    ! Tapering function
    IF (.NOT.longrange) THEN
       ft = taper_eval(taper,rij,error)
       dft= dtaper_eval(taper,rij,error)
       rep_d(:) = dft*rep(:)+ft*rep_d(:)
       rep(:) = ft*rep(:)
    END IF

    IF ( shortrange ) THEN
       CALL dterep_ana ( sepi, sepj, rij, rep_lr, rep_d_lr, taper, shortrange=.FALSE.,&
            do_longrange=.TRUE., error=error)
       ! Correct Integral and Derivatives values
        rep  =   rep -  rep_lr
       rep_d = rep_d - rep_d_lr
    END IF

    ! Possibly debug 2el 2cent integrals and derivatives
    IF (debug_this_module) THEN
       CALL check_dterep_ana(sepi, sepj, rij, rep, rep_d, shortrange,taper=taper, error=error)
    END IF

  END SUBROUTINE dterep_ana

! *****************************************************************************
!> \brief Calculates the derivative and the value of two-electron repulsion 
!>      integrals and the nuclear attraction integrals w.r.t. |r| - sp core only
!> \param on input rij     = interatomic distance
!>               sepi    = paramters of atom i
!>               sepj    = paramters of atom j
!>     on output ri      = array of two-electron repulsion integrals
!>              dri      = array of derivatives of two-electron repulsion integrals
!>        The two-centre repulsion integrals (over local coordinates) are
!>        stored as follows (where p-sigma = O,  and p-pi = P and P* )
!>        (SS/SS)=1,      (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
!>        (SO/SO)=6,      (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
!>        (SS/OO)=11,     (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
!>        (OO/OO)=16,     (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
!>        (PP/P*P*)=21,   (P*P/P*P)=22.
!> \note
!>      Analytical version - Analytical evaluation of gradients
!>      Teodoro Laino - Zurich University 04.2007
!>      routine adapted from mopac7 (repp)
!>      vector version written by Ernest R. Davidson, Indiana University
!> \par History
!>      04.2007 created [tlaino]
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE dterep_sp_ana ( sepi, sepj, rij, ri, dri, taper, longrange, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(22), INTENT(INOUT)   :: ri, dri
    TYPE(taper_type), POINTER                :: taper
    LOGICAL, INTENT(IN)                      :: longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dterep_sp_ana', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: si, sj
    REAL(dp) :: ade, adi, adj, adq, aed, aee, aeq, ami, amj, aqd, aqe, aqi, &
      aqj, aqq, axx, da, db, ddi, ddj, ddxdx, ddxqxz, ddzdz, ddze, ddzqxx, &
      ddzqzz, dedz, dee, deqxx, deqzz, dqxxdz, dqxxe, dqxxqxx, dqxxqyy, &
      dqxxqzz, dqxzdx, dqxzqxz, dqzzdz, dqzze, dqzzqxx, dqzzqzz, drsq, dwww, &
      dxdx, dxqxz, dxxx, dyyy, dzdz, dze, dzqxx, dzqzz, dzzz, edz, ee, eqxx, &
      eqzz, qa, qb, qqi, qqj, qxxdz, qxxe, qxxqxx, qxxqyy, qxxqzz, qxzdx, &
      qxzqxz, qzzdz, qzze, qzzqxx, qzzqzz, r, rsq, www, xxx, yyy, zi, zj, zzz
    REAL(dp), DIMENSION(72)                  :: arg, darg, dsqr, sqr

    r=rij
    si = (sepi%natorb > 1)
    sj = (sepj%natorb > 1)
    zi = sepi%zeff
    zj = sepj%zeff
    ddi = sepi%cs(2)
    ddj = sepj%cs(2)
    qqi = sepi%cs(3)
    qqj = sepj%cs(3)

    !  hydrogen - hydrogen  (SS/SS)
    ami  = sepi%ko(1)
    amj  = sepj%ko(1)
    aee  = ami + amj
    aee  = aee * aee
    IF (longrange) aee = 0.0_dp
    rsq  = r*r
    drsq = 2.0_dp*r
    arg(1)  = rsq + aee
    darg(1) = drsq
    sqr(1)  = pxy(1)/SQRT(arg(1))
    dsqr(1) = -(0.5_dp/arg(1))*darg(1)*sqr(1)
    ee      = sqr(1)
    dee     = dsqr(1)
    ri(1)   = ee 
    dri(1)  = dee
    IF (si.OR.sj) THEN
       ! In case (i) is an heavy atom
       IF (si) THEN
          da=ddi
          qa=qqi
          IF (longrange) THEN
             ade = 0.0_dp
             aqe = 0.0_dp
          ELSE
             adi = sepi%ko(2)
             aqi = sepi%ko(3)
             ade = adi + amj
             ade = ade * ade
             aqe = aqi + amj
             aqe = aqe * aqe
          END IF
          xxx       = r+da
          arg(2)    = xxx*xxx + ade
          darg(2)   = 2.0_dp*xxx
          xxx       = r-da
          arg(3)    = xxx*xxx + ade
          darg(3)   = 2.0_dp*xxx
          xxx       = r+qa
          arg(4)    = xxx*xxx + aqe
          darg(4)   = 2.0_dp*xxx
          xxx       = r-qa
          arg(5)    = xxx*xxx + aqe
          darg(5)   = 2.0_dp*xxx
          arg(6)    = rsq + aqe
          darg(6)   = drsq
          arg(7)    = arg(6) + qa*qa
          darg(7)   = darg(6)
          sqr(2:7)  = pxy(2:7)/SQRT(arg(2:7))
          dsqr(2:7) = -(0.5_dp/arg(2:7))*darg(2:7)*sqr(2:7)
          dze      = sqr(2) + sqr(3)
          qzze     = sqr(4) + sqr(5) + sqr(6)
          qxxe     = sqr(7) + sqr(6)
          ri(2)    = dze
          ri(3)    = ee + qzze
          ri(4)    = ee + qxxe
          ddze     = dsqr(2) + dsqr(3)
          dqzze    = dsqr(4) + dsqr(5) + dsqr(6)
          dqxxe    = dsqr(7) + dsqr(6)
          dri(2)   = ddze
          dri(3)   = dee + dqzze
          dri(4)   = dee + dqxxe
       END IF

       ! In case (j) is an heavy atom
       IF (sj) THEN
          db=ddj
          qb=qqj
          IF (longrange) THEN
             aed = 0.0_dp
             aeq = 0.0_dp
          ELSE
             adj = sepj%ko(2)
             aqj = sepj%ko(3)
             aed = ami + adj
             aed = aed * aed
             aeq = ami + aqj
             aeq = aeq * aeq
          END IF
          xxx = r-db
          arg(8) = xxx*xxx + aed
          darg(8)= 2.0_dp*xxx
          xxx = r+db
          arg(9) = xxx*xxx + aed
          darg(9)= 2.0_dp*xxx
          xxx = r-qb
          arg(10) = xxx*xxx + aeq
          darg(10)= 2.0_dp*xxx
          xxx = r+qb
          arg(11) = xxx*xxx + aeq
          darg(11)= 2.0_dp*xxx
          arg(12) = rsq + aeq
          darg(12)= drsq
          arg(13) = arg(12) + qb*qb
          darg(13)= darg(12)
          sqr(8:13) = pxy(8:13)/SQRT(arg(8:13))
          dsqr(8:13) = -(0.5_dp/arg(8:13))*darg(8:13)*sqr(8:13)
          edz       = sqr(8) + sqr(9)
          eqzz      = sqr(10) + sqr(11) + sqr(12)
          eqxx      = sqr(12) + sqr(13)
          ri(5)     = edz
          ri(11)    = ee + eqzz
          ri(12)    = ee + eqxx
          dedz      = dsqr(8)  + dsqr(9)
          deqzz     = dsqr(10) + dsqr(11) + dsqr(12)
          deqxx     = dsqr(12) + dsqr(13)
          dri(5)    = dedz
          dri(11)   = dee + deqzz
          dri(12)   = dee + deqxx
       END IF

       ! In case both (i) and (j) are heavy atoms
       IF (si.AND.sj) THEN
          IF (longrange) THEN
             axx = 0.0_dp
             adq = 0.0_dp
             aqd = 0.0_dp
             aqq = 0.0_dp
          ELSE
             axx = adi + adj
             axx = axx * axx
             adq = adi + aqj
             adq = adq * adq
             aqd = aqi + adj
             aqd = aqd * aqd
             aqq = aqi + aqj
             aqq = aqq * aqq
          END IF
          xxx = da-db
          arg(14) = rsq + axx + xxx*xxx
          darg(14)= drsq
          xxx = da+db
          arg(15) = rsq + axx + xxx*xxx
          darg(15)= drsq
          xxx = r + da - db
          arg(16) = xxx*xxx + axx
          darg(16)= 2.0_dp* xxx
          xxx = r - da + db
          arg(17) = xxx*xxx + axx
          darg(17)= 2.0_dp* xxx
          xxx = r - da - db
          arg(18) = xxx*xxx + axx
          darg(18)= 2.0_dp* xxx
          xxx = r + da + db
          arg(19) = xxx*xxx + axx
          darg(19)= 2.0_dp* xxx
          xxx = r + da
          arg(20) = xxx*xxx + adq
          darg(20)= 2.0_dp* xxx
          arg(21) = arg(20) + qb*qb
          darg(21)= darg(20)
          xxx = r - da
          arg(22) = xxx*xxx + adq
          darg(22)= 2.0_dp* xxx
          arg(23) = arg(22) + qb*qb
          darg(23)= darg(22)
          xxx = r - db
          arg(24) = xxx*xxx + aqd
          darg(24)= 2.0_dp* xxx
          arg(25) = arg(24) + qa*qa
          darg(25)= darg(24)
          xxx = r + db
          arg(26) = xxx*xxx + aqd
          darg(26)= 2.0_dp* xxx
          arg(27) = arg(26) + qa*qa
          darg(27)= darg(26)
          xxx = r + da - qb
          arg(28) = xxx*xxx + adq
          darg(28)= 2.0_dp*xxx
          xxx = r - da - qb
          arg(29) = xxx*xxx + adq
          darg(29)= 2.0_dp* xxx
          xxx = r + da + qb
          arg(30) = xxx*xxx + adq
          darg(30)= 2.0_dp* xxx
          xxx = r - da + qb
          arg(31) = xxx*xxx + adq
          darg(31)= 2.0_dp* xxx
          xxx = r + qa - db
          arg(32) = xxx*xxx + aqd
          darg(32)= 2.0_dp* xxx
          xxx = r + qa + db
          arg(33) = xxx*xxx + aqd
          darg(33)= 2.0_dp* xxx
          xxx = r - qa - db
          arg(34) = xxx*xxx + aqd
          darg(34)= 2.0_dp* xxx
          xxx = r - qa + db
          arg(35) = xxx*xxx + aqd
          darg(35)= 2.0_dp* xxx
          arg(36) = rsq + aqq
          darg(36)= drsq
          xxx = qa - qb
          arg(37) = arg(36) + xxx*xxx
          darg(37)= darg(36)
          xxx = qa + qb
          arg(38) = arg(36) + xxx*xxx
          darg(38)= darg(36)
          arg(39) = arg(36) + qa*qa
          darg(39)= darg(36)
          arg(40) = arg(36) + qb*qb
          darg(40)= darg(36)
          arg(41) = arg(39) + qb*qb
          darg(41)= darg(39)
          xxx = r - qb
          arg(42) = xxx*xxx + aqq
          darg(42)= 2.0_dp*xxx
          arg(43) = arg(42) + qa*qa
          darg(43)= darg(42)
          xxx = r + qb
          arg(44) = xxx*xxx + aqq
          darg(44)= 2.0_dp*xxx
          arg(45) = arg(44) + qa*qa
          darg(45)= darg(44)
          xxx = r + qa
          arg(46) = xxx*xxx + aqq
          darg(46)= 2.0_dp*xxx
          arg(47) = arg(46) + qb*qb
          darg(47)= darg(46)
          xxx = r - qa
          arg(48) = xxx*xxx + aqq
          darg(48)= 2.0_dp*xxx
          arg(49) = arg(48) + qb*qb
          darg(49)= darg(48)
          xxx = r + qa - qb
          arg(50) = xxx*xxx + aqq
          darg(50)= 2.0_dp*xxx
          xxx = r + qa + qb
          arg(51) = xxx*xxx + aqq
          darg(51)= 2.0_dp*xxx
          xxx = r - qa - qb
          arg(52) = xxx*xxx + aqq
          darg(52)= 2.0_dp*xxx
          xxx = r - qa + qb
          arg(53) = xxx*xxx + aqq
          darg(53)= 2.0_dp*xxx
          qa=qa/2.0_dp
          qb=qb/2.0_dp
          xxx = da - qb
          dxxx= 0.0_dp
          xxx = xxx*xxx
          yyy = r - qb
          dyyy= 2.0_dp*yyy
          yyy = yyy*yyy
          zzz = da + qb
          dzzz= 0.0_dp
          zzz = zzz*zzz
          www = r + qb
          dwww= 2.0_dp*www
          www = www*www
          arg(54) =  xxx +  yyy + adq
          darg(54)= dxxx + dyyy 
          arg(55) =  xxx +  www + adq
          darg(55)= dxxx + dwww
          arg(56) =  zzz +  yyy + adq
          darg(56)= dzzz + dyyy
          arg(57) =  zzz +  www + adq
          darg(57)= dzzz + dwww
          xxx = qa - db
          dxxx= 0.0_dp
          xxx = xxx*xxx
          yyy = qa + db
          dyyy= 0.0_dp
          yyy = yyy*yyy
          zzz = r + qa
          dzzz= 2.0_dp*zzz
          zzz = zzz*zzz
          www = r - qa
          dwww= 2.0_dp*www
          www = www*www
          arg(58) =  zzz +  xxx + aqd
          darg(58)= dzzz + dxxx
          arg(59) =  www +  xxx + aqd
          darg(59)= dwww + dxxx
          arg(60) =  zzz +  yyy + aqd
          darg(60)= dzzz + dyyy
          arg(61) =  www +  yyy + aqd
          darg(61)= dwww + dyyy
          xxx = qa - qb
          xxx = xxx*xxx
          arg(62) =  arg(36) + 2._dp*xxx
          darg(62)= darg(36)
          yyy = qa + qb
          yyy = yyy*yyy
          arg(63) =  arg(36) + 2._dp*yyy
          darg(63)= darg(36)
          arg(64) =  arg(36) + 2._dp*(qa*qa+qb*qb)
          darg(64)= darg(36)
          zzz = r + qa - qb
          dzzz= 2.0_dp*zzz
          zzz = zzz*zzz
          arg(65) = zzz + xxx + aqq
          darg(65)= dzzz
          arg(66) = zzz + yyy + aqq
          darg(66)= dzzz
          zzz = r + qa + qb
          dzzz= 2.0_dp*zzz
          zzz = zzz*zzz
          arg(67) = zzz + xxx + aqq
          darg(67)= dzzz
          arg(68) = zzz + yyy + aqq
          darg(68)= dzzz
          zzz = r - qa - qb
          dzzz= 2.0_dp*zzz
          zzz = zzz*zzz
          arg(69) = zzz + xxx + aqq
          darg(69)= dzzz
          arg(70) = zzz + yyy + aqq
          darg(70)= dzzz
          zzz = r - qa + qb
          dzzz= 2.0_dp*zzz
          zzz = zzz*zzz
          arg(71) = zzz + xxx + aqq
          darg(71)= dzzz
          arg(72) = zzz + yyy + aqq
          darg(72)= dzzz
          sqr(14:72)  = pxy(14:72)/SQRT(arg(14:72))
          dsqr(14:72) = -(0.5_dp/arg(14:72))*darg(14:72)*sqr(14:72)
          dxdx   = sqr(14) + sqr(15)
          dzdz   = sqr(16) + sqr(17) + sqr(18) + sqr(19)
          dzqxx  = sqr(20) + sqr(21) + sqr(22) + sqr(23)
          qxxdz  = sqr(24) + sqr(25) + sqr(26) + sqr(27)
          dzqzz  = sqr(28) + sqr(29) + sqr(30) + sqr(31) + sqr(22) + sqr(20)
          qzzdz  = sqr(32) + sqr(33) + sqr(34) + sqr(35) + sqr(24) + sqr(26)
          qxxqxx = sqr(37) + sqr(38) + sqr(39) + sqr(40) + sqr(36)
          qxxqyy = sqr(41) + sqr(39) + sqr(40) + sqr(36)
          qxxqzz = sqr(43) + sqr(45) + sqr(42) + sqr(44) + sqr(39) + sqr(36)
          qzzqxx = sqr(47) + sqr(49) + sqr(46) + sqr(48) + sqr(40) + sqr(36)
          qzzqzz = sqr(50) + sqr(51) + sqr(52) + sqr(53) + sqr(48) + sqr(46) + sqr(42) + sqr(44) + sqr(36)
          dxqxz  = sqr(54) + sqr(55) + sqr(56) + sqr(57)
          qxzdx  = sqr(58) + sqr(59) + sqr(60) + sqr(61)
          qxzqxz = sqr(65) + sqr(67) + sqr(69) + sqr(71) + sqr(66) + sqr(68) + sqr(70) + sqr(72)
          ddxdx   = dsqr(14) + dsqr(15)
          ddzdz   = dsqr(16) + dsqr(17) + dsqr(18) + dsqr(19)
          ddzqxx  = dsqr(20) + dsqr(21) + dsqr(22) + dsqr(23)
          dqxxdz  = dsqr(24) + dsqr(25) + dsqr(26) + dsqr(27)
          ddzqzz  = dsqr(28) + dsqr(29) + dsqr(30) + dsqr(31) + dsqr(22) + dsqr(20)
          dqzzdz  = dsqr(32) + dsqr(33) + dsqr(34) + dsqr(35) + dsqr(24) + dsqr(26)
          dqxxqxx = dsqr(37) + dsqr(38) + dsqr(39) + dsqr(40) + dsqr(36)
          dqxxqyy = dsqr(41) + dsqr(39) + dsqr(40) + dsqr(36)
          dqxxqzz = dsqr(43) + dsqr(45) + dsqr(42) + dsqr(44) + dsqr(39) + dsqr(36)
          dqzzqxx = dsqr(47) + dsqr(49) + dsqr(46) + dsqr(48) + dsqr(40) + dsqr(36)
          dqzzqzz = dsqr(50) + dsqr(51) + dsqr(52) + dsqr(53) + dsqr(48) + dsqr(46) + dsqr(42) + dsqr(44) + dsqr(36)
          ddxqxz  = dsqr(54) + dsqr(55) + dsqr(56) + dsqr(57)
          dqxzdx  = dsqr(58) + dsqr(59) + dsqr(60) + dsqr(61)
          dqxzqxz = dsqr(65) + dsqr(67) + dsqr(69) + dsqr(71) + dsqr(66) + dsqr(68) + dsqr(70) + dsqr(72)
          ri(6)  = dzdz
          ri(7)  = dxdx
          ri(8)  =  edz - qzzdz
          ri(9)  =  edz - qxxdz
          ri(10) =      - qxzdx
          ri(13) =  dze - dzqzz
          ri(14) =  dze - dzqxx
          ri(15) =      - dxqxz
          ri(16) = ee + eqzz + qzze + qzzqzz
          ri(17) = ee + eqzz + qxxe + qxxqzz
          ri(18) = ee + eqxx + qzze + qzzqxx
          ri(19) = ee + eqxx + qxxe + qxxqxx
          ri(20) =                    qxzqxz
          ri(21) = ee + eqxx + qxxe + qxxqyy
          ri(22) = 0.5_dp * (qxxqxx - qxxqyy)
          dri(6)  = ddzdz
          dri(7)  = ddxdx
          dri(8)  =  dedz - dqzzdz
          dri(9)  =  dedz - dqxxdz
          dri(10) =       - dqxzdx
          dri(13) =  ddze - ddzqzz
          dri(14) =  ddze - ddzqxx
          dri(15) =       - ddxqxz
          dri(16) = dee + deqzz + dqzze + dqzzqzz
          dri(17) = dee + deqzz + dqxxe + dqxxqzz
          dri(18) = dee + deqxx + dqzze + dqzzqxx
          dri(19) = dee + deqxx + dqxxe + dqxxqxx
          dri(20) =                       dqxzqxz
          dri(21) = dee + deqxx + dqxxe + dqxxqyy
          dri(22) = 0.5_dp  * ( dqxxqxx - dqxxqyy )
       END IF
    END IF

  END SUBROUTINE dterep_sp_ana

END MODULE semi_empirical_int_ana
