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
!>      Teodoro Laino (05.2008) [tlaino] - University of Zurich : 
! *****************************************************************************
MODULE semi_empirical_int_ana

  USE f77_blas
  USE input_constants,                 ONLY: do_method_am1,&
                                             do_method_pchg,&
                                             do_method_pdg,&
                                             do_method_pm3,&
                                             do_method_pm6,&
                                             do_method_undef
  USE kinds,                           ONLY: dp
  USE physcon,                         ONLY: angstrom,&
                                             evolt
  USE semi_empirical_int_arrays,       ONLY: &
       fac_x_to_z, ijkl_ind, ijkl_sym, inddd, inddp, indexa, indexb, indpp, &
       int2c_type, l_index, map_x_to_z, rij_threshold
  USE semi_empirical_int_num,          ONLY: nucint_d_num,&
                                             nucint_sp_num,&
                                             terep_d_num,&
                                             terep_sp_num
  USE semi_empirical_int_utils,        ONLY: d_ijkl_d,&
                                             d_ijkl_sp,&
                                             dtaper_eval,&
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
!>      Teodoro Laino [tlaino] - University of Zurich 04.2008 : removed the core-core part
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
       CALL rotmat (sepi, sepj, rijv, rij, ij_matrix, do_derivatives=lgrad, do_invert=invert,&
            error=error)

       IF (lgrad) THEN
          drij(:) = rijv/rij
          ! Possibly Invert Frame
          IF (invert) THEN
             xtmp=drij(3)
             drij(3)=drij(1)
             drij(1)=xtmp
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
       IF (invert.AND.l_e1b)  CALL invert_integral  (sepi, sepi, int1el=e1b, error=error)
       IF (invert.AND.l_e2a)  CALL invert_integral  (sepj, sepj, int1el=e2a, error=error)

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
       IF (invert.AND.l_de1b) CALL invert_derivative(sepi, sepi, dint1el=de1b, error=error)
       IF (invert.AND.l_de2a) CALL invert_derivative(sepj, sepj, dint1el=de2a, error=error)
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
       ! Compute Integrals in diatomic frame
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
               poi = sepi%ko(9)
             END IF
             IF ( ABS(sepj%rho) > EPSILON(0.0_dp) ) THEN
               poj = sepj%rho
             ELSE
               poj = sepj%ko(9)
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
                enuc_loc=scale*(1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))+9.28_dp*EXP(-5.98_dp*rija))
                IF (l_denuc) THEN
                   denuc_loc = dscale * (1._dp+2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))+9.28_dp*EXP(-5.98_dp*rija))&
                        -scale*2._dp*xab*EXP(-aab*(rija+0.0003_dp*rija**6))*(aab+6.0_dp*0.0003_dp*rija**5)*drija &
                        -scale*9.28_dp*EXP(-5.98_dp*rija)*5.98_dp*drija
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
!>      05.2008 New driver for integral invertion (supports d-orbitals)
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE invert_integral(sepi, sepj, int1el, int2el, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), DIMENSION(:), INTENT(INOUT), &
      OPTIONAL                               :: int1el, int2el
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'invert_integral', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: fdim, gind, gknd, i, imap, &
                                                ind, j, jmap, jnd, k, kmap, &
                                                knd, l, lmap, lnd, ndim, &
                                                sdim, tdim, tind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: ifac, jfac, kfac, lfac
    REAL(KIND=dp), DIMENSION(2025)           :: tmp2el
    REAL(KIND=dp), DIMENSION(45)             :: tmp1el

    failure = .FALSE.
    ! One-electron integral
    IF (PRESENT(int1el)) THEN
       fdim   = (sepi%natorb+1)*sepi%natorb/2
       ndim   = 0
       tmp1el(1:fdim) = 0.0_dp
       DO i = 1, sepi%natorb             
       DO j = 1, i
          ndim = ndim + 1
          
          ! Get the integral in the original frame (along z)
          DO ind = 1, 2
             imap = map_x_to_z(ind,i)
             IF (imap==0) CYCLE
             ifac = fac_x_to_z(ind,i)
             DO jnd = 1, 2
                jmap = map_x_to_z(jnd,j)
                IF (jmap==0) CYCLE
                jfac = fac_x_to_z(jnd,j)
                gind = indexb(imap,jmap)

                tmp1el(ndim) = tmp1el(ndim) + ifac*jfac*int1el(gind)
             END DO
          END DO
       END DO
       END DO
       int1el(1:fdim) = tmp1el(1:fdim)
    END IF

    ! Two electron integrals
    IF (PRESENT(int2el)) THEN
       sdim   = (sepi%natorb+1)*sepi%natorb/2
       tdim   = (sepj%natorb+1)*sepj%natorb/2
       fdim   = sdim*tdim
       ndim   = 0
       tmp2el(1:fdim) = 0.0_dp
       DO i = 1, sepi%natorb             
       DO j = 1, i
          DO k = 1, sepj%natorb             
          DO l = 1, k
             ndim = ndim + 1

             ! Get the integral in the original frame (along z)
             DO ind = 1, 2
                imap = map_x_to_z(ind,i)
                IF (imap==0) CYCLE
                ifac = fac_x_to_z(ind,i)
                DO jnd = 1, 2
                   jmap = map_x_to_z(jnd,j)
                   IF (jmap==0) CYCLE
                   jfac = fac_x_to_z(jnd,j)
                   gind = indexb(imap,jmap)
                   
             ! Get the integral in the original frame (along z)
             DO knd = 1, 2
                kmap = map_x_to_z(knd,k)
                IF (kmap==0) CYCLE
                kfac = fac_x_to_z(knd,k)
                DO lnd = 1, 2
                   lmap = map_x_to_z(lnd,l)
                   IF (lmap==0) CYCLE
                   lfac = fac_x_to_z(lnd,l)
                   gknd = indexb(kmap,lmap)
                   
                   tind = (gind-1)*tdim+gknd
                   tmp2el(ndim) = tmp2el(ndim) + ifac*jfac*lfac*kfac*int2el(tind)
                END DO
             END DO
                       
                END DO
             END DO
                   
          END DO
          END DO
       END DO
       END DO
       int2el(1:fdim) = tmp2el(1:fdim)
    END IF
  END SUBROUTINE invert_integral
  
! *****************************************************************************
!> \brief Exploits inversion symmetry to avoid divergence
!> \par History
!>      04.2007 created [tlaino]
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE invert_derivative(sepi, sepj, dint1el, dint2el, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: dint1el, dint2el
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'invert_derivative', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, m
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: tmp

    failure = .FALSE.
    ! Integral part
    DO i=1,3
       IF (PRESENT(dint1el)) THEN
          CALL invert_integral(sepi, sepj, int1el=dint1el(:,i), error=error)
       END IF
       IF (PRESENT(dint2el)) THEN
          CALL invert_integral(sepi, sepj, int2el=dint2el(:,i), error=error)
       END IF
    END DO

    ! Derivatives part
    IF (PRESENT(dint1el)) THEN
       DO m=1,SIZE(dint1el,1)
          tmp=dint1el(m,3)
          dint1el(m,3)=dint1el(m,1)
          dint1el(m,1)=tmp
       END DO
    END IF
    IF (PRESENT(dint2el)) THEN
       DO m=1,SIZE(dint2el,1)
          tmp=dint2el(m,3)
          dint2el(m,3)=dint2el(m,1)
          dint2el(m,1)=tmp
       END DO
    END IF
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
    
    ! Value of the integrals for sp shell
    CALL  nucint_sp_num(sepi, sepj, rij, ssss,  core, itype, longrange, error=error)
    ! Integrals derivatives for sp shell
    CALL dnucint_sp_ana(sepi, sepj, rij,dssss, dcore, itype, longrange, error=error)

    IF ((sepi%dorb.OR.sepj%dorb).AND.l_core) THEN
       ! Compute the contribution from d-orbitals
       CALL  nucint_d_num(sepi, sepj, rij, core, itype, longrange, error=error)
       ! Compute the contribution from d-orbitals
       CALL dnucint_d_ana(sepi, sepj, rij,dcore, itype, longrange, error=error)
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
!> \brief Calculates the nuclear attraction integrals and derivatives for sp basis
!> \param on input rij     = interatomic distance
!>                sepi    = paramters of atom i
!>                sepj    = paramters of atom j
!>                itype   = type of semi_empirical model
!>                          extension to the original routine to compute qm/mm integrals
!>      on output core    = derivative of 4 X 2 array of electron-core attraction integrals
!>                ssss    = derivative of (ssss) integral
!>                          derivatives are intended w.r.t. rij
!>         The storage of the nuclear attraction integrals  core(kl/ij) iS
!>         (SS/)=1,   (SP/)=2,   (PP/)=3,   (P+P+/)=4,   (SD/)=5,
!>         (DP/)=6,   (DD/)=7,   (D+P+)=8,  (D+D+/)=9,   (D#D#)=10
!>
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
!>      05.2008 Teodoro Laino [tlaino] - University of Zurich: major rewriting
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE dnucint_sp_ana ( sepi, sepj, rij, dssss, dcore, itype, longrange, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), INTENT(INOUT), OPTIONAL        :: dssss
    REAL(dp), DIMENSION(10, 2), &
      INTENT(INOUT), OPTIONAL                :: dcore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dnucint_sp_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ij, kl
    LOGICAL                                  :: failure, l_core, l_ssss, si, &
                                                sj

    failure   = .FALSE.
    l_core  = PRESENT(dcore)
    l_ssss  = PRESENT(dssss)

    IF (.NOT.(l_core.OR.l_ssss)) RETURN
    si  = (sepi%natorb > 1)
    sj  = (sepj%natorb > 1)

    ij = indexa(1, 1)
    IF (l_ssss) THEN
       ! Store the value for the derivative of <S  S  | S  S  > (Used for computing the core-core interactions)
       dssss = d_ijkl_sp (sepi, sepj, ij, ij, 0, 0, 0, 0, -1, rij, longrange, itype, error)
    END IF

    IF (l_core) THEN
       !     <S  S  | S  S  >
       kl = indexa(1, 1)
       dcore(1,1)    = d_ijkl_sp (sepi, sepj, kl, ij, 0, 0, 0, 0, 1, rij, longrange, itype, error) * sepj%zeff
       IF (si) THEN
          !  <S  P  | S  S  >
          kl = indexa(2, 1)
          dcore(2,1) = d_ijkl_sp (sepi, sepj, kl, ij, 0, 1, 0, 0, 1, rij, longrange, itype, error) * sepj%zeff
          !  <P  P  | S  S  >
          kl = indexa(2, 2)
          dcore(3,1) = d_ijkl_sp (sepi, sepj, kl, ij, 1, 1, 0, 0, 1, rij, longrange, itype, error) * sepj%zeff
          !  <P+ P+ | S  S  >
          kl = indexa(3, 3)
          dcore(4,1) = d_ijkl_sp (sepi, sepj, kl, ij, 1, 1, 0, 0, 1, rij, longrange, itype, error) * sepj%zeff
       END IF

       !     <S  S  | S  S  >
       kl = indexa(1, 1)
       dcore(1,2)    = d_ijkl_sp (sepi, sepj, ij, kl, 0, 0, 0, 0, 2, rij, longrange, itype, error) * sepi%zeff
       IF (sj) THEN
          !  <S  S  | S  P  >
          kl = indexa(2, 1)
          dcore(2,2) = d_ijkl_sp (sepi, sepj, ij, kl, 0, 0, 0, 1, 2, rij, longrange, itype, error) * sepi%zeff
          !  <S  S  | P  P  >
          kl = indexa(2, 2)
          dcore(3,2) = d_ijkl_sp (sepi, sepj, ij, kl, 0, 0, 1, 1, 2, rij, longrange, itype, error) * sepi%zeff
          !  <S  S  | P+ P+ >
          kl = indexa(3, 3)
          dcore(4,2) = d_ijkl_sp (sepi, sepj, ij, kl, 0, 0, 1, 1, 2, rij, longrange, itype, error) * sepi%zeff
       END IF
    END IF

  END SUBROUTINE dnucint_sp_ana

! *****************************************************************************
!> \brief Calculates the analytical derivative of the nuclear attraction 
!>        integrals involving d orbitals
!> \param on input rij   = interatomic distance
!>                sepi   = paramters of atom i
!>                sepj   = paramters of atom j
!>                itype  = type of semi_empirical model
!>                         extension to the original routine to compute qm/mm integrals
!>        short_range    = requests only the evaluation of the short range part
!>          longrange    = requests only the evaluation of the long range part
!>      on output core   = 4 X 2 array of electron-core attraction integrals
!>
!>         The storage of the nuclear attraction integrals  core(kl/ij) iS
!>         (SS/)=1,   (SP/)=2,   (PP/)=3,   (P+P+/)=4,   (SD/)=5,
!>         (DP/)=6,   (DD/)=7,   (D+P+)=8,  (D+D+/)=9,   (D#D#)=10
!>
!>         where ij=1 if the orbitals centred on atom i,  =2 if on atom j.
!> \author
!>      Teodoro Laino (05.2008) [tlaino] - University of Zurich: created
! *****************************************************************************
  SUBROUTINE dnucint_d_ana( sepi, sepj, rij, dcore, itype, longrange, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(10, 2), &
      INTENT(INOUT)                          :: dcore
    INTEGER, INTENT(IN)                      :: itype
    LOGICAL, INTENT(IN)                      :: longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dnucint_d_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij, kl
    LOGICAL                                  :: failure

    failure = .FALSE.
    ! Check if d-orbitals are present
    IF (sepi%dorb.OR.sepj%dorb) THEN
       DO i = 5, 10
          dcore(i, 1) = 0.0_dp
          dcore(i, 2) = 0.0_dp
       END DO
       !
       ij = indexa(1, 1)
       IF (sepj%dorb) THEN
          !  <S S | D S>
          kl = indexa(5, 1)
          dcore(5, 2) = d_ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 0, 1, rij, longrange, itype, error) * sepi%zeff
          !  <S S | D P >
          kl = indexa(5, 2)
          dcore(6, 2) = d_ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 1, 1, rij, longrange, itype, error) * sepi%zeff
          !  <S S | D D >
          kl = indexa(5, 5)
          dcore(7, 2) = d_ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, rij, longrange, itype, error) * sepi%zeff
          !  <S S | D+P+>
          kl = indexa(6, 3)
          dcore(8, 2) = d_ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 1, 1, rij, longrange, itype, error) * sepi%zeff
          !  <S S | D+D+>
          kl = indexa(6, 6)
          dcore(9, 2) = d_ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, rij, longrange, itype, error) * sepi%zeff
          !  <S S | D#D#>
          kl = indexa(8, 8)
          dcore(10, 2)= d_ijkl_d (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, rij, longrange, itype, error) * sepi%zeff
       END IF
       IF (sepi%dorb) THEN
          !  <D S | S S>
          kl = indexa(5, 1)
          dcore(5, 1) = d_ijkl_d (sepi, sepj, kl, ij, 2, 0, 0, 0, 2, rij, longrange, itype, error) * sepj%zeff
          !  <D P | S S >
          kl = indexa(5, 2)
          dcore(6, 1) = d_ijkl_d (sepi, sepj, kl, ij, 2, 1, 0, 0, 2, rij, longrange, itype, error) * sepj%zeff
          !  <D D | S S >
          kl = indexa(5, 5)
          dcore(7, 1) = d_ijkl_d (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, rij, longrange, itype, error) * sepj%zeff
          !  <D+P+| S S >
          kl = indexa(6, 3)
          dcore(8, 1) = d_ijkl_d (sepi, sepj, kl, ij, 2, 1, 0, 0, 2, rij, longrange, itype, error) * sepj%zeff
          !  <D+D+| S S >
          kl = indexa(6, 6)
          dcore(9, 1) = d_ijkl_d (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, rij, longrange, itype, error) * sepj%zeff
          !  <D#D#| S S >
          kl = indexa(8, 8)
          dcore(10, 1)= d_ijkl_d (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, rij, longrange, itype, error) * sepj%zeff
       END IF

    END IF
  END SUBROUTINE dnucint_d_ana

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
       CALL rotmat (sepi, sepj, rijv, rij, ij_matrix, do_derivatives=lgrad, do_invert=invert, error=error)          

       ! Compute integrals in diatomic frame as well their derivatives (if requested)
       CALL dterep_ana(sepi,sepj,rij,rep,rep_d,taper,shortrange,error=error)

       IF (lgrad) THEN
          drij(:) = rijv/rij
          ! Possibly Invert Frame
          IF (invert) THEN
             xtmp=drij(3)
             drij(3)=drij(1)
             drij(1)=xtmp
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
             IF (invert) CALL invert_integral(sepi, sepj, int2el=w, error=error)
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
             IF (invert) CALL invert_derivative(sepi, sepj, dint2el=dw, error=error)
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

    LOGICAL                                  :: failure, longrange
    REAL(dp), DIMENSION(491)                 :: rep_d_lr, rep_lr
    REAL(KIND=dp)                            :: dft, ft

    failure = .FALSE.
    longrange = .FALSE.
    IF (PRESENT(do_longrange))  longrange = do_longrange
    rep   = 0.0_dp
    rep_d = 0.0_dp

    ! Integral Values for sp shells only
    CALL  terep_sp_num(sepi, sepj, rij,   rep(1:34), longrange, error=error)
    ! Integral Derivatives
    CALL dterep_sp_ana(sepi, sepj, rij, rep_d(1:34), longrange, error=error)

    IF (sepi%dorb.OR.sepj%dorb) THEN
       ! Compute the contribution from d-orbitals
       CALL  terep_d_num (sepi, sepj, rij, rep,   longrange, error=error )
       ! Compute the derivatives from d-orbitals
       CALL dterep_d_ana (sepi, sepj, rij, rep_d, longrange, error=error )
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
!>      integrals and the nuclear attraction integrals w.r.t. |r| - sp shells only
!> \param on input rij     = interatomic distance
!>               sepi    = paramters of atom i
!>               sepj    = paramters of atom j
!>     on output ri      = array of two-electron repulsion integrals
!>              dri      = array of derivatives of two-electron repulsion integrals
!> \note
!>      Analytical version - Analytical evaluation of gradients
!>      Teodoro Laino - Zurich University 04.2007
!>      routine adapted from mopac7 (repp)
!>      vector version written by Ernest R. Davidson, Indiana University
!> \par History
!>      04.2007 created [tlaino]
!>      Teodoro Laino (03.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integrals
!>      05.2008 Teodoro Laino [tlaino] - University of Zurich: major rewriting
!> \author Teodoro Laino - Zurich University
! *****************************************************************************
  SUBROUTINE dterep_sp_ana ( sepi, sepj, rij, drep, longrange, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(34), INTENT(INOUT)   :: drep
    LOGICAL, INTENT(IN)                      :: longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dterep_sp_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij, j, k, kl, l, lasti, &
                                                lastj, li, lj, lk, ll, nold, &
                                                numb
    LOGICAL                                  :: failure

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
                      drep(numb) =  drep( nold)
                   ELSE IF (nold < 0) THEN
                      drep(numb) = -drep(-nold)
                   ELSE IF (nold == 0) THEN
                      drep(numb) = d_ijkl_sp (sepi, sepj, ij, kl, li, lj, lk, ll, 0, rij,&
                                              longrange, do_method_undef, error)
                   END IF
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE dterep_sp_ana

! *****************************************************************************
!> \brief Calculates the derivatives of the two-electron repulsion integrals - d shell only
!> \param on input rij   = interatomic distance
!>               sepi    = paramters of atom i
!>               sepj    = paramters of atom j
!>           longrange   = requests only the evaluation of the long range part
!>
!>     on output rep     = array of two-electron repulsion integrals
!>
!> \par History
!>      Teodoro Laino (05.2008) [tlaino] - University of Zurich : new driver
!>                 for computing integral derivatives for d-orbitals
! *****************************************************************************
  SUBROUTINE dterep_d_ana (sepi, sepj, rij, drep, longrange, error )
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    REAL(dp), INTENT(IN)                     :: rij
    REAL(dp), DIMENSION(491), INTENT(INOUT)  :: drep
    LOGICAL, INTENT(IN)                      :: longrange
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dterep_d_ana', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij, j, k, kl, l, lasti, &
                                                lastj, li, lj, lk, ll, nold, &
                                                numb
    LOGICAL                                  :: failure

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
                      drep(numb) =  drep( nold)
                   ELSE IF (nold <-34) THEN
                      drep(numb) = -drep(-nold)
                   ELSE IF (nold == 0) THEN
                      drep(numb) = d_ijkl_d (sepi, sepj, ij, kl, li, lj, lk, ll, 0, rij,&
                                             longrange, do_method_undef, error)
                   END IF
                END IF
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE dterep_d_ana

END MODULE semi_empirical_int_ana
