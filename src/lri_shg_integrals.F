!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates integral matrices for LRIGPW method using solid harmonic
!>        Gaussians (SHG)
!>        lri : local resolution of the identity
!> \par History
!>      created [06.2015]
!> \author Dorothea Golze  
! *****************************************************************************
MODULE lri_shg_integrals
  USE ai_shg_overlap,                  ONLY: get_W_matrix,&
                                             get_real_scaled_solid_harmonic,&
                                             overlap_shg_ab,&
                                             overlap_shg_aba,&
                                             overlap_shg_abb,&
                                             s_overlap_ab,&
                                             s_overlap_abx
  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE kinds,                           ONLY: dp
  USE lri_environment_types,           ONLY: lri_environment_type
  USE mathconstants,                   ONLY: fac,&
                                             ifac,&
                                             pi
  USE orbital_pointers,                ONLY: nsoset
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_shg_integrals'

  PUBLIC :: precalc_angular_shg_part, lri_int_ab_shg, lri_int_aba_shg,&
            lri_int_abb_shg

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief precalculates the angular part of the SHG integrals for the matrices
!>        (fa,fb), (a,b), (a,b,fa) and (b,fb,a)
!> \param oba orbital basis on a
!> \param obb orbital basis on b
!> \param fba aux basis on a
!> \param fbb aux basis on b
!> \param rab distance vector between a and b
!> \param Waux_mat matrix storing angular-dependent part
!> \param error ...
! *****************************************************************************
  SUBROUTINE precalc_angular_shg_part(oba,obb,fba,fbb,rab,Waux_mat,error)

    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fba, fbb
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: Waux_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'precalc_angular_shg_part', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, isize, j, k, la_max, laa_max, lb_max, lbb_max, &
      lca_max, lcb_max, li_max, lj_max, lmax, mdim(4), size_int(4,2), stat, &
      temp
    INTEGER, DIMENSION(:), POINTER           :: li_max_all
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rab2
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: Rc, Rs

    CALL timeset(routineN,handle)
    NULLIFY(li_max_all) 

    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)

    !*** 1 Waux_mat of size (li_max,lj_max) for elements
    !                    i        j
    !    [aab]    --> (laa_max, lb_max)
    !    [bba]    --> (lbb_max, la_max) --> use for [abb]
    !    [ab] ri  --> (lca_max, lcb_max)
    !    [ab] orb --> (la_max , lb_max)

    la_max   = MAXVAL(oba%lmax)
    lb_max   = MAXVAL(obb%lmax)
    lca_max  = MAXVAL(fba%lmax)
    lcb_max  = MAXVAL(fbb%lmax)

    laa_max = la_max + lca_max
    lbb_max = lb_max + lcb_max
    li_max = MAX(laa_max,lbb_max)
    lj_max = MAX(la_max, lb_max, lcb_max)
    lmax = li_max
 
    ALLOCATE(li_max_all(0:lj_max),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(Rc(0:lmax,-2*lmax:2*lmax),Rs(0:lmax,-2*lmax:2*lmax),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    Rc=0._dp
    Rs=0._dp
    mdim(1) = nsoset(li_max) + 1
    mdim(2) = nsoset(lj_max) + 1
    mdim(3) = li_max+lj_max  + 1
    mdim(4) = 4
    ALLOCATE(Waux_mat(mdim(1),mdim(2),mdim(3),mdim(4)),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    !Waux_mat = 0._dp .. takes time

    !*** Waux_mat (li_max,lj_max) contains elements not needed,
    !*** make indixing so that only required ones are computed
    !*** li_max_all(j) --> li_max dependent on j 
    size_int(1,:) = (/laa_max,lb_max/)
    size_int(2,:) = (/lbb_max,la_max/)
    size_int(3,:) = (/lca_max,lcb_max/)
    size_int(4,:) = (/la_max,lb_max/)

    li_max_all(:) = 0
    DO isize = 1,4
      i = size_int(isize,1)
      j = size_int(isize,2)
      k = li_max_all(j)
      IF(k < i) li_max_all(j) = i
    ENDDO
    temp =  li_max_all(lj_max)
    DO j = lj_max,0,-1
      IF(li_max_all(j) < temp) THEN
        li_max_all(j) = temp
      ELSE
        temp = li_max_all(j)
      ENDIF
    ENDDO

    !*** -rab, since Eq. in Ref. use Ra-Rb, not Rb-Ra
    CALL get_real_scaled_solid_harmonic(Rc,Rs,lmax,-rab,rab2)
    CALL get_W_matrix(li_max_all,lj_max,lmax,Rc,Rs,Waux_mat)

    DEALLOCATE(Rc,Rs,li_max_all,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE precalc_angular_shg_part

! *****************************************************************************
!> \brief calculate integrals (a,b,fa)
!> \param lri_env ...
!> \param abaint integral (a,b,fa)
!> \param rab distance vector between A and B
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fba auxiliary basis set at center A
!> \param scona_shg contraction matrix for orb bas on A
!> \param sconb_shg contraction matrix for orb bas on B
!> \param sconca_shg contraction matrix for aux bas on A
!> \param Waux_mat angular-depedent part
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_int_aba_shg(lri_env,abaint,rab,oba,obb,fba,scona_shg,sconb_shg,&
                             sconca_shg,Waux_mat,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: abaint
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fba
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: scona_shg, sconb_shg, &
                                                sconca_shg
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: Waux_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_aba_shg', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ids, ids_start, iil, il, ipgfa, iset, ishella, jpgfb, &
      jset, jshellb, kpgfa, kset, kshella, la_max_set, laa, laa_max_set, lai, &
      lak, lb_max_set, lbj, lca_max_set, nds, nl, nl_set, npgfa_set, &
      npgfb_set, npgfca_set, nseta, nsetb, nsetca, nshella_set, nshellb_set, &
      nshellca_set, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, lb_max, lca_max, &
                                                npgfa, npgfb, npgfca, &
                                                nshella, nshellb, nshellca
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfca, la, lb, lca
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: a, aif, bjf, const, consta, &
                                                constca, dab, prefac, rab2, &
                                                shg_fac(5)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b, &
                                                set_radius_ca
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta, zetb, zetca
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: swork, swork_cont

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, lb_max, lca_max, npgfa, npgfb, npgfca) 
    NULLIFY(first_sgfa, first_sgfb, first_sgfca, set_radius_a, set_radius_b,&
            set_radius_ca, zeta, zetb, zetca) 

    ! basis ikind
    first_sgfa    =>  oba%first_sgf
    la_max        =>  oba%lmax
    la            =>  oba%l
    npgfa         =>  oba%npgf
    nshella       =>  oba%nshell
    nseta         =   oba%nset
    set_radius_a  =>  oba%set_radius
    zeta          =>  oba%zet
    ! basis jkind
    first_sgfb    =>  obb%first_sgf
    lb_max        =>  obb%lmax
    lb            =>  obb%l
    npgfb         =>  obb%npgf
    nshellb       =>  obb%nshell
    nsetb         =   obb%nset
    set_radius_b  =>  obb%set_radius
    zetb          =>  obb%zet

    ! basis RI A 
    first_sgfca    =>  fba%first_sgf
    lca_max        =>  fba%lmax
    lca            =>  fba%l
    npgfca         =>  fba%npgf
    nshellca       =>  fba%nshell
    nsetca         =   fba%nset
    set_radius_ca  =>  fba%set_radius
    zetca          =>  fba%zet

    dab = SQRT( SUM(rab**2) )
    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)

    la_max_set   = MAXVAL(la_max)
    lb_max_set   = MAXVAL(lb_max)
    lca_max_set  = MAXVAL(lca_max)
    npgfa_set    = MAXVAL(npgfa)
    npgfb_set    = MAXVAL(npgfb)
    npgfca_set   = MAXVAL(npgfca)
    nshella_set  = MAXVAL(nshella)
    nshellb_set  = MAXVAL(nshellb)
    nshellca_set = MAXVAL(nshellca)

    laa_max_set = la_max_set+lca_max_set

    ! allocate some work storage....
    nds = laa_max_set + lb_max_set + 1
    nl_set = INT((laa_max_set)/2)
    ALLOCATE(swork(npgfa_set,npgfb_set,npgfca_set,0:nl_set,nds),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(swork_cont(nshella_set,nshellb_set,nshellca_set,0:nl_set,nds),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO iset=1,nseta

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          DO kset=1,nsetca

             IF (set_radius_b(jset) + set_radius_ca(kset) < dab) CYCLE

             CALL s_overlap_abx(la_max(iset),npgfa(iset),zeta(:,iset),&
                                lb_max(jset),npgfb(jset),zetb(:,jset),&
                                lca_max(kset),npgfca(kset),zetca(:,kset),&
                                rab,swork,calc_aba=.TRUE.,error=error)

             swork_cont = 0.0_dp
             DO ishella = 1, nshella(iset)
                lai = la(ishella,iset)
                consta=SQRT(1.0_dp/REAL(2*lai+1,dp))
                DO jshellb = 1, nshellb(jset)
                   lbj = lb(jshellb,jset)
                   DO kshella = 1, nshellca(kset)
                      lak = lca(kshella,kset)
                      constca=SQRT(1.0_dp/REAL(2*lak+1,dp))
                      nl = INT((lai + lak)/2) 
                      IF(lai == 0 .OR. lak == 0) nl = 0
                      DO il = 0, nl
                         laa = lai + lak - 2*il
                         nds = laa + lbj + 1
                         ids_start = nds - MIN(laa,lbj) 
                         const = consta * constca * 2.0_dp* SQRT(pi*REAL(2*laa+1,dp))
                         DO iil = 1,il
                           shg_fac(iil) = fac(laa+iil-1) * ifac(laa) * REAL(laa,dp)&
                                          * fac(il)/fac(il-iil)/fac(iil) 
                         ENDDO
                         DO ipgfa = 1, npgfa(iset)
                            DO jpgfb = 1, npgfb(jset)
                               bjf  = 1.0_dp/((2._dp*zetb(jpgfb,jset))**lbj)
                               DO kpgfa = 1, npgfca(kset)
                                  a = zeta(ipgfa,iset)+zetca(kpgfa,kset)
                                  aif  = 1.0_dp/((2._dp*a)**laa)
                                  prefac = const*aif*bjf&
                                           *scona_shg(ipgfa,ishella,iset)&
                                           *sconb_shg(jpgfb,jshellb,jset)&
                                           *sconca_shg(kpgfa,kshella,kset)
                                  DO ids = ids_start, nds 
                                     IF(il==0) THEN
                                       swork_cont(ishella,jshellb,kshella,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshella,il,ids)&
                                                 + prefac * swork(ipgfa,jpgfb,kpgfa,il,ids)
                                     ELSEIF(il==1) THEN
                                       swork_cont(ishella,jshellb,kshella,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshella,il,ids)&
                                                 + prefac * (swork(ipgfa,jpgfb,kpgfa,il,ids)&
                                                  + shg_fac(il) * swork(ipgfa,jpgfb,kpgfa,il-1,ids)/a) 
                                     ELSEIF(il==2) THEN
                                       swork_cont(ishella,jshellb,kshella,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshella,il,ids)&
                                                 +prefac*(swork(ipgfa,jpgfb,kpgfa,il,ids)&
                                                  + shg_fac(il-1) * swork(ipgfa,jpgfb,kpgfa,il-1,ids)/a &
                                                  + shg_fac(il) * swork(ipgfa,jpgfb,kpgfa,il-2,ids)/(a**2)) 
                                     ELSEIF(il==3) THEN
                                       swork_cont(ishella,jshellb,kshella,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshella,il,ids)&
                                                 +prefac*(swork(ipgfa,jpgfb,kpgfa,il,ids)&
                                                   + shg_fac(il-2) * swork(ipgfa,jpgfb,kpgfa,il-1,ids)/a &
                                                   + shg_fac(il-1) * swork(ipgfa,jpgfb,kpgfa,il-2,ids)/(a**2)& 
                                                   + shg_fac(il) * swork(ipgfa,jpgfb,kpgfa,il-3,ids)/(a**3))
                                     ELSEIF(il==4) THEN
                                       swork_cont(ishella,jshellb,kshella,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshella,il,ids)&
                                                 +prefac*(swork(ipgfa,jpgfb,kpgfa,il,ids)&
                                                   + shg_fac(il-3) * swork(ipgfa,jpgfb,kpgfa,il-1,ids)/a&
                                                   + shg_fac(il-2) * swork(ipgfa,jpgfb,kpgfa,il-2,ids)/(a**2) &
                                                   + shg_fac(il-1) * swork(ipgfa,jpgfb,kpgfa,il-3,ids)/(a**3)& 
                                                   + shg_fac(il) * swork(ipgfa,jpgfb,kpgfa,il-4,ids)/(a**4)) 
                                     ELSEIF(il==5) THEN
                                       swork_cont(ishella,jshellb,kshella,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshella,il,ids)&
                                                 +prefac*(swork(ipgfa,jpgfb,kpgfa,il,ids)&
                                                   + shg_fac(il-4) * swork(ipgfa,jpgfb,kpgfa,il-1,ids)/a&
                                                   + shg_fac(il-3) * swork(ipgfa,jpgfb,kpgfa,il-2,ids)/a**2&
                                                   + shg_fac(il-2) * swork(ipgfa,jpgfb,kpgfa,il-3,ids)/(a**3) &
                                                   + shg_fac(il-1) * swork(ipgfa,jpgfb,kpgfa,il-4,ids)/(a**4)& 
                                                   + shg_fac(il) * swork(ipgfa,jpgfb,kpgfa,il-5,ids)/(a**5))
                                     ELSE 
                                      CALL cp_unimplemented_error(fromWhere=routineP, &
                                           message="SHG integrals not implemented when l quantum number"//&
                                           " of orbital and ri basis larger than 11", &
                                           error=error, error_level=cp_failure_level)
                                     ENDIF
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             CALL overlap_shg_aba(la(:,iset),first_sgfa(:,iset),nshella(iset),&
                                   lb(:,jset),first_sgfb(:,jset),nshellb(jset),&
                                   lca(:,kset),first_sgfca(:,kset),nshellca(kset),&
                                   swork_cont,abaint,Waux_mat,lri_env%cg_shg%cg_coeff,&
                                   lri_env%cg_shg%cg_none0_list,lri_env%cg_shg%ncg_none0)

          END DO
       END DO
    END DO

    DEALLOCATE(swork_cont,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(swork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE lri_int_aba_shg

! *****************************************************************************
!> \brief calculate integrals (a,b,fb)
!> \param lri_env ...
!> \param abbint integral (a,b,fb)
!> \param rab distance vector between A and B
!> \param oba orbital basis at center A
!> \param obb orbital basis at center B
!> \param fbb auxiliary basis set at center B
!> \param scona_shg contraction matrix for orb bas on A
!> \param sconb_shg contraction matrix for orb bas on B
!> \param sconcb_shg contraction matrix for aux bas on B
!> \param Waux_mat angular-dependent part
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_int_abb_shg(lri_env,abbint,rab,oba,obb,fbb,scona_shg,sconb_shg,&
                             sconcb_shg,Waux_mat,error)

    TYPE(lri_environment_type), POINTER      :: lri_env
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: abbint
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    TYPE(gto_basis_set_type), POINTER        :: oba, obb, fbb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: scona_shg, sconb_shg, &
                                                sconcb_shg
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: Waux_mat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_abb_shg', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ids, ids_start, iil, il, ipgfa, iset, ishella, jpgfb, &
      jset, jshellb, kpgfb, kset, kshellb, la_max_set, lai, lb_max_set, lbb, &
      lbb_max_set, lbj, lbk, lcb_max_set, nds, nl, nl_set, npgfa_set, &
      npgfb_set, npgfcb_set, nseta, nsetb, nsetcb, nshella_set, nshellb_set, &
      nshellcb_set, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, lb_max, lcb_max, &
                                                npgfa, npgfb, npgfcb, &
                                                nshella, nshellb, nshellcb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, &
                                                first_sgfcb, la, lb, lcb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aif, b, bjf, const, constb, &
                                                constcb, dab, prefac, rab2, &
                                                shg_fac(5)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b, &
                                                set_radius_cb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta, zetb, zetcb
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), POINTER      :: swork, swork_cont

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, lb_max, lcb_max, npgfa, npgfb, npgfcb) 
    NULLIFY(first_sgfa, first_sgfb, first_sgfcb, set_radius_a, set_radius_b,&
            set_radius_cb, zeta, zetb, zetcb) 

    ! basis ikind
    first_sgfa    =>  oba%first_sgf
    la_max        =>  oba%lmax
    la            =>  oba%l
    npgfa         =>  oba%npgf
    nshella       =>  oba%nshell
    nseta         =   oba%nset
    set_radius_a  =>  oba%set_radius
    zeta          =>  oba%zet
    ! basis jkind
    first_sgfb    =>  obb%first_sgf
    lb_max        =>  obb%lmax
    lb            =>  obb%l
    npgfb         =>  obb%npgf
    nshellb       =>  obb%nshell
    nsetb         =   obb%nset
    set_radius_b  =>  obb%set_radius
    zetb          =>  obb%zet

    ! basis RI on B
    first_sgfcb    =>  fbb%first_sgf
    lcb_max        =>  fbb%lmax
    lcb            =>  fbb%l
    npgfcb         =>  fbb%npgf
    nshellcb       =>  fbb%nshell
    nsetcb         =   fbb%nset
    set_radius_cb  =>  fbb%set_radius
    zetcb          =>  fbb%zet

    dab = SQRT( SUM(rab**2) )
    rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)

    la_max_set   = MAXVAL(la_max)
    lb_max_set   = MAXVAL(lb_max)
    lcb_max_set  = MAXVAL(lcb_max)
    npgfa_set    = MAXVAL(npgfa)
    npgfb_set    = MAXVAL(npgfb)
    npgfcb_set   = MAXVAL(npgfcb)
    nshella_set  = MAXVAL(nshella)
    nshellb_set  = MAXVAL(nshellb)
    nshellcb_set = MAXVAL(nshellcb)

    lbb_max_set = lb_max_set + lcb_max_set

    ! allocate some work storage....
    nds = la_max_set + lbb_max_set + 1
    nl_set = INT((lbb_max_set)/2)
    ALLOCATE(swork(1:npgfa_set,1:npgfb_set,1:npgfcb_set,0:nl_set,1:nds),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(swork_cont(nshella_set,nshellb_set,nshellcb_set,0:nl_set,nds),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO iset=1,nseta

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE

          DO kset=1,nsetcb

             IF (set_radius_a(iset) + set_radius_cb(kset) < dab) CYCLE

             CALL s_overlap_abx(la_max(iset),npgfa(iset),zeta(:,iset),&
                                lb_max(jset),npgfb(jset),zetb(:,jset),&
                                lcb_max(kset),npgfcb(kset),zetcb(:,kset),&
                                rab,swork,calc_aba=.FALSE.,error=error)

             swork_cont = 0.0_dp
             DO ishella = 1, nshella(iset)
                lai = la(ishella,iset)
                DO jshellb = 1, nshellb(jset)
                   lbj = lb(jshellb,jset)
                   constb=SQRT(1.0_dp/REAL(2*lbj+1,dp))
                   DO kshellb = 1, nshellcb(kset)
                      lbk = lcb(kshellb,kset)
                      constcb=SQRT(1.0_dp/REAL(2*lbk+1,dp))
                      nl = INT((lbj + lbk)/2) 
                      IF(lbj == 0 .OR. lbk == 0) nl = 0
                      DO il = 0, nl
                         lbb = lbj + lbk - 2*il
                         nds = lai + lbb + 1
                         ids_start = nds - MIN(lai,lbb)
                         const = constb * constcb * 2.0_dp* SQRT(pi*REAL(2*lbb+1,dp))
                         DO iil = 1,il
                           shg_fac(iil) = fac(lbb+iil-1) * ifac(lbb) * REAL(lbb,dp)&
                                          * fac(il)/fac(il-iil)/fac(iil) 
                         ENDDO
                         DO ipgfa = 1, npgfa(iset)
                            aif  = 1.0_dp/((2._dp*zeta(ipgfa,iset))**lai)
                            DO jpgfb = 1, npgfb(jset)
                               DO kpgfb = 1, npgfcb(kset)
                                  b = zetb(jpgfb,jset)+zetcb(kpgfb,kset)
                                  bjf  = 1.0_dp/((2._dp*b)**lbb)
                                  prefac = const*aif*bjf&
                                           *scona_shg(ipgfa,ishella,iset)&
                                           *sconb_shg(jpgfb,jshellb,jset)&
                                           *sconcb_shg(kpgfb,kshellb,kset)
                                  DO ids = ids_start, nds 
                                     IF(il==0) THEN
                                       swork_cont(ishella,jshellb,kshellb,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshellb,il,ids)&
                                                 + prefac * swork(ipgfa,jpgfb,kpgfb,il,ids)
                                     ELSEIF(il==1) THEN
                                       swork_cont(ishella,jshellb,kshellb,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshellb,il,ids)&
                                                 + prefac * (swork(ipgfa,jpgfb,kpgfb,il,ids)&
                                                  + shg_fac(il) * swork(ipgfa,jpgfb,kpgfb,il-1,ids)/b) 
                                     ELSEIF(il==2) THEN
                                       swork_cont(ishella,jshellb,kshellb,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshellb,il,ids)&
                                                 +prefac*(swork(ipgfa,jpgfb,kpgfb,il,ids)&
                                                  + shg_fac(il-1) * swork(ipgfa,jpgfb,kpgfb,il-1,ids)/b &
                                                  + shg_fac(il) * swork(ipgfa,jpgfb,kpgfb,il-2,ids)/(b**2)) 
                                     ELSEIF(il==3) THEN
                                       swork_cont(ishella,jshellb,kshellb,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshellb,il,ids)&
                                                 +prefac*(swork(ipgfa,jpgfb,kpgfb,il,ids)&
                                                   + shg_fac(il-2) * swork(ipgfa,jpgfb,kpgfb,il-1,ids)/b &
                                                   + shg_fac(il-1) * swork(ipgfa,jpgfb,kpgfb,il-2,ids)/(b**2)& 
                                                   + shg_fac(il) * swork(ipgfa,jpgfb,kpgfb,il-3,ids)/(b**3))
                                     ELSEIF(il==4) THEN
                                       swork_cont(ishella,jshellb,kshellb,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshellb,il,ids)&
                                                 +prefac*(swork(ipgfa,jpgfb,kpgfb,il,ids)&
                                                   + shg_fac(il-3) * swork(ipgfa,jpgfb,kpgfb,il-1,ids)/b&
                                                   + shg_fac(il-2) * swork(ipgfa,jpgfb,kpgfb,il-2,ids)/(b**2) &
                                                   + shg_fac(il-1) * swork(ipgfa,jpgfb,kpgfb,il-3,ids)/(b**3)& 
                                                   + shg_fac(il) * swork(ipgfa,jpgfb,kpgfb,il-4,ids)/(b**4)) 
                                     ELSEIF(il==5) THEN
                                       swork_cont(ishella,jshellb,kshellb,il,ids) = &
                                                 swork_cont(ishella,jshellb,kshellb,il,ids)&
                                                 +prefac*(swork(ipgfa,jpgfb,kpgfb,il,ids)&
                                                   + shg_fac(il-4) * swork(ipgfa,jpgfb,kpgfb,il-1,ids)/b&
                                                   + shg_fac(il-3) * swork(ipgfa,jpgfb,kpgfb,il-2,ids)/b**2&
                                                   + shg_fac(il-2) * swork(ipgfa,jpgfb,kpgfb,il-3,ids)/(b**3) &
                                                   + shg_fac(il-1) * swork(ipgfa,jpgfb,kpgfb,il-4,ids)/(b**4)& 
                                                   + shg_fac(il) * swork(ipgfa,jpgfb,kpgfb,il-5,ids)/(b**5))
                                     ELSE 
                                      CALL cp_unimplemented_error(fromWhere=routineP, &
                                           message="SHG integrals not implemented when l quantum number"//&
                                           " of orbital and ri basis larger than 11", &
                                           error=error, error_level=cp_failure_level)
                                     ENDIF
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
             CALL overlap_shg_abb(la(:,iset),first_sgfa(:,iset),nshella(iset),&
                                    lb(:,jset),first_sgfb(:,jset),nshellb(jset),&
                                    lcb(:,kset),first_sgfcb(:,kset),nshellcb(kset),&
                                    swork_cont,abbint,Waux_mat,lri_env%cg_shg%cg_coeff,&
                                    lri_env%cg_shg%cg_none0_list,lri_env%cg_shg%ncg_none0)

          END DO
       END DO
    END DO

    DEALLOCATE(swork_cont,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(swork,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE lri_int_abb_shg

! *****************************************************************************
!> \brief calculate overlap integrals (a,b)
!> \param sab integral (a,b)
!> \param rab distance vector
!> \param fba basis at center A
!> \param fbb basis at center B
!> \param scona_shg contraction matrix A
!> \param sconb_shg contraxtion matrix B
!> \param Waux_mat angular-dependent part
!> \param error ...
! *****************************************************************************
  SUBROUTINE lri_int_ab_shg(sab,rab,fba,fbb,scona_shg,sconb_shg,Waux_mat,error)

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: sab
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    TYPE(gto_basis_set_type), POINTER        :: fba, fbb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: scona_shg, sconb_shg
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: Waux_mat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lri_int_ab_shg', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ids, ids_start, ipgfa, iset, ishella, jpgfb, jset, &
      jshellb, la_max_set, lai, lb_max_set, lbj, nds, npgfa_set, npgfb_set, &
      nseta, nsetb, nshella_set, nshellb_set, stat
    INTEGER, DIMENSION(:), POINTER           :: la_max, lb_max, npgfa, npgfb, &
                                                nshella, nshellb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa, first_sgfb, la, lb
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: aif, bjf, dab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a, set_radius_b
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: zeta, zetb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: swork, swork_cont

    failure = .FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(la_max, lb_max, npgfa, npgfb, first_sgfa, first_sgfb, set_radius_a,&
             set_radius_b, zeta, zetb)

    ! basis ikind
    first_sgfa   =>  fba%first_sgf
    la_max       =>  fba%lmax
    la           =>  fba%l
    npgfa        =>  fba%npgf
    nseta        =   fba%nset
    set_radius_a =>  fba%set_radius
    zeta         =>  fba%zet
    nshella      =>  fba%nshell
    ! basis jkind
    first_sgfb   =>  fbb%first_sgf
    lb_max       =>  fbb%lmax
    lb           =>  fbb%l
    npgfb        =>  fbb%npgf
    nsetb        =   fbb%nset
    set_radius_b =>  fbb%set_radius
    zetb         =>  fbb%zet
    nshellb      =>  fbb%nshell

    dab = SQRT( SUM(rab**2) )

    la_max_set = MAXVAL(la_max)
    lb_max_set = MAXVAL(lb_max)

    ! allocate some work matrices
    npgfa_set = MAXVAL(npgfa)
    npgfb_set = MAXVAL(npgfb)
    nshella_set = MAXVAL(nshella)
    nshellb_set = MAXVAL(nshellb)
    ALLOCATE(swork(npgfa_set,npgfb_set,la_max_set+lb_max_set+1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(swork_cont(nshella_set,nshellb_set,la_max_set+lb_max_set+1),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO iset=1,nseta

       DO jset=1,nsetb

          IF (set_radius_a(iset) + set_radius_b(jset) < dab) CYCLE


          CALL s_overlap_ab(la_max(iset),npgfa(iset),zeta(:,iset),&
                            lb_max(jset),npgfb(jset),zetb(:,jset),&
                            rab,swork)

          swork_cont = 0.0_dp
          DO ishella = 1, nshella(iset)
             lai = la(ishella,iset)
             DO jshellb = 1, nshellb(jset)
                lbj = lb(jshellb,jset)
                nds = lai + lbj + 1
                ids_start = nds - MIN(lai,lbj)
                DO ipgfa = 1, npgfa(iset)
                   aif  = 1.0_dp/((2._dp*zeta(ipgfa,iset))**lai)
                   DO jpgfb = 1, npgfb(jset)
                      bjf  = 1.0_dp/((2._dp*zetb(jpgfb,jset))**lbj)
                      DO ids =ids_start, nds 
                         swork_cont(ishella,jshellb,ids) = swork_cont(ishella,jshellb,ids)&
                                                          +scona_shg(ipgfa,ishella,iset)&
                                                          *sconb_shg(jpgfb,jshellb,jset)&
                                                          *aif*bjf*swork(ipgfa,jpgfb,ids) 
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          CALL overlap_shg_ab(la(:,iset),first_sgfa(:,iset),nshella(iset),&
                              lb(:,jset),first_sgfb(:,jset),nshellb(jset),&
                              swork_cont,sab,Waux_mat)
       END DO
    END DO

    DEALLOCATE(swork,swork_cont,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE lri_int_ab_shg

END MODULE lri_shg_integrals
