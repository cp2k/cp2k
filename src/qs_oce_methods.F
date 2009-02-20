!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the construction of the coefficients
!>      for the expansion  of the atomic
!>      densities rho1_hard and rho1_soft in terms of primitive spherical gaussians.
!> \par History
!>      05-2004 created
!> \author MI
! *****************************************************************************
MODULE qs_oce_methods

  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE block_p_types,                   ONLY: block_p_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: dfac,&
                                             pi
  USE orbital_pointers,                ONLY: indco,&
                                             init_orbital_pointers,&
                                             nco,&
                                             ncoset,&
                                             nso
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next
  USE qs_util,                         ONLY: exp_radius
  USE sap_kind_types,                  ONLY: clist_type,&
                                             sap_int_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_oce_methods'

! *** Public subroutines ***

  PUBLIC :: build_oce_matrices,proj_blk

CONTAINS

! *****************************************************************************
  SUBROUTINE build_oce_block(oceh,oces,atom_ka,atom_kb,rab,dab,nder,sgf_list,nsgf_cnt,sgf_soft_only,&
       eps_fit,error)

    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces
    TYPE(atomic_kind_type), POINTER          :: atom_ka, atom_kb
    REAL(dp), DIMENSION(3)                   :: rab
    REAL(dp)                                 :: dab
    INTEGER, INTENT(IN)                      :: nder
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sgf_list
    INTEGER, INTENT(OUT)                     :: nsgf_cnt
    LOGICAL, INTENT(OUT)                     :: sgf_soft_only
    REAL(dp), INTENT(IN)                     :: eps_fit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'build_oce_block', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_col, ic, ico, ider, igau, ip, ipgf, iprjc, iprjs, is, &
      iset, isgfb, isgfb_cnt, iso, isp, istat, ithread, jc, jset, lds, lm, &
      lpoint, lprj, lsgfb, lsgfb_cnt, lshell, lx, ly, lz, m, m1, maxcob, &
      maxder, maxlb, maxlprj, maxnprja, maxsoa, n, ncob, np_car, np_sph, &
      nseta, nsetb, nsoatot, nthread, ntotsgfb, omp_get_max_threads, &
      omp_get_thread_num, sgf_hard_only
    INTEGER, DIMENSION(:), POINTER           :: fp_cara, fp_spha, lb_max, &
                                                lb_min, npgfb, nprjla, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: calculate_forces, failure, &
                                                paw_atom_a, paw_atom_b
    REAL(dp)                                 :: hard_radius_a, hard_radius_b, &
                                                radius, rcprja
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: s
    REAL(dp), DIMENSION(:), POINTER          :: set_radius_b, zisominb
    REAL(dp), DIMENSION(:, :), POINTER :: cprj_s, local_oce_h, local_oce_s, &
      ovc, ovs, rpgfb, rzetprja, spa_sb, spa_tmp, sphi_b, zetb, zetprja
    REAL(dp), DIMENSION(:, :, :), POINTER    :: sdpa_tmp, spa_sbt, spa_tmpt
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: sdpa_tmpt
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_a, paw_proj_b

   failure = .FALSE.
   NULLIFY(orb_basis_a,paw_proj_a)
   CALL get_atomic_kind(atomic_kind=atom_ka,orb_basis_set=orb_basis_a,&
                        paw_proj_set=paw_proj_a,paw_atom=paw_atom_a,&
                        hard_radius=hard_radius_a)

   NULLIFY(orb_basis_b,paw_proj_b)
   CALL get_atomic_kind(atomic_kind=atom_kb,orb_basis_set=orb_basis_b,&
                         paw_proj_set=paw_proj_b,paw_atom=paw_atom_b,&
                         hard_radius=hard_radius_b)

   IF(.NOT. paw_atom_a) RETURN

   NULLIFY(cprj_s,nprjla,fp_cara,fp_spha,rzetprja,zetprja)
   CALL get_paw_proj_set(paw_proj_set=paw_proj_a,cprj_s=cprj_s,maxl=maxlprj, &
                         nprj=nprjla,ncgauprj=np_car,nsgauprj=np_sph,rcprj=rcprja, &
                         first_prj=fp_cara,first_prjs=fp_spha,&
                         rzetprj=rzetprja,zetprj=zetprja)

   NULLIFY(first_sgfb,lb_max,lb_min,npgfb,nsgfb,rpgfb,sphi_b,set_radius_b,zetb,zisominb)
   CALL get_gto_basis_set(gto_basis_set=orb_basis_b,nset=nsetb, nsgf=ntotsgfb, &
                          set_radius=set_radius_b,lmax=lb_max,lmin=lb_min, &
                          npgf=npgfb,nsgf_set=nsgfb,pgf_radius=rpgfb, &
                          sphi=sphi_b,zet=zetb,first_sgf=first_sgfb, &
                          maxco=maxcob,maxl=maxlb)

   CALL get_gto_basis_set(gto_basis_set=orb_basis_a,nset=nseta,maxso=maxsoa)

!  *** Add the block ab ***

    NULLIFY(spa_sbt,spa_tmpt, sdpa_tmpt)
    NULLIFY(spa_tmp,ovc,ovs,sdpa_tmp)

    nthread = 1
!$  nthread = omp_get_max_threads()

    maxder = ncoset(nder)
    nsoatot = maxsoa*nseta
    maxnprja = SIZE(zetprja,1)

    calculate_forces = .FALSE.
    IF(nder > 0) THEN
      calculate_forces = .TRUE.
    ENDIF

    lm = MAX(maxlb, maxlprj)
    lds = ncoset(lm+nder+1)

    ALLOCATE(s(lds,lds,ncoset(nder+1)),STAT= istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(spa_sbt(np_car,ntotsgfb,0:nthread-1), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    spa_sbt(1:np_car,1:ntotsgfb,0:nthread-1) =0.0_dp
    ALLOCATE(spa_tmpt(maxnprja*ncoset(maxlprj),maxcob*maxder,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    spa_tmpt(1:maxnprja*ncoset(maxlprj),1:maxcob*maxder,0:nthread-1) =0.0_dp
    ALLOCATE(sdpa_tmpt(maxnprja*ncoset(maxlprj+1),maxcob,4,0:nthread-1),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    sdpa_tmpt(1:maxnprja*ncoset(maxlprj+1),1:maxcob,1:4,0:nthread-1) =0.0_dp
    ALLOCATE(ovc(np_car,maxcob*nsetb*maxder),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ovc = 0.0_dp
    ALLOCATE(ovs(np_sph,maxcob*nsetb*maxder),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ovs = 0.0_dp

    ithread = 0
!$  ithread = omp_get_thread_num()
    spa_sb => spa_sbt(:,:,ithread)
    spa_tmp => spa_tmpt(:,:,ithread)
    sdpa_tmp => sdpa_tmpt(:,:,:,ithread)
    m1 = 0
    nsgf_cnt = 0
    isgfb_cnt = 1
    sgf_hard_only = 0
    DO jset = 1,nsetb
       !
       ! Set the contribution list
       IF(hard_radius_a + set_radius_b(jset) >= dab) THEN

          ! check if this function is hard
          radius = exp_radius(lb_max(jset),MAXVAL(zetb(1:npgfb(jset),jset)),eps_fit,1.0_dp)
          IF(radius.LE.hard_radius_b) sgf_hard_only = sgf_hard_only + 1

          isgfb = first_sgfb(1,jset)
          lsgfb = isgfb - 1 + nsgfb(jset)
          DO jc = isgfb,lsgfb
             nsgf_cnt = nsgf_cnt + 1
             sgf_list(nsgf_cnt) = jc
          ENDDO

!     ***integral between proj of iatom and primitives of jatom
!     *** Calculate the primitives overlap ***
        spa_tmp = 0.0_dp
        spa_sb  = 0.0_dp
        sdpa_tmp = 0.0_dp

        s = 0.0_dp
        ncob = npgfb(jset)*ncoset(lb_max(jset))

        isgfb = first_sgfb(1,jset)
        lsgfb = isgfb - 1 + nsgfb(jset)

        lsgfb_cnt = isgfb_cnt - 1 + nsgfb(jset)

        DO lprj = 0,maxlprj
          IF (calculate_forces) THEN
            CALL overlap(lprj,lprj,nprjla(lprj),&
                       rzetprja(:,lprj),zetprja(:,lprj),&
                       lb_max(jset),lb_min(jset),npgfb(jset),&
                       rpgfb(:,jset),zetb(:,jset),&
                       -rab,dab,spa_tmp,&
                       0,.FALSE.,s,lds,sdpa_tmp)

            DO ipgf = 1,nprjla(lprj)
              n = (ipgf -1)*ncoset(lprj)
              lpoint = ncoset(lprj -1)+1+ n
              m = fp_cara(lprj)+(ipgf -1)*nco(lprj)
              is = 0
              isp =  0
              DO ip = 1,npgfb(jset)
                 ic = (ip-1)*ncoset(lb_max(jset))
                 DO jc = ncoset(lb_min(jset)-1)+1,ncoset(lb_max(jset))
                   ovc(m:m-1+nco(lprj),jc+ic+isp+m1) =&
                     spa_tmp(lpoint:lpoint-1+nco(lprj),ic+jc+is)
                 ENDDO
              ENDDO
            ENDDO
            DO ider =2,maxder
              isp =  (ider-1)*maxcob*nsetb
              DO ipgf = 1,nprjla(lprj)
                n = (ipgf -1)*ncoset(lprj)
                lpoint = ncoset(lprj -1)+1+ n
                m = fp_cara(lprj)+(ipgf -1)*nco(lprj)
                DO ip = 1,npgfb(jset)
                   ic = (ip-1)*ncoset(lb_max(jset))
                   DO jc = ncoset(lb_min(jset)-1)+1,ncoset(lb_max(jset))
                     ovc(m:m-1+nco(lprj),jc+ic+isp+m1) =&
                       sdpa_tmp(lpoint:lpoint-1+nco(lprj),ic+jc,ider)
                   ENDDO
                ENDDO
              ENDDO
            ENDDO

          ELSE
            CALL overlap(lprj,lprj,nprjla(lprj),&
                       rzetprja(:,lprj),zetprja(:,lprj),&
                       lb_max(jset),lb_min(jset),npgfb(jset),&
                       rpgfb(:,jset),zetb(:,jset),&
                       -rab,dab,spa_tmp,&
                       nder,.FALSE.,s,lds)

            DO ipgf = 1,nprjla(lprj)
              n = (ipgf -1)*ncoset(lprj)
              lpoint = ncoset(lprj -1)+1+ n
              m = fp_cara(lprj)+(ipgf -1)*nco(lprj)

              is = 0
              isp =  0
              DO ip = 1,npgfb(jset)
                 ic = (ip-1)*ncoset(lb_max(jset))
                 DO jc = ncoset(lb_min(jset)-1)+1,ncoset(lb_max(jset))
                   ovc(m:m-1+nco(lprj),jc+ic+m1) =&
                     spa_tmp(lpoint:lpoint-1+nco(lprj),ic+jc)
                 ENDDO
              ENDDO

            ENDDO

         END IF
      ENDDO

        n = ncoset(lb_max(jset))
        DO ider =1,maxder
          is = (ider-1)*maxcob
          isp =  (ider-1)*maxcob*nsetb
          DO ipgf = 1,npgfb(jset)
            DO lshell = lb_min(jset),lb_max(jset)
              DO ic = 1,nco(lshell)
                igau = ic + ncoset(lshell-1) + n*(ipgf-1) + m1 + isp
                DO lprj = 0,maxlprj
                  DO ip = 1,nprjla(lprj)
                    DO iso = 1,nso(lprj)
                      iprjs = fp_spha(lprj) + iso - 1 + nso(lprj)*(ip-1)
                      DO ico = 1,nco(lprj)
                        lx = indco(1,ico+ncoset(lprj-1))
                        ly = indco(2,ico+ncoset(lprj-1))
                        lz = indco(3,ico+ncoset(lprj-1))
                        iprjc = fp_cara(lprj) + ico - 1 + nco(lprj)*(ip-1)
                        ovs(iprjs,igau) = ovs(iprjs,igau) + &
                           orbtramat(lprj)%c2s(iso,ico) * ovc(iprjc,igau) /&
                           SQRT((4.0_dp*pi)/dfac(2*lprj+1)*&
                           dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))
                      ENDDO ! ico
                    ENDDO ! iso
                  ENDDO ! ip
                ENDDO ! lprj
              ENDDO ! ic
            ENDDO ! lshell
          ENDDO ! ipgf
        ENDDO ! ider

        IF(paw_atom_b) THEN
          CALL get_paw_proj_set(paw_proj_set=paw_proj_b,zisomin=zisominb)
          DO ipgf = 1,npgfb(jset)
            DO lshell = lb_min(jset),lb_max(jset)
              IF(zetb(ipgf,jset) >= zisominb(lshell)) THEN
                igau = n*(ipgf-1) + ncoset(lshell - 1)
                DO ider =1,maxder
                  is = maxcob*(ider-1)
                  isp =  (ider-1)*maxcob*nsetb
                  ovs(:,igau+1+isp+m1:igau+nco(lshell)+isp+m1) = 0.0_dp
                ENDDO
              ENDIF
            ENDDO
          END DO
        ENDIF

!       *** Contraction step (integrals and derivatives)
        IF(dab <= 1.0E-8_dp)THEN
          iset = jset
          n = maxsoa*(iset-1)

          NULLIFY(local_oce_h,local_oce_s)
          CALL get_paw_proj_set(paw_proj_set=paw_proj_a, &
                                local_oce_sphi_h=local_oce_h, &
                                local_oce_sphi_s=local_oce_s)
          oceh(1)%block(n+1:n+maxsoa,isgfb:lsgfb) = &
                    local_oce_h(1:maxsoa,isgfb:lsgfb)
          oces(1)%block(n+1:n+maxsoa,isgfb:lsgfb) = &
                    local_oce_s(1:maxsoa,isgfb:lsgfb)

          DO ider =2,maxder
            first_col = (ider - 1)*maxcob*nsetb + 1 + m1
            CALL dgemm("N","N",np_sph,nsgfb(jset),ncob,&
                      1.0_dp,ovs(1,first_col),SIZE(ovs,1),&
                      paw_proj_a%sphi_h(1,isgfb),SIZE(paw_proj_a%sphi_h,1),&
                      0.0_dp,spa_sb(1,isgfb),SIZE(spa_sb,1))

            CALL dgemm("T","N",nsoatot,nsgfb(jset),np_sph,&
                      1.0_dp,cprj_s(1,1),SIZE(cprj_s,1),&
                      spa_sb(1,isgfb),SIZE(spa_sb,1),&
                      1.0_dp,oceh(ider)%block(1,isgfb),SIZE(oceh(ider)%block,1))

            CALL dgemm("N","N",np_sph,nsgfb(jset),ncob,&
                      1.0_dp,ovs(1,first_col),SIZE(ovs,1),&
                      paw_proj_a%sphi_s(1,isgfb),SIZE(paw_proj_a%sphi_s,1),&
                      0.0_dp,spa_sb(1,isgfb),SIZE(spa_sb,1))

            CALL dgemm("T","N",nsoatot,nsgfb(jset),np_sph,&
                      1.0_dp,cprj_s(1,1),SIZE(cprj_s,1),&
                      spa_sb(1,isgfb),SIZE(spa_sb,1),&
                      1.0_dp,oces(ider)%block(1,isgfb),SIZE(oces(ider)%block,1))
          ENDDO
        ELSE
          DO ider = 1,maxder
            first_col = (ider - 1)*maxcob*nsetb + 1 + m1
            CALL dgemm("N","N",np_sph,nsgfb(jset),ncob,&
                      1.0_dp,ovs(1,first_col),SIZE(ovs,1),&
                      sphi_b(1,isgfb),SIZE(sphi_b,1),&
                      0.0_dp,spa_sb(1,isgfb),SIZE(spa_sb,1))

            CALL dgemm("T","N",nsoatot,nsgfb(jset),np_sph,&
                      1.0_dp,cprj_s(1,1),SIZE(cprj_s,1),&
                      spa_sb(1,isgfb),SIZE(spa_sb,1),&
                      1.0_dp,oceh(ider)%block(1,isgfb_cnt),SIZE(oceh(ider)%block,1))

            oces(ider)%block(1:nsoatot,isgfb_cnt:lsgfb_cnt) = &
                      oceh(ider)%block(1:nsoatot,isgfb_cnt:lsgfb_cnt)
          ENDDO
        END IF
        isgfb_cnt = isgfb_cnt + nsgfb(jset)
      END IF ! radius
      m1 = m1 + maxcob
    ENDDO  !jset

    ! check if the screened functions are all soft
    sgf_soft_only = .FALSE.
    IF(sgf_hard_only.EQ.0) sgf_soft_only = .TRUE.

    DEALLOCATE(s,spa_sbt,spa_tmpt,sdpa_tmpt,ovc,ovs ,STAT = istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE build_oce_block

! *****************************************************************************
!> \brief Set up the sparse matrix for the coefficients of one center expansions
!>      This routine uses the same logic as the nonlocal pseudopotential
!> \param intac: TYPE that holds the integrals (a=basis; c=projector)
!> \par History
!>      02.2009 created
!> \author jgh
! *****************************************************************************
  SUBROUTINE build_oce_matrices(intac, calculate_forces, nder,&
                    atomic_kind_set, particle_set, sap_oce, eps_fit, error)

    TYPE(sap_int_type), DIMENSION(:), &
      POINTER                                :: intac
    LOGICAL, INTENT(IN)                      :: calculate_forces
    INTEGER                                  :: nder
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sap_oce
    REAL(dp), INTENT(IN)                     :: eps_fit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_oce_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iac, ikind, ilist, jkind, &
      jneighbor, maxco, maxder, maxl, maxlgto, maxlprj, maxprj, maxsgf, &
      maxsoa, maxsob, mlprj, natom, ncoa_sum, nkind, nlist, nneighbor, &
      np_car, np_sph, nseta, nsetb, nsgf_cnt, nsgfa, nsobtot, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sgf_list
    INTEGER, DIMENSION(3)                    :: cell_b
    INTEGER, DIMENSION(:), POINTER           :: fp_car, fp_sph, la_max, &
                                                la_min, npgfa, nprjla, &
                                                nsgf_seta
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, paw_atom_b, &
                                                sgf_soft_only
    REAL(KIND=dp)                            :: dab, rcprj
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: cprj, rpgfa, rzetprj, sphi_a, &
                                                zeta, zetb
    TYPE(atomic_kind_type), POINTER          :: at_a, at_b, atomic_kind
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces
    TYPE(clist_type), POINTER                :: clist
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_paw, orb_basis_set
    TYPE(neighbor_list_type), POINTER        :: sap_oce_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sap_oce_neighbor_node
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_b

    failure = .FALSE.
    IF (calculate_forces) THEN
      CALL timeset(routineN//" (forces)",handle)
    ELSE
      CALL timeset(routineN,handle)
    ENDIF

    IF ( ASSOCIATED(sap_oce) ) THEN

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

      maxder = ncoset(nder)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               maxco=maxco,&
                               maxlgto=maxlgto,&
                               maxlprj=maxlprj,&
                               maxco_proj=maxprj,&
                               maxsgf=maxsgf)

      maxl = MAX(maxlgto,maxlprj)
      CALL init_orbital_pointers(maxl+nder+1)

      ALLOCATE(oceh(maxder),oces(maxder),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      DO i=1,nkind*nkind
        NULLIFY(intac(i)%alist)
        intac(i)%nalist=0
      END DO

      !calculate the overlap integrals <a|p>
      DO ikind=1,nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_set)

         IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                first_sgf=first_sgfa,&
                                lmax=la_max,&
                                lmin=la_min,&
                                nco_sum=ncoa_sum,&
                                maxso=maxsoa,&
                                npgf=npgfa,&
                                nset=nseta,&
                                nsgf=nsgfa,&
                                nsgf_set=nsgf_seta,&
                                pgf_radius=rpgfa,&
                                set_radius=set_radius_a,&
                                sphi=sphi_a,&
                                zet=zeta)

         DO jkind=1,nkind
           atomic_kind => atomic_kind_set(jkind)

           iac = ikind + nkind*(jkind - 1)
           IF (.NOT.ASSOCIATED(sap_oce(iac)%neighbor_list_set)) CYCLE

           NULLIFY(paw_proj_b)
           CALL get_atomic_kind(atomic_kind=atomic_kind,paw_proj_set=paw_proj_b,paw_atom=paw_atom_b)
           IF (.NOT.paw_atom_b) CYCLE

           CALL get_atomic_kind(atomic_kind=atomic_kind,orb_basis_set=orb_basis_paw)
           IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE
           CALL get_gto_basis_set(gto_basis_set=orb_basis_paw,maxso=maxsob,nset=nsetb)
           nsobtot = maxsob*nsetb

           CALL get_paw_proj_set(paw_proj_set=paw_proj_b,&
                                 cprj=cprj,&
                                 maxl=mlprj,&
                                 nprj=nprjla,&
                                 ncgauprj=np_car,&
                                 nsgauprj=np_sph,&
                                 rcprj=rcprj, &
                                 first_prj=fp_car,&
                                 first_prjs=fp_sph,&
                                 rzetprj=rzetprj,&
                                 zetprj=zetb)

           CALL get_neighbor_list_set(neighbor_list_set=sap_oce(iac)%neighbor_list_set,nlist=nlist)
           NULLIFY(sap_oce_neighbor_list)
           intac(iac)%a_kind = ikind
           intac(iac)%p_kind = jkind
           intac(iac)%nalist = nlist
           ALLOCATE(intac(iac)%alist(nlist),STAT=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

           DO ilist = 1, nlist
              IF ( .NOT. ASSOCIATED(sap_oce_neighbor_list) ) THEN
                sap_oce_neighbor_list => first_list(sap_oce(iac)%neighbor_list_set)
              ELSE
                sap_oce_neighbor_list => next(sap_oce_neighbor_list)
              END IF
              CALL get_neighbor_list(neighbor_list=sap_oce_neighbor_list,atom=atom_a,nnode=nneighbor)
              intac(iac)%alist(ilist)%aatom  = atom_a
              intac(iac)%alist(ilist)%nclist = nneighbor 
              ALLOCATE(intac(iac)%alist(ilist)%clist(nneighbor),STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

              sap_oce_neighbor_node => first_node(sap_oce_neighbor_list)
              DO jneighbor=1,nneighbor
                 CALL get_neighbor_node(neighbor_node=sap_oce_neighbor_node,neighbor=atom_b,r=rab,cell=cell_b)
                 dab = SQRT(SUM(rab*rab))
                 clist => intac(iac)%alist(ilist)%clist(jneighbor)
                 clist%catom = atom_b
                 clist%cell  = cell_b
                 clist%rac   = rab
                 clist%nsgf_cnt = 0
                 clist%maxac    = 0.0_dp
                 clist%maxach   = 0.0_dp
                 NULLIFY(clist%acint,clist%achint,clist%sgf_list)

                 DO i=1,maxder
                    ALLOCATE(oceh(i)%block(nsobtot,nsgfa),oces(i)%block(nsobtot,nsgfa),STAT=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                    oceh(i)%block=0._dp
                    oces(i)%block=0._dp
                 END DO
                 ALLOCATE(sgf_list(nsgfa),STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

                 at_a => atomic_kind_set(jkind)
                 at_b => atomic_kind_set(ikind)

                 CALL build_oce_block(oceh,oces,at_a,at_b,rab,dab,nder,sgf_list,nsgf_cnt,&
                      sgf_soft_only,eps_fit,error)

                 clist%sgf_soft_only = sgf_soft_only
                 clist%nsgf_cnt = nsgf_cnt
                 IF(nsgf_cnt.GT.0) THEN

                    ALLOCATE(clist%acint(nsgf_cnt,nsobtot,maxder),clist%sgf_list(nsgf_cnt),STAT=stat)
                    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                    clist%acint(:,:,:)=0._dp
                    clist%sgf_list(:)=HUGE(0)

                    IF( atom_a == atom_b .AND. ALL(cell_b == 0) ) THEN
                       CPPostcondition(nsgf_cnt.EQ.nsgfa,cp_failure_level,routineP,error,failure)
                       ! *** Special case: A=B
                       ALLOCATE(clist%achint(nsgfa,nsobtot,maxder),STAT=stat)
                       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                       clist%achint=0._dp
                       clist%acint(1:nsgfa,1:nsobtot,1) = TRANSPOSE(oces(1)%block(1:nsobtot,1:nsgfa))
                       clist%achint(1:nsgfa,1:nsobtot,1) = TRANSPOSE(oceh(1)%block(1:nsobtot,1:nsgfa))
                    ELSE
                       DO i=1,maxder
                          clist%acint(1:nsgf_cnt,1:nsobtot,i) = TRANSPOSE(oces(i)%block(1:nsobtot,1:nsgf_cnt))
                       END DO
                    END IF

                    clist%maxac=MAXVAL(ABS(clist%acint(:,:,1)))
                    clist%maxach=0._dp
                    clist%sgf_list(1:nsgf_cnt) = sgf_list(1:nsgf_cnt)
                 ENDIF

                 sap_oce_neighbor_node => next(sap_oce_neighbor_node)

                 DO i=1,maxder
                   DEALLOCATE(oceh(i)%block,oces(i)%block,STAT=stat)
                   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
                 END DO
                 DEALLOCATE(sgf_list,STAT=stat)
                 CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
              END DO
           END DO
        END DO
     END DO

     DEALLOCATE(oceh,oces,STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

   END IF 

   CALL timestop(handle)
    
 END SUBROUTINE build_oce_matrices

!==========================================================================================================

! *****************************************************************************
!> \brief Project a matrix block onto the local atomic functions.
!>        
!> \par History
!>      02.2009 created
!> \author vw
! *****************************************************************************
  SUBROUTINE proj_blk(h_a,s_a,na,h_b,s_b,nb,blk,ldb,proj_h,proj_s,nso,buf1,buf2,fac,distab)

    INTEGER                                  :: na
    REAL(dp)                                 :: s_a(na,*), h_a(na,*)
    INTEGER                                  :: nb
    REAL(dp)                                 :: s_b(nb,*), h_b(nb,*)
    INTEGER                                  :: ldb
    REAL(dp)                                 :: blk(ldb,*)
    INTEGER                                  :: nso
    REAL(dp)                                 :: proj_s(nso,*), proj_h(nso,*), &
                                                buf1(*), buf2(*), fac
    LOGICAL                                  :: distab

!

    IF(na.EQ.0.OR.nb.EQ.0.OR.nso.EQ.0) RETURN
    !
    ! handle special cases
    IF(na.EQ.1.AND.nb.EQ.1) THEN
       !
       ! hard
       CALL dger(nso,nso,fac*blk(1,1),h_a(1,1),1,h_b(1,1),1,proj_h(1,1),nso)
       !
       ! soft
       CALL dger(nso,nso,fac*blk(1,1),s_a(1,1),1,s_b(1,1),1,proj_s(1,1),nso)
    ELSE
       IF (distab) THEN
          !
          ! hard
          CALL dgemm('N','N',na,nso,nb,fac,blk(1,1),ldb,h_b(1,1),nb,0.0_dp,buf1(1),na)
          CALL dgemm('T','N',nso,nso,na,1.0_dp,h_a(1,1),na,buf1(1),na,0.0_dp,buf2(1),nso)
          CALL daxpy(nso*nso,1.0_dp,buf2(1),1,proj_h(1,1),1)
          !
          ! soft
          CALL daxpy(nso*nso,1.0_dp,buf2(1),1,proj_s(1,1),1)
       ELSE
          !
          ! hard
          CALL dgemm('N','N',na,nso,nb,fac,blk(1,1),ldb,h_b(1,1),nb,0.0_dp,buf1(1),na)
          CALL dgemm('T','N',nso,nso,na,1.0_dp,h_a(1,1),na,buf1(1),na,1.0_dp,proj_h(1,1),nso)
          !
          ! soft
          CALL dgemm('N','N',na,nso,nb,fac,blk(1,1),ldb,s_b(1,1),nb,0.0_dp,buf1(1),na)
          CALL dgemm('T','N',nso,nso,na,1.0_dp,s_a(1,1),na,buf1(1),na,1.0_dp,proj_s(1,1),nso)
       ENDIF
    ENDIF
    !
  END SUBROUTINE proj_blk
 

END MODULE qs_oce_methods
