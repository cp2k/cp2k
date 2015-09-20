!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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

  USE ai_overlap,                      ONLY: overlap
  USE ao_util,                         ONLY: exp_radius
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
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE sap_kind_types,                  ONLY: clist_type,&
                                             sap_int_type,&
                                             sap_sort
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_oce_methods'

! *** Public subroutines ***

  PUBLIC :: build_oce_matrices, proj_blk, prj_scatter, prj_gather

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param oces ...
!> \param atom_ka ...
!> \param atom_kb ...
!> \param rab ...
!> \param nder ...
!> \param sgf_list ...
!> \param nsgf_cnt ...
!> \param sgf_soft_only ...
!> \param eps_fit ...
! *****************************************************************************
  SUBROUTINE build_oce_block(oces,atom_ka,atom_kb,rab,nder,sgf_list,nsgf_cnt,sgf_soft_only,&
       eps_fit)


    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oces
    TYPE(qs_kind_type), POINTER              :: atom_ka, atom_kb
    REAL(dp), DIMENSION(3)                   :: rab
    INTEGER, INTENT(IN)                      :: nder
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sgf_list
    INTEGER, INTENT(OUT)                     :: nsgf_cnt
    LOGICAL, INTENT(OUT)                     :: sgf_soft_only
    REAL(dp), INTENT(IN)                     :: eps_fit

    CHARACTER(len=*), PARAMETER :: routineN = 'build_oce_block', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_col, ic, ico, ider, ig1, igau, ip, ipgf, is, isgfb, &
      isgfb_cnt, iso, isp, istat, jc, jset, lds, lm, lpoint, lprj, lsgfb, &
      lsgfb_cnt, lshell, lx, ly, lz, m, m1, maxcob, maxder, maxlb, maxlprj, &
      maxnprja, maxsoa, msab, n, ncob, np_car, np_sph, nsatbas, nseta, nsetb, &
      nsoatot, ntotsgfb, sgf_hard_only
    INTEGER, DIMENSION(:), POINTER           :: fp_cara, fp_spha, lb_max, &
                                                lb_min, npgfb, nprjla, nsgfb
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfb
    LOGICAL                                  :: calculate_forces, failure, &
                                                paw_atom_a, paw_atom_b
    REAL(dp)                                 :: dab, hard_radius_a, &
                                                hard_radius_b, radius, rcprja
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: c2s, ovs, spa_sb
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: s
    REAL(dp), DIMENSION(:), POINTER          :: set_radius_b, zisominb
    REAL(dp), DIMENSION(:, :), POINTER       :: cprj_s, csprj, rpgfb, &
                                                rzetprja, spa_tmp, sphi_b, &
                                                zetb, zetprja
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a, orb_basis_b
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_a, paw_proj_b

   failure = .FALSE.
   NULLIFY(orb_basis_a,paw_proj_a)
   CALL get_qs_kind(qs_kind=atom_ka,basis_set=orb_basis_a,&
                        paw_proj_set=paw_proj_a,paw_atom=paw_atom_a,&
                        hard_radius=hard_radius_a)

   NULLIFY(orb_basis_b,paw_proj_b)
   CALL get_qs_kind(qs_kind=atom_kb,basis_set=orb_basis_b,&
                         paw_proj_set=paw_proj_b,paw_atom=paw_atom_b,&
                         hard_radius=hard_radius_b)

   IF(.NOT. paw_atom_a) RETURN

   NULLIFY(cprj_s,nprjla,fp_cara,fp_spha,rzetprja,zetprja)
   CALL get_paw_proj_set(paw_proj_set=paw_proj_a,cprj_s=cprj_s,csprj=csprj,maxl=maxlprj, &
                         nprj=nprjla,ncgauprj=np_car,nsgauprj=np_sph,nsatbas=nsatbas,rcprj=rcprja, &
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
    dab = SQRT(SUM(rab*rab))

    maxder = ncoset(nder)
    nsoatot = maxsoa*nseta
    maxnprja = SIZE(zetprja,1)

    calculate_forces = .FALSE.
    IF(nder > 0) THEN
      calculate_forces = .TRUE.
    ENDIF

    lm = MAX(maxlb, maxlprj)
    lds = ncoset(lm+nder+1)
    msab = MAX(maxnprja*ncoset(maxlprj),maxcob)

    ALLOCATE(c2s(lds,lds),STAT= istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(s(lds,lds,ncoset(nder+1)),STAT= istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(spa_sb(np_car,ntotsgfb), STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(spa_tmp(msab,msab*maxder),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(ovs(np_sph,maxcob*nsetb*maxder),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)

    m1 = 0
    nsgf_cnt = 0
    isgfb_cnt = 1
    sgf_hard_only = 0
    DO jset = 1,nsetb
       !
       ! Set the contribution list
       IF(hard_radius_a + set_radius_b(jset) >= dab) THEN
         isgfb = first_sgfb(1,jset)
         lsgfb = isgfb - 1 + nsgfb(jset)
         DO jc = isgfb,lsgfb
           nsgf_cnt = nsgf_cnt + 1
           sgf_list(nsgf_cnt) = jc
         ENDDO

         ! check if this function is hard
         radius = exp_radius(lb_max(jset),MAXVAL(zetb(1:npgfb(jset),jset)),eps_fit,1.0_dp)
         IF(radius.LE.hard_radius_b) sgf_hard_only = sgf_hard_only + 1

         ! ***integral between proj of iatom and primitives of jatom
         ! *** Calculate the primitives overlap ***
         spa_tmp = 0.0_dp
         ovs = 0.0_dp
         s = 0.0_dp
         ncob = npgfb(jset)*ncoset(lb_max(jset))
         isgfb = first_sgfb(1,jset)
         lsgfb = isgfb - 1 + nsgfb(jset)

         lsgfb_cnt = isgfb_cnt - 1 + nsgfb(jset)

         DO lprj = 0,maxlprj
           CALL overlap(lprj,lprj,nprjla(lprj),&
                       rzetprja(:,lprj),zetprja(:,lprj),&
                       lb_max(jset),lb_min(jset),npgfb(jset),&
                       rpgfb(:,jset),zetb(:,jset),&
                       -rab,dab,spa_tmp,&
                       nder,.TRUE.,s,lds)
           DO iso = 1,nso(lprj)
             DO ico = 1,nco(lprj)
               lx = indco(1,ico+ncoset(lprj-1))
               ly = indco(2,ico+ncoset(lprj-1))
               lz = indco(3,ico+ncoset(lprj-1))
               c2s(iso,ico) = orbtramat(lprj)%c2s(iso,ico)/SQRT((4.0_dp*pi)/dfac(2*lprj+1)*&
                              dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1))
             ENDDO
           ENDDO
           DO ider=1,maxder
             is = (ider-1)*SIZE(spa_tmp,1)
             isp =  (ider-1)*maxcob*nsetb
             DO ipgf = 1,nprjla(lprj)
               lpoint = ncoset(lprj -1)+1+ (ipgf -1)*ncoset(lprj)
               m = fp_spha(lprj)+(ipgf -1)*nso(lprj)
               DO ip = 1,npgfb(jset)
                 ic = (ip-1)*ncoset(lb_max(jset))
                 igau = isp+ic+m1+ncoset(lb_min(jset)-1)+1
                 ig1 = is+ic+ncoset(lb_min(jset)-1)+1
                 n = ncoset(lb_max(jset)) - ncoset(lb_min(jset)-1)
                 CALL DGEMM("N","N",nso(lprj),n,nco(lprj),1._dp,c2s,lds,spa_tmp(lpoint,ig1),msab,&
                            0._dp,ovs(m,igau),np_sph)
               ENDDO
             ENDDO
           ENDDO
         ENDDO

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

         ! *** Contraction step (integrals and derivatives)
         DO ider = 1,maxder
          first_col = (ider - 1)*maxcob*nsetb + 1 + m1
          CALL dgemm("N","N",np_sph,nsgfb(jset),ncob,&
                    1.0_dp,ovs(1,first_col),SIZE(ovs,1),&
                    sphi_b(1,isgfb),SIZE(sphi_b,1),&
                    0.0_dp,spa_sb(1,isgfb),SIZE(spa_sb,1))

          CALL dgemm("T","N",nsatbas,nsgfb(jset),np_sph,&
                    1.0_dp,csprj(1,1),SIZE(csprj,1),&
                    spa_sb(1,isgfb),SIZE(spa_sb,1),&
                    1.0_dp,oces(ider)%block(1,isgfb_cnt),SIZE(oces(ider)%block,1))
        ENDDO
        isgfb_cnt = isgfb_cnt + nsgfb(jset)
      END IF ! radius
      m1 = m1 + maxcob
    ENDDO  !jset

    ! check if the screened functions are all soft
    sgf_soft_only = .FALSE.
    IF(sgf_hard_only.EQ.0) sgf_soft_only = .TRUE.

    DEALLOCATE(c2s,s,spa_sb,spa_tmp,ovs,STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)

  END SUBROUTINE build_oce_block

! *****************************************************************************
!> \brief ...
!> \param oceh ...
!> \param oces ...
!> \param atom_ka ...
!> \param sgf_list ...
!> \param nsgf_cnt ...
! *****************************************************************************
  SUBROUTINE build_oce_block_local(oceh,oces,atom_ka,sgf_list,nsgf_cnt)

    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces
    TYPE(qs_kind_type), POINTER              :: atom_ka
    INTEGER, DIMENSION(:), INTENT(OUT)       :: sgf_list
    INTEGER, INTENT(OUT)                     :: nsgf_cnt

    CHARACTER(len=*), PARAMETER :: routineN = 'build_oce_block_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iset, isgfa, j, jc, lsgfa, &
                                                maxlprj, maxsoa, n, nsatbas, &
                                                nseta, nsgfa, stat
    INTEGER, DIMENSION(:), POINTER           :: n2oindex, nsgf_seta
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, paw_atom_a
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: prjloc_h, prjloc_s
    REAL(dp), DIMENSION(:, :), POINTER       :: local_oce_h, local_oce_s
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_a
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_a

    failure = .FALSE.
    NULLIFY(orb_basis_a,paw_proj_a)
    CALL get_qs_kind(qs_kind=atom_ka,basis_set=orb_basis_a,&
                         paw_proj_set=paw_proj_a,paw_atom=paw_atom_a)

    IF(.NOT. paw_atom_a) RETURN

    CALL get_paw_proj_set(paw_proj_set=paw_proj_a,maxl=maxlprj,nsatbas=nsatbas,&
                          n2oindex=n2oindex)

    CALL get_gto_basis_set(gto_basis_set=orb_basis_a,first_sgf=first_sgfa,&
                           nsgf=nsgfa,nsgf_set=nsgf_seta,nset=nseta,maxso=maxsoa)

    NULLIFY(local_oce_h,local_oce_s)
    CALL get_paw_proj_set(paw_proj_set=paw_proj_a, &
                          local_oce_sphi_h=local_oce_h, &
                          local_oce_sphi_s=local_oce_s)

    ALLOCATE(prjloc_h(nseta*maxsoa,nsgfa),prjloc_s(nseta*maxsoa,nsgfa),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)
    prjloc_h = 0._dp
    prjloc_s = 0._dp

    nsgf_cnt = 0
    DO iset = 1,nseta
      isgfa = first_sgfa(1,iset)
      lsgfa = isgfa - 1 + nsgf_seta(iset)
      DO jc = isgfa,lsgfa
        nsgf_cnt = nsgf_cnt + 1
        sgf_list(nsgf_cnt) = jc
      ENDDO
      n = maxsoa*(iset-1)

      prjloc_h(n+1:n+maxsoa,isgfa:lsgfa) = local_oce_h(1:maxsoa,isgfa:lsgfa)
      prjloc_s(n+1:n+maxsoa,isgfa:lsgfa) = local_oce_s(1:maxsoa,isgfa:lsgfa)
    ENDDO

    DO i=1,nsgfa
      DO j=1,nsatbas
        jc=n2oindex(j)
        oceh(1)%block(j,i) = prjloc_h(jc,i)
        oces(1)%block(j,i) = prjloc_s(jc,i)
      END DO
    END DO

    DEALLOCATE(prjloc_h,prjloc_s,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,failure)

  END SUBROUTINE build_oce_block_local


! *****************************************************************************
!> \brief Set up the sparse matrix for the coefficients of one center expansions
!>      This routine uses the same logic as the nonlocal pseudopotential
!> \param intac TYPE that holds the integrals (a=basis; c=projector)
!> \param calculate_forces ...
!> \param nder ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param sap_oce ...
!> \param eps_fit ...
!> \par History
!>      02.2009 created
!> \author jgh
! *****************************************************************************
  SUBROUTINE build_oce_matrices(intac, calculate_forces, nder,&
                    qs_kind_set, particle_set, sap_oce, eps_fit)

    TYPE(sap_int_type), DIMENSION(:), &
      POINTER                                :: intac
    LOGICAL, INTENT(IN)                      :: calculate_forces
    INTEGER                                  :: nder
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sap_oce
    REAL(dp), INTENT(IN)                     :: eps_fit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'build_oce_matrices', &
      routineP = moduleN//':'//routineN

    INTEGER :: atom_a, atom_b, handle, i, iac, ikind, ilist, jkind, &
      jneighbor, ldai, ldsab, maxco, maxder, maxl, maxlgto, maxlprj, maxprj, &
      maxsgf, maxsoa, maxsob, mlprj, natom, ncoa_sum, nkind, nlist, &
      nneighbor, nsatbas, nseta, nsetb, nsgf_cnt, nsgfa, nsobtot, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: sgf_list
    INTEGER, DIMENSION(3)                    :: cell_b
    INTEGER, DIMENSION(:), POINTER           :: fp_car, fp_sph, la_max, &
                                                la_min, npgfa, nprjla, &
                                                nsgf_seta
    INTEGER, DIMENSION(:, :), POINTER        :: first_sgfa
    LOGICAL                                  :: failure, local, paw_atom_b, &
                                                sgf_soft_only
    REAL(KIND=dp)                            :: dab, rcprj
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: sab, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: ai_work
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:), POINTER     :: set_radius_a
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rpgfa, rzetprj, sphi_a, zeta, &
                                                zetb
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces
    TYPE(clist_type), POINTER                :: clist
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_paw, orb_basis_set
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(paw_proj_set_type), POINTER         :: paw_proj_b
    TYPE(qs_kind_type), POINTER              :: at_a, at_b, qs_kind

    failure = .FALSE.
    IF (calculate_forces) THEN
      CALL timeset(routineN//"_forces",handle)
    ELSE
      CALL timeset(routineN,handle)
    ENDIF

    IF ( ASSOCIATED(sap_oce) ) THEN

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      maxder = ncoset(nder)

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set,&
                               maxco=maxco,&
                               maxlgto=maxlgto,&
                               maxlprj=maxlprj,&
                               maxco_proj=maxprj,&
                               maxsgf=maxsgf)

      maxl = MAX(maxlgto,maxlprj)
      CALL init_orbital_pointers(maxl+nder+1)

      ldsab = MAX(maxco,ncoset(maxlprj),maxsgf,maxprj)
      ldai = ncoset(maxl+nder+1)
      ALLOCATE(sab(ldsab,ldsab*maxder),work(ldsab,ldsab*maxder),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,failure)
      sab=0.0_dp
      ALLOCATE (ai_work(ldai,ldai,ncoset(nder+1)),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,failure)
      ai_work=0.0_dp

      ALLOCATE(oceh(maxder),oces(maxder),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,failure)

      DO i=1,nkind*nkind
        NULLIFY(intac(i)%alist,intac(i)%asort,intac(i)%aindex)
        intac(i)%nalist=0
      END DO

      !calculate the overlap integrals <a|p>
      CALL neighbor_list_iterator_create(nl_iterator,sap_oce)
      DO WHILE (neighbor_list_iterate(nl_iterator)==0)
         CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,iatom=atom_a,jatom=atom_b,&
              nlist=nlist,ilist=ilist,nnode=nneighbor,inode=jneighbor,cell=cell_b,r=rab)

         iac = ikind + nkind*(jkind - 1)
         dab = SQRT(SUM(rab*rab))

         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind,basis_set=orb_basis_set)

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

         qs_kind => qs_kind_set(jkind)

         NULLIFY(paw_proj_b)
         CALL get_qs_kind(qs_kind=qs_kind,paw_proj_set=paw_proj_b,paw_atom=paw_atom_b)
         IF (.NOT.paw_atom_b) CYCLE

         CALL get_qs_kind(qs_kind=qs_kind,basis_set=orb_basis_paw)
         IF (.NOT.ASSOCIATED(orb_basis_paw)) CYCLE
         CALL get_gto_basis_set(gto_basis_set=orb_basis_paw,maxso=maxsob,nset=nsetb)

         CALL get_paw_proj_set(paw_proj_set=paw_proj_b,&
                               maxl=mlprj,&
                               nprj=nprjla,&
                               nsatbas=nsatbas,&
                               rcprj=rcprj, &
                               first_prj=fp_car,&
                               first_prjs=fp_sph,&
                               rzetprj=rzetprj,&
                               zetprj=zetb)

         nsobtot = nsatbas

         IF(.NOT.ASSOCIATED(intac(iac)%alist)) THEN
            intac(iac)%a_kind = ikind
            intac(iac)%p_kind = jkind
            intac(iac)%nalist = nlist
            ALLOCATE(intac(iac)%alist(nlist),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,failure)
            DO i=1,nlist
               NULLIFY(intac(iac)%alist(i)%clist)
               intac(iac)%alist(i)%aatom  = 0
               intac(iac)%alist(i)%nclist = 0
            END DO
         END IF
         IF(.NOT.ASSOCIATED(intac(iac)%alist(ilist)%clist)) THEN
            intac(iac)%alist(ilist)%aatom  = atom_a
            intac(iac)%alist(ilist)%nclist = nneighbor
            ALLOCATE(intac(iac)%alist(ilist)%clist(nneighbor),STAT=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,failure)
         END IF

         clist => intac(iac)%alist(ilist)%clist(jneighbor)
         clist%catom = atom_b
         clist%cell  = cell_b
         clist%rac   = rab
         clist%nsgf_cnt = 0
         clist%maxac    = 0.0_dp
         clist%maxach   = 0.0_dp
         NULLIFY(clist%acint,clist%achint,clist%sgf_list)

         ALLOCATE(sgf_list(nsgfa),STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,failure)

         at_a => qs_kind_set(jkind)
         at_b => qs_kind_set(ikind)

         local = (atom_a == atom_b .AND. ALL(cell_b == 0))

         IF( local ) THEN
            DO i=1,maxder
               ALLOCATE(oceh(i)%block(nsobtot,nsgfa),oces(i)%block(nsobtot,nsgfa),STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,failure)
               oceh(i)%block=0._dp
               oces(i)%block=0._dp
            END DO
            CALL build_oce_block_local(oceh,oces,at_a,sgf_list,nsgf_cnt)
            clist%nsgf_cnt = nsgf_cnt
            clist%sgf_soft_only = .FALSE.
            IF(nsgf_cnt > 0) THEN
               ALLOCATE(clist%acint(nsgf_cnt,nsobtot,maxder),clist%sgf_list(nsgf_cnt),STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,failure)
               clist%acint(:,:,:)=0._dp
               clist%sgf_list(:)=HUGE(0)
               CPPostcondition(nsgf_cnt == nsgfa,cp_failure_level,routineP,failure)
               ! *** Special case: A=B
               ALLOCATE(clist%achint(nsgfa,nsobtot,maxder),STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,failure)
               clist%achint=0._dp
               clist%acint(1:nsgfa,1:nsobtot,1) = TRANSPOSE(oces(1)%block(1:nsobtot,1:nsgfa))
               clist%achint(1:nsgfa,1:nsobtot,1) = TRANSPOSE(oceh(1)%block(1:nsobtot,1:nsgfa))
               clist%maxac=MAXVAL(ABS(clist%acint(:,:,1)))
               clist%maxach=0._dp
               clist%sgf_list(1:nsgf_cnt) = sgf_list(1:nsgf_cnt)
            ENDIF
            DO i=1,maxder
              DEALLOCATE(oceh(i)%block,oces(i)%block,STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,failure)
            END DO
         ELSE
            DO i=1,maxder
               ALLOCATE(oces(i)%block(nsobtot,nsgfa),STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,failure)
               oces(i)%block=0._dp
            END DO
            CALL build_oce_block(oces,at_a,at_b,rab,nder,sgf_list,nsgf_cnt,sgf_soft_only,eps_fit)
            clist%nsgf_cnt = nsgf_cnt
            clist%sgf_soft_only = sgf_soft_only
            IF(nsgf_cnt > 0) THEN
               ALLOCATE(clist%acint(nsgf_cnt,nsobtot,maxder),clist%sgf_list(nsgf_cnt),STAT=stat)
               CPPostcondition(stat==0,cp_failure_level,routineP,failure)
               clist%acint(:,:,:)=0._dp
               clist%sgf_list(:)=HUGE(0)
               DO i=1,maxder
                  clist%acint(1:nsgf_cnt,1:nsobtot,i) = TRANSPOSE(oces(i)%block(1:nsobtot,1:nsgf_cnt))
               END DO
               clist%maxac=MAXVAL(ABS(clist%acint(:,:,1)))
               clist%maxach=0._dp
               clist%sgf_list(1:nsgf_cnt) = sgf_list(1:nsgf_cnt)
            ENDIF
            DO i=1,maxder
              DEALLOCATE(oces(i)%block,STAT=stat)
              CPPostcondition(stat==0,cp_failure_level,routineP,failure)
            END DO
         END IF

         DEALLOCATE(sgf_list,STAT=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,failure)

     END DO
     CALL neighbor_list_iterator_release(nl_iterator)

     DEALLOCATE(sab,work,ai_work,STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,failure)
     DEALLOCATE(oceh,oces,STAT=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,failure)

     ! set up sort index
     CALL sap_sort(intac)

   END IF

   CALL timestop(handle)

 END SUBROUTINE build_oce_matrices

! *****************************************************************************
!> \brief Project a matrix block onto the local atomic functions.
!>
!> \param h_a ...
!> \param s_a ...
!> \param na ...
!> \param h_b ...
!> \param s_b ...
!> \param nb ...
!> \param blk ...
!> \param ldb ...
!> \param proj_h ...
!> \param proj_s ...
!> \param nso ...
!> \param buf1 ...
!> \param buf2 ...
!> \param fac ...
!> \param distab ...
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

! *****************************************************************************
!> \brief ...
!> \param ain matrix in old indexing
!> \param aout matrix in new compressed indexing
!> \param atom ...
! *****************************************************************************
  SUBROUTINE prj_gather(ain,aout,atom)

    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: ain
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: aout
    TYPE(qs_kind_type), INTENT(IN)           :: atom

    CHARACTER(len=*), PARAMETER :: routineN = 'prj_gather', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ip, j, jp, nbas
    INTEGER, DIMENSION(:), POINTER           :: n2oindex
    LOGICAL                                  :: paw_atom
    TYPE(paw_proj_set_type), POINTER         :: paw_proj

    NULLIFY(paw_proj)
    CALL get_qs_kind(qs_kind=atom,paw_proj_set=paw_proj,paw_atom=paw_atom)
    NULLIFY(n2oindex)
    CALL get_paw_proj_set(paw_proj_set=paw_proj,nsatbas=nbas,n2oindex=n2oindex)

    DO i=1,nbas
      ip = n2oindex(i)
      DO j=1,nbas
        jp = n2oindex(j)
        aout(j,i) = ain(jp,ip)
      END DO
    END DO

  END SUBROUTINE prj_gather

! *****************************************************************************
!> \brief ...
!> \param ain  matrix in new compressed indexing
!> \param aout matrix in old indexing (addup)
!> \param atom ...
! *****************************************************************************
  SUBROUTINE prj_scatter(ain,aout,atom)

    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: ain
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: aout
    TYPE(qs_kind_type), INTENT(IN)           :: atom

    CHARACTER(len=*), PARAMETER :: routineN = 'prj_scatter', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ip, j, jp, nbas
    INTEGER, DIMENSION(:), POINTER           :: n2oindex
    LOGICAL                                  :: paw_atom
    TYPE(paw_proj_set_type), POINTER         :: paw_proj

    NULLIFY(paw_proj)
    CALL get_qs_kind(qs_kind=atom,paw_proj_set=paw_proj,paw_atom=paw_atom)
    NULLIFY(n2oindex)
    CALL get_paw_proj_set(paw_proj_set=paw_proj,nsatbas=nbas,n2oindex=n2oindex)

    DO i=1,nbas
      ip = n2oindex(i)
      DO j=1,nbas
        jp = n2oindex(j)
        aout(jp,ip) = aout(jp,ip) + ain(j,i)
      END DO
    END DO

  END SUBROUTINE prj_scatter

END MODULE qs_oce_methods
