!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_oce_methods
!!
!!   NAME
!!     qs_oce_methods
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE qs_oce_methods

  USE ai_overlap_new,                  ONLY: overlap
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE kinds,                           ONLY: dp,dp_size
  USE mathconstants,                   ONLY: dfac,&
                                             pi
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: indco,&
                                             nco,&
                                             ncoset,&
                                             nso
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_neighbor_list_types,          ONLY: block_p_type,&
                                             find_neighbor_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next,&
                                             qlist_type
  USE qs_oce_types,                    ONLY: oce_matrix_type
  USE sparse_matrix_types,             ONLY: add_1d_block_node,&
                                             get_block_node
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

  USE util,                            ONLY: locate

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_oce_methods"

! *** Public subroutines ***

  PUBLIC :: set_up_oce

!!***
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE build_oce_block(oceh,oces,atom_ka,atom_kb,rab,dab,nder)

    TYPE(block_p_type), DIMENSION(:),POINTER   :: oceh, oces
    TYPE(atomic_kind_type), POINTER            :: atom_ka, atom_kb
    INTEGER,INTENT(IN)                         :: nder
    REAL(dp)                                   :: dab
    REAL(dp), DIMENSION(3)                     :: rab


!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine_name = "build_oce_block"

!   *** Local variables ***

    TYPE(gto_basis_set_type), POINTER          :: orb_basis_a, orb_basis_b

    TYPE(paw_proj_set_type),  POINTER          :: paw_proj_a, paw_proj_b

    INTEGER, DIMENSION(:), POINTER             :: lb_max, lb_min, nprjla, &
                                                  nsgfb, npgfb
    INTEGER, DIMENSION(:), POINTER             :: fp_cara,&
                                                  fp_spha
    INTEGER, DIMENSION(:,:), POINTER           :: first_sgfb
    REAL(dp)                                   :: rcprja
    REAL(dp), DIMENSION(:), POINTER            :: set_radius_b, zisominb
    REAL(dp), DIMENSION(:,:), POINTER          :: cprj_s, local_oce_s, local_oce_h, &
                                                  rpgfb, rzetprja, &
                                                  sphi_b, zetb, zetprja
    REAL(dp), DIMENSION(:,:), POINTER          :: spa_sb,spa_tmp,ovc,ovs
    REAL(dp), DIMENSION(:,:,:), POINTER        :: spa_sbt,spa_tmpt,sdpa_tmp
    REAL(dp), DIMENSION(:,:,:,:), POINTER      :: sdpa_tmpt
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    :: s

    INTEGER  :: first_col, i, ic, ico, ider, igau, ip, ipgf, iprjc, iprjs, is, iset, &
                isgfb, iso, isp, istat, jc, jset, lds, lm, lpoint, lprj, lsgfb, lshell, &
                lx,ly,lz, m, m1, maxcob,maxder, maxlb, maxlprj, maxnprja,maxsoa, &
                n, ncob, np_car,np_sph, nseta, nsetb, nsoatot, ntotsgfb
    INTEGER :: ithread,nthread, omp_get_max_threads,omp_get_thread_num
    LOGICAL :: paw_atom_a,paw_atom_b,calculate_forces

#if defined(__AIX) && defined(__INIT_WITH_NANS)
    REAL(dp) :: nans_dp
    nans_dp = z"7FF0000000000001"
#endif

!   ---------------------------------------------------------------------------

   NULLIFY(orb_basis_a,paw_proj_a)
   CALL get_atomic_kind(atomic_kind=atom_ka,orb_basis_set=orb_basis_a,&
                        paw_proj_set=paw_proj_a,paw_atom=paw_atom_a)

   NULLIFY(orb_basis_b,paw_proj_b)
   CALL get_atomic_kind(atomic_kind=atom_kb,orb_basis_set=orb_basis_b,&
                         paw_proj_set=paw_proj_b,paw_atom=paw_atom_b)

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

    IF(nder > 0) THEN
      calculate_forces = .TRUE.
    ENDIF

    lm = MAX(maxlb, maxlprj)
    lds = ncoset(lm+nder+1)
    ALLOCATE(s(lds,lds,ncoset(nder+1)),STAT= istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "s",lds*ncoset(nder+1)*dp_size)

    CALL reallocate(spa_sbt,1,np_car,1,ntotsgfb,0,nthread-1)
    CALL reallocate(spa_tmpt,1,maxnprja*ncoset(maxlprj),1,maxcob*maxder,0,nthread-1)
    CALL reallocate(sdpa_tmpt,1,maxnprja*ncoset(maxlprj+1),1,maxcob,1,4,0,nthread-1)
    CALL reallocate(ovc,1,np_car,1,maxcob*nsetb*maxder)
    CALL reallocate(ovs,1,np_sph,1,maxcob*nsetb*maxder)

    ithread = 0
!$  ithread = omp_get_thread_num()
    spa_sb => spa_sbt(:,:,ithread)
    spa_tmp => spa_tmpt(:,:,ithread)
    sdpa_tmp => sdpa_tmpt(:,:,:,ithread)

    m1 = 0
    DO jset = 1,nsetb
      IF (rcprja + set_radius_b(jset) >= dab) THEN

!     ***integral between proj of iatom and primitives of jatom
!     *** Calculate the primitives overlap ***
        spa_tmp = 0.0_dp
        spa_sb  = 0.0_dp
        sdpa_tmp = 0.0_dp

        s = 0.0_dp
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        isgfb = first_sgfb(1,jset)
        lsgfb = isgfb - 1 + nsgfb(jset)

        DO lprj = 0,maxlprj
          IF (calculate_forces) THEN
            CALL overlap(lprj,lprj,nprjla(lprj),&
                       rzetprja(:,lprj),zetprja(:,lprj),&
                       lb_max(jset),lb_min(jset),npgfb(jset),&
                       rpgfb(:,jset),zetb(:,jset),&
                       -rab,dab,spa_tmp,&
                       0,.false.,s,lds,sdpa_tmp)

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
                       nder,.false.,s,lds)

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
                      1.0_dp,oceh(ider)%block(1,isgfb),SIZE(oceh(ider)%block,1))

            oces(ider)%block(1:nsoatot,isgfb:lsgfb) = &
                      oceh(ider)%block(1:nsoatot,isgfb:lsgfb)
          ENDDO
        END IF 
      END IF ! radius
      m1 = m1 + maxcob
    ENDDO  !jset
    
    DEALLOCATE(s,spa_sbt,spa_tmpt,sdpa_tmpt,ovc,ovs ,STAT = istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "s,spa_sbt,spa_tmpt,sdpa_tmpt,ovc,ovs")

  END SUBROUTINE build_oce_block

! *****************************************************************************

  SUBROUTINE set_up_oce(qs_env,oce,ikind,iatom,jkind,jatom,rab,&
                        oceh,oces,nder)

!   Purpose: set up the sparse matrix for the coefficients  of one center expansions

!   History: - Creation (5-02-04 MI)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(oce_matrix_type), POINTER           :: oce
    INTEGER, INTENT(IN)                      :: ikind, iatom, jkind, jatom
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rab
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces
    INTEGER, INTENT(IN)                      :: nder

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "set_up_oce"

    INTEGER                                  :: handle, i, ibc, j, katom, kkind, &
                                                kkneighbor, kneighbor, ndist, nkind, &
                                                nneighbor
    INTEGER, DIMENSION(3)                    :: cell_c
    LOGICAL                                  :: add_oce,dist_done,paw_atom
    REAL(KIND=dp)                            :: dac, dbc
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc, rdist
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: at_a, at_b, atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sac_oce, sbc_oce
    TYPE(neighbor_list_type), POINTER        :: sbc_oce_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sbc_oce_neighbor_node
    TYPE(qlist_type), DIMENSION(:, :, :), &
      POINTER                                :: tmp_neigh
    TYPE(qlist_type), POINTER                :: sac_oce_neighbor

!   ***************************************************************************

    CALL timeset(routine_name,"I"," ",handle)

    NULLIFY(atomic_kind_set)
    NULLIFY(sac_oce)
    NULLIFY(sbc_oce)

    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    sac_oce=sac_oce,&
                    sbc_oce=sbc_oce)
    nkind = SIZE(atomic_kind_set)

    DO kkind=1,nkind

       IF(.NOT.ASSOCIATED(oce%voce(kkind,iatom)%neighbor)) CYCLE

       ibc = jkind + nkind*(kkind - 1)

       IF (.NOT.ASSOCIATED(sbc_oce(ibc)%neighbor_list_set)) CYCLE
       NULLIFY(sbc_oce_neighbor_list)

       sbc_oce_neighbor_list =>&
             find_neighbor_list(neighbor_list_set=&
             sbc_oce(ibc)%neighbor_list_set,&
             atom=jatom)
       IF (.NOT.ASSOCIATED(sbc_oce_neighbor_list)) CYCLE

       atomic_kind => atomic_kind_set(kkind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            orb_basis_set=orb_basis_set,&
                            paw_atom=paw_atom)

       IF(.NOT. paw_atom)  CYCLE

       sbc_oce_neighbor_list =>&
         find_neighbor_list(neighbor_list_set=&
                         sbc_oce(ibc)%neighbor_list_set,atom=jatom)

       CALL get_neighbor_list(neighbor_list=sbc_oce_neighbor_list,&
                                     nnode=nneighbor)

       sbc_oce_neighbor_node => first_node(sbc_oce_neighbor_list)

       knodeloop: DO kneighbor=1,nneighbor
          CALL get_neighbor_node(neighbor_node=sbc_oce_neighbor_node,&
                            neighbor=katom,&
                            cell=cell_c,&
                            r=rbc)
          dbc = SQRT(rbc(1)*rbc(1)+ rbc(2)*rbc(2) + rbc(3)*rbc(3))

          tmp_neigh => oce%voce(kkind,iatom)%neighbor
          DO i=1, 3
            IF ( (cell_c(i) < LBOUND(tmp_neigh,i)) .OR. &
                 (cell_c(i) > UBOUND(tmp_neigh,i))) THEN
              sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
              CYCLE knodeloop
            END IF
          END DO

          sac_oce_neighbor => oce%voce(kkind,iatom)%neighbor(cell_c(1),&
                                                            cell_c(2),&
                                                            cell_c(3))
          IF (sac_oce_neighbor%n == 0) THEN
              sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
              CYCLE
          END IF

!         *** Locate operator atom in the sac_oce neighbor list ***

           kkneighbor = locate(sac_oce_neighbor%list,katom)

           IF (kkneighbor == 0) THEN
               sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
               CYCLE
           END IF

           rac(:) = sac_oce_neighbor%r(:,kkneighbor)
           dac = SQRT(rac(1)*rac(1)+ rac(2)*rac(2) + rac(3)*rac(3))

!          *** Check whether this block is already done
           add_oce = .FALSE.

           IF(oce%couple(katom,iatom)%block_created) THEN
             dist_done = .FALSE.
             DO i = 1,oce%couple(katom,iatom)%ndist
               rdist(1:3) = oce%couple(katom,iatom)%rdist_old(1:3,i)
               IF(ABS(rdist(1)-rac(1))< 1.E-8_dp .AND.&
                    ABS(rdist(2)-rac(2))< 1.E-8_dp .AND.&
                    ABS(rdist(3)-rac(3))< 1.E-8_dp ) THEN
                 dist_done = .TRUE.

                 GOTO 10
               END IF
             END DO
             add_oce = .TRUE.
           END IF

           IF (add_oce) THEN
             DO i=1,SIZE(oce%hardo)
                NULLIFY (oceh(i)%block)
                CALL get_block_node(matrix=oce%hardo(i)%matrix,&
                                    block_row=katom,&
                                    block_col=iatom,&
                                    BLOCK=oceh(i)%block)
                NULLIFY (oces(i)%block)
                CALL get_block_node(matrix=oce%softo(i)%matrix,&
                                    block_row=katom,&
                                    block_col=iatom,&
                                    BLOCK=oces(i)%block)
             END DO
           ELSE
             DO i=1,SIZE(oce%hardo)
                NULLIFY (oceh(i)%block)
                CALL add_1d_block_node(matrix=oce%hardo(i)%matrix,&
                                 block_row=katom,&
                                 block_col=iatom,&
                                 BLOCK=oceh(i)%block)
                NULLIFY (oces(i)%block)
                CALL add_1d_block_node(matrix=oce%softo(i)%matrix,&
                                 block_row=katom,&
                                 block_col=iatom,&
                                 BLOCK=oces(i)%block)
             ENDDO
           END IF
           at_a => atomic_kind_set(kkind)
           at_b => atomic_kind_set(ikind)

           CALL build_oce_block(oceh,oces, &
                                at_a,&
                                at_b,&
                                rac,dac,&
                                nder)

!         *** This block has been done ***
          oce%block_added(katom,iatom) = .TRUE.
          oce%couple(katom,iatom)%block_created = .TRUE.
          oce%couple(katom,iatom)%ndist =oce%couple(katom,iatom)%ndist + 1
          ndist = oce%couple(katom,iatom)%ndist
          CALL reallocate(oce%couple(katom,iatom)%rdist_old,1,3,1,ndist)
          oce%couple(katom,iatom)%rdist_old(1:3,ndist) = rac(1:3) 

10         CONTINUE
           IF(jatom == iatom) THEN
             sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
             CYCLE
           ENDIF
!          *** Check whether this block is already done
           add_oce = .FALSE.

           IF(oce%couple(katom,jatom)%block_created) THEN
             dist_done = .FALSE.
             DO i = 1,oce%couple(katom,jatom)%ndist
               rdist(1:3) = oce%couple(katom,jatom)%rdist_old(1:3,i)
               IF(ABS(rdist(1)-rbc(1))< 1.E-8_dp .AND.&
                    ABS(rdist(2)-rbc(2))< 1.E-8_dp .AND.&
                    ABS(rdist(3)-rbc(3))< 1.E-8_dp ) THEN
                 dist_done = .TRUE.

                 GOTO 20
               END IF
             END DO
             add_oce = .TRUE.
           END IF

           IF(add_oce) THEN
             DO i=1,SIZE(oce%hardo)
                NULLIFY (oceh(i)%block)
                CALL get_block_node(matrix=oce%hardo(i)%matrix,&
                                block_row=katom,&
                                block_col=jatom,&
                                BLOCK=oceh(i)%block)
                NULLIFY (oces(i)%block)
                CALL get_block_node(matrix=oce%softo(i)%matrix,&
                                block_row=katom,&
                                block_col=jatom,&
                                BLOCK=oces(i)%block)
             END DO
           ELSE
             DO i=1,SIZE(oce%hardo)
                NULLIFY (oceh(i)%block)
                CALL add_1d_block_node(matrix=oce%hardo(i)%matrix,&
                               block_row=katom,&
                               block_col=jatom,&
                               BLOCK=oceh(i)%block)
                NULLIFY (oces(i)%block)
                CALL add_1d_block_node(matrix=oce%softo(i)%matrix,&
                               block_row=katom,&
                               block_col=jatom,&
                               BLOCK=oces(i)%block)
             ENDDO
           END IF

           at_a => atomic_kind_set(kkind)
           at_b => atomic_kind_set(jkind)


           CALL build_oce_block(oceh,oces, &
                              at_a,&
                              at_b,&
                              rbc,dbc,&
                              nder)
!         *** This block has been done ***
          oce%block_added(katom,jatom) = .TRUE.
          oce%couple(katom,jatom)%block_created = .TRUE.
          oce%couple(katom,jatom)%ndist =oce%couple(katom,jatom)%ndist + 1
          ndist = oce%couple(katom,jatom)%ndist
          CALL reallocate(oce%couple(katom,jatom)%rdist_old,1,3,1,ndist)
          oce%couple(katom,jatom)%rdist_old(1:3,ndist) = rbc(1:3)

20         CONTINUE
           sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)

       ENDDO  knodeloop ! kneighbor
    ENDDO  ! kkind
 
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE set_up_oce

! *****************************************************************************
END MODULE qs_oce_methods
