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
  USE kinds,                           ONLY: dp
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
  USE sparse_matrix_types,             ONLY: add_1d_block_node
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

  SUBROUTINE build_oce_block(oceh,oces,atom_ka,atom_kb,rab,dab,nder,ai_work,ldai)


    TYPE(block_p_type), DIMENSION(:),POINTER   :: oceh, oces
    TYPE(atomic_kind_type), POINTER            :: atom_ka, atom_kb
    INTEGER,INTENT(IN)                         :: ldai, nder
    REAL(dp)                                   :: dab
    REAL(dp), DIMENSION(3)                     :: rab
    REAL(dp), DIMENSION(ldai,ldai,*), &
                       INTENT(INOUT)           :: ai_work


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
    REAL(dp), DIMENSION(:,:,:), POINTER        :: spa_bbt,spa_sbt,spa_tmpt
    REAL(dp), DIMENSION(:,:), POINTER          :: spa_bb,spa_sb,spa_tmp,ovc,ovs

    INTEGER  :: first_col, i, ider, ipgf, iset, isgfb, jset, lpoint, lprj, lsgfb, lshell, &
                m, maxsoa, maxcob, ic,jc,is,ip, m1,ico,iso,iprjc,iprjs,igau,lx,ly,lz,&
                maxder, maxlb, maxlprj, maxnprja, n, nsoatot, ncob, np_car,np_sph, nseta, nsetb, ntotsgfb

    INTEGER :: ithread,nthread
    INTEGER :: omp_get_max_threads,omp_get_thread_num
    LOGICAL :: paw_atom_a,paw_atom_b,calculate_forces

#if defined(__AIX) && defined(__INIT_WITH_NANS)
    REAL(dp) :: nans_dp
    nans_dp = z"7FF0000000000001"
#endif

!   ---------------------------------------------------------------------------


   CALL get_atomic_kind(atomic_kind=atom_ka,orb_basis_set=orb_basis_a,&
                        paw_proj_set=paw_proj_a,paw_atom=paw_atom_a)

   CALL get_atomic_kind(atomic_kind=atom_kb,orb_basis_set=orb_basis_b,&
                         paw_proj_set=paw_proj_b,paw_atom=paw_atom_b)

   IF(.NOT. paw_atom_a) RETURN

   CALL get_paw_proj_set(paw_proj_set=paw_proj_a,cprj_s=cprj_s,maxl=maxlprj, &
                         nprj=nprjla,ncgauprj=np_car,nsgauprj=np_sph,rcprj=rcprja, & 
                         first_prj=fp_cara,first_prjs=fp_spha,&
                         rzetprj=rzetprja,zetprj=zetprja)

   CALL get_gto_basis_set(gto_basis_set=orb_basis_b,nset=nsetb, nsgf=ntotsgfb, &
                          set_radius=set_radius_b,lmax=lb_max,lmin=lb_min, &
                          npgf=npgfb,nsgf_set=nsgfb,pgf_radius=rpgfb, &
                          sphi=sphi_b,zet=zetb,first_sgf=first_sgfb, &
                          maxco=maxcob,maxl=maxlb)

   CALL get_gto_basis_set(gto_basis_set=orb_basis_a,nset=nseta,maxso=maxsoa)

!  *** Add the block ab ***

    NULLIFY(spa_sbt)
    NULLIFY(spa_tmp,ovc,ovs)

    nthread = 1
!$  nthread = omp_get_max_threads()

    maxder = ncoset(nder)
    nsoatot = maxsoa*nseta
    maxnprja = SIZE(zetprja,1)

    IF(nder > 0) THEN
      calculate_forces = .TRUE.
    ENDIF

    CALL reallocate(spa_sbt,1,np_car,1,ntotsgfb,0,nthread-1)
    CALL reallocate(spa_tmpt,1,maxnprja*ncoset(maxlprj),1,maxcob*maxder,0,nthread-1)
    CALL reallocate(ovc,1,np_car,1,maxcob*nsetb*maxder)
    CALL reallocate(ovs,1,np_sph,1,maxcob*nsetb*maxder)

    ithread = 0
!$  ithread = omp_get_thread_num()
    spa_sb => spa_sbt(:,:,ithread)
    spa_tmp => spa_tmpt(:,:,ithread)

    m1 = 0
    DO jset = 1,nsetb
      IF (rcprja + set_radius_b(jset) >= dab) THEN

!     ***integral between proj of iatom and primitives of jatom
!     *** Calculate the primitives overlap ***
        spa_tmp = 0.0_dp
        spa_sb  = 0.0_dp
#if defined(__AIX) && defined(__INIT_WITH_NANS)
        ai_work = nans_dp
#endif
        ncob = npgfb(jset)*ncoset(lb_max(jset))
        isgfb = first_sgfb(1,jset)
        lsgfb = isgfb - 1 + nsgfb(jset)



        DO lprj = 0,maxlprj
          CALL overlap(lprj,lprj,nprjla(lprj),&
                       rzetprja(:,lprj),zetprja(:,lprj),&
                       lb_max(jset),lb_min(jset),npgfb(jset),&
                       rpgfb(:,jset),zetb(:,jset),&
                       rab,dab,spa_tmp,&
                       nder,calculate_forces,ai_work,ldai)
          DO ipgf = 1,nprjla(lprj)
            n = (ipgf -1)*ncoset(lprj)
            lpoint = ncoset(lprj -1)+1+ n
            m = fp_cara(lprj)+(ipgf -1)*nco(lprj)
            DO ider =1,maxder
            is = (ider-1)*maxcob
            DO ip = 1,npgfb(jset)
               ic = (ip-1)*ncoset(lb_max(jset))
               DO jc = ncoset(lb_min(jset)-1)+1,ncoset(lb_max(jset)) 
                 ovc(m:m-1+nco(lprj),jc+ic+is+m1) = spa_tmp(lpoint:lpoint-1+nco(lprj),ic+jc+is)
               ENDDO
            ENDDO
            ENDDO
          ENDDO
        ENDDO

        n = ncoset(lb_max(jset))
        DO ider =1,maxder
        is = (ider-1)*maxcob
        DO ipgf = 1,npgfb(jset)
          DO lshell = lb_min(jset),lb_max(jset)
            DO ic = 1,nco(lshell)
              igau = ic + ncoset(lshell-1) + n*(ipgf-1) + m1 + is
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
                  ovs(:,igau+1+is:igau+nco(lshell)+is) = 0.0_dp
                ENDDO
              ENDIF
            ENDDO 
          END DO
        ENDIF

!       *** Contraction step (integrals and derivatives)
        IF(dab <= 1.0E-8_dp)THEN

          iset = jset
          n = maxsoa*(iset-1) 

          CALL get_paw_proj_set(paw_proj_set=paw_proj_a, &
                                local_oce_sphi_h=local_oce_h, & 
                                local_oce_sphi_s=local_oce_s)
          oceh(1)%block(n+1:n+maxsoa,isgfb:lsgfb) = &
                    local_oce_h(1:maxsoa,isgfb:lsgfb)
          oces(1)%block(n+1:n+maxsoa,isgfb:lsgfb) = &
                    local_oce_s(1:maxsoa,isgfb:lsgfb)

          DO ider =2,maxder
            first_col = (ider - 1)*maxcob + 1 + m1
            CALL dgemm("N","N",np_sph,nsgfb(jset),ncob,&
                      1.0_dp,ovs(1,first_col),SIZE(ovs,1),&
                      sphi_b(1,isgfb),SIZE(sphi_b,1),&
                      0.0_dp,spa_sb(1,isgfb),SIZE(spa_sb,1))

            CALL dgemm("T","N",nsoatot,nsgfb(jset),np_sph,&
                      1.0_dp,cprj_s(1,1),SIZE(cprj_s,1),&
                      spa_sb(1,isgfb),SIZE(spa_sb,1),&
                      0.0_dp,oceh(ider)%block(1,isgfb),SIZE(oceh(ider)%block,1)) 

            CALL dgemm("N","N",np_sph,nsgfb(jset),ncob,&
                      1.0_dp,ovs(1,first_col),SIZE(ovs,1),&
                      local_oce_s(1,isgfb),SIZE(local_oce_s,1),&
                      0.0_dp,spa_sb(1,isgfb),SIZE(spa_sb,1))

            CALL dgemm("T","N",nsoatot,nsgfb(jset),np_sph,&
                      1.0_dp,cprj_s(1,1),SIZE(cprj_s,1),&
                      spa_sb(1,isgfb),SIZE(spa_sb,1),&
                      0.0_dp,oceh(ider)%block(1,isgfb),SIZE(oces(ider)%block,1))
          ENDDO
        ELSE       
          DO ider = 1,maxder
            first_col = (ider - 1)*maxcob + 1 + m1
            CALL dgemm("N","N",np_sph,nsgfb(jset),ncob,&
                      1.0_dp,ovs(1,first_col),SIZE(ovs,1),&
                      sphi_b(1,isgfb),SIZE(sphi_b,1),&
                      0.0_dp,spa_sb(1,isgfb),SIZE(spa_sb,1))

            CALL dgemm("T","N",nsoatot,nsgfb(jset),np_sph,&
                      1.0_dp,cprj_s(1,1),SIZE(cprj_s,1),&
                      spa_sb(1,isgfb),SIZE(spa_sb,1),&
                      0.0_dp,oceh(ider)%block(1,isgfb),SIZE(oceh(ider)%block,1))

            oces(ider)%block(1:nsoatot,isgfb:lsgfb) = &
                      oceh(ider)%block(1:nsoatot,isgfb:lsgfb)

          ENDDO
        END IF 
      END IF ! radius
      m1 = m1 + maxcob
    ENDDO  !jset
    
! stop 'oce'

    DEALLOCATE(spa_sbt)

  END SUBROUTINE build_oce_block

! *****************************************************************************

  SUBROUTINE set_up_oce(qs_env,oce,ikind,iatom,jkind,jatom,rab,&
                        oceh,oces,nder,ai_work,ldai)

!   Purpose: set up the sparse matrix for the coefficients  of one center expansions

!   History: - Creation (5-02-04 MI)

!   ***************************************************************************

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(oce_matrix_type), POINTER           :: oce
    INTEGER, INTENT(IN)                      :: ikind, iatom, jkind, jatom
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rab
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces
    INTEGER, INTENT(IN)                      :: nder, ldai
    REAL(KIND=dp), DIMENSION(ldai, ldai, *)  :: ai_work

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "set_up_oce"

    INTEGER                                  :: i, ibc, katom, kkind, &
                                                kkneighbor, kneighbor, nkind, &
                                                nneighbor
    INTEGER, DIMENSION(3)                    :: cell_c
    LOGICAL                                  :: paw_atom
    REAL(KIND=dp)                            :: dac, dbc
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc
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

       sbc_oce_neighbor_list =>&
             find_neighbor_list(neighbor_list_set=&
             sbc_oce(ibc)%neighbor_list_set,&
             atom=jatom)

       IF (.NOT.ASSOCIATED(sbc_oce(ibc)%neighbor_list_set)) CYCLE

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

           at_a => atomic_kind_set(kkind)
           at_b => atomic_kind_set(ikind)

           CALL build_oce_block(oceh,oces, &
                                at_a,&
                                at_b,&
                                rac,dac,&
                                nder,ai_work,ldai)

           IF(jatom == iatom) THEN
             sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
             CYCLE
           ENDIF
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

           at_a => atomic_kind_set(kkind)
           at_b => atomic_kind_set(jkind)
           CALL build_oce_block(oceh,oces, &
                              at_a,&
                              at_b,&
                              rbc,dbc,&
                              nder,ai_work,ldai)
 
           sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)

       ENDDO  knodeloop ! kneighbor
    ENDDO  ! kkind

  END SUBROUTINE set_up_oce

! *****************************************************************************
END MODULE qs_oce_methods
