!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****h* cp2k/qs_oce_methods
!!
!!   NAME
!!     qs_oce_methods
!!
!!   FUNCTION
!!     Routines for the construction of the coefficients
!!     for the expansion  of the atomic
!!     densities rho1_hard and rho1_soft in terms of primitive spherical gaussians.
!!
!!   AUTHOR
!!     MI
!!   
!!   MODIFICATION HISTORY
!!     05-2004 created
!!
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
                                             reduced_3c_list_type
  USE qs_oce_types,                    ONLY: oce_matrix_type,&
                                             dist_list_4oce_type
  USE sparse_matrix_types,             ONLY: add_1d_block_node,&
                                             get_block_node
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

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

    calculate_forces = .FALSE.
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
!!****f*  qs_oce_methods/set_up_oce  [1.0] *
!!
!!   NAME
!!     set_up_oce
!!
!!   FUNCTION
!!     Set up the sparse matrix for the coefficients  of one center expansions
!!     To this end, it uses the same pair lists employed for the construction
!!     of the density matriax. Here for each pair iatom jatom, all the katom
!!     which are close to them both are taken into account.
!!
!!   NOTES
!!     - The OCE blocks are dimensioned:
!!       number of primitives on katom times of contractions on iatom/jatom
!!
!!   INPUTS
!!     - atomic_kind_set: all the information of the atomic kinds
!!     - reduced_3c_oce: special 3 center lists that give the neighbors to iatom and jatom
!!     - oce: sparse matrix of the one center coefficients
!!     - ikind, iatom, jkind, jatom : indexes and kinds of the atoms of the pair
!!     - oceh,oces : blocks of the oce matrix for the hard and soft densities
!!     - nder: it is =/ 0 if the forces are required
!!     - dist_list: tool used to know which triplets have been already taken into account
!!                  (to avoid double counting in presence of PBC)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     MI
!!
!!   MODIFICATION HISTORY
!!     02.2004 created 
!!
!!*** **********************************************************************
  SUBROUTINE set_up_oce(atomic_kind_set,reduced_3c_oce,oce,&
                        ikind,iatom,jkind,jatom,rab,&
                        oceh,oces,nder,dist_list)


    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce
    TYPE(oce_matrix_type), POINTER           :: oce
    INTEGER, INTENT(IN)                      :: ikind, iatom, jkind, jatom
    REAL(dp), DIMENSION(3), INTENT(IN)       :: rab
    TYPE(block_p_type), DIMENSION(:), &
      POINTER                                :: oceh, oces
    INTEGER, INTENT(IN)                      :: nder
    TYPE(dist_list_4oce_type), DIMENSION(:),&
      POINTER                                :: dist_list

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "set_up_oce"

    INTEGER                                  :: handle, i, j, katom, kkind, &
                                                kneighbor, ndist, nkind, &
                                                nneighbor
    LOGICAL                                  :: add_oce,dist_done,paw_atom
    REAL(KIND=dp)                            :: dac, dbc
    REAL(KIND=dp), DIMENSION(3)              :: rac, rbc, rdist
    TYPE(atomic_kind_type), POINTER          :: at_a, at_b, atomic_kind

!   ---------------------------------------------------------------------------

    CALL timeset(routine_name,"I"," ",handle)

    nkind = SIZE(atomic_kind_set)

    DO kkind=1,nkind

       atomic_kind => atomic_kind_set(kkind)

       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            paw_atom=paw_atom)

       IF(.NOT. paw_atom)  CYCLE

        nneighbor = reduced_3c_oce(kkind)%nnode
        DO kneighbor = 1,nneighbor

           katom = reduced_3c_oce(kkind)%index_atom(kneighbor)
           rac(1:3) = reduced_3c_oce(kkind)%rac(1:3,kneighbor)
           dac = SQRT(reduced_3c_oce(kkind)%rac2(kneighbor))
           rbc(1:3) = reduced_3c_oce(kkind)%rbc(1:3,kneighbor)
           dbc = SQRT(reduced_3c_oce(kkind)%rbc2(kneighbor))

!          *** Check whether this block is already done
           add_oce = .FALSE.

           IF(oce%couple(katom,iatom)%block_created) THEN
             dist_done = .FALSE.
             DO i = 1,dist_list(katom)%ndist
               rdist(1:3) = dist_list(katom)%rdist(1:3,i)
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
          oce%couple(katom,iatom)%block_created = .TRUE.
          dist_list(katom)%ndist =dist_list(katom)%ndist + 1
          ndist = dist_list(katom)%ndist
          CALL reallocate(dist_list(katom)%rdist,1,3,1,ndist)
          dist_list(katom)%rdist(1:3,ndist) = rac(1:3) 

10         CONTINUE
           IF(jatom == iatom) THEN
             CYCLE
           ENDIF
!          *** Check whether this block is already done
           add_oce = .FALSE.

           IF(oce%couple(katom,jatom)%block_created) THEN
             dist_done = .FALSE.
             DO i = 1,dist_list(katom)%ndist
               rdist(1:3) = dist_list(katom)%rdist(1:3,i)
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
          oce%couple(katom,jatom)%block_created = .TRUE.
          dist_list(katom)%ndist =dist_list(katom)%ndist + 1
          ndist = dist_list(katom)%ndist
          CALL reallocate(dist_list(katom)%rdist,1,3,1,ndist)
          dist_list(katom)%rdist(1:3,ndist) = rbc(1:3)

20         CONTINUE

       ENDDO   ! kneighbor
    ENDDO  ! kkind
 
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE set_up_oce

! *****************************************************************************
END MODULE qs_oce_methods
