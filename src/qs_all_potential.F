!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 -2004 CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_all_potential [1.0] *
!!
!!   NAME
!!     qs_all_potential
!!
!!   FUNCTION
!!     Calculation of Core Hamiltonian contributions due to all-electron potential
!!     The Coulomb integrals are calculated over Cartesian Gaussian-type functions
!!     (electron repulsion integrals, ERIs).
!!
!!   AUTHOR
!!
!!   LITERATURE
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE qs_all_potential 

  USE ai_verfc,                        ONLY: verfc
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential
  USE kinds,                           ONLY: dp,dp_size,&
                                             int_size
  USE memory_utilities,                ONLY: reallocate
  USE orbital_pointers,                ONLY: coset,indco,ncoset

  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_node_type,&
                                             next,&
                                             qlist_type,&
                                             reduced_3c_list_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE util,                            ONLY: locate

  IMPLICIT NONE

  PRIVATE


! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_all_potential"

! *** Define vall type
  TYPE vall_type
    TYPE(qlist_type), DIMENSION(:,:,:), POINTER :: neighbor
    REAL(KIND = dp)                             :: alpha_c,&
                                                   core_charge,&
                                                   core_radius
    INTEGER                                     :: zeta_c
    REAL(KIND=dp), DIMENSION(:,:), POINTER      :: force_c
  END TYPE vall_type

! Public Types

  PUBLIC :: vall_type

! Public Subroutine
 
  PUBLIC :: all_integrals, radii_3c_orbxall

!******************************************************************************
 
  CONTAINS
 
!******************************************************************************
  SUBROUTINE  all_integrals(la_max,la_min,npgfa,rpgfa,zeta,&
                            lb_max,lb_min,npgfb,rpgfb,zetb,&
                            rab,dab,jkind,jatom,nkind,sbc_3c,&
                            reduced_3c_rho0,vall,hab,nder,pab,f_a,f_b)

    INTEGER, INTENT(IN)                        :: la_max, la_min, npgfa
    REAL(dp), DIMENSION(:), INTENT(IN)         :: rpgfa, zeta
    INTEGER, INTENT(IN)                        :: lb_max, lb_min, npgfb
    REAL(dp), DIMENSION(:), INTENT(IN)         :: rpgfb, zetb
    REAL(dp), DIMENSION(3), INTENT(IN)         :: rab
    REAL(dp), INTENT(IN)                       :: dab
    INTEGER, INTENT(IN)                        :: jkind,jatom,nkind
    TYPE(neighbor_list_set_p_type), DIMENSION(:),&
              POINTER                          :: sbc_3c
    TYPE(reduced_3c_list_type), DIMENSION(:),&
              POINTER                          :: reduced_3c_rho0
    TYPE(vall_type), DIMENSION(:), POINTER     :: vall
    REAL(dp), DIMENSION(:, :), INTENT(INOUT)   :: hab
    INTEGER, INTENT(IN)                        :: nder
    REAL(dp), DIMENSION(:, :), INTENT(IN), &
             OPTIONAL                          :: pab
    REAL(dp), DIMENSION(3), INTENT(OUT), &
             OPTIONAL                          :: f_a, f_b

    CHARACTER(LEN=*), PARAMETER :: routine_name = "all_integrals"
    REAL(dp), PARAMETER         :: EPS = 1.E-8_dp

    TYPE(neighbor_list_type),     POINTER      :: sbc_3c_neighbor_list
    TYPE(neighbor_node_type),     POINTER      :: sbc_3c_neighbor_node
    TYPE(qlist_type), DIMENSION(:, :, :), &
                                       POINTER :: tmp_neighbor
    TYPE(qlist_type),             POINTER      :: sac_3c_neighbor
    INTEGER  :: bc, i, istat, ic_a, ic_ap,&
                ic_am, ic_b, ic_bm, icoa, icoa_plus, icoam_plus, icoap_plus, icob, &
                icobm, ipgfa, ipgfb,  katom, kkind, kneighbor, knode, na_plus,&
                na, nap, nb, nb_plus, ncoa, ncob, nnode_c 
    INTEGER, DIMENSION(3)                      :: cell_c, la, lap, lam, lb, lbm
    REAL(dp) :: cerf, der_intab_a, der_intab_b, force_a, force_b, &
                rpgfc, zc, zetc, rab2, rac(3), &
                rac2, rbc(3), rbc2, zax2, zbx2
    REAL(dp), DIMENSION(:),   ALLOCATABLE      :: ff
    REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: Int_ab_sum, Int_ab_plus
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    :: vverf, vvnuc

!   ---------------------------------------------------------------------------

    ncoa = npgfa*ncoset(la_max)
    ncob = npgfb*ncoset(lb_max)
    rab2 = dab*dab

    ALLOCATE(Int_ab_sum(ncoa,ncob),STAT = istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routine_name,module_name,__LINE__,"Int_ab_sum",ncoa*ncob*dp_size)
    END IF
    Int_ab_sum = 0.0_dp

    IF (nder > 0) THEN
      na_plus = npgfa*ncoset(la_max+nder)
      nb_plus = npgfb*ncoset(lb_max)
      ALLOCATE(Int_ab_plus(na_plus,nb_plus),STAT=istat)
      IF (istat /= 0) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,"Int_ab_plus",&
                  na_plus*nb_plus*dp_size)
      END IF
      f_a = 0.0_dp
      f_b = 0.0_dp
    END IF

    ALLOCATE(vverf(ncoset(la_max+nder),ncoset(lb_max),la_max+nder+lb_max+1),STAT=istat)
    ALLOCATE(vvnuc(ncoset(la_max+nder),ncoset(lb_max),la_max+nder+lb_max+1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "vverf,vvnuc",2*(ncoset(la_max+nder)*ncoset(lb_max)*(la_max+nder+lb_max+1))*dp_size)
    ALLOCATE(ff(0:la_max+nder+lb_max),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "ff",(1+la_max+nder+lb_max)*dp_size)

    DO kkind = 1,nkind

      IF (.NOT. ASSOCIATED(vall(kkind)%neighbor)) CYCLE

      zetc = vall(kkind)%alpha_c
      zc = REAL(vall(kkind)%zeta_c,dp)
      cerf = vall(kkind)%core_charge
      rpgfc = vall(kkind)%core_radius

      IF (zc == 0.0_dp) CYCLE

      bc = jkind + nkind*(kkind - 1)

       nnode_c = reduced_3c_rho0(kkind)%nnode
       CALL reallocate(vall(kkind)%force_c,1,3,1,nnode_c)

       DO knode=1,nnode_c

         katom = reduced_3c_rho0(kkind)%index_atom(knode)
         rac(1:3) = reduced_3c_rho0(kkind)%rac(1:3,knode)
         rac2 = reduced_3c_rho0(kkind)%rac2(knode)
         rbc(1:3) = reduced_3c_rho0(kkind)%rbc(1:3,knode)
         rbc2 = reduced_3c_rho0(kkind)%rbc2(knode)

        IF(ABS(rac(1)-(rab(1) + rbc(1)))   > EPS .OR. &
             ABS(rac(2)-(rab(2) + rbc(2))) > EPS .OR. &
             ABS(rac(3)-(rab(3) + rbc(3))) > EPS ) THEN
           CYCLE 
        END IF

        IF(nder > 0) THEN

          Int_ab_plus = 0.0_dp
          CALL verfc(la_max+nder,npgfa,zeta,rpgfa,la_min,&
                     lb_max,npgfb,zetb,rpgfb,lb_min,&
                     zetc,rpgfc,zc,cerf,&
                     rab,rab2,rac,rac2,rbc,rbc2,Int_ab_sum,&
                     vverf,vvnuc,ff(0:),&
                     maxder=nder,vabc_plus=Int_ab_plus)

          DO i = 1,3
            force_a = 0.0_dp
            force_b = 0.0_dp

            na = ncoset(la_max)
            nap = ncoset(la_max+nder)
            nb = ncoset(lb_max)

            DO ipgfa = 1,npgfa
              zax2 = zeta(ipgfa)*2.0_dp
              DO ic_a = ncoset(la_min-1)+1,ncoset(la_max) 
                la(1:3) = indco(1:3,ic_a)
                lap(1:3) = la(1:3)
                lap(i) = la(i) + 1
                ic_ap = coset(lap(1),lap(2),lap(3))
                lam(1:3) = la(1:3)
                lam(i) = la(i) - 1
                ic_am = coset(lam(1),lam(2),lam(3))

                icoa =  ic_a + (ipgfa-1)*na
                icoa_plus = ic_a  + (ipgfa-1)*nap
                icoap_plus = ic_ap + (ipgfa-1)*nap
                icoam_plus = ic_am + (ipgfa-1)*nap

                DO ipgfb = 1,npgfb
                  zbx2 = zetb(ipgfb)*2.0_dp
                  DO ic_b = ncoset(lb_min-1)+1,ncoset(lb_max)
                    lb(1:3) = indco(1:3,ic_b)
                    lbm(1:3) = lb(1:3)
                    lbm(i) = lb(i) - 1
                    ic_bm = coset(lbm(1),lbm(2),lbm(3))

                    icob = ic_b + (ipgfb-1)*nb
                    icobm = ic_bm + (ipgfb-1)*nb

                    IF(lam(i)<0) THEN
                      der_intab_a = -zax2*Int_ab_plus(icoap_plus,icob)
                    ELSE
                      der_intab_a = -zax2*Int_ab_plus(icoap_plus,icob)+&
                              REAL(la(i),dp)*Int_ab_plus(icoam_plus,icob)

                    END IF

                    force_a = force_a - pab(icoa,icob)*der_intab_a

                    IF(lbm(i)<0) THEN
                      der_intab_b = -zbx2*(Int_ab_plus(icoap_plus,icob)-&
                                    rab(i)*Int_ab_plus(icoa_plus,icob))
                    ELSE
                      der_intab_b = -zbx2*(Int_ab_plus(icoap_plus,icob)-&
                                    rab(i)*Int_ab_plus(icoa_plus,icob))+&
                              REAL(lb(i),dp)*Int_ab_plus(icoa_plus,icobm)
                    END IF

                    force_b = force_b - pab(icoa,icob)*der_intab_b

                  END DO  ! ic_b
                END DO  ! ipgfb
              END DO  ! ic_a
            END DO  ! ipgfa

            f_a(i) = f_a(i) + force_a
            f_b(i) = f_b(i) + force_b
            vall(kkind)%force_c(i,knode) = vall(kkind)%force_c(i,knode) -&
                                        (force_a+force_b)
          END DO  ! i
        ELSE

          CALL verfc(la_max,npgfa,zeta,rpgfa,la_min,&
                     lb_max,npgfb,zetb,rpgfb,lb_min,&
                     zetc,rpgfc,zc,cerf,&
                     rab,rab2,rac,rac2,rbc,rbc2,Int_ab_sum,&
                     vverf,vvnuc,ff(0:))

        END IF

      END DO  ! knode

    END DO  ! kkind

    hab(1:ncoa,1:ncob) = hab(1:ncoa,1:ncob) + Int_ab_sum(1:ncoa,1:ncob) 

    IF (nder > 0) THEN
      DEALLOCATE(Int_ab_plus,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "Int_ab_plus")
    END IF
    DEALLOCATE(vverf,vvnuc,ff,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "vverf,vvnuc,ff")
    END IF
    DEALLOCATE(Int_ab_sum,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "Int_ab_sum")
    END IF

  END SUBROUTINE all_integrals

!******************************************************************************

 SUBROUTINE radii_3c_orbxall(atomic_kind_set,all_potential)

    TYPE(atomic_kind_type), DIMENSION(:),&
                                  POINTER      :: atomic_kind_set
    TYPE(all_potential_type), POINTER          :: all_potential

    CHARACTER(LEN=*), PARAMETER :: routine_name = "radii_3c_orbxall"
    REAL(dp), PARAMETER         :: EPS = 1.E-6_dp

    TYPE(atomic_kind_type),       POINTER      :: atom_kind
    TYPE(gto_basis_set_type),     POINTER      :: orb_basis

    INTEGER :: i, ic, ico, ico0, icoa, ikind, ipgf, iset, is, istat, &
               j, l0, la, lx, ly, lz, maxlgto, n, ncoc, ncoa, nkind, nseta, zeta_c
    INTEGER, DIMENSION(:),        POINTER      :: la_max, la_min, npgfa 

    REAL(dp) :: cerf, intmax0, r(3,0:60), r2(0:60), rab(3), rab2, rac(3), rac2, ri, &
                rmax, rpgfc, zc, zetc
    REAL(dp), DIMENSION(:),       POINTER      :: orb_rad
    REAL(dp), DIMENSION(:,:),     POINTER      :: set_rad
    REAL(dp), DIMENSION(:,:),     POINTER      :: rpgfa, sphi_a, zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: r_cut_set 
    REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: int_tmp 
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    :: vverf, vvnuc
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: ff

    NULLIFY(orb_rad,set_rad)

    nkind = SIZE(atomic_kind_set,1)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxlgto=maxlgto)

    ALLOCATE(vverf(ncoset(maxlgto),ncoset(maxlgto),1+maxlgto+maxlgto),STAT=istat)
    ALLOCATE(vvnuc(ncoset(maxlgto),ncoset(maxlgto),1+maxlgto+maxlgto),STAT=istat)
    ALLOCATE(ff(0:maxlgto+maxlgto),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "vverf,vvnuc,ff")

    ncoc = 1
    CALL get_potential(potential=all_potential,&
              alpha_core_charge=zetc,&
              zeff=zeta_c,&
              ccore_charge=cerf,&
              core_charge_radius=rpgfc)
 
    zc = REAL(zeta_c,dp)

    orb_rad => all_potential%orb_radius_3c
    set_rad => all_potential%set_radius_3c

    rmax = 20.0_dp
    do i=0,60
      ri=REAL(i,dp)*rmax/60._dp+0.01_dp
      r2(i)=ri**2
      r(:,i)=(/ri,0._dp,0._dp/)
    enddo
    rab(1:3) = 0.0_dp
    rab2 = 0.0_dp

    DO ikind = 1,nkind
      NULLIFY(atom_kind,orb_basis,la_max,la_min,npgfa,rpgfa,zeta)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis)
      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             lmax=la_max,&
                             lmin=la_min,&
                             npgf=npgfa,&
                             nset=nseta,&
                             pgf_radius=rpgfa,&
                             sphi=sphi_a,&
                             zet=zeta)

      orb_rad(ikind) = 0.0_dp
      set_rad(:,ikind) = 0.0_dp
      ALLOCATE(r_cut_set(nseta),STAT=istat)


      DO iset = 1,nseta

        ncoa = npgfa(iset)*ncoset(la_max(iset))
        r_cut_set = 0.0_dp

        ALLOCATE(int_tmp(ncoa,ncoa),STAT=istat)        

        DO i = 10,60 
          rac(1:3) = r(1:3,i)
          rac2 = r2(i)

          int_tmp = 0.0_dp
          intmax0 = 0.0_dp
          CALL verfc(la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                     la_max(iset),npgfa(iset),zeta(:,iset),rpgfa(:,iset),la_min(iset),&
                     zetc,rpgfc,zc,cerf,&
                     rab,rab2,rac,rac2,rac,rac2,int_tmp,vverf,vvnuc,ff(0:))

          DO ipgf = 1,npgfa(iset)
             n=(ipgf-1)*ncoset(la_max(iset))
             DO ico = 1,ncoset(la_max(iset))
               icoa = ico + n          

               intmax0 = MAX(intmax0,ABS(int_tmp(icoa,icoa)))
             END DO
          END DO

          IF( intmax0 < eps ) THEN
            r_cut_set(iset) = r(1,i)
            EXIT
          END IF
  
        END DO  ! i

        orb_rad(ikind) = MAX(orb_rad(ikind) ,r_cut_set(iset) ) 
        set_rad(iset,ikind) = r_cut_set(iset)*r_cut_set(iset)  

        DEALLOCATE(int_tmp,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "int_tmp")

      END DO  ! iset
      
      DEALLOCATE(r_cut_set,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "r_cut_set")

!  write(*,'(I4,5f12.5)') ikind , orb_rad(ikind),set_rad(1:nseta,ikind)
    END DO  ! ikind

    DEALLOCATE(vverf,vvnuc,ff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "vverf,vvnuc,ff")

  END SUBROUTINE radii_3c_orbxall

!******************************************************************************

END MODULE qs_all_potential

!******************************************************************************
