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
                                             qlist_type
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
    INTEGER, DIMENSION(:), POINTER              :: index_c
  END TYPE vall_type

! Public Types

  PUBLIC :: vall_type

! Public Subroutine
 
  PUBLIC :: all_integrals

!******************************************************************************
 
  CONTAINS
 
!******************************************************************************
  SUBROUTINE  all_integrals(la_max,la_min,npgfa,rpgfa,zeta,&
                            lb_max,lb_min,npgfb,rpgfb,zetb,&
                            rab,dab,jkind,jatom,nkind,sbc_3c,&
                            vall,hab,nder,pab,f_a,f_b)

    INTEGER, INTENT(IN)                        :: la_max, la_min, npgfa
    REAL(dp), DIMENSION(:), INTENT(IN)         :: rpgfa, zeta
    INTEGER, INTENT(IN)                        :: lb_max, lb_min, npgfb
    REAL(dp), DIMENSION(:), INTENT(IN)         :: rpgfb, zetb
    REAL(dp), DIMENSION(3), INTENT(IN)         :: rab
    REAL(dp), INTENT(IN)                       :: dab
    INTEGER, INTENT(IN)                        :: jkind,jatom,nkind
    TYPE(neighbor_list_set_p_type), DIMENSION(:),&
              POINTER                          :: sbc_3c
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
    REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: Int_ab_sum, Int_ab_plus

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
    ELSE
      ALLOCATE(Int_ab_plus(1,1),STAT=istat)
      IF (istat /= 0) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,"Int_ab_plus",1*dp_size)
      END IF
    END IF

    DO kkind = 1,nkind

      IF (.NOT. ASSOCIATED(vall(kkind)%neighbor)) CYCLE

      zetc = vall(kkind)%alpha_c
      zc = REAL(vall(kkind)%zeta_c,dp)
      cerf = vall(kkind)%core_charge
      rpgfc = vall(kkind)%core_radius

      IF (zc == 0.0_dp) CYCLE

      bc = jkind + nkind*(kkind - 1)

      IF (.NOT. ASSOCIATED(sbc_3c(bc)%neighbor_list_set)) CYCLE
      NULLIFY(sbc_3c_neighbor_list)

      sbc_3c_neighbor_list => find_neighbor_list(neighbor_list_set=&
           sbc_3c(bc)%neighbor_list_set,atom=jatom)

      CALL get_neighbor_list(neighbor_list=sbc_3c_neighbor_list, &
           nnode = nnode_c)

!      IF (nder > 0) THEN
        CALL reallocate(vall(kkind)%force_c,1,3,1,nnode_c)
        CALL reallocate(vall(kkind)%index_c,1,nnode_c)
!      END IF

      sbc_3c_neighbor_node => first_node(sbc_3c_neighbor_list)

      knodeloop: DO knode=1,nnode_c

        CALL get_neighbor_node(neighbor_node=sbc_3c_neighbor_node,&
             neighbor=katom,&
                cell=cell_c,&
                       r=rbc)

        rbc2= rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3) 

        vall(kkind)%index_c(knode) = katom

        tmp_neighbor => vall(kkind)%neighbor
        DO i=1, 3
          IF ( (cell_c(i) < LBOUND(tmp_neighbor,i)) .OR. &
                (cell_c(i) > UBOUND(tmp_neighbor,i))) THEN
            sbc_3c_neighbor_node => next(sbc_3c_neighbor_node)
            CYCLE knodeloop
          END IF
        END DO
        sac_3c_neighbor => vall(kkind)%neighbor(cell_c(1),&
                                                cell_c(2),&
                                                cell_c(3))
        IF(sac_3c_neighbor%n == 0) THEN
          sbc_3c_neighbor_node => next(sbc_3c_neighbor_node)
          CYCLE knodeloop
        END IF

        kneighbor = locate(sac_3c_neighbor%list,katom)

        IF (kneighbor == 0) THEN
           sbc_3c_neighbor_node => next(sbc_3c_neighbor_node)
           CYCLE knodeloop
        END IF

        rac(:) = sac_3c_neighbor%r(:,kneighbor)
        rac2= rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3) 
                
        IF(ABS(rac(1)-(rab(1) + rbc(1)))   > EPS .OR. &
             ABS(rac(2)-(rab(2) + rbc(2))) > EPS .OR. &
             ABS(rac(3)-(rab(3) + rbc(3))) > EPS ) THEN
           sbc_3c_neighbor_node => next(sbc_3c_neighbor_node)
           CYCLE knodeloop
        END IF

        IF(nder > 0) THEN

          Int_ab_plus = 0.0_dp
          CALL verfc(la_max+nder,npgfa,zeta,rpgfa,la_min,&
                     lb_max,npgfb,zetb,rpgfb,lb_min,&
                     zetc,rpgfc,zc,cerf,&
                     rab,rab2,rac,rac2,rbc,rbc2,Int_ab_sum,&
                     nder,Int_ab_plus)

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

          Int_ab_plus = 0.0_dp

          CALL verfc(la_max,npgfa,zeta,rpgfa,la_min,&
                     lb_max,npgfb,zetb,rpgfb,lb_min,&
                     zetc,rpgfc,zc,cerf,&
                     rab,rab2,rac,rac2,rbc,rbc2,Int_ab_sum,&
                     nder,Int_ab_plus)

        END IF

        sbc_3c_neighbor_node => next(sbc_3c_neighbor_node)
      END DO knodeloop  ! knode

    END DO  ! kkind

    hab(1:ncoa,1:ncob) = hab(1:ncoa,1:ncob) + Int_ab_sum(1:ncoa,1:ncob) 

    DEALLOCATE(Int_ab_sum,Int_ab_plus,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "Int_ab_sum,Int_ab_plus")
    END IF

  END SUBROUTINE all_integrals

END MODULE qs_all_potential
!******************************************************************************
