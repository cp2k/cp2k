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
  USE orbital_pointers,                ONLY: ncoset
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

    TYPE(neighbor_list_type),     POINTER      :: sbc_3c_neighbor_list
    TYPE(neighbor_node_type),     POINTER      :: sbc_3c_neighbor_node
    TYPE(qlist_type), DIMENSION(:, :, :), &
                                       POINTER :: tmp_neighbor
    TYPE(qlist_type),             POINTER      :: sac_3c_neighbor
    INTEGER  :: bc, i, istat, katom, kkind, kneighbor, knode, na_plus, nb_plus, ncoa, ncob, nnode_c 
    INTEGER, DIMENSION(3)                      :: cell_c
    REAL(dp) :: cerf, rpgfc, zc, zetc, rab2, rac(3), rac2, rbc(3), rbc2
    REAL(dp), DIMENSION(:,:), ALLOCATABLE      :: Int_ab_sum, Int_ab_plus
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    :: Int_dab_sum, Int_adb_sum

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
      ALLOCATE(Int_dab_sum(ncoa,ncob,3),STAT=istat)
      ALLOCATE(Int_adb_sum(ncoa,ncob,3),STAT=istat)
      IF (istat /= 0) THEN
         CALL stop_memory(routine_name,module_name,__LINE__,"Int_dab_plus,Int_adb_plus",&
                  2*3*ncoa*ncob*dp_size)
      END IF
      na_plus = npgfa*ncoset(la_max+nder)
      nb_plus = npgfb*ncoset(lb_max+nder)
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
      NULLIFY(vall(kkind)%force_c)
      NULLIFY(vall(kkind)%index_c)
      IF (zc == 0.0_dp) CYCLE

      bc = jkind + nkind*(kkind - 1)

      IF (.NOT. ASSOCIATED(sbc_3c(bc)%neighbor_list_set)) CYCLE
      NULLIFY(sbc_3c_neighbor_list)

      sbc_3c_neighbor_list => find_neighbor_list(neighbor_list_set=&
           sbc_3c(bc)%neighbor_list_set,atom=jatom)

      CALL get_neighbor_list(neighbor_list=sbc_3c_neighbor_list, &
           nnode = nnode_c)

      CALL reallocate(vall(kkind)%force_c,1,3,1,nnode_c)
      CALL reallocate(vall(kkind)%index_c,1,nnode_c)

      sbc_3c_neighbor_node => first_node(sbc_3c_neighbor_list)

      knodeloop: DO knode=1,nnode_c

        CALL get_neighbor_node(neighbor_node=sbc_3c_neighbor_node,&
             neighbor=katom,&
                cell=cell_c,&
                       r=rbc)

        rbc2= rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3) 
        rac(:) = rab(:) + rbc(:)
        rac2= rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3) 

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

        IF(nder > 0) THEN

          Int_ab_plus = 0.0_dp
          CALL verfc(la_max+nder,npgfa,zeta,rpgfa,la_min,&
                     lb_max+nder,npgfb,zetb,rpgfb,lb_min,&
                     zetc,rpgfc,zc,cerf,&
                     rab,rab2,rac,rac2,rbc,rbc2,Int_ab_sum,&
                     nder,Int_ab_plus)

        ELSE

          Int_ab_plus = 0.0_dp

          CALL verfc(la_max,npgfa,zeta,rpgfa,la_min,&
                     lb_max,npgfb,zetb,rpgfb,lb_min,&
                     zetc,rpgfc,zc,cerf,&
                     rab,rab2,rac,rac2,rbc,rbc2,Int_ab_sum,&
                     nder,Int_ab_plus)

        END IF

      END DO knodeloop  ! knode

    END DO  ! kkind

    hab(1:ncoa,1:ncob) = hab(1:ncoa,1:ncob) + Int_ab_sum(1:ncoa,1:ncob) 

    DEALLOCATE(Int_ab_sum,Int_ab_plus,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "Int_ab_sum,Int_ab_plus")
    END IF
    IF(nder > 0) THEN
      DEALLOCATE(Int_dab_sum,Int_adb_sum,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,&
                       "Int_dab_sum,Int_adb_sum")
      END IF
    END IF

  END SUBROUTINE all_integrals

END MODULE qs_all_potential
!******************************************************************************
