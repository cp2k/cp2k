!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002-1004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_ks_atom [1.0] *
!!
!!   NAME
!!     qs_ks_atom
!!
!!   FUNCTION
!!     routines that build the Kohn-Sham matrix  contributions
!!     coming from local atomic densities
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!
!!
!!
!!   SOURCE
!****************************************************************************
MODULE qs_ks_atom
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level
  USE hartree_local_types,             ONLY: ecoul_1center_type,&
                                             get_ecoul_1c
  USE kinds,                           ONLY: dp, int_size
  USE memory_utilities,                ONLY: reallocate
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next,&
                                             qlist_type
  USE qs_oce_types,                    ONLY: oce_matrix_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: trace_r_AxB
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='qs_ks_atom'

  PUBLIC :: update_ks_atom

!******************************************************************************

CONTAINS

!******************************************************************************

  SUBROUTINE update_ks_atom(qs_env,ks,p,ispin,forces,error)

    TYPE(qs_environment_type), INTENT(IN)    :: qs_env
    TYPE(real_matrix_p_type), INTENT(INOUT)  :: ks
    TYPE(real_matrix_p_type), INTENT(IN)     :: p
    INTEGER                                  :: ispin
    LOGICAL, INTENT(IN)                      :: forces
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER              :: routine = "update_ks_atom"

    INTEGER :: i, ia_kind, iab, iatom, ibc, ikind, ilist, inode, &
               istat, ja_kind, jatom, jkind, ka_kind, katom, &
               kkind, kneighbor, knode, last_jatom, max_gau, max_nsgf, &
               maxsoc, natom, nkind, nlist, nnode, nnode_c, &
               nsetc, nsgfa, nsgfb, nsoctot
    INTEGER, DIMENSION(3)                    :: cell_c
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind
    LOGICAL                                  :: paw_atom, new_atom_b
    REAL(dp), DIMENSION(3)                   :: rbc
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C_coeff_h, C_coeff_s
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block, p_block
    REAL(dp), DIMENSION(3,3)                 :: force_tmp
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_oce, sbc_oce
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local, &
                                                sbc_oce_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node, &
                                                sbc_oce_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(ecoul_1center_type), DIMENSION(:), & 
                                 POINTER     :: ecoul_1c
    TYPE(qlist_type), DIMENSION(:, :, :), &
      POINTER                                :: tmp_neigh
    TYPE(qlist_type), POINTER                :: sac_oce_neighbor
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom
    TYPE(rho_atom_type), POINTER             :: rho_at

!******************************************************************************

    NULLIFY(atomic_kind_set,force,oce,rho_atom,sab_orb,sac_oce,sbc_oce)
    NULLIFY(orb_basis,ecoul_1c)

    CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           force=force,&
           oce=oce, &
           rho_atom_set=rho_atom,&
           sab_orb=sab_orb,&
           sac_oce=sac_oce,&
           sbc_oce=sbc_oce,&
           ecoul_1c=ecoul_1c)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxsgf=max_nsgf,&
                             maxgtops=max_gau,natom=natom)
    IF(forces) THEN
      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"atom_of_kind",natom*int_size)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)
    END IF


    ALLOCATE (C_coeff_h(max_gau,max_nsgf))
    ALLOCATE (C_coeff_s(max_gau,max_nsgf))

    nkind = SIZE(atomic_kind_set)

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis,&
                           paw_atom=paw_atom)
      CALL  get_gto_basis_set(gto_basis_set=orb_basis,&
                                nsgf=nsgfa)

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)

        CALL get_atomic_kind(atomic_kind=atomic_kind, &
                             orb_basis_set=orb_basis)
        CALL  get_gto_basis_set(gto_basis_set=orb_basis,&
                                nsgf=nsgfb)

        iab = ikind + nkind*(jkind - 1)

        IF (.NOT.ASSOCIATED(sab_orb(iab)%neighbor_list_set)) CYCLE

        CALL get_neighbor_list_set(neighbor_list_set=&
                                   sab_orb(iab)%neighbor_list_set,&
                                   nlist=nlist)

        NULLIFY ( sab_orb_neighbor_list )

        DO ilist=1,nlist

          IF ( .NOT. ASSOCIATED(sab_orb_neighbor_list) ) THEN
            sab_orb_neighbor_list => first_list(sab_orb(iab)%neighbor_list_set)
          ELSE
            sab_orb_neighbor_list => next(sab_orb_neighbor_list)
          END IF
          sab_orb_neighbor_list_local => sab_orb_neighbor_list

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom,&
                                 nnode=nnode)

          last_jatom = 0

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode
            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
            END IF

            IF (.NOT. new_atom_b) THEN
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

            NULLIFY(h_block, p_block)
            IF (iatom <= jatom) THEN
              CALL get_block_node(matrix=ks%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=h_block)
            ELSE
              CALL get_block_node(matrix=ks%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  BLOCK=h_block)
            END IF

            IF(forces) THEN
              IF (iatom <= jatom) THEN
                CALL get_block_node(matrix=p%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=p_block)
              ELSE
                CALL get_block_node(matrix=p%matrix,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    BLOCK=p_block)
              END IF
            END IF

            DO kkind=1,nkind

              IF(.NOT.ASSOCIATED(oce%voce(kkind,iatom)%neighbor)) CYCLE

              ibc = jkind + nkind*(kkind - 1)

              IF (.NOT.ASSOCIATED(sbc_oce(ibc)%neighbor_list_set)) CYCLE
              NULLIFY(sbc_oce_neighbor_list)

              sbc_oce_neighbor_list =>&
                         find_neighbor_list(neighbor_list_set=&
                         sbc_oce(ibc)%neighbor_list_set,&
                         atom=jatom)
              IF(.NOT. ASSOCIATED(sbc_oce_neighbor_list)) CYCLE

              atomic_kind => atomic_kind_set(kkind)

              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   orb_basis_set=orb_basis,&
                                   paw_atom=paw_atom)

              CALL  get_gto_basis_set(gto_basis_set=orb_basis,&
                                      nset =nsetc,maxso=maxsoc)
              nsoctot = maxsoc*nsetc

              IF(.NOT. paw_atom)  CYCLE

              CALL get_neighbor_list(neighbor_list=sbc_oce_neighbor_list,&
                                     nnode=nnode_c)
              sbc_oce_neighbor_node => first_node(sbc_oce_neighbor_list)

              knodeloop: DO knode=1,nnode_c

                CALL get_neighbor_node(neighbor_node=sbc_oce_neighbor_node,&
                          neighbor=katom,&
                          cell=cell_c,&
                          r=rbc)
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

!               *** Locate operator atom in the sac_ppnl neighbor list ***
                kneighbor = locate(sac_oce_neighbor%list,katom)
                IF (kneighbor == 0) THEN
                    sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
                    CYCLE
                END IF
                rho_at => rho_atom(katom)

                CALL add_vhxca_to_ks(h_block,oce,rho_at,ecoul_1c,ispin, &
                                    iatom,nsgfa,jatom,nsgfb,katom,nsoctot,&
                                    C_coeff_h,C_coeff_s,max_gau)

                IF(forces)THEN
                  ia_kind = atom_of_kind(iatom)
                  ja_kind = atom_of_kind(jatom)
                  ka_kind = atom_of_kind(katom)
                  IF(iatom <= jatom) THEN
                    rho_at => rho_atom(katom)
                    force_tmp(1:3,1:3) = 0.0_dp 
                    CALL add_vhxca_forces(p_block,oce,rho_at,force_tmp,ecoul_1c,&
                                    ispin,iatom,nsgfa,jatom,nsgfb,katom,nsoctot,&
                                    C_coeff_h,C_coeff_s,max_gau)
                    force(kkind)%vhxc_atom(1:3,ka_kind) = &
                                          force(kkind)%vhxc_atom(1:3,ka_kind) + &
                                          force_tmp(1:3,3) 
                    force(ikind)%vhxc_atom(1:3,ia_kind) = &
                                         force(ikind)%vhxc_atom(1:3,ia_kind) + &
                                         force_tmp(1:3,1) 
                    force(jkind)%vhxc_atom(1:3,ja_kind) = &
                                         force(jkind)%vhxc_atom(1:3,ja_kind) + &
                                         force_tmp(1:3,2) 
                  ELSE
                    rho_at => rho_atom(katom)
                    force_tmp(1:3,1:3) = 0.0_dp 
                    CALL add_vhxca_forces(p_block,oce,rho_at,force_tmp,ecoul_1c,&
                                    ispin,jatom,nsgfb,iatom,nsgfa,katom,nsoctot,&
                                    C_coeff_h,C_coeff_s,max_gau)
                    force(kkind)%vhxc_atom(1:3,ka_kind) = &
                                          force(kkind)%vhxc_atom(1:3,ka_kind) + &
                                          force_tmp(1:3,3) 
                    force(ikind)%vhxc_atom(1:3,ia_kind) = &
                                          force(ikind)%vhxc_atom(1:3,ia_kind) + &
                                          force_tmp(1:3,2) 
                    force(jkind)%vhxc_atom(1:3,ja_kind) = &
                                          force(jkind)%vhxc_atom(1:3,ja_kind) + &
                                          force_tmp(1:3,1) 
                  END IF
                END IF

                sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
              ENDDO knodeloop    ! knode
            ENDDO    ! kkind

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          ENDDO    ! inode
        ENDDO    ! ilist 
      ENDDO    ! jkind
    ENDDO    ! ikind

    DEALLOCATE (C_coeff_h)
    DEALLOCATE (C_coeff_s)
    IF(forces) DEALLOCATE (atom_of_kind)


  END SUBROUTINE update_ks_atom
!******************************************************************************

  SUBROUTINE add_vhxca_to_ks(h_block,oce,rho_atom,ecoul_1c,ispin,ia,na,ja,nb,&
                             ka,nsp,C_h,C_s,ldc)


    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(ecoul_1center_type), DIMENSION(:), & 
                                 POINTER     :: ecoul_1c
    INTEGER, INTENT(IN)                      :: ispin, ia, na, ja, nb, ka, &
                                                nsp, ldc
    REAL(KIND=dp), DIMENSION(ldc, nb)        :: C_s, C_h

    CHARACTER(LEN=*), PARAMETER              :: routine = "add_vhxca_to_ks"

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: c_block, IntC_h, IntC_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_xc_h, int_xc_s
    REAL(dp), DIMENSION(:,:), POINTER        :: int_h_h, int_h_s, &
                                                int_h_0h, int_h_0s,&
                                                int_h_2c, int_h_3c
    REAL(dp), DIMENSION(:,:), POINTER        :: int_hard, int_soft
    INTEGER                                  :: is1, is2

!******************************************************************************

    NULLIFY(int_hard,int_soft,intC_h,intC_s)
    NULLIFY(int_xc_h, int_xc_s,int_h_h, int_h_s, int_h_0h, int_h_0s)
    NULLIFY(int_h_2c, int_h_3c)

    CALL reallocate(int_hard,1,nsp,1,nsp)
    CALL reallocate(int_soft,1,nsp,1,nsp)

    CALL get_rho_atom(rho_atom=rho_atom,&
                     ga_Vxc_gb_h=int_xc_h,ga_Vxc_gb_s=int_xc_s)

    CALL get_ecoul_1c(ecoul_1c=ecoul_1c,iatom=ka,aVh1b_h=int_h_h,&
                      aVh1b_s=int_h_s,aVh1b_0=int_h_0h,&
                      aVhsb_0s=int_h_0s, &
                      aVh2c0b=int_h_2c, aVh3c0b=int_h_3c)

    DO is1 = 1,nsp
      DO is2 = is1,nsp
        int_hard(is2,is1) = int_xc_h(ispin)%r_coef(is2,is1) + &
                         int_h_h(is2,is1) - int_h_0h(is2,is1) +&
                         int_h_0s(is2,is1) + &
                         int_h_2c(is2,is1) + int_h_3c(is2,is1)
        int_soft(is2,is1) = int_xc_s(ispin)%r_coef(is2,is1) + &
                         int_h_s(is2,is1) - int_h_0h(is2,is1) +&
                         int_h_0s(is2,is1) + &
                         int_h_2c(is2,is1) + int_h_3c(is2,is1)
        int_hard(is1,is2) = int_hard(is2,is1)
        int_soft(is1,is2) = int_soft(is2,is1)
      END DO
    END DO

    IF (ia <= ja) THEN
      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=c_block)
      C_h = 0.0_dp
      C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

! write(*,*) 'C_h_b ',nsp,nb
!  do i = 1,nsp
!    write(*,'(i3,20f10.5)') i,C_h(i,1:nb)
!  enddo


      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=c_block)
      C_s = 0.0_dp
      C_s(1:nsp,1:nb)=c_block(1:nsp,1:nb)

      CALL reallocate(intC_h,1,nsp,1,nb)
      CALL reallocate(intC_s,1,nsp,1,nb)

      CALL dgemm("N","N",nsp,nb,nsp,&
               1.0_dp,int_hard,nsp,&
               C_h,ldc,0.0_dp,&
               intC_h,nsp)

      CALL dgemm("N","N",nsp,nb,nsp,&
               1.0_dp,int_soft,nsp,&
               C_s,ldc,0.0_dp,&
               intC_s,nsp)

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=c_block)
      C_h = 0.0_dp
      C_h(1:nsp,1:na)=c_block(1:nsp,1:na)

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=c_block)
      C_s = 0.0_dp
      C_s(1:nsp,1:na)=c_block(1:nsp,1:na)

      CALL dgemm("T","N",na,nb,nsp,&
               1.0_dp,C_h,ldc,&
               intC_h,nsp,1.0_dp,&
               h_block,SIZE(h_block,1))
       
      CALL dgemm("T","N",na,nb,nsp,&
               -1.0_dp,C_s,ldc,&
               intC_s,nsp,1.0_dp,&
               h_block,SIZE(h_block,1))

    ELSE

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=c_block)
      C_h = 0.0_dp
      C_h(1:nsp,1:na)=c_block(1:nsp,1:na)

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=c_block)
      C_s = 0.0_dp
      C_s(1:nsp,1:na)=c_block(1:nsp,1:na)

      CALL reallocate(intC_h,1,nsp,1,na)
      CALL reallocate(intC_s,1,nsp,1,na)

      CALL dgemm("N","N",nsp,na,nsp,&
               1.0_dp,int_hard,nsp,&
               C_h,ldc,0.0_dp,&
               intC_h,nsp)

      CALL dgemm("N","N",nsp,na,nsp,&
               1.0_dp,int_soft,nsp,&
               C_s,ldc,0.0_dp,&
               intC_s,nsp)

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=c_block)
      C_h = 0.0_dp
      C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=c_block)
      C_s = 0.0_dp
      C_s(1:nsp,1:nb)=c_block(1:nsp,1:nb)

      CALL dgemm("T","N",nb,na,nsp,&
               1.0_dp,C_h,ldc,&
               intC_h,nsp,1.0_dp,&
               h_block,SIZE(h_block,1))

      CALL dgemm("T","N",nb,na,nsp,&
               -1.0_dp,C_s,ldc,&
               intC_s,nsp,1.0_dp,&
               h_block,SIZE(h_block,1))

    END IF

  END SUBROUTINE add_vhxca_to_ks

! *****************************************************************************

  SUBROUTINE add_vhxca_forces(p_block,oce,rho_atom,force,ecoul_1c,ispin,&
                            ia,na,ja,nb,ka,nsp,C_h,C_s,ldc)

! *****************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(3,3)                 :: force
    TYPE(ecoul_1center_type), DIMENSION(:), & 
                                 POINTER     :: ecoul_1c
    INTEGER, INTENT(IN)                      :: ispin, ia, na, ja, nb, ka, &
                                                nsp, ldc
    REAL(KIND=dp), DIMENSION(ldc, nb)        :: C_s, C_h

    CHARACTER(LEN=*), PARAMETER              :: routine = "add_vhxca_forces"

    INTEGER                                  :: dir
    REAL(KIND=dp)                            :: ieqj, trace
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: c_block, dCPC_h, dCPC_s, &
                                                p_blockt, PC_h, PC_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_xc_h, int_xc_s
    REAL(dp), DIMENSION(:,:), POINTER        :: int_h_h, int_h_s, &
                                                int_h_0h, int_h_0s,&
                                                int_h_2c, int_h_3c
    REAL(dp), DIMENSION(:,:), POINTER        :: int_hard, int_soft
    INTEGER                                  :: is1, is2

! *****************************************************************************

!   *** factor 2 because only half of the pairs with ia =/ ja are  considered
    ieqj = 2.0_dp
!   iatom == jatom
    IF(ia == ja) ieqj = 1.0_dp

    NULLIFY(int_hard, int_soft, PC_h, PC_s, dCPC_h, dCPC_s)
    NULLIFY(int_xc_h, int_xc_s,int_h_h, int_h_s, int_h_0h, int_h_0s)
    NULLIFY(int_h_2c, int_h_3c)

    CALL reallocate(int_hard,1,nsp,1,nsp)
    CALL reallocate(int_soft,1,nsp,1,nsp)

    CALL get_rho_atom(rho_atom=rho_atom,&
                     ga_Vxc_gb_h=int_xc_h,ga_Vxc_gb_s=int_xc_s)

    CALL get_ecoul_1c(ecoul_1c=ecoul_1c,iatom=ka,aVh1b_h=int_h_h,&
                      aVh1b_s=int_h_s,aVh1b_0=int_h_0h,aVhsb_0s=int_h_0s,&
                      aVh2c0b=int_h_2c,aVh3c0b=int_h_3c)

    DO is1 = 1,nsp
      DO is2 = is1,nsp
        int_hard(is2,is1) = int_xc_h(ispin)%r_coef(is2,is1) + &
                         int_h_h(is2,is1) - int_h_0h(is2,is1) + &
                         int_h_0s(is2,is1) + &
                         int_h_2c(is2,is1) + int_h_3c(is2,is1)
        int_soft(is2,is1) = int_xc_s(ispin)%r_coef(is2,is1) + &
                         int_h_s(is2,is1) - int_h_0h(is2,is1) + &
                         int_h_0s(is2,is1) + &
                         int_h_2c(is2,is1) + int_h_3c(is2,is1)
        int_hard(is1,is2) = int_hard(is2,is1)
        int_soft(is1,is2) = int_soft(is2,is1)
      END DO
    END DO

! i-k contributions

    NULLIFY(c_block)
    CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=c_block)
    C_h = 0.0_dp
    C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

    NULLIFY(c_block)
    CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=c_block)
    C_s = 0.0_dp
    C_s(1:nsp,1:nb)=c_block(1:nsp,1:nb)

    CALL reallocate(PC_h,1,na,1,nsp)
    CALL reallocate(PC_s,1,na,1,nsp)

    CALL dgemm("N","T",na,nsp,nb,&
               1.0_dp,p_block,SIZE(p_block,1),&
               C_h,ldc,0.0_dp,&
               PC_h,na)

    CALL dgemm("N","N",na,nsp,nb,&
               1.0_dp,p_block,SIZE(p_block,1),&
               C_s,ldc,0.0_dp,&
               PC_s,na)

    DO dir = 2,4
       NULLIFY(c_block)
       CALL get_block_node(matrix=oce%hardo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=c_block)
       C_h = 0.0_dp
       C_h(1:nsp,1:na)=c_block(1:nsp,1:na)

       NULLIFY(c_block)
       CALL get_block_node(matrix=oce%softo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=c_block)
       C_s = 0.0_dp
       C_s(1:nsp,1:na)=c_block(1:nsp,1:na)

       CALL reallocate(dCPC_h,1,nsp,1,nsp)
       CALL reallocate(dCPC_s,1,nsp,1,nsp)


       CALL dgemm("N","N",nsp,nsp,na,&
                  1.0_dp,C_h,ldc,&
                  PC_h,na,0.0_dp,&
                  dCPC_h,nsp)

       trace = trace_r_AxB(dCPC_h,nsp,int_hard,nsp,nsp,nsp)

!       force(dir-1,3) = force(dir-1,3) - 2.0_dp * trace
!       force(dir-1,1) = force(dir-1,1) + 2.0_dp * trace
       force(dir-1,3) = force(dir-1,3) - ieqj * trace
       force(dir-1,1) = force(dir-1,1) + ieqj * trace

       CALL dgemm("N","N",nsp,nsp,na,&
                  1.0_dp,C_s,ldc,&
                  PC_s,na,0.0_dp,&
                  dCPC_s,nsp)

       trace = trace_r_AxB(dCPC_s,nsp,int_soft,nsp,nsp,nsp)

!       force(dir-1,3) = force(dir-1,3) + 2.0_dp * trace
!       force(dir-1,1) = force(dir-1,1) - 2.0_dp * trace
       force(dir-1,3) = force(dir-1,3) +  ieqj * trace
       force(dir-1,1) = force(dir-1,1) -  ieqj * trace

    END DO

! j-k contributioms

    NULLIFY( PC_h, PC_s, dCPC_h, dCPC_s)

    p_blockt = TRANSPOSE(p_block)

    NULLIFY(c_block)
    CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=c_block)
    C_h = 0.0_dp
    C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

    NULLIFY(c_block)
    CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        BLOCK=c_block)
    C_s = 0.0_dp
    C_s(1:nsp,1:na)=c_block(1:nsp,1:na)

    CALL reallocate(PC_h,1,nb,1,nsp)
    CALL reallocate(PC_s,1,nb,1,nsp)

    CALL dgemm("N","T",nb,nsp,na,&
               1.0_dp,p_blockt,SIZE(p_blockt,1),&
               C_h,ldc,0.0_dp,&
               PC_h,nb)

    CALL dgemm("N","N",nb,nsp,na,&
               1.0_dp,p_blockt,SIZE(p_blockt,1),&
               C_s,ldc,0.0_dp,&
               PC_s,nb)

    DO dir = 2,4
       NULLIFY(c_block)
       CALL get_block_node(matrix=oce%hardo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=c_block)
       C_h = 0.0_dp
       C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

       NULLIFY(c_block)
       CALL get_block_node(matrix=oce%softo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        BLOCK=c_block)
       C_s = 0.0_dp
       C_s(1:nsp,1:nb)=c_block(1:nsp,1:nb)

       CALL reallocate(dCPC_h,1,nsp,1,nsp)
       CALL reallocate(dCPC_s,1,nsp,1,nsp)

       CALL dgemm("N","N",nsp,nsp,nb,&
                  1.0_dp,C_h,ldc,&
                  PC_h,na,0.0_dp,&
                  dCPC_h,nsp)

       trace = trace_r_AxB(dCPC_h,nsp,int_hard,nsp,nsp,nsp)

!       force(dir-1,3) = force(dir-1,3) - 2.0_dp * trace
!       force(dir-1,2) = force(dir-1,2) + 2.0_dp * trace
       force(dir-1,3) = force(dir-1,3) -  ieqj * trace
       force(dir-1,2) = force(dir-1,2) +  ieqj * trace

       CALL dgemm("N","N",nsp,nsp,nb,&
                  1.0_dp,C_s,ldc,&
                  PC_s,na,0.0_dp,&
                  dCPC_s,nsp)

       trace = trace_r_AxB(dCPC_s,nsp,int_soft,nsp,nsp,nsp)

!       force(dir-1,3) = force(dir-1,3) + 2.0_dp * trace
!       force(dir-1,2) = force(dir-1,2) - 2.0_dp * trace
       force(dir-1,3) = force(dir-1,3) +  ieqj * trace
       force(dir-1,2) = force(dir-1,2) -  ieqj * trace

    END DO

  END SUBROUTINE add_vhxca_forces

! *****************************************************************************

END MODULE qs_ks_atom
