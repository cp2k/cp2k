!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002-1004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_vxc_atom [1.0] *
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
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             list_type, &
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next
  USE qs_oce_types,                    ONLY: oce_matrix_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: trace_r_AxB
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type,&
                                             transfer_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='qs_ks_atom'

  PUBLIC :: update_ks_vxc_atom

!******************************************************************************

CONTAINS

!******************************************************************************

  SUBROUTINE update_ks_vxc_atom(qs_env,ks,p,ispin,forces,error)

    TYPE(real_matrix_p_type), INTENT(INOUT)  :: ks
    TYPE(real_matrix_p_type), INTENT(IN)     :: p
    TYPE(qs_environment_type), INTENT(IN)    :: qs_env
    LOGICAL, INTENT(IN)                      :: forces

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine = "update_ks_vxc_atom"

    TYPE(atomic_kind_type), DIMENSION(:),&
                            POINTER          :: atomic_kind_set
    TYPE(atomic_kind_type),  POINTER         :: atomic_kind
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(gto_basis_set_type), POINTER        :: orb_basis

    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(rho_atom_type), DIMENSION(:), &
                                   POINTER   :: rho_atom
    TYPE(rho_atom_type), POINTER             :: rho_at
    TYPE(neighbor_list_set_p_type), &
               DIMENSION(:),  POINTER        :: sab_orb, sac_oce, sbc_oce

    TYPE(list_type), POINTER                 :: sac_oce_neighbor
    TYPE(neighbor_list_type), POINTER        :: sac_oce_neighbor_list, &
                                                sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local,&
                                                sbc_oce_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sac_oce_neighbor_node, &
                                                sab_orb_neighbor_node, &
                                                sbc_oce_neighbor_node
    REAL(dp), DIMENSION(3)                   :: rbc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: c_block,h_block, &
                                                p_block
    REAL(KIND = dp), DIMENSION(:,:), &
                                ALLOCATABLE  :: C_coeff_h,&
                                                C_coeff_s
 
    LOGICAL                                  :: paw_atom

    INTEGER     ::  i,iab, iac, iatom, ibc, ikind, ilist, inode, ispin, istat, jatom, &
                    jkind, katom, kkind, kneighbor, knode, max_gau, max_nsgf, &
                    maxsoc, natom, nsoctot, nkind, nlist, nnode, &
                    nneighbor, nnode_c, nsgfa, nsgfb, nsgfc, nsetc
    INTEGER, DIMENSION(3)                    :: cell_c

!******************************************************************************

    NULLIFY(atomic_kind_set,force,oce,rho_atom,sab_orb,sac_oce,sbc_oce)
    NULLIFY(orb_basis)

    CALL get_qs_env(qs_env=qs_env,&
           atomic_kind_set=atomic_kind_set,&
           force=force,&
           oce=oce, &
           rho_atom_set=rho_atom,&
           sab_orb=sab_orb,&
           sac_oce=sac_oce,&
           sbc_oce=sbc_oce)


    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,maxsgf=max_nsgf,&
                             maxgtops=max_gau,natom=natom)

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
          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode
            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom)

            IF (iatom <= jatom) THEN
              CALL get_block_node(matrix=ks%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=h_block)
            ELSE
              CALL get_block_node(matrix=ks%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  block=h_block)
            END IF
!dbg
!  write(*,*) 'iat ', iatom,jatom
!  write(*,*) 'size ', size(h_block,1),size(h_block,2)
!  do i = 1,size(h_block,1)
!    write(*,'(i3,20f10.5)') i,h_block(i,:)
!  enddo
!  stop 'h_b'
            IF(forces) THEN
              IF (iatom <= jatom) THEN
                NULLIFY ( p_block )
                CALL get_block_node(matrix=p%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  block=p_block)
              ELSE
                NULLIFY ( p_block )
                CALL get_block_node(matrix=p%matrix,&
                                    block_row=jatom,&
                                    block_col=iatom,&
                                    block=p_block)
              END IF
            END IF

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
                DO i=1, 3
                   IF ( (cell_c(i) < LBOUND(oce%voce(kkind,iatom)%neighbor,i)) .OR. &
                         (cell_c(i) > UBOUND(oce%voce(kkind,iatom)%neighbor,i))) THEN
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
!dbg
!   write(*,*) 'katom ', katom
                CALL add_vxca_to_ks(h_block,oce,rho_at,ispin, &
                                    iatom,nsgfa,jatom,nsgfb,katom,nsoctot,&
                                    C_coeff_h,C_coeff_s,max_gau)

                IF(forces)THEN
                  IF(iatom <= jatom) THEN
                    rho_at => rho_atom(katom)
                    CALL add_vxca_forces(p_block,oce,rho_at,force(kkind),ispin,&
                                    iatom,nsgfa,jatom,nsgfb,katom,nsoctot,&
                                    C_coeff_h,C_coeff_s,max_gau)
                  ELSE
                    rho_at => rho_atom(katom)
                    CALL add_vxca_forces(p_block,oce,rho_at,force(kkind),ispin,&
                                    jatom,nsgfb,iatom,nsgfa,katom,nsoctot,&
                                    C_coeff_h,C_coeff_s,max_gau)
                  END IF
                END IF

                sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
              ENDDO knodeloop    ! knode
            ENDDO    ! kkind
!dbg
!  write(*,*) 'ks dopo vxc atom'
!  do i = 1,size(h_block,1)
!    write(*,'(i3,20f10.5)') i,h_block(i,:)
!  enddo

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          ENDDO    ! inode
        ENDDO    ! ilist 
      ENDDO    ! jkind
    ENDDO    ! ikind
!  stop 'ks up'
    DEALLOCATE (C_coeff_h)
    DEALLOCATE (C_coeff_s)


  END SUBROUTINE update_ks_vxc_atom
!******************************************************************************

  SUBROUTINE add_vxca_to_ks(h_block,oce,rho_atom,ispin,ia,na,ja,nb,ka,nsp,C_h,C_s,ldc)


    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_block
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(KIND = dp), DIMENSION(ldc,nb)       :: C_h,C_s ! work arrays
    INTEGER, INTENT(IN)                      :: ia,&  ! i atom index 
                                                ispin, & ! spin index
                                                na, & ! # spherical contr. iatom
                                                ja, & ! j atom index
                                                ldc, & ! dimension of work arrays
                                                nb, & ! # spherical contr. jatom
                                                ka, & ! k atom index 
                                                nsp ! # spherical primitives katom

    CHARACTER(LEN=*), PARAMETER :: routine = "add_vxca_to_ks"

    TYPE(rho_atom_coeff), DIMENSION(:),POINTER :: int_h, int_s
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: c_block, IntC_h, IntC_s
!dbg
!    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: h_test
    INTEGER :: i
!******************************************************************************

    NULLIFY(int_h,int_s,intC_h,intC_s)

    CALL get_rho_atom(rho_atom=rho_atom,&
                     ga_Vxc_gb_h=int_h,ga_Vxc_gb_s=int_s)

    IF (ia <= ja) THEN
      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        block=c_block)
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
                        block=c_block)
      C_s = 0.0_dp
      C_s(1:nsp,1:nb)=c_block(1:nsp,1:nb)

! write(*,*) 'C_s_b ',nsp,nb
!  do i = 1,nsp
!    write(*,'(i3,20f10.5)') i,C_s(i,1:nb)
!  enddo


      call reallocate(intC_h,1,nsp,1,nb)
      call reallocate(intC_s,1,nsp,1,nb)

      call dgemm("N","N",nsp,nb,nsp,&
               1.0_dp,int_h(ispin)%r_coef,nsp,&
               C_h,ldc,0.0_dp,&
               intC_h,nsp)

      call dgemm("N","N",nsp,nb,nsp,&
               1.0_dp,int_s(ispin)%r_coef,nsp,&
               C_s,ldc,0.0_dp,&
               intC_s,nsp)

! write(*,*) 'intC_s_b ',nsp,nb
!  do i = 1,nsp
!    write(*,'(i3,20f10.5)') i,intC_s(i,1:nb)
!  enddo

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        block=c_block)
      C_h = 0.0_dp
      C_h(1:nsp,1:na)=c_block(1:nsp,1:na)

!write(*,*) 'C_h_a ',nsp,na
!  do i = 1,nsp
!    write(*,'(i3,20f10.5)') i,C_h(i,1:na)
!  enddo

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        block=c_block)
      C_s = 0.0_dp
      C_s(1:nsp,1:na)=c_block(1:nsp,1:na)
!write(*,*) 'C_s_a ',nsp,na
!  do i = 1,nsp
!    write(*,'(i3,20f10.5)') i,C_s(i,1:na)
!  enddo

      call dgemm("T","N",na,nb,nsp,&
               1.0_dp,C_h,ldc,&
               intC_h,nsp,1.0_dp,&
               h_block,SIZE(h_block,1))
       
      call dgemm("T","N",na,nb,nsp,&
               -1.0_dp,C_s,ldc,&
               intC_s,nsp,1.0_dp,&
               h_block,SIZE(h_block,1))

    ELSE

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        block=c_block)
      C_h = 0.0_dp
      C_h(1:nsp,1:na)=c_block(1:nsp,1:na)

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        block=c_block)
      C_s = 0.0_dp
      C_s(1:nsp,1:na)=c_block(1:nsp,1:na)

      call reallocate(intC_h,1,nsp,1,na)
      call reallocate(intC_s,1,nsp,1,na)

      call dgemm("N","N",nsp,na,nsp,&
               1.0_dp,int_h(ispin)%r_coef,nsp,&
               C_h,ldc,0.0_dp,&
               intC_h,nsp)

      call dgemm("N","N",nsp,na,nsp,&
               1.0_dp,int_s(ispin)%r_coef,nsp,&
               C_s,ldc,0.0_dp,&
               intC_s,nsp)

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        block=c_block)
      C_h = 0.0_dp
      C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

      NULLIFY(c_block)
      CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        block=c_block)
      C_s = 0.0_dp
      C_s(1:nsp,1:nb)=c_block(1:nsp,1:nb)

      call dgemm("T","N",nb,na,nsp,&
               1.0_dp,C_h,ldc,&
               intC_h,nsp,1.0_dp,&
               h_block,SIZE(h_block,1))

      call dgemm("T","N",nb,na,nsp,&
               -1.0_dp,C_s,ldc,&
               intC_s,nsp,1.0_dp,&
               h_block,SIZE(h_block,1))

    END IF

  END SUBROUTINE add_vxca_to_ks

! *****************************************************************************

  SUBROUTINE add_vxca_forces(p_block,oce,rho_atom,force,ispin,&
                            ia,na,ja,nb,ka,nsp,C_h,C_s,ldc)

! *****************************************************************************

    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(KIND = dp), DIMENSION(ldc,nb)        :: C_h,C_s
    TYPE(qs_force_type)                      :: force
    INTEGER, INTENT(IN)                      :: ia,ispin,na,ja,ldc,nb,ka,nsp

    CHARACTER(LEN=*), PARAMETER :: routine = "add_vxca_forces"

    TYPE(rho_atom_coeff), DIMENSION(:), POINTER :: int_h, int_s
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: c_block, PC_h, PC_s, &
                                                dCPC_h, dCPC_s, p_blockt
    REAL(KIND=dp)                            :: trace
    INTEGER                                  :: dir

! *****************************************************************************

    NULLIFY(int_h, int_s, PC_h, PC_s, dCPC_h, dCPC_s)

    CALL get_rho_atom(rho_atom=rho_atom,&
                     ga_Vxc_gb_h=int_h,ga_Vxc_gb_s=int_s)


! i-k contributions

    NULLIFY(c_block)
    CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        block=c_block)
    C_h = 0.0_dp
    C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

    NULLIFY(c_block)
    CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        block=c_block)
    C_s = 0.0_dp
    C_s(1:nsp,1:nb)=c_block(1:nsp,1:nb)

    call reallocate(PC_h,1,na,1,nsp)
    call reallocate(PC_s,1,na,1,nsp)

    call dgemm("N","T",na,nsp,nb,&
               1.0_dp,p_block,SIZE(p_block,1),&
               C_h,ldc,0.0_dp,&
               PC_h,na)

    call dgemm("N","N",na,nsp,nb,&
               1.0_dp,p_block,SIZE(p_block,1),&
               C_s,ldc,0.0_dp,&
               PC_s,na)

    DO dir = 2,4
       NULLIFY(c_block)
       CALL get_block_node(matrix=oce%hardo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        block=c_block)
       C_h = 0.0_dp
       C_h(1:nsp,1:na)=c_block(1:nsp,1:na)

       NULLIFY(c_block)
       CALL get_block_node(matrix=oce%softo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        block=c_block)
       C_s = 0.0_dp
       C_s(1:nsp,1:na)=c_block(1:nsp,1:na)

       call reallocate(dCPC_h,1,nsp,1,nsp)
       call reallocate(dCPC_s,1,nsp,1,nsp)


       call dgemm("N","N",nsp,nsp,na,&
                  1.0_dp,C_h,ldc,&
                  PC_h,na,0.0_dp,&
                  dCPC_h,nsp)

       trace = trace_r_AxB(dCPC_h,nsp,int_h(ispin)%r_coef,nsp,nsp,nsp)

       force%vxc_atom(dir,ka) = force%vxc_atom(dir,ka) + 2.0_dp * trace
       force%vxc_atom(dir,ia) = force%vxc_atom(dir,ia) - 2.0_dp * trace

       call dgemm("N","N",nsp,nsp,na,&
                  1.0_dp,C_s,ldc,&
                  PC_s,na,0.0_dp,&
                  dCPC_s,nsp)

       trace = trace_r_AxB(dCPC_s,nsp,int_s(ispin)%r_coef,nsp,nsp,nsp)

       force%vxc_atom(dir,ka) = force%vxc_atom(dir,ka) - 2.0_dp * trace
       force%vxc_atom(dir,ia) = force%vxc_atom(dir,ia) + 2.0_dp * trace

    END DO

! j-k contributioms

    NULLIFY( PC_h, PC_s, dCPC_h, dCPC_s)

    p_blockt = transpose(p_block)

    NULLIFY(c_block)
    CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        block=c_block)
    C_h = 0.0_dp
    C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

    NULLIFY(c_block)
    CALL get_block_node(matrix=oce%softo(1)%matrix,&
                        block_row=ka,&
                        block_col=ia,&
                        block=c_block)
    C_s = 0.0_dp
    C_s(1:nsp,1:na)=c_block(1:nsp,1:na)

    call reallocate(PC_h,1,nb,1,nsp)
    call reallocate(PC_s,1,nb,1,nsp)

    call dgemm("N","T",nb,nsp,na,&
               1.0_dp,p_blockt,SIZE(p_block,1),&
               C_h,ldc,0.0_dp,&
               PC_h,nb)

    call dgemm("N","N",nb,nsp,na,&
               1.0_dp,p_blockt,SIZE(p_block,1),&
               C_s,ldc,0.0_dp,&
               PC_s,nb)

    DO dir = 2,4
       NULLIFY(c_block)
       CALL get_block_node(matrix=oce%hardo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        block=c_block)
       C_h = 0.0_dp
       C_h(1:nsp,1:nb)=c_block(1:nsp,1:nb)

       NULLIFY(c_block)
       CALL get_block_node(matrix=oce%softo(dir)%matrix,&
                        block_row=ka,&
                        block_col=ja,&
                        block=c_block)
       C_s = 0.0_dp
       C_s(1:nsp,1:nb)=c_block(1:nsp,1:nb)

       call reallocate(dCPC_h,1,nsp,1,nsp)
       call reallocate(dCPC_s,1,nsp,1,nsp)

       call dgemm("N","N",nsp,nsp,nb,&
                  1.0_dp,C_h,ldc,&
                  PC_h,na,0.0_dp,&
                  dCPC_h,nsp)

       trace = trace_r_AxB(dCPC_h,nsp,int_h(ispin)%r_coef,nsp,nsp,nsp)

       force%vxc_atom(dir,ka) = force%vxc_atom(dir,ka) + 2.0_dp * trace
       force%vxc_atom(dir,ja) = force%vxc_atom(dir,ja) - 2.0_dp * trace

       call dgemm("N","N",nsp,nsp,nb,&
                  1.0_dp,C_s,ldc,&
                  PC_s,na,0.0_dp,&
                  dCPC_s,nsp)

       trace = trace_r_AxB(dCPC_s,nsp,int_s(ispin)%r_coef,nsp,nsp,nsp)

       force%vxc_atom(dir,ka) = force%vxc_atom(dir,ka) - 2.0_dp * trace
       force%vxc_atom(dir,ja) = force%vxc_atom(dir,ja) + 2.0_dp * trace

    END DO


  END SUBROUTINE add_vxca_forces

! *****************************************************************************

END MODULE qs_ks_atom
