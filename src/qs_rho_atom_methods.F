!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_rho_atom_methods [1.0] *
!!
!!   NAME
!!     qs_rho_atom_methods
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     
!!
!!   SOURCE
!******************************************************************************

MODULE qs_rho_atom_methods

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             gapw_control_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                             init_lebedev_grids,&
                                             lebedev_grid
  USE mathconstants,                   ONLY: pi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             ncoset,&
                                             nsoset
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_grid_atom,                    ONLY: create_grid_atom,&
                                             grid_atom_type
  USE qs_harmonics_atom,               ONLY: create_harmonics_atom,&
                                             harmonics_atom_type
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
  USE qs_rho_atom_types,               ONLY: allocate_rho_atom_set,&
                                             rho_atom_type
  USE qs_util,                         ONLY: symmetric_add
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                             clebsch_gordon_init
  USE termination,                     ONLY: stop_memory
  USE util,                            ONLY: locate

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "qs_rho_atom_methods"

! *** Public subroutines ***

  PUBLIC :: calculate_rho_atom, &
            calculate_rho_atom_coeff, &
            init_rho_atom

! *****************************************************************************

 CONTAINS

! *****************************************************************************
  SUBROUTINE calculate_rho_atom(qs_env,atom_kind)

    TYPE(qs_environment_type), INTENT(inout) :: qs_env
    TYPE(atomic_kind_type), POINTER          :: atom_kind

    CHARACTER(len=*), PARAMETER :: routine_name = 'calculate_rho_atom', &
      routineP = module_name//':'//routine_name

    INTEGER :: i, i1, i2, iat, iatom, icg, ico1, ico1_first, ico1_last, ico2, &
      ico2_first, ico2_last, ipgf1, ipgf2, iset1, iset2, iso, iso1, &
      iso1_first, iso1_last, iso2, iso2_first, iso2_last, istat, l, lmax12, &
      lmin12, m1, m1s, m2, m2s, max_s_harm, maxcg, maxco, maxl, maxso, n1, &
      n1s, n2, n2s, natom, nr, nset, nspins, size1, size2
    INTEGER, DIMENSION(:), POINTER           :: atom_list, lmax, lmin, npgf
    REAL(dp)                                 :: rtot1, rtot2
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: g1, g2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: dgg, gg
    REAL(dp), DIMENSION(:, :), POINTER       :: coeff, CPCH_sphere, &
                                                CPCS_sphere, zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set

!dbg
!    REAL(dp), DIMENSION(:,:), POINTER     :: sab
!    REAL(dp), DIMENSION(3)     :: rab
!    INTEGER, DIMENSION(:), POINTER        :: nsgf_s
!    INTEGER, DIMENSION(:,:), POINTER     :: first_sgf
!    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ai_work
!    REAL(dp), DIMENSION(:,:), POINTER     :: rpgf
! dbg
! *****************************************************************************

    NULLIFY(rho_atom_set)
    NULLIFY(dft_control)
    NULLIFY(orb_basis)
    NULLIFY(grid_atom)
    NULLIFY(harmonics)
    NULLIFY(atom_list)

    NULLIFY(lmin,lmax,npgf,zet,my_CG,coeff,CPCH_sphere,CPCS_sphere)

    CALL get_qs_env(qs_env=qs_env,rho_atom_set=rho_atom_set,&
                    dft_control=dft_control)

    nspins = dft_control%nspins

    CALL get_atomic_kind(atomic_kind=atom_kind, orb_basis_set=orb_basis,&
                         atom_list=atom_list,harmonics=harmonics,&
                         grid_atom=grid_atom,natom=natom)
! Set density coefficient to zero befor the calculation
    DO iat = 1,natom
      iatom = atom_list(iat)
      DO i = 1,nspins
        rho_atom_set(iatom)%rho_rad_h(i)%r_coef = 0.0_dp
        rho_atom_set(iatom)%rho_rad_s(i)%r_coef = 0.0_dp
        rho_atom_set(iatom)%drho_rad_h(i)%r_coef = 0.0_dp
        rho_atom_set(iatom)%drho_rad_s(i)%r_coef = 0.0_dp
      ENDDO
    ENDDO

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxco=maxco,maxl=maxl,npgf=npgf,nset=nset,zet=zet,&
                          maxso=maxso)

    max_s_harm = harmonics%max_s_harm
    maxcg = harmonics%maxcg
    nr = grid_atom%nr

    my_CG => harmonics%my_CG

    CALL reallocate(CPCH_sphere,1,nsoset(maxl),1,nsoset(maxl))
    CALL reallocate(CPCS_sphere,1,nsoset(maxl),1,nsoset(maxl))

!    ldai = ncoset(maxl + 1)
!    ALLOCATE(ai_work(ldai,ldai,ncoset(1)))
!    CALL reallocate(sab,1,maxco,1,maxco)
!    rab = 0.0_dp
!    dab = 0.0_dp
 
    ALLOCATE(g1(nr),g2(nr),gg(nr,0:2*maxl),dgg(nr,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "g1,g2,gg,dgg")

!******************************************************************************

    g1 = 0.0_dp
    g2 = 0.0_dp
    m1 = 0
    m1s = 0

    DO iset1 = 1,nset
      n1 = ncoset(lmax(iset1))
      n1s = nsoset(lmax(iset1))
      DO ipgf1  = 1,npgf(iset1)
        ico1_first = ncoset(lmin(iset1)-1)+1+n1*(ipgf1-1)+m1
        ico1_last  = ncoset(lmax(iset1))+n1*(ipgf1-1)+m1
        iso1_first = nsoset(lmin(iset1)-1)+1+n1s*(ipgf1-1)+m1s
        iso1_last  = nsoset(lmax(iset1))+n1s*(ipgf1-1)+m1s
        size1 = iso1_last - iso1_first + 1
        ico1 = ncoset(lmin(iset1)-1)+1
        i1 = nsoset(lmin(iset1)-1)+1

        g1(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))
        m2 = 0
        m2s = 0
        DO iset2 = 1,nset
          n2 = ncoset(lmax(iset2))
          n2s=nsoset(lmax(iset2))
          DO ipgf2 = 1,npgf(iset2)
            ico2_first = ncoset(lmin(iset2)-1)+1+n2*(ipgf2-1)+m2
            ico2_last  = ncoset(lmax(iset2))+n2*(ipgf2-1)+m2
            iso2_first = nsoset(lmin(iset2)-1)+1+n2s*(ipgf2-1)+m2s
            iso2_last  = nsoset(lmax(iset2))+n2s*(ipgf2-1)+m2s
            size2 = iso2_last - iso2_first + 1
            ico2 = ncoset(lmin(iset2)-1)+1
            i2 = nsoset(lmin(iset2)-1)+1

            g2(1:nr) = EXP(-zet(ipgf2,iset2)*grid_atom%rad2(1:nr))
            lmin12 = lmin(iset1)+lmin(iset2)
            lmax12 = lmax(iset1)+lmax(iset2)

            gg = 0.0_dp
            dgg = 0.0_dp
            IF (lmin12 == 0) THEN
              gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
            ELSE
              gg(1:nr,lmin12) = grid_atom%rad(1:nr)**lmin12*g1(1:nr)*g2(1:nr)
              dgg(1:nr,lmin12) = REAL(lmin12,dp)/grid_atom%rad(1:nr)*&
                                gg(1:nr,lmin12)
            END IF

            DO l=lmin12+1,lmax12
              dgg(1:nr,l) = REAL(l,dp)*gg(1:nr,l-1)
              gg(1:nr,l) = grid_atom%rad(1:nr)*gg(1:nr,l-1)
              dgg(1:nr,l-1) = dgg(1:nr,l-1) -2.0_dp*(zet(ipgf1,iset1)+ &
                             zet(ipgf2,iset2))*gg(1:nr,l)

            END DO
            dgg(1:nr,lmax12) = dgg(1:nr,lmax12) -2.0_dp*(zet(ipgf1,iset1)+ &
                             zet(ipgf2,iset2))*grid_atom%rad(1:nr)*&
                             gg(1:nr,lmax12)

            DO iat = 1,natom
              iatom = atom_list(iat)

              DO i=1,nspins
                CPCH_sphere = 0.0_dp
                CPCS_sphere = 0.0_dp
                coeff => rho_atom_set(iatom)%cpc_h(i)%r_coef
                CPCH_sphere(i1:i1+size1-1,i2:i2+size2-1) = &
                  coeff(iso1_first:iso1_last,iso2_first:iso2_last)
                coeff => rho_atom_set(iatom)%cpc_s(i)%r_coef
                CPCS_sphere(i1:i1+size1-1,i2:i2+size2-1) = &
                  coeff(iso1_first:iso1_last,iso2_first:iso2_last)


!**********************************************************************************
!dbg
!        do iso1 = nsoset(lmin(iset1)-1)+1,nsoset(lmax(iset1))
!          do iso2 = nsoset(lmin(iset2)-1)+1,nsoset(lmax(iset2))
!             l = indso(1,iso1) + indso(1,iso2)
!!             if(iso1<5 .AND. iso2<5)then
!             do is2 = 1,nr
!               do is1 = 1,grid_atom%ng_sphere
!                  rtot1 = rtot1 + gg(is2,l)*CPCH_sphere(iso1,iso2)*&
!                         harmonics%slm(is1,iso1)*harmonics%slm(is1,iso2)*&
!                         grid_atom%wr(is2)*grid_atom%wa(is1)
!                  rtot2 = rtot2 + gg(is2,l)*CPCS_sphere(iso1,iso2)*&
!                         harmonics%slm(is1,iso1)*harmonics%slm(is1,iso2)*&
!                         grid_atom%wr(is2)*grid_atom%wa(is1)
!               enddo
!             enddo
!!             endif
!          enddo
!        enddo
!***********************************************************************************
                DO iso = 1,max_s_harm
                  DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                    iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                    iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)

                    l = indso(1,iso1) + indso(1,iso2)
                    rho_atom_set(iatom)%rho_rad_h(i)%r_coef(1:nr,iso) =&
                      rho_atom_set(iatom)%rho_rad_h(i)%r_coef(1:nr,iso) +&
                      gg(1:nr,l)*CPCH_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                    rho_atom_set(iatom)%drho_rad_h(i)%r_coef(1:nr,iso) =&
                      rho_atom_set(iatom)%drho_rad_h(i)%r_coef(1:nr,iso) +&
                      dgg(1:nr,l)*CPCH_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)


                    rho_atom_set(iatom)%rho_rad_s(i)%r_coef(1:nr,iso) =&
                      rho_atom_set(iatom)%rho_rad_s(i)%r_coef(1:nr,iso) +&
                      gg(1:nr,l)*CPCS_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                    rho_atom_set(iatom)%drho_rad_s(i)%r_coef(1:nr,iso) =&
                      rho_atom_set(iatom)%drho_rad_s(i)%r_coef(1:nr,iso) +&
                      dgg(1:nr,l)*CPCS_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                  ENDDO ! icg
                ENDDO  ! iso     

              ENDDO ! i
            ENDDO  ! iat
        
          ENDDO ! ipgf2
          m2 = m2+maxco
          m2s = m2s+maxso

        ENDDO  ! iset2
      ENDDO   ! ipgf1
      m1 = m1 + maxco
      m1s = m1s+maxso
    ENDDO    ! iset1
!    write(*,*) 'rho 1 2 A', rtot1,rtot2
    rtot1 = 0.0_dp
    rtot2 = 0.0_dp
!dbg
!********************************************************************************
!           DO iso = 1,max_s_harm
!             do is2 = 1,nr
!               do is1 = 1,grid_atom%ng_sphere
!                  rtot1 = rtot1 + rho_atom_set(3)%rho_rad_h(1)%r_coef(is2,iso)*&
!                         harmonics%slm(is1,iso)*&
!                         grid_atom%wr(is2)*grid_atom%wa(is1)
!                   rtot2 = rtot2 + rho_atom_set(3)%rho_rad_s(1)%r_coef(is2,iso)*&
!                         harmonics%slm(is1,iso)*&
!                         grid_atom%wr(is2)*grid_atom%wa(is1)
!               enddo
!             enddo
!          ENDDO   
!    write(*,*) 'rho 1 2 B', rtot1,rtot2
!********************************************************************************

    DEALLOCATE(g1,g2,gg,dgg,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routine_name,module_name,__LINE__,&
                         "g1,g2,gg,dgg")
    ENDIF

  END SUBROUTINE calculate_rho_atom

! *****************************************************************************

  SUBROUTINE calculate_rho_atom_coeff(globenv,qs_env,rho_ao)

    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao

    CHARACTER(LEN=*), PARAMETER :: routine_name = "calculate_rho_atom_coeff"

    INTEGER :: i, iab, iatom, ibc, icol, ikind, ilist, inode, irow, ispin, &
      jatom, jkind, katom, kkind, kneighbor, knode, max_gau, max_nsgf, &
      maxsoc, natom, nkind, nlist, nnode, nnode_c, nsetc, nsgfa, nsgfb, &
      nsoctot, nspins
    INTEGER, DIMENSION(3)                    :: cell_c
    LOGICAL                                  :: paw_atom
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C_coeff_h_a, C_coeff_h_b, &
                                                C_coeff_s_a, C_coeff_s_b, PC
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: coef_aux_h, coef_aux_s, &
                                                p_coeff
    REAL(KIND=dp), DIMENSION(3)              :: rab, rbc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: c_block, p_block, p_blockt
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb, sac_oce, sbc_oce
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local, &
                                                sbc_oce_neighbor_list
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node, &
                                                sbc_oce_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(qlist_type), DIMENSION(:, :, :), &
      POINTER                                :: tmp_neigh
    TYPE(qlist_type), POINTER                :: sac_oce_neighbor
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

! *****************************************************************************

    NULLIFY(dft_control)
    NULLIFY(atomic_kind_set)
    NULLIFY(oce)
    NULLIFY(sab_orb)
    NULLIFY(sac_oce)
    NULLIFY(sbc_oce)
    NULLIFY(rho_atom)
    NULLIFY(orb_basis_set)
    NULLIFY(c_block,p_block)


    CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            atomic_kind_set=atomic_kind_set,&
            oce=oce,sab_orb=sab_orb,&
            sac_oce=sac_oce,&
            sbc_oce=sbc_oce,&
            rho_atom_set=rho_atom)

    nspins=dft_control%nspins

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=max_nsgf,&
                             maxgtops=max_gau,natom=natom)

    ALLOCATE (C_coeff_h_a(max_gau,max_nsgf))
    ALLOCATE (C_coeff_s_a(max_gau,max_nsgf))
    ALLOCATE (C_coeff_h_b(max_gau,max_nsgf))
    ALLOCATE (C_coeff_s_b(max_gau,max_nsgf))
    ALLOCATE (p_coeff(max_nsgf,max_nsgf,nspins))
    ALLOCATE (PC(max_nsgf,max_gau))
    ALLOCATE (coef_aux_h(max_gau,max_gau,nspins))
    ALLOCATE (coef_aux_s(max_gau,max_gau,nspins))

    DO i = 1,natom
      DO ispin = 1, nspins
        rho_atom(i)%cpc_h(ispin)%r_coef = 0.0_dp
        rho_atom(i)%cpc_s(ispin)%r_coef = 0.0_dp
      ENDDO
    ENDDO
    nkind = SIZE(atomic_kind_set)



    DO ikind = 1 ,nkind
      atomic_kind => atomic_kind_set(ikind)

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set)
      IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
                             nsgf=nsgfa)

      DO jkind=1,nkind

        atomic_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atomic_kind,&
                             orb_basis_set=orb_basis_set)

        IF (.NOT.ASSOCIATED(orb_basis_set)) CYCLE

        CALL get_gto_basis_set(gto_basis_set=orb_basis_set,&
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
                                   neighbor=jatom,&
                                   r=rab)

            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
               DO ispin = 1,nspins
                  NULLIFY(p_block)
                  CALL get_block_node(matrix=rho_ao(ispin)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=p_block)
                  p_coeff(1:nsgfa,1:nsgfb,ispin) = p_block(1:nsgfa,1:nsgfb)
               ENDDO
            ELSE
               irow = jatom
               icol = iatom
               DO ispin = 1,nspins
                  NULLIFY(p_blockt)
                  NULLIFY(p_block)
                  CALL get_block_node(matrix=rho_ao(ispin)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=p_blockt)

                  p_coeff(1:nsgfa,1:nsgfb,ispin)=TRANSPOSE( p_blockt )

               ENDDO
            ENDIF
            
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

              CALL  get_gto_basis_set(gto_basis_set=orb_basis_set,&
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

!               *** Locate operator atom in the sac_oce neighbor list ***

                kneighbor = locate(sac_oce_neighbor%list,katom)

                IF (kneighbor == 0) THEN
                    sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
                    CYCLE
                END IF

                 NULLIFY(c_block)
                 CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=iatom,&
                                    BLOCK=c_block)
                 C_coeff_h_a = 0.0_dp
                 C_coeff_h_a(1:nsoctot,1:nsgfa)=c_block(1:nsoctot,1:nsgfa)

                 NULLIFY(c_block)
                 CALL get_block_node(matrix=oce%softo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=iatom,&
                                    BLOCK=c_block)
                 C_coeff_s_a = 0.0_dp
                 C_coeff_s_a(1:nsoctot,1:nsgfa)=c_block(1:nsoctot,1:nsgfa)

                 IF(iatom == jatom) THEN
!        if(katom == iatom) then
                   DO ispin = 1,nspins
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,p_coeff(1,1,ispin),&
                                SIZE(p_coeff,1),C_coeff_h_a,SIZE(C_coeff_h_a,1),&
                                0.0_dp,PC,max_nsgf)

                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_h_a,&
                                max_gau,PC,max_nsgf,1.0_dp,&
                                rho_atom(katom)%cpc_h(ispin)%r_coef,nsoctot)

                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,p_coeff(1,1,ispin),&
                                SIZE(p_coeff,1),C_coeff_s_a,SIZE(C_coeff_s_a,1),&
                                0.0_dp,PC,max_nsgf)

                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_s_a,&
                                max_gau,PC,max_nsgf,1.0_dp,&
                                rho_atom(katom)%cpc_s(ispin)%r_coef,nsoctot)
                   END DO
!         endif
!      write(*,*) iatom,jatom,katom
!     do i = 1,nsgfa
!     write(*,'(i3,14f10.5)') i,p_coeff(i,1:nsgfb,1)
!     enddo
!!
!     do i = 1,nsoctot
!     write(*,'(i3,14f10.5)') i,C_coeff_h_a(i,1:nsgfa)
!     enddo
!      if(katom ==1 ) then
!      write(*,*) iatom,jatom
!     do i = 1,nsoctot
!     write(*,'(i3,32f8.3)') i,rho_atom(katom)%cpc_h(1)%r_coef(i,1:nsoctot)
!     enddo
!    endif
!      stop 'katom'
                 ELSE
                    NULLIFY(c_block)
                    CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=jatom,&
                                    BLOCK=c_block)
                    C_coeff_h_b = 0.0_dp
                    C_coeff_h_b(1:nsoctot,1:nsgfb)=c_block(1:nsoctot,1:nsgfb)


                    NULLIFY(c_block)
                    CALL get_block_node(matrix=oce%softo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=jatom,&
                                    BLOCK=c_block)
                    C_coeff_s_b = 0.0_dp
                    C_coeff_s_b(1:nsoctot,1:nsgfb)=c_block(1:nsoctot,1:nsgfb)


                    DO ispin = 1,nspins
                      CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,&
                                 1.0_dp,p_coeff(1,1,ispin),SIZE(p_coeff,1),&
                                 C_coeff_h_b,SIZE(C_coeff_h_b,1),&
                                 0.0_dp,PC,max_nsgf)


                      CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 1.0_dp,C_coeff_h_a,max_gau,&
                                 PC,max_nsgf,0.0_dp,&
                                 coef_aux_h(1,1,ispin),max_gau)

                      CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                                 p_coeff(1,1,ispin),SIZE(p_coeff,1),&
                                 C_coeff_s_b,SIZE(C_coeff_h_b,1),&
                                 0.0_dp,PC,max_nsgf)

                      CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,&
                                 C_coeff_s_a,max_gau,PC,&
                                 max_nsgf,0.0_dp,&
                                 coef_aux_s(1,1,ispin),max_gau)

                      CALL symmetric_add(rho_atom(katom)%cpc_h(ispin)%r_coef,nsoctot,&
                                       coef_aux_h(1,1,ispin),max_gau,nsoctot)
                      CALL symmetric_add(rho_atom(katom)%cpc_s(ispin)%r_coef,nsoctot,&
                                       coef_aux_s(1,1,ispin),max_gau,nsoctot)
                    END DO
!     do i = 1,nsoctot
!     write(*,'(i3,32f8.3)') i,rho_atom(katom)%cpc_h(1)%r_coef(i,1:nsoctot)
!     enddo
                 ENDIF
10   CONTINUE
                 sbc_oce_neighbor_node => next(sbc_oce_neighbor_node)
              ENDDO  knodeloop ! knode
            ENDDO  ! kkind
          sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          ENDDO  ! inode
        ENDDO  ! ilist
      ENDDO  ! jkind
    ENDDO ! ikind

!     do i = 1,32
!     write(*,'(i3,16f10.5)') i,rho_atom(1)%cpc_h(1)%r_coef(i,1:16)
!     enddo
!  stop 'coef'
#ifdef PARALLEL
    DO iatom = 1,natom
       DO ispin = 1,nspins
         mp_sum(rho_atom(iatom)%cpc_h(ispin)%r_coef,globenv%group)
         mp_sum(rho_atom(iatom)%cpc_s(ispin)%r_coef,globenv%group)
       ENDDO
    END DO
#endif

    DEALLOCATE (C_coeff_h_a)
    DEALLOCATE (C_coeff_s_a)
    DEALLOCATE (C_coeff_h_b)
    DEALLOCATE (C_coeff_s_b)
    DEALLOCATE (p_coeff)
    DEALLOCATE (PC)
    DEALLOCATE (coef_aux_h)
    DEALLOCATE (coef_aux_s)

  END SUBROUTINE calculate_rho_atom_coeff

! *****************************************************************************

  SUBROUTINE init_rho_atom(qs_env,globenv)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type)            :: globenv

    CHARACTER(LEN=*), PARAMETER              :: routine_name = "init_rho_atom"

    INTEGER :: i, ia, iat, iatom, ikind, il, iso, iso1, iso2, l1, l1l2, l2, &
      lc1, lc2, ll, llmax, lmax_sphere, lp, m1, m2, max_s_harm, max_s_set, &
      maxco, maxlgto, maxso, mm, mp, na, nat, natom, ncotot, nkind, nr, nset, &
      nsotot, nspins, quadrature
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: paw_atom
    REAL(dp)                                 :: cosia
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rga
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gapw_control_type), POINTER         :: gapwc
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set

! *****************************************************************************

    NULLIFY(rho_atom_set)
    NULLIFY(atomic_kind_set)
    NULLIFY(atomic_kind)
    NULLIFY(dft_control)
    NULLIFY(gapwc)
    NULLIFY(orb_basis_set)
    NULLIFY(my_CG, atomic_kind, grid_atom, harmonics,atom_list)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                    dft_control=dft_control) 

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             natom = natom)

    nspins = dft_control%nspins

    CALL allocate_rho_atom_set(rho_atom_set,natom)

    nkind = SIZE(atomic_kind_set)
    gapwc => dft_control%qs_control%gapw_control
    lmax_sphere = gapwc%lmax_sphere

    llmax = MAX(lmax_sphere, 2*maxlgto)
    max_s_harm = nsoset(llmax)
    max_s_set = nsoset(maxlgto)

!   *** allocate calculate the CG coefficients up to the maxl *** 
    CALL clebsch_gordon_init ( llmax )
    CALL reallocate(my_CG,1,max_s_set,1,max_s_set,1,max_s_harm)

    ALLOCATE(rga(llmax,2))
    DO lc1 = 0,maxlgto
      DO  iso1=nsoset(lc1-1)+1,nsoset(lc1)
        l1 = indso(1,iso1)
        m1 = indso(2,iso1)
        DO lc2 = 0,maxlgto
          DO iso2 = nsoset(lc2-1)+1,nsoset(lc2) 
            l2 = indso(1,iso2)
            m2 = indso(2,iso2)
            CALL clebsch_gordon (l1,m1,l2,m2,rga)
            IF(l1+l2 > llmax) THEN
              l1l2 = llmax
            ELSE
              l1l2 = l1+l2
            END IF
            mp = m1 + m2
            mm = m1 - m2
             IF ( m1*m2 < 0  .OR. (m1*m2==0 .AND. (m1<0 .OR. m2<0))) THEN
               mp = -ABS(mp)
               mm = -ABS(mm)
             ELSE
               mp = ABS(mp)
               mm = ABS(mm)
             END IF
            DO  lp=MOD(l1+l2,2),l1l2,2
               il = lp/2 + 1
               IF( ABS(mp) <= lp) THEN
               IF(mp >= 0) THEN
                 iso = nsoset(lp-1) + lp + 1 + mp
               ELSE
                 iso = nsoset(lp-1) + lp + 1 - ABS(mp)
               END IF
               my_CG(iso1,iso2,iso) =  rga(il,1)
               ENDIF
               IF(mp /= mm .AND. ABS(mm) <= lp) THEN
               IF(mm >= 0) THEN
                 iso = nsoset(lp-1) + lp + 1 + mm
               ELSE
                  iso = nsoset(lp-1) + lp + 1 - ABS(mm)
               END IF
               my_CG(iso1,iso2,iso) =  rga(il,2)
               ENDIF
             END DO
          ENDDO  ! iso2
        ENDDO  ! lc2
      ENDDO    ! iso1
    ENDDO      ! lc1
    DEALLOCATE(rga)
 
!   *** initialize the Lebedev grids ***
    CALL init_lebedev_grids
    quadrature=gapwc%quadrature

    DO ikind = 1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set, &
                           atom_list=atom_list,natom=nat, &
                           paw_atom=paw_atom,&
                           grid_atom=grid_atom,&
                           harmonics=harmonics,&
                           ngrid_rad=nr,ngrid_ang=na)

      IF(.NOT. paw_atom) CYCLE

!     *** determine the Lebedev grid for this kind ***
      ll =  get_number_of_lebedev_grid(n=na)
      na = lebedev_grid(ll)%n
      grid_atom%ng_sphere = na
      grid_atom%nr = nr

!     *** calculate the radial grid ***
      CALL create_grid_atom(grid_atom,nr,na,quadrature,globenv)

      grid_atom%wa(1:na) = 4._dp*pi*lebedev_grid(ll)%w(1:na)

      DO ia = 1,na
        cosia = lebedev_grid(ll)%r(3,ia)
        grid_atom%cos_teta(ia) = cosia 
        grid_atom%sin_teta(ia) =  SQRT(1-cosia*cosia)
      ENDDO

!     *** calculate the spherical harmonics on the grid ***
      CALL create_harmonics_atom(harmonics,orb_basis_set,&
                   my_CG,na,llmax,max_s_harm,ll,globenv)


      CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                             maxco=maxco,maxso=maxso,nset=nset)
      nsotot = maxso * nset 
      ncotot = maxco * nset
      DO iat = 1,nat
        iatom = atom_list(iat) 
!       *** allocate the radial density for each LM,for each atom ***
 
        ALLOCATE (rho_atom_set(iatom)%cpc_h(nspins))
        ALLOCATE (rho_atom_set(iatom)%cpc_s(nspins))
        ALLOCATE (rho_atom_set(iatom)%rho_rad_h(nspins))
        ALLOCATE (rho_atom_set(iatom)%rho_rad_s(nspins))
        ALLOCATE (rho_atom_set(iatom)%drho_rad_h(nspins))
        ALLOCATE (rho_atom_set(iatom)%drho_rad_s(nspins))
        ALLOCATE (rho_atom_set(iatom)%ga_Vxc_gb_h(nspins))
        ALLOCATE (rho_atom_set(iatom)%ga_Vxc_gb_s(nspins))

        DO i = 1,nspins
          CALL reallocate(rho_atom_set(iatom)%cpc_h(i)%r_coef,&
                                      1,nsotot,1,nsotot)
          CALL reallocate(rho_atom_set(iatom)%cpc_s(i)%r_coef,&
                                      1,nsotot,1,nsotot)
          ALLOCATE(rho_atom_set(iatom)%rho_rad_h(i)%r_coef(&
                                      1:nr,1:max_s_harm))
          ALLOCATE(rho_atom_set(iatom)%rho_rad_s(i)%r_coef(& 
                                      1:nr,1:max_s_harm))

          CALL reallocate(rho_atom_set(iatom)%drho_rad_h(i)%r_coef,&
                                      1,nr,1,max_s_harm)
          CALL reallocate(rho_atom_set(iatom)%drho_rad_s(i)%r_coef,& 
                                      1,nr,1,max_s_harm)

          CALL reallocate(rho_atom_set(iatom)%ga_Vxc_gb_h(i)%r_coef,&
                                      1,nsotot,1,nsotot)
          CALL reallocate(rho_atom_set(iatom)%ga_Vxc_gb_s(i)%r_coef,&
                                      1,nsotot,1,nsotot)
        END DO

      END DO ! iat
    END DO

    CALL set_qs_env(qs_env=qs_env, rho_atom_set=rho_atom_set)

  END SUBROUTINE init_rho_atom

! *****************************************************************************

END MODULE qs_rho_atom_methods
