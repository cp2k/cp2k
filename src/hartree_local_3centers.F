!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/hartree_local_3centers
!!
!!   NAME
!!     hartree_local_3centers
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE  hartree_local_3centers

  USE ai_coulomb,                      ONLY: coulomb3
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,dp_size,&
                                             int_size
  USE hartree_local_types,             ONLY: ecoul_1center_type,&
                                             get_ecoul_1c
  USE mathconstants,                   ONLY: dfac,fourpi,pi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: coset,&
                                             indco,&
                                             nco,&
                                             ncoset,&
                                             nso,&
                                             nsoset
  USE orbital_transformation_matrices, ONLY: orbtramat
  USE particle_types,                  ONLY: particle_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type
  USE qs_neighbor_list_types,          ONLY: find_neighbor_list,&
                                             first_list,&
                                             first_node,&
                                             get_neighbor_list,&
                                             get_neighbor_list_set,&
                                             get_neighbor_node,&
                                             neighbor_list_p_type,&
                                             neighbor_list_set_p_type,&
                                             neighbor_list_set_type,&
                                             neighbor_list_type,&
                                             neighbor_node_type,&
                                             next,&
                                             qlist_type
  USE qs_oce_types,                    ONLY: vtriple_type
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
  USE qs_util,                         ONLY: trace_r_AxB
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: locate

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "hartree_local_3centers"

! Public Subroutine
 
  PUBLIC :: integrate_a_vhg0_b

!******************************************************************************
 
  CONTAINS
 
!******************************************************************************
  SUBROUTINE integrate_a_vhg0_b(qs_env,h,p,energy,just_energy,&
                                calculate_forces,error)

    TYPE(qs_energy_type), INTENT(INOUT)        :: energy
    TYPE(real_matrix_p_type), INTENT(IN)       :: p
    TYPE(real_matrix_p_type), INTENT(INOUT)    :: h
    TYPE(qs_environment_type), POINTER         :: qs_env
    LOGICAL, INTENT(IN)                        :: calculate_forces,&
                                                  just_energy
    TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine_name = "integrate_a_vhg0_b"

    TYPE(atomic_kind_type), DIMENSION(:),&
                                  POINTER      :: atomic_kind_set
    TYPE(atomic_kind_type),       POINTER      :: atom_kind
    TYPE(cell_type),              POINTER      :: cell
    TYPE(ecoul_1center_type), DIMENSION(:), & 
                                 POINTER       :: ecoul_1c
    TYPE(gto_basis_set_type),     POINTER      :: orb_basis
    TYPE(neighbor_list_set_p_type), &
                         DIMENSION(:), POINTER :: sab_orb,sbc_3c
    TYPE(neighbor_list_type),     POINTER      :: sab_orb_neighbor_list,&
                                                  sab_orb_neighbor_list_local,&
                                                  sbc_3c_neighbor_list
    TYPE(neighbor_node_type),     POINTER      :: sab_orb_neighbor_node,&
                                                  sbc_3c_neighbor_node
    TYPE(cp_para_env_type), POINTER            :: para_env
    TYPE(particle_type), DIMENSION(:), &
                                     POINTER   :: particle_set
    TYPE(qlist_type), DIMENSION(:, :, :), &
                                       POINTER :: tmp_neighbor
    TYPE(qlist_type),             POINTER      :: sac_3c_neighbor 
    TYPE(qs_force_type), DIMENSION(:), POINTER :: force
    TYPE(rho0_mpole_type),        POINTER      :: rho0_mpole 
    TYPE(vtriple_type), DIMENSION(:,:),&
                                  POINTER      :: VH_3c_list
    REAL(dp) :: alpha_core, c1, dab, ieqj, ra(3), rab(3), rab2, rac(3), rac2, rb(3), &
                rbc(3), rbc2,  rc(3), rpgf0_h, rpgf0_s, trace, z0_h, z0_s, zhard 
    REAL(dp), DIMENSION(:),       POINTER      :: ng0_h, ng0_s
    REAL(dp), DIMENSION(:),       POINTER      :: normhard
    REAL(dp), DIMENSION(:),   ALLOCATABLE      :: sum_ab_sph
    REAL(dp), DIMENSION(:),       POINTER      :: QN_h, QN_s,Qlm
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE    :: sum_ab_car
    REAL(dp), DIMENSION(:,:),     POINTER      :: h_block,p_block
    REAL(dp), DIMENSION(:,:),     POINTER      :: rpgfa,rpgfb,&
                                                  sphi_a,sphi_b,&
                                                  zeta,zetb
    REAL(dp), DIMENSION(:,:),     POINTER      :: CPC_ab
    REAL(dp), DIMENSION(:,:),     POINTER      :: hab, work
    REAL(dp), DIMENSION(:,:),     POINTER      :: Int_ab_sum
    REAL(dp), DIMENSION(:,:),     POINTER      :: aVh3c0b
    REAL(dp), DIMENSION(:,:,:),   POINTER      :: Qlm_gg
    REAL(dp), DIMENSION(:,:,:),   POINTER      :: Int_abc
    REAL(dp), DIMENSION(:,:,:), POINTER        :: Int_dab_sum, Int_adb_sum

    INTEGER :: ab, bc, handle, iatom,  i_atom, i, ic, ico, icol, idir, ifs1, ifs2, ikind,&
               ilist, inode, irow, iset1, iset2, is, iso, istat, jatom, j_atom, jkind, &
               katom, kat, k_atom, kkind, kneighbor, knode, l0_kkind, last_jatom, lmax_0, &
               lshell, lx, ly, lz, m1, m2, maxch0, maxco_global, maxcoa, maxcob, maxder, &
               maxsoc, maxsgf, maxsgf_set, nkind, nat, natom, nchan0_car,&
               nchan0_sph, nco1, nco2, ncotota, ncototb, nlist, nnode, nnode_c,&
               ns1, ns2, nsgfa, nsgfb, nseta, nsetb, nsetc, nsototc, nz1, nz2
    INTEGER, DIMENSION(3)                      :: cell_c
    INTEGER, DIMENSION(:),        POINTER      :: atom_list,la_max,la_min,&
                                                  lb_max,lb_min,&
                                                  npgfa,npgfb,nsgfa_set,nsgfb_set
    INTEGER, DIMENSION(:),        ALLOCATABLE  :: atom_of_kind
    INTEGER, DIMENSION(:,:),      POINTER      :: first_sgfa,first_sgfb
    LOGICAL                                    :: new_atom_b, paw_atom
    INTEGER  :: ic1, ic2, ip1,ip2, n1, n2, nn1, nn2,jc
    REAL(dp) ::  e3c
!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind_set,cell,ecoul_1c)
    NULLIFY( para_env,particle_set,rho0_mpole)
    NULLIFY(sab_orb,sbc_3c,VH_3c_list)

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                    cell=cell,ecoul_1c=ecoul_1c,&
                    para_env=para_env,& 
                    particle_set=particle_set,&
                    rho0_mpole=rho0_mpole,&
                    VH_3c_list=VH_3c_list,&
                    sab_orb=sab_orb,sbc_3c=sbc_3c)

    nkind = SIZE(atomic_kind_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco_global,&
                             maxsgf=maxsgf, maxsgf_set=maxsgf_set,&
                             natom=natom)

    maxder = 0
    IF (calculate_forces) THEN
 
      CALL timeset(routine_name//" (forces)","I"," ",handle)

      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                       "atom_of_kind",natom*int_size)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      CALL get_qs_env(qs_env=qs_env, force=force)

      maxder = 1
    ELSE
     
      CALL timeset(routine_name,"I"," ",handle)

    END IF

    NULLIFY(ng0_h,ng0_s)

    CALL get_rho0_mpole(rho0_mpole=rho0_mpole,zet0_h=z0_h,zet0_s=z0_s,&
                        norm_g0l_h=ng0_h,norm_g0l_s=ng0_s,lmax_0=lmax_0)

    maxch0 = ncoset(lmax_0)
    NULLIFY(normhard)
    CALL reallocate(normhard,0,lmax_0)

    ALLOCATE(sum_ab_car(maxch0,natom),sum_ab_sph(maxch0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                           "sum_ab_car,sum_ab_sph",natom*maxch0+maxch0)
    sum_ab_car = 0.0_dp

    DO ikind=1,nkind

      NULLIFY(atom_kind,orb_basis)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,soft_basis_set=orb_basis)

      IF (.NOT.ASSOCIATED(orb_basis)) CYCLE

      NULLIFY(first_sgfa,la_max,la_min,npgfa,nsgfa_set,rpgfa,sphi_a,zeta)
      CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             first_sgf=first_sgfa,&
                             lmax=la_max,&
                             lmin=la_min,&
                             maxco=maxcoa,&
                             npgf=npgfa,&
                             nset=nseta,&
                             nsgf=nsgfa,&
                             nsgf_set=nsgfa_set,&
                             pgf_radius=rpgfa,&
                             sphi=sphi_a,&
                             zet=zeta)
      ncotota = maxcoa * nseta

      DO jkind = 1,nkind
        NULLIFY(atom_kind,orb_basis)
        atom_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atom_kind,soft_basis_set=orb_basis)

        IF (.NOT.ASSOCIATED(orb_basis)) CYCLE
        NULLIFY(first_sgfb,lb_max,lb_min,npgfb,nsgfb_set,rpgfb,sphi_b,zetb)
        CALL get_gto_basis_set(gto_basis_set=orb_basis,&
                             first_sgf=first_sgfb,&
                             lmax=lb_max,&
                             lmin=lb_min,&
                             maxco=maxcob,&
                             npgf=npgfb,&
                             nset=nsetb,&
                             nsgf=nsgfb,&
                             nsgf_set=nsgfb_set,&
                             pgf_radius=rpgfb,&
                             sphi=sphi_b,&
                             zet=zetb)
        ncototb = maxcob * nsetb

        ab = ikind + nkind*(jkind - 1)

        NULLIFY(sab_orb_neighbor_list)
        IF (ASSOCIATED(sab_orb(ab)%neighbor_list_set)) THEN

           CALL get_neighbor_list_set(neighbor_list_set=&
                                      sab_orb(ab)%neighbor_list_set,&
                                      nlist=nlist)
           sab_orb_neighbor_list => first_list(sab_orb(ab)%neighbor_list_set)
        ELSE
           nlist=0
        END IF
        
        NULLIFY(work)
        CALL reallocate(work,1,maxco_global,1,maxsgf_set)
        DO ilist = 1,nlist

          sab_orb_neighbor_list_local => sab_orb_neighbor_list

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom, nnode=nnode)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          last_jatom = 0
          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

!DBG
!      WRITE(*,*) 'pair atoms', iatom, jatom
            rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = sqrt(rab2)

            rb(:) = pbc(particle_set(jatom)%r,cell)

            IF(iatom == jatom) THEN
              ieqj = 1.0_dp
            ELSE
              ieqj = 2.0_dp
            END IF 

            NULLIFY(CPC_ab,hab,h_block, p_block)
            NULLIFY(Int_ab_sum,Int_dab_sum,Int_adb_sum,Int_abc)
            IF (iatom <= jatom) THEN
!             Retrieve KS and density matrix for this couple of atoms
              CALL get_block_node(matrix=h%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=h_block)

              CALL get_block_node(matrix=p%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=p_block)

              CALL reallocate(CPC_ab,1,ncotota,1,ncototb)
              CALL reallocate(hab,1,maxcoa,1,maxcob)

!             Build the CPC coefficients that should then multiplied to the 
!             3-centers integrals of primitives to get the energy and the aVh3c0b
!             integrals for atom C (the one where the g0 are centered)

              m1 = 0
              DO iset1 = 1,nseta
                ifs1 = first_sgfa(1,iset1)
                nco1 = npgfa(iset1)*ncoset(la_max(iset1))
                ns1  = nsgfa_set(iset1)
                m2 = 0
                DO iset2 = 1,nsetb
                  ifs2 = first_sgfb(1,iset2)
                  nco2 = npgfb(iset2)*ncoset(lb_max(iset2))
                  ns2  = nsgfb_set(iset2)

                 hab = 0.0_dp
                  CALL dgemm("N","N",nco1,ns2,ns1,&
                              1.0_dp,sphi_a(1,ifs1),maxcoa,&
                              p_block(ifs1,ifs2),nsgfa,&
                              0.0_dp,work(1,1),maxco_global)
                  CALL dgemm("N","T",nco1,nco2,ns2,&
                             1.0_dp,work(1,1),maxco_global,&
                             sphi_b(1,ifs2), maxcob,&
                             0.0_dp,hab(1,1),maxcoa)

                 DO ip1 = 1,npgfa(iset1)
                   ic1 = ncoset(la_min(iset1)-1)+1
                   n1 = ncoset(la_max(iset1))*(ip1-1)
                   nn1 = ncoset(la_max(iset1))*ip1
                   DO ip2 = 1,npgfb(iset2)
                     ic2 = ncoset(lb_min(iset2)-1)+1
                     n2 =  ncoset(lb_max(iset2))*(ip2-1)
                     nn2 = ncoset(lb_max(iset2))*ip2
                     CPC_ab(ic1+n1+m1:nn1+m1,ic2+n2+m2:nn2+m2) = hab(ic1+n1:nn1,ic2+n2:nn2)
                   END DO
                 END DO

                 m2 = m2 + maxcob
                END DO  ! iset2
                m1 = m1 + maxcoa
              END DO  ! iset1

!             Allocate array for sum of integrals
              CALL reallocate(Int_ab_sum,1,maxcoa,1,maxcob)
              IF(calculate_forces) THEN
                CALL reallocate(Int_dab_sum,1,maxcoa,1,maxcob,1,3)
                CALL reallocate(Int_adb_sum,1,maxcoa,1,maxcob,1,3)
              END IF
!             Allocate array for the 3-centers integrals
              CALL reallocate(Int_abc,1,ncotota,1,ncototb,1,maxch0)

            ELSE
              CALL get_block_node(matrix=h%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  BLOCK=h_block)

              CALL get_block_node(matrix=p%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  BLOCK=p_block)

              CALL reallocate(CPC_ab,1,ncototb,1,ncotota)
              CALL reallocate(hab,1,maxcob,1,maxcoa)

              m1 = 0
              DO iset1 = 1,nsetb
                ifs1 = first_sgfb(1,iset1)
                nco1 = npgfb(iset1)*ncoset(lb_max(iset1))
                ns1  = nsgfb_set(iset1)
                m2 = 0
                DO iset2 = 1,nseta
                  ifs2 = first_sgfa(1,iset2)
                  nco2 = npgfa(iset2)*ncoset(la_max(iset2))
                  ns2  = nsgfa_set(iset2)
                  hab = 0.0_dp
                  CALL dgemm("N","N",nco1,ns2,ns1,&
                              1.0_dp,sphi_b(1,ifs1),maxcob,&
                              p_block(ifs1,ifs2),nsgfb,&
                              0.0_dp,work(1,1),maxco_global)
                  CALL dgemm("N","T",nco1,nco2,ns2,&
                             1.0_dp,work(1,1),maxco_global,&
                             sphi_a(1,ifs2), maxcoa,&
                             0.0_dp,hab(1,1),maxcob)

                 DO ip1 = 1,npgfb(iset1)
                   ic1 = ncoset(lb_min(iset1)-1)+1
                   n1 = ncoset(lb_max(iset1))*(ip1-1)
                   nn1 = ncoset(lb_max(iset1))*ip1
                   DO ip2 = 1,npgfa(iset2)
                     ic2 = ncoset(la_min(iset2)-1)+1
                     n2 =  ncoset(la_max(iset2))*(ip2-1)
                     nn2 = ncoset(la_max(iset2))*ip2
                     CPC_ab(ic1+n1+m1:nn1+m1,ic2+n2+m2:nn2+m2) = hab(ic1+n1:nn1,ic2+n2:nn2)
                   END DO
                 END DO
                 m2 = m2 + maxcoa
                END DO  ! iset2
                m1 = m1 + maxcob
              END DO  ! iset1
!
!             Allocate array for sum of integrals
              CALL reallocate(Int_ab_sum,1,maxcob,1,maxcoa)
              IF(calculate_forces) THEN
                CALL reallocate(Int_dab_sum,1,maxcob,1,maxcoa,1,3)
                CALL reallocate(Int_adb_sum,1,maxcob,1,maxcoa,1,3)
              END IF
!             Allocate array for the 3-centers integrals
              CALL reallocate(Int_abc,1,ncototb,1,ncotota,1,maxch0)
            END IF

            DO kkind=1,nkind
              IF(.NOT.ASSOCIATED(VH_3c_list(kkind,iatom)%neighbor)) CYCLE

              bc = jkind + nkind*(kkind - 1)

              IF (.NOT.ASSOCIATED(sbc_3c(bc)%neighbor_list_set)) CYCLE
              NULLIFY(sbc_3c_neighbor_list)

              sbc_3c_neighbor_list =>&
                         find_neighbor_list(neighbor_list_set=&
                         sbc_3c(bc)%neighbor_list_set,&
                         atom=jatom)

              IF (.NOT.ASSOCIATED(sbc_3c_neighbor_list)) CYCLE

              CALL get_rho0_mpole(rho0_mpole=rho0_mpole,ikind=kkind,&
                                  l0_ikind=l0_kkind,&
                                  rpgf0_h=rpgf0_h,rpgf0_s=rpgf0_s)

              nchan0_car = ncoset(l0_kkind)

              NULLIFY(atom_kind,orb_basis)
              atom_kind => atomic_kind_set(kkind)
              CALL get_atomic_kind(atomic_kind=atom_kind,&
                                   paw_atom=paw_atom,&
                                   alpha_core_charge=alpha_core)

              CALL get_neighbor_list(neighbor_list=sbc_3c_neighbor_list,&
                                     nnode=nnode_c)
              IF(paw_atom)THEN
                normhard(0:l0_kkind) = ng0_h(0:l0_kkind)
                zhard = z0_h
              ELSE
                c1 = alpha_core/pi
                c1 = c1*c1*c1
                normhard(0) = SQRT(c1) 
                zhard =  alpha_core
              END IF

              sbc_3c_neighbor_node => first_node(sbc_3c_neighbor_list)

              knodeloop: DO knode=1,nnode_c

                CALL get_neighbor_node(neighbor_node=sbc_3c_neighbor_node,&
                          neighbor=katom,&
                          cell=cell_c,&
                          r=rbc)

                rbc2= rbc(1)*rbc(1) + rbc(2)*rbc(2) + rbc(3)*rbc(3) 
                rc(:) = pbc(particle_set(katom)%r,cell)
                rac(:) = rab(:) + rbc(:)
                rac2= rac(1)*rac(1) + rac(2)*rac(2) + rac(3)*rac(3) 

                tmp_neighbor => VH_3c_list(kkind,iatom)%neighbor
                DO i=1, 3
                   IF ( (cell_c(i) < LBOUND(tmp_neighbor,i)) .OR. &
                         (cell_c(i) > UBOUND(tmp_neighbor,i))) THEN
                      sbc_3c_neighbor_node => next(sbc_3c_neighbor_node)
                      CYCLE knodeloop
                   END IF
                END DO
                sac_3c_neighbor => VH_3c_list(kkind,iatom)%neighbor(cell_c(1),&
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

                Int_abc = 0.0_dp
                NULLIFY(QN_s,QN_h,Qlm)
                CALL  get_rho0_mpole(rho0_mpole=rho0_mpole, iat=katom,&
                                     Qlm_car_s=QN_s, Qlm_car_h=QN_h, Qlm_tot=Qlm)

                IF(.NOT. paw_atom) THEN
                  QN_h(1) = QN_h(1)*normhard(0)/ng0_h(0)
                ENDIF

                IF(iatom <= jatom) THEN

                  m1 = 0
                  DO iset1 = 1,nseta
                    ifs1 = first_sgfa(1,iset1)
                    nz1  = npgfa(iset1)
                    nco1 = npgfa(iset1)*ncoset(la_max(iset1))
                    ns1  = nsgfa_set(iset1)
                    m2 = 0
                    DO iset2 = 1,nsetb
                      ifs2 = first_sgfb(1,iset2)
                      nz2  = npgfb(iset2)
                      nco2 = npgfb(iset2)*ncoset(lb_max(iset2))
                      ns2  = nsgfb_set(iset2)

                      Int_ab_sum = 0.0_dp
                      CALL sumup_3c_integrals(&
                           la_max(iset1),la_min(iset1),nco1,m1,&
                           rpgfa(:,iset1),zeta(:,iset1),nz1,&
                           lb_max(iset2),lb_min(iset2),nco2,m2,&
                           rpgfb(:,iset2),zetb(:,iset2),nz2,&
                           l0_kkind,normhard(0:),ng0_s(0:),rpgf0_h,rpgf0_s,&
                           QN_h,-QN_s,Qlm,zhard,z0_s,&
                           rab,rab2,rac,rac2,rbc,rbc2,&
                           Int_ab_sum,Int_abc,maxder,&
                           Int_dab_sum,Int_adb_sum)

                      CALL dgemm("N","N",nco1,ns2,nco2,&
                           1.0_dp,Int_ab_sum(1,1),maxcoa,&
                           sphi_b(1,ifs2),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                      CALL dgemm("T","N",ns1,ns2,nco1,&
                            1.0_dp,sphi_a(1,ifs1),SIZE(sphi_a,1),&
                           work(1,1),SIZE(work,1),&
                           1.0_dp,h_block(ifs1,ifs2),SIZE(h_block,1))

                      hab(1:nco1,1:nco2) = CPC_ab(1+m1:m1+nco1,1+m2:m2+nco2)  

                      trace = trace_r_AxB(hab,maxcoa,Int_ab_sum,maxcoa,nco1,nco2)
                      energy%hartree_3c = energy%hartree_3c + ieqj*trace 

                      IF(calculate_forces) THEN
                        i_atom = atom_of_kind(iatom)
                        j_atom = atom_of_kind(jatom)
                        k_atom = atom_of_kind(katom)
                        DO idir = 1,3 
                          trace = trace_r_AxB(hab,maxcoa,&
                                  Int_dab_sum(:,:,idir),maxcoa,nco1,nco2)
                          force(ikind)%hartree_3c(idir,i_atom) = &
                                  force(ikind)%hartree_3c(idir,i_atom) - ieqj*trace
                          force(kkind)%hartree_3c(idir,k_atom) = &
                                  force(kkind)%hartree_3c(idir,k_atom) + ieqj*trace

                          trace = trace_r_AxB(hab,maxcoa,&
                                  Int_adb_sum(:,:,idir),maxcoa,nco1,nco2)
                          force(jkind)%hartree_3c(idir,j_atom) = &
                                  force(jkind)%hartree_3c(idir,j_atom) - ieqj*trace
                          force(kkind)%hartree_3c(idir,k_atom) = &
                                  force(kkind)%hartree_3c(idir,k_atom) + ieqj*trace
                        END DO 
                      END IF

                      m2 = m2 + maxcob
                    END DO  ! iset2
                    m1 = m1 + maxcoa
                  END DO  ! iset1

                ELSE

                  m1 = 0
                  DO iset1 = 1,nsetb
                    ifs1 = first_sgfb(1,iset1)
                    nz1  = npgfb(iset1)
                    nco1 = npgfb(iset1)*ncoset(lb_max(iset1))
                    ns1  = nsgfb_set(iset1)
                    m2 = 0
                    DO iset2 = 1,nseta
                      ifs2 = first_sgfa(1,iset2)
                      nz2  = npgfa(iset2)
                      nco2 = npgfa(iset2)*ncoset(la_max(iset2))
                      ns2  = nsgfa_set(iset2)
                      Int_ab_sum = 0.0_dp
                      CALL sumup_3c_integrals(&
                           lb_max(iset1),lb_min(iset1),nco1,m1,&
                           rpgfb(:,iset1),zetb(:,iset1),nz1,&
                           la_max(iset2),la_min(iset2),nco2,m2,&
                           rpgfa(:,iset2),zeta(:,iset2),nz2,&
                           l0_kkind,normhard(0:),ng0_s(0:),rpgf0_h,rpgf0_s,&
                           QN_h,-QN_s,Qlm,zhard,z0_s,&
                           -rab,rab2,rbc,rbc2,rac,rac2,&
                           Int_ab_sum,Int_abc,maxder,&
                           Int_adb_sum,Int_dab_sum)

                      CALL dgemm("N","N",nco1,ns2,nco2,&
                           1.0_dp,Int_ab_sum(1,1),maxcob,&
                           sphi_a(1,ifs2),SIZE(sphi_a,1),&
                           0.0_dp,work(1,1),SIZE(work,1))
                      CALL dgemm("T","N",ns1,ns2,nco1,&
                            1.0_dp,sphi_b(1,ifs1),SIZE(sphi_b,1),&
                           work(1,1),SIZE(work,1),&
                           1.0_dp,h_block(ifs1,ifs2),SIZE(h_block,1))

                      hab(1:nco1,1:nco2) = CPC_ab(1+m1:m1+nco1,1+m2:m2+nco2)  
                      trace = trace_r_AxB(hab,maxcob,Int_ab_sum,maxcob,nco1,nco2)
                      energy%hartree_3c = energy%hartree_3c + 2.0_dp*trace 

                      IF(calculate_forces) THEN
                        i_atom = atom_of_kind(iatom)
                        j_atom = atom_of_kind(jatom)
                        k_atom = atom_of_kind(katom)
                        DO idir = 1,3 
                          trace = trace_r_AxB(hab,maxcob,&
                                  Int_dab_sum(:,:,idir),maxcob,nco1,nco2)
                          force(ikind)%hartree_3c(idir,i_atom) = &
                                  force(ikind)%hartree_3c(idir,i_atom) - 2._dp*trace
                          force(kkind)%hartree_3c(idir,k_atom) = &
                                  force(kkind)%hartree_3c(idir,k_atom) + 2._dp*trace

                          trace = trace_r_AxB(hab,maxcob,&
                                  Int_adb_sum(:,:,idir),maxcob,nco1,nco2)
                          force(jkind)%hartree_3c(idir,j_atom) = &
                                  force(jkind)%hartree_3c(idir,j_atom) - 2._dp*trace
                          force(kkind)%hartree_3c(idir,k_atom) = &
                                  force(kkind)%hartree_3c(idir,k_atom) + 2._dp*trace
                        END DO 
                      END IF
                      m2 = m2 + maxcoa
                    END DO  ! iset2
                    m1 = m1 + maxcob
                  END DO  ! iset1

                END IF 

                ns1 = SIZE(CPC_ab,1)
                ns2 = SIZE(CPC_ab,2)

                DO ico = 1,nchan0_car
                  sum_ab_car(ico,katom) = sum_ab_car(ico,katom) + &
                                  ieqj*trace_r_AxB(CPC_ab,ns1,Int_abc(:,:,ico),&
                                  ns1,ns1,ns2)
                END DO  ! ico

                sbc_3c_neighbor_node => next(sbc_3c_neighbor_node)
              END DO knodeloop  ! knode
            END DO  !kkind
            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          END DO  ! inode
          sab_orb_neighbor_list => next(sab_orb_neighbor_list)
        END DO  ! ilist 

      END DO  ! jkind
    END DO  ! ikind

    CALL mp_sum(energy%hartree_3c,para_env%group)
    e3c = 0.0_dp

    IF(just_energy) GOTO 100
!   One center contributions from the 3-centers term
!   These contributions are added to the KS matrix together with the other
!   one center contributions, in the module qs_ks_atom
    DO kkind = 1,nkind

       NULLIFY(atom_kind,orb_basis,atom_list)
       atom_kind => atomic_kind_set(kkind)
       CALL get_atomic_kind(atomic_kind=atom_kind,&
               orb_basis_set=orb_basis,&
               atom_list=atom_list,natom=nat,paw_atom=paw_atom)

       CALL  get_gto_basis_set(gto_basis_set=orb_basis,&
                               nset =nsetc,maxso=maxsoc)
       nsototc = maxsoc * nsetc
       IF(.NOT. paw_atom) CYCLE

       NULLIFY(Qlm_gg)
       CALL get_rho0_mpole(rho0_mpole=rho0_mpole,ikind=kkind,&
                           l0_ikind=l0_kkind,Qlm_gg=Qlm_gg)

       nchan0_car = ncoset(l0_kkind)
       nchan0_sph = nsoset(l0_kkind)

       DO kat = 1,nat
         katom = atom_list(kat)

         NULLIFY(aVh3c0b)
         CALL get_ecoul_1c(ecoul_1c=ecoul_1c,iatom=katom,aVh3c0b=aVh3c0b)
         aVh3c0b = 0.0_dp

!        IF PARALLEL sum up the contributions to sum_ab_car coming from
!        ab pairs handled by different processors
         CALL mp_sum(sum_ab_car(:,katom),para_env%group)

         sum_ab_sph = 0.0_dp
         DO lshell = 0,l0_kkind
           DO is = 1,nso(lshell)
             iso = nsoset(lshell-1) + is 
             DO ic = 1,nco(lshell)
                ico = ncoset(lshell-1) + ic
                lx = indco(1,ico)
                ly = indco(2,ico)
                lz = indco(3,ico)
                sum_ab_sph(iso) =  sum_ab_sph(iso) + &
                   orbtramat(lshell)%c2s(is,ic)*sum_ab_car(ico,katom)/&
                   SQRT(fourpi*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)/&
                   dfac(2*lshell+1))
             END DO  ! ico
           END DO  ! iso
         END DO  ! lshell

         DO iso = 1,nchan0_sph
           DO is = 1,nsototc
              aVh3c0b(1:nsototc,is) = aVh3c0b(1:nsototc,is) + &
                            sum_ab_sph(iso)*Qlm_gg(1:nsototc,is,iso)
           END DO
         END DO 
       END DO  ! kat
    END DO  ! kkind

100 CONTINUE
    DEALLOCATE(sum_ab_car,sum_ab_sph,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,&
                       "sum_ab_car,sum_ab_sph")
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_a_vhg0_b

!******************************************************************************

  SUBROUTINE sumup_3c_integrals(l1_max,l1_min,nco1,m1,rp1,z1,nz1,&
                                l2_max,l2_min,nco2,m2,rp2,z2,nz2,&
                                l0_max,ng0_h,ng0_s,rp0_h,rp0_s,QN_h,QN_s,Qlm,z0_h,z0_s,&
                                r12,r122,r13,r132,r23,r232,vabc,Int_abc,maxder,&
                                vdabc,vadbc)

    INTEGER, INTENT(IN)                :: l1_max, l1_min, m1, nco1, nz1
    INTEGER, INTENT(IN)                :: l2_max, l2_min, m2, nco2, nz2
    INTEGER, INTENT(IN)                :: l0_max, maxder
    REAL(dp), INTENT(IN)               :: rp0_h, rp0_s, z0_h, z0_s
    REAL(dp), DIMENSION(:), INTENT(IN) :: rp1, z1
    REAL(dp), DIMENSION(:), INTENT(IN) :: rp2, z2
    REAL(dp), DIMENSION(0:),INTENT(IN) :: ng0_h, ng0_s
    REAL(dp), DIMENSION(:), INTENT(IN) :: Qlm
    REAL(dp), DIMENSION(:), INTENT(IN) :: QN_h, QN_s
    REAL(dp), INTENT(IN)               :: r122, r132, r232
    REAL(dp), DIMENSION(3), INTENT(IN) :: r12, r13, r23
    REAL(dp), DIMENSION(:, :), &
      INTENT(INOUT)                    :: vabc
    REAL(dp), DIMENSION(:,:,:), &
       POINTER                         :: Int_abc,vdabc,vadbc

    CHARACTER(LEN=*), PARAMETER   :: routine_name = "sumup_3c_integrals"

    INTEGER     ::  ico, ico1, ico1m, ico1p, ico2, ico2m, ico2p, &
                    idir, idx_1, idx_2, idx_1_L, idx_1m_L, idx_1p_L,&
                    idx_2_L, idx_2m_L, ipgf1, ipgf2, istat, l, n1,&
                    n1_plus, n1L, n2, n2_plus, n2L, nch0, npgf1, npgf2

    INTEGER, DIMENSION(3)              :: l1,l1m,l1p,l2,l2m
    REAL(dp)                           :: z1x2, z2x2

    REAL(dp), DIMENSION(:,:),&
                    ALLOCATABLE        :: vabc_plus
    REAL(dp), DIMENSION(:,:,:),&
                    ALLOCATABLE        :: int_tmp

!   ---------------------------------------------------------------------------

    nch0 = ncoset(l0_max)
    ALLOCATE(int_tmp(nco1,nco2,nch0),STAT=istat)
    int_tmp = 0.0_dp
    IF (istat /= 0) THEN
      CALL stop_memory(routine_name,module_name,__LINE__,"int_tmp",nco1*nco2*nch0*dp_size)
    END IF

    IF(maxder > 0) THEN
      IF (.NOT.ASSOCIATED(vdabc) .OR. .NOT.ASSOCIATED(vadbc)) THEN
        CALL stop_program(routine_name,&
             "The vadbc and vdabc pointers are not associated ")
      END IF

      n1_plus = nz1*ncoset(l1_max+1)
      n2_plus = nz2*ncoset(l2_max+1)
      ALLOCATE(vabc_plus(n1_plus,n2_plus),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"vabc_plus",&
                         n1_plus*n2_plus*dp_size)
      END IF

      vabc_plus = 0.0_dp

      CALL coulomb3(l1_max+maxder,nz1,z1,rp1,l1_min,&
                  l2_max+maxder,nz2,z2,rp2,l2_min,&
                  l0_max,z0_h,rp0_h,0,QN_h,&
                  r12,r122,r13,r132,r23,r232,vabc,int_tmp,&
                  maxder,vabc_plus)
    ELSE

      ALLOCATE(vabc_plus(1,1),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"vabc_plus",1*dp_size)
      END IF
      vabc_plus = 0.0_dp

      CALL coulomb3(l1_max,nz1,z1,rp1,l1_min,&
                  l2_max,nz2,z2,rp2,l2_min,&
                  l0_max,z0_h,rp0_h,0,QN_h,&
                  r12,r122,r13,r132,r23,r232,vabc,int_tmp,&
                  maxder,vabc_plus)
    END IF

    DO ico = 1,nch0
       l = indco(1,ico)+indco(2,ico)+indco(3,ico)

       Int_abc(1+m1:nco1+m1,1+m2:nco2+m2,ico) = &
                     ng0_h(l) * int_tmp(1:nco1,1:nco2,ico)
    END DO 
    
    IF(maxder > 0) THEN
      CALL coulomb3(l1_max+maxder,nz1,z1,rp1,l1_min,&
                  l2_max+maxder,nz2,z2,rp2,l2_min,&
                  l0_max,z0_s,rp0_s,0,QN_s,&
                  r12,r122,r13,r132,r23,r232,vabc,int_tmp,&
                  maxder,vabc_plus)
    ELSE
      CALL coulomb3(l1_max,nz1,z1,rp1,l1_min,&
                  l2_max,nz2,z2,rp2,l2_min,&
                  l0_max,z0_s,rp0_s,0,QN_s,&
                  r12,r122,r13,r132,r23,r232,vabc,int_tmp,&
                  maxder,vabc_plus)

    END IF

    DO ico = 1,nch0
       l = indco(1,ico)+indco(2,ico)+indco(3,ico)

       Int_abc(1+m1:nco1+m1,1+m2:nco2+m2,ico) = &
                    Int_abc(1+m1:nco1+m1,1+m2:nco2+m2,ico) - &
                    ng0_s(l) * int_tmp(1:nco1,1:nco2,ico)

    END DO

    IF(maxder > 0) THEN

      npgf1 = nz1
      npgf2 = nz2
      n1 = ncoset(l1_max)
      n1L = ncoset(l1_max+maxder)
      n2 = ncoset(l2_max)
      n2L = ncoset(l2_max+maxder)
      DO idir = 1,3
        DO ipgf1 = 1,npgf1
          z1x2 = z1(ipgf1)*2._dp
          DO ico1 = ncoset(l1_min-1)+1,ncoset(l1_max)
            l1(1:3) = indco(1:3,ico1)
            l1p(1:3) = l1(1:3)
            l1p(idir) =  l1(idir) + 1
            ico1p = coset(l1p(1),l1p(2),l1p(3))
            l1m(1:3) = l1(1:3)
            l1m(idir) =  l1(idir) - 1

            idx_1 = ico1 + (ipgf1-1)*n1
            idx_1_L = ico1 + (ipgf1-1)*n1L
            idx_1p_L = ico1p + (ipgf1-1)*n1L
 
            DO ipgf2 = 1,npgf2
              z2x2 = z2(ipgf1)*2._dp
              DO ico2 = ncoset(l2_min-1)+1,ncoset(l2_max)

                l2(1:3) = indco(1:3,ico2)
                l2m(1:3) = l2(1:3)
                l2m(idir) =  l2(idir) - 1
                ico2m = coset(l2m(1),l2m(2),l2m(3))

                idx_2 = ico2 + (ipgf2-1)*n2
                idx_2_L = ico2 + (ipgf2-1)*n2L

                IF(l1m(idir) < 0) THEN
                   vdabc(idx_1,idx_2,idir) = &
                   z1x2*vabc_plus(idx_1p_L,idx_2_L)
                ELSE
                  ico1m = coset(l1m(1),l1m(2),l1m(3))
                  idx_1m_L = ico1m + (ipgf1-1)*n1L
                  vdabc(idx_1,idx_2,idir) = &
                     z1x2*vabc_plus(idx_1p_L,idx_2_L)&
                    -l1(idir)*vabc_plus(idx_1m_L,idx_2_L)
                END IF

                IF(l2m(idir) < 0) THEN
                  vadbc(idx_1,idx_2,idir) = &
                     z2x2*(vabc_plus(idx_1p_L,idx_2_L)+&
                         r12(idir)*vabc(idx_1,idx_2))
                ELSE
                  ico2m = coset(l2m(1),l2m(2),l2m(3))
                  idx_2m_L = ico2m + (ipgf2-1)*n2L
                  vadbc(idx_1,idx_2,idir) = &
                     z2x2*(vabc_plus(idx_1p_L,idx_2_L)+&
                         r12(idir)*vabc(idx_1,idx_2))&
                    -l2(idir)*vabc_plus(idx_1_L,idx_2m_L)
                END IF

              END DO  ! ico2
            END DO  ! ipf2

          END DO  ! ico1
        END DO  ! ipgf1

      END DO  

    END IF

    DEALLOCATE(int_tmp, STAT=istat)
    DEALLOCATE(vabc_plus,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routine_name,module_name,__LINE__,"int_tmp,vabc_plus")
    END IF

  END SUBROUTINE sumup_3c_integrals

!******************************************************************************
END MODULE hartree_local_3centers
!******************************************************************************
