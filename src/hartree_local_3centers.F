!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/hartree_local_3centers
!!
!!   NAME
!!     hartree_local_3centers
!!
!!   FUNCTION
!!     Perform the calculation of the hartree 3 center terms for the GAPW method
!!
!!   AUTHOR
!!      MI
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
  USE input_constants,                 ONLY: tddfpt_singlet
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,dp_size,&
                                             int_size
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
                                             reduced_3c_list_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: trace_r_AxB, trace_r_AxB_new 
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit!, &
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: moduleN = "hartree_local_3centers"

! Public Subroutine
 
  PUBLIC :: integrate_a_vhg0_b

!******************************************************************************
 
  CONTAINS
 
!!****f* hartree_local_3centers/integrate_a_vhg0_b [1.0] *
!!
!!   NAME
!!     integrate_a_vhg0_b
!!
!!   FUNCTION
!!     Analytic calculation of the 3 center integrals for the GAPW 3 center terms
!!     This terms give the Hartree contribution due to the interaction of the soft density
!!     and the (rho0_hard-rho0_soft) density, i.e. the difference between the hard and
!!     soft compensation charges. (rho0_hard-rho0_soft) is a sum over atomic 
!!     contributions, which should be very shor ranged. 
!!     For this terms, properly constructed 3 center lists are used, in order to restrict
!!     the calculation only to those triplets of atoms with overlapping functions.
!!     For atom A and atom B the only soft primitives of the basis set are taken into 
!!     account, for the atom C, instead, the g0 primitives (see compensation charges)
!!     are used in the integrals.
!!
!!   NOTE
!!     For large systems the 3 center lists can become huge and this part of
!!     the calculation is the bottleneck. It needs optimization for better performance
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     MI
!!
!!
!!*** **********************************************************************

  SUBROUTINE integrate_a_vhg0_b(qs_env,h_mat,p_mat,energy_hartree_3c,just_energy,&
                                calculate_forces,tddft,p_env,error)

    REAL(kind=dp), INTENT(INOUT)               :: energy_hartree_3c
    TYPE(real_matrix_p_type), DIMENSION(:), &
      INTENT(IN)                               :: p_mat
    TYPE(real_matrix_p_type), DIMENSION(:), &
      INTENT(INOUT)                            :: h_mat
    TYPE(qs_environment_type), POINTER         :: qs_env
    LOGICAL, INTENT(IN)                        :: calculate_forces,&
                                                  just_energy
    LOGICAL, INTENT(IN), OPTIONAL              :: tddft
    TYPE(qs_p_env_type), POINTER, OPTIONAL     :: p_env
    TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "integrate_a_vhg0_b"
    REAL(dp), PARAMETER         :: EPS_a_vhg0_b = 1.E-8_dp ! other parameter for input ?

    TYPE(atomic_kind_type), DIMENSION(:),&
                                  POINTER      :: atomic_kind_set
    TYPE(atomic_kind_type),       POINTER      :: atom_kind
    TYPE(cell_type),              POINTER      :: cell
    TYPE(gto_basis_set_type),     POINTER      :: orb_basis
    TYPE(neighbor_list_set_p_type), &
                         DIMENSION(:), POINTER :: sab_orb,sbc_3c
    TYPE(neighbor_list_type),     POINTER      :: sab_orb_neighbor_list,&
                                                  sab_orb_neighbor_list_local
    TYPE(neighbor_node_type),     POINTER      :: sab_orb_neighbor_node
    TYPE(cp_para_env_type), POINTER            :: para_env
    TYPE(particle_type), DIMENSION(:), &
                                     POINTER   :: particle_set
    TYPE(qs_force_type), DIMENSION(:), POINTER :: force
    TYPE(reduced_3c_list_type), DIMENSION(:), &
                                  POINTER      :: reduced_3c_rho0
    TYPE(rho0_mpole_type),        POINTER      :: rho0_mpole 
    TYPE(rho_atom_type), DIMENSION(:), &
                                     POINTER   :: rho_atom_set
    TYPE(rho_atom_type),       POINTER         :: rho_atom
    TYPE(rho_atom_coeff), DIMENSION(:), &
        POINTER                                :: int_local_h,int_local_s
    REAL(dp) :: alpha_core, c1, dab, fac, ieqj, rab(3), rab2, rac(3), rac2, &
                rbc(3), rbc2, r_V0, trace, z0_h, z0_s, zhard , e3
    REAL(dp), DIMENSION(:),       POINTER      :: ng0_h, ng0_s
    REAL(dp), DIMENSION(:),   ALLOCATABLE      :: normhard, ff
    REAL(dp), DIMENSION(:),   ALLOCATABLE      :: sum_ab_sph
    REAL(dp), DIMENSION(:),       POINTER      :: QN_h, QN_s,Qlm
    REAL(dp), DIMENSION(:),   ALLOCATABLE      :: QN_h_c
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE    ::  work
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE  :: sum_ab_car
    REAL(dp), DIMENSION(:,:),     POINTER      :: h_block,p_block
    REAL(dp), DIMENSION(:,:),     POINTER      :: h_block2,p_block2
    REAL(dp), DIMENSION(:,:),     POINTER      :: rpgfa,rpgfb,&
                                                  sphi_a,sphi_b,&
                                                  zeta,zetb
    REAL(dp), DIMENSION(:,:),     POINTER      :: set_radius_3c
    REAL(dp), DIMENSION(:,:),     POINTER      :: hab
    REAL(dp), DIMENSION(:,:,:),    POINTER      :: CPC_ab
    REAL(dp), DIMENSION(:,:),     POINTER      :: Int_ab_sum
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE    :: int_tmpp
    REAL(dp), DIMENSION(:,:,:),   POINTER      :: Qlm_gg
    REAL(dp), DIMENSION(:,:,:),   POINTER      :: Int_abc
    REAL(dp), DIMENSION(:,:,:), POINTER        :: Int_dab_sum, Int_adb_sum
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE  :: vv

    INTEGER :: ab, bo(2), handle, iatom,  i_atom, i, ic, ico, icol, &
               idir, ifs1, ifs2, ikind, ilist, inode, irow, iset1, iset2,&
               is, iso, istat, ispin, jatom, j_atom, jkind, katom, kat, k_atom, &
               kkind, knode, l0_kkind, lmax_0, lshell, lx, ly, lz,&
               m1, m2, maxch0, maxco_global, maxcoa, &
               maxcob, maxder, maxlgto, maxlgto1, &
               maxsoc, maxsgf, maxsgf_set, mepos, nat, natom, nchan0_car,&
               nchan0_sph, nco1, nco2, ncotota, ncotota_old, &
               ncototb, ncototb_old, nkind, nlist, nnode, nnode_c,&
               ns1, ns2, nsgfa, nsgfb, nseta, nsetb, nsetc, &
               nsototc, nspins, num_pe, nz1, nz2
    INTEGER, DIMENSION(:),        POINTER      :: atom_list,la_max,la_min,&
                                                  lb_max,lb_min,&
                                                  npgfa,npgfb,nsgfa_set,nsgfb_set
    INTEGER, DIMENSION(:),        ALLOCATABLE  :: atom_of_kind
    INTEGER, DIMENSION(:,:),      POINTER      :: first_sgfa,first_sgfb
    LOGICAL                                    :: paw_atom, my_tddft
    INTEGER  :: ic1, ic1_min_m1, ic1_max, ic2, ic2_min_m1, ic2_max, &
                ip1,ip2, ldcpc1, ldcpc2, ldhab1, ldhab2, ldintab1, ldintab2, &
                n1, n2, nn1, nn2,jc
    REAL(KIND=dp) :: zeff

!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind_set,cell)
    NULLIFY( para_env,particle_set,rho0_mpole)
    NULLIFY(sab_orb,reduced_3c_rho0)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft

    ! TRIPLETS don't need a coulomb contribution
    IF (my_tddft) THEN
       IF ((qs_env%dft_control%nspins == 1) &
         .AND. (qs_env%dft_control%tddfpt_control%res_etype /= tddfpt_singlet)) RETURN
    END IF

    fac = 1.0_dp
    nspins = qs_env%dft_control%nspins 

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    para_env=para_env,& 
                    particle_set=particle_set,&
                    rho0_mpole=rho0_mpole,&
                    rho_atom_set=rho_atom_set,&
                    sab_orb=sab_orb)

    NULLIFY(ng0_h,ng0_s)
    CALL get_rho0_mpole(rho0_mpole=rho0_mpole,zet0_h=z0_h,zet0_s=z0_s,&
                        norm_g0l_h=ng0_h,norm_g0l_s=ng0_s,lmax_0=lmax_0,&
                        set_radius_3c=set_radius_3c)
    IF(z0_h==z0_s) RETURN

    IF (my_tddft) THEN
       rho0_mpole   => p_env%local_rho_set%rho0_mpole
       rho_atom_set => p_env%local_rho_set%rho_atom_set
    END IF

    nkind = SIZE(atomic_kind_set,1)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxco=maxco_global,maxlgto=maxlgto,&
                             maxsgf=maxsgf, maxsgf_set=maxsgf_set,&
                             natom=natom)

    maxder = 0
    IF (calculate_forces) THEN
 
      CALL timeset(routineN//" (forces)","I"," ",handle)

      ALLOCATE (atom_of_kind(natom),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "atom_of_kind",natom*int_size)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                               atom_of_kind=atom_of_kind)

      CALL get_qs_env(qs_env=qs_env, force=force)

      maxder = 1
      maxlgto1 = maxlgto + maxder
    ELSE
     
      CALL timeset(routineN,"I"," ",handle)

      maxlgto1 = maxlgto
    END IF



    maxch0 = ncoset(lmax_0)

!   Allocate  work arrays
    ALLOCATE(work(maxco_global,maxsgf_set),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "work",maxco_global*maxsgf_set)

    ALLOCATE(normhard(0:lmax_0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "normhard",(1+lmax_0)*dp_size)
    normhard(0:) = 0.0_dp
    ALLOCATE(QN_h_c(ncoset(lmax_0)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "QN_h_c",ncoset(lmax_0)*dp_size)
    QN_h_c(:) = 0.0_dp

    ALLOCATE(sum_ab_car(maxch0,natom),sum_ab_sph(maxch0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "sum_ab_car,sum_ab_sph",(natom*maxch0+maxch0)*dp_size)
    sum_ab_car = 0.0_dp

    NULLIFY(CPC_ab,hab)
    ALLOCATE(CPC_ab(1,1,nspins),STAT=istat)
    ldcpc1 = 1
    ldcpc2 = 1
    ALLOCATE(hab(1,1),STAT=istat)
    ldhab1 = 1
    ldhab2 = 1
    NULLIFY(Int_ab_sum,Int_dab_sum,Int_adb_sum,Int_abc)
    ALLOCATE(Int_ab_sum(1,1))
    ldintab1 = 1
    ldintab2 = 1
    CPC_ab = 0.0_dp
    hab    = 0.0_dp
    Int_ab_sum = 0.0_dp
    IF(calculate_forces) THEN
      ALLOCATE(Int_dab_sum(1,1,3),Int_adb_sum(1,1,3))
      Int_dab_sum=0.0_dp
      Int_adb_sum=0.0_dp
    END IF

    ALLOCATE(int_tmpp(maxco_global,maxco_global,maxch0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "int_tmpp",(maxco_global*maxco_global*maxch0)*dp_size)
    
    ALLOCATE(vv(ncoset(maxlgto1),ncoset(maxlgto),ncoset(lmax_0),&
         maxlgto1+maxlgto+lmax_0+1),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "vv",ncoset(maxlgto1)**2*ncoset(lmax_0)*(maxlgto1+maxlgto+lmax_0+1)*dp_size)
    ALLOCATE(ff(0:maxlgto1+maxlgto+lmax_0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "ff",(1+maxlgto1+maxlgto+lmax_0)*dp_size)

    ncotota_old = 0
    ncototb_old = 0

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
        
        DO ilist = 1,nlist

          sab_orb_neighbor_list_local => sab_orb_neighbor_list

          CALL get_neighbor_list(neighbor_list=sab_orb_neighbor_list_local,&
                                 atom=iatom, nnode=nnode)

          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab,&
                                   reduced_3c_rho0=reduced_3c_rho0)

            rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
            dab = SQRT(rab2)

            IF(iatom == jatom) THEN
              ieqj = 1.0_dp
            ELSE
              ieqj = 2.0_dp
            END IF 

            NULLIFY(h_block, p_block, p_block2, h_block2)
            IF (iatom <= jatom) THEN
!             Retrieve KS and density matrix for this couple of atoms
              CALL get_block_node(matrix=h_mat(1)%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=h_block)

              CALL get_block_node(matrix=p_mat(1)%matrix,&
                                  block_row=iatom,&
                                  block_col=jatom,&
                                  BLOCK=p_block)

              IF(nspins == 2) THEN
                CALL get_block_node(matrix=h_mat(2)%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    BLOCK=h_block2)

                CALL get_block_node(matrix=p_mat(2)%matrix,&
                                    block_row=iatom,&
                                    block_col=jatom,&
                                    BLOCK=p_block2)
              END IF

!             IF(ncotota .gt. ncotota_old .OR. ncototb .GT. ncototb_old) then
!                CALL reallocate(CPC_ab,1,ncotota,1,ncototb)
              CALL reallocate(CPC_ab,1,ncotota,1,ncototb,1,nspins)
              ldcpc1 = ncotota
              ldcpc2 = ncototb
!              END IF
              CPC_ab = 0.0_dp

!              IF(maxcoa .GT. SIZE(hab,1) .OR. maxcob .GT. SIZE(hab,2)) THEN
              IF(maxcoa .GT. ldhab1 .OR. maxcob .GT. ldhab2) THEN
                CALL reallocate(hab,1,maxcoa,1,maxcob)
                ldhab1 = maxcoa
                ldhab2 = maxcob
              ELSE
                hab = 0.0_dp
              END IF

!             Build the CPC coefficients that should then multiplied to the 
!             3-centers integrals of primitives to get the energy and the aVh3c0b
!             integrals for atom C (the one where the g0 are centered)

              m1 = 0
              DO iset1 = 1,nseta
                ifs1 = first_sgfa(1,iset1)
                nco1 = npgfa(iset1)*ncoset(la_max(iset1))
                ns1  = nsgfa_set(iset1)
                m2 = 0
                ic1_min_m1 = ncoset(la_min(iset1)-1)
                ic1_max    = ncoset(la_max(iset1))
                DO iset2 = 1,nsetb
                  ifs2 = first_sgfb(1,iset2)
                  nco2 = npgfb(iset2)*ncoset(lb_max(iset2))
                  ns2  = nsgfb_set(iset2)
                  ic2_min_m1 = ncoset(lb_min(iset2)-1)
                  ic2_max    = ncoset(lb_max(iset2))
 
                  CALL dgemm("N","N",nco1,ns2,ns1,&
                              1.0_dp,sphi_a(1,ifs1),maxcoa,&
                              p_block(ifs1,ifs2),nsgfa,&
                              0.0_dp,work(1,1),maxco_global)
                  CALL dgemm("N","T",nco1,nco2,ns2,&
                             1.0_dp,work(1,1),maxco_global,&
                             sphi_b(1,ifs2), maxcob,&
                             0.0_dp,hab(1,1),ldhab1)


                  DO ip1 = 1,npgfa(iset1)
                    ic1 = ic1_min_m1+1
                    n1 = ic1_max*(ip1-1)
                    nn1 = ic1_max*ip1
                    DO ip2 = 1,npgfb(iset2)
                      ic2 = ic2_min_m1+1
                      n2 =  ic2_max*(ip2-1)
                      nn2 = ic2_max*ip2
                      CPC_ab(ic1+n1+m1:nn1+m1,ic2+n2+m2:nn2+m2,1) = &
                                       hab(ic1+n1:nn1,ic2+n2:nn2)
                    END DO
                  END DO

                  IF(nspins == 2) THEN
                    CALL dgemm("N","N",nco1,ns2,ns1,&
                                1.0_dp,sphi_a(1,ifs1),maxcoa,&
                                p_block2(ifs1,ifs2),nsgfa,&
                                0.0_dp,work(1,1),maxco_global)
                    CALL dgemm("N","T",nco1,nco2,ns2,&
                               1.0_dp,work(1,1),maxco_global,&
                               sphi_b(1,ifs2), maxcob,&
                               0.0_dp,hab(1,1),ldhab1)


                    DO ip1 = 1,npgfa(iset1)
                      ic1 = ic1_min_m1+1
                      n1 = ic1_max*(ip1-1)
                      nn1 = ic1_max*ip1
                      DO ip2 = 1,npgfb(iset2)
                        ic2 = ic2_min_m1+1
                        n2 =  ic2_max*(ip2-1)
                        nn2 = ic2_max*ip2
                        CPC_ab(ic1+n1+m1:nn1+m1,ic2+n2+m2:nn2+m2,2) = &
                                         hab(ic1+n1:nn1,ic2+n2:nn2)
                      END DO
                    END DO

                  END IF

                 m2 = m2 + maxcob
                END DO  ! iset2
                m1 = m1 + maxcoa
              END DO  ! iset1

!             Allocate array for sum of integrals

              IF(calculate_forces) THEN
                 IF(maxcoa .GT. ldintab1 .OR. maxcob .GT. ldintab2) THEN
                    CALL reallocate(Int_ab_sum,1,maxcoa,1,maxcob)
                    CALL reallocate(Int_dab_sum,1,maxcoa,1,maxcob,1,3)
                    CALL reallocate(Int_adb_sum,1,maxcoa,1,maxcob,1,3)
                    ldintab1 = maxcoa
                    ldintab2 = maxcob
                 ELSE
                   Int_dab_sum=0.0_dp
                   Int_adb_sum=0.0_dp
                 END IF
              ELSE
                 IF(maxcoa .GT. ldintab1 .OR. maxcob .GT. ldintab2) THEN
                    CALL reallocate(Int_ab_sum,1,maxcoa,1,maxcob)
                    ldintab1 = maxcoa
                    ldintab2 = maxcob
                 END IF
              END IF
!             Allocate array for the 3-centers integrals
              CALL reallocate(Int_abc,1,ncotota,1,ncototb,1,maxch0)

            ELSE

              CALL get_block_node(matrix=h_mat(1)%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  BLOCK=h_block)

              CALL get_block_node(matrix=p_mat(1)%matrix,&
                                  block_row=jatom,&
                                  block_col=iatom,&
                                  BLOCK=p_block)

              IF(nspins == 2) THEN
                 CALL get_block_node(matrix=h_mat(2)%matrix,&
                                     block_row=jatom,&
                                     block_col=iatom,&
                                     BLOCK=h_block2)
 
                 CALL get_block_node(matrix=p_mat(2)%matrix,&
                                     block_row=jatom,&
                                     block_col=iatom,&
                                     BLOCK=p_block2)

              END IF

!              IF(ncototb .gt. ncototb_old .OR. ncotota .GT. ncotota_old) then
              CALL reallocate(CPC_ab,1,ncototb,1,ncotota,1,nspins)
              ldcpc1 = ncototb
              ldcpc2 = ncotota
!              END IF
               CPC_ab = 0.0_dp

              IF(maxcob .GT. ldhab1 .OR. maxcoa .GT. ldhab2) THEN
                CALL reallocate(hab,1,maxcob,1,maxcoa)
                ldhab1 = maxcob
                ldhab2 = maxcoa
              ELSE
                hab = 0.0_dp
              END IF

              m1 = 0
              DO iset1 = 1,nsetb
                ifs1 = first_sgfb(1,iset1)
                nco1 = npgfb(iset1)*ncoset(lb_max(iset1))
                ns1  = nsgfb_set(iset1)
                m2 = 0
                ic1_min_m1 = ncoset(lb_min(iset1)-1)
                ic1_max    = ncoset(lb_max(iset1))
                DO iset2 = 1,nseta
                  ifs2 = first_sgfa(1,iset2)
                  nco2 = npgfa(iset2)*ncoset(la_max(iset2))
                  ns2  = nsgfa_set(iset2)
                  ic2_min_m1 = ncoset(la_min(iset2)-1)
                  ic2_max    = ncoset(la_max(iset2))

                  CALL dgemm("N","N",nco1,ns2,ns1,&
                             1.0_dp,sphi_b(1,ifs1),maxcob,&
                             p_block(ifs1,ifs2),nsgfb,&
                             0.0_dp,work(1,1),maxco_global)
                  CALL dgemm("N","T",nco1,nco2,ns2,&
                             1.0_dp,work(1,1),maxco_global,&
                             sphi_a(1,ifs2), maxcoa,&
                             0.0_dp,hab(1,1),ldhab1)
!                             0.0_dp,hab(1,1),SIZE(hab,1))

                  DO ip1 = 1,npgfb(iset1)
                    ic1 = ic1_min_m1+1
                    n1 = ic1_max*(ip1-1)
                    nn1 = ic1_max*ip1
                    DO ip2 = 1,npgfa(iset2)
                      ic2 = ic2_min_m1+1
                      n2 =  ic2_max*(ip2-1)
                      nn2 = ic2_max*ip2
                      CPC_ab(ic1+n1+m1:nn1+m1,ic2+n2+m2:nn2+m2,1) = &
                                       hab(ic1+n1:nn1,ic2+n2:nn2)
                    END DO
                  END DO

                  IF (nspins==2) THEN

                    CALL dgemm("N","N",nco1,ns2,ns1,&
                               1.0_dp,sphi_b(1,ifs1),maxcob,&
                               p_block2(ifs1,ifs2),nsgfb,&
                               0.0_dp,work(1,1),maxco_global)
                    CALL dgemm("N","T",nco1,nco2,ns2,&
                               1.0_dp,work(1,1),maxco_global,&
                               sphi_a(1,ifs2), maxcoa,&
                               0.0_dp,hab(1,1),ldhab1)

                    DO ip1 = 1,npgfb(iset1)
                      ic1 = ic1_min_m1+1
                      n1 = ic1_max*(ip1-1)
                      nn1 = ic1_max*ip1
                      DO ip2 = 1,npgfa(iset2)
                        ic2 = ic2_min_m1+1
                        n2 =  ic2_max*(ip2-1)
                        nn2 = ic2_max*ip2
                        CPC_ab(ic1+n1+m1:nn1+m1,ic2+n2+m2:nn2+m2,2) = &
                                       hab(ic1+n1:nn1,ic2+n2:nn2)
                      END DO
                    END DO

                  END IF

                  m2 = m2 + maxcoa
                END DO  ! iset2
                m1 = m1 + maxcob
              END DO  ! iset1
!
!             Allocate array for sum of integrals

              IF(calculate_forces) THEN
                IF(maxcob .GT. ldintab1 .OR. maxcoa .GT. ldintab2) THEN
                  CALL reallocate(Int_ab_sum,1,maxcob,1,maxcoa)
                  CALL reallocate(Int_dab_sum,1,maxcob,1,maxcoa,1,3)
                  CALL reallocate(Int_adb_sum,1,maxcob,1,maxcoa,1,3)
                  ldintab1 = maxcob
                  ldintab2 = maxcoa
                ELSE
                  Int_dab_sum=0.0_dp
                  Int_adb_sum=0.0_dp
                END IF
              ELSE
                IF(maxcob .GT. ldintab1 .OR. maxcoa .GT. ldintab2) THEN
                  CALL reallocate(Int_ab_sum,1,maxcob,1,maxcoa)
                  ldintab1 = maxcob
                  ldintab2 = maxcoa
                END IF
              END IF
!             Allocate array for the 3-centers integrals
              CALL reallocate(Int_abc,1,ncototb,1,ncotota,1,maxch0)
            END IF

            DO kkind=1,nkind

              CALL get_rho0_mpole(rho0_mpole=rho0_mpole,ikind=kkind,&
                                  l0_ikind=l0_kkind, r_V0hmV0s=r_V0)

              nchan0_car = ncoset(l0_kkind)

              NULLIFY(atom_kind,orb_basis)
              atom_kind => atomic_kind_set(kkind)
              CALL get_atomic_kind(atomic_kind=atom_kind,&
                                   paw_atom=paw_atom,zeff=zeff,&
                                   alpha_core_charge=alpha_core)

              IF(paw_atom)THEN
                normhard(0:l0_kkind) = ng0_h(0:l0_kkind)
                zhard = z0_h
              ELSE
                c1 = alpha_core/pi
                c1 = c1*c1*c1
                normhard(:) = 0.0_dp
                normhard(0) = SQRT(c1) 
                zhard =  alpha_core
              END IF

              nnode_c = reduced_3c_rho0(kkind)%nnode
              DO knode=1,nnode_c

                katom = reduced_3c_rho0(kkind)%index_atom(knode)
                rac(1:3) = reduced_3c_rho0(kkind)%rac(1:3,knode)
                rac2 = reduced_3c_rho0(kkind)%rac2(knode)
                rbc(1:3) = reduced_3c_rho0(kkind)%rbc(1:3,knode)
                rbc2 = reduced_3c_rho0(kkind)%rbc2(knode)

                IF(ABS(rac(1)-(rab(1) + rbc(1)))   > EPS_a_vhg0_b .OR. &
                     ABS(rac(2)-(rab(2) + rbc(2))) > EPS_a_vhg0_b .OR. &
                     ABS(rac(3)-(rab(3) + rbc(3))) > EPS_a_vhg0_b ) THEN
                    CYCLE
                END IF
                
                Int_abc = 0.0_dp
                NULLIFY(QN_s,QN_h,Qlm)
                CALL  get_rho0_mpole(rho0_mpole=rho0_mpole, iat=katom,&
                                     Qlm_car_s=QN_s, Qlm_car_h=QN_h, Qlm_tot=Qlm)

                IF(paw_atom) THEN
                   QN_h_c(1:ncoset(l0_kkind)) = QN_h(1:ncoset(l0_kkind))
                ELSE
                  QN_h_c(:) = 0.0_dp
                  QN_h_c(1) = -normhard(0)*zeff
                ENDIF

                IF(iatom <= jatom) THEN

                  m1 = 0
                  DO iset1 = 1,nseta

!                    IF(rac2<=set_radius_3c(iset1,ikind)) THEN

                    ifs1 = first_sgfa(1,iset1)
                    nz1  = npgfa(iset1)
                    nco1 = npgfa(iset1)*ncoset(la_max(iset1))
                    ns1  = nsgfa_set(iset1)
                    m2 = 0
                    DO iset2 = 1,nsetb
!                      IF(rbc2<=set_radius_3c(iset2,jkind)) THEN

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
                           l0_kkind,normhard(0:),ng0_s(0:),r_V0,r_V0,&
                           QN_h_c,-QN_s,Qlm,zhard,z0_s,&
                           rab,rab2,rac,rac2,rbc,rbc2,&
                           Int_ab_sum,Int_abc,maxder,&
                           Int_dab_sum,Int_adb_sum,int_tmpp,vv,ff(0:))

                      CALL dgemm("N","N",nco1,ns2,nco2,&
                           1.0_dp,Int_ab_sum(1,1),ldintab1,&
                           sphi_b(1,ifs2),SIZE(sphi_b,1),&
                           0.0_dp,work(1,1),maxco_global)
                      
                      CALL dgemm("T","N",ns1,ns2,nco1,&
                            fac,sphi_a(1,ifs1),SIZE(sphi_a,1),&
                           work(1,1),maxco_global,&
                           1.0_dp,h_block(ifs1,ifs2),SIZE(h_block,1))
!
                      IF(nspins == 2) THEN
                        CALL dgemm("T","N",ns1,ns2,nco1,&
                              fac,sphi_a(1,ifs1),SIZE(sphi_a,1),&
                             work(1,1),maxco_global,&
                             1.0_dp,h_block2(ifs1,ifs2),SIZE(h_block,1))
                      END IF

                      DO ispin = 1,nspins
                        trace = trace_r_AxB_new(CPC_ab(:,:,ispin),ldcpc1,m1,m2,&
                                          Int_ab_sum,ldintab1,0,0,nco1,nco2)
                        energy_hartree_3c = energy_hartree_3c + ieqj*trace 
                      ENDDO

                      IF(calculate_forces) THEN

                        i_atom = atom_of_kind(iatom)
                        j_atom = atom_of_kind(jatom)
                        k_atom = atom_of_kind(katom)
                        DO idir = 1,3 
                          DO ispin = 1,nspins
                            trace = trace_r_AxB_new(CPC_ab(:,:,ispin),ldcpc1,m1,m2,&
                                    Int_dab_sum(:,:,idir),ldintab1,0,0,nco1,nco2)

                            force(ikind)%hartree_3c(idir,i_atom) = &
                                    force(ikind)%hartree_3c(idir,i_atom) - ieqj*trace
                            force(kkind)%hartree_3c(idir,k_atom) = &
                                    force(kkind)%hartree_3c(idir,k_atom) + ieqj*trace

                            trace = trace_r_AxB_new(CPC_ab(:,:,ispin),ldcpc1,m1,m2,&
                                    Int_adb_sum(:,:,idir),ldintab1,0,0,nco1,nco2)

                            force(jkind)%hartree_3c(idir,j_atom) = &
                                    force(jkind)%hartree_3c(idir,j_atom) - ieqj*trace
                            force(kkind)%hartree_3c(idir,k_atom) = &
                                    force(kkind)%hartree_3c(idir,k_atom) + ieqj*trace
                          END DO
                        END DO 

                      END IF

!                      END IF  ! set_radius_3c

                      m2 = m2 + maxcob
                    END DO  ! iset2

!                    END IF  ! set_radius_3c

                    m1 = m1 + maxcoa
                  END DO  ! iset1

                ELSE

                  m1 = 0
                  DO iset1 = 1,nsetb

!                    IF(rbc2<=set_radius_3c(iset1,jkind)) THEN

                    ifs1 = first_sgfb(1,iset1)
                    nz1  = npgfb(iset1)
                    nco1 = npgfb(iset1)*ncoset(lb_max(iset1))
                    ns1  = nsgfb_set(iset1)
                    m2 = 0
                    DO iset2 = 1,nseta

!                      IF(rac2<=set_radius_3c(iset2,ikind)) THEN

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
                           l0_kkind,normhard(0:),ng0_s(0:),r_V0,r_V0,&
                           QN_h_c,-QN_s,Qlm,zhard,z0_s,&
                           -rab,rab2,rbc,rbc2,rac,rac2,&
                           Int_ab_sum,Int_abc,maxder,&
                           Int_adb_sum,Int_dab_sum,int_tmpp,vv,ff(0:))

                      CALL dgemm("N","N",nco1,ns2,nco2,&
                           1.0_dp,Int_ab_sum(1,1),ldintab1,&
                           sphi_a(1,ifs2),SIZE(sphi_a,1),&
                           0.0_dp,work(1,1),maxco_global)
                      CALL dgemm("T","N",ns1,ns2,nco1,&
                            fac,sphi_b(1,ifs1),SIZE(sphi_b,1),&
                           work(1,1),maxco_global,&
                           1.0_dp,h_block(ifs1,ifs2),SIZE(h_block,1))

                      IF(nspins == 2) THEN
                        CALL dgemm("T","N",ns1,ns2,nco1,&
                              fac,sphi_b(1,ifs1),SIZE(sphi_b,1),&
                             work(1,1),maxco_global,&
                             1.0_dp,h_block2(ifs1,ifs2),SIZE(h_block,1))
                      END IF

 
                      DO ispin = 1,nspins
                        trace = trace_r_AxB_new(CPC_ab(:,:,ispin),ldcpc1,m1,m2,&
                                            Int_ab_sum,ldintab1,0,0,nco1,nco2)

                        energy_hartree_3c = energy_hartree_3c + 2.0_dp*trace 
                      END DO

                      IF(calculate_forces) THEN

                        i_atom = atom_of_kind(iatom)
                        j_atom = atom_of_kind(jatom)
                        k_atom = atom_of_kind(katom)
                        DO idir = 1,3 
                          DO ispin = 1,nspins
                            trace = trace_r_AxB_new(CPC_ab(:,:,ispin),ldcpc1,m1,m2,&
                                    Int_dab_sum(:,:,idir),ldintab1,0,0,nco1,nco2)

                            force(ikind)%hartree_3c(idir,i_atom) = &
                                    force(ikind)%hartree_3c(idir,i_atom) - 2._dp*trace
                            force(kkind)%hartree_3c(idir,k_atom) = &
                                    force(kkind)%hartree_3c(idir,k_atom) + 2._dp*trace

                            trace = trace_r_AxB_new(CPC_ab(:,:,ispin),ldcpc1,m1,m2,&
                                    Int_adb_sum(:,:,idir),ldintab1,0,0,nco1,nco2)
  
                            force(jkind)%hartree_3c(idir,j_atom) = &
                                    force(jkind)%hartree_3c(idir,j_atom) - 2._dp*trace
                            force(kkind)%hartree_3c(idir,k_atom) = &
                                    force(kkind)%hartree_3c(idir,k_atom) + 2._dp*trace
                          END DO
                        END DO 
                      END IF

!                      END IF  ! set_radius_3c

                      m2 = m2 + maxcoa
                    END DO  ! iset2

!                    END IF  ! set_radius_3c

                    m1 = m1 + maxcob
                  END DO  ! iset1

                END IF 

                DO ispin = 1,nspins
                  DO ico = 1,nchan0_car
                    sum_ab_car(ico,katom) = sum_ab_car(ico,katom) + &
                                 ieqj*trace_r_AxB(CPC_ab(:,:,ispin),&
                                            ldcpc1,Int_abc(:,:,ico),&
                                            ldcpc1,ldcpc1,ldcpc2)
                  END DO  ! ico
                END DO

              END DO  ! knode
            END DO  !kkind

            sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          END DO  ! inode
          sab_orb_neighbor_list => next(sab_orb_neighbor_list)
        END DO  ! ilist 

      END DO  ! jkind
    END DO  ! ikind

!   IF PARALLEL sum up the contributions to sum_ab_car coming from
!   ab pairs handled by different processors
    CALL mp_sum(energy_hartree_3c,para_env%group)
    CALL mp_sum(sum_ab_car,para_env%group)

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

!      Distribute the atoms of this kind
       num_pe = para_env%num_pe
       mepos  = para_env%mepos
       bo = get_limit( nat, num_pe, mepos )

       DO kat = bo(1), bo(2) !1,nat
         katom = atom_list(kat)

         rho_atom => rho_atom_set(katom)
         CALL get_rho_atom(rho_atom=rho_atom,&
              ga_Vlocal_gb_h=int_local_h,ga_Vlocal_gb_s=int_local_s)

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

             DO ispin = 1,nspins

                int_local_h(ispin)%r_coef(1:nsototc,is) = &
                          int_local_h(ispin)%r_coef(1:nsototc,is) + &
                          sum_ab_sph(iso)*Qlm_gg(1:nsototc,is,iso)
                int_local_s(ispin)%r_coef(1:nsototc,is) = &
                          int_local_s(ispin)%r_coef(1:nsototc,is) + &
                          sum_ab_sph(iso)*Qlm_gg(1:nsototc,is,iso)
             END DO
           END DO
         END DO 


       END DO  ! kat

    END DO  ! kkind

100 CONTINUE

    IF(ASSOCIATED(Int_ab_sum)) THEN
       DEALLOCATE(Int_ab_sum,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "Int_ab_sum")
    END IF
    IF(ASSOCIATED(Int_abc)) THEN
       DEALLOCATE(Int_abc,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "Int_abc")

    END IF
    IF(calculate_forces) THEN
       DEALLOCATE (atom_of_kind,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
            "atom_of_kind")
       IF(ASSOCIATED(Int_dab_sum)) THEN
         DEALLOCATE(Int_dab_sum,Int_adb_sum,STAT=istat)
         IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                      "Int_dab_sum,Int_adb_sum")
       END IF
    END IF
    IF(ASSOCIATED(CPC_ab)) THEN
      DEALLOCATE(CPC_ab, STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
          "CPC_ab") 
    END IF
    IF(ASSOCIATED(hab)) THEN
      DEALLOCATE(hab, STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
          "hab") 
    END IF

    DEALLOCATE(int_tmpp,vv,ff,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
           "int_tmpp,vv,ff")
    END IF

    DEALLOCATE(work,sum_ab_car,sum_ab_sph,STAT=istat)
    IF (istat /= 0) THEN
      CALL stop_memory(routineN,moduleN,__LINE__,&
           "work,sum_ab_car,sum_ab_sph")
    END IF
    DEALLOCATE(normhard,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "normhard")
    DEALLOCATE(QN_h_c,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "QN_h_c")

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE integrate_a_vhg0_b

!******************************************************************************
!!****f* hartree_local_3centers/sumup_3c_integrals [1.0] *
!!
!!   NAME
!!     sumup_3c_integrals
!!
!!   FUNCTION
!!     Driver for the calculation of the 3 center integrals.
!!     It already performes the difference (rho0_h_C-rho0_s_C)rho_soft
!!
!!   NOTE
!!     For large systems the 3 center lists can become huge and this part of
!!     the calculation is the bottleneck. This routine is called too many times
!!     It needs optimization for better performance
!!
!!   AUTHOR
!!     MI
!!
!!
!!*** **********************************************************************


  SUBROUTINE sumup_3c_integrals(l1_max,l1_min,nco1,m1,rp1,z1,nz1,&
                                l2_max,l2_min,nco2,m2,rp2,z2,nz2,&
                                l0_max,ng0_h,ng0_s,rp0_h,rp0_s,QN_h,QN_s,Qlm,z0_h,z0_s,&
                                r12,r122,r13,r132,r23,r232,vabc,Int_abc,maxder,&
                                vdabc,vadbc,int_tmp,vv,ff)

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
    REAL(dp), DIMENSION(0:)            :: ff
    REAL(dp), DIMENSION(:, :), &
      INTENT(INOUT)                    :: vabc
    REAL(dp), DIMENSION(:,:,:)         :: int_tmp
    REAL(dp), DIMENSION(:,:,:), &
       POINTER                         :: Int_abc,vdabc,vadbc
    REAL(dp), DIMENSION(:,:,:,:)       :: vv

    CHARACTER(LEN=*), PARAMETER   :: routineN = "sumup_3c_integrals"

    INTEGER     ::  handle,ico, ico1, ico1m, ico1p, ico2, ico2m, &
                    idir, idx_1, idx_2, idx_1_L, idx_1m_L, idx_1p_L,&
                    idx_2m, ipgf1, ipgf2, istat, l, n1,&
                    n1_plus, n1L, n2, n2_plus, nch0, npgf1, npgf2

    INTEGER, DIMENSION(3)              :: l1,l1m,l1p,l2,l2m
    REAL(dp)                           :: z1x2, z2x2

    REAL(dp), DIMENSION(:,:),&
                    ALLOCATABLE        :: vabc_plus

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I"," ",handle)
    nch0 = ncoset(l0_max)

    IF(maxder > 0) THEN
      IF (.NOT.ASSOCIATED(vdabc) .OR. .NOT.ASSOCIATED(vadbc)) THEN
        CALL stop_program(routineN,&
             "The vadbc and vdabc pointers are not associated ")
      END IF

      n1_plus = nz1*ncoset(l1_max+1)
      n2_plus = nz2*ncoset(l2_max)
      ALLOCATE(vabc_plus(n1_plus,n2_plus),STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"vabc_plus",&
                         n1_plus*n2_plus*dp_size)
      END IF

      vabc_plus = 0.0_dp

      CALL coulomb3(l1_max+maxder,nz1,z1,rp1,l1_min,&
                  l2_max,nz2,z2,rp2,l2_min,&
                  l0_max,z0_h,rp0_h,0,QN_h,&
                  r12,r122,r13,r132,r23,r232,vabc,int_tmp,&
                  vv,ff(0:),maxder=maxder,vabc_plus=vabc_plus)
    ELSE


      CALL coulomb3(l1_max,nz1,z1,rp1,l1_min,&
                  l2_max,nz2,z2,rp2,l2_min,&
                  l0_max,z0_h,rp0_h,0,QN_h,&
                  r12,r122,r13,r132,r23,r232,vabc,int_tmp,&
                  vv,ff(0:))
    END IF

    DO ico = 1,nch0
       l = indco(1,ico)+indco(2,ico)+indco(3,ico)

!       Int_abc(1+m1:nco1+m1,1+m2:nco2+m2,ico) = &
!                     ng0_h(l) * int_tmp(1:nco1,1:nco2,ico)
       DO ico1 = 1,nco2
         CALL DSCAL(nco1,ng0_h(l),int_tmp(1,ico1,ico),1)
         CALL DCOPY(nco1,int_tmp(1,ico1,ico),1,Int_abc(m1+1,m2+ico1,ico),1)
       END DO
    END DO 
    
    IF(maxder > 0) THEN
      CALL coulomb3(l1_max+maxder,nz1,z1,rp1,l1_min,&
                  l2_max,nz2,z2,rp2,l2_min,&
                  l0_max,z0_s,rp0_s,0,QN_s,&
                  r12,r122,r13,r132,r23,r232,vabc,int_tmp,&
                  vv,ff(0:),maxder=maxder,vabc_plus=vabc_plus)
    ELSE
      CALL coulomb3(l1_max,nz1,z1,rp1,l1_min,&
                  l2_max,nz2,z2,rp2,l2_min,&
                  l0_max,z0_s,rp0_s,0,QN_s,&
                  r12,r122,r13,r132,r23,r232,vabc,int_tmp,&
                  vv,ff(0:))

    END IF

    DO ico = 1,nch0
       l = indco(1,ico)+indco(2,ico)+indco(3,ico)

!       Int_abc(1+m1:nco1+m1,1+m2:nco2+m2,ico) = &
!                    Int_abc(1+m1:nco1+m1,1+m2:nco2+m2,ico) - &
!                    ng0_s(l) * int_tmp(1:nco1,1:nco2,ico)

       DO ico1 = 1,nco2
         CALL DAXPY(nco1,-ng0_s(l),int_tmp(1,ico1,ico),1,Int_abc(m1+1,m2+ico1,ico),1)
       END DO
    END DO

    IF(maxder > 0) THEN

      npgf1 = nz1
      npgf2 = nz2
      n1 = ncoset(l1_max)
      n1L = ncoset(l1_max+maxder)
      n2 = ncoset(l2_max)
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
              z2x2 = z2(ipgf2)*2._dp
              DO ico2 = ncoset(l2_min-1)+1,ncoset(l2_max)

                l2(1:3) = indco(1:3,ico2)
                l2m(1:3) = l2(1:3)
                l2m(idir) =  l2(idir) - 1
                ico2m = coset(l2m(1),l2m(2),l2m(3))

                idx_2 = ico2 + (ipgf2-1)*n2

                IF(l1m(idir) < 0) THEN
                   vdabc(idx_1,idx_2,idir) = &
                   -z1x2*vabc_plus(idx_1p_L,idx_2)
                ELSE
                  ico1m = coset(l1m(1),l1m(2),l1m(3))
                  idx_1m_L = ico1m + (ipgf1-1)*n1L
                  vdabc(idx_1,idx_2,idir) = &
                     -z1x2*vabc_plus(idx_1p_L,idx_2)&
                    +REAL(l1(idir),dp)*vabc_plus(idx_1m_L,idx_2)
                END IF

                IF(l2m(idir) < 0) THEN
                  vadbc(idx_1,idx_2,idir) = &
                    -z2x2*(vabc_plus(idx_1p_L,idx_2)-&
                         r12(idir)*vabc(idx_1,idx_2))
                ELSE
                  ico2m = coset(l2m(1),l2m(2),l2m(3))
                  idx_2m = ico2m + (ipgf2-1)*n2
                  vadbc(idx_1,idx_2,idir) = &
                     -z2x2*(vabc_plus(idx_1p_L,idx_2)-&
                         r12(idir)*vabc(idx_1,idx_2))&
                    +REAL(l2(idir),dp)*vabc_plus(idx_1_L,idx_2m)
                END IF

              END DO  ! ico2
            END DO  ! ipf2

          END DO  ! ico1
        END DO  ! ipgf1

      END DO  

    END IF

    IF(maxder > 0) THEN
      DEALLOCATE(vabc_plus,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routineN,moduleN,__LINE__,"vabc_plus")
      END IF
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE sumup_3c_integrals

!******************************************************************************
END MODULE hartree_local_3centers
!******************************************************************************
