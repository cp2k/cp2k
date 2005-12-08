!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/hartree_local_2centers
!!
!!   NAME
!!     hartree_local_2centers
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   SOURCE
!******************************************************************************

MODULE  hartree_local_2centers

  USE ai_coulomb,                      ONLY: coulomb2
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
                                             next
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                             rho0_mpole_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: trace_r_AxB
  USE simulation_cell,                 ONLY: cell_type,&
                                             pbc
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop 
  USE util,                            ONLY: get_limit

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(LEN=*), PARAMETER :: module_name = "hartree_local_3centers"

! Public Subroutine
 
  PUBLIC :: integrate_vhg0_g0

!******************************************************************************
 
  CONTAINS
 
!******************************************************************************

    SUBROUTINE integrate_vhg0_g0(qs_env,energy_hartree_2c,just_energy,calculate_forces, &
                               tddft, p_env)

    REAL(kind=dp), INTENT(INOUT)               :: energy_hartree_2c
    TYPE(qs_environment_type), POINTER         :: qs_env
    LOGICAL, INTENT(IN)                        :: calculate_forces,&
                                                  just_energy
    LOGICAL, INTENT(IN), OPTIONAL              :: tddft
    TYPE(qs_p_env_type), POINTER, OPTIONAL     :: p_env

    CHARACTER(LEN=*), PARAMETER :: routine_name = "integrate_vhg0_g0"

    TYPE(atomic_kind_type), DIMENSION(:),&
                                  POINTER      :: atomic_kind_set
    TYPE(atomic_kind_type),       POINTER      :: atom_kind
    TYPE(cell_type),              POINTER      :: cell
    TYPE(gto_basis_set_type),     POINTER      :: orb_basis
    TYPE(neighbor_list_set_p_type), &
                         DIMENSION(:), POINTER :: sab_2c
    TYPE(neighbor_list_type),     POINTER      :: sab_2c_neighbor_list,&
                                                  sab_2c_neighbor_list_local
    TYPE(neighbor_node_type),     POINTER      :: sab_2c_neighbor_node
    TYPE(cp_para_env_type), POINTER            :: para_env
    TYPE(particle_type), DIMENSION(:), &
                                     POINTER   :: particle_set
    TYPE(qs_force_type), DIMENSION(:), POINTER :: force
    TYPE(rho0_mpole_type),        POINTER      :: rho0_mpole
    TYPE(rho_atom_type), DIMENSION(:), &
                                     POINTER   :: rho_atom_set
    TYPE(rho_atom_type),       POINTER         :: rho_atom
    TYPE(rho_atom_coeff), DIMENSION(:), &
        POINTER                                :: int_local_h,int_local_s

    REAL(dp) :: alpha_core, c1, ene_2c, force_a(3), ieqj, ra(3), rab(3), rab2, rb(3),&
                r0cut_a, r0cut_b,  rpgf0_s_a, rpgf0_s_b, z0_h, z0_h_a, z0_h_b, z0_s
    REAL(dp), DIMENSION(:),       POINTER      :: ng0_h, ng0_s
    REAL(dp), DIMENSION(:),       ALLOCATABLE  :: nhard_a, nhard_b
    REAL(dp), DIMENSION(:),       POINTER      :: QN_h,QN_s
    REAL(dp), DIMENSION(:),       ALLOCATABLE  :: QN_h_a,QN_s_a,&
                                                  QN_h_b,QN_s_b 
    REAL(dp), DIMENSION(:),       ALLOCATABLE  :: Isum_a_sph
    REAL(dp), DIMENSION(:),       ALLOCATABLE  :: ff
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE  :: Isum_a_car
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE  :: vab
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE  :: vv
    REAL(dp), DIMENSION(:,:,:),   POINTER      :: Qlm_gg

    INTEGER :: ab, bo(2), handle, i_atom, iat, iatom, ic, ico, idir, ikind, ilist, &
               inode, is, iso, ispin, istat, j_atom, jatom, jkind, l0_a, l0_b, &
               lmax_0, lmax_01, lshell, lx, ly, lz, maxch0, maxder, maxso, &
               mepos, nat, natom, nch0_sph, nkind, nlist, nnode, nset, nsotot, num_pe
    REAL(KIND=dp) ::  zeff_a, zeff_b
    INTEGER, DIMENSION(:),        POINTER      :: atom_list
    INTEGER, DIMENSION(:),        ALLOCATABLE  :: atom_of_kind
    LOGICAL                                    :: paw_atom_a, paw_atom_b, my_tddft

!   ---------------------------------------------------------------------------

    NULLIFY(atomic_kind_set,cell)
    NULLIFY(para_env,particle_set,rho0_mpole)
    NULLIFY(rho_atom_set,sab_2c)

    my_tddft = .FALSE.
    IF (PRESENT(tddft)) my_tddft = tddft
    IF (my_tddft) THEN
       IF (qs_env%dft_control%nspins==1 .AND. &
            qs_env%dft_control%tddfpt_control%res_etype/=tddfpt_singlet) RETURN
    END IF

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                    cell=cell,&
                    para_env=para_env,& 
                    particle_set=particle_set,&
                    rho0_mpole=rho0_mpole,&
                    rho_atom_set=rho_atom_set,&
                    sab_2c=sab_2c)

    NULLIFY(ng0_h,ng0_s)
    CALL get_rho0_mpole(rho0_mpole=rho0_mpole,zet0_h=z0_h,zet0_s=z0_s,&
                        norm_g0l_h=ng0_h,norm_g0l_s=ng0_s,lmax_0=lmax_0)

    IF(z0_h==z0_s) THEN
      RETURN
    END IF

    IF (my_tddft) THEN
       rho0_mpole   => p_env%local_rho_set%rho0_mpole
       rho_atom_set => p_env%local_rho_set%rho_atom_set
    END IF

    nkind = SIZE(atomic_kind_set,1)
    natom = SIZE(particle_set,1)

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


    IF (calculate_forces) THEN
      lmax_01 = lmax_0 + maxder
    ELSE
      lmax_01 = lmax_0
    END IF

    maxch0 = ncoset(lmax_0)
    ALLOCATE(nhard_a(0:lmax_0),nhard_b(0:lmax_0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "nhard_a,nhard_b",(2+2*lmax_0)*dp_size)
    nhard_a(0:) = 0.0_dp
    nhard_b(0:) = 0.0_dp
    ALLOCATE(QN_h_a(maxch0),QN_h_b(maxch0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "QN_h_a,QN_h_b",(2*maxch0)*dp_size)
    QN_h_a(:) = 0.0_dp
    QN_h_b(:) = 0.0_dp
    ALLOCATE(QN_s_a(maxch0),QN_s_b(maxch0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "QN_s_a,QN_s_b",(2*maxch0)*dp_size)
    QN_s_a(:) = 0.0_dp
    QN_s_b(:) = 0.0_dp
    ALLOCATE(Isum_a_car(maxch0,natom),Isum_a_sph(maxch0),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "Isum_a_car, Isum_a_sph",(natom*maxch0+maxch0)*dp_size)
    Isum_a_car = 0.0_dp

    ALLOCATE(vab(maxch0,maxch0),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "vab",maxch0*maxch0*dp_size)
    ALLOCATE(vv(ncoset(lmax_01),maxch0,lmax_01+lmax_0+1),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "vv",maxch0*maxch0*(lmax_01+lmax_0+1)*dp_size)
    ALLOCATE(ff(0:lmax_01+lmax_0),STAT = istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "ff", (lmax_0+lmax_0+1)*dp_size)

    DO ikind = 1,nkind

      NULLIFY(atom_kind)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,&
                           paw_atom=paw_atom_a,zeff=zeff_a,&
                           alpha_core_charge=alpha_core)

      CALL get_rho0_mpole(rho0_mpole=rho0_mpole,ikind=ikind,&
                              l0_ikind=l0_a,r_V0hmV0s=r0cut_a,&
                              rpgf0_s=rpgf0_s_a)
      nhard_a = 0.0_dp
      IF(paw_atom_a) THEN
        nhard_a(0:l0_a) = ng0_h(0:l0_a)
        z0_h_a = z0_h
      ELSE
        c1 = alpha_core/pi
        c1 = c1*c1*c1
        nhard_a(0) = SQRT(c1)
        z0_h_a = alpha_core
      END IF

      DO jkind = 1,nkind

        ab = ikind + nkind*(jkind-1)

        NULLIFY(sab_2c_neighbor_list)
        IF(ASSOCIATED(sab_2c(ab)%neighbor_list_set)) THEN
          CALL get_neighbor_list_set(neighbor_list_set=&
               sab_2c(ab)%neighbor_list_set,nlist=nlist)
 
          sab_2c_neighbor_list => first_list(sab_2c(ab)%neighbor_list_set)
        ELSE
          nlist = 0
        END IF

        NULLIFY(atom_kind)
        atom_kind => atomic_kind_set(jkind)
        CALL get_atomic_kind(atomic_kind=atom_kind,&
                            paw_atom=paw_atom_b,zeff=zeff_b,&
                            alpha_core_charge=alpha_core)

        CALL get_rho0_mpole(rho0_mpole=rho0_mpole,ikind=jkind,&
                            l0_ikind=l0_b,r_V0hmV0s=r0cut_b,&
                            rpgf0_s=rpgf0_s_b)

        nhard_b = 0.0_dp
        IF(paw_atom_b) THEN
          nhard_b(0:l0_b) = ng0_h(0:l0_b)
          z0_h_b = z0_h
        ELSE
          c1 = alpha_core/pi
          c1 = c1*c1*c1
          nhard_b(0) = SQRT(c1)
          z0_h_b = alpha_core
        END IF

        DO ilist = 1,nlist

          sab_2c_neighbor_list_local => sab_2c_neighbor_list
          CALL get_neighbor_list(neighbor_list=sab_2c_neighbor_list_local,&
                                 atom=iatom,nnode=nnode)

          ra(:) = pbc(particle_set(iatom)%r,cell)

          sab_2c_neighbor_node => first_node(sab_2c_neighbor_list_local)

          NULLIFY(QN_h,QN_s)
          CALL  get_rho0_mpole(rho0_mpole=rho0_mpole, iat=iatom,&
                               Qlm_car_h=QN_h,Qlm_car_s=QN_s)

          IF(paw_atom_a) THEN
             QN_h_a(1:ncoset(l0_a)) = QN_h(1:ncoset(l0_a))
          ELSE
             QN_h_a(1) = -nhard_a(0)*zeff_a
          ENDIF
          QN_s_a(1:ncoset(l0_a)) = QN_s(1:ncoset(l0_a))

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_2c_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab)
            rab2  = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

            rb(:) = pbc(particle_set(jatom)%r,cell)

            IF(iatom == jatom) THEN
              ieqj = 1.0_dp
            ELSE
              ieqj = 2.0_dp
            END IF 

            NULLIFY(QN_s,QN_h)
            CALL  get_rho0_mpole(rho0_mpole=rho0_mpole, iat=jatom,&
                                 Qlm_car_s=QN_s, Qlm_car_h=QN_h)

            IF(paw_atom_b) THEN
               QN_h_b(1:ncoset(l0_b)) = QN_h(1:ncoset(l0_b))
            ELSE
               QN_h_b(1) = -nhard_b(0)*zeff_b
            ENDIF
            QN_s_b(1:ncoset(l0_b)) = QN_s(1:ncoset(l0_b))

            IF(calculate_forces) THEN
              CALL sumup_2c_integrals(z0_h_a,z0_h_b,z0_s,&
                                      nhard_a(0:), nhard_b(0:),ng0_s(0:),&
                                      l0_a,r0cut_a,QN_h_a,QN_s_a,&
                                      l0_b,r0cut_b,QN_h_b,QN_s_b,&
!                       l0_a,rpgf0_s_a,QN_h_a,QN_s_a,&
!                       l0_b,rpgf0_s_b,QN_h_b,QN_s_b,&
                                      rab,rab2,ieqj,ene_2c,Isum_a_car(:,iatom),&
                                      Isum_a_car(:,jatom),maxder,force_a,&
                                      vab,vv,ff(0:))

!             ***  sum up in the force arrey  ***
              i_atom = atom_of_kind(iatom)
              j_atom = atom_of_kind(jatom)
              DO idir = 1,3
                force(ikind)%hartree_2c(idir,i_atom) = &
                  force(ikind)%hartree_2c(idir,i_atom) - force_a(idir)
                force(jkind)%hartree_2c(idir,j_atom) = &
                  force(jkind)%hartree_2c(idir,j_atom) + force_a(idir)
              END DO 

            ELSE

              CALL sumup_2c_integrals(z0_h_a,z0_h_b,z0_s,&
                                      nhard_a(0:), nhard_b(0:),ng0_s(0:),&
                                      l0_a,r0cut_a,QN_h_a,QN_s_a,&
                                      l0_b,r0cut_b,QN_h_b,QN_s_b,&
!                       l0_a,rpgf0_s_a,QN_h_a,QN_s_a,&
!                       l0_b,rpgf0_s_b,QN_h_b,QN_s_b,&
                                      rab,rab2,ieqj,ene_2c,Isum_a_car(:,iatom),&
                                      Isum_a_car(:,jatom),maxder,force_a,&
                                      vab,vv,ff(0:))
            END IF 

!           *** sum up in the energy of hartree local from 2 centers terms
            energy_hartree_2c = energy_hartree_2c + 0.5_dp*ieqj*ene_2c
            sab_2c_neighbor_node => next(sab_2c_neighbor_node)
          END DO  ! inode

          sab_2c_neighbor_list => next(sab_2c_neighbor_list)
        END DO  ! ilist

      END DO  ! jkind

    END DO  ! ikind

!   *** IF PARALLEL sum up the contributions to aVh2c0b(iatom) coming from
!   *** pair-interactions calculated on different processors
    CALL mp_sum(energy_hartree_2c,para_env%group)
    CALL mp_sum(Isum_a_car,para_env%group)

    IF(just_energy) GOTO 100
!   One center contributions from the 2-centers term
!   These contributions are added to the KS matrix together with the other
!   one center contributions, in the module qs_ks_atom
 ene_2c = 0.0_dp
    DO ikind = 1,nkind

       NULLIFY(atom_kind,orb_basis,atom_list)
       atom_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atom_kind,&
                            orb_basis_set=orb_basis,&
                            atom_list=atom_list,natom=nat,&
                            paw_atom=paw_atom_a)

       CALL  get_gto_basis_set(gto_basis_set=orb_basis,&
                                 nset =nset,maxso=maxso)
       nsotot = maxso * nset

       IF(.NOT. paw_atom_a) CYCLE

       NULLIFY(Qlm_gg)
       CALL get_rho0_mpole(rho0_mpole=rho0_mpole,ikind=ikind,&
                           l0_ikind=l0_a,Qlm_gg=Qlm_gg)

       nch0_sph = nsoset(l0_a)

!      Distribute the atoms of this kind
       num_pe = para_env%num_pe
       mepos  = para_env%mepos
       bo = get_limit( nat, num_pe, mepos )

       DO iat = bo(1), bo(2) !1,nat
         iatom = atom_list(iat)

         rho_atom => rho_atom_set(iatom)
         CALL get_rho_atom(rho_atom=rho_atom,&
              ga_Vlocal_gb_h=int_local_h,ga_Vlocal_gb_s=int_local_s)

         Isum_a_sph = 0.0_dp
         DO lshell = 0,l0_a
           DO is = 1,nso(lshell)
             iso = nsoset(lshell-1) + is 
             DO ic = 1,nco(lshell)
                ico = ncoset(lshell-1) + ic
                lx = indco(1,ico)
                ly = indco(2,ico)
                lz = indco(3,ico)
                Isum_a_sph(iso) =  Isum_a_sph(iso) + &
                   orbtramat(lshell)%c2s(is,ic)*Isum_a_car(ico,iatom)/&
                   SQRT(fourpi*dfac(2*lx-1)*dfac(2*ly-1)*dfac(2*lz-1)/&
                   dfac(2*lshell+1))
             END DO  ! ico
           END DO  ! iso
         END DO  ! lshell

         DO iso = 1,nch0_sph
           DO is = 1,nsotot

              DO ispin = 1,qs_env%dft_control%nspins
                int_local_h(ispin)%r_coef(1:nsotot,is) = &
                          int_local_h(ispin)%r_coef(1:nsotot,is) + &
                          Isum_a_sph(iso)*Qlm_gg(1:nsotot,is,iso)
                            
                int_local_s(ispin)%r_coef(1:nsotot,is) = &
                          int_local_s(ispin)%r_coef(1:nsotot,is) + &
                          Isum_a_sph(iso)*Qlm_gg(1:nsotot,is,iso)
              END DO
 
           END DO
         END DO

       END DO  ! iat
    END DO  ! ikind

100 CONTINUE

    DEALLOCATE(vab,vv,ff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
          "vab.vv.ff")

    DEALLOCATE(nhard_a, nhard_b, QN_h_a, QN_h_b, QN_s_a, QN_s_b, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "nhard_a, nhard_b, QN_h_a, QN_h_b, QN_s_a, QN_s_b")
    DEALLOCATE(Isum_a_car, Isum_a_sph, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
         "Isum_a_car, Isum_a_sph")

    IF (calculate_forces) THEN
       DEALLOCATE(atom_of_kind,STAT=istat)
       IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
            "atom_of_kind")
    END IF

    CALL timestop(0.0_dp,handle) 

  END SUBROUTINE integrate_vhg0_g0

!******************************************************************************
 
  SUBROUTINE sumup_2c_integrals(z0_h_a,z0_h_b,z0_s,ng0_h_a,ng0_h_b,ng0_s,&
                       l0_a,r0_a,QN_h_a,QN_s_a,&
                       l0_b,r0_b,QN_h_b,QN_s_b,&
                       rab,rab2,ieqj,ene_2c,Isum_a,Isum_b,maxder,force_a,&
                       vab,vv,ff) 
!   ---------------------------------------------------------------------------
    INTEGER, INTENT(IN)                   :: l0_a, l0_b
    INTEGER, INTENT(IN)                   :: maxder
    REAL(dp), INTENT(IN)                  :: z0_h_a, z0_h_b ,z0_s
    REAL(dp), INTENT(IN)                  :: r0_a
    REAL(dp), INTENT(IN)                  :: r0_b 
    REAL(dp), INTENT(IN)                  :: ieqj,rab2
    REAL(dp)                              :: ene_2c
    REAL(dp), DIMENSION(3), INTENT(IN)    :: rab
    REAL(dp), DIMENSION(3), INTENT(OUT)   :: force_a
    REAL(dp), DIMENSION(0:), INTENT(IN)   :: ng0_h_a, ng0_h_b, ng0_s
    REAL(dp), DIMENSION(:), INTENT(IN)    :: QN_h_a,QN_s_a
    REAL(dp), DIMENSION(:), INTENT(IN)    :: QN_h_b,QN_s_b
    REAL(dp), DIMENSION(:), INTENT(INOUT) :: Isum_a, Isum_b
    REAL(dp), DIMENSION(0:)               :: ff
    REAL(dp), DIMENSION(:,:)              :: vab
    REAL(dp), DIMENSION(:,:,:)            :: vv

    CHARACTER(LEN=*), PARAMETER           :: routine_name = "sumup_2c_integrals"

    INTEGER  :: icoa, icoap, icoam, icob, idir, l, la(3), lam(3), lap(3), &
                istat, nch_a, nch_a_p, nch_b
    REAL(dp) :: sumb, zx2, zx2_a
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: I_tmp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vab_plus
    REAL(dp), DIMENSION(1)                :: zet_a,zet_b,rad_a,rad_b

!   ---------------------------------------------------------------------------

    nch_a = ncoset(l0_a)
    nch_b = ncoset(l0_b)

    IF(maxder > 0) THEN

      nch_a_p = ncoset(l0_a+1)
      ALLOCATE(vab_plus(nch_a_p,nch_b), I_tmp(nch_a_p), STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
           "vab_plus,I_tmp",nch_a*nch_b*dp_size+nch_a_p*dp_size)

      ene_2c = 0.0_dp
      force_a(1:3) = 0.0_dp 
      zet_a(1)  = z0_h_a
      zet_b(1)  = z0_h_b
      rad_a(1) = r0_a
      rad_b(1) = r0_b      
      zx2_a = z0_h_a*2.0_dp
      CALL coulomb2(l0_a+1,1,zet_a,rad_a,0,l0_b,1,zet_b,rad_b,0,&
                    rab,rab2,vab,vv,ff(0:),maxder=maxder,vac_plus=vab_plus)

      I_tmp = 0.0_dp
      DO icoa = 1,nch_a
        l = indco(1,icoa)+indco(2,icoa)+indco(3,icoa)

        DO icob = 1,nch_b
          I_tmp(icoa) = I_tmp(icoa) + QN_h_b(icob)*vab(icoa,icob)
        END DO

        ene_2c = ene_2c + QN_h_a(icoa) * I_tmp(icoa)
        Isum_a(icoa) = Isum_a(icoa) +  ng0_h_a(l) * I_tmp(icoa)

      END DO
      IF(ieqj > 1.0_dp) THEN
        DO icob = 1,nch_b
          l = indco(1,icob)+indco(2,icob)+indco(3,icob)
          sumb = 0.0_dp
          DO icoa = 1,nch_a
            sumb = sumb + QN_h_a(icoa)*vab(icoa,icob)
          END DO
          Isum_b(icob) = Isum_b(icob) + ng0_h_b(l)*sumb
        END DO
      END IF
      DO icoa = ncoset(l0_a)+1,ncoset(l0_a+1)
        DO icob = 1,nch_b
          I_tmp(icoa) = I_tmp(icoa) + QN_h_b(icob)*vab_plus(icoa,icob)
        END DO 
      END DO
      
      DO idir = 1,3
        DO icoa = 1,nch_a
          la(1:3) = indco(1:3,icoa)
          lap(1:3) = la(1:3)
          lam(1:3) = la(1:3)
          lap(idir) = la(idir) + 1
          lam(idir) = la(idir) - 1
          icoap = coset(lap(1),lap(2),lap(3))

          IF(lam(idir) < 0) THEN

            force_a(idir) = force_a(idir) - QN_h_a(icoa)*& 
                                   zx2_a* I_tmp(icoap)
          ELSE
            icoam = coset(lam(1),lam(2),lam(3))

            force_a(idir) = force_a(idir) - QN_h_a(icoa)*& 
                   (zx2_a* I_tmp(icoap)-la(idir)*I_tmp(icoam))
          END IF

        END DO
      END DO  ! idir

      zet_a(1)  = z0_s
      zet_b(1)  = z0_s
      rad_a(1) = r0_a
      rad_b(1) = r0_b      
      zx2 = z0_s * 2.0_dp
      CALL coulomb2(l0_a+1,1,zet_a,rad_a,0,l0_b,1,zet_b,rad_b,0,&
                    rab,rab2,vab,vv,ff(0:),maxder=maxder,vac_plus=vab_plus)

      I_tmp = 0.0_dp
      DO icoa = 1,nch_a
        l = indco(1,icoa)+indco(2,icoa)+indco(3,icoa)

        DO icob = 1,nch_b
          I_tmp(icoa) = I_tmp(icoa) + QN_s_b(icob)*vab(icoa,icob)
        END DO

        ene_2c = ene_2c - QN_s_a(icoa) * I_tmp(icoa)

        Isum_a(icoa) = Isum_a(icoa) -  ng0_s(l) * I_tmp(icoa)

      END DO
      IF(ieqj > 1.0_dp) THEN
        DO icob = 1,nch_b
          l = indco(1,icob)+indco(2,icob)+indco(3,icob)
          sumb = 0.0_dp
          DO icoa = 1,nch_a
            sumb = sumb + QN_s_a(icoa)*vab(icoa,icob)
          END DO
          Isum_b(icob) = Isum_b(icob) - ng0_s(l)*sumb
        END DO
      END IF
      DO icoa = ncoset(l0_a)+1,ncoset(l0_a+1)
        DO icob = 1,nch_b
          I_tmp(icoa) = I_tmp(icoa) + QN_s_b(icob)*vab_plus(icoa,icob)
        END DO 
      END DO

      DO idir = 1,3
        DO icoa = 1,nch_a
          la(1:3) = indco(1:3,icoa)
          lap(1:3) = la(1:3)
          lam(1:3) = la(1:3)
          lap(idir) = la(idir) + 1
          lam(idir) = la(idir) - 1
          icoap = coset(lap(1),lap(2),lap(3))

          IF(lam(idir) < 0) THEN

            force_a(idir) = force_a(idir) + &
                             QN_s_a(icoa) * zx2 * I_tmp(icoap)
          ELSE
            icoam = coset(lam(1),lam(2),lam(3))

            force_a(idir) = force_a(idir) + &
               QN_s_a(icoa)*(zx2*I_tmp(icoap) - la(idir)*I_tmp(icoam))
          END IF

        END DO
      END DO  ! idir

      DEALLOCATE(vab_plus, I_tmp, STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
           "vab_plus,I_tmp")


    ELSE

      ene_2c = 0.0_dp

      zet_a(1)  = z0_h_a
      zet_b(1)  = z0_h_b
      rad_a(1) = r0_a
      rad_b(1) = r0_b      
      CALL coulomb2(l0_a,1,zet_a,rad_a,0,l0_b,1,zet_b,rad_b,0,&
                    rab,rab2,vab,vv,ff(0:))

      DO icoa = 1,nch_a
        l = indco(1,icoa)+indco(2,icoa)+indco(3,icoa)

        sumb = 0.0_dp
        DO icob = 1,nch_b
          sumb = sumb + QN_h_b(icob)*vab(icoa,icob)
        END DO

        ene_2c = ene_2c + QN_h_a(icoa) * sumb

        Isum_a(icoa) = Isum_a(icoa) +  ng0_h_a(l)*sumb
      END DO
      IF(ieqj > 1.0_dp) THEN
        DO icob = 1,nch_b
          l = indco(1,icob)+indco(2,icob)+indco(3,icob)
          sumb = 0.0_dp
          DO icoa = 1,nch_a
            sumb = sumb + QN_h_a(icoa)*vab(icoa,icob)
          END DO
          Isum_b(icob) = Isum_b(icob) + ng0_h_b(l)*sumb
        END DO
      END IF

      zet_a(1)  = z0_s
      zet_b(1)  = z0_s
      rad_a(1) = r0_a
      rad_b(1) = r0_b      
      CALL coulomb2(l0_a,1,zet_a,rad_a,0,l0_b,1,zet_b,rad_b,0,&
                    rab,rab2,vab,vv,ff(0:))

      DO icoa = 1,nch_a
        l = indco(1,icoa)+indco(2,icoa)+indco(3,icoa)

        sumb = 0.0_dp
        DO icob = 1,nch_b
          sumb = sumb + QN_s_b(icob)*vab(icoa,icob)
        END DO

        ene_2c = ene_2c - QN_s_a(icoa) * sumb

        Isum_a(icoa) = Isum_a(icoa) -  ng0_s(l)*sumb
      END DO
      IF(ieqj > 1.0_dp) THEN
        DO icob = 1,nch_b
          l = indco(1,icob)+indco(2,icob)+indco(3,icob)
          sumb = 0.0_dp
          DO icoa = 1,nch_a
            sumb = sumb + QN_s_a(icoa)*vab(icoa,icob)
          END DO
          Isum_b(icob) = Isum_b(icob) - ng0_s(l)*sumb
        END DO
      END IF

    END IF

  END SUBROUTINE sumup_2c_integrals

!******************************************************************************

END MODULE hartree_local_2centers

!******************************************************************************
