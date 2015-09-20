!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!
! *****************************************************************************
MODULE qs_rho_atom_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_p_type,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             gapw_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_get_block_p,&
                                             cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erf
  USE kinds,                           ONLY: dp
  USE kpoint_types,                    ONLY: get_kpoint_info,&
                                             kpoint_type
  USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                             get_number_of_lebedev_grid,&
                                             init_lebedev_grids,&
                                             lebedev_grid
  USE mathconstants,                   ONLY: fourpi,&
                                             pi
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE paw_proj_set_types,              ONLY: get_paw_proj_set,&
                                             paw_proj_set_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_grid_atom,                    ONLY: create_grid_atom,&
                                             grid_atom_type
  USE qs_harmonics_atom,               ONLY: create_harmonics_atom,&
                                             get_maxl_CG,&
                                             get_none0_cg_list,&
                                             harmonics_atom_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             get_qs_kind_set,&
                                             qs_kind_type
  USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                             neighbor_list_iterate,&
                                             neighbor_list_iterator_create,&
                                             neighbor_list_iterator_p_type,&
                                             neighbor_list_iterator_release,&
                                             neighbor_list_set_p_type
  USE qs_oce_methods,                  ONLY: proj_blk
  USE qs_oce_types,                    ONLY: oce_matrix_type
  USE qs_rho_atom_types,               ONLY: allocate_rho_atom_set,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE sap_kind_types,                  ONLY: alist_pre_align_blk,&
                                             alist_type,&
                                             get_alist
  USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                             clebsch_gordon_deallocate,&
                                             clebsch_gordon_init
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  USE whittaker,                       ONLY: whittaker_c0a,&
                                             whittaker_ci
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_atom_methods'

! *** Public subroutines ***

  PUBLIC :: allocate_rho_atom_internals, &
            calculate_rho_atom, &
            calculate_rho_atom_coeff, &
            init_rho_atom

 CONTAINS

! *****************************************************************************
!> \brief ...
!> \param para_env ...
!> \param rho_atom_set ...
!> \param qs_kind ...
!> \param atom_list ...
!> \param grid_atom ...
!> \param natom ...
!> \param nspins ...
!> \param tot_rho1_h ...
!> \param tot_rho1_s ...
! *****************************************************************************
  SUBROUTINE calculate_rho_atom(para_env,rho_atom_set,qs_kind,atom_list,&
                                grid_atom, natom, nspins, tot_rho1_h, tot_rho1_s)

    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(qs_kind_type), INTENT(IN)           :: qs_kind
    INTEGER, DIMENSION(:), INTENT(IN)        :: atom_list
    TYPE(grid_atom_type), INTENT(IN)         :: grid_atom
    INTEGER, INTENT(IN)                      :: natom, nspins
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: tot_rho1_h, tot_rho1_s

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: damax_iso_not0_local, dbmax_iso_not0_local, dmax_iso_not0, &
      handle, i, i1, i2, iat, iatom, icg, ipgf1, ipgf2, ir, iset1, iset2, &
      iso, iso1, iso1_first, iso1_last, iso2, iso2_first, iso2_last, istat, &
      j, l, l1, l2, l_iso, l_sub, l_sum, lmax12, lmax_expansion, lmin12, m1, &
      m1s, m2, m2s, max_iso_not0, max_iso_not0_local, max_s_harm, maxl, &
      maxso, mepos, n1s, n2s, nr, nset, num_pe, size1, size2
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cg_n_list, dacg_n_list, &
                                                dbcg_n_list
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cg_list, dacg_list, dbcg_list
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf, o2nindex
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: done_vgg
    REAL(dp)                                 :: c1, c2, rho_h, rho_s, &
                                                root_zet12, zet12
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: dd, erf_zet12, g1, g2, gg0, &
                                                int1, int2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: CPCH_sphere, CPCS_sphere, &
                                                dgg, gg, gg_lm1
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: vgg
    REAL(dp), DIMENSION(:, :), POINTER       :: coeff, zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG, my_dCG
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: my_CG_dxyz
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(paw_proj_set_type), POINTER         :: paw_proj

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(orb_basis)
    NULLIFY(harmonics)
    NULLIFY(lmin,lmax,npgf,zet,my_CG,my_dCG,my_CG_dxyz,coeff)

    CALL get_qs_kind(qs_kind,basis_set=orb_basis,&
                         paw_proj_set=paw_proj,harmonics=harmonics)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxl=maxl,npgf=npgf,nset=nset,zet=zet,&
                           maxso=maxso)

    CALL get_paw_proj_set(paw_proj_set=paw_proj,o2nindex=o2nindex)

    max_iso_not0 = harmonics%max_iso_not0
    dmax_iso_not0 = harmonics%dmax_iso_not0
    max_s_harm = harmonics%max_s_harm

    nr           = grid_atom%nr
    lmax_expansion = indso(1,max_iso_not0)
!   Distribute the atoms of this kind
    num_pe = para_env%num_pe
    mepos  = para_env%mepos
    bo = get_limit( natom, num_pe, mepos )

    my_CG        => harmonics%my_CG
    my_dCG       => harmonics%my_dCG
    my_CG_dxyz   => harmonics%my_CG_dxyz

    ALLOCATE(CPCH_sphere(nsoset(maxl),nsoset(maxl)),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(CPCS_sphere(nsoset(maxl),nsoset(maxl)),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(g1(nr),g2(nr),gg0(nr),gg(nr,0:2*maxl),dgg(nr,0:2*maxl),&
             gg_lm1(nr,0:2*maxl),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(erf_zet12(nr),vgg(nr,0:2*maxl,0:indso(1,max_iso_not0)),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(done_vgg(0:2*maxl,0:indso(1,max_iso_not0)),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(dd(nr),int1(nr),int2(nr),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    ALLOCATE(cg_list(2,nsoset(maxl)**2,max_s_harm),cg_n_list(max_s_harm),&
             dacg_list(2,nsoset(maxl)**2,max_s_harm),dacg_n_list(max_s_harm),&
             dbcg_list(2,nsoset(maxl)**2,max_s_harm),dbcg_n_list(max_s_harm),&
             STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)

    DO iat = bo(1),bo(2)
      iatom = atom_list(iat)
      DO i=1,nspins
        IF(.NOT. ASSOCIATED(rho_atom_set(iatom)%rho_rad_h(i)%r_coef)) THEN
            CALL allocate_rho_atom_rad(rho_atom_set,iatom,i,nr,max_iso_not0)
        ELSE
            CALL set2zero_rho_atom_rad(rho_atom_set,iatom,i)
        ENDIF
      END DO
    END DO

    m1s = 0
    DO iset1 = 1,nset
        m2s = 0
        DO iset2 = 1,nset

           CALL get_none0_cg_list(my_CG,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                  max_s_harm,lmax_expansion,cg_list,cg_n_list,max_iso_not0_local)
           CPPrecondition(max_iso_not0_local.LE.max_iso_not0,cp_failure_level,routineP,failure)
           CALL get_none0_cg_list(my_CG_dxyz,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                  max_s_harm,lmax_expansion,dacg_list,dacg_n_list,damax_iso_not0_local)
           CALL get_none0_cg_list(my_dCG,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                  max_s_harm,lmax_expansion,dbcg_list,dbcg_n_list,dbmax_iso_not0_local)

           n1s = nsoset(lmax(iset1))
           DO ipgf1  = 1,npgf(iset1)
              iso1_first = nsoset(lmin(iset1)-1)+1+n1s*(ipgf1-1)+m1s
              iso1_last  = nsoset(lmax(iset1))+n1s*(ipgf1-1)+m1s
              size1 = iso1_last - iso1_first + 1
              iso1_first = o2nindex(iso1_first)
              iso1_last  = o2nindex(iso1_last)
              i1 = iso1_last - iso1_first + 1
              CPPrecondition(size1==i1,cp_failure_level,routineP,failure)
              i1 = nsoset(lmin(iset1)-1)+1

              g1(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))

              n2s=nsoset(lmax(iset2))
              DO ipgf2 = 1,npgf(iset2)
                 iso2_first = nsoset(lmin(iset2)-1)+1+n2s*(ipgf2-1)+m2s
                 iso2_last  = nsoset(lmax(iset2))+n2s*(ipgf2-1)+m2s
                 size2 = iso2_last - iso2_first + 1
                 iso2_first = o2nindex(iso2_first)
                 iso2_last  = o2nindex(iso2_last)
                 i2 = iso2_last - iso2_first + 1
                 CPPrecondition(size2==i2,cp_failure_level,routineP,failure)
                 i2 = nsoset(lmin(iset2)-1)+1

                 g2(1:nr) = EXP(-zet(ipgf2,iset2)*grid_atom%rad2(1:nr))
                 lmin12 = lmin(iset1)+lmin(iset2)
                 lmax12 = lmax(iset1)+lmax(iset2)

                 zet12 = zet(ipgf1,iset1)+zet(ipgf2,iset2)
                 root_zet12 = SQRT(zet(ipgf1,iset1)+zet(ipgf2,iset2))
                 DO ir = 1,nr
                    erf_zet12(ir) = erf(root_zet12*grid_atom%rad(ir))
                 END DO

                 gg = 0.0_dp
                 dgg = 0.0_dp
                 gg_lm1 = 0.0_dp
                 vgg = 0.0_dp
                 done_vgg = .FALSE.
                 ! reduce the number of terms in the expansion local densities
                 IF(lmin12 .LE. lmax_expansion) THEN
                    IF (lmin12 == 0) THEN
                       gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
                       gg_lm1(1:nr,lmin12) = 0.0_dp
                       gg0(1:nr) = gg(1:nr,lmin12)
                    ELSE
                       gg0(1:nr) = g1(1:nr)*g2(1:nr)
                       gg(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12)*g1(1:nr)*g2(1:nr)
                       gg_lm1(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12-1)*g1(1:nr)*g2(1:nr)
                    END IF

                    ! reduce the number of terms in the expansion local densities
                    IF(lmax12 .GT. lmax_expansion) lmax12 = lmax_expansion

                    DO l=lmin12+1,lmax12
                       gg(1:nr,l) = grid_atom%rad(1:nr)*gg(1:nr,l-1)
                       gg_lm1(1:nr,l) = gg(1:nr,l-1)
                       dgg(1:nr,l-1) = dgg(1:nr,l-1) -2.0_dp*(zet(ipgf1,iset1)+ &
                            zet(ipgf2,iset2))*gg(1:nr,l)

                    END DO
                    dgg(1:nr,lmax12) = dgg(1:nr,lmax12) -2.0_dp*(zet(ipgf1,iset1)+ &
                         zet(ipgf2,iset2))*grid_atom%rad(1:nr)*&
                         gg(1:nr,lmax12)

                    c2 = SQRT(pi*pi*pi/(zet12*zet12*zet12))

                    DO iso = 1,max_iso_not0_local
                       l_iso = indso(1,iso)
                       c1 = fourpi/(2._dp*REAL(l_iso,dp)+1._dp)
                       DO icg = 1,cg_n_list(iso)
                          iso1 = cg_list(1,icg,iso)
                          iso2 = cg_list(2,icg,iso)

                          l = indso(1,iso1) + indso(1,iso2)
                          CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,failure)
                          IF(.NOT. failure) THEN

                             IF(done_vgg(l,l_iso)) CYCLE
                             L_sum = l + l_iso
                             L_sub = l - l_iso

                             IF(l_sum == 0) THEN
                                vgg(1:nr,l,l_iso) = erf_zet12(1:nr)* grid_atom%oorad2l(1:nr,1) * c2
                             ELSE
                                CALL whittaker_c0a(int1,grid_atom%rad,gg0,erf_zet12,zet12,l,l_iso,nr)
                                CALL whittaker_ci(int2,grid_atom%rad,gg0,zet12,L_sub,nr)

                                DO ir = 1,nr
                                   int2(ir) = grid_atom%rad2l(ir,l_iso)*int2(ir)
                                   vgg(ir,l,l_iso) = c1 * (int1(ir) + int2(ir))
                                END DO
                             END IF
                             done_vgg(l,l_iso) = .TRUE.
                          END IF ! failure

                       END DO
                    END DO
                 END IF ! lmax_expansion

                 DO iat = bo(1),bo(2)
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

                       DO iso = 1,max_iso_not0_local
                          l_iso = indso(1,iso)
                          DO icg = 1,cg_n_list(iso)
                             iso1 = cg_list(1,icg,iso)
                             iso2 = cg_list(2,icg,iso)

                             l1 = indso(1,iso1)
                             l2 = indso(1,iso2)
                             m1 = indso(2,iso1)
                             m2 = indso(2,iso2)

                             dd(1:nr) = REAL(l1*l2,dp)/grid_atom%rad2(1:nr) &
                                  - 2._dp*REAL(l1,dp)*zet(ipgf2,iset2) &
                                  - 2._dp*REAL(l2,dp)*zet(ipgf1,iset1) &
                                  + 4._dp*zet(ipgf1,iset1)*zet(ipgf2,iset2)*grid_atom%rad2(1:nr)
                             l = indso(1,iso1) + indso(1,iso2)
                             CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,failure)
                             IF(.NOT. failure) THEN
                                rho_atom_set(iatom)%rho_rad_h(i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%rho_rad_h(i)%r_coef(1:nr,iso) +&
                                     gg(1:nr,l)*CPCH_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                                rho_atom_set(iatom)%rho_rad_s(i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%rho_rad_s(i)%r_coef(1:nr,iso) +&
                                     gg(1:nr,l)*CPCS_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                                rho_atom_set(iatom)%drho_rad_h(i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%drho_rad_h(i)%r_coef(1:nr,iso) +&
                                     dgg(1:nr,l)*CPCH_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                                rho_atom_set(iatom)%drho_rad_s(i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%drho_rad_s(i)%r_coef(1:nr,iso) +&
                                     dgg(1:nr,l)*CPCS_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                                rho_atom_set(iatom)%vrho_rad_h(i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%vrho_rad_h(i)%r_coef(1:nr,iso) +&
                                     vgg(1:nr,l,l_iso)*CPCH_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                                rho_atom_set(iatom)%vrho_rad_s(i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%vrho_rad_s(i)%r_coef(1:nr,iso) +&
                                     vgg(1:nr,l,l_iso)*CPCS_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                                rho_atom_set(iatom)%trho_rad_h(1,i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%trho_rad_h(1,i)%r_coef(1:nr,iso) +&
                                     0.5_dp*gg(1:nr,l)*dd(1:nr)*CPCH_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                                rho_atom_set(iatom)%trho_rad_s(1,i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%trho_rad_s(1,i)%r_coef(1:nr,iso) +&
                                     0.5_dp*gg(1:nr,l)*dd(1:nr)*CPCS_sphere(iso1,iso2)*my_CG(iso1,iso2,iso)

                                rho_atom_set(iatom)%trho_rad_h(3,i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%trho_rad_h(3,i)%r_coef(1:nr,iso) +&
                                     0.5_dp*REAL(m1*m2,dp)*gg(1:nr,l)*CPCH_sphere(iso1,iso2)*&
                                     my_CG(iso1,iso2,iso)/grid_atom%rad2(1:nr)

                                rho_atom_set(iatom)%trho_rad_s(3,i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%trho_rad_s(3,i)%r_coef(1:nr,iso) +&
                                     0.5_dp*REAL(m1*m2,dp)*gg(1:nr,l)*CPCS_sphere(iso1,iso2)*&
                                     my_CG(iso1,iso2,iso)/grid_atom%rad2(1:nr)

                             END IF !failure

                          ENDDO  ! icg

                       ENDDO ! iso

                       DO iso = 1,max_iso_not0 !damax_iso_not0_local
                          l_iso = indso(1,iso)
                          DO icg = 1,dacg_n_list(iso)
                             iso1 = dacg_list(1,icg,iso)
                             iso2 = dacg_list(2,icg,iso)
                             l = indso(1,iso1) + indso(1,iso2)
                             CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,failure)
                             IF(.NOT. failure) THEN
                                DO j = 1,3
                                   rho_atom_set(iatom)%rho_rad_h_d(j,i)%r_coef(1:nr,iso) =&
                                        rho_atom_set(iatom)%rho_rad_h_d(j,i)%r_coef(1:nr,iso) +&
                                        gg_lm1(1:nr,l)*CPCH_sphere(iso1,iso2)*my_CG_dxyz(j,iso1,iso2,iso)

                                   rho_atom_set(iatom)%rho_rad_s_d(j,i)%r_coef(1:nr,iso) =&
                                        rho_atom_set(iatom)%rho_rad_s_d(j,i)%r_coef(1:nr,iso) +&
                                        gg_lm1(1:nr,l)*CPCS_sphere(iso1,iso2)*my_CG_dxyz(j,iso1,iso2,iso)
                                END DO
                             END IF !failure
                          END DO  ! icg

                       ENDDO  ! iso

                       !Expansion for tau functionals (polar angle term)
                       !We ignore the contributions iso > max_iso_not0 (should be dmax_iso_not0)

                       DO iso = 1,max_iso_not0!dbmax_iso_not0_local
                          l_iso = indso(1,iso)
                          DO icg = 1,dbcg_n_list(iso)
                             iso1 = dbcg_list(1,icg,iso)
                             iso2 = dbcg_list(2,icg,iso)

                             l = indso(1,iso1) + indso(1,iso2)
                             CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,failure)
                             IF(.NOT.failure) THEN
                                rho_atom_set(iatom)%trho_rad_h(2,i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%trho_rad_h(2,i)%r_coef(1:nr,iso) +&
                                     0.5_dp*gg(1:nr,l)*CPCH_sphere(iso1,iso2)*my_dCG(iso1,iso2,iso)/grid_atom%rad2(1:nr)

                                rho_atom_set(iatom)%trho_rad_s(2,i)%r_coef(1:nr,iso) =&
                                     rho_atom_set(iatom)%trho_rad_s(2,i)%r_coef(1:nr,iso) +&
                                     0.5_dp*gg(1:nr,l)*CPCS_sphere(iso1,iso2)*my_dCG(iso1,iso2,iso)/grid_atom%rad2(1:nr)
                             END IF !failure
                          END DO  ! icg
                       ENDDO  ! iso

                    ENDDO ! i
                 ENDDO  ! iat

              ENDDO ! ipgf2
           ENDDO   ! ipgf1
           m2s = m2s+maxso
        ENDDO  ! iset2
        m1s = m1s+maxso
    ENDDO    ! iset1

    DO iat = bo(1),bo(2)
      iatom = atom_list(iat)

      DO i = 1,nspins

        DO iso = 1,max_iso_not0
          rho_s = 0.0_dp
          rho_h = 0.0_dp
          DO ir = 1,nr
            rho_h = rho_h + &
                    rho_atom_set(iatom)%rho_rad_h(i)%r_coef(ir,iso)*grid_atom%wr(ir)
            rho_s = rho_s + &
                    rho_atom_set(iatom)%rho_rad_s(i)%r_coef(ir,iso)*grid_atom%wr(ir)
          END DO  ! ir
          tot_rho1_h(i) = tot_rho1_h(i) + rho_h*harmonics%slm_int(iso)
          tot_rho1_s(i) = tot_rho1_s(i) + rho_s*harmonics%slm_int(iso)
        END DO  ! iso

      END DO  ! ispin

    END DO  ! iat

    DEALLOCATE(CPCH_sphere,CPCS_sphere,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    DEALLOCATE(g1,g2,gg0,gg,gg_lm1,dgg,vgg,done_vgg,erf_zet12,int1,int2,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    DEALLOCATE(cg_list,cg_n_list,dacg_list,dacg_n_list,dbcg_list,dbcg_n_list,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)

    CALL timestop(handle)

  END SUBROUTINE calculate_rho_atom

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho_ao ...
!> \param rho_atom_set ...
! *****************************************************************************
  SUBROUTINE calculate_rho_atom_coeff(qs_env, rho_ao, rho_atom_set)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: rho_ao
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_atom_set

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_atom_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2), handle, i, iac, iatom, ibc, icol, ikind, img, irow, &
      ispin, istat, jatom, jkind, kac, katom, kbc, kkind, max_gau, max_nsgf, &
      maxsoc, mepos, n_cont_a, n_cont_b, nat_kind, natom, nimages, nkind, &
      nsatbas, nsetc, nsgfa, nsgfb, nsoctot, nspins, num_pe
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: kind_of
    INTEGER, DIMENSION(3)                    :: cell_b
    INTEGER, DIMENSION(:), POINTER           :: a_list, list_a, list_b
    INTEGER, DIMENSION(:, :, :), POINTER     :: cell_to_index
    LOGICAL                                  :: dista, distab, distb, &
                                                failure, found, paw_atom
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: p_matrix
    REAL(KIND=dp)                            :: eps_cpc, factor, pmax
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: CPC, PC1
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: C_coeff_hh_a, C_coeff_hh_b, &
                                                C_coeff_ss_a, C_coeff_ss_b, &
                                                r_coef_h, r_coef_s
    TYPE(alist_type), POINTER                :: alist_ac, alist_bc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_p_type), &
      DIMENSION(:), POINTER                  :: basis_set_list
    TYPE(gto_basis_set_type), POINTER        :: basis_set_a, basis_set_b, &
                                                orb_basis_set
    TYPE(kpoint_type), POINTER               :: kpoints
    TYPE(neighbor_list_iterator_p_type), &
      DIMENSION(:), POINTER                  :: nl_iterator
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(paw_proj_set_type), POINTER         :: paw_proj
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: p_block_spin
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(dft_control)
    NULLIFY(qs_kind_set)
    NULLIFY(oce)
    NULLIFY(sab_orb)
    NULLIFY(rho_atom)
    NULLIFY(orb_basis_set)

    CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            atomic_kind_set=atomic_kind_set,&
            qs_kind_set=qs_kind_set,&
            para_env=para_env,&
            oce=oce,sab_orb=sab_orb,&
            rho_atom_set=rho_atom)

    eps_cpc = dft_control%qs_control%gapw_control%eps_cpc

    IF (PRESENT(rho_atom_set)) rho_atom => rho_atom_set

    nspins = dft_control%nspins
    nimages = dft_control%nimages

    NULLIFY(cell_to_index)
    IF (nimages > 1) THEN
       CALL get_qs_env(qs_env=qs_env,kpoints=kpoints)
       CALL get_kpoint_info(kpoint=kpoints,cell_to_index=cell_to_index)
    END IF

    ALLOCATE(p_block_spin(nspins),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    DO ispin = 1,nspins
      NULLIFY(p_block_spin(ispin)%r_coef)
    END DO

    CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
    CALL get_qs_kind_set(qs_kind_set, maxsgf=max_nsgf, maxgtops=max_gau)

    ALLOCATE (PC1(max_nsgf*max_gau),CPC(max_gau*max_gau),p_matrix(max_nsgf,max_nsgf),&
              STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)

    nkind = SIZE(atomic_kind_set)
!   Inizialize to 0 the CPC coefficients and the local density arrays
    DO ikind = 1, nkind
      CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=a_list,natom=nat_kind)
      CALL get_qs_kind(qs_kind_set(ikind),paw_atom=paw_atom)

      IF(.NOT. paw_atom) CYCLE
      DO i = 1,nat_kind
        iatom = a_list(i)
        DO ispin = 1, nspins
          rho_atom(iatom)%cpc_h(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%cpc_s(ispin)%r_coef = 0.0_dp
        ENDDO  ! ispin
      ENDDO  ! i

      num_pe = para_env%num_pe
      mepos = para_env%mepos
      bo = get_limit( nat_kind, num_pe, mepos )
      DO i = bo(1), bo(2)
        iatom = a_list(i)
        DO ispin = 1, nspins
          rho_atom(iatom)%ga_Vlocal_gb_h(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%ga_Vlocal_gb_s(ispin)%r_coef = 0.0_dp
        ENDDO  ! ispin
      ENDDO  ! i
    END DO ! ikind

    ALLOCATE (basis_set_list(nkind),STAT=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,failure)
    DO ikind=1,nkind
      CALL get_qs_kind(qs_kind_set(ikind),basis_set=basis_set_a)
      IF (ASSOCIATED(basis_set_a)) THEN
        basis_set_list(ikind)%gto_basis_set => basis_set_a
      ELSE
        NULLIFY(basis_set_list(ikind)%gto_basis_set)
      END IF
    END DO
    CALL neighbor_list_iterator_create(nl_iterator,sab_orb)
    DO WHILE (neighbor_list_iterate(nl_iterator)==0)
       CALL get_iterator_info(nl_iterator,ikind=ikind,jkind=jkind,&
                              iatom=iatom,jatom=jatom,cell=cell_b,r=rab)
       basis_set_a => basis_set_list(ikind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_a)) CYCLE
       basis_set_b => basis_set_list(jkind)%gto_basis_set
       IF (.NOT.ASSOCIATED(basis_set_b)) CYCLE
       nsgfa = basis_set_a%nsgf
       nsgfb = basis_set_b%nsgf

       pmax = 0._dp
       IF (iatom <= jatom) THEN
          irow = iatom
          icol = jatom
       ELSE
          irow = jatom
          icol = iatom
       ENDIF

       IF(nimages > 1) THEN
          img = cell_to_index(cell_b(1),cell_b(2),cell_b(3))
          CPPostcondition(img > 0,cp_failure_level,routineP,failure)
       ELSE
          img = 1
       END IF

       DO ispin = 1,nspins
          NULLIFY(p_block_spin(ispin)%r_coef)
          CALL cp_dbcsr_get_block_p(matrix=rho_ao(ispin,img)%matrix,&
               row=irow,col=icol,BLOCK=p_block_spin(ispin)%r_coef,&
               found=found)
          pmax = pmax+MAXVAL(ABS(p_block_spin(ispin)%r_coef))
       ENDDO

       DO kkind = 1 ,nkind
         CALL get_qs_kind(qs_kind_set(kkind),&
                          basis_set=orb_basis_set,&
                          paw_proj_set=paw_proj,paw_atom=paw_atom)

         CALL get_gto_basis_set(gto_basis_set=orb_basis_set,nset =nsetc,maxso=maxsoc)

         IF(.NOT. paw_atom)  CYCLE

         CALL get_paw_proj_set(paw_proj_set=paw_proj,nsatbas=nsatbas)
         nsoctot = nsatbas

         iac = ikind + nkind*(kkind - 1)
         ibc = jkind + nkind*(kkind - 1)

         IF (.NOT.ASSOCIATED(oce%intac(iac)%alist)) CYCLE
         IF (.NOT.ASSOCIATED(oce%intac(ibc)%alist)) CYCLE
         CALL get_alist(oce%intac(iac), alist_ac, iatom)
         CALL get_alist(oce%intac(ibc), alist_bc, jatom)
         IF (.NOT.ASSOCIATED(alist_ac)) CYCLE
         IF (.NOT.ASSOCIATED(alist_bc)) CYCLE

         DO kac=1,alist_ac%nclist
           DO kbc=1,alist_bc%nclist
             IF ( alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom ) CYCLE
             IF ( ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0) ) THEN
               IF(pmax*alist_bc%clist(kbc)%maxac*alist_ac%clist(kac)%maxac < eps_cpc) CYCLE

               n_cont_a = alist_ac%clist(kac)%nsgf_cnt
               n_cont_b = alist_bc%clist(kbc)%nsgf_cnt
               IF(n_cont_a.EQ.0.OR.n_cont_b.EQ.0) CYCLE
               list_a => alist_ac%clist(kac)%sgf_list
               list_b => alist_bc%clist(kbc)%sgf_list

               katom = alist_ac%clist(kac)%catom

               IF(iatom==katom .AND. ALL(alist_ac%clist(kac)%cell == 0)) THEN
                 C_coeff_hh_a => alist_ac%clist(kac)%achint(:,:,1)
                 C_coeff_ss_a => alist_ac%clist(kac)%acint(:,:,1)
                 dista=.FALSE.
               ELSE
                 C_coeff_hh_a => alist_ac%clist(kac)%acint(:,:,1)
                 C_coeff_ss_a => alist_ac%clist(kac)%acint(:,:,1)
                 dista=.TRUE.
               END IF
               IF(jatom==katom .AND. ALL(alist_bc%clist(kbc)%cell == 0)) THEN
                 C_coeff_hh_b => alist_bc%clist(kbc)%achint(:,:,1)
                 C_coeff_ss_b => alist_bc%clist(kbc)%acint(:,:,1)
                 distb=.FALSE.
               ELSE
                 C_coeff_hh_b => alist_bc%clist(kbc)%acint(:,:,1)
                 C_coeff_ss_b => alist_bc%clist(kbc)%acint(:,:,1)
                 distb=.TRUE.
               END IF

               distab = dista.AND.distb

               DO ispin = 1,nspins

                  IF (iatom <= jatom) THEN
                     CALL alist_pre_align_blk(p_block_spin(ispin)%r_coef,&
                          SIZE(p_block_spin(ispin)%r_coef,1),p_matrix,SIZE(p_matrix,1),&
                          list_a,n_cont_a,list_b,n_cont_b)
                  ELSE
                     CALL alist_pre_align_blk(p_block_spin(ispin)%r_coef,&
                          SIZE(p_block_spin(ispin)%r_coef,1),p_matrix,SIZE(p_matrix,1),&
                          list_b,n_cont_b,list_a,n_cont_a)
                  ENDIF

                  factor = 1.0_dp
                  IF(iatom == jatom) factor = 0.5_dp

                  r_coef_h => rho_atom(katom)%cpc_h(ispin)%r_coef
                  r_coef_s => rho_atom(katom)%cpc_s(ispin)%r_coef

                  IF(iatom <= jatom) THEN
                     CALL proj_blk(C_coeff_hh_a,C_coeff_ss_a,n_cont_a,&
                                   C_coeff_hh_b,C_coeff_ss_b,n_cont_b,&
                                   p_matrix,max_nsgf,r_coef_h,r_coef_s,nsoctot,&
                                   PC1,CPC,factor,distab)
                  ELSE
                     CALL proj_blk(C_coeff_hh_b,C_coeff_ss_b,n_cont_b,&
                                   C_coeff_hh_a,C_coeff_ss_a,n_cont_a,&
                                   p_matrix,max_nsgf,r_coef_h,r_coef_s,nsoctot,&
                                   PC1,CPC,factor,distab)
                  ENDIF

               END DO
               EXIT !search loop over jatom-katom list
             END IF
           END DO
         END DO
       END DO
    END DO
    CALL neighbor_list_iterator_release(nl_iterator)

    ALLOCATE (kind_of(natom),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    CALL get_atomic_kind_set(atomic_kind_set,kind_of=kind_of)

    DO iatom = 1,natom
      ikind = kind_of(iatom)

       DO ispin = 1,nspins
         IF(ASSOCIATED(rho_atom(iatom)%cpc_h(ispin)%r_coef)) THEN
           CALL mp_sum(rho_atom(iatom)%cpc_h(ispin)%r_coef, para_env%group)
           CALL mp_sum(rho_atom(iatom)%cpc_s(ispin)%r_coef, para_env%group)
           r_coef_h => rho_atom(iatom)%cpc_h(ispin)%r_coef
           r_coef_s => rho_atom(iatom)%cpc_s(ispin)%r_coef
           r_coef_h(:,:) = r_coef_h(:,:)+TRANSPOSE(r_coef_h(:,:))
           r_coef_s(:,:) = r_coef_s(:,:)+TRANSPOSE(r_coef_s(:,:))
         END IF
       ENDDO

    END DO

    DEALLOCATE(PC1,CPC,p_block_spin,p_matrix,kind_of,basis_set_list,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)

    CALL timestop(handle)

  END SUBROUTINE calculate_rho_atom_coeff

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gapw_control ...
! *****************************************************************************
  SUBROUTINE init_rho_atom(qs_env,gapw_control)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(gapw_control_type), POINTER         :: gapw_control

    CHARACTER(len=*), PARAMETER :: routineN = 'init_rho_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ikind, il, iso, iso1, iso2, istat, l1, l1l2, l2, la, &
      lc1, lc2, lcleb, ll, llmax, lmax_sphere, lp, m1, m2, max_s_harm, &
      max_s_set, maxl, maxlgto, maxs, mm, mp, na, nat, natom, nr, nspins, &
      quadrature
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, paw_atom
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rga
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(rho_atom_set)
    NULLIFY(qs_kind_set, atomic_kind_set)
    NULLIFY(dft_control)
    NULLIFY(orb_basis_set)
    NULLIFY(my_CG, grid_atom, harmonics,atom_list)

    CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set,&
                    dft_control=dft_control)
    CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
    CALL get_qs_kind_set(qs_kind_set,maxlgto=maxlgto)

    nspins = dft_control%nspins

    lmax_sphere = gapw_control%lmax_sphere

    llmax = MIN(lmax_sphere, 2*maxlgto)
    max_s_harm = nsoset(llmax)
    max_s_set = nsoset(maxlgto)

    lcleb = MAX(llmax,2*maxlgto,1)

!   *** allocate calculate the CG coefficients up to the maxl ***
    CALL clebsch_gordon_init ( lcleb )
    CALL reallocate(my_CG,1,max_s_set,1,max_s_set,1,max_s_harm)

    ALLOCATE(rga(lcleb,2),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
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
    DEALLOCATE(rga,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,failure)
    CALL clebsch_gordon_deallocate()

!   *** initialize the Lebedev grids ***
    CALL init_lebedev_grids()
    quadrature=gapw_control%quadrature

    DO ikind = 1, SIZE(atomic_kind_set)
      CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=atom_list,natom=nat)
      CALL get_qs_kind(qs_kind_set(ikind),&
                       paw_atom=paw_atom,&
                       grid_atom=grid_atom,&
                       harmonics=harmonics,&
                       ngrid_rad=nr,ngrid_ang=na)

!     *** determine the Lebedev grid for this kind ***

      ll = get_number_of_lebedev_grid(n=na)
      na = lebedev_grid(ll)%n
      la = lebedev_grid(ll)%l
      grid_atom%ng_sphere = na
      grid_atom%nr = nr

      IF (llmax > la) THEN
        WRITE(*,'(/,72("*"))')
        WRITE(*,'(T2,A,T66,I4)') &
          "WARNING: the lebedev grid is built for angular momentum l up to ",la,&
          "         the max l of spherical harmonics is larger, l_max = ", llmax,&
          "         good integration is guaranteed only for l <= ",la
        WRITE(*,'(72("*"),/)')
      END IF

!     *** calculate the radial grid ***
      CALL create_grid_atom(grid_atom,nr,na,llmax,ll,quadrature)

!     *** calculate the spherical harmonics on the grid ***

      CALL get_qs_kind(qs_kind_set(ikind),basis_set=orb_basis_set)
      CALL get_gto_basis_set(gto_basis_set=orb_basis_set,maxl=maxl)
      maxs=nsoset(maxl)
      CALL create_harmonics_atom(harmonics,&
                   my_CG,na,llmax,maxs,max_s_harm,ll,grid_atom%wa,&
                   grid_atom%pol,grid_atom%cos_pol,&
                   grid_atom%cotan_azi)
      CALL get_maxl_CG(harmonics,orb_basis_set,llmax,max_s_harm)

    END DO

    CALL deallocate_lebedev_grids()
    DEALLOCATE(my_CG)

    CALL allocate_rho_atom_internals(qs_env, rho_atom_set)

    CALL set_qs_env(qs_env=qs_env, rho_atom_set=rho_atom_set)

    CALL timestop(handle)

  END SUBROUTINE init_rho_atom

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho_atom_set ...
! *****************************************************************************
  SUBROUTINE allocate_rho_atom_internals(qs_env,rho_atom_set)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rho_atom_internals', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2), handle, iat, iatom, ikind, ispin, istat, j, &
      max_iso_not0, mepos, nat, natom, nr, nsatbas, nsotot, nspins, num_pe
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, paw_atom
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(paw_proj_set_type), POINTER         :: paw_proj
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN,handle)

    NULLIFY(qs_kind_set,para_env,rho_atom_set,dft_control)
    failure = .FALSE.

    CALL get_qs_env(qs_env,&
                    qs_kind_set=qs_kind_set, &
                    atomic_kind_set=atomic_kind_set,&
                    para_env=para_env,&
                    dft_control=dft_control)

    CALL get_atomic_kind_set(atomic_kind_set, natom=natom)

    nspins = dft_control%nspins

    CALL allocate_rho_atom_set(rho_atom_set,natom)

    DO ikind = 1, SIZE(atomic_kind_set)

       NULLIFY(atom_list, harmonics)
       CALL get_atomic_kind(atomic_kind_set(ikind),atom_list=atom_list,natom=nat)
       CALL get_qs_kind(qs_kind_set(ikind),&
                        paw_proj_set=paw_proj,&
                        paw_atom=paw_atom,&
                        harmonics=harmonics,&
                        ngrid_rad=nr)

       IF (paw_atom) THEN
         CALL get_paw_proj_set(paw_proj_set=paw_proj,nsatbas=nsatbas,nsotot=nsotot)
       END IF

       max_iso_not0 = harmonics%max_iso_not0
       DO iat = 1,nat
          iatom = atom_list(iat)
          !       *** allocate the radial density for each LM,for each atom ***

          ALLOCATE (rho_atom_set(iatom)%rho_rad_h(nspins),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,failure)
          ALLOCATE (rho_atom_set(iatom)%rho_rad_s(nspins),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,failure)
          ALLOCATE (rho_atom_set(iatom)%vrho_rad_h(nspins),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,failure)
          ALLOCATE (rho_atom_set(iatom)%vrho_rad_s(nspins),STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,failure)

          istat = 0
          DO ispin = 1,nspins

             NULLIFY(rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef, &
                     rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef, &
                     rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef, &
                     rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef)
          END DO

          ALLOCATE (rho_atom_set(iatom)%cpc_h(nspins), &
                    rho_atom_set(iatom)%cpc_s(nspins), &
                    rho_atom_set(iatom)%drho_rad_h(nspins), &
                    rho_atom_set(iatom)%drho_rad_s(nspins), &
                    rho_atom_set(iatom)%trho_rad_h(3,nspins), &
                    rho_atom_set(iatom)%trho_rad_s(3,nspins), &
                    rho_atom_set(iatom)%rho_rad_h_d(3,nspins), &
                    rho_atom_set(iatom)%rho_rad_s_d(3,nspins), &
                    STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,failure)

          IF (paw_atom) THEN
             DO ispin = 1,nspins
                NULLIFY(rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef, &
                        rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef, &
                        rho_atom_set(iatom)%trho_rad_h(1,ispin)%r_coef, &
                        rho_atom_set(iatom)%trho_rad_h(2,ispin)%r_coef, &
                        rho_atom_set(iatom)%trho_rad_h(3,ispin)%r_coef, &
                        rho_atom_set(iatom)%trho_rad_s(1,ispin)%r_coef, &
                        rho_atom_set(iatom)%trho_rad_s(2,ispin)%r_coef, &
                        rho_atom_set(iatom)%trho_rad_s(3,ispin)%r_coef)
                ALLOCATE(rho_atom_set(iatom)%cpc_h(ispin)%r_coef(1:nsatbas,1:nsatbas), &
                         rho_atom_set(iatom)%cpc_s(ispin)%r_coef(1:nsatbas,1:nsatbas),&
                         STAT=istat)
                CPPostcondition(istat==0,cp_failure_level,routineP,failure)

                rho_atom_set(iatom)%cpc_h(ispin)%r_coef = 0.0_dp
                rho_atom_set(iatom)%cpc_s(ispin)%r_coef = 0.0_dp
             END DO

          ELSE
             DO ispin = 1,nspins
                NULLIFY(rho_atom_set(iatom)%cpc_h(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%cpc_s(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef)

                DO j = 1,3
                  NULLIFY(rho_atom_set(iatom)%trho_rad_h(j,ispin)%r_coef)
                  NULLIFY(rho_atom_set(iatom)%trho_rad_s(j,ispin)%r_coef)
                  NULLIFY(rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef)
                  NULLIFY(rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef)
                END DO
             END DO
          ENDIF

       END DO ! iat

       num_pe = para_env%num_pe
       mepos = para_env%mepos
       bo = get_limit( nat, num_pe, mepos )
       DO iat = bo(1), bo(2)
          iatom = atom_list(iat)
          ALLOCATE (rho_atom_set(iatom)%ga_Vlocal_gb_h(nspins), &
                    rho_atom_set(iatom)%ga_Vlocal_gb_s(nspins), &
                    STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,failure)
          IF (paw_atom) THEN
             DO ispin = 1,nspins
                NULLIFY(rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef)

                CALL reallocate(rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef,&
                                1,nsotot,1,nsotot)
                CALL reallocate(rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef,&
                                1,nsotot,1,nsotot)

                rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef = 0.0_dp
                rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef = 0.0_dp
             END DO
          ELSE
             DO ispin = 1,nspins
                NULLIFY(rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef)
             END DO
          ENDIF

       END DO ! iat

    END DO

    CALL timestop(handle)

  END SUBROUTINE allocate_rho_atom_internals

! *****************************************************************************
!> \brief ...
!> \param rho_atom_set ...
!> \param iatom ...
!> \param ispin ...
!> \param nr ...
!> \param max_iso_not0 ...
! *****************************************************************************
  SUBROUTINE allocate_rho_atom_rad(rho_atom_set,iatom,ispin,nr,max_iso_not0)

    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    INTEGER, INTENT(IN)                      :: iatom, ispin, nr, max_iso_not0

    CHARACTER(len=*), PARAMETER :: routineN = 'allocate_rho_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, j
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure = .FALSE.

         ALLOCATE(rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef(1:nr,1:max_iso_not0), &
                  rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef(1:nr,1:max_iso_not0), &
                  rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef(1:nr,1:max_iso_not0), &
                  rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef(1:nr,1:max_iso_not0), &
                  STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,failure)

         rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef=0.0_dp
         rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef=0.0_dp
         rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef=0.0_dp
         rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef=0.0_dp

         ALLOCATE(rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef(nr,max_iso_not0),&
                  rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef(nr,max_iso_not0),&
                  STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,failure)
         rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef = 0.0_dp
         rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef = 0.0_dp

         DO j = 1,3
           ALLOCATE(rho_atom_set(iatom)%trho_rad_h(j,ispin)%r_coef(nr,max_iso_not0),&
                    rho_atom_set(iatom)%trho_rad_s(j,ispin)%r_coef(nr,max_iso_not0),&
                    STAT=istat)
           CPPrecondition(istat==0,cp_failure_level,routineP,failure)
           rho_atom_set(iatom)%trho_rad_h(j,ispin)%r_coef = 0.0_dp
           rho_atom_set(iatom)%trho_rad_s(j,ispin)%r_coef = 0.0_dp
         END DO

         DO j = 1,3
            NULLIFY(rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef, &
                    rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef)
            ALLOCATE(rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef(nr,max_iso_not0), &
                     rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef(nr,max_iso_not0), &
                     STAT=istat)
            CPPrecondition(istat==0,cp_failure_level,routineP,failure)
            rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef = 0.0_dp
            rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef = 0.0_dp
         END DO

    CALL timestop(handle)

  END SUBROUTINE  allocate_rho_atom_rad

! *****************************************************************************
!> \brief ...
!> \param rho_atom_set ...
!> \param iatom ...
!> \param ispin ...
! *****************************************************************************
  SUBROUTINE set2zero_rho_atom_rad(rho_atom_set,iatom,ispin)

!    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    INTEGER, INTENT(IN)                      :: iatom, ispin

    CHARACTER(len=*), PARAMETER :: routineN = 'set2zero_rho_atom_rad', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j
    LOGICAL                                  :: failure

    failure = .FALSE.

      rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef = 0.0_dp

      rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef = 0.0_dp

      rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef = 0.0_dp
      rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef = 0.0_dp

      DO j = 1,3
        rho_atom_set(iatom)%trho_rad_h(j,ispin)%r_coef = 0.0_dp
        rho_atom_set(iatom)%trho_rad_s(j,ispin)%r_coef = 0.0_dp
        rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef = 0.0_dp
        rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef = 0.0_dp
      END DO

  END SUBROUTINE set2zero_rho_atom_rad

END MODULE qs_rho_atom_methods
