
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
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE erf_fn,                          ONLY: erf
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
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
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_grid_atom,                    ONLY: create_grid_atom,&
                                             grid_atom_type
  USE qs_harmonics_atom,               ONLY: create_harmonics_atom,&
                                             harmonics_atom_type
  USE qs_neighbor_list_types,          ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_type, &
       neighbor_node_type, next, reduced_3c_list_type
  USE qs_oce_types,                    ONLY: oce_matrix_type
  USE qs_rho_atom_types,               ONLY: allocate_rho_atom_set,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: symmetric_add
  USE sparse_matrix_types,             ONLY: get_block_node,&
                                             real_matrix_p_type
  USE spherical_harmonics,             ONLY: clebsch_gordon,&
                                             clebsch_gordon_deallocate,&
                                             clebsch_gordon_init
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
!,locate
  USE whittaker,                       ONLY: whittaker_c0,&
                                             whittaker_ci
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters (only in this module)

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_rho_atom_methods'

! *** Public subroutines ***

  PUBLIC :: allocate_rho_atom_internals, &
            calculate_rho_atom, &
            calculate_rho_atom_coeff, &
            init_rho_atom

! *****************************************************************************

 CONTAINS

! *****************************************************************************
  SUBROUTINE calculate_rho_atom(para_env,rho_atom_set,atom_kind,atom_list,&
                                grid_atom,natom,nspins,tot_rho1_h,tot_rho1_s,error)

    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    INTEGER, DIMENSION(:), INTENT(IN)        :: atom_list
    TYPE(grid_atom_type), INTENT(IN)         :: grid_atom
    INTEGER, INTENT(IN)                      :: natom, nspins
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: tot_rho1_h, tot_rho1_s
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, i1, i2, iat, iatom, icg, ipgf1, ipgf2, ir, iset1, &
      iset2, iso, iso1, iso1_first, iso1_last, iso2, iso2_first, iso2_last, &
      istat, j, l, l_iso, l_sub, l_sum, lmax12, lmax_expansion, lmin12, m1s, &
      m2s, max_iso_not0, maxcg, maxco, maxl, maxso, mepos, n1s, n2s, nr, &
      nset, num_pe, size1, size2
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf
    LOGICAL                                  :: failure
    LOGICAL, ALLOCATABLE, DIMENSION(:, :)    :: done_vgg
    REAL(dp)                                 :: c1, c2, rho_h, rho_s, &
                                                root_zet12, zet12
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: erf_zet12, g1, g2, gg0, int1, &
                                                int2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: CPCH_sphere, CPCS_sphere, &
                                                dgg, gg
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: vgg
    REAL(dp), DIMENSION(:, :), POINTER       :: coeff, zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: my_CG_dxyz
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics

! *****************************************************************************

    CALL timeset("calculate_rho_atom","I"," ",handle)

    failure = .FALSE.
    NULLIFY(orb_basis)
    NULLIFY(harmonics)
    NULLIFY(lmin,lmax,npgf,zet,my_CG,my_CG_dxyz,coeff)

    CALL get_atomic_kind(atomic_kind=atom_kind, orb_basis_set=orb_basis,&
                         harmonics=harmonics)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxco=maxco,maxl=maxl,npgf=npgf,nset=nset,zet=zet,&
                          maxso=maxso)

    max_iso_not0 = harmonics%max_iso_not0
    maxcg        = harmonics%maxcg
    nr           = grid_atom%nr
    lmax_expansion = indso(1,max_iso_not0)
!   Distribute the atoms of this kind
    num_pe = para_env%num_pe
    mepos  = para_env%mepos
    bo = get_limit( natom, num_pe, mepos )

    my_CG        => harmonics%my_CG
    my_CG_dxyz   => harmonics%my_CG_dxyz

    ALLOCATE(CPCH_sphere(nsoset(maxl),nsoset(maxl)),STAT=istat)
    ALLOCATE(CPCS_sphere(nsoset(maxl),nsoset(maxl)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "CPCH_sphere,CPCS_sphere",2*nsoset(maxl)*nsoset(maxl)*dp_size)
    ALLOCATE(g1(nr),g2(nr),gg0(nr),gg(nr,0:2*maxl),dgg(nr,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "g1,g2,gg0,gg,dgg",3*nr+2*nr*(2*maxl+1)*dp_size)
    ALLOCATE(erf_zet12(nr),vgg(nr,0:2*maxl,0:indso(1,max_iso_not0)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "erf_zet12,vgg",nr+nr*(2*maxl+1)*(2*maxl+1)*dp_size)
    ALLOCATE(done_vgg(0:2*maxl,0:indso(1,max_iso_not0)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "done_vgg",(2*maxl+1)*(2*maxl+1)*int_size)
    ALLOCATE(int1(nr),int2(nr),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "int1,int2",2*nr*dp_size)

!******************************************************************************

    m1s = 0

    DO iset1 = 1,nset
      n1s = nsoset(lmax(iset1))
      DO ipgf1  = 1,npgf(iset1)

        iso1_first = nsoset(lmin(iset1)-1)+1+n1s*(ipgf1-1)+m1s
        iso1_last  = nsoset(lmax(iset1))+n1s*(ipgf1-1)+m1s
        size1 = iso1_last - iso1_first + 1
        i1 = nsoset(lmin(iset1)-1)+1

        g1(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))
        m2s = 0
        DO iset2 = 1,nset
          n2s=nsoset(lmax(iset2))
          DO ipgf2 = 1,npgf(iset2)
            iso2_first = nsoset(lmin(iset2)-1)+1+n2s*(ipgf2-1)+m2s
            iso2_last  = nsoset(lmax(iset2))+n2s*(ipgf2-1)+m2s
            size2 = iso2_last - iso2_first + 1
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
            vgg = 0.0_dp
            done_vgg = .FALSE.
!test reduce expansion local densities
       IF(lmin12 .LE. lmax_expansion) THEN
            IF (lmin12 == 0) THEN
              gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
              gg0(1:nr) = gg(1:nr,lmin12)
            ELSE
              gg0(1:nr) = g1(1:nr)*g2(1:nr)
              gg(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12)*g1(1:nr)*g2(1:nr)
            END IF
           
!test reduce expansion local densities
            IF(lmax12 .GT. lmax_expansion) lmax12 = lmax_expansion

            DO l=lmin12+1,lmax12
              gg(1:nr,l) = grid_atom%rad(1:nr)*gg(1:nr,l-1)
              dgg(1:nr,l-1) = dgg(1:nr,l-1) -2.0_dp*(zet(ipgf1,iset1)+ &
                             zet(ipgf2,iset2))*gg(1:nr,l)

            END DO
            dgg(1:nr,lmax12) = dgg(1:nr,lmax12) -2.0_dp*(zet(ipgf1,iset1)+ &
                             zet(ipgf2,iset2))*grid_atom%rad(1:nr)*&
                             gg(1:nr,lmax12)

            c2 = SQRT(pi*pi*pi/(zet12*zet12*zet12))
            DO iso = 1,max_iso_not0
               l_iso = indso(1,iso)
               c1 = fourpi/(2._dp*REAL(l_iso,dp)+1._dp)

               DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                 iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                 iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)

                 l = indso(1,iso1) + indso(1,iso2)
!test reduce expansion local densities
                    CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)
               IF(.NOT. failure) THEN

                 IF(done_vgg(l,l_iso)) CYCLE
                 L_sum = l + l_iso
                 L_sub = l - l_iso

                 IF(l_sum == 0) THEN
                  vgg(1:nr,l,l_iso) = erf_zet12(1:nr)* grid_atom%oorad2l(1:nr,1) * c2
                 ELSE
                   CALL whittaker_c0(int1,grid_atom%rad,gg0,erf_zet12,zet12,L_sum,nr)
                   CALL whittaker_ci(int2,grid_atom%rad,gg0,zet12,L_sub,nr)

                   DO ir = 1,nr
                     int1(ir) = grid_atom%oorad2l(ir,l_iso+1)*int1(ir)
                     int2(ir) = grid_atom%rad2l(ir,l_iso)*int2(ir)
                     vgg(ir,l,l_iso) = c1 * (int1(ir) + int2(ir))
                   END DO
                 END IF 
                 done_vgg(l,l_iso) = .TRUE.
!test reduce expansion local densities
               END IF ! failure

               END DO
            END DO
!test reduce expansion local densities
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

                DO iso = 1,max_iso_not0
                  l_iso = indso(1,iso)
                   DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                    iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                    iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)

                    l = indso(1,iso1) + indso(1,iso2)
!test reduce expansion local densities
                    CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)
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
!test reduce expansion local densities
                    END IF !failure

                  ENDDO  ! icg
                  DO icg = 1,harmonics%ncg_dxyz(iso,iset1,iset2)
                    iso1 = harmonics%ind_cg_a_dxyz(icg,iso,iset1,iset2)
                    iso2 = harmonics%ind_cg_b_dxyz(icg,iso,iset1,iset2)

                    l = indso(1,iso1) + indso(1,iso2)
!test reduce expansion local densities
                    CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)
                    IF(.NOT. failure) THEN
                    DO j = 1,3
                      rho_atom_set(iatom)%rho_rad_h_d(j,i)%r_coef(1:nr,iso) =&
                         rho_atom_set(iatom)%rho_rad_h_d(j,i)%r_coef(1:nr,iso) +&
                         gg(1:nr,l)*CPCH_sphere(iso1,iso2)*my_CG_dxyz(j,iso1,iso2,iso)

                      rho_atom_set(iatom)%rho_rad_s_d(j,i)%r_coef(1:nr,iso) =&
                         rho_atom_set(iatom)%rho_rad_s_d(j,i)%r_coef(1:nr,iso) +&
                         gg(1:nr,l)*CPCS_sphere(iso1,iso2)*my_CG_dxyz(j,iso1,iso2,iso)
                    END DO
!test reduce expansion local densities
                    END IF !failure
                  END DO  ! icg

                ENDDO  ! iso 

             ENDDO ! i
            ENDDO  ! iat
        
          ENDDO ! ipgf2
          m2s = m2s+maxso

        ENDDO  ! iset2
      ENDDO   ! ipgf1
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
                     rho_atom_set(iatom)%rho_rad_h(i)%r_coef(ir,iso)*&
                     grid_atom%wr(ir)
                rho_s = rho_s + &
                     rho_atom_set(iatom)%rho_rad_s(i)%r_coef(ir,iso)*&
                     grid_atom%wr(ir)
           !dbg
    !         write(17+iso+100*iat,'(2(I0,1x),3f14.6)') ir , iso, grid_atom%rad(ir), &
    !             rho_atom_set(iatom)%rho_rad_h(i)%r_coef(ir,iso), &
    !              rho_atom_set(iatom)%rho_rad_s(i)%r_coef(ir,iso)
           !dbg
             END DO  ! ir
             tot_rho1_h(i) = tot_rho1_h(i) + rho_h*harmonics%slm_int(iso)
             tot_rho1_s(i) = tot_rho1_s(i) + rho_s*harmonics%slm_int(iso)
          END DO  ! iso

       END DO  ! ispin

    END DO  ! iat

    DEALLOCATE(CPCH_sphere,CPCS_sphere,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
            "CPCH_sphere,CPCS_sphere")
    ENDIF
    DEALLOCATE(g1,g2,gg0,gg,dgg,vgg,done_vgg,erf_zet12,int1,int2,STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,&
            "g1,g2,gg0,gg,dgg,vgg,done_vgg,erf_zet12,int1,int2")
    ENDIF

   CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_rho_atom

! *****************************************************************************

  SUBROUTINE calculate_rho_atom_coeff(qs_env, rho_ao, rho_atom_set)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: rho_ao
    TYPE(rho_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: rho_atom_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_rho_atom_coeff'

    INTEGER :: handle, i, iab, iatom, icol, ikind, ilist, inode, irow, ispin, &
      istat, j, jatom, jkind, katom, kkind, kneighbor, last_jatom, &
      last_katom, max_gau, max_nsgf, maxsoc, nat_kind, natom, nkind, nlist, &
      nnode, nnode_c, nsetc, nsgfa, nsgfb, nsoctot, nspins
    INTEGER, DIMENSION(:), POINTER           :: a_list
    LOGICAL                                  :: new_atom_b, new_atom_c, &
                                                paw_atom
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: PC
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: coef_aux_h, coef_aux_s
    REAL(KIND=dp), DIMENSION(3)              :: rab
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: C_coeff_hh_a, C_coeff_hh_b, &
                                                C_coeff_ss_a, C_coeff_ss_b
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_orb
    TYPE(neighbor_list_type), POINTER        :: sab_orb_neighbor_list, &
                                                sab_orb_neighbor_list_local
    TYPE(neighbor_node_type), POINTER        :: sab_orb_neighbor_node
    TYPE(oce_matrix_type), POINTER           :: oce
    TYPE(reduced_3c_list_type), &
      DIMENSION(:), POINTER                  :: reduced_3c_oce
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: p_block_spin
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom

!
! *****************************************************************************

    CALL timeset("calculate_rho_atom_coeff","I"," ",handle)

    NULLIFY(dft_control)
    NULLIFY(atomic_kind_set)
    NULLIFY(oce)
    NULLIFY(sab_orb)
    NULLIFY(rho_atom)
    NULLIFY(orb_basis_set)
    NULLIFY(reduced_3c_oce)

    CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            atomic_kind_set=atomic_kind_set,&
            oce=oce,sab_orb=sab_orb,&
            rho_atom_set=rho_atom)

    IF (PRESENT(rho_atom_set)) rho_atom => rho_atom_set

    nspins=dft_control%nspins

    ALLOCATE(p_block_spin(nspins),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                         "p_block_spin")
    DO ispin = 1,nspins
      NULLIFY(p_block_spin(ispin)%r_coef)
    END DO

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxsgf=max_nsgf,&
                             maxgtops=max_gau,natom=natom)

    ALLOCATE (PC(max_nsgf,max_gau),STAT=istat)
    ALLOCATE (coef_aux_h(max_gau,max_gau,nspins),STAT=istat)
    ALLOCATE (coef_aux_s(max_gau,max_gau,nspins),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "PC,coef_aux_h,coef_aux_s",3*nspins*max_gau*max_gau*dp_size)

    nkind = SIZE(atomic_kind_set)
!   Inizialize to 0 the CPC coefficients and the local density arrays
    DO ikind = 1 ,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=a_list,natom=nat_kind,&
                           paw_atom=paw_atom)

      IF(.NOT. paw_atom) CYCLE
      DO i = 1,nat_kind
        iatom = a_list(i)
        DO ispin = 1, nspins
          rho_atom(iatom)%cpc_h(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%cpc_s(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%rho_rad_h(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%rho_rad_s(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%drho_rad_h(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%drho_rad_s(ispin)%r_coef = 0.0_dp
          DO j = 1,3
             rho_atom(iatom)%rho_rad_h_d(j,ispin)%r_coef = 0.0_dp
             rho_atom(iatom)%rho_rad_s_d(j,ispin)%r_coef = 0.0_dp 
          END DO
          rho_atom(iatom)%vrho_rad_h(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%vrho_rad_s(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%ga_Vlocal_gb_h(ispin)%r_coef = 0.0_dp
          rho_atom(iatom)%ga_Vlocal_gb_s(ispin)%r_coef = 0.0_dp
        ENDDO  ! ispin
      ENDDO  ! i
    END DO ! ikind

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
          last_jatom = 0
          sab_orb_neighbor_node => first_node(sab_orb_neighbor_list_local)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=sab_orb_neighbor_node,&
                                   neighbor=jatom,&
                                   r=rab,&
                                   reduced_3c_oce=reduced_3c_oce)

            IF (jatom /= last_jatom) THEN
              new_atom_b = .TRUE.
              last_jatom = jatom
            ELSE
              new_atom_b = .FALSE.
              sab_orb_neighbor_node => next(sab_orb_neighbor_node)
              CYCLE
            END IF

            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom

               DO ispin = 1,nspins
                  NULLIFY(p_block_spin(ispin)%r_coef)
                  CALL get_block_node(matrix=rho_ao(ispin)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=p_block_spin(ispin)%r_coef)
               ENDDO

            ELSE
               irow = jatom
               icol = iatom
               DO ispin = 1,nspins
                  NULLIFY(p_block_spin(ispin)%r_coef)
                  CALL get_block_node(matrix=rho_ao(ispin)%matrix,&
                                    block_row=irow,&
                                    block_col=icol,&
                                    BLOCK=p_block_spin(ispin)%r_coef)

               ENDDO
            ENDIF
            
            DO kkind=1,nkind

              atomic_kind => atomic_kind_set(kkind)

              CALL get_atomic_kind(atomic_kind=atomic_kind,&
                                   orb_basis_set=orb_basis_set,&
                                   paw_atom=paw_atom)

              CALL  get_gto_basis_set(gto_basis_set=orb_basis_set,&
                                      nset =nsetc,maxso=maxsoc)

              nsoctot = maxsoc*nsetc

              IF(.NOT. paw_atom)  CYCLE

              last_katom = 0

              nnode_c = reduced_3c_oce(kkind)%nnode
              DO kneighbor = 1,nnode_c

                katom = reduced_3c_oce(kkind)%index_atom(kneighbor)

                IF (katom /= last_katom) THEN
                   new_atom_c = .TRUE.
                   last_katom = katom
                ELSE
                   new_atom_c = .FALSE.
                   CYCLE
                END IF

                NULLIFY(C_coeff_hh_a)
                CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=iatom,&
                                    BLOCK=C_coeff_hh_a)

                NULLIFY(C_coeff_ss_a)
                CALL get_block_node(matrix=oce%softo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=iatom,&
                                    BLOCK=C_coeff_ss_a)

                IF(iatom == jatom) THEN

                   DO ispin = 1,nspins
                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                                p_block_spin(ispin)%r_coef,&
                                SIZE(p_block_spin(ispin)%r_coef,1),&
                                C_coeff_hh_a,SIZE(C_coeff_hh_a,1),&
                                0.0_dp,PC,max_nsgf)

                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_hh_a,&
                                SIZE(C_coeff_hh_a,1),PC,max_nsgf,1.0_dp,&
                                rho_atom(katom)%cpc_h(ispin)%r_coef,nsoctot)

                     CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                                p_block_spin(ispin)%r_coef,&
                                SIZE(p_block_spin(ispin)%r_coef,1),&
                                C_coeff_ss_a,SIZE(C_coeff_ss_a,1),&
                                0.0_dp,PC,max_nsgf)

                     CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,C_coeff_ss_a,&
                                SIZE(C_coeff_ss_a,1),PC,max_nsgf,1.0_dp,&
                                rho_atom(katom)%cpc_s(ispin)%r_coef,nsoctot)
                   END DO

                ELSEIF(iatom < jatom) THEN
                   NULLIFY(C_coeff_hh_b)
                   CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=jatom,&
                                    BLOCK=C_coeff_hh_b)

                   NULLIFY(C_coeff_ss_b)
                   CALL get_block_node(matrix=oce%softo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=jatom,&
                                    BLOCK=C_coeff_ss_b)

                   DO ispin = 1,nspins
                      CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                                 p_block_spin(ispin)%r_coef,&
                                 SIZE(p_block_spin(ispin)%r_coef,1),&
                                 C_coeff_hh_b,SIZE(C_coeff_hh_b,1),&
                                 0.0_dp,PC,max_nsgf)


                      CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,&
                                 1.0_dp,C_coeff_hh_a,SIZE(C_coeff_hh_a,1),&
                                 PC,max_nsgf,0.0_dp,&
                                 coef_aux_h(1,1,ispin),max_gau)

                      CALL dgemm('N','T',nsgfa,nsoctot,nsgfb,1.0_dp,&
                                 p_block_spin(ispin)%r_coef,&
                                 SIZE(p_block_spin(ispin)%r_coef,1),&
                                 C_coeff_ss_b,SIZE(C_coeff_ss_b,1),&
                                 0.0_dp,PC,max_nsgf)

                      CALL dgemm('N','N',nsoctot,nsoctot,nsgfa,1.0_dp,&
                                 C_coeff_ss_a,SIZE(C_coeff_ss_a,1),PC,&
                                 max_nsgf,0.0_dp,&
                                 coef_aux_s(1,1,ispin),max_gau)

                      CALL symmetric_add(rho_atom(katom)%cpc_h(ispin)%r_coef,nsoctot,&
                                       coef_aux_h(1,1,ispin),max_gau,nsoctot)
                      CALL symmetric_add(rho_atom(katom)%cpc_s(ispin)%r_coef,nsoctot,&
                                       coef_aux_s(1,1,ispin),max_gau,nsoctot)
                   END DO

                ELSEIF(iatom > jatom) THEN
                   NULLIFY(C_coeff_hh_b)
                   CALL get_block_node(matrix=oce%hardo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=jatom,&
                                    BLOCK=C_coeff_hh_b)

                   NULLIFY(C_coeff_ss_b)
                   CALL get_block_node(matrix=oce%softo(1)%matrix,&
                                    block_row=katom,&
                                    block_col=jatom,&
                                    BLOCK=C_coeff_ss_b)

                   DO ispin = 1,nspins
                      CALL dgemm('N','T',nsgfb,nsoctot,nsgfa,1.0_dp,&
                                 p_block_spin(ispin)%r_coef,&
                                 SIZE(p_block_spin(ispin)%r_coef,1),&
                                 C_coeff_hh_a,SIZE(C_coeff_hh_a,1),&
                                 0.0_dp,PC,max_nsgf)


                      CALL dgemm('N','N',nsoctot,nsoctot,nsgfb,&
                                 1.0_dp,C_coeff_hh_b,SIZE(C_coeff_hh_b,1),&
                                 PC,max_nsgf,0.0_dp,&
                                 coef_aux_h(1,1,ispin),max_gau)

                      CALL dgemm('N','T',nsgfb,nsoctot,nsgfa,1.0_dp,&
                                 p_block_spin(ispin)%r_coef,&
                                 SIZE(p_block_spin(ispin)%r_coef,1),&
                                 C_coeff_ss_a,SIZE(C_coeff_ss_a,1),&
                                 0.0_dp,PC,max_nsgf)

                      CALL dgemm('N','N',nsoctot,nsoctot,nsgfb,1.0_dp,&
                                 C_coeff_ss_b,SIZE(C_coeff_ss_b,1),PC,&
                                 max_nsgf,0.0_dp,&
                                 coef_aux_s(1,1,ispin),max_gau)

                      CALL symmetric_add(rho_atom(katom)%cpc_h(ispin)%r_coef,nsoctot,&
                                       coef_aux_h(1,1,ispin),max_gau,nsoctot)
                      CALL symmetric_add(rho_atom(katom)%cpc_s(ispin)%r_coef,nsoctot,&
                                       coef_aux_s(1,1,ispin),max_gau,nsoctot)
                   END DO

                ENDIF
              ENDDO  ! kneighbor
            ENDDO  ! kkind
          sab_orb_neighbor_node => next(sab_orb_neighbor_node)
          ENDDO  ! inode
        ENDDO  ! ilist
      ENDDO  ! jkind
    ENDDO ! ikind

    DO iatom = 1,natom

       DO ispin = 1,nspins
         IF(ASSOCIATED(rho_atom(iatom)%cpc_h(ispin)%r_coef)) THEN
           CALL mp_sum(rho_atom(iatom)%cpc_h(ispin)%r_coef,qs_env%para_env%group)
           CALL mp_sum(rho_atom(iatom)%cpc_s(ispin)%r_coef,qs_env%para_env%group)
         END IF
       ENDDO

    END DO

    DEALLOCATE (PC,coef_aux_h,coef_aux_s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
        "PC,coef_aux_h,coef_aux_s")
    DEALLOCATE (p_block_spin, STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
        "p_block_spin")

   CALL timestop(0.0_dp,handle)

  END SUBROUTINE calculate_rho_atom_coeff

! *****************************************************************************

  SUBROUTINE init_rho_atom(qs_env,gapw_control)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(gapw_control_type), POINTER         :: gapw_control

    CHARACTER(LEN=*), PARAMETER              :: routineN = 'init_rho_atom'

    INTEGER :: ia, ikind, il, iso, iso1, iso2, istat, l1, l1l2, l2, lc1, lc2, &
      lcleb, ll, llmax, lmax_sphere, lp, m1, m2, max_s_harm, max_s_set, &
      maxlgto, mm, mp, na, nat, natom, nkind, nr, nspins, quadrature
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: paw_atom
    REAL(dp)                                 :: cosia, phi
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: rga
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set

! -----------------------------------------------------------------------------

    NULLIFY(rho_atom_set)
    NULLIFY(atomic_kind_set)
    NULLIFY(atomic_kind)
    NULLIFY(dft_control)
    NULLIFY(orb_basis_set)
    NULLIFY(my_CG, atomic_kind, grid_atom, harmonics,atom_list)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                    dft_control=dft_control) 
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxlgto=maxlgto,&
                             natom = natom)

    nspins = dft_control%nspins

    nkind = SIZE(atomic_kind_set)
    lmax_sphere = gapw_control%lmax_sphere

    llmax = MIN(lmax_sphere, 2*maxlgto)
    max_s_harm = nsoset(llmax)
    max_s_set = nsoset(maxlgto)

    lcleb = MAX(llmax,2*maxlgto,1)

!   *** allocate calculate the CG coefficients up to the maxl *** 
    CALL clebsch_gordon_init ( lcleb )
    CALL reallocate(my_CG,1,max_s_set,1,max_s_set,1,max_s_harm)

    ALLOCATE(rga(lcleb,2),STAT=istat)
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
    CALL clebsch_gordon_deallocate 

!   *** initialize the Lebedev grids ***
    CALL init_lebedev_grids
    quadrature=gapw_control%quadrature

    DO ikind = 1,nkind
      atomic_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           orb_basis_set=orb_basis_set, &
                           atom_list=atom_list,natom=nat, &
                           paw_atom=paw_atom,&
                           grid_atom=grid_atom,&
                           harmonics=harmonics,&
                           ngrid_rad=nr,ngrid_ang=na)

!     *** determine the Lebedev grid for this kind ***

      ll =  get_number_of_lebedev_grid(n=na)
      na = lebedev_grid(ll)%n
      grid_atom%ng_sphere = na
      grid_atom%nr = nr
  
      IF(llmax>ll) THEN
        WRITE(6,'(/,72("*"))')
        WRITE(6,'(T2,A,T66,I4)') & 
          "WARNING: the lebedev grid is built for angular momentum l up to ",ll,&
          "            the max l of spherical harmonics is larger, l_max = ", llmax,&
          "                   good integration is guaranteed only for l <= ", ll
        WRITE(6,'(72("*"),/)')
      END IF

!     *** calculate the radial grid ***
      CALL create_grid_atom(grid_atom,nr,na,llmax,quadrature)

      grid_atom%wa(1:na) = 4._dp*pi*lebedev_grid(ll)%w(1:na)

      DO ia = 1,na
        cosia = lebedev_grid(ll)%r(3,ia)
        grid_atom%cos_phi(ia)   =   cosia
        grid_atom%theta(ia)     =  &
               ATAN2( lebedev_grid(ll)%r(2,ia) , lebedev_grid(ll)%r(1,ia) )
        grid_atom%cos_theta(ia) = COS( grid_atom%theta(ia) )
        phi =  ACOS(cosia)
        IF( grid_atom%sin_phi(ia) < 0.0_dp ) phi = - phi
        grid_atom%phi(ia) = phi
        grid_atom%sin_phi(ia)  = SIN(grid_atom%phi(ia)) 
        IF(ABS(grid_atom%sin_phi(ia)) > EPSILON(1.0_dp)) THEN
          grid_atom%cotan_phi(ia) = &
                grid_atom%cos_phi(ia)/grid_atom%sin_phi(ia)
        ELSE
          grid_atom%cotan_phi(ia) = 0.0_dp
        END IF

        grid_atom%sin_theta(ia) = SIN( grid_atom%theta(ia) ) 
        IF(ABS(grid_atom%sin_theta(ia)) > EPSILON(1.0_dp)) THEN
          grid_atom%usin_theta(ia)  = 1.0_dp/grid_atom%sin_theta(ia)
        ELSE
          grid_atom%usin_theta(ia) = 1.0_dp
        END IF

      ENDDO

!     *** calculate the spherical harmonics on the grid ***

      CALL create_harmonics_atom(harmonics, orb_basis_set,&
                   my_CG,na,llmax,max_s_harm,ll,grid_atom%wa,&
                   grid_atom%phi,grid_atom%cos_phi,&
                   grid_atom%theta,grid_atom%cos_theta,grid_atom%cotan_phi)

    END DO

    CALL deallocate_lebedev_grids
    DEALLOCATE(my_CG)

    CALL allocate_rho_atom_internals(qs_env, rho_atom_set)

    CALL set_qs_env(qs_env=qs_env, rho_atom_set=rho_atom_set)

  END SUBROUTINE init_rho_atom

! *****************************************************************************

  SUBROUTINE allocate_rho_atom_internals(qs_env, rho_atom_set)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set

    CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_rho_atom_internals'

    INTEGER :: iat, iatom, ikind, ispin, istat, j, max_iso_not0, maxso, nat, &
      natom, nkind, nr, nset, nsotot, nspins, tmpstat
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: paw_atom
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    TYPE(harmonics_atom_type), POINTER       :: harmonics

! --------------------------------------------------------------------------

    NULLIFY(atomic_kind_set)
    NULLIFY(rho_atom_set)

    CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom = natom)

    nkind = SIZE(atomic_kind_set)
    nspins = qs_env%dft_control%nspins

    CALL allocate_rho_atom_set(rho_atom_set,natom)

    DO ikind = 1, nkind

       NULLIFY(atom_list, atomic_kind, harmonics)

       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            orb_basis_set=orb_basis_set, &
                            atom_list=atom_list,natom=nat, &
                            paw_atom=paw_atom,&
                            harmonics=harmonics,&
                            ngrid_rad=nr)

       CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                              maxso=maxso, nset=nset)
       nsotot = maxso * nset 
       max_iso_not0 = harmonics%max_iso_not0
       DO iat = 1,nat
          iatom = atom_list(iat) 
          !       *** allocate the radial density for each LM,for each atom ***

          ALLOCATE (rho_atom_set(iatom)%rho_rad_h(nspins),STAT=istat)
          ALLOCATE (rho_atom_set(iatom)%rho_rad_s(nspins),STAT=istat)
          ALLOCATE (rho_atom_set(iatom)%vrho_rad_h(nspins),STAT=istat)
          ALLOCATE (rho_atom_set(iatom)%vrho_rad_s(nspins),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
              "rho_rad_h,rho_rad_s,vrho_rad_h,vrho_rad_s",4*nspins)

          istat = 0
          DO ispin = 1,nspins

             NULLIFY(rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef, &
                     rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef, &
                     rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef, &
                     rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef)

             ALLOCATE(rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef(1:nr,1:max_iso_not0), &
                      rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef(1:nr,1:max_iso_not0), &
                      rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef(1:nr,1:max_iso_not0), &
                      rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef(1:nr,1:max_iso_not0), &
                      STAT=tmpstat)
             istat = istat + tmpstat

             rho_atom_set(iatom)%rho_rad_h(ispin)%r_coef=0.0_dp
             rho_atom_set(iatom)%rho_rad_s(ispin)%r_coef=0.0_dp
             rho_atom_set(iatom)%vrho_rad_h(ispin)%r_coef=0.0_dp
             rho_atom_set(iatom)%vrho_rad_s(ispin)%r_coef=0.0_dp

          END DO
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
              "r_coef of rho_rad_h,rho_rad_s,vrho_rad_h,vrho_rad_s",&
              4*nr*max_iso_not0*nspins*dp_size)

          ALLOCATE (rho_atom_set(iatom)%cpc_h(nspins), &
                    rho_atom_set(iatom)%cpc_s(nspins), &
                    rho_atom_set(iatom)%drho_rad_h(nspins), &
                    rho_atom_set(iatom)%drho_rad_s(nspins), &
                    rho_atom_set(iatom)%rho_rad_h_d(3,nspins), &
                    rho_atom_set(iatom)%rho_rad_s_d(3,nspins), &
                    STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                               "cpc_h,cpc_s,drho_rad_h,drho_rad_s,rho_rad_h_d,rho_rad_s_d")


          ALLOCATE (rho_atom_set(iatom)%ga_Vlocal_gb_h(nspins), &
                    rho_atom_set(iatom)%ga_Vlocal_gb_s(nspins), &
                    STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "ga_Vlocal_gb_h,ga_Vlocal_gb_s")

          IF (paw_atom) THEN 
             istat = 0
             DO ispin = 1,nspins
                NULLIFY(rho_atom_set(iatom)%cpc_h(ispin)%r_coef, &
                        rho_atom_set(iatom)%cpc_s(ispin)%r_coef, &
                        rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef, &
                        rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef)
                ALLOCATE(rho_atom_set(iatom)%cpc_h(ispin)%r_coef(1:nsotot,1:nsotot), &
                         rho_atom_set(iatom)%cpc_s(ispin)%r_coef(1:nsotot,1:nsotot), &
                         rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef(nr,max_iso_not0), &
                         rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef(nr,max_iso_not0), &
                         STAT=tmpstat)
                istat = istat + tmpstat

                rho_atom_set(iatom)%cpc_h(ispin)%r_coef = 0.0_dp
                rho_atom_set(iatom)%cpc_s(ispin)%r_coef = 0.0_dp
                rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef = 0.0_dp
                rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef = 0.0_dp

                DO j = 1,3
                   NULLIFY(rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef, &
                           rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef)
                   ALLOCATE(rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef(nr,max_iso_not0), &
                            rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef(nr,max_iso_not0), &
                            STAT=tmpstat)
                   istat = istat + tmpstat
                   rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef = 0.0_dp
                   rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef = 0.0_dp
                END DO

                NULLIFY(rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef)

                CALL reallocate(rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef,&
                                1,nsotot,1,nsotot)
                CALL reallocate(rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef,&
                                1,nsotot,1,nsotot)

                rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef = 0.0_dp
                rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef = 0.0_dp
             END DO

             IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                 " r_coef of cpc_h,cpc_s,drho_rad_h,drho_rad_s,rho_rad_h_d,"//&
                 "rho_rad_s_d",(5*nr*max_iso_not0+4*nsotot*nsotot)*nspins*dp_size)
          ELSE
             DO ispin = 1,nspins
                NULLIFY(rho_atom_set(iatom)%cpc_h(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%cpc_s(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%drho_rad_h(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%drho_rad_s(ispin)%r_coef)

                DO j = 1,3
                   NULLIFY(rho_atom_set(iatom)%rho_rad_h_d(j,ispin)%r_coef)
                   NULLIFY(rho_atom_set(iatom)%rho_rad_s_d(j,ispin)%r_coef)
                END DO
                NULLIFY(rho_atom_set(iatom)%ga_Vlocal_gb_h(ispin)%r_coef)
                NULLIFY(rho_atom_set(iatom)%ga_Vlocal_gb_s(ispin)%r_coef)
             END DO
          ENDIF

       END DO ! iat

    END DO

  END SUBROUTINE allocate_rho_atom_internals

END MODULE qs_rho_atom_methods
