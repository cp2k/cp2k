!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002 -2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

!!****h* cp2k/qs_vxc_atom [1.0] *
!!
!!   NAME
!!     qs_vxc_atom
!!
!!   FUNCTION
!!     routines that build the integrals of the Vxc potential calculated  
!!     for the atomic density in the basis set of spherical primitives 
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
MODULE qs_vxc_atom
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,nco,nso,&
                                             ncoset,nsoset
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: transform_s2c
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_atom,                         ONLY: vxc_of_r
  USE xc_functionals,                  ONLY: xc_calculate_lda,&
                                             xc_calculate_lsd,&
                                             xc_get_pot_size
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_release

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='qs_vxc_atom'

  PUBLIC :: calculate_vxc_atom


!******************************************************************************

  INTERFACE calculate_gaVxcgb
    MODULE PROCEDURE gaVxcgb_noGC,gaVxcgb_GC
  END INTERFACE

!******************************************************************************

CONTAINS

!******************************************************************************

 SUBROUTINE calculate_vxc_atom(qs_env,error)

!   ***************************************************************************

    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    TYPE(qs_environment_type), INTENT(INOUT) :: qs_env

    CHARACTER(LEN=*), PARAMETER         :: routine_name = "calculate_vxc_atom",&
                routineP = module_name//':'//routine_name

    TYPE(atomic_kind_type), DIMENSION(:),&
                            POINTER     :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER     :: atom_kind
    TYPE(dft_control_type), POINTER     :: dft_control
    TYPE(grid_atom_type), POINTER       :: grid_atom
    TYPE(harmonics_atom_type), POINTER  :: harmonics
    TYPE(pw_env_type), POINTER          :: pw_env
    TYPE(pw_pool_type), POINTER         :: pw_pool
    TYPE(qs_energy_type), POINTER       :: energy
    TYPE(rho_atom_type), DIMENSION(:), &
                              POINTER   :: rho_atom_set
    TYPE(rho_atom_type),  POINTER       :: rho_atom
    TYPE(rho_atom_coeff),DIMENSION(:),& 
                  POINTER               :: r_h,r_s,dr_h,dr_s
    TYPE(xc_derivative_set_type), &
                             POINTER    :: derivative_set

    REAL(dp), DIMENSION(:), POINTER     :: weight 
    REAL(dp), DIMENSION(:,:), &
                              POINTER   :: rho_h,rho_s
    REAL(dp), DIMENSION(:,:,:), &
                              POINTER   :: drho_h,drho_s,&
                                           vxc_h,vxc_s,&
                                           vxg_h,vxg_s
    LOGICAL                             :: cross_term,&
                                           gradient_f,lsd,paw_atom
    REAL(dp)                            :: exc_h,exc_s,rad,urad,usin_ia,rtot
    REAL(dp), ALLOCATABLE, &
                        DIMENSION(:,:)  :: locpot
    INTEGER                             :: i,ia,iat,iatom, ierr,&
                                           ikind,ipot,ir,&
                                           iso,ispin,lp,na,natom,nkind,&
                                           nr,npot,nspins,order 
    INTEGER, DIMENSION(:), POINTER      :: atom_list

!   ***************************************************************************

    NULLIFY(dft_control)
    NULLIFY(energy)
    NULLIFY(atomic_kind_set,atom_kind)
    NULLIFY(pw_env)
    NULLIFY(rho_atom_set,rho_atom)
    NULLIFY(pw_pool)
    NULLIFY(atom_kind)
    NULLIFY(grid_atom,harmonics)
    NULLIFY(atom_list)
    NULLIFY(derivative_set)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    energy=energy,&
                    atomic_kind_set=atomic_kind_set,&
                    pw_env=pw_env,&
                    rho_atom_set=rho_atom_set)

   nspins = dft_control%nspins
   gradient_f=ANY(dft_control%xc_control%gradient_functionals)
   cross_term = ANY (dft_control%xc_control%crossterms) .AND. nspins>1

   lsd = (nspins==2)
   order = 1
   npot = xc_get_pot_size(lsd,gradient_f , cross_term, order)
   
   CALL pw_env_get(pw_env,auxbas_pw_pool=pw_pool,error=error)

   CALL xc_dset_create(derivative_set, pw_pool=pw_pool, error=error)

   nkind = SIZE(atomic_kind_set)

!  Put to 0 the energy contribution from the one center XC terms 
   energy%exc1 = 0.0_dp

!  Here starts the loop over all the atoms
   DO ikind = 1,nkind

     atom_kind => atomic_kind_set(ikind)
     CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                        natom=natom,paw_atom=paw_atom,&
                        harmonics=harmonics,grid_atom=grid_atom)
     IF(.NOT. paw_atom) CYCLE

     NULLIFY(r_h,r_s,dr_h,dr_s)
     NULLIFY(rho_h,drho_h,rho_s,drho_s,weight)
     NULLIFY(vxc_h,vxc_s,vxg_h,vxg_s)

     nr = grid_atom%nr
     na = grid_atom%ng_sphere

     CALL reallocate(rho_h,1,na,1,nspins)
     CALL reallocate(rho_s,1,na,1,nspins)
     CALL reallocate(weight,1,na)
     CALL reallocate(vxc_h,1,na,1,nr,1,nspins)
     CALL reallocate(vxc_s,1,na,1,nr,1,nspins)
  
     ALLOCATE(derivative_set%data(na, 0:npot-1), stat=ierr)
!     CPPostconditionNoFail(ierr==0,cp_failure_level,routineP,error)

     IF (gradient_f) THEN
        ALLOCATE(locpot(na,npot), STAT=ierr)
!        CPPostconditionNoFail(ierr==0,cp_failure_level,routineP,error)
        CALL reallocate(drho_h,1,na,1,4,1,nspins)
        CALL reallocate(drho_s,1,na,1,4,1,nspins)
        CALL reallocate(vxg_h,1,na,1,nr,1,nspins)
        CALL reallocate(vxg_s,1,na,1,nr,1,nspins)
     END IF

     DO iat = 1,natom
       iatom = atom_list(iat)

       rho_atom => rho_atom_set(iatom)
       IF(gradient_f) THEN
         CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,&
                           rho_rad_s=r_s,drho_rad_h=dr_h,&
                           drho_rad_s=dr_s)
       ELSE
         CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,&
                           rho_rad_s=r_s)
       END IF

       rho_atom%exc_h = 0.0_dp
       rho_atom%exc_s = 0.0_dp
         
       rtot = 0.0_dp
       DO ir = 1,nr
         weight(1:na) = grid_atom%wr(ir)*grid_atom%wa(1:na)
         rad = grid_atom%rad(ir)
         urad = 1.0_dp/rad
         rho_h = 0.0_dp
         rho_s = 0.0_dp
!        calculate rho and drho on the points of the sphere
         IF(gradient_f) THEN
           DO ispin = 1,nspins
             DO iso = 1,harmonics%max_iso_not0
                DO ia = 1,na
                  usin_ia = 1.0_dp/grid_atom%sin_teta(ia) 
                  rho_h(ia,ispin) = rho_h(ia,ispin) + & 
                       r_h(ispin)%r_coef(ir,iso)*&
                       harmonics%slm(ia,iso)
                  rho_s(ia,ispin) = rho_s(ia,ispin) + & 
                       r_s(ispin)%r_coef(ir,iso)*&
                       harmonics%slm(ia,iso)

!                 *** Derivatives in spherical coordinates ***
!                 *** d/dr
                  drho_h(ia,1,ispin) = drho_h(ia,1,ispin) + &
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)
                  drho_s(ia,1,ispin) = drho_s(ia,1,ispin) + &
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

!                 *** d/d teta
                  drho_h(ia,2,ispin) = drho_h(ia,2,ispin) + &
                                       r_h(ispin)%r_coef(ir,iso)*&
                                       urad*&
                                       harmonics%dslm(1,ia,iso) 
                  drho_s(ia,2,ispin) = drho_s(ia,2,ispin) + &
                                       r_s(ispin)%r_coef(ir,iso)*&
                                       urad*&
                                       harmonics%dslm(1,ia,iso)

!                 *** d/d phi
                  drho_h(ia,3,ispin) = drho_h(ia,3,ispin) + &
                                       r_h(ispin)%r_coef(ir,iso)*&
                                       urad*usin_ia*&
                                       harmonics%dslm(2,ia,iso)
                  drho_s(ia,3,ispin) = drho_s(ia,3,ispin) + &
                                       r_s(ispin)%r_coef(ir,iso)*&
                                       urad*usin_ia*&
                                       harmonics%dslm(2,ia,iso)

                  drho_h(ia,4,ispin) = SQRT(&
                               drho_h(ia,1,ispin)*drho_h(ia,1,ispin)+& 
                               drho_h(ia,2,ispin)*drho_h(ia,2,ispin)+&
                               drho_h(ia,3,ispin)*drho_h(ia,3,ispin))

                  drho_s(ia,4,ispin) = SQRT(&
                               drho_s(ia,1,ispin)*drho_s(ia,1,ispin)+& 
                               drho_s(ia,2,ispin)*drho_s(ia,2,ispin)+&
                               drho_s(ia,3,ispin)*drho_s(ia,3,ispin))
                END DO
             END DO
           END DO
         ELSE         
           DO ispin = 1,nspins
             DO iso = 1,harmonics%max_iso_not0
                DO ia = 1,na
                  rho_h(ia,ispin) = rho_h(ia,ispin) + & 
                       r_h(ispin)%r_coef(ir,iso)*&
                       harmonics%slm(ia,iso)
                  rho_s(ia,ispin) = rho_s(ia,ispin) + & 
                       r_s(ispin)%r_coef(ir,iso)*&
                       harmonics%slm(ia,iso)
                END DO
             END DO
           END DO
!                DO ia = 1,na
!                  rtot = rtot + rho_h(ia,1)*weight(ia)
!                ENDDO
         END IF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        *** hard atom density ***
         CALL dcopy(na*npot, 0.0_dp, 0, derivative_set%data, 1)

         IF(gradient_f) THEN
           CALL vxc_of_r(dft_control%xc_control,derivative_set,nspins,na,rho_h,&
                         drho=drho_h,locpot=locpot,error=error)
         ELSE
           CALL vxc_of_r(dft_control%xc_control,derivative_set,nspins,na,rho_h)
         ENDIF

!        *** Sum up the exc contribution  to exc_h ***
         exc_h = SUM ( derivative_set%data(:,0)*weight(:) )
         rho_atom%exc_h = rho_atom%exc_h + exc_h

!        ***  vxc_h(ia,ir) 
         DO ispin = 1,nspins
           DO ia =1,na
             vxc_h(ia,ir,ispin) =  derivative_set%data(ia,ispin)*weight(ia)
           END DO
         END DO

!        *** if gradient_f :  vxg_h(ia,ir)
         IF (gradient_f) THEN
           DO ispin = 1,nspins
             ipot = nspins + ispin
             DO ia =1,na
               vxg_h(ia,ir,ispin) =  derivative_set%data(ia,ipot)*weight(ia)
             END DO
           END DO
         END IF

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        *** soft atom density ***
         CALL dcopy(na*npot, 0.0_dp, 0, derivative_set%data, 1)

         IF(gradient_f) THEN
           CALL vxc_of_r(dft_control%xc_control,derivative_set,nspins,na,rho_s,&
                         drho=drho_s,locpot=locpot,error=error)
         ELSE
           CALL vxc_of_r(dft_control%xc_control,derivative_set,nspins,na,rho_s)
         ENDIF

!        *** Sum up the exc contribution  to exc_s ***
         exc_s = SUM ( derivative_set%data(:,0)*weight(:) )
         rho_atom%exc_s = rho_atom%exc_s + exc_s
!        ***  vxc_s(ia,ir) 
         DO ispin = 1,nspins
           DO ia =1,na
             vxc_s(ia,ir,ispin) =  derivative_set%data(ia,ispin)*weight(ia)
           END DO
         END DO

!        *** if gradient_f :  vxg_s(ia,ir)
         IF (gradient_f) THEN
           DO ispin = 1,nspins
             ipot = nspins + ispin
             DO ia =1,na
               vxg_s(ia,ir,ispin) =  derivative_set%data(ia,ipot)*weight(ia)
             END DO
           END DO
         END IF

       END DO  !nr 

!      *** Add contributions to the exc energy 
       energy%exc1 = energy%exc1 + rho_atom%exc_h - &
                                 rho_atom%exc_s

!      *** Integration to get the matrix elements relative to the vxc_atom ***
!      *** here the products with the primitives is done: gaVxcgb ***
!      *** internal transformation to get the integral in cartesian Gaussians ***

       atom_kind => atomic_kind_set(ikind)
       IF(gradient_f) THEN
         CALL calculate_gaVxcgb(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,&
                                rho_atom,nspins)
       ELSE
         CALL calculate_gaVxcgb(vxc_h,vxc_s,atom_kind,rho_atom,nspins)
       ENDIF

       NULLIFY(r_h,r_s,dr_h,dr_s)
     END DO ! iat

     DEALLOCATE(derivative_set%data, stat=ierr)
     IF (gradient_f) THEN
        DEALLOCATE(locpot,stat=ierr)
     END IF 

   ENDDO ! ikind

   CALL xc_dset_release(derivative_set, error=error)

 END SUBROUTINE calculate_vxc_atom

! *****************************************************************************

 SUBROUTINE  gaVxcgb_noGC(vxc_h,vxc_s,atom_kind,rho_atom,nspins)

! *****************************************************************************

    TYPE(atomic_kind_type), POINTER       :: atom_kind
    TYPE(rho_atom_type),  POINTER         :: rho_atom

    REAL(dp), DIMENSION(:,:,:),POINTER    :: vxc_h,vxc_s
    REAL(dp), DIMENSION(:,:), POINTER     :: vx
    INTEGER, INTENT(IN)                   :: nspins

    CHARACTER(LEN=*), PARAMETER           :: routine_name = "gaVxcgb_noGC"

    TYPE(gto_basis_set_type), POINTER     :: orb_basis

    TYPE(grid_atom_type), POINTER         :: grid_atom
    TYPE(harmonics_atom_type), &
                           POINTER        :: harmonics
    TYPE(rho_atom_coeff),  DIMENSION(:), &
                           POINTER        :: int_h, int_s
    REAL(dp), DIMENSION(:,:,:), POINTER   :: my_CG

    REAL(dp), DIMENSION(:), ALLOCATABLE   :: g1,g2
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gg,&
                                             gVg_h, gVg_s,&
                                             matso_h, matso_s,&
                                             matco_h,matco_s 
    REAL(dp), DIMENSION(:,:), POINTER     :: zet
    INTEGER, DIMENSION(:), POINTER        :: lmin,lmax,npgf

    INTEGER                               :: ia,ic,icg,ico1,&
                                             ico1_last,ico2,&
                                             ico2_last,ipgf1,ipgf2,iset1,&
                                             iset2,iso,iso1,iso2,&
                                             ispin,istat,l,ld,&
                                             lmin12,lmax12,maxso, maxl, m1,&
                                             m2,n1,n2,na, ngau1,ngau2,nngau1,nr, nset,&
                                             size1

! *****************************************************************************

    NULLIFY(lmin,lmax,npgf,zet,my_CG,harmonics,grid_atom) 

    CALL get_atomic_kind(atomic_kind=atom_kind, orb_basis_set=orb_basis,&
                         harmonics=harmonics,&
                         grid_atom=grid_atom)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxso=maxso,maxl=maxl,npgf=npgf,&
                           nset=nset,zet=zet)

    nr = grid_atom%nr
    na = grid_atom%ng_sphere 
    my_CG => harmonics%my_CG

    ALLOCATE(g1(nr),g2(nr),gg(nr,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "g1,g2,gg")

    ALLOCATE(gVg_h(na,0:2*maxl),gVg_s(na,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gVg_h,gVg_s")

    ALLOCATE(matso_h(nsoset(maxl),nsoset(maxl)),&
             matso_s(nsoset(maxl),nsoset(maxl)),STAT=istat)
    ALLOCATE(matco_h(ncoset(maxl),ncoset(maxl)),&
             matco_s(ncoset(maxl),ncoset(maxl)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "matco,matso")

    NULLIFY(int_h,int_s)

    CALL get_rho_atom(rho_atom=rho_atom,ga_Vxc_gb_h=int_h,ga_Vxc_gb_s=int_s)

    g1 = 0.0_dp
    g2 = 0.0_dp
    m1 = 0
    DO iset1 = 1,nset
      n1 = nsoset(lmax(iset1))
      m2 = 0
      DO iset2 = 1,nset

        n2 = nsoset(lmax(iset2)) 
        DO ipgf1  = 1,npgf(iset1)
          ngau1  = n1*(ipgf1-1)+m1 
          size1  = nsoset(lmax(iset1))-nsoset(lmin(iset1)-1)
          nngau1 = nsoset(lmin(iset1)-1)+ngau1

          g1(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))
          DO ipgf2 = 1,npgf(iset2)
            ngau2 = n2*(ipgf2-1)+m2

            g2(1:nr) = EXP(-zet(ipgf2,iset2)*grid_atom%rad2(1:nr))
            lmin12 = lmin(iset1)+lmin(iset2)
            lmax12 = lmax(iset1)+lmax(iset2)

            gg = 0.0_dp
            IF (lmin12 == 0) THEN
              gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
            ELSE
              gg(1:nr,lmin12) = grid_atom%rad(1:nr)**lmin12*g1(1:nr)*g2(1:nr)
            END IF
            DO l=lmin12+1,lmax12
              gg(1:nr,l) = grid_atom%rad(1:nr)*gg(:,l-1)
            END DO

            DO ispin=1,nspins
              ld = lmax12+1
              vx => vxc_h(:,:,ispin)
              CALL dgemm('N','N',na,ld,nr,1.0_dp,vx(1:na,1:nr),na,&
                          gg(1:nr,0:lmax12),nr,0.0_dp,gVg_h(1:na,0:lmax12),na)
              vx => vxc_s(:,:,ispin)
              CALL dgemm('N','N',na,ld,nr,1.0_dp,vx(1:na,1:nr),na,&
                          gg(1:nr,0:lmax12),nr,0.0_dp,gVg_s(1:na,0:lmax12),na)

              matso_h = 0.0_dp
              matso_s = 0.0_dp
              DO iso = 1,harmonics%max_iso_not0
                DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                  iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                  iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                  l = indso(1,iso1) + indso(1,iso2)
                  DO ia = 1,na
                    matso_h(iso1,iso2) = matso_h(iso1,iso2) + & 
                            gVg_h(ia,l)*&
                            my_CG(iso1,iso2,iso)*&
                            harmonics%slm(ia,iso)
                    matso_s(iso1,iso2) = matso_s(iso1,iso2) + & 
                            gVg_s(ia,l)*&
                            my_CG(iso1,iso2,iso)*&
                            harmonics%slm(ia,iso)
                  END DO
                END DO
              END DO

!             Write in the global matrix
              DO ic = nsoset(lmin(iset2)-1)+1,nsoset(lmax(iset2))
                iso1 = nsoset(lmin(iset1)-1)+1
                iso2 = ngau2+ic
                CALL dcopy(size1,matso_h(iso1,ic),1,&
                           int_h(ispin)%r_coef(nngau1+1,iso2),1)
                CALL dcopy(size1,matso_s(iso1,ic),1,&
                           int_s(ispin)%r_coef(nngau1+1,iso2),1)
               END DO

            END DO ! ispin

          END DO ! ipfg2
        END DO ! ipfg1
        m2 = m2 + maxso
      END DO ! iset2
      m1 = m1 + maxso
    END DO  ! iset1

!  write(*,*) 'int'
!  do ic = 1,maxso*nset
!     write(*,'(i3,32f10.5)') ic,int_h(1)%r_coef(ic,:)
!  enddo
!    stop 'gVg'
   DEALLOCATE(g1,g2,gg,matco_h,matco_s,matso_h,matso_s,gVg_s,gVg_h)

 END SUBROUTINE gaVxcgb_noGC

! *****************************************************************************

 SUBROUTINE  gaVxcgb_GC(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,rho_atom,nspins)

! *****************************************************************************

    TYPE(atomic_kind_type), POINTER       :: atom_kind
    TYPE(rho_atom_type),  POINTER         :: rho_atom

    REAL(dp), DIMENSION(:,:,:),POINTER    :: vxc_h,vxc_s,&
                                             vxg_h,vxg_s
    INTEGER, INTENT(IN)                   :: nspins

    CHARACTER(LEN=*), PARAMETER           :: routine_name = "gaVxcgb_GC"

    TYPE(gto_basis_set_type), POINTER     :: orb_basis

    TYPE(grid_atom_type), POINTER         :: grid_atom
    TYPE(harmonics_atom_type),&
                           POINTER        :: harmonics
    TYPE(rho_atom_coeff),  DIMENSION(:), &
                           POINTER        :: int_h, int_s
    TYPE(rho_atom_coeff),DIMENSION(:),& 
                  POINTER                 :: r_h,r_s,dr_h,dr_s
    REAL(dp)                              :: rad,urad,usin_ia
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: g1,g2
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gg,dgg,&
                                             gVXCg_h, gVXCg_s,&
                                             matso_h, matso_s,&
                                             matco_h,matco_s 
    REAL(dp), DIMENSION(:,:), POINTER     :: zet
    REAL(dp), DIMENSION(:,:,:), POINTER   :: my_CG

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: gVXGg_h, gVXGg_s
    REAL(dp), DIMENSION(:,:,:),ALLOCATABLE :: drho_h,drho_s
    INTEGER, DIMENSION(:), POINTER        :: lmin,lmax,npgf

    INTEGER                               :: ia,ic,icg,ico1,&
                                             ico1_last,ico2,&
                                             ico2_last,ipgf1,ipgf2,ir,iset1,&
                                             iset2,iso,iso1,iso2,&
                                             ispin,istat,l,&
                                             lmin12,lmax12,maxso, maxl, m1,&
                                             m2,n1,n2,na, ngau1, ngau2, nngau1, nr, nset,&
                                             size1

! *****************************************************************************

    NULLIFY(lmin,lmax,npgf,zet,my_CG,harmonics,grid_atom,r_h,r_s,dr_h,dr_s) 

    CALL get_atomic_kind(atomic_kind=atom_kind, orb_basis_set=orb_basis,&
                         harmonics=harmonics,&
                         grid_atom=grid_atom)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxso=maxso,maxl=maxl,npgf=npgf,&
                           nset=nset,zet=zet)

    nr = grid_atom%nr
    na = grid_atom%ng_sphere 
    my_CG => harmonics%my_CG

    ALLOCATE(g1(nr),g2(nr),gg(nr,2*maxl),dgg(nr,2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "g1,g2,gg,dgg")

    ALLOCATE(gVXCg_h(na,2*maxl),gVXCg_s(na,2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gVXCg_h,gVXCg_s")
    ALLOCATE(gVXGg_h(3,na,2*maxl),gVXGg_s(3,na,2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gVXCg_h,gVXCg_s")
    Allocate(drho_h(3,nr,na),drho_s(3,nr,na),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                      "drho_h,drho_s")

    NULLIFY(int_h,int_s)

    ALLOCATE(matso_h(nsoset(maxl),nsoset(maxl)),&
             matso_s(nsoset(maxl),nsoset(maxl)),STAT=istat)
    ALLOCATE(matco_h(ncoset(maxl),ncoset(maxl)),&
             matco_s(ncoset(maxl),ncoset(maxl)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "matco,matso")

    CALL get_rho_atom(rho_atom=rho_atom,&
                      rho_rad_h=r_h,&
                      rho_rad_s=r_s,&
                      drho_rad_h=dr_h,&
                      drho_rad_s=dr_s,& 
                      ga_Vxc_gb_h=int_h,ga_Vxc_gb_s=int_s)

    DO ispin=1,nspins
    drho_h = 0.0_dp
    drho_s = 0.0_dp
!   calculate rho and drho on the points of the sphere
    DO iso = 1,harmonics%max_iso_not0
       DO ia = 1,na
         usin_ia = 1.0_dp/grid_atom%sin_teta(ia) 
         DO ir = 1,nr
             rad = grid_atom%rad(ir)
             urad = 1.0_dp/rad
!            *** Derivatives in spherical coordinates ***
!            *** d/dr
             drho_h(1,ir,ia) = drho_h(1,ir,ia) + &
                                     dr_h(ispin)%r_coef(ir,iso)*&
                                     harmonics%slm(ia,iso)
             drho_s(1,ir,ia) = drho_s(1,ir,ir) + &
                                     dr_s(ispin)%r_coef(ir,iso)*&
                                     harmonics%slm(ia,iso)

!            *** d/d teta
             drho_h(2,ir,ia) = drho_h(2,ir,ia) + &
                                     r_h(ispin)%r_coef(ir,iso)*&
                                     urad*&
                                     harmonics%dslm(1,ia,iso) 
             drho_s(2,ir,ia) = drho_s(2,ir,ia) + &
                                     r_s(ispin)%r_coef(ir,iso)*&
                                     urad*&
                                     harmonics%dslm(1,ia,iso)

!            *** d/d phi
             drho_h(3,ir,ia) = drho_h(3,ir,ia) + &
                                       r_h(ispin)%r_coef(ir,iso)*&
                                       urad*usin_ia*&
                                       harmonics%dslm(2,ia,iso)
             drho_s(3,ir,ia) = drho_s(3,ir,ia) + &
                                       r_s(ispin)%r_coef(ir,iso)*&
                                       urad*usin_ia*&
                                       harmonics%dslm(2,ia,iso)
        END DO            ! ir
      END DO     ! ia
    END DO       ! iso

    g1 = 0.0_dp
    g2 = 0.0_dp
    m1 = 0
    DO iset1 = 1,nset
      n1 = nsoset(lmax(iset1))
      m2 = 0
      DO iset2 = 1,nset

        n2 = nsoset(lmax(iset2)) 
        DO ipgf1  = 1,npgf(iset1)
          ngau1  = n1*(ipgf1-1)+m1
          size1  = nsoset(lmax(iset1))-nsoset(lmin(iset1)-1)
          nngau1 = nsoset(lmin(iset1)-1)+ngau1

          g1(1:nr) = EXP(-zet(ipgf1,iset1)*grid_atom%rad2(1:nr))
          DO ipgf2 = 1,npgf(iset2)
            ngau2 = n2*(ipgf2-1)+m2

            g2(1:nr) = EXP(-zet(ipgf2,iset2)*grid_atom%rad2(1:nr))
            lmin12 = lmin(iset1)+lmin(iset2)
            lmax12 = lmax(iset1)+lmax(iset2)

            gg = 0.0_dp
            dgg = 0.0_dp
            IF (lmin12 == 0) THEN
              gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
            ELSE
              gg(1:nr,lmin12) = grid_atom%rad(1:nr)**lmin12*g1(1:nr)*g2(1:nr)
            END IF
            DO l=lmin12+1,lmax12
              gg(1:nr,l) = grid_atom%rad(1:nr)*gg(:,l-1)
            END DO

              gVXCg_h = 0.0_dp
              gVXCg_s = 0.0_dp
              DO l = lmin12,lmax12
                DO ia = 1,na
                  DO ir = 1,nr
                     gVXCg_h(ia,l) = gVXCg_h(ia,l) + &
                                    gg(ir,l) * vxc_h(ia,ir,ispin)
                     gVXCg_s(ia,l) = gVXCg_s(ia,l) + &
                                    gg(ir,l) * vxc_s(ia,ir,ispin)

                     urad = 1.0_dp/grid_atom%rad(ir)
                     usin_ia = urad/grid_atom%sin_teta(ia)

                     gVXGg_h(1,ia,l) = gVXGg_h(1,ia,l) + &
                                    drho_s(1,ir,ia)*vxg_s(ia,ir,ispin)*&
                                    dgg(ir,l)*0.5_dp
                     gVXGg_h(2,ia,l) = gVXGg_h(2,ia,l) + &
                                    drho_s(2,ir,ia)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*0.5_dp*urad
                     gVXGg_h(3,ia,l) = gVXGg_h(3,ia,l) + &
                                    drho_s(3,ir,ia)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*0.5_dp*usin_ia

                     gVXGg_s(1,ia,l) = gVXGg_s(1,ia,l) + &
                                    drho_s(1,ir,ia)*vxg_s(ia,ir,ispin)*&
                                    dgg(ir,l)*0.5_dp
                     gVXGg_s(2,ia,l) = gVXGg_s(2,ia,l) + &
                                    drho_s(2,ir,ia)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*0.5_dp*urad
                     gVXGg_s(3,ia,l) = gVXGg_s(3,ia,l) + &
                                    drho_s(3,ir,ia)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*0.5_dp*usin_ia
                  ENDDO ! ir
                ENDDO ! ia
              ENDDO  ! l

              matso_h = 0.0_dp
              matso_s = 0.0_dp
              DO iso = 1,harmonics%max_iso_not0
                DO icg = 1,harmonics%ncg(iso,iset1,iset2)
                  iso1 = harmonics%ind_cg_a(icg,iso,iset1,iset2)
                  iso2 = harmonics%ind_cg_b(icg,iso,iset1,iset2)
                  l = indso(1,iso1) + indso(1,iso2)
                  DO ia = 1,na
                    matso_h(iso1,iso2) = matso_h(iso1,iso2) + & 
                            ((gVXCg_h(ia,l)+gVXGg_h(1,ia,l))*&
                            harmonics%slm(ia,iso)+&
                            gVXGg_h(2,ia,l)*harmonics%dslm(1,ia,iso)+&
                            gVXGg_h(3,ia,l)*harmonics%dslm(2,ia,iso))*&
                            my_CG(iso1,iso2,iso)
                    matso_s(iso1,iso2) = matso_s(iso1,iso2) + & 
                            ((gVXCg_s(ia,l)+gVXGg_h(1,ia,l))*&
                            harmonics%slm(ia,iso)+&
                            gVXGg_s(2,ia,l)*harmonics%dslm(1,ia,iso)+&
                            gVXGg_s(3,ia,l)*harmonics%dslm(2,ia,iso))*&
                            my_CG(iso1,iso2,iso)
                  END DO
                END DO
              END DO

!             Write in the global matrix
              DO ic = nsoset(lmin(iset2)-1)+1,nsoset(lmax(iset2))
                iso1 = nsoset(lmin(iset1)-1)+1
                iso2 = ngau2+ic
                CALL dcopy(size1,matso_h(iso1,ic),1,&
                           int_h(ispin)%r_coef(nngau1+1,iso2),1)
                CALL dcopy(size1,matso_s(iso1,ic),1,&
                           int_s(ispin)%r_coef(nngau1+1,iso2),1)
               END DO

          END DO ! ipfg2
        END DO ! ipfg1
        m2 = m2 + maxso
      END DO ! iset2
      m1 = m1 + maxso
    END DO  ! iset1
    END DO ! ispin

 END SUBROUTINE gaVxcgb_GC

! *****************************************************************************

END MODULE qs_vxc_atom

