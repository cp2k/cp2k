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
  USE kinds,                           ONLY: dp,dp_size
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
      OPTIONAL                          :: error
    TYPE(qs_environment_type), POINTER  :: qs_env

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
                  POINTER               :: r_h,r_s
    TYPE(rho_atom_coeff),DIMENSION(:),&
                  POINTER               :: dr_h,dr_s
    TYPE(rho_atom_coeff),DIMENSION(:,:),&
                  POINTER               :: r_h_d,r_s_d
    TYPE(xc_derivative_set_type), &
                             POINTER    :: derivative_set

    REAL(dp), DIMENSION(:), POINTER     :: weight 
    REAL(dp), DIMENSION(:,:), &
                              POINTER   :: rho_h,rho_s
    REAL(dp), DIMENSION(:,:), &
                              POINTER   :: oodrho_h,oodrho_s
    REAL(dp), DIMENSION(:,:,:), &
                              POINTER   :: vxc_h,vxc_s,&
                                           vxg_h,vxg_s
    REAL(dp), DIMENSION(:,:,:,:), &
                              POINTER   :: drho_h,drho_s
    LOGICAL                             :: cross_term,&
                                           gradient_f,lsd,paw_atom
    REAL(dp)                            :: exc_h,exc_s,rad,urad,usin_ia,rtot
    REAL(dp), ALLOCATABLE, &
                        DIMENSION(:,:)  :: locpot
    INTEGER                             :: handle,i,ia,iat,iatom, ierr,&
                                           ikind,ipot,ir,&
                                           iso,ispin,lp,na,natom,nkind,&
                                           nr,npot,nspins,order 
    INTEGER, DIMENSION(:), POINTER      :: atom_list

!   ***************************************************************************

    CALL timeset("calculate_vxc_atom","I","",handle)

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

     NULLIFY(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
     NULLIFY(rho_h,drho_h,rho_s,drho_s,weight,oodrho_h,oodrho_s)
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
        CALL reallocate(drho_h,1,4,1,na,1,nr,1,nspins)
        CALL reallocate(drho_s,1,4,1,na,1,nr,1,nspins)
        CALL reallocate(oodrho_h,1,na,1,nspins)
        CALL reallocate(oodrho_s,1,na,1,nspins)
        CALL reallocate(vxg_h,1,na,1,nr,1,nspins)
        CALL reallocate(vxg_s,1,na,1,nr,1,nspins)
     END IF

     DO iat = 1,natom
       iatom = atom_list(iat)

       rho_atom => rho_atom_set(iatom)
       IF(gradient_f) THEN
         CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,&
                           rho_rad_s=r_s,drho_rad_h=dr_h,&
                           drho_rad_s=dr_s,rho_rad_h_d=r_h_d,&
                           rho_rad_s_d=r_s_d)
         drho_h = 0.0_dp
         drho_s = 0.0_dp
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
         urad = grid_atom%oorad2l(ir,1)
         rho_h = 0.0_dp
         rho_s = 0.0_dp
!        calculate rho and drho on the points of the sphere
         IF(gradient_f) THEN
           oodrho_h = 0.0_dp
           oodrho_s = 0.0_dp
           DO ispin = 1,nspins
             DO iso = 1,harmonics%max_iso_not0
                DO ia = 1,na
                  usin_ia = grid_atom%usin_theta(ia) 
                  rho_h(ia,ispin) = rho_h(ia,ispin) + & 
                       r_h(ispin)%r_coef(ir,iso)*&
                       harmonics%slm(ia,iso)
                  rho_s(ia,ispin) = rho_s(ia,ispin) + & 
                       r_s(ispin)%r_coef(ir,iso)*&
                       harmonics%slm(ia,iso)

!dbg
!                 *** Derivatives in spherical coordinates ***
!                 *** d/dr
!                  drho_h(1,ia,ir,ispin) = drho_h(1,ia,irispin) + &
!                                       dr_h(ispin)%r_coef(ir,iso)*&
!                                       harmonics%slm(ia,iso)
!                  drho_s(1,ia,ir,ispin) = drho_s(1,ia,ir,ispin) + &
!                                       dr_s(ispin)%r_coef(ir,iso)*&
 !                                      harmonics%slm(ia,iso)
!
!                 *** d/d theta
!                  drho_h(2,ia,ir,ispin) = drho_h(2,ia,ir,ispin) + &
!                                       r_h(ispin)%r_coef(ir,iso)*&
!                                       urad*&
!                                       harmonics%dslm(1,ia,iso) 
!                  drho_s(2,ia,ir,ispin) = drho_s(2,ia,ir,ispin) + &
!                                       r_s(ispin)%r_coef(ir,iso)*&
!                                       urad*&
!                                       harmonics%dslm(1,ia,iso)

!                 *** d/d phi
!                  drho_h(3,ia,ir,ispin) = drho_h(3,ia,ir,ispin) + &
!                                       r_h(ispin)%r_coef(ir,iso)*&
!                                       urad*usin_ia*&
!                                       harmonics%dslm(2,ia,iso)
!                  drho_s(3,ia,ir,ispin) = drho_s(3,ia,ir,ispin) + &
!                                       r_s(ispin)%r_coef(ir,iso)*&
!                                       urad*usin_ia*&
!                                       harmonics%dslm(2,ia,iso)


                  drho_h(1,ia,ir,ispin) = drho_h(1,ia,ir,ispin) + & 
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(1,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_h_d(1,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)
!                                       urad*r_h(ispin)%r_coef(ir,iso)*&
!                                       (1.0_dp-harmonics%a(1,ia)*harmonics%a(1,ia))*&
!                                       harmonics%dslm_dxyz(1,ia,iso)
                  drho_h(2,ia,ir,ispin) = drho_h(2,ia,ir,ispin) + & 
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(2,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_h_d(2,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

!                                       urad*r_h(ispin)%r_coef(ir,iso)*&
!                                       (1.0_dp-harmonics%a(2,ia)*harmonics%a(2,ia))*&
!                                       harmonics%dslm_dxyz(2,ia,iso)
                  drho_h(3,ia,ir,ispin) = drho_h(3,ia,ir,ispin) + & 
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(3,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_h_d(3,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

!                                       urad*r_h(ispin)%r_coef(ir,iso)*&
!                                       (1.0_dp-harmonics%a(3,ia)*harmonics%a(3,ia))*&
!                                       harmonics%dslm_dxyz(3,ia,iso)

                  drho_s(1,ia,ir,ispin) = drho_s(1,ia,ir,ispin) + & 
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(1,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_s_d(1,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

!                                       urad*r_s(ispin)%r_coef(ir,iso)*&
!                                       (1.0_dp-harmonics%a(1,ia)*harmonics%a(1,ia))*&
!                                       harmonics%dslm_dxyz(1,ia,iso)
                  drho_s(2,ia,ir,ispin) = drho_s(2,ia,ir,ispin) + & 
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(2,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_s_d(2,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

!                                       urad*r_s(ispin)%r_coef(ir,iso)*&
!                                       (1.0_dp-harmonics%a(2,ia)*harmonics%a(2,ia))*&
!                                       harmonics%dslm_dxyz(2,ia,iso)
                  drho_s(3,ia,ir,ispin) = drho_s(3,ia,ir,ispin) + & 
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(3,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_s_d(3,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

!                                       urad*r_s(ispin)%r_coef(ir,iso)*&
!                                       (1.0_dp-harmonics%a(3,ia)*harmonics%a(3,ia))*&
!                                       harmonics%dslm_dxyz(3,ia,iso)



                  drho_h(4,ia,ir,ispin) = SQRT(&
                               drho_h(1,ia,ir,ispin)*drho_h(1,ia,ir,ispin)+& 
                               drho_h(2,ia,ir,ispin)*drho_h(2,ia,ir,ispin)+&
                               drho_h(3,ia,ir,ispin)*drho_h(3,ia,ir,ispin))
                  oodrho_h(ia,ispin) = -1.0_dp/MIN(-drho_h(4,ia,ir,ispin),-1.E-12_dp)
                  drho_s(4,ia,ir,ispin) = SQRT(&
                               drho_s(1,ia,ir,ispin)*drho_s(1,ia,ir,ispin)+& 
                               drho_s(2,ia,ir,ispin)*drho_s(2,ia,ir,ispin)+&
                               drho_s(3,ia,ir,ispin)*drho_s(3,ia,ir,ispin))
                  oodrho_s(ia,ispin) = -1.0_dp/MIN(-drho_s(4,ia,ir,ispin),-1.E-12_dp)
                END DO
             END DO
!  WRITE(*,*) 'iatom ', iatom, 'ir=10', grid_atom%rad(ir)
!  DO ia = 1,na
!    WRITE(50,'(I4,11f16.10)') ia,grid_atom%theta(ia),grid_atom%phi(ia),&
!             harmonics%slm(ia,1:9)   
!    WRITE(51,'(I4,11f16.10)') ia,grid_atom%theta(ia),grid_atom%phi(ia),&
!             harmonics%dslm(1,ia,1:9)   
!    WRITE(52,'(I4,11f16.10)') ia,grid_atom%theta(ia),grid_atom%phi(ia),&
!             harmonics%dslm(2,ia,1:9)   
! if(iatom == 3) then
!    WRITE(51,'(I4,7f16.10)') ia,grid_atom%rad(ir),rho_h(ia,ispin),&
!                  drho_h(ia,1,ispin),drho_h(ia,2,ispin),drho_h(ia,3,ispin),drho_h(ia,4,ispin)
! end if
!  end do 
!    WRITE(46,*)
           END DO
!  stop 'gra'
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
         END IF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        *** hard atom density ***
         CALL dcopy(na*npot, 0.0_dp, 0, derivative_set%data, 1)

         IF(gradient_f) THEN
           CALL vxc_of_r(dft_control%xc_control,derivative_set,nspins,na,rho_h,&
                         drho=drho_h,ir=ir,locpot=locpot,error=error)
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
               vxg_h(ia,ir,ispin) =  derivative_set%data(ia,ipot)*weight(ia)*oodrho_h(ia,ispin)
             END DO
           END DO
         END IF

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        *** soft atom density ***
         CALL dcopy(na*npot, 0.0_dp, 0, derivative_set%data, 1)

         IF(gradient_f) THEN
           CALL vxc_of_r(dft_control%xc_control,derivative_set,nspins,na,rho_s,&
                         drho=drho_s,ir=ir,locpot=locpot,error=error)
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
               vxg_s(ia,ir,ispin) =  derivative_set%data(ia,ipot)*weight(ia)*oodrho_s(ia,ispin)
             END DO
           END DO
         END IF

       END DO  !nr 

!      *** Add contributions to the exc energy 
       energy%exc1 = energy%exc1 + rho_atom%exc_h - &
                                 rho_atom%exc_s

!     WRITE(*,*) rho_atom%exc_h, rho_atom%exc_s, energy%exc1

!      *** Integration to get the matrix elements relative to the vxc_atom ***
!      *** here the products with the primitives is done: gaVxcgb ***
!      *** internal transformation to get the integral in cartesian Gaussians ***

       atom_kind => atomic_kind_set(ikind)
       IF(gradient_f) THEN
         CALL calculate_gaVxcgb(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,&
                                rho_atom,drho_h,drho_s,nspins)
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
!   STOP 'exc'

   CALL timestop(0.0_dp,handle)
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
                                             matso_h, matso_s
    REAL(dp), DIMENSION(:,:), POINTER     :: zet
    INTEGER, DIMENSION(:), POINTER        :: lmin,lmax,npgf

    INTEGER                               :: ia,ic,icg,&
                                             ipgf1,ipgf2,iset1,&
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
                                     "g1,g2,gg",(2*nr+2*maxl*nr)*dp_size)

    ALLOCATE(gVg_h(na,0:2*maxl),gVg_s(na,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gVg_h,gVg_s",4*na*maxl*dp_size)

    ALLOCATE(matso_h(nsoset(maxl),nsoset(maxl)),&
             matso_s(nsoset(maxl),nsoset(maxl)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                    "matso_h,matso_s",2*nsoset(maxl)*nsoset(maxl)*dp_size)

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
              gg(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12)*g1(1:nr)*g2(1:nr)
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

   DEALLOCATE(g1,g2,gg,matso_h,matso_s,gVg_s,gVg_h,STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "g1,g2,gg,matso_h,matso_s,gVg_s,gVg_h")
 END SUBROUTINE gaVxcgb_noGC

! *****************************************************************************

 SUBROUTINE  gaVxcgb_GC(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,rho_atom,drho_h,drho_s,nspins)

! *****************************************************************************

    TYPE(atomic_kind_type), POINTER       :: atom_kind
    TYPE(rho_atom_type),  POINTER         :: rho_atom

    REAL(dp), DIMENSION(:,:,:),POINTER    :: vxc_h,vxc_s,&
                                             vxg_h,vxg_s
    REAL(dp), DIMENSION(:,:,:,:), &
                                POINTER   :: drho_h,drho_s
    INTEGER, INTENT(IN)                   :: nspins

    CHARACTER(LEN=*), PARAMETER           :: routine_name = "gaVxcgb_GC"

    TYPE(gto_basis_set_type), POINTER     :: orb_basis

    TYPE(grid_atom_type), POINTER         :: grid_atom
    TYPE(harmonics_atom_type),&
                           POINTER        :: harmonics
    TYPE(rho_atom_coeff),  DIMENSION(:), &
                           POINTER        :: int_h, int_s
    REAL(dp)                              :: rad,urad,usin_ia
    REAL(dp), DIMENSION(:), ALLOCATABLE   :: g1,g2
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gg,dgg,&
                                             gVXCg_h, gVXCg_s,&
                                             matso_h, matso_s

    REAL(dp), DIMENSION(:,:), POINTER     :: zet
    REAL(dp), DIMENSION(:,:,:), POINTER   :: my_CG
    REAL(dp), DIMENSION(:,:,:,:), POINTER   :: my_CG_dxyz

    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: gVXGg_h, gVXGg_s
    INTEGER, DIMENSION(:), POINTER        :: lmin,lmax,npgf

    INTEGER                               :: ia,ic,icg,&
                                             ipgf1,ipgf2,ir,iset1,&
                                             iset2,iso,iso1,iso2,&
                                             ispin,istat,l,&
                                             lmin12,lmax12,maxso, maxl, m1,&
                                             m2,n1,n2,na, ngau1, ngau2, nngau1, nr, nset,&
                                             size1

! *****************************************************************************

    NULLIFY(lmin,lmax,npgf,zet,my_CG,my_CG_dxyz,harmonics,grid_atom) 

    CALL get_atomic_kind(atomic_kind=atom_kind, orb_basis_set=orb_basis,&
                         harmonics=harmonics,&
                         grid_atom=grid_atom)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxso=maxso,maxl=maxl,npgf=npgf,&
                           nset=nset,zet=zet)

    nr = grid_atom%nr
    na = grid_atom%ng_sphere 
    my_CG => harmonics%my_CG
    my_CG_dxyz => harmonics%my_CG_dxyz

    ALLOCATE(g1(nr),g2(nr),gg(nr,0:2*maxl),dgg(nr,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "g1,g2,gg,dgg",(2*nr+4*nr*maxl)*dp_size)

    ALLOCATE(gVXCg_h(na,0:2*maxl),gVXCg_s(na,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gVXCg_h,gVXCg_s",4*na*maxl*dp_size)
    ALLOCATE(gVXGg_h(3,na,0:2*maxl),gVXGg_s(3,na,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                     "gVXCg_h,gVXCg_s",12*na*maxl*dp_size)

    NULLIFY(int_h,int_s)

    ALLOCATE(matso_h(nsoset(maxl),nsoset(maxl)),&
             matso_s(nsoset(maxl),nsoset(maxl)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "matso_h,matso_s",2*nsoset(maxl)*nsoset(maxl)*dp_size)

    CALL get_rho_atom(rho_atom=rho_atom,&
                      ga_Vxc_gb_h=int_h,ga_Vxc_gb_s=int_s)

    DO ispin=1,nspins

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
              gg(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12)*g1(1:nr)*g2(1:nr)
!              dgg(1:nr,lmin12) = REAL(lmin12,dp)*grid_atom%oorad2l(1:nr,1)*&
!                                gg(1:nr,lmin12)
            END IF
            DO l=lmin12+1,lmax12
!              dgg(1:nr,l) = REAL(l,dp)*gg(1:nr,l-1)
              gg(1:nr,l) = grid_atom%rad(1:nr)*gg(:,l-1)
              dgg(1:nr,l-1) = dgg(1:nr,l-1) -2.0_dp*(zet(ipgf1,iset1)+ &
                             zet(ipgf2,iset2))*gg(1:nr,l)
            END DO
            dgg(1:nr,lmax12) = dgg(1:nr,lmax12) -2.0_dp*(zet(ipgf1,iset1)+ &
                             zet(ipgf2,iset2))*grid_atom%rad(1:nr)*&
                             gg(1:nr,lmax12)

              gVXCg_h = 0.0_dp
              gVXCg_s = 0.0_dp
              gVXGg_h = 0.0_dp
              gVXGg_s = 0.0_dp
              DO l = lmin12,lmax12
                DO ia = 1,na
                  DO ir = 1,nr
                     gVXCg_h(ia,l) = gVXCg_h(ia,l) + &
                                    gg(ir,l) * vxc_h(ia,ir,ispin)+&
                                    dgg(ir,l)*vxg_h(ia,ir,ispin)*&
                                   (drho_h(1,ia,ir,ispin)*harmonics%a(1,ia)+&
                                    drho_h(2,ia,ir,ispin)*harmonics%a(2,ia)+&
                                    drho_h(3,ia,ir,ispin)*harmonics%a(3,ia))

                     gVXCg_s(ia,l) = gVXCg_s(ia,l) + &
                                    gg(ir,l) * vxc_s(ia,ir,ispin)+&
                                    dgg(ir,l)*vxg_s(ia,ir,ispin)*&
                                   (drho_s(1,ia,ir,ispin)*harmonics%a(1,ia)+&
                                    drho_s(2,ia,ir,ispin)*harmonics%a(2,ia)+&
                                    drho_s(3,ia,ir,ispin)*harmonics%a(3,ia))

                     urad = grid_atom%oorad2l(ir,1)


                     gVXGg_h(1,ia,l) = gVXGg_h(1,ia,l) + &
                                    drho_h(1,ia,ir,ispin)*vxg_h(ia,ir,ispin)*&
                                    gg(ir,l)*urad
                     gVXGg_h(2,ia,l) = gVXGg_h(2,ia,l) + &
                                    drho_h(2,ia,ir,ispin)*vxg_h(ia,ir,ispin)*&
                                    gg(ir,l)*urad
                     gVXGg_h(3,ia,l) = gVXGg_h(3,ia,l) + &
                                    drho_h(3,ia,ir,ispin)*vxg_h(ia,ir,ispin)*&
                                    gg(ir,l)*urad

                     gVXGg_s(1,ia,l) = gVXGg_s(1,ia,l) + &
                                    drho_s(1,ia,ir,ispin)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*urad
                     gVXGg_s(2,ia,l) = gVXGg_s(2,ia,l) + &
                                    drho_s(2,ia,ir,ispin)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*urad
                     gVXGg_s(3,ia,l) = gVXGg_s(3,ia,l) + &
                                    drho_s(3,ia,ir,ispin)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*urad
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
                                         gVXCg_h(ia,l)*&
                                         harmonics%slm(ia,iso)*&
                                         my_CG(iso1,iso2,iso)
                    matso_s(iso1,iso2) = matso_s(iso1,iso2) + & 
                                         gVXCg_s(ia,l)*&
                                         harmonics%slm(ia,iso)*&
                                         my_CG(iso1,iso2,iso)
                  END DO
                END DO
                DO icg = 1,harmonics%ncg_dxyz(iso,iset1,iset2)
                  iso1 = harmonics%ind_cg_a_dxyz(icg,iso,iset1,iset2)
                  iso2 = harmonics%ind_cg_b_dxyz(icg,iso,iset1,iset2)
                  l = indso(1,iso1) + indso(1,iso2)
                  DO ia = 1,na
                    matso_h(iso1,iso2) = matso_h(iso1,iso2) + &
                            (gVXGg_h(1,ia,l)*my_CG_dxyz(1,iso1,iso2,iso)+&
                             gVXGg_h(2,ia,l)*my_CG_dxyz(2,iso1,iso2,iso)+&
                             gVXGg_h(3,ia,l)*my_CG_dxyz(3,iso1,iso2,iso))*&
                             harmonics%slm(ia,iso)

                    matso_s(iso1,iso2) = matso_s(iso1,iso2) + &
                            (gVXGg_s(1,ia,l)*my_CG_dxyz(1,iso1,iso2,iso)+&
                             gVXGg_s(2,ia,l)*my_CG_dxyz(2,iso1,iso2,iso)+&
                             gVXGg_s(3,ia,l)*my_CG_dxyz(3,iso1,iso2,iso))*&
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

          END DO ! ipfg2
        END DO ! ipfg1
        m2 = m2 + maxso
      END DO ! iset2
      m1 = m1 + maxso
    END DO  ! iset1
    END DO ! ispin

    DEALLOCATE(g1,g2,gg,dgg,matso_h,matso_s,gVXCg_h,gVXCg_s,gVXGg_h,gVXGg_s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                     "g1,g2,gg,matso_h,matso_s,gVXCg_s,gVXCg_h,gVXGg_h,gVXGg_s")

 END SUBROUTINE gaVxcgb_GC

! *****************************************************************************

END MODULE qs_vxc_atom

