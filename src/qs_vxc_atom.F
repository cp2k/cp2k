!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002 -2004  CP2K developers group                           !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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
  USE cp_control_types,                ONLY: dft_control_type,&
                                             xc_control_type
  USE cp_error_handling,               ONLY: cp_error_type, &
                                             cp_a_l,&
                                             cp_assert
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_to_string
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,dp_size
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,nco,nso,&
                                             ncoset,nsoset
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: harmonics_atom_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE qs_util,                         ONLY: transform_s2c
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  USE xc_atom,                         ONLY: vxc_of_r,&
                                             vxc_of_r_new,&
                                             xc_2nd_deriv_of_r,&
                                             xc_rho_set_atom_update
  USE xc_functionals,                  ONLY: xc_calculate_lda,&
                                             xc_calculate_lsd,&
                                             xc_get_pot_size
  USE xc_derivatives,                  ONLY: xc_functional_get_info
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_release
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_setall,&
                                             xc_rho_cflags_type,&
                                             xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='qs_vxc_atom'

  PUBLIC :: calculate_vxc_atom, calculate_xc_2nd_deriv_atom


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
    TYPE(cp_para_env_type), POINTER     :: para_env
    TYPE(dft_control_type), POINTER     :: dft_control
    TYPE(grid_atom_type), POINTER       :: grid_atom
    TYPE(harmonics_atom_type), POINTER  :: harmonics
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
    TYPE(xc_control_type), POINTER      :: xc_control
    TYPE(xc_derivative_set_type), &
                             POINTER    :: deriv_set
    TYPE(xc_rho_cflags_type)            :: needs
    TYPE(xc_rho_set_type), POINTER      :: rho_set

    REAL(dp), DIMENSION(:), POINTER     :: weight 
    REAL(dp), DIMENSION(:,:), &
                              POINTER   :: rho_h,rho_s
    REAL(dp), DIMENSION(:,:), &
                              POINTER   :: oodrho_h,oodrho_s
    REAL(dp), DIMENSION(:,:,:), &
                              POINTER   :: vxc_h,vxc_s,&
                                           vxg_h,vxg_s
    REAL(dp), DIMENSION(:,:), POINTER   :: vxgc_h,vxgc_s
    REAL(dp), DIMENSION(:,:,:,:), &
                              POINTER   :: drho_h,drho_s
    LOGICAL                             :: cross_term,&
                                           gradient_f,lsd,paw_atom
    REAL(dp)                            :: exc_h,exc_s,rad,urad,usin_ia,rtot
    INTEGER                             :: bo(2),handle,i,ia,iat,iatom, ierr,&
                                           ifunc,ikind,ipot,ir,iso,ispin,&
                                           istat,lp,mepos,na,natom,nkind,&
                                           nr,nspins,num_pe 
    INTEGER, DIMENSION(:), POINTER      :: atom_list
    INTEGER, DIMENSION(2, 3)            :: bounds


!   ***************************************************************************

    CALL timeset("calculate_vxc_atom","I","",handle)

    NULLIFY(dft_control,xc_control,para_env)
    NULLIFY(energy)
    NULLIFY(atomic_kind_set,atom_kind)
    NULLIFY(rho_atom_set,rho_atom)
    NULLIFY(atom_kind)
    NULLIFY(grid_atom,harmonics)
    NULLIFY(atom_list)
    NULLIFY(deriv_set)
    NULLIFY(rho_set)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    energy=energy,&
                    atomic_kind_set=atomic_kind_set,&
                    rho_atom_set=rho_atom_set)

   nspins = dft_control%nspins
   gradient_f=ANY(dft_control%xc_control%gradient_functionals)
   lsd = (nspins==2)

   nkind = SIZE(atomic_kind_set)

!  Put to 0 the energy contribution from the one center XC terms 
   energy%exc1 = 0.0_dp

!  Nullify some pointers fro work-arrays
   NULLIFY(rho_h,drho_h,rho_s,drho_s,weight,oodrho_h,oodrho_s)
   NULLIFY(vxc_h,vxc_s,vxg_h,vxg_s,vxgc_h,vxgc_s)

!  Here starts the loop over all the atoms
   DO ikind = 1,nkind

     atom_kind => atomic_kind_set(ikind)
     CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                        natom=natom,paw_atom=paw_atom,&
                        harmonics=harmonics,grid_atom=grid_atom)
     IF(.NOT. paw_atom) CYCLE


     nr = grid_atom%nr
     na = grid_atom%ng_sphere

!    Prepare the structures needed to calculate and store the xc derivatives
     xc_control => dft_control%xc_control
!    Array dimension: here anly one dimensional arrays are used,
!    i.e. only the first column of deriv_data is read.
!    The other to dimensions  are set to size equal 1
     bounds(1:2,1:3) = 1
     bounds(2,1) = na

     ! create a place where to put the derivatives
     CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)
     ! create the place where to store the argument for the functionals
     CALL xc_rho_set_create(rho_set,bounds,rho_cutoff=xc_control%density_cut,&
            drho_cutoff=xc_control%gradient_cut,tau_cutoff=xc_control%tau_cut,&
            error=error)
     ! assesses the argument needed by the functionals
     CALL xc_rho_cflags_setall(needs,.FALSE.,error=error)
     DO ifunc=1,SIZE(xc_control%functionals)
        CALL xc_functional_get_info(xc_control%functionals(ifunc),lsd=lsd,&
             needs=needs, error=error)
     END DO
     ! allocate the required 3d arrays where to store rho and drho
     CALL xc_rho_set_atom_update(rho_set,needs,nspins,na)
     
     cross_term = .FALSE.
     IF(lsd) THEN
       cross_term = needs%norm_drho .OR. needs%drhoa_drhob
     END IF

     CALL reallocate(rho_h,1,na,1,nspins)
     CALL reallocate(rho_s,1,na,1,nspins)
     CALL reallocate(weight,1,na)
     CALL reallocate(vxc_h,1,na,1,nr,1,nspins)
     CALL reallocate(vxc_s,1,na,1,nr,1,nspins)
  
     IF (gradient_f) THEN
        CALL reallocate(drho_h,1,4,1,na,1,nr,1,nspins)
        CALL reallocate(drho_s,1,4,1,na,1,nr,1,nspins)
        CALL reallocate(oodrho_h,1,na,1,nspins)
        CALL reallocate(oodrho_s,1,na,1,nspins)
        CALL reallocate(vxg_h,1,na,1,nr,1,nspins)
        CALL reallocate(vxg_s,1,na,1,nr,1,nspins)
        IF(cross_term) THEN
          CALL reallocate(vxgc_h,1,na,1,nr)
          CALL reallocate(vxgc_s,1,na,1,nr)
        END IF
     END IF

!    Distribute the atoms of this kind
     num_pe = para_env%num_pe
     mepos  = para_env%mepos
     bo = get_limit( natom, num_pe, mepos )

     DO iat = 1,natom
       iatom = atom_list(iat)
       DO ispin = 1,nspins
         rho_atom_set(iatom)%ga_Vxc_gb_h(ispin)%r_coef = 0.0_dp
         rho_atom_set(iatom)%ga_Vxc_gb_s(ispin)%r_coef = 0.0_dp
       END DO
       rho_atom_set(iatom)%exc_h = 0.0_dp
       rho_atom_set(iatom)%exc_s = 0.0_dp
     END DO

     DO iat = bo(1),bo(2) !1,natom
       iatom = atom_list(iat)

       rho_atom => rho_atom_set(iatom)
       IF(gradient_f) THEN
         NULLIFY(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
         CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,&
                           rho_rad_s=r_s,drho_rad_h=dr_h,&
                           drho_rad_s=dr_s,rho_rad_h_d=r_h_d,&
                           rho_rad_s_d=r_s_d)
         drho_h = 0.0_dp
         drho_s = 0.0_dp
       ELSE
         NULLIFY(r_h,r_s)
         CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,&
                           rho_rad_s=r_s)
       END IF

       rtot = 0.0_dp

       DO ir = 1,nr

         CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f,&
                               xc_control%gradient_cut, &
                               ir, r_h, r_s, rho_h, rho_s, &
                               dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s,  &
                               oodrho_h, oodrho_s, error)

         weight = grid_atom%wr(ir)*grid_atom%wa(1:na)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        *** hard atom density ***

!        Use the new routine to calculate the xc derivatives
!        rho_h and drho_h are given for the radial point ir and
!        all the na points on the angular grid
         CALL vxc_of_r_new(xc_control = xc_control,&
                           rho_set=rho_set,deriv_set=deriv_set,deriv_order=1,needs=needs,&
                           rho=rho_h,drho=drho_h,oodrho=oodrho_h,w=weight,&
                           lsd=lsd,na=na,ir_pnt=ir,&
                           exc=exc_h,vxc=vxc_h,vxg=vxg_h,vxg_cross=vxgc_h,error=error)

!        *** Sum up the exc contribution  to exc_h ***
         rho_atom%exc_h = rho_atom%exc_h + exc_h

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!        *** soft atom density ***
!        Use the new routine to calculate the xc derivatives
         CALL vxc_of_r_new(xc_control=xc_control,&
                           rho_set=rho_set,deriv_set=deriv_set,deriv_order=1,needs=needs,&
                           rho=rho_s,drho=drho_s,oodrho=oodrho_s,w=weight,&
                           lsd=lsd,na=na,ir_pnt=ir,&
                           exc=exc_s,vxc=vxc_s,vxg=vxg_s,vxg_cross=vxgc_s,error=error)

!        *** Sum up the exc contribution  to exc_s ***
         rho_atom%exc_s = rho_atom%exc_s + exc_s

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       END DO  !nr 

!      *** Add contributions to the exc energy 
       energy%exc1 = energy%exc1 + rho_atom%exc_h - rho_atom%exc_s

!     WRITE(*,*) rho_atom%exc_h, rho_atom%exc_s, energy%exc1

!      *** Integration to get the matrix elements relative to the vxc_atom ***
!      *** here the products with the primitives is done: gaVxcgb ***
!      *** internal transformation to get the integral in cartesian Gaussians ***

       atom_kind => atomic_kind_set(ikind)
       IF(cross_term) THEN
         CALL calculate_gaVxcgb(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,&
                                    rho_atom,drho_h,drho_s,nspins,&
                                        vxgc_h=vxgc_h,vxgc_s=vxgc_s)
       ELSE IF(gradient_f) THEN
         CALL calculate_gaVxcgb(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,&
                                      rho_atom,drho_h,drho_s,nspins)
       ELSE
         CALL calculate_gaVxcgb(vxc_h,vxc_s,atom_kind,rho_atom,nspins)
       ENDIF

       NULLIFY(r_h,r_s,dr_h,dr_s)
     END DO ! iat

!    Release the xc structure used to store the xc derivatives
     CALL xc_dset_release(deriv_set, error=error)
     CALL xc_rho_set_release(rho_set,error=error)

   ENDDO ! ikind

   natom = SIZE(rho_atom_set,1)
   DO iatom = 1,natom
     DO ispin = 1,nspins
       CALL mp_sum(rho_atom_set(iatom)%ga_Vxc_gb_h(ispin)%r_coef,para_env%group)
       CALL mp_sum(rho_atom_set(iatom)%ga_Vxc_gb_s(ispin)%r_coef,para_env%group)
     END DO
     CALL mp_sum(rho_atom_set(iatom)%exc_h,para_env%group)
     CALL mp_sum(rho_atom_set(iatom)%exc_s,para_env%group)
   END DO
   CALL mp_sum(energy%exc1,para_env%group)

   DEALLOCATE(rho_h,rho_s,weight,vxc_h,vxc_s,STAT=istat)
   IF (istat /= 0)  CALL stop_memory(routine_name,module_name,__LINE__,&
        "rho_h,rho_s,weight,vxc_h,vxc_s")
   IF (gradient_f) THEN
      DEALLOCATE(drho_h,drho_s,oodrho_h,oodrho_s,vxg_h,vxg_s,STAT=istat)
      IF (istat /= 0)   CALL stop_memory(routine_name,module_name,__LINE__,&
           "drho_h,drho_s,oodrho_h,oodrho_s,vxg_h,vxg_s")
      IF(cross_term) THEN
         DEALLOCATE(vxgc_h,vxgc_s,STAT=istat)
         IF (istat /= 0)   CALL stop_memory(routine_name,module_name,__LINE__,&
              "vxgc_h,vxgc_s")
      END IF
   END IF

   CALL timestop(0.0_dp,handle)

 END SUBROUTINE calculate_vxc_atom

! *****************************************************************************

  SUBROUTINE calc_rho_angular(grid_atom, harmonics, nspins, grad_func, grad_cut, &
                              ir, r_h, r_s, rho_h, rho_s, &
                              dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, &
                              oodrho_h, oodrho_s, error)

   ! ARGUMENTS
   TYPE(grid_atom_type), POINTER       :: grid_atom
   TYPE(harmonics_atom_type), POINTER  :: harmonics
   INTEGER, INTENT(IN)                 :: nspins
   LOGICAL, INTENT(IN)                 :: grad_func
   REAL(KIND=dp), INTENT(IN)           :: grad_cut
   INTEGER, INTENT(IN)                 :: ir
   TYPE(rho_atom_coeff),DIMENSION(:),& 
        POINTER                        :: r_h,r_s
   REAL(KIND=dp), DIMENSION(:,:), &
        POINTER                        :: rho_h, rho_s
   TYPE(rho_atom_coeff),DIMENSION(:),&
        POINTER                        :: dr_h,dr_s
   TYPE(rho_atom_coeff),DIMENSION(:,:),&
        POINTER                        :: r_h_d,r_s_d
   REAL(KIND=dp), DIMENSION(:,:,:,:), &
        POINTER                        :: drho_h,drho_s
   REAL(KIND=dp), DIMENSION(:,:), &
        POINTER                        :: oodrho_h, oodrho_s
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                       :: error

   ! LOCALS
   INTEGER                             :: asize, ispin, iso
   INTEGER                             :: ia, na
   LOGICAL                             :: failure
   REAL(KIND=dp)                       :: rad, urad, usin_ia

   ! PARAMETERS
   CHARACTER(LEN=*), PARAMETER         :: routineN = "calculate_rho_angular", &
                                          routineP = module_name//"/"//routineN

   failure = .FALSE.
   CPPostcondition(ASSOCIATED(r_h),cp_failure_level,routineP,error,failure)
   CPPostcondition(ASSOCIATED(r_s),cp_failure_level,routineP,error,failure)
   CPPostcondition(ASSOCIATED(rho_h),cp_failure_level,routineP,error,failure)
   CPPostcondition(ASSOCIATED(rho_s),cp_failure_level,routineP,error,failure)
   IF (grad_func) THEN
      CPPostcondition(ASSOCIATED(dr_h),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(dr_s),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(r_h_d),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(r_s_d),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(drho_h),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(drho_s),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(oodrho_h),cp_failure_level,routineP,error,failure)
      CPPostcondition(ASSOCIATED(oodrho_s),cp_failure_level,routineP,error,failure)
   END IF
   
   IF (failure) RETURN

   na = grid_atom%ng_sphere
   rad = grid_atom%rad(ir)
   urad = grid_atom%oorad2l(ir,1)
   asize = SIZE(rho_h,1)*SIZE(rho_h,2)
   CALL DCOPY(asize, 0.0_dp, 0, rho_h, 1)
   CALL DCOPY(asize, 0.0_dp, 0, rho_s, 1)

   DO ispin = 1,nspins
      DO iso = 1,harmonics%max_iso_not0
         DO ia = 1,na
            rho_h(ia,ispin) = rho_h(ia,ispin) + & 
                              r_h(ispin)%r_coef(ir,iso)*&
                              harmonics%slm(ia,iso)
            rho_s(ia,ispin) = rho_s(ia,ispin) + & 
                              r_s(ispin)%r_coef(ir,iso)*&
                              harmonics%slm(ia,iso)
         END DO  ! ia
      END DO  ! iso
   END DO  ! ispin

   IF(grad_func) THEN
      asize = SIZE(oodrho_h,1)*SIZE(oodrho_h,2)
      CALL DCOPY(asize, 0.0_dp, 0, oodrho_h, 1)
      CALL DCOPY(asize, 0.0_dp, 0, oodrho_s, 1)
      DO ispin = 1,nspins
         DO iso = 1,harmonics%max_iso_not0
            DO ia = 1,na
               usin_ia = grid_atom%usin_theta(ia) 
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

!              icomponents of the gradient of rho1 hard 
               drho_h(1,ia,ir,ispin) = drho_h(1,ia,ir,ispin) + & 
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(1,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_h_d(1,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_h(2,ia,ir,ispin) = drho_h(2,ia,ir,ispin) + & 
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(2,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_h_d(2,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_h(3,ia,ir,ispin) = drho_h(3,ia,ir,ispin) + & 
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(3,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_h_d(3,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

!              components of the gradient of rho1 soft 
               drho_s(1,ia,ir,ispin) = drho_s(1,ia,ir,ispin) + & 
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(1,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_s_d(1,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_s(2,ia,ir,ispin) = drho_s(2,ia,ir,ispin) + & 
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(2,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_s_d(2,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_s(3,ia,ir,ispin) = drho_s(3,ia,ir,ispin) + & 
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(3,ia)*harmonics%slm(ia,iso) +&
                                       urad*r_s_d(3,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)


               drho_h(4,ia,ir,ispin) = SQRT(&
                    drho_h(1,ia,ir,ispin)*drho_h(1,ia,ir,ispin)+& 
                    drho_h(2,ia,ir,ispin)*drho_h(2,ia,ir,ispin)+&
                    drho_h(3,ia,ir,ispin)*drho_h(3,ia,ir,ispin))
               IF(drho_h(4,ia,ir,ispin) > grad_cut ) THEN
                  oodrho_h(ia,ispin) = 1.0_dp/drho_h(4,ia,ir,ispin)
               ELSE
                  oodrho_h(ia,ispin) =  0.0_dp
               END IF
               
               drho_s(4,ia,ir,ispin) = SQRT(&
                    drho_s(1,ia,ir,ispin)*drho_s(1,ia,ir,ispin)+& 
                    drho_s(2,ia,ir,ispin)*drho_s(2,ia,ir,ispin)+&
                    drho_s(3,ia,ir,ispin)*drho_s(3,ia,ir,ispin))
               IF(drho_h(4,ia,ir,ispin) > grad_cut ) THEN
                  oodrho_s(ia,ispin) = 1.0_dp/drho_s(4,ia,ir,ispin)
               ELSE
                  oodrho_s(ia,ispin) =  0.0_dp
               END IF

            END DO  ! ia
         END DO  ! iso
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
      END DO  ! ispin
   END IF

 END SUBROUTINE calc_rho_angular


! *****************************************************************************

 SUBROUTINE calculate_xc_2nd_deriv_atom(p_env, qs_env, error)

   ! ARGUMENTS
   TYPE(qs_p_env_type), POINTER        :: p_env
   TYPE(qs_environment_type), POINTER  :: qs_env
   TYPE(cp_error_type), INTENT(inout), &
        OPTIONAL                       :: error

   ! LOCALS
   TYPE(atomic_kind_type), POINTER     :: atom_kind
   TYPE(atomic_kind_type), DIMENSION(:),&
        POINTER                        :: atomic_kind_set
   TYPE(xc_derivative_set_type), &
        POINTER                        :: deriv_set
   TYPE(dft_control_type), POINTER     :: dft_control
   TYPE(qs_energy_type), POINTER       :: energy
   TYPE(grid_atom_type), POINTER       :: grid_atom
   TYPE(harmonics_atom_type), POINTER  :: harmonics
   TYPE(xc_rho_cflags_type)            :: needs
   TYPE(cp_para_env_type), POINTER     :: para_env
   TYPE(rho_atom_type), DIMENSION(:), &
        POINTER                        :: rho_atom_set
   TYPE(rho_atom_coeff),DIMENSION(:),& 
        POINTER                        :: r_h, r_s, dr_h, dr_s
   TYPE(rho_atom_coeff),DIMENSION(:,:),&
        POINTER                        :: r_h_d, r_s_d
   TYPE(rho_atom_type),  POINTER       :: rho_atom
   TYPE(xc_rho_set_type), POINTER      :: rho_set
   TYPE(xc_control_type), POINTER      :: xc_control
   INTEGER                             :: i, iso, stat, asize
   INTEGER                             :: atom, iatom, natom
   INTEGER                             :: ikind, nkind
   INTEGER                             :: ispin, nspins
   INTEGER                             :: ia, na, ir, nr
   INTEGER, DIMENSION(:), POINTER      :: atom_list
   INTEGER, DIMENSION(2)               :: local_loop_limit
   INTEGER, DIMENSION(2,3)             :: bounds
   LOGICAL                             :: cross_term, gradient_functional, failure
   LOGICAL                             :: lsd, paw_atom
   REAL(KIND=dp)                       :: rtot, rad, urad, usin_ia, exc
   REAL(KIND=dp), DIMENSION(:), POINTER:: weight 
   REAL(KIND=dp), DIMENSION(:,:), &
        POINTER                        :: rho_h, rho_s, oodrho_h, oodrho_s, &
                                          vxgc_h, vxgc_s
   REAL(KIND=dp), DIMENSION(:,:,:,:), &
        POINTER                        :: drho_h, drho_s
   REAL(KIND=dp), DIMENSION(:,:,:), &
        POINTER                        :: vxc_h, vxc_s, vxg_h, vxg_s

   ! PARAMETERS
   CHARACTER(LEN=*), PARAMETER         :: routineN = "calculate_xc_2nd_deriv_atom", &
                                          routineP = module_name//"/"//routineN

   failure = .FALSE.
   NULLIFY(atomic_kind_set, dft_control, energy, para_env, xc_control)
   NULLIFY(rho_h, rho_s, drho_h, drho_s, oodrho_h, oodrho_s, weight)
   NULLIFY(vxgc_h, vxgc_s, vxc_h, vxc_s, vxg_h, vxg_s)
   CALL get_qs_env(qs_env=qs_env,&
                   dft_control=dft_control,&
                   para_env=para_env,&
                   energy=energy,&
                   atomic_kind_set=atomic_kind_set)
   rho_atom_set => p_env%local_rho_set%rho_atom_set
   xc_control => dft_control%tddfpt_control%xc_control

   nspins = dft_control%nspins
   gradient_functional = ANY(xc_control%gradient_functionals)
   lsd = (nspins==2)
   ! TDDFT needs an unrestricted 2nd derivative in any case
   IF (dft_control%do_tddfpt_calculation) lsd = .TRUE. 
   nkind = SIZE(atomic_kind_set)

!  Here starts the loop over all the atoms
   DO ikind = 1,nkind
      
      NULLIFY(atom_kind, atom_list, harmonics, grid_atom)
      atom_kind => atomic_kind_set(ikind)
      CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                           natom=natom,paw_atom=paw_atom,&
                           harmonics=harmonics,grid_atom=grid_atom)
      IF (.NOT.paw_atom) CYCLE
      
      nr = grid_atom%nr
      na = grid_atom%ng_sphere

      ! Array dimension: here anly one dimensional arrays are used,
      ! i.e. only the first column of deriv_data is read.
      ! The other to dimensions  are set to size equal 1.
      bounds(1:2,1:3) = 1
      bounds(2,1) = na

      NULLIFY(deriv_set, rho_set)
      CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)
      CALL xc_rho_set_create(rho_set,bounds,rho_cutoff=xc_control%density_cut,&
                             drho_cutoff=xc_control%gradient_cut,tau_cutoff=xc_control%tau_cut,&
                             error=error)

      ! assesses the argument needed by the functionals
      CALL xc_rho_cflags_setall(needs,.FALSE.,error=error)
      DO i=1,SIZE(xc_control%functionals)
         CALL xc_functional_get_info(xc_control%functionals(i), lsd=lsd, &
                                     needs=needs, error=error)
      END DO  
      ! basic components
      IF (lsd) THEN
         needs%rho_spin=.TRUE.
         needs%tau_spin=needs%tau_spin.OR.needs%tau
      ELSE
         needs%rho=.TRUE.
      END IF
      IF (needs%drho .OR. needs%drho_spin .OR. &
          needs%norm_drho .OR. needs%norm_drho_spin) THEN
         IF (lsd) THEN
            needs%drho_spin = .TRUE.
         ELSE
            needs%drho = .TRUE.
         END IF
      END IF

      ! allocate the required 3d arrays where to store rho and drho
      IF (.NOT.dft_control%do_tddfpt_calculation) THEN
         CALL xc_rho_set_atom_update(rho_set, needs, nspins, na)
      ELSE
         CALL xc_rho_set_atom_update(rho_set, needs, 2, na)
      END IF

      !! Will I need to change this??
      cross_term = .FALSE.
      IF (lsd) THEN
         cross_term = needs%norm_drho .OR. needs%drhoa_drhob
      END IF
      
      ALLOCATE(rho_h(1:na,1:nspins), &
               rho_s(1:na,1:nspins), stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(weight(1:na), stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      IF (gradient_functional) THEN
         ALLOCATE(drho_h(1:4,1:na,1:nr,1:nspins), &
                  drho_s(1:4,1:na,1:nr,1:nspins), stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(oodrho_h(1:na,1:nspins), &
                  oodrho_s(1:na,1:nspins), stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      END IF

      DO iatom = 1,natom
         atom = atom_list(iatom)
         DO ispin = 1,nspins
            rho_atom_set(atom)%ga_Vxc_gb_h(ispin)%r_coef = 0.0_dp
            rho_atom_set(atom)%ga_Vxc_gb_s(ispin)%r_coef = 0.0_dp
         END DO
         rho_atom_set(atom)%exc_h = 0.0_dp
         rho_atom_set(atom)%exc_s = 0.0_dp
      END DO

      ! parallelization
      local_loop_limit = get_limit( natom, para_env%num_pe, para_env%mepos )

      DO iatom = local_loop_limit(1), local_loop_limit(2) !1,natom
         atom = atom_list(iatom)
         
         rho_atom => rho_atom_set(atom)
         NULLIFY(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
         IF(gradient_functional) THEN
            CALL get_rho_atom(rho_atom=rho_atom, &
                              rho_rad_h=r_h, rho_rad_s=r_s, &
                              drho_rad_h=dr_h, drho_rad_s=dr_s, &
                              rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
            asize = SIZE(drho_h,1)*SIZE(drho_h,2)*SIZE(drho_h,3)*SIZE(drho_h,4)
            CALL DCOPY(asize, 0.0_dp, 0, drho_h, 1)
            CALL DCOPY(asize, 0.0_dp, 0, drho_s, 1)
         ELSE
            CALL get_rho_atom(rho_atom=rho_atom, &
                              rho_rad_h=r_h, rho_rad_s=r_s)
         END IF

      END DO

      rtot = 0.0_dp
      
      DO ir = 1,nr

         CALL calc_rho_angular(grid_atom, harmonics, nspins, &
                     gradient_functional, xc_control%gradient_cut, &
                     ir, r_h, r_s, rho_h, rho_s, &
                     dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, oodrho_h, oodrho_s, &
                     error)
         weight(1:na) = grid_atom%wr(ir)*grid_atom%wa(1:na)

         ! hard atom density 
         CALL xc_2nd_deriv_of_r(xc_control = xc_control,&
                                rho_set=rho_set,deriv_set=deriv_set,deriv_order=2,needs=needs,&
                                rho=rho_h,drho=drho_h,oodrho=oodrho_h,w=weight,&
                                lsd=lsd,nspins=nspins,na=na,ir_pnt=ir,&
                                exc=exc,vxc=vxc_h,vxg=vxg_h,vxg_cross=vxgc_h,error=error)

!        *** Sum up the exc contribution  to exc_h ***
         rho_atom%exc_h = rho_atom%exc_h + exc

         ! atom density
         CALL xc_2nd_deriv_of_r(xc_control=xc_control,&
                                rho_set=rho_set,deriv_set=deriv_set,deriv_order=2,needs=needs,&
                                rho=rho_s,drho=drho_s,oodrho=oodrho_s,w=weight,&
                                lsd=lsd,nspins=nspins,na=na,ir_pnt=ir,&
                                exc=exc,vxc=vxc_s,vxg=vxg_s,vxg_cross=vxgc_s,error=error)

!        *** Sum up the exc contribution  to exc_s ***
         rho_atom%exc_s = rho_atom%exc_s + exc

      END DO

      ! some cleanup
      DEALLOCATE(rho_h, rho_s, weight, stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      NULLIFY(rho_h, rho_s, weight)
      IF (gradient_functional) THEN
         DEALLOCATE(drho_h, drho_s, oodrho_h, oodrho_s, stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         NULLIFY(drho_h, drho_s, oodrho_h, oodrho_s)
      END IF
   END DO

 END SUBROUTINE calculate_xc_2nd_deriv_atom

! *****************************************************************************

 SUBROUTINE  gaVxcgb_noGC(vxc_h,vxc_s,atom_kind,rho_atom,nspins)


    TYPE(atomic_kind_type), POINTER       :: atom_kind
    TYPE(rho_atom_type),  POINTER         :: rho_atom

    REAL(dp), DIMENSION(:,:,:),POINTER    :: vxc_h,vxc_s
    INTEGER, INTENT(IN)                   :: nspins

    CHARACTER(LEN=*), PARAMETER           :: routine_name = "gaVxcgb_noGC"

    TYPE(gto_basis_set_type), POINTER     :: orb_basis

    TYPE(grid_atom_type), POINTER         :: grid_atom
    TYPE(harmonics_atom_type), &
                           POINTER        :: harmonics
    TYPE(rho_atom_coeff),  DIMENSION(:), &
                           POINTER        :: int_h, int_s

    REAL(dp), DIMENSION(:), ALLOCATABLE   :: g1,g2
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: gg,&
                                             gVg_h, gVg_s,&
                                             matso_h, matso_s
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vx
    REAL(dp), DIMENSION(:,:), POINTER     :: zet
    REAL(dp), DIMENSION(:,:,:), POINTER   :: my_CG
    INTEGER, DIMENSION(:), POINTER        :: lmin,lmax,npgf

    INTEGER                               :: ia,ic,icg,&
                                             ipgf1,ipgf2,ir,iset1,&
                                             iset2,iso,iso1,iso2,&
                                             ispin,istat,l,ld,&
                                             lmin12,lmax12,maxso, maxl, m1,&
                                             m2,n1,n2,na, ngau1,ngau2,nngau1,nr, nset,&
                                             size1

! -----------------------------------------------------------------------------

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

    ALLOCATE(vx(na,nr), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                                                     "vx",na*nr*dp_size)

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
              DO ir =1,nr
                vx(1:na,ir) = vxc_h(1:na,ir,ispin)
              END DO
              CALL dgemm('N','N',na,ld,nr,1.0_dp,vx(1:na,1:nr),na,&
                          gg(1:nr,0:lmax12),nr,0.0_dp,gVg_h(1:na,0:lmax12),na)
              DO ir =1,nr
                vx(1:na,ir) = vxc_s(1:na,ir,ispin)
              END DO
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

   DEALLOCATE(g1,g2,gg,matso_h,matso_s,gVg_s,gVg_h,vx,STAT=istat)
   IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,&
                       "g1,g2,gg,matso_h,matso_s,gVg_s,gVg_h")
 END SUBROUTINE gaVxcgb_noGC

! *****************************************************************************

 SUBROUTINE  gaVxcgb_GC(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,rho_atom,drho_h,&
                        drho_s,nspins,vxgc_h,vxgc_s)


    TYPE(atomic_kind_type), POINTER       :: atom_kind
    TYPE(rho_atom_type),  POINTER         :: rho_atom

    REAL(dp), DIMENSION(:,:,:),POINTER    :: vxc_h,vxc_s,&
                                             vxg_h,vxg_s
    REAL(dp), DIMENSION(:,:,:,:), &
                                POINTER   :: drho_h,drho_s
    INTEGER, INTENT(IN)                   :: nspins
    REAL(dp), DIMENSION(:,:), POINTER,&
              OPTIONAL                    :: vxgc_h, vxgc_s

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
    LOGICAL                               :: cross_term

! -----------------------------------------------------------------------------

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

    cross_term = .FALSE.
    IF(PRESENT(vxgc_h)) cross_term = .TRUE.

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

              IF(cross_term) THEN

!             Cross Term
              DO l = lmin12,lmax12
                DO ia = 1,na
                  DO ir = 1,nr
                     gVXCg_h(ia,l) = gVXCg_h(ia,l) + &
                            gg(ir,l) * vxc_h(ia,ir,ispin)+&
                            dgg(ir,l)*vxg_h(ia,ir,ispin)*&
                           (drho_h(1,ia,ir,ispin)*harmonics%a(1,ia)+&
                            drho_h(2,ia,ir,ispin)*harmonics%a(2,ia)+&
                            drho_h(3,ia,ir,ispin)*harmonics%a(3,ia))+&
                            dgg(ir,l)*vxgc_h(ia,ir)*&
                           ((drho_h(1,ia,ir,1)+drho_h(1,ia,ir,2))*harmonics%a(1,ia)+&
                            (drho_h(2,ia,ir,1)+drho_h(2,ia,ir,2))*harmonics%a(2,ia)+&
                            (drho_h(3,ia,ir,1)+drho_h(3,ia,ir,2))*harmonics%a(3,ia))

                     gVXCg_s(ia,l) = gVXCg_s(ia,l) + &
                            gg(ir,l) * vxc_s(ia,ir,ispin)+&
                            dgg(ir,l)*vxg_s(ia,ir,ispin)*&
                           (drho_s(1,ia,ir,ispin)*harmonics%a(1,ia)+&
                            drho_s(2,ia,ir,ispin)*harmonics%a(2,ia)+&
                            drho_s(3,ia,ir,ispin)*harmonics%a(3,ia))+&
                            dgg(ir,l)*vxgc_s(ia,ir)*&
                           ((drho_s(1,ia,ir,1)+drho_s(1,ia,ir,2))*harmonics%a(1,ia)+&
                            (drho_s(2,ia,ir,1)+drho_s(2,ia,ir,2))*harmonics%a(2,ia)+&
                            (drho_s(3,ia,ir,1)+drho_s(3,ia,ir,2))*harmonics%a(3,ia))

                     urad = grid_atom%oorad2l(ir,1)


                     gVXGg_h(1,ia,l) = gVXGg_h(1,ia,l) + &
                                    drho_h(1,ia,ir,ispin)*vxg_h(ia,ir,ispin)*&
                                    gg(ir,l)*urad + &
                                    (drho_h(1,ia,ir,1)+drho_h(1,ia,ir,2))*&
                                    vxgc_h(ia,ir)*gg(ir,l)*urad
                     gVXGg_h(2,ia,l) = gVXGg_h(2,ia,l) + &
                                    drho_h(2,ia,ir,ispin)*vxg_h(ia,ir,ispin)*&
                                    gg(ir,l)*urad + &
                                    (drho_h(2,ia,ir,1)+drho_h(2,ia,ir,2))*&
                                    vxgc_h(ia,ir)*gg(ir,l)*urad
                     gVXGg_h(3,ia,l) = gVXGg_h(3,ia,l) + &
                                    drho_h(3,ia,ir,ispin)*vxg_h(ia,ir,ispin)*&
                                    gg(ir,l)*urad + &
                                    (drho_h(3,ia,ir,1)+drho_h(3,ia,ir,2))*&
                                    vxgc_h(ia,ir)*gg(ir,l)*urad

                     gVXGg_s(1,ia,l) = gVXGg_s(1,ia,l) + &
                                    drho_s(1,ia,ir,ispin)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*urad + &
                                    (drho_s(1,ia,ir,1)+drho_s(1,ia,ir,2))*&
                                    vxgc_s(ia,ir)*gg(ir,l)*urad
                     gVXGg_s(2,ia,l) = gVXGg_s(2,ia,l) + &
                                    drho_s(2,ia,ir,ispin)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*urad + &
                                    (drho_s(2,ia,ir,1)+drho_s(2,ia,ir,2))*&
                                    vxgc_s(ia,ir)*gg(ir,l)*urad
                     gVXGg_s(3,ia,l) = gVXGg_s(3,ia,l) + &
                                    drho_s(3,ia,ir,ispin)*vxg_s(ia,ir,ispin)*&
                                    gg(ir,l)*urad + &
                                    (drho_s(3,ia,ir,1)+drho_s(3,ia,ir,2))*&
                                    vxgc_s(ia,ir)*gg(ir,l)*urad
                  ENDDO ! ir
                ENDDO ! ia
              ENDDO  ! l

              ELSE

!             Without Cross Term
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
              END IF

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

