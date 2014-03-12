!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the integrals of the Vxc potential calculated
!>      for the atomic density in the basis set of spherical primitives
! *****************************************************************************
MODULE qs_vxc_atom

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                             gto_basis_set_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: tddfpt_excitations,&
                                             tddfpt_triplet,&
                                             xc_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: indso,&
                                             nsoset
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_harmonics_atom,               ONLY: get_none0_cg_list,&
                                             harmonics_atom_type
  USE qs_linres_types,                 ONLY: nablavks_atom_type
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                             rho_atom_coeff,&
                                             rho_atom_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  USE xc_atom,                         ONLY: fill_rho_set,&
                                             vxc_of_r_new,&
                                             xc_2nd_deriv_of_r,&
                                             xc_rho_set_atom_update
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_release,&
                                             xc_dset_zero_all
  USE xc_derivatives,                  ONLY: xc_functionals_get_needs
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vxc_atom'

  PUBLIC :: calculate_vxc_atom,&
            calculate_xc_2nd_deriv_atom,&
            calc_rho_angular,&
            gaVxcgb_noGC

CONTAINS

! *****************************************************************************
  SUBROUTINE calculate_vxc_atom(qs_env,energy_only,error,gradient_atom_set,adiabatic_rescale_factor)

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: energy_only
    TYPE(cp_error_type), INTENT(inout)       :: error
    TYPE(nablavks_atom_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: gradient_atom_set
    REAL(dp), INTENT(IN), OPTIONAL           :: adiabatic_rescale_factor

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_vxc_atom', &
      routineP = moduleN//':'//routineN

    INTEGER :: bo(2), handle, ia, iat, iatom, idir, ikind, ir, ispin, istat, &
      mepos, myfun, na, natom, nkind, nr, nspins, num_pe
    INTEGER, DIMENSION(2, 3)                 :: bounds
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: epr_xc, failure, gradient_f, &
                                                lsd, paw_atom, tau_f
    REAL(dp)                                 :: density_cut, exc_h, exc_s, &
                                                gradient_cut, &
                                                my_adiabatic_rescale_factor, &
                                                tau_cut
    REAL(dp), DIMENSION(:, :), POINTER       :: rho_h, rho_s, tau_h, tau_s, &
                                                weight
    REAL(dp), DIMENSION(:, :, :), POINTER    :: vtau_h, vtau_s, vxc_h, vxc_s
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: drho_h, drho_s, vxg_h, vxg_s
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: dr_h, dr_s, r_h, r_s
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      POINTER                                :: r_h_d, r_s_d, trho_h, trho_s
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(section_vals_type), POINTER         :: input, xc_fun_section, &
                                                xc_section
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho_set_h, rho_set_s

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (atom_kind)
    NULLIFY (atom_list)
    NULLIFY (atomic_kind_set)
    NULLIFY (deriv_set)
    NULLIFY (energy)
    NULLIFY (grid_atom)
    NULLIFY (harmonics)
    NULLIFY (input)
    NULLIFY (para_env)
    NULLIFY (rho_atom)
    NULLIFY (rho_atom_set)
    NULLIFY (rho_set_h)
    NULLIFY (rho_set_s)

    epr_xc = .FALSE.
    IF(PRESENT(gradient_atom_set)) THEN
       epr_xc = .TRUE.
    END IF

    IF(PRESENT(adiabatic_rescale_factor)) THEN
      my_adiabatic_rescale_factor = adiabatic_rescale_factor
    ELSE
      my_adiabatic_rescale_factor = 1.0_dp
    END IF

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    energy=energy,&
                    atomic_kind_set=atomic_kind_set,&
                    input=input,&
                    rho_atom_set=rho_atom_set,error=error)

    IF(epr_xc) THEN
       xc_section => section_vals_get_subs_vals(input,&
                     "PROPERTIES%LINRES%EPR%PRINT%G_TENSOR%XC",error=error)
    ELSE
       xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
    END IF
    xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",&
         error=error)
    CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",&
                               i_val=myfun,error=error)

    IF(myfun == xc_none) THEN
      energy%exc1 = 0.0_dp
      rho_atom_set(:)%exc_h = 0.0_dp
      rho_atom_set(:)%exc_s = 0.0_dp
    ELSE
      CALL section_vals_val_get(xc_section,"DENSITY_CUTOFF",&
           r_val=density_cut,error=error)
      CALL section_vals_val_get(xc_section,"GRADIENT_CUTOFF",&
           r_val=gradient_cut,error=error)
      CALL section_vals_val_get(xc_section,"TAU_CUTOFF",&
           r_val=tau_cut,error=error)

      lsd = dft_control%lsd
      nspins = dft_control%nspins
      needs = xc_functionals_get_needs(xc_fun_section,&
                                       lsd=lsd,&
                                       add_basic_components=.TRUE.,&
                                       error=error)

      ! whatever the xc, if epr_xc, drho_spin is needed
      IF(epr_xc) needs%drho_spin = .TRUE.

      gradient_f = (needs%drho.OR.needs%drho_spin)
      tau_f = (needs%tau.OR.needs%tau_spin)

      nkind = SIZE(atomic_kind_set)

      ! Initialize energy contribution from the one center XC terms to zero
      energy%exc1 = 0.0_dp

      ! Nullify some pointers for work-arrays
      NULLIFY (rho_h,drho_h,rho_s,drho_s,weight)
      NULLIFY (vxc_h,vxc_s,vxg_h,vxg_s)
      NULLIFY (tau_h,tau_s)
      NULLIFY (vtau_h,vtau_s)

      ! Here starts the loop over all the atoms

      DO ikind = 1,nkind

        atom_kind => atomic_kind_set(ikind)
        CALL get_atomic_kind(atomic_kind=atom_kind,atom_list=atom_list,&
                             natom=natom,paw_atom=paw_atom,&
                             harmonics=harmonics,grid_atom=grid_atom)

        IF (.NOT.paw_atom) CYCLE

        nr = grid_atom%nr
        na = grid_atom%ng_sphere

        ! Prepare the structures needed to calculate and store the xc derivatives

        ! Array dimension: here anly one dimensional arrays are used,
        ! i.e. only the first column of deriv_data is read.
        ! The other to dimensions  are set to size equal 1
        bounds(1:2,1:3) = 1
        bounds(2,1) = na
        bounds(2,2) = nr

        ! create a place where to put the derivatives
        CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)
        ! create the place where to store the argument for the functionals
        CALL xc_rho_set_create(rho_set_h,bounds,rho_cutoff=density_cut,&
                               drho_cutoff=gradient_cut,tau_cutoff=tau_cut,&
                               error=error)
        CALL xc_rho_set_create(rho_set_s,bounds,rho_cutoff=density_cut,&
                               drho_cutoff=gradient_cut,tau_cutoff=tau_cut,&
                               error=error)

        ! allocate the required 3d arrays where to store rho and drho
        CALL xc_rho_set_atom_update(rho_set_h,needs,nspins,bounds)
        CALL xc_rho_set_atom_update(rho_set_s,needs,nspins,bounds)

        CALL reallocate(rho_h,1,na,1,nspins)
        CALL reallocate(rho_s,1,na,1,nspins)
        weight => grid_atom%weight
        CALL reallocate(vxc_h,1,na,1,nr,1,nspins)
        CALL reallocate(vxc_s,1,na,1,nr,1,nspins)

        IF (gradient_f) THEN
          CALL reallocate(drho_h,1,4,1,na,1,nr,1,nspins)
          CALL reallocate(drho_s,1,4,1,na,1,nr,1,nspins)
          CALL reallocate(vxg_h,1,3,1,na,1,nr,1,nspins)
          CALL reallocate(vxg_s,1,3,1,na,1,nr,1,nspins)
        END IF

        IF (tau_f) THEN
          CALL reallocate(tau_h,1,na,1,nspins)
          CALL reallocate(tau_s,1,na,1,nspins)
          CALL reallocate(vtau_h,1,na,1,nr,1,nspins)
          CALL reallocate(vtau_s,1,na,1,nr,1,nspins)
        END IF

        ! Distribute the atoms of this kind

        num_pe = para_env%num_pe
        mepos  = para_env%mepos
        bo = get_limit( natom, num_pe, mepos )

        DO iat = bo(1),bo(2)
          iatom = atom_list(iat)

          rho_atom_set(iatom)%exc_h = 0.0_dp
          rho_atom_set(iatom)%exc_s = 0.0_dp

          rho_atom => rho_atom_set(iatom)
          IF (gradient_f) THEN
            NULLIFY(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
            CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,&
                              rho_rad_s=r_s,drho_rad_h=dr_h,&
                              drho_rad_s=dr_s,rho_rad_h_d=r_h_d,&
                              rho_rad_s_d=r_s_d)
            drho_h = 0.0_dp
            drho_s = 0.0_dp
          ELSE
            NULLIFY(r_h,r_s)
            CALL get_rho_atom(rho_atom=rho_atom,rho_rad_h=r_h,rho_rad_s=r_s)
          END IF
          IF (tau_f) THEN
            NULLIFY(trho_h,trho_s)
            CALL get_rho_atom(rho_atom=rho_atom,trho_rad_h=trho_h,trho_rad_s=trho_s)
          END IF

          DO ir = 1,nr
            CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f,&
                                  ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, &
                                  r_h_d, r_s_d, drho_h, drho_s, error=error)
            IF (tau_f) THEN
              CALL calc_tau_angular(grid_atom, harmonics, nspins, ir, &
                                  trho_h, trho_s, tau_h, tau_s, error=error)
            END IF
            CALL fill_rho_set(rho_set_h,lsd,nspins,needs,rho_h,drho_h,tau_h,na,ir,error=error)
            CALL fill_rho_set(rho_set_s,lsd,nspins,needs,rho_s,drho_s,tau_s,na,ir,error=error)
          END DO

          !-------------------!
          ! hard atom density !
          !-------------------!
          CALL xc_dset_zero_all(deriv_set, error)
          CALL vxc_of_r_new(xc_fun_section, rho_set_h, deriv_set, 1, needs, weight, &
                            lsd, na, nr, exc_h, vxc_h, vxg_h, vtau_h, energy_only=energy_only, &
                            epr_xc=epr_xc, adiabatic_rescale_factor=my_adiabatic_rescale_factor, error=error)
          rho_atom%exc_h = rho_atom%exc_h + exc_h

          !-------------------!
          ! soft atom density !
          !-------------------!
          CALL xc_dset_zero_all(deriv_set, error)
          CALL vxc_of_r_new(xc_fun_section, rho_set_s, deriv_set, 1, needs, weight, &
                            lsd, na, nr, exc_s, vxc_s, vxg_s,  vtau_s, energy_only=energy_only, &
                            epr_xc=epr_xc, adiabatic_rescale_factor=my_adiabatic_rescale_factor, error=error)
          rho_atom%exc_s = rho_atom%exc_s + exc_s

          IF (epr_xc) THEN
             DO ispin = 1,nspins
                DO idir = 1,3
                   DO ir = 1,nr
                      DO ia = 1,na
                         gradient_atom_set(iatom)%nablavks_vec_rad_h(idir,ispin)%r_coef(ir,ia) = &
                            gradient_atom_set(iatom)%nablavks_vec_rad_h(idir,ispin)%r_coef(ir,ia) &
                            + vxg_h(idir,ia,ir,ispin)
                         gradient_atom_set(iatom)%nablavks_vec_rad_s(idir,ispin)%r_coef(ir,ia) = &
                            gradient_atom_set(iatom)%nablavks_vec_rad_s(idir,ispin)%r_coef(ir,ia) &
                            + vxg_s(idir,ia,ir,ispin)
                      END DO ! ia
                   END DO ! ir
                END DO ! idir
             END DO ! ispin
          END IF

          ! Add contributions to the exc energy

          energy%exc1 = energy%exc1 + rho_atom%exc_h - rho_atom%exc_s

          ! Integration to get the matrix elements relative to the vxc_atom
          ! here the products with the primitives is done: gaVxcgb
          ! internal transformation to get the integral in cartesian Gaussians

          IF (.NOT.energy_only) THEN
            atom_kind => atomic_kind_set(ikind)
            IF (gradient_f) THEN
              CALL gaVxcgb_GC(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,&
                              rho_atom,drho_h,drho_s,nspins,error=error)
            ELSE
              CALL gaVxcgb_noGC(vxc_h,vxc_s,atom_kind,rho_atom,nspins,error=error)
            END IF
            IF (tau_f) THEN
              CALL dgaVtaudgb(vtau_h,vtau_s,atom_kind,rho_atom,nspins,error=error)
            END IF
          END IF  ! energy_only
          NULLIFY (r_h,r_s,dr_h,dr_s)
        END DO ! iat

        ! Release the xc structure used to store the xc derivatives
        CALL xc_dset_release(deriv_set, error=error)
        CALL xc_rho_set_release(rho_set_h,error=error)
        CALL xc_rho_set_release(rho_set_s,error=error)

      END DO ! ikind

      CALL mp_sum(energy%exc1,para_env%group)

      IF (ASSOCIATED(rho_h)) THEN
        DEALLOCATE (rho_h,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(rho_s)) THEN
        DEALLOCATE (rho_s,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(vxc_h)) THEN
        DEALLOCATE (vxc_h,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (ASSOCIATED(vxc_s)) THEN
        DEALLOCATE (vxc_s,STAT=istat)
        CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
      END IF

      IF (gradient_f) THEN

        IF (ASSOCIATED(drho_h)) THEN
          DEALLOCATE (drho_h,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF (ASSOCIATED(drho_s)) THEN
          DEALLOCATE (drho_s,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF (ASSOCIATED(vxg_h)) THEN
          DEALLOCATE (vxg_h,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF (ASSOCIATED(vxg_s)) THEN
          DEALLOCATE (vxg_s,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

      END IF

      IF (tau_f) THEN

        IF (ASSOCIATED(tau_h)) THEN
          DEALLOCATE (tau_h,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF (ASSOCIATED(tau_s)) THEN
          DEALLOCATE (tau_s,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF (ASSOCIATED(vtau_h)) THEN
          DEALLOCATE (vtau_h,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

        IF (ASSOCIATED(vtau_s)) THEN
          DEALLOCATE (vtau_s,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
        END IF

      END IF

    END IF !xc_none

    CALL timestop(handle)

  END SUBROUTINE calculate_vxc_atom

! *****************************************************************************
  SUBROUTINE calculate_xc_2nd_deriv_atom(p_env,qs_env,xc_section,do_tddft,do_triplet,error)

    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: xc_section
    LOGICAL, INTENT(IN), OPTIONAL            :: do_tddft, do_triplet
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_xc_2nd_deriv_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: atom, excitations, handle, &
                                                iatom, ikind, ir, na, natom, &
                                                nkind, nr, nspins, res_etype, &
                                                stat
    INTEGER, DIMENSION(2)                    :: local_loop_limit
    INTEGER, DIMENSION(2, 3)                 :: bounds
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, gradient_functional, &
                                                lsd, lsd_singlets, &
                                                my_do_tddft, paw_atom, tau_f
    REAL(KIND=dp)                            :: density_cut, exc, &
                                                gradient_cut, rtot, tau_cut
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rho1_h, rho1_s, rho_h, rho_s, &
                                                tau1_h, tau1_s, tau_h, tau_s, &
                                                weight
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: vxc_h, vxc_s
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: drho1_h, drho1_s, drho_h, &
                                                drho_s, vxg_h, vxg_s
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: dr1_h, dr1_s, dr_h, dr_s, &
                                                r1_h, r1_s, r_h, r_s
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      POINTER                                :: r1_h_d, r1_s_d, r_h_d, r_s_d
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho1_atom_set, rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho1_atom, rho_atom
    TYPE(section_vals_type), POINTER         :: input, xc_fun_section
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho1_set_h, rho1_set_s, &
                                                rho_set_h, rho_set_s

! -------------------------------------------------------------------------

   CALL timeset(routineN,handle)

   failure = .FALSE.

   NULLIFY(atomic_kind_set, energy, para_env)
   NULLIFY(rho_h, rho_s, drho_h, drho_s, weight)
   NULLIFY(rho1_h, rho1_s, drho1_h, drho1_s)
   NULLIFY(vxc_h, vxc_s, vxg_h, vxg_s)
   NULLIFY(tau_h, tau_s, tau1_h, tau1_s)

   my_do_tddft = .TRUE.
   IF(PRESENT(do_tddft)) my_do_tddft = do_tddft

   CALL get_qs_env(qs_env=qs_env,&
                   para_env=para_env,&
                   energy=energy,&
                   input=input,&
                   atomic_kind_set=atomic_kind_set,&
                   rho_atom_set=rho_atom_set,&
                   error=error)

   rho1_atom_set => p_env%local_rho_set%rho_atom_set

    CALL section_vals_val_get(input,"DFT%LSD",l_val=lsd,error=error)
    CALL section_vals_val_get(xc_section,"DENSITY_CUTOFF",&
         r_val=density_cut,error=error)
    CALL section_vals_val_get(xc_section,"GRADIENT_CUTOFF",&
         r_val=gradient_cut,error=error)
    CALL section_vals_val_get(xc_section,"TAU_CUTOFF",&
         r_val=tau_cut,error=error)
    IF(my_do_tddft) THEN
       CALL section_vals_val_get(input,"DFT%EXCITATIONS",&
            i_val=excitations,error=error)
       CALL section_vals_val_get(input,"DFT%TDDFPT%LSD_SINGLETS",&
            l_val=lsd_singlets,error=error)
       CALL section_vals_val_get(input,"DFT%TDDFPT%RES_ETYPE",&
            i_val=res_etype,error=error)
    ENDIF
    xc_fun_section => section_vals_get_subs_vals(xc_section,&
         "XC_FUNCTIONAL",error=error)
    IF (lsd) THEN
       nspins=2
    ELSE
       nspins=1
    END IF
    needs=xc_functionals_get_needs(xc_fun_section,lsd=lsd,&
         add_basic_components=.TRUE., error=error)
    gradient_functional = needs%drho .OR. needs%drho_spin
    tau_f = (needs%tau.OR.needs%tau_spin)
    IF ( tau_f ) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
              message="Tau functionals not implemented for GAPW 2nd derivatives", &
              error=error, error_level=cp_failure_level)
    END IF

    IF(my_do_tddft) THEN
       IF (excitations==tddfpt_excitations) THEN
          IF (nspins == 1 .AND. (lsd_singlets .OR. res_etype == tddfpt_triplet)) THEN
             lsd = .TRUE.
          END IF
       END IF
    ELSE
       IF(PRESENT(do_triplet)) THEN
          IF(nspins == 1 .AND. do_triplet) THEN
             lsd = .TRUE.
          ENDIF
       ENDIF
    ENDIF

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
      bounds(2,2) = nr

      NULLIFY(deriv_set, rho_set_h, rho_set_s, rho1_set_h, rho1_set_s)
      CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)
      CALL xc_rho_set_create(rho_set_h,bounds,rho_cutoff=density_cut,&
                             drho_cutoff=gradient_cut,tau_cutoff=tau_cut,&
                             error=error)
      CALL xc_rho_set_create(rho_set_s,bounds,rho_cutoff=density_cut,&
                             drho_cutoff=gradient_cut,tau_cutoff=tau_cut,&
                             error=error)
      CALL xc_rho_set_create(rho1_set_h,bounds,rho_cutoff=density_cut,&
                             drho_cutoff=gradient_cut,tau_cutoff=tau_cut,&
                             error=error)
      CALL xc_rho_set_create(rho1_set_s,bounds,rho_cutoff=density_cut,&
                             drho_cutoff=gradient_cut,tau_cutoff=tau_cut,&
                             error=error)

      ! allocate the required 3d arrays where to store rho and drho
      IF (nspins == 1 .AND. .NOT. lsd) THEN
         CALL xc_rho_set_atom_update(rho_set_h, needs, 1, bounds)
         CALL xc_rho_set_atom_update(rho1_set_h, needs, 1, bounds)
         CALL xc_rho_set_atom_update(rho_set_s, needs, 1, bounds)
         CALL xc_rho_set_atom_update(rho1_set_s, needs, 1, bounds)
      ELSE
         CALL xc_rho_set_atom_update(rho_set_h, needs, 2, bounds)
         CALL xc_rho_set_atom_update(rho1_set_h, needs, 2, bounds)
         CALL xc_rho_set_atom_update(rho_set_s, needs, 2, bounds)
         CALL xc_rho_set_atom_update(rho1_set_s, needs, 2, bounds)
      END IF

      ALLOCATE (rho_h(1:na,1:nspins),rho1_h(1:na,1:nspins), &
                rho_s(1:na,1:nspins),rho1_s(1:na,1:nspins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ALLOCATE (vxc_h(1:na,1:nr,1:nspins), vxc_s(1:na,1:nr,1:nspins),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      vxc_h = 0.0_dp
      vxc_s = 0.0_dp

      weight => grid_atom%weight

      IF (gradient_functional) THEN
         ALLOCATE(drho_h(1:4,1:na,1:nr,1:nspins), drho1_h(1:4,1:na,1:nr,1:nspins), &
                  drho_s(1:4,1:na,1:nr,1:nspins), drho1_s(1:4,1:na,1:nr,1:nspins), stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

         ALLOCATE(vxg_h(1:3,1:na,1:nr,1:nspins), vxg_s(1:3,1:na,1:nr,1:nspins), stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         vxg_h = 0.0_dp
         vxg_s = 0.0_dp
      END IF

      ! parallelization
      local_loop_limit = get_limit( natom, para_env%num_pe, para_env%mepos )

      DO iatom = local_loop_limit(1), local_loop_limit(2) !1,natom
         atom = atom_list(iatom)

         rho_atom_set(atom)%exc_h = 0.0_dp
         rho_atom_set(atom)%exc_s = 0.0_dp
         rho1_atom_set(atom)%exc_h = 0.0_dp
         rho1_atom_set(atom)%exc_s = 0.0_dp

         rho_atom => rho_atom_set(atom)
         rho1_atom => rho1_atom_set(atom)
         NULLIFY(r_h,r_s,dr_h,dr_s,r_h_d,r_s_d)
         NULLIFY(r1_h,r1_s,dr1_h,dr1_s,r1_h_d,r1_s_d)
         IF(gradient_functional) THEN
            CALL get_rho_atom(rho_atom=rho_atom, &
                              rho_rad_h=r_h, rho_rad_s=r_s, &
                              drho_rad_h=dr_h, drho_rad_s=dr_s, &
                              rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
            CALL get_rho_atom(rho_atom=rho1_atom, &
                              rho_rad_h=r1_h, rho_rad_s=r1_s, &
                              drho_rad_h=dr1_h, drho_rad_s=dr1_s, &
                              rho_rad_h_d=r1_h_d, rho_rad_s_d=r1_s_d)
            drho_h = 0.0_dp ; drho_s = 0.0_dp
            drho1_h = 0.0_dp ; drho1_s = 0.0_dp
         ELSE
            CALL get_rho_atom(rho_atom=rho_atom, &
                              rho_rad_h=r_h, rho_rad_s=r_s)
            CALL get_rho_atom(rho_atom=rho1_atom, &
                              rho_rad_h=r1_h, rho_rad_s=r1_s)
         END IF

         rtot = 0.0_dp

         DO ir = 1,nr

            CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_functional, &
                                  ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, r_h_d, r_s_d, &
                                  drho_h, drho_s, error)
            CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_functional, &
                                  ir, r1_h, r1_s, rho1_h, rho1_s, dr1_h, dr1_s, r1_h_d, r1_s_d, &
                                  drho1_h, drho1_s, error)

            CALL fill_rho_set(rho_set_h,lsd,nspins,needs,rho_h,drho_h,tau_h,na,ir,error)
            CALL fill_rho_set(rho1_set_h,lsd,nspins,needs,rho1_h,drho1_h,tau1_h,na,ir,error)
            CALL fill_rho_set(rho_set_s,lsd,nspins,needs,rho_s,drho_s,tau_s,na,ir,error)
            CALL fill_rho_set(rho1_set_s,lsd,nspins,needs,rho1_s,drho1_s,tau1_s,na,ir,error)

         END DO

         CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                rho_set=rho_set_h, rho1_set=rho1_set_h, &
                                deriv_set=deriv_set, needs=needs, &
                                w=weight, exc=exc, vxc=vxc_h, vxg=vxg_h, error=error)
         CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                rho_set=rho_set_s, rho1_set=rho1_set_s, &
                                deriv_set=deriv_set, needs=needs, &
                                w=weight, exc=exc, vxc=vxc_s, vxg=vxg_s, error=error)

         atom_kind => atomic_kind_set(ikind)

         IF(gradient_functional) THEN
            CALL gaVxcgb_GC(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,&
                            rho1_atom, drho_h, drho_s, nspins,error=error)
         ELSE
            CALL gaVxcgb_noGC(vxc_h, vxc_s, atom_kind, &
                              rho1_atom, nspins,error=error)
         ENDIF

         NULLIFY(r_h,r_s,dr_h,dr_s)

      END DO

      ! some cleanup
      DEALLOCATE(rho_h, rho_s, rho1_h, rho1_s, vxc_h, vxc_s, stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      NULLIFY(rho_h, rho_s, weight)
      IF (gradient_functional) THEN
         DEALLOCATE(drho_h, drho_s, vxg_h, vxg_s, stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(drho1_h, drho1_s, stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         NULLIFY(drho_h, drho_s, drho1_h, drho1_s, vxg_h, vxc_s)
      END IF

      CALL xc_dset_release(deriv_set, error=error)
      CALL xc_rho_set_release(rho_set_h,error=error)
      CALL xc_rho_set_release(rho1_set_h,error=error)
      CALL xc_rho_set_release(rho_set_s,error=error)
      CALL xc_rho_set_release(rho1_set_s,error=error)

   END DO

   CALL timestop(handle)

 END SUBROUTINE calculate_xc_2nd_deriv_atom

! *****************************************************************************
  SUBROUTINE calc_rho_angular(grid_atom, harmonics, nspins, grad_func, &
                              ir, r_h, r_s, rho_h, rho_s, &
                              dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s, error)

    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    INTEGER, INTENT(IN)                      :: nspins
    LOGICAL, INTENT(IN)                      :: grad_func
    INTEGER, INTENT(IN)                      :: ir
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: r_h, r_s
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rho_h, rho_s
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: dr_h, dr_s
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      POINTER                                :: r_h_d, r_s_d
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      POINTER                                :: drho_h, drho_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_rho_angular', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iso, ispin, na
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rad, urad

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
   END IF

   IF (failure) RETURN

   na = grid_atom%ng_sphere
   rad = grid_atom%rad(ir)
   urad = grid_atom%oorad2l(ir,1)
   rho_h = 0.0_dp ; rho_s = 0.0_dp
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
      DO ispin = 1,nspins
         DO iso = 1,harmonics%max_iso_not0
            DO ia = 1,na

!              components of the gradient of rho1 hard
               drho_h(1,ia,ir,ispin) = drho_h(1,ia,ir,ispin) + &
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(1,ia)*harmonics%slm(ia,iso) +&
                                       r_h_d(1,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_h(2,ia,ir,ispin) = drho_h(2,ia,ir,ispin) + &
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(2,ia)*harmonics%slm(ia,iso) +&
                                       r_h_d(2,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_h(3,ia,ir,ispin) = drho_h(3,ia,ir,ispin) + &
                                       dr_h(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(3,ia)*harmonics%slm(ia,iso) +&
                                       r_h_d(3,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

!              components of the gradient of rho1 soft
               drho_s(1,ia,ir,ispin) = drho_s(1,ia,ir,ispin) + &
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(1,ia)*harmonics%slm(ia,iso) +&
                                       r_s_d(1,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_s(2,ia,ir,ispin) = drho_s(2,ia,ir,ispin) + &
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(2,ia)*harmonics%slm(ia,iso) +&
                                       r_s_d(2,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_s(3,ia,ir,ispin) = drho_s(3,ia,ir,ispin) + &
                                       dr_s(ispin)%r_coef(ir,iso)*&
                                       harmonics%a(3,ia)*harmonics%slm(ia,iso) +&
                                      r_s_d(3,ispin)%r_coef(ir,iso)*&
                                       harmonics%slm(ia,iso)

               drho_h(4,ia,ir,ispin) = SQRT(&
                    drho_h(1,ia,ir,ispin)*drho_h(1,ia,ir,ispin)+&
                    drho_h(2,ia,ir,ispin)*drho_h(2,ia,ir,ispin)+&
                    drho_h(3,ia,ir,ispin)*drho_h(3,ia,ir,ispin))

               drho_s(4,ia,ir,ispin) = SQRT(&
                    drho_s(1,ia,ir,ispin)*drho_s(1,ia,ir,ispin)+&
                    drho_s(2,ia,ir,ispin)*drho_s(2,ia,ir,ispin)+&
                    drho_s(3,ia,ir,ispin)*drho_s(3,ia,ir,ispin))

            END DO  ! ia
         END DO  ! iso
      END DO  ! ispin
   END IF

 END SUBROUTINE calc_rho_angular
! *****************************************************************************
  SUBROUTINE calc_tau_angular(grid_atom, harmonics, nspins, ir, &
                              trho_h, trho_s, tau_h, tau_s, error)

    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    INTEGER, INTENT(IN)                      :: nspins, ir
    TYPE(rho_atom_coeff), DIMENSION(:, :), &
      POINTER                                :: trho_h, trho_s
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tau_h, tau_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_tau_angular', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ia, iso, ispin, na
    LOGICAL                                  :: failure

    failure = .FALSE.
    CPPostcondition(ASSOCIATED(trho_h),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(trho_s),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(tau_h),cp_failure_level,routineP,error,failure)
    CPPostcondition(ASSOCIATED(tau_s),cp_failure_level,routineP,error,failure)

    IF (failure) RETURN

    na = grid_atom%ng_sphere
    tau_h = 0.0_dp
    tau_s = 0.0_dp

    DO ispin = 1,nspins
      DO iso = 1,harmonics%max_iso_not0
        DO ia = 1,na
          tau_h(ia,ispin) = tau_h(ia,ispin) + &
                            trho_h(1,ispin)%r_coef(ir,iso)*harmonics%slm(ia,iso)
          tau_h(ia,ispin) = tau_h(ia,ispin) + &
            trho_h(3,ispin)%r_coef(ir,iso)*harmonics%slm(ia,iso)*grid_atom%usin_pol(ia)**2
          tau_s(ia,ispin) = tau_s(ia,ispin) + &
                            trho_s(1,ispin)%r_coef(ir,iso)*harmonics%slm(ia,iso)
          tau_s(ia,ispin) = tau_s(ia,ispin) + &
            trho_s(3,ispin)%r_coef(ir,iso)*harmonics%slm(ia,iso)*grid_atom%usin_pol(ia)**2
        END DO  ! ia
      END DO  ! iso
      DO iso = 1,harmonics%max_iso_not0
        DO ia = 1,na
          tau_h(ia,ispin) = tau_h(ia,ispin) + &
                            trho_h(2,ispin)%r_coef(ir,iso)*harmonics%slm(ia,iso)
          tau_s(ia,ispin) = tau_s(ia,ispin) + &
                            trho_s(2,ispin)%r_coef(ir,iso)*harmonics%slm(ia,iso)
       END DO  ! ia
     END DO  ! iso
   END DO  ! ispin

  END SUBROUTINE calc_tau_angular

! *****************************************************************************
  SUBROUTINE  gaVxcgb_noGC(vxc_h,vxc_s,atom_kind,rho_atom,nspins,error)

    REAL(dp), DIMENSION(:, :, :), POINTER    :: vxc_h, vxc_s
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(rho_atom_type), POINTER             :: rho_atom
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gaVxcgb_noGC', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, ia, ic, icg, ipgf1, ipgf2, ir, iset1, iset2, iso, &
      iso1, iso2, ispin, istat, l, ld, lmax12, lmax_expansion, lmin12, m1, &
      m2, max_iso_not0, max_iso_not0_local, max_s_harm, maxl, maxso, n1, n2, &
      na, ngau1, ngau2, nngau1, nr, nset, size1
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cg_n_list
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cg_list
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: g1, g2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: gg, gVg_h, gVg_s, matso_h, &
                                                matso_s, vx
    REAL(dp), DIMENSION(:, :), POINTER       :: zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_hh, int_ss

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY (lmin,lmax,npgf,zet,my_CG,harmonics,grid_atom)

    CALL get_atomic_kind(atomic_kind=atom_kind, orb_basis_set=orb_basis,&
                         harmonics=harmonics,&
                         grid_atom=grid_atom)

    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxso=maxso,maxl=maxl,npgf=npgf,&
                           nset=nset,zet=zet)

    nr = grid_atom%nr
    na = grid_atom%ng_sphere
    my_CG => harmonics%my_CG
    max_iso_not0 = harmonics%max_iso_not0
    lmax_expansion = indso(1,max_iso_not0)
    max_s_harm = harmonics%max_s_harm

    ALLOCATE (g1(nr),g2(nr),gg(nr,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "g1,g2,gg",(2*nr+2*maxl*nr)*dp_size)

    ALLOCATE(gVg_h(na,0:2*maxl),gVg_s(na,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gVg_h,gVg_s",4*na*maxl*dp_size)

    ALLOCATE(matso_h(nsoset(maxl),nsoset(maxl)),&
             matso_s(nsoset(maxl),nsoset(maxl)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "matso_h,matso_s",&
                                     2*nsoset(maxl)*nsoset(maxl)*dp_size)

    ALLOCATE(vx(na,nr), STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "vx",na*nr*dp_size)

    ALLOCATE(cg_list(2,nsoset(maxl)**2,max_s_harm),cg_n_list(max_s_harm),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    NULLIFY (int_hh,int_ss)

    CALL get_rho_atom(rho_atom=rho_atom,ga_Vlocal_gb_h=int_hh,ga_Vlocal_gb_s=int_ss)

    g1 = 0.0_dp
    g2 = 0.0_dp
    m1 = 0
    DO iset1 = 1,nset
      n1 = nsoset(lmax(iset1))
      m2 = 0
      DO iset2 = 1,nset
         CALL get_none0_cg_list(my_CG,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                max_s_harm,lmax_expansion,cg_list,cg_n_list,max_iso_not0_local,error)
         CPPrecondition(max_iso_not0_local.LE.max_iso_not0,cp_failure_level,routineP,error,failure)

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

            ! reduce expansion local densities
            IF(lmin12 .LE. lmax_expansion) THEN

            gg = 0.0_dp
            IF (lmin12 == 0) THEN
              gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
            ELSE
              gg(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12)*g1(1:nr)*g2(1:nr)
            END IF

            ! limit the expansion of the local densities to a max L
            IF(lmax12 .GT. lmax_expansion) lmax12 = lmax_expansion

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
              DO iso = 1,max_iso_not0_local
                 DO icg = 1,cg_n_list(iso)
                    iso1 = cg_list(1,icg,iso)
                    iso2 = cg_list(2,icg,iso)
                    l = indso(1,iso1) + indso(1,iso2)

                    CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)
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
                CALL daxpy(size1,1.0_dp,matso_h(iso1,ic),1,&
                           int_hh(ispin)%r_coef(nngau1+1,iso2),1)
                CALL daxpy(size1,1.0_dp,matso_s(iso1,ic),1,&
                           int_ss(ispin)%r_coef(nngau1+1,iso2),1)
               END DO

            END DO ! ispin

            END IF ! lmax_expansion

          END DO ! ipfg2
        END DO ! ipfg1
        m2 = m2 + maxso
      END DO ! iset2
      m1 = m1 + maxso
    END DO  ! iset1

   DEALLOCATE(g1,g2,gg,matso_h,matso_s,gVg_s,gVg_h,vx,STAT=istat)
   IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "g1,g2,gg,matso_h,matso_s,gVg_s,gVg_h")

   DEALLOCATE(cg_list,cg_n_list,STAT=istat)
   CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

   CALL timestop(handle)

  END SUBROUTINE gaVxcgb_noGC

! *****************************************************************************
  SUBROUTINE  gaVxcgb_GC(vxc_h,vxc_s,vxg_h,vxg_s,atom_kind,rho_atom,drho_h,&
                         drho_s,nspins,error)

    REAL(dp), DIMENSION(:, :, :), POINTER    :: vxc_h, vxc_s
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: vxg_h, vxg_s
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(rho_atom_type), POINTER             :: rho_atom
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: drho_h, drho_s
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'gaVxcgb_GC', &
      routineP = moduleN//':'//routineN

    INTEGER :: dmax_iso_not0, dmax_iso_not0_local, handle, ia, ic, icg, &
      ipgf1, ipgf2, ir, iset1, iset2, iso, iso1, iso2, ispin, istat, l, &
      lmax12, lmax_expansion, lmin12, m1, m2, max_iso_not0, &
      max_iso_not0_local, max_s_harm, maxl, maxso, n1, n2, na, ngau1, ngau2, &
      nngau1, nr, nset, size1
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cg_n_list, dcg_n_list
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cg_list, dcg_list
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf
    LOGICAL                                  :: failure
    REAL(dp)                                 :: urad
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: g1, g2
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: dgg, gg, gVXCg_h, gVXCg_s, &
                                                matso_h, matso_s
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: gVXGg_h, gVXGg_s
    REAL(dp), DIMENSION(:, :), POINTER       :: zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG
    REAL(dp), DIMENSION(:, :, :, :), POINTER :: my_CG_dxyz
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_hh, int_ss

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

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
    max_iso_not0 = harmonics%max_iso_not0
    dmax_iso_not0 = harmonics%dmax_iso_not0
    lmax_expansion = indso(1,max_iso_not0)
    max_s_harm = harmonics%max_s_harm

    ALLOCATE(g1(nr),g2(nr),gg(nr,0:2*maxl),dgg(nr,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "g1,g2,gg,dgg",(2*nr+4*nr*maxl)*dp_size)

    ALLOCATE(gVXCg_h(na,0:2*maxl),gVXCg_s(na,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gVXCg_h,gVXCg_s",4*na*maxl*dp_size)
    ALLOCATE(gVXGg_h(3,na,0:2*maxl),gVXGg_s(3,na,0:2*maxl),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gVXCg_h,gVXCg_s",12*na*maxl*dp_size)
    ALLOCATE(cg_list(2,nsoset(maxl)**2,max_s_harm),cg_n_list(max_s_harm),&
             dcg_list(2,nsoset(maxl)**2,max_s_harm),dcg_n_list(max_s_harm),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(int_hh,int_ss)

    ALLOCATE(matso_h(nsoset(maxl),nsoset(maxl)),&
             matso_s(nsoset(maxl),nsoset(maxl)),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                       "matso_h,matso_s",2*nsoset(maxl)*nsoset(maxl)*dp_size)

    CALL get_rho_atom(rho_atom=rho_atom, ga_Vlocal_gb_h=int_hh,ga_Vlocal_gb_s=int_ss)

    DO ispin=1,nspins

    g1 = 0.0_dp
    g2 = 0.0_dp
    m1 = 0
    DO iset1 = 1,nset
       n1 = nsoset(lmax(iset1))
       m2 = 0
       DO iset2 = 1,nset
          CALL get_none0_cg_list(my_CG,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                 max_s_harm,lmax_expansion,cg_list,cg_n_list,max_iso_not0_local,error)
          CPPrecondition(max_iso_not0_local.LE.max_iso_not0,cp_failure_level,routineP,error,failure)
          CALL get_none0_cg_list(my_CG_dxyz,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                 max_s_harm,lmax_expansion,dcg_list,dcg_n_list,dmax_iso_not0_local,error)
          !CPPrecondition(dmax_iso_not0_local.LE.dmax_iso_not0,cp_failure_level,routineP,error,failure)

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

                !test reduce expansion local densities
                IF(lmin12 .LE. lmax_expansion) THEN

                   gg = 0.0_dp
                   dgg = 0.0_dp

                   IF (lmin12 == 0) THEN
                      gg(1:nr,lmin12) = g1(1:nr)*g2(1:nr)
                   ELSE
                      gg(1:nr,lmin12) = grid_atom%rad2l(1:nr,lmin12)*g1(1:nr)*g2(1:nr)
                   END IF

                   !test reduce expansion local densities
                   IF(lmax12 .GT. lmax_expansion) lmax12 = lmax_expansion

                   DO l=lmin12+1,lmax12
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

              ! Cross Term
              DO l = lmin12,lmax12
                 DO ia = 1,na
                    DO ir = 1,nr
                       gVXCg_h(ia,l) = gVXCg_h(ia,l) + &
                            gg(ir,l) * vxc_h(ia,ir,ispin)+&
                            dgg(ir,l)*&
                            (vxg_h(1,ia,ir,ispin)*harmonics%a(1,ia)+&
                            vxg_h(2,ia,ir,ispin)*harmonics%a(2,ia)+&
                            vxg_h(3,ia,ir,ispin)*harmonics%a(3,ia))

                       gVXCg_s(ia,l) = gVXCg_s(ia,l) + &
                            gg(ir,l) * vxc_s(ia,ir,ispin)+&
                            dgg(ir,l)*&
                            (vxg_s(1,ia,ir,ispin)*harmonics%a(1,ia)+&
                            vxg_s(2,ia,ir,ispin)*harmonics%a(2,ia)+&
                            vxg_s(3,ia,ir,ispin)*harmonics%a(3,ia))

                       urad = grid_atom%oorad2l(ir,1)

                       gVXGg_h(1,ia,l) = gVXGg_h(1,ia,l) + &
                            vxg_h(1,ia,ir,ispin)*&
                            gg(ir,l)*urad

                       gVXGg_h(2,ia,l) = gVXGg_h(2,ia,l) + &
                            vxg_h(2,ia,ir,ispin)*&
                            gg(ir,l)*urad

                       gVXGg_h(3,ia,l) = gVXGg_h(3,ia,l) + &
                            vxg_h(3,ia,ir,ispin)*&
                            gg(ir,l)*urad

                       gVXGg_s(1,ia,l) = gVXGg_s(1,ia,l) + &
                            vxg_s(1,ia,ir,ispin)*&
                            gg(ir,l)*urad

                       gVXGg_s(2,ia,l) = gVXGg_s(2,ia,l) + &
                            vxg_s(2,ia,ir,ispin)*&
                            gg(ir,l)*urad

                       gVXGg_s(3,ia,l) = gVXGg_s(3,ia,l) + &
                            vxg_s(3,ia,ir,ispin)*&
                            gg(ir,l)*urad

                    ENDDO ! ir
                 ENDDO ! ia
              ENDDO  ! l

              matso_h = 0.0_dp
              matso_s = 0.0_dp
              DO iso = 1,max_iso_not0_local
                 DO icg = 1,cg_n_list(iso)
                    iso1 = cg_list(1,icg,iso)
                    iso2 = cg_list(2,icg,iso)

                    l = indso(1,iso1) + indso(1,iso2)

                    !test reduce expansion local densities
                    CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)
                    IF(.NOT. failure) THEN

                       DO ia = 1,na
                          matso_h(iso1,iso2) = matso_h(iso1,iso2) + &
                                               gVXCg_h(ia,l)*&
                                               harmonics%slm(ia,iso)*&
                                               my_CG(iso1,iso2,iso)
                          matso_s(iso1,iso2) = matso_s(iso1,iso2) + &
                                               gVXCg_s(ia,l)*&
                                               harmonics%slm(ia,iso)*&
                                               my_CG(iso1,iso2,iso)
                       END DO ! ia

                    !test reduce expansion local densities
                    END IF ! failure

                 END DO

              END DO ! iso

              DO iso = 1,dmax_iso_not0_local
                 DO icg = 1,dcg_n_list(iso)
                    iso1 = dcg_list(1,icg,iso)
                    iso2 = dcg_list(2,icg,iso)

                    l = indso(1,iso1) + indso(1,iso2)
                    !test reduce expansion local densities
                    CPPrecondition(l<=lmax_expansion,cp_failure_level,routineP,error,failure)
                    IF(.NOT. failure) THEN

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

                       END DO ! ia

                    !test reduce expansion local densities
                    END IF ! failure

                 END DO ! icg
              END DO ! iso
              !test reduce expansion local densities
           END IF ! lmax_expansion

!             Write in the global matrix
              DO ic = nsoset(lmin(iset2)-1)+1,nsoset(lmax(iset2))
                iso1 = nsoset(lmin(iset1)-1)+1
                iso2 = ngau2+ic
                CALL daxpy(size1,1.0_dp,matso_h(iso1,ic),1,&
                           int_hh(ispin)%r_coef(nngau1+1,iso2),1)
                CALL daxpy(size1,1.0_dp,matso_s(iso1,ic),1,&
                           int_ss(ispin)%r_coef(nngau1+1,iso2),1)
               END DO

          END DO ! ipfg2
        END DO ! ipfg1
        m2 = m2 + maxso
      END DO ! iset2
      m1 = m1 + maxso
    END DO  ! iset1
    END DO ! ispin

    DEALLOCATE(g1,g2,gg,dgg,matso_h,matso_s,gVXCg_h,gVXCg_s,gVXGg_h,gVXGg_s,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
         "g1,g2,gg,matso_h,matso_s,gVXCg_s,gVXCg_h,gVXGg_h,gVXGg_s")
    DEALLOCATE(cg_list,cg_n_list,dcg_list,dcg_n_list,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE gaVxcgb_GC

! *****************************************************************************
  SUBROUTINE  dgaVtaudgb(vtau_h,vtau_s,atom_kind,rho_atom,nspins,error)

    REAL(dp), DIMENSION(:, :, :), POINTER    :: vtau_h, vtau_s
    TYPE(atomic_kind_type), POINTER          :: atom_kind
    TYPE(rho_atom_type), POINTER             :: rho_atom
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dgaVtaudgb', &
      routineP = moduleN//':'//routineN

    INTEGER :: dmax_iso_not0_local, handle, ic, icg, ipgf1, ipgf2, iset1, &
      iset2, iso, iso1, iso2, ispin, istat, l1, l2, lmax12, lmax_expansion, &
      lmin12, m1, m2, max_iso_not0, max_iso_not0_local, max_s_harm, maxiso, &
      maxl, maxso, mm1, mm2, n1, n2, na, ngau1, ngau2, nngau1, nr, nset, size1
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: cg_n_list, dcg_n_list
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: cg_list, dcg_list
    INTEGER, DIMENSION(:), POINTER           :: lmax, lmin, npgf
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: dd, g1, g2, gg, gr, ww
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: dgr, dvth, dvts, matso_h, &
                                                matso_s, vth, vts
    REAL(dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: fgr
    REAL(dp), DIMENSION(:, :), POINTER       :: zet
    REAL(dp), DIMENSION(:, :, :), POINTER    :: my_CG, my_dCG
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(gto_basis_set_type), POINTER        :: orb_basis
    TYPE(harmonics_atom_type), POINTER       :: harmonics
    TYPE(rho_atom_coeff), DIMENSION(:), &
      POINTER                                :: int_hh, int_ss

! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.

    NULLIFY(harmonics,grid_atom)
    CALL get_atomic_kind(atomic_kind=atom_kind,orb_basis_set=orb_basis,&
                         harmonics=harmonics,grid_atom=grid_atom)

    NULLIFY(lmin,lmax,npgf,zet)
    CALL get_gto_basis_set(gto_basis_set=orb_basis,lmax=lmax,lmin=lmin,&
                           maxso=maxso,maxl=maxl,npgf=npgf,&
                           nset=nset,zet=zet)

    my_CG => harmonics%my_CG
    my_dCG => harmonics%my_dCG
    max_iso_not0 = harmonics%max_iso_not0
    lmax_expansion = indso(1,max_iso_not0)
    max_s_harm = harmonics%max_s_harm

    nr = grid_atom%nr
    na = grid_atom%ng_sphere
    ALLOCATE(g1(nr),g2(nr),gg(nr),dd(nr),gr(nr),ww(na),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(fgr(nr,0:maxl,0:maxl),dgr(nr,0:lmax_expansion),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(cg_list(2,nsoset(maxl)**2,max_s_harm),cg_n_list(max_s_harm),&
             dcg_list(2,nsoset(maxl)**2,max_s_harm),dcg_n_list(max_s_harm),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    ALLOCATE(matso_h(nsoset(maxl),nsoset(maxl)),&
             matso_s(nsoset(maxl),nsoset(maxl)),STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    NULLIFY(int_hh,int_ss)
    CALL get_rho_atom(rho_atom=rho_atom,ga_Vlocal_gb_h=int_hh,ga_Vlocal_gb_s=int_ss)

    vtau_h = 0.5_dp*vtau_h
    vtau_s = 0.5_dp*vtau_s

    DO ispin=1,nspins
      maxiso = SIZE(harmonics%slm,2)
      ALLOCATE(vth(nr,maxiso),vts(nr,maxiso),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      ALLOCATE(dvth(nr,maxiso),dvts(nr,maxiso),STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
      DO iso=1,maxiso
         ww(1:na) = harmonics%slm(1:na,iso)
         vth(1:nr,iso) = MATMUL(ww(1:na),vtau_h(1:na,1:nr,ispin))
         vts(1:nr,iso) = MATMUL(ww(1:na),vtau_s(1:na,1:nr,ispin))
         ww(1:na) = harmonics%slm(1:na,iso)*grid_atom%usin_pol(1:na)**2
         dvth(1:nr,iso) = MATMUL(ww(1:na),vtau_h(1:na,1:nr,ispin))
         dvts(1:nr,iso) = MATMUL(ww(1:na),vtau_s(1:na,1:nr,ispin))
      END DO
      m1 = 0
      DO iset1 = 1,nset
        n1 = nsoset(lmax(iset1))
        m2 = 0
        DO iset2 = 1,nset
           CALL get_none0_cg_list(my_CG,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                  max_s_harm,lmax_expansion,cg_list,cg_n_list,max_iso_not0_local,error)
           CPPrecondition(max_iso_not0_local.LE.max_iso_not0,cp_failure_level,routineP,error,failure)
           CALL get_none0_cg_list(my_dCG,lmin(iset1),lmax(iset1),lmin(iset2),lmax(iset2),&
                                  max_s_harm,lmax_expansion,dcg_list,dcg_n_list,dmax_iso_not0_local,error)

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

              IF(lmin12 .LE. lmax_expansion) THEN
                 fgr = 0._dp
                 dgr = 0._dp
                 DO l1=0,maxl
                    DO l2=0,maxl
                       IF (l1+l2 > lmax_expansion) CYCLE
                       IF (l1+l2 > 0) THEN
                          gg(1:nr) = g1(1:nr)*g2(1:nr) * grid_atom%rad2l(1:nr,l1+l2)
                       ELSE
                          gg(1:nr) = g1(1:nr)*g2(1:nr)
                       END IF
                       dd(1:nr) = REAL(l1*l2,dp)/grid_atom%rad2(1:nr) - &
                            2._dp*REAL(l1,dp)*zet(ipgf2,iset2) - &
                            2._dp*REAL(l2,dp)*zet(ipgf1,iset1) + &
                            4._dp*zet(ipgf1,iset1)*zet(ipgf2,iset2)*grid_atom%rad2(1:nr)
                       fgr(1:nr,l1,l2) = dd(1:nr) * gg(1:nr)
                       dgr(1:nr,l1+l2) = gg(1:nr)/grid_atom%rad2(1:nr)
                    END DO
                 END DO
              END IF

              matso_h = 0.0_dp
              matso_s = 0.0_dp
              IF(lmin12 .LE. lmax_expansion) THEN
                 DO iso = 1,max_iso_not0_local
                    DO icg = 1,cg_n_list(iso)
                       iso1 = cg_list(1,icg,iso)
                       iso2 = cg_list(2,icg,iso)

                       l1 = indso(1,iso1)
                       l2 = indso(1,iso2)
                       mm1 = indso(2,iso1)
                       mm2 = indso(2,iso2)

                       IF (l1+l2 > lmax_expansion) CYCLE

                       matso_h(iso1,iso2) = matso_h(iso1,iso2) + my_CG(iso1,iso2,iso)*&
                            DOT_PRODUCT(vth(1:nr,ispin),fgr(1:nr,l1,l2))
                       matso_s(iso1,iso2) = matso_s(iso1,iso2) + my_CG(iso1,iso2,iso)*&
                            DOT_PRODUCT(vts(1:nr,ispin),fgr(1:nr,l1,l2))
                       ! d azimuthal
!deb
IF(mm1*mm2 /= 0) THEN
  IF(ABS(indso(2,iso)) .LE. 1) THEN
     WRITE(*,*) " PROBLEM ::::: ",l1,mm1,l2,mm2,indso(1,iso),indso(2,iso)
  ENDIF
ENDIF
!deb
                       matso_h(iso1,iso2) = matso_h(iso1,iso2) + REAL(mm1*mm2,dp)*my_CG(iso1,iso2,iso)*&
                            DOT_PRODUCT(dvth(1:nr,ispin),dgr(1:nr,l1+l2))
                       matso_s(iso1,iso2) = matso_s(iso1,iso2) + REAL(mm1*mm2,dp)*my_CG(iso1,iso2,iso)*&
                            DOT_PRODUCT(dvts(1:nr,ispin),dgr(1:nr,l1+l2))
                    END DO
                 END DO
              END IF
              ! d polar
              DO iso = 1,dmax_iso_not0_local
                 DO icg = 1,dcg_n_list(iso)
                    iso1 = dcg_list(1,icg,iso)
                    iso2 = dcg_list(2,icg,iso)

                    l1 = indso(1,iso1)
                    l2 = indso(1,iso2)
                    matso_h(iso1,iso2) = matso_h(iso1,iso2) + my_dCG(iso1,iso2,iso)*&
                         DOT_PRODUCT(vth(1:nr,ispin),dgr(1:nr,l1+l2))
                    matso_s(iso1,iso2) = matso_s(iso1,iso2) + my_dCG(iso1,iso2,iso)*&
                         DOT_PRODUCT(vts(1:nr,ispin),dgr(1:nr,l1+l2))
                 END DO
              END DO

              ! Write in the global matrix
              DO ic = nsoset(lmin(iset2)-1)+1,nsoset(lmax(iset2))
                iso1 = nsoset(lmin(iset1)-1)+1
                iso2 = ngau2+ic
                CALL daxpy(size1,1.0_dp,matso_h(iso1,ic),1,int_hh(ispin)%r_coef(nngau1+1,iso2),1)
                CALL daxpy(size1,1.0_dp,matso_s(iso1,ic),1,int_ss(ispin)%r_coef(nngau1+1,iso2),1)
              END DO

            END DO ! ipfg2
          END DO ! ipfg1
          m2 = m2 + maxso
        END DO ! iset2
        m1 = m1 + maxso
      END DO  ! iset1
      DEALLOCATE(vth,vts,dvth,dvts,STAT=istat)
      CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    END DO ! ispin

    DEALLOCATE(matso_h,matso_s,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(g1,g2,gg,dd,gr,ww,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(fgr,dgr,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(cg_list,cg_n_list,dcg_list,dcg_n_list,STAT=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    vtau_h = 2._dp*vtau_h
    vtau_s = 2._dp*vtau_s

    CALL timestop(handle)

  END SUBROUTINE dgaVtaudgb

END MODULE qs_vxc_atom
