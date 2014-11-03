!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that build the integrals of the Vxc potential calculated
!>      for the atomic code
! *****************************************************************************
MODULE atom_xc
  USE atom_types,                      ONLY: GTH_PSEUDO,&
                                             atom_type,&
                                             opmat_type
! ZMP added coulomb_potential_numeric
  USE atom_utils,                      ONLY: atom_core_density,&
                                             atom_density,&
                                             coulomb_potential_numeric,&
                                             integrate_grid,&
                                             numpot_matrix
  USE cp_files,                        ONLY: close_file,&
                                             get_unit_number,&
                                             open_file
  USE input_constants,                 ONLY: xc_none
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: fourpi
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_atom,                         ONLY: xc_rho_set_atom_update
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_create,&
                                             xc_dset_get_derivative,&
                                             xc_dset_release,&
                                             xc_dset_zero_all
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_derivatives,                  ONLY: xc_functionals_eval,&
                                             xc_functionals_get_needs
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_xc'

! ZMP added public subroutines calculate_atom_zmp, calculate_atom_ext_vxc 
  PUBLIC :: calculate_atom_vxc_lda, calculate_atom_vxc_lsd, &
            calculate_atom_zmp, calculate_atom_ext_vxc

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ZMP subroutine for the calculation of the effective constraint 
!>        potential in the atomic code.
!> 
!> \param ext_density ...
!> \param atom ...
!> \param lprint ...
!> \param xcmat ...
!> \param error ...
!> \author D. Varsano [daniele.varsano@nano.cnr.it] 
!>
! *****************************************************************************
  SUBROUTINE calculate_atom_zmp(ext_density,atom,lprint,xcmat,error)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: ext_density
    TYPE(atom_type), INTENT(INOUT)           :: atom
    LOGICAL                                  :: lprint
    TYPE(opmat_type), OPTIONAL, POINTER      :: xcmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atom_zmp', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extunit, handle, ierr, ir, &
                                                nr, z
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: int1, int2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: deltarho, rho_dum, vxc, vxc1, &
                                                vxc2
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rho

    CALL timeset(routineN,handle)
    failure = .FALSE.

    nr = atom%basis%grid%nr
    z  = atom%z
    ALLOCATE(rho(nr,1),vxc(nr),vxc1(nr),vxc2(nr),rho_dum(nr),deltarho(nr),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    CALL atom_density(rho(:,1),atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
    !Vxc1 
    int1=integrate_grid(ext_density,atom%basis%grid)
    int1=fourpi*int1
    CALL coulomb_potential_numeric(vxc1,rho(:,1),atom%basis%grid,error)

    vxc1 = -vxc1/z

    !Vxc2
    rho_dum=rho(:,1)*int1/z
    deltarho=rho_dum-ext_density
    int2=integrate_grid(deltarho,atom%basis%grid)
    CALL coulomb_potential_numeric(vxc2,deltarho,atom%basis%grid,error)

    vxc2=vxc2*atom%lambda

    !Vxc
    vxc=vxc1+vxc2

    atom%energy%exc = fourpi*integrate_grid(vxc,rho(:,1),atom%basis%grid)
    atom%rho_diff_integral = fourpi*int2

    IF (lprint) THEN
      extunit=get_unit_number()
      CALL open_file(file_name="linear_potentials.dat",file_status="UNKNOWN",&
                     file_form="FORMATTED", file_action="WRITE",&
                     unit_number=extunit)

      WRITE(extunit,fmt='("#",T10,"R[bohr]",T36,"Rho[au]",T61,"DRho[au]",T86,"V_XC[au]",T111,"V_XC1[au]",T136,"V_XC2[au]")')
      DO ir=1,nr
        WRITE(extunit,fmt='(T1,E24.15,T26,E24.15,T51,E24.15,T76,E24.15,T101,E24.15,T126,E24.15)') &
                            atom%basis%grid%rad(ir), rho(ir,1), deltarho(ir), vxc(ir), vxc1(ir), vxc2(ir)
      ENDDO
      CALL close_file(unit_number=extunit)
    ENDIF

    IF(PRESENT(xcmat)) CALL numpot_matrix(xcmat%op,vxc,atom%basis,0,error)

    DEALLOCATE(rho,vxc,vxc1,vxc2,rho_dum,deltarho,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE calculate_atom_zmp
! *****************************************************************************
!> \brief ZMP subroutine for the scf external density from the external v_xc 
!> 
!> \param vxc ...
!> \param atom ...
!> \param lprint ...
!> \param xcmat ...
!> \param error ...
!> \author D. Varsano [daniele.varsano@nano.cnr.it] 
!>
! *****************************************************************************
  SUBROUTINE calculate_atom_ext_vxc(vxc,atom,lprint,xcmat,error)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: vxc
    TYPE(atom_type), INTENT(INOUT)           :: atom
    LOGICAL, INTENT(in)                      :: lprint
    TYPE(opmat_type), INTENT(inout), &
      OPTIONAL, POINTER                      :: xcmat
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atom_ext_vxc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: extunit, handle, ierr, ir, nr
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: rho

    CALL timeset(routineN,handle)
    failure = .FALSE.
    nr = atom%basis%grid%nr

    ALLOCATE(rho(nr,1),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    CALL atom_density(rho(:,1),atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)

    IF (lprint) THEN
      extunit=get_unit_number()
      CALL open_file(file_name="linear_potentials.dat",file_status="UNKNOWN",&
                     file_form="FORMATTED", file_action="WRITE",&
                     unit_number=extunit)

      WRITE(extunit,fmt='("#",T10,"R[bohr]",T36,"Rho[au]",T61,"V_XC[au]")')
      DO ir=1,nr
        WRITE(extunit,fmt='(T1,E24.15,T26,E24.15,T51,E24.15)') &
                            atom%basis%grid%rad(ir), rho(ir,1), vxc(ir)
      ENDDO
      CALL close_file(unit_number=extunit)
    ENDIF

    atom%energy%exc = fourpi*integrate_grid(vxc,rho(:,1),atom%basis%grid)
    IF(PRESENT(xcmat)) CALL numpot_matrix(xcmat%op,vxc,atom%basis,0,error)

    DEALLOCATE(rho,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE calculate_atom_ext_vxc
! *****************************************************************************
!> \brief ...
!> \param xcmat ...
!> \param atom ...
!> \param xc_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_atom_vxc_lda(xcmat,atom,xc_section,error)
    TYPE(opmat_type), POINTER                :: xcmat
    TYPE(atom_type), INTENT(INOUT)           :: atom
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atom_vxc_lda', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: deriv_order, handle, i, ierr, &
                                                l, myfun, n1, n2, n3, nr, &
                                                nspins
    INTEGER, DIMENSION(2, 3)                 :: bounds
    LOGICAL                                  :: failure, lsd, nlcc
    REAL(KIND=dp)                            :: density_cut, gradient_cut, &
                                                tau_cut
    REAL(KIND=dp), DIMENSION(:), POINTER     :: exc, vxc
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: drho, rho, tau
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: taumat, xcpot
    TYPE(section_vals_type), POINTER         :: xc_fun_section
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho_set

    failure = .FALSE.
! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    nlcc=.FALSE.
    IF ( atom%potential%ppot_type==GTH_PSEUDO ) THEN
      nlcc = atom%potential%gth_pot%nlcc
    END IF

    IF ( ASSOCIATED(xc_section) ) THEN
      NULLIFY(rho_set)

      xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",error=error)
      CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",i_val=myfun,error=error)

      IF(myfun == xc_none) THEN
        atom%energy%exc = 0._dp
      ELSE
        CALL section_vals_val_get(xc_section,"DENSITY_CUTOFF",r_val=density_cut,error=error)
        CALL section_vals_val_get(xc_section,"GRADIENT_CUTOFF",r_val=gradient_cut,error=error)
        CALL section_vals_val_get(xc_section,"TAU_CUTOFF",r_val=tau_cut,error=error)

        lsd = .FALSE.
        nspins = 1
        needs = xc_functionals_get_needs(xc_fun_section,lsd=lsd,add_basic_components=.FALSE.,error=error)

        ! Prepare the structures needed to calculate and store the xc derivatives

        ! Array dimension: here anly one dimensional arrays are used,
        ! i.e. only the first column of deriv_data is read.
        ! The other to dimensions  are set to size equal 1
        nr = atom%basis%grid%nr
        bounds(1:2,1:3) = 1
        bounds(2,1) = nr

        ! create a place where to put the derivatives
        NULLIFY(deriv_set)
        CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)
        ! create the place where to store the argument for the functionals
        CALL xc_rho_set_create(rho_set,bounds,rho_cutoff=density_cut,&
               drho_cutoff=gradient_cut,tau_cutoff=tau_cut,error=error)
        ! allocate the required 3d arrays where to store rho and drho
        CALL xc_rho_set_atom_update(rho_set,needs,nspins,bounds)

        NULLIFY(rho,drho,tau)
        IF ( needs%rho ) THEN
          ALLOCATE(rho(nr,1),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          CALL atom_density(rho(:,1),atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
          IF ( nlcc ) THEN
            CALL atom_core_density(rho(:,1),atom%potential%gth_pot,typ="RHO",rr=atom%basis%grid%rad,error=error)
          END IF
        END IF
        IF ( needs%norm_drho ) THEN
          ALLOCATE(drho(nr,1),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          CALL atom_density(drho(:,1),atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="DER",error=error)
          IF ( nlcc ) THEN
            CALL atom_core_density(drho(:,1),atom%potential%gth_pot,typ="DER",rr=atom%basis%grid%rad,error=error)
          END IF
        END IF
        IF ( needs%tau ) THEN
          ALLOCATE(tau(nr,1),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          CALL atom_density(tau(:,1),atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,&
                            typ="KIN",rr=atom%basis%grid%rad2,error=error)
        END IF

        CALL fill_rho_set(rho_set,nspins,needs,rho,drho,tau,nr,error=error)

        CALL xc_dset_zero_all(deriv_set, error)

        deriv_order = 1
        CALL xc_functionals_eval(xc_fun_section,lsd=lsd,rho_set=rho_set,deriv_set=deriv_set,&
                                 deriv_order=deriv_order,error=error)

        ! Integration to get the matrix elements and energy
        deriv => xc_dset_get_derivative(deriv_set,"",allocate_deriv=.FALSE., error=error)
        CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
        atom%energy%exc = fourpi*integrate_grid(xcpot(:,1,1),atom%basis%grid)

        ! dump grid density and xcpot (xc energy?)
        IF (.FALSE.) THEN
           OPEN(UNIT=17,FILE="atom.dat")
           DO i=1,SIZE(xcpot(:,1,1))
              WRITE(17,*) atom%basis%grid%rad(i),rho(i,1),xcpot(i,1,1)
           ENDDO
           CLOSE(UNIT=17)
        ENDIF

        IF ( needs%rho ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(rho)",allocate_deriv=.FALSE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
          CALL numpot_matrix(xcmat%op,xcpot(:,1,1),atom%basis,0,error)
          DEALLOCATE(rho,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        END IF
        IF ( needs%norm_drho ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",allocate_deriv=.FALSE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
          CALL numpot_matrix(xcmat%op,xcpot(:,1,1),atom%basis,1,error)
          DEALLOCATE(drho,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        END IF
        IF ( needs%tau ) THEN
          deriv => xc_dset_get_derivative(deriv_set,"(tau)",allocate_deriv=.FALSE.,error=error)
          CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
          n1 = SIZE(xcmat%op,1)
          n2 = SIZE(xcmat%op,2)
          n3 = SIZE(xcmat%op,3)
          ALLOCATE(taumat(n1,n2,0:n3-1),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          taumat = 0._dp

          xcpot(:,1,1) = 0.5_dp * xcpot(:,1,1)
          CALL numpot_matrix(xcmat%op,xcpot(:,1,1),atom%basis,2,error)
          xcpot(:,1,1) = xcpot(:,1,1)/atom%basis%grid%rad2(:)
          CALL numpot_matrix(taumat,xcpot(:,1,1),atom%basis,0,error)
          DO l=0,3
            xcmat%op(:,:,l) = xcmat%op(:,:,l) + REAL(l*(l+1),dp)*taumat(:,:,l)
          END DO

          DEALLOCATE(tau,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          DEALLOCATE(taumat,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        END IF

        ! Release the xc structure used to store the xc derivatives
        CALL xc_dset_release(deriv_set, error=error)
        CALL xc_rho_set_release(rho_set,error=error)

      END IF !xc_none

    ELSE

      ! we don't have an xc_section, use a default setup
      nr = atom%basis%grid%nr
      ALLOCATE(rho(nr,1),exc(nr),vxc(nr),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

      CALL atom_density(rho(:,1),atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
      IF ( nlcc ) THEN
        CALL atom_core_density(rho(:,1),atom%potential%gth_pot,typ="RHO",rr=atom%basis%grid%rad,error=error)
      END IF
      CALL lda_pade(rho(:,1),exc,vxc,error)

      atom%energy%exc = fourpi*integrate_grid(exc,atom%basis%grid)
      CALL numpot_matrix(xcmat%op,vxc,atom%basis,0,error)

      DEALLOCATE(rho,exc,vxc,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_atom_vxc_lda

! *****************************************************************************
!> \brief ...
!> \param xcmata ...
!> \param xcmatb ...
!> \param atom ...
!> \param xc_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_atom_vxc_lsd(xcmata,xcmatb,atom,xc_section,error)
    TYPE(opmat_type), POINTER                :: xcmata, xcmatb
    TYPE(atom_type), INTENT(INOUT)           :: atom
    TYPE(section_vals_type), POINTER         :: xc_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atom_vxc_lsd', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: deriv_order, handle, ierr, l, &
                                                myfun, n1, n2, n3, nr, nspins
    INTEGER, DIMENSION(2, 3)                 :: bounds
    LOGICAL                                  :: failure, lsd
    REAL(KIND=dp)                            :: density_cut, gradient_cut, &
                                                tau_cut
    REAL(KIND=dp), DIMENSION(:), POINTER     :: exc, vxca, vxcb
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: drho, rho, tau
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: taumat, xcpot
    TYPE(section_vals_type), POINTER         :: xc_fun_section
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_derivative_type), POINTER        :: deriv
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho_set

    failure = .FALSE.
! -------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    IF ( ASSOCIATED(xc_section) ) THEN
       NULLIFY(rho_set)

       xc_fun_section => section_vals_get_subs_vals(xc_section,"XC_FUNCTIONAL",error=error)
       CALL section_vals_val_get(xc_fun_section,"_SECTION_PARAMETERS_",i_val=myfun,error=error)

       IF(myfun == xc_none) THEN
         atom%energy%exc = 0._dp
       ELSE
         CALL section_vals_val_get(xc_section,"DENSITY_CUTOFF",r_val=density_cut,error=error)
         CALL section_vals_val_get(xc_section,"GRADIENT_CUTOFF",r_val=gradient_cut,error=error)
         CALL section_vals_val_get(xc_section,"TAU_CUTOFF",r_val=tau_cut,error=error)

         lsd = .TRUE.
         nspins = 2
         needs = xc_functionals_get_needs(xc_fun_section,lsd=lsd,add_basic_components=.FALSE.,error=error)

         ! Prepare the structures needed to calculate and store the xc derivatives

         ! Array dimension: here anly one dimensional arrays are used,
         ! i.e. only the first column of deriv_data is read.
         ! The other to dimensions  are set to size equal 1
         nr = atom%basis%grid%nr
         bounds(1:2,1:3) = 1
         bounds(2,1) = nr

         ! create a place where to put the derivatives
         NULLIFY(deriv_set)
         CALL xc_dset_create(deriv_set, local_bounds=bounds, error=error)
         ! create the place where to store the argument for the functionals
         CALL xc_rho_set_create(rho_set,bounds,rho_cutoff=density_cut,&
                drho_cutoff=gradient_cut,tau_cutoff=tau_cut,error=error)
         ! allocate the required 3d arrays where to store rho and drho
         CALL xc_rho_set_atom_update(rho_set,needs,nspins,bounds)

         NULLIFY(rho,drho,tau)
         IF ( needs%rho_spin ) THEN
           ALLOCATE(rho(nr,2),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           CALL atom_density(rho(:,1),atom%orbitals%pmata,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
           CALL atom_density(rho(:,2),atom%orbitals%pmatb,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
         END IF
         IF ( needs%norm_drho_spin ) THEN
           ALLOCATE(drho(nr,2),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           CALL atom_density(drho(:,1),atom%orbitals%pmata,atom%basis,atom%state%maxl_occ,typ="DER",error=error)
           CALL atom_density(drho(:,2),atom%orbitals%pmatb,atom%basis,atom%state%maxl_occ,typ="DER",error=error)
         END IF
         IF ( needs%tau_spin ) THEN
           ALLOCATE(tau(nr,2),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           CALL atom_density(tau(:,1),atom%orbitals%pmata,atom%basis,atom%state%maxl_occ,&
                             typ="KIN",rr=atom%basis%grid%rad2,error=error)
           CALL atom_density(tau(:,2),atom%orbitals%pmatb,atom%basis,atom%state%maxl_occ,&
                             typ="KIN",rr=atom%basis%grid%rad2,error=error)
         END IF

         CALL fill_rho_set(rho_set,nspins,needs,rho,drho,tau,nr,error=error)

         CALL xc_dset_zero_all(deriv_set, error)

         deriv_order = 1
         CALL xc_functionals_eval(xc_fun_section,lsd=lsd,rho_set=rho_set,deriv_set=deriv_set,&
                                  deriv_order=deriv_order,error=error)

         ! Integration to get the matrix elements and energy
         deriv => xc_dset_get_derivative(deriv_set,"",allocate_deriv=.FALSE., error=error)
         CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
         atom%energy%exc = fourpi*integrate_grid(xcpot(:,1,1),atom%basis%grid)

         IF ( needs%rho_spin ) THEN
           deriv => xc_dset_get_derivative(deriv_set,"(rhoa)",allocate_deriv=.FALSE.,error=error)
           CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
           CALL numpot_matrix(xcmata%op,xcpot(:,1,1),atom%basis,0,error)
           deriv => xc_dset_get_derivative(deriv_set,"(rhob)",allocate_deriv=.FALSE.,error=error)
           CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
           CALL numpot_matrix(xcmatb%op,xcpot(:,1,1),atom%basis,0,error)
           DEALLOCATE(rho,STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         END IF
         IF ( needs%norm_drho_spin ) THEN
           ! drhoa
           NULLIFY(deriv)
           deriv => xc_dset_get_derivative(deriv_set,"(norm_drhoa)",allocate_deriv=.FALSE.,error=error)
           CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
           CALL numpot_matrix(xcmata%op,xcpot(:,1,1),atom%basis,1,error)
           ! drhob
           NULLIFY(deriv)
           deriv => xc_dset_get_derivative(deriv_set,"(norm_drhob)",allocate_deriv=.FALSE.,error=error)
           CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
           CALL numpot_matrix(xcmatb%op,xcpot(:,1,1),atom%basis,1,error)
           ! Cross Terms
           NULLIFY(deriv)
           deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",error=error)
           IF(ASSOCIATED(deriv)) THEN
             CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
             CALL numpot_matrix(xcmata%op,xcpot(:,1,1),atom%basis,1,error)
             CALL numpot_matrix(xcmatb%op,xcpot(:,1,1),atom%basis,1,error)
           END IF
           DEALLOCATE(drho,STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         END IF
         IF ( needs%tau_spin ) THEN
           n1 = SIZE(xcmata%op,1)
           n2 = SIZE(xcmata%op,2)
           n3 = SIZE(xcmata%op,3)
           ALLOCATE(taumat(n1,n2,0:n3-1),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

           deriv => xc_dset_get_derivative(deriv_set,"(tau_a)",allocate_deriv=.FALSE.,error=error)
           CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
           taumat = 0._dp
           xcpot(:,1,1) = 0.5_dp * xcpot(:,1,1)
           CALL numpot_matrix(xcmata%op,xcpot(:,1,1),atom%basis,2,error)
           xcpot(:,1,1) = xcpot(:,1,1)/atom%basis%grid%rad2(:)
           CALL numpot_matrix(taumat,xcpot(:,1,1),atom%basis,0,error)
           DO l=0,3
             xcmata%op(:,:,l) = xcmata%op(:,:,l) + REAL(l*(l+1),dp)*taumat(:,:,l)
           END DO

           deriv => xc_dset_get_derivative(deriv_set,"(tau_b)",allocate_deriv=.FALSE.,error=error)
           CALL xc_derivative_get(deriv,deriv_data=xcpot,error=error)
           taumat = 0._dp
           xcpot(:,1,1) = 0.5_dp * xcpot(:,1,1)
           CALL numpot_matrix(xcmatb%op,xcpot(:,1,1),atom%basis,2,error)
           xcpot(:,1,1) = xcpot(:,1,1)/atom%basis%grid%rad2(:)
           CALL numpot_matrix(taumat,xcpot(:,1,1),atom%basis,0,error)
           DO l=0,3
             xcmatb%op(:,:,l) = xcmatb%op(:,:,l) + REAL(l*(l+1),dp)*taumat(:,:,l)
           END DO

           DEALLOCATE(tau,STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           DEALLOCATE(taumat,STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         END IF

         ! Release the xc structure used to store the xc derivatives
         CALL xc_dset_release(deriv_set, error=error)
         CALL xc_rho_set_release(rho_set,error=error)

       END IF !xc_none

     ELSE

      ! we don't have an xc_section, use a default setup
      nr = atom%basis%grid%nr
      ALLOCATE(rho(nr,2),exc(nr),vxca(nr),vxcb(nr),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

      CALL atom_density(rho(:,1),atom%orbitals%pmata,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
      CALL atom_density(rho(:,2),atom%orbitals%pmatb,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
      CALL lsd_pade(rho(:,1),rho(:,2),exc,vxca,vxcb,error)

      atom%energy%exc = fourpi*integrate_grid(exc,atom%basis%grid)
      CALL numpot_matrix(xcmata%op,vxca,atom%basis,0,error)
      CALL numpot_matrix(xcmatb%op,vxcb,atom%basis,0,error)

      DEALLOCATE(rho,exc,vxca,vxcb,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE calculate_atom_vxc_lsd

! *****************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param nspins ...
!> \param needs ...
!> \param rho ...
!> \param drho ...
!> \param tau ...
!> \param na ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE fill_rho_set(rho_set,nspins,needs,rho,drho,tau,na,error)

    TYPE(xc_rho_set_type), POINTER           :: rho_set
    INTEGER, INTENT(IN)                      :: nspins
    TYPE(xc_rho_cflags_type), INTENT(in)     :: needs
    REAL(dp), DIMENSION(:, :), POINTER       :: rho, drho, tau
    INTEGER, INTENT(IN)                      :: na
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fill_rho_set', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: f13 = (1.0_dp/3.0_dp)

    INTEGER                                  :: ia
    LOGICAL                                  :: failure

    failure = .FALSE.

    SELECT CASE(nspins)
    CASE(1)
       CPPrecondition(.NOT.needs%rho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%norm_drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drhoa_drhob,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%rho_spin_1_3,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%tau_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drho,cp_failure_level,routineP,error,failure)
       ! Give rho to 1/3
       IF (needs%rho_1_3) THEN
          DO ia = 1,na
             rho_set%rho_1_3(ia,1,1) = MAX(rho(ia,1),0.0_dp)**f13
          END DO
          rho_set%owns%rho_1_3=.TRUE.
          rho_set%has%rho_1_3=.TRUE.
       END IF
       ! Give the density
       IF (needs%rho) THEN
          DO ia = 1,na
             rho_set%rho(ia,1,1) = rho(ia,1)
          END DO
          rho_set%owns%rho=.TRUE.
          rho_set%has%rho=.TRUE.
       END IF
       ! Give the norm of the gradient of the density
       IF (needs%norm_drho) THEN
          DO ia = 1,na
             rho_set%norm_drho(ia,1,1) = drho(ia,1)
          END DO
          rho_set%owns%norm_drho=.TRUE.
          rho_set%has%norm_drho=.TRUE.
       END IF
    CASE(2)
       CPPrecondition(.NOT.needs%drho,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drhoa_drhob,cp_failure_level,routineP,error,failure)
       ! Give the total density
       IF (needs%rho) THEN
          DO ia = 1,na
             rho_set%rho(ia,1,1) = rho(ia,1) + rho(ia,2)
          END DO
          rho_set%owns%rho=.TRUE.
          rho_set%has%rho=.TRUE.
       END IF
       ! Give the norm of the total gradient of the density
       IF (needs%norm_drho) THEN
          DO ia = 1,na
            rho_set%norm_drho(ia,1,1) = drho(ia,1) + drho(ia,2)
          END DO
          rho_set%owns%norm_drho=.TRUE.
          rho_set%has%norm_drho=.TRUE.
       END IF
       ! Give rho_spin
       IF (needs%rho_spin) THEN
         DO ia = 1,na
            rho_set%rhoa(ia,1,1) = rho(ia,1)
            rho_set%rhob(ia,1,1) = rho(ia,2)
         END DO
         rho_set%owns%rho_spin=.TRUE.
         rho_set%has%rho_spin=.TRUE.
       END IF
       ! Give rho_spin to 1/3
       IF (needs%rho_spin_1_3) THEN
          DO ia = 1,na
             rho_set%rhoa_1_3(ia,1,1) = MAX(rho(ia,1),0.0_dp)**f13
             rho_set%rhob_1_3(ia,1,1) = MAX(rho(ia,2),0.0_dp)**f13
          END DO
          rho_set%owns%rho_1_3=.TRUE.
          rho_set%has%rho_1_3=.TRUE.
       END IF
       ! Give the norm of the gradient of rhoa and of rhob separatedly
       IF (needs%norm_drho_spin) THEN
          DO ia = 1,na
             rho_set%norm_drhoa(ia,1,1) = drho(ia,1)
             rho_set%norm_drhob(ia,1,1) = drho(ia,2)
          END DO
          rho_set%owns%norm_drho_spin=.TRUE.
          rho_set%has%norm_drho_spin=.TRUE.
       END IF
       !
    END SELECT

    ! tau part
    IF (needs%tau) THEN
       IF (nspins==2) THEN
          DO ia = 1,na
             rho_set%tau(ia,1,1) = tau(ia,1)+tau(ia,2)
          END DO
          rho_set%owns%tau=.TRUE.
          rho_set%has%tau=.TRUE.
       ELSE
          DO ia = 1,na
             rho_set%tau(ia,1,1) = tau(ia,1)
          END DO
          rho_set%owns%tau=.TRUE.
          rho_set%has%tau=.TRUE.
       END IF
    END IF
    IF (needs%tau_spin) THEN
       CPPrecondition(nspins==2,cp_failure_level,routineP,error,failure)
       DO ia = 1,na
          rho_set%tau_a(ia,1,1) = tau(ia,1)
          rho_set%tau_b(ia,1,1) = tau(ia,2)
       END DO
       rho_set%owns%tau_spin=.TRUE.
       rho_set%has%tau_spin=.TRUE.
    END IF

  END SUBROUTINE fill_rho_set
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param exc ...
!> \param vxc ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE lda_pade(rho,exc,vxc,error)

    REAL(dp), DIMENSION(:)                   :: rho, exc, vxc
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lda_pade', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: a0 = 0.4581652932831429E+0_dp, &
      a1 = 0.2217058676663745E+1_dp, a2 = 0.7405551735357053E+0_dp, &
      a3 = 0.1968227878617998E-1_dp, b1 = 1.0000000000000000E+0_dp, &
      b2 = 0.4504130959426697E+1_dp, b3 = 0.1110667363742916E+1_dp, &
      b4 = 0.2359291751427506E-1_dp, f13 = (1.0_dp/3.0_dp), &
      rsfac = 0.6203504908994000166680065_dp

    INTEGER                                  :: i, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: depade, dpv, dq, epade, p, q, &
                                                rs

    failure = .FALSE.

    n = SIZE(rho)
    exc(1:n) = 0._dp
    vxc(1:n) = 0._dp

    DO i=1,n
      IF ( rho(i) > 1.e-20_dp ) THEN
         rs = rsfac * rho(i)**(-f13)
         p = a0 + (a1 + (a2 + a3*rs)*rs)*rs
         q = (b1 + (b2 + (b3 + b4*rs)*rs)*rs)*rs
         epade = -p/q

         dpv = a1 + (2.0_dp*a2 + 3.0_dp*a3*rs)*rs
         dq = b1 + (2.0_dp*b2 + (3.0_dp*b3 + 4.0_dp*b4*rs)*rs)*rs
         depade = f13 * rs * (dpv*q - p*dq) / (q*q)

         exc(i) = epade * rho(i)
         vxc(i) = epade + depade
      END IF
    END DO

  END SUBROUTINE lda_pade
! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param rhob ...
!> \param exc ...
!> \param vxca ...
!> \param vxcb ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE lsd_pade ( rhoa, rhob, exc, vxca, vxcb, error )

    REAL(dp), DIMENSION(:)                   :: rhoa, rhob, exc, vxca, vxcb
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'lsd_pade', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: a0 = 0.4581652932831429E+0_dp, &
      a1 = 0.2217058676663745E+1_dp, a2 = 0.7405551735357053E+0_dp, &
      a3 = 0.1968227878617998E-1_dp, b1 = 1.0000000000000000E+0_dp, &
      b2 = 0.4504130959426697E+1_dp, b3 = 0.1110667363742916E+1_dp, &
      b4 = 0.2359291751427506E-1_dp, da0 = 0.119086804055547E+0_dp, &
      da1 = 0.6157402568883345E+0_dp, da2 = 0.1574201515892867E+0_dp, &
      da3 = 0.3532336663397157E-2_dp, db1 = 0.0000000000000000E+0_dp, &
      db2 = 0.2673612973836267E+0_dp, db3 = 0.2052004607777787E+0_dp, &
      db4 = 0.4200005045691381E-2_dp, f13 = (1.0_dp/3.0_dp), &
      f43 = (4.0_dp/3.0_dp)
    REAL(KIND=dp), PARAMETER :: fxfac = 1.923661050931536319759455_dp , &
      rsfac = 0.6203504908994000166680065_dp

    INTEGER                                  :: i, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dc, dpv, dq, dr, dx, fa0, &
                                                fa1, fa2, fa3, fb1, fb2, fb3, &
                                                fb4, fx1, fx2, p, q, rhoab, &
                                                rs, x, xp, xq

! 1/(2^(4/3) - 2)

    failure = .FALSE.

    n = SIZE(rhoa)
    exc(1:n) = 0._dp
    vxca(1:n) = 0._dp
    vxcb(1:n) = 0._dp

    DO i=1,n
      rhoab = rhoa(i) + rhob(i)
      IF ( rhoab > 1.e-20_dp ) THEN
         rs = rsfac * rhoab**(-f13)

         x = (rhoa(i) - rhob(i)) / rhoab
         IF ( x < -1.0_dp ) THEN
            fx1 = 1.0_dp
            fx2 = -f43*fxfac*2.0_dp**f13
         ELSE IF ( x > 1.0_dp ) THEN
            fx1 = 1.0_dp
            fx2 = f43*fxfac*2.0_dp**f13
         ELSE
            fx1 = ( (1.0_dp+x)**f43 + (1.0_dp-x)**f43 - 2.0_dp ) * fxfac
            fx2 = ( (1.0_dp+x)**f13 - (1.0_dp-x)**f13 ) * fxfac * f43
         END IF

         fa0 = a0 + fx1*da0
         fa1 = a1 + fx1*da1
         fa2 = a2 + fx1*da2
         fa3 = a3 + fx1*da3
         fb1 = b1 + fx1*db1
         fb2 = b2 + fx1*db2
         fb3 = b3 + fx1*db3
         fb4 = b4 + fx1*db4

         p = fa0 + (fa1 + (fa2 + fa3*rs)*rs)*rs
         q = (fb1 + (fb2 + (fb3 + fb4*rs)*rs)*rs)*rs
         dpv = fa1 + (2.0_dp*fa2 + 3.0_dp*fa3*rs)*rs
         dq = fb1 + (2.0_dp*fb2 + (3.0_dp*fb3 + &
                     4.0_dp*fb4*rs)*rs)*rs
         xp = da0 + (da1 + (da2 + da3*rs)*rs)*rs
         xq = (db1 + (db2 + (db3 + db4*rs)*rs)*rs)*rs

         dr = (dpv*q - p*dq)/(q*q)
         dx = 2.0_dp*(xp*q - p*xq)/(q*q)*fx2/rhoab
         dc = f13*rs*dr - p/q

         exc(i) =  -p/q*rhoab
         vxca(i) = dc - dx*rhob(i)
         vxcb(i) = dc + dx*rhoa(i)
      END IF
   END DO

  END SUBROUTINE lsd_pade

! *****************************************************************************

END MODULE atom_xc
