
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/xc_atom
!!
!!   NAME
!!     xc_atom
!!     
!!
!!
!!   SOURCE
!****************************************************************************
MODULE xc_atom 

  USE cp_control_types,                ONLY: xc_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_linked_list_xc_deriv,         ONLY: cp_sll_xc_deriv_next,&
                                             cp_sll_xc_deriv_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_functionals,                  ONLY: xc_calculate_lda,&
                                             xc_calculate_lsd,&
                                             xc_get_pot_size,&
                                             xc_get_order,&
                                             xc_set
  USE xc_derivative_desc,              ONLY: MAX_DERIVATIVE_DESC_LENGTH
  USE xc_derivative_set_types,         ONLY: xc_dset_get_derivative,&
                                             xc_derivative_set_type
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_derivatives,                  ONLY: xc_functional_eval
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_type,&
                                             xc_rho_set_get,&
                                             xc_rho_set_type

 IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PRIVATE, PARAMETER :: module_name='xc_atom'

  PUBLIC :: vxc_of_r, vxc_of_r_new, xc_rho_set_atom_update

!******************************************************************************

CONTAINS

!******************************************************************************

  SUBROUTINE vxc_of_r(xc_control,derivative_set,nspins,na,rho,drho,ir,locpot,error)

    TYPE(xc_control_type)               :: xc_control
    TYPE(xc_derivative_set_type), &
                             POINTER    :: derivative_set
    INTEGER                             :: nspins,na
    REAL(dp), DIMENSION(:,:), &
                             POINTER    :: rho
    INTEGER, INTENT(IN), OPTIONAL       :: ir
    REAL(dp), DIMENSION(:,:,:,:), &
                     OPTIONAL, POINTER  :: drho
    REAL(dp), DIMENSION(:,:), OPTIONAL :: locpot
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                          :: error

    CHARACTER(LEN=*), PARAMETER         :: routine_name = "vxc_of_r"

    REAL(dp), DIMENSION(:), POINTER     :: rhoa,rhob,drhoa,drhob
    CHARACTER(LEN=40), DIMENSION(3)     :: funct
    LOGICAL, DIMENSION(3)               :: grad_funct
    LOGICAL                             :: docalc( 3 ),cross_term,&
                                           gradient_f,lsd
    INTEGER  :: i,ia,ifun,ind(100),istat,lp,npot,order

!******************************************************************************

   gradient_f=PRESENT(drho)
   cross_term = ANY (xc_control%crossterms) .AND. nspins>1

   lsd = (nspins==2)
   order = 1
   npot = xc_get_pot_size(lsd,gradient_f , cross_term, order)

   funct(:) = xc_control%functionals(:)
   grad_funct(:) = xc_control%gradient_functionals(:)

   docalc ( : ) = ( INDEX ( funct(:), "NONE" ) == 0 )
   CALL xc_set(xc_control%density_cut, xc_control%gradient_cut)
   IF(lsd) THEN
      ALLOCATE(rhoa(na),rhob(na),STAT=istat)
      DO ia = 1,na
        rhoa(ia) = rho(ia,1)
        rhob(ia) = rho(ia,2)
      END DO
      IF(PRESENT(drho)) THEN
        ALLOCATE(drhoa(na),drhob(na),STAT=istat)
        DO ia = 1,na
          drhoa(ia) = drho(4,ia,ir,1) 
          drhob(ia) = drho(4,ia,ir,2)
        END DO
      END IF
      DO ifun=1, 3
        IF (docalc(ifun)) THEN
          IF (gradient_f) THEN
            CALL dcopy(na*npot, 0.0_dp, 0, locpot, 1)
            lp = xc_get_pot_size(.FALSE., grad_funct(ifun), .FALSE., order)
            CALL xc_calculate_lsd(derivative_set, funct(ifun),&
                                  rhoa,rhob,locpot, order,&
                                  drhoa=drhoa,drhob=drhob,error=error)
            CALL xc_get_order(.FALSE., grad_funct(ifun), .FALSE., &
                              gradient_f, .FALSE., order, ind(1:lp))
            DO i = 1, lp
               CALL daxpy(na,1._dp, locpot(1,i), 1,& 
                          derivative_set%data(1,ind(i)-1), 1)
            END DO
          ELSE
            CALL xc_calculate_lsd(derivative_set, funct(ifun),&
                                  rhoa,rhob,&
                                  derivative_set%data, order,error=error)
          END IF
        END IF
      END DO
      DEALLOCATE(rhoa,rhob,STAT=istat)
      IF(PRESENT(drho)) THEN
        DEALLOCATE(drhoa,drhob,STAT=istat)
      END IF
   ELSE
      ALLOCATE(rhoa(na),STAT=istat)
      DO ia = 1,na
        rhoa(ia) = rho(ia,1)
      END DO
      IF(PRESENT(drho)) THEN
        ALLOCATE(drhoa(na),STAT=istat) 
        DO ia = 1,na
          drhoa(ia) = drho(4,ia,ir,1)
        END DO
      ENDIF
      DO ifun=1, 3
        IF (docalc(ifun)) THEN
          IF (gradient_f) THEN
            CALL dcopy(na*npot, 0.0_dp, 0, locpot, 1)
            lp = xc_get_pot_size(.FALSE., grad_funct(ifun), .FALSE., order)
            CALL xc_calculate_lda(derivative_set, functional=funct(ifun),&
                                 rho=rhoa, drho=drhoa, &
                                  pot=locpot, order=order, error=error)
            CALL xc_get_order(.FALSE., grad_funct(ifun), .FALSE., &
                              gradient_f, .FALSE., order, ind(1:lp))
            DO i = 1, lp
              CALL daxpy(na,1._dp, locpot(1,i), 1,& 
                         derivative_set%data(1,ind(i)-1), 1)
            END DO
          ELSE

            CALL xc_calculate_lda(derivative_set, functional=funct(ifun),rho=rhoa,&
                                  pot=derivative_set%data, order=order, error=error)

          END IF
        END IF
      END DO
      DEALLOCATE (rhoa,STAT=istat)
      IF(PRESENT(drho)) THEN
        DEALLOCATE(drhoa,STAT=istat)
      END IF
    END IF

  END SUBROUTINE vxc_of_r

!****************************************************************************


  SUBROUTINE vxc_of_r_new(xc_control,rho_set,deriv_set,deriv_order,needs,rho,drho,oodrho,w,&
                          lsd,na,ir_pnt,exc,vxc,vxg,vxg_cross,error)

! This routine updates rho_set by giving to it the rho and drho that are needed.
! Since for the local densities rho1_h and rho1_s local grids are used it is not possible
! to call xc_rho_set_update. 
! As input of this routine one gets rho and drho on a one dimensional grid. 
! The grid is the angular grid corresponding to a given point ir_pnt on the radial grid.
! The derivatives are calculated on this one dimensional grid, the results are stored in 
! exc, vxc(1:na,ir_pnt,ispin), vxg(1:na,ir_pnt,ispin), vxg_cross(1:na,ir_pnt,ispin)
! Afterwords the arrays containing the derivatives are put to zero so that the routine
! can safely be called for the next radial point ir_pnt

    TYPE(xc_control_type), POINTER                :: xc_control
    TYPE(xc_rho_set_type), POINTER                :: rho_set
    TYPE(xc_derivative_set_type), POINTER         :: deriv_set
    INTEGER, INTENT(in)                           :: deriv_order
    TYPE(xc_rho_cflags_type), INTENT(in)          :: needs
    REAL(dp), DIMENSION(:,:),  POINTER            :: rho 
    REAL(dp), DIMENSION(:,:,:,:),  POINTER        :: drho 
    REAL(dp), DIMENSION(:,:),  POINTER            :: oodrho 
    REAL(dp), DIMENSION(:),  POINTER              :: w 
    LOGICAL, INTENT(IN)                           :: lsd
    INTEGER, INTENT(IN)                           :: na,ir_pnt
    REAL(dp)                                      :: exc
    REAL(dp), DIMENSION(:,:,:), POINTER           :: vxc,vxg
    REAL(dp), DIMENSION(:,:), POINTER             :: vxg_cross
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                    :: error

    CHARACTER(LEN=*), PARAMETER                   :: routine_name = "vxc_of_r_new",&
                                         routineP = module_name//':'//routine_name


    TYPE(cp_sll_xc_deriv_type), POINTER           :: pos
    TYPE(xc_derivative_type), POINTER             :: deriv_att
    REAL(KIND=dp), PARAMETER                      :: f13 = (1.0_dp/3.0_dp)
    REAL(dp)                                      :: drho_cutoff
    REAL(dp), DIMENSION(:,:,:), POINTER           :: deriv_data, norm_drho
    INTEGER :: handle, i, ia, idir, ispin, n_deriv, nspins, order
    LOGICAL                                       :: failure, gradient_f
    CHARACTER(len=MAX_DERIVATIVE_DESC_LENGTH)     :: desc
    CHARACTER(len=30), DIMENSION(2), PARAMETER    :: &
      norm_drho_spin_name = (/ "(norm_drhoa)", "(norm_drhob)" /)


!   ***************************************************************************

    CALL timeset(routine_name,"I"," ",handle)
    failure = .FALSE.

    nspins = 1
    IF(lsd) nspins = 2

    gradient_f=(needs%drho_spin.OR.needs%norm_drho_spin.OR.&
                needs%drhoa_drhob.OR.needs%drho.OR.needs%norm_drho)

    CALL fill_rho_set(rho_set,lsd,.FALSE.,needs,rho,drho,oodrho,w,na,ir_pnt,error)

!  Calculate the derivatives
    DO i=1,SIZE(xc_control%functionals)
       CALL xc_functional_eval(xc_control%functionals(i), &
            lsd=lsd,&
            rho_set=rho_set, &
            deriv_set=deriv_set,&
            deriv_order=deriv_order, &
            error=error)
    END DO

   NULLIFY (deriv_data)
!  EXC energy
   deriv_att => xc_dset_get_derivative(deriv_set,"", error=error)
   CPPrecondition(ASSOCIATED(deriv_att),cp_failure_level,routineP,error,failure)
   CALL xc_derivative_get(deriv_att,deriv_data=deriv_data, error=error)
   exc = SUM( deriv_data(:,1,1)*w(:))
   deriv_data = 0.0_dp   ! put to zero
   NULLIFY (deriv_data)
!  Derivative with respect to the density
   IF(lsd) THEN
     deriv_att => xc_dset_get_derivative(deriv_set, "(rhoa)")
     IF (ASSOCIATED(deriv_att)) THEN
       CALL xc_derivative_get(deriv_att,deriv_data=deriv_data, error=error)
       DO ia = 1,na
         vxc(ia,ir_pnt,1) = deriv_data(ia,1,1)*w(ia)
       END DO
       NULLIFY (deriv_data)
     END IF
     deriv_att => xc_dset_get_derivative(deriv_set, "(rhob)")
     IF (ASSOCIATED(deriv_att)) THEN
       CALL xc_derivative_get(deriv_att,deriv_data=deriv_data, error=error)
       DO ia = 1,na
         vxc(ia,ir_pnt,2) = deriv_data(ia,1,1)*w(ia)
       END DO
       NULLIFY(deriv_data)
     END IF
   ELSE
     deriv_att => xc_dset_get_derivative(deriv_set, "(rho)")
     IF (ASSOCIATED(deriv_att)) THEN
       CALL xc_derivative_get(deriv_att,deriv_data=deriv_data, error=error)
       DO ia = 1,na
         vxc(ia,ir_pnt,1) = deriv_data(ia,1,1)*w(ia)
       END DO
       NULLIFY (deriv_data)
     END IF
   END IF

   deriv_att => xc_dset_get_derivative(deriv_set, "(rho)")
   IF(lsd .AND. ASSOCIATED(deriv_att)) THEN
     CALL xc_derivative_get(deriv_att,deriv_data=deriv_data, error=error)
     DO ia = 1,na
       vxc(ia,ir_pnt,1) = vxc(ia,ir_pnt,1) + deriv_data(ia,1,1)*w(ia)
       vxc(ia,ir_pnt,2) = vxc(ia,ir_pnt,2) + deriv_data(ia,1,1)*w(ia)
     END DO
     NULLIFY (deriv_data)
   END IF

!  Derivatives with respect to the gradient
   IF (lsd) THEN
     DO ispin = 1,nspins
       deriv_att => xc_dset_get_derivative(deriv_set, norm_drho_spin_name(ispin))
       IF (ASSOCIATED(deriv_att)) THEN
          CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,error=error)
          DO ia = 1,na
            vxg(ia,ir_pnt,ispin) = deriv_data(ia,1,1)*w(ia)*oodrho(ia,ispin)
          END DO
          NULLIFY (deriv_data)
       END IF
     END DO
   ELSE
     deriv_att => xc_dset_get_derivative(deriv_set,"(norm_drho)")
     IF (ASSOCIATED(deriv_att)) THEN
       CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,error=error)
       DO ia = 1,na
         vxg(ia,ir_pnt,1) = deriv_data(ia,1,1)*w(ia)*oodrho(ia,1)
       END DO
       NULLIFY (deriv_data)
     END IF
   END IF

!  Cross Terms
   deriv_att => xc_dset_get_derivative(deriv_set, "(norm_drho)")
   IF(lsd .AND. ASSOCIATED(deriv_att)) THEN
     CALL xc_derivative_get(deriv_att,deriv_data=deriv_data,error=error)
     CALL xc_rho_set_get(rho_set=rho_set,norm_drho=norm_drho,drho_cutoff=drho_cutoff)
     CPPrecondition(ASSOCIATED(norm_drho),cp_failure_level,routineP,error,failure)
     CPPrecondition(ASSOCIATED(vxg_cross),cp_failure_level,routineP,error,failure)
     DO ia = 1,na
      vxg_cross(ia,ir_pnt) = deriv_data(ia,1,1)*w(ia)/MAX(norm_drho(ia,1,1),drho_cutoff)
     END DO
     NULLIFY (deriv_data)
   END IF

!  Put to zero the arrays for the derivatives to be ready for the next call to this routin
   n_deriv = 0
   pos => deriv_set%derivs
   DO WHILE (cp_sll_xc_deriv_next(pos,el_att=deriv_att, error=error)) 
     CALL xc_derivative_get(deriv_att,order=order,desc=desc,deriv_data=deriv_data)
     IF(order==1) THEN
       IF(lsd) THEN
         SELECT CASE(desc)
           CASE("(rho)","(rhoa)","(rhob)","(norm_drho)","(norm_drhoa)",&
                "(norm_drhob)","(drhoa_drhob)")
             n_deriv = n_deriv + 1
             deriv_data = 0.0_dp
           CASE("(tau)","(tau_a)","(tau_b)")
             n_deriv = n_deriv + 1
             deriv_data = 0.0_dp
         END SELECT
       ELSE
         SELECT CASE(desc)
           CASE("(rho)","(norm_drho)")
             n_deriv = n_deriv + 1
             deriv_data = 0.0_dp
           CASE("(tau)")
             n_deriv = n_deriv + 1
             deriv_data = 0.0_dp
         END SELECT
       END IF
     END IF
   END DO

   CALL timestop(0.0_dp,handle)

  END SUBROUTINE vxc_of_r_new

!   ***************************************************************************

  SUBROUTINE xc_rho_set_atom_update(rho_set,needs,nspins,na)

!   This routine allocates the storage arrays for rho and drho
!   In calculate_vxc_atom this is called once for each atomic_kind, 
!   After the loop over all the atoms of the kind and over all the points
!   of the radial grid for each atom, rho_set is deallocated.
!   Within the same kind, at each new point on the radial grid, the rho_set
!   arrays rho and drho are overwritten.

    TYPE(xc_rho_set_type), POINTER      :: rho_set
    TYPE(xc_rho_cflags_type)            :: needs
    INTEGER, INTENT(IN)                 :: nspins, na

    CHARACTER(LEN=*), PARAMETER                   :: routine_name = "vxc_of_r_new",&
                                         routineP = module_name//':'//routine_name

    INTEGER :: idir
!   ---------------------------------------------------------------------------

    SELECT CASE(nspins)
    CASE(1)
!     What is this for?
      IF (needs%rho_1_3) THEN
        NULLIFY(rho_set%rho_1_3)
        CALL reallocate(rho_set%rho_1_3,1,na,1,1,1,1)
        rho_set%owns%rho_1_3=.TRUE.
        rho_set%has%rho_1_3= .FALSE.
      END IF
!     Allocate the storage space for the density 
      IF (needs%rho) THEN
        NULLIFY(rho_set%rho)
        CALL reallocate(rho_set%rho,1,na,1,1,1,1)
        rho_set%owns%rho=.TRUE.
        rho_set%has%rho=.FALSE.
      END IF
!     Allocate the storage space for  the norm of the gradient of the density
      IF (needs%norm_drho) THEN
        NULLIFY(rho_set%norm_drho)
        CALL reallocate(rho_set%norm_drho,1,na,1,1,1,1) 
        rho_set%owns%norm_drho=.TRUE.
        rho_set%has%norm_drho=.FALSE.
      END IF
!     Allocate the storage space for the three components of the gradient of the density
      IF (needs%drho) THEN
        DO idir = 1,3
          NULLIFY (rho_set%drho(idir)%array)
          CALL reallocate(rho_set%drho(idir)%array,1,na,1,1,1,1)
        END DO
        rho_set%owns%drho=.TRUE.
        rho_set%has%drho=.FALSE.
      END IF
    CASE(2)
!     Allocate the storage space for the total density
      IF (needs%rho) THEN
     ! this should never be the case unless you use LDA functionals with LSD
        NULLIFY(rho_set%rho)
        CALL reallocate(rho_set%rho,1,na,1,1,1,1)
        rho_set%owns%rho=.TRUE.
        rho_set%has%rho=.FALSE.
      END IF
!     What is this for?
      IF (needs%rho_1_3) THEN
        NULLIFY(rho_set%rho_1_3)
        CALL reallocate(rho_set%rho_1_3,1,na,1,1,1,1) 
        rho_set%owns%rho_1_3=.TRUE.
        rho_set%has%rho_1_3=.FALSE.
      END IF
!     What is this for?
      IF (needs%rho_spin_1_3) THEN
        NULLIFY(rho_set%rhoa_1_3,rho_set%rhob_1_3)
        CALL reallocate(rho_set%rhoa_1_3,1,na,1,1,1,1)
        CALL reallocate(rho_set%rhob_1_3,1,na,1,1,1,1) 
        rho_set%owns%rho_spin_1_3=.TRUE.
        rho_set%has%rho_spin_1_3=.FALSE.
      END IF
!     Allocate the storage space for the spin densities rhoa and rhob
      IF (needs%rho_spin) THEN
        NULLIFY(rho_set%rhoa,rho_set%rhob)
        CALL reallocate(rho_set%rhoa,1,na,1,1,1,1)
        CALL reallocate(rho_set%rhob,1,na,1,1,1,1)
        rho_set%owns%rho_spin=.TRUE.
        rho_set%has%rho_spin=.FALSE.
      END IF
!     Allocate the storage space for the norm of the gradient of the total density
      IF (needs%norm_drho) THEN
        NULLIFY(rho_set%norm_drho)
        CALL reallocate(rho_set%norm_drho,1,na,1,1,1,1)
        rho_set%owns%norm_drho=.TRUE.
        rho_set%has%norm_drho=.FALSE.
      END IF
!     Allocate the storage space for the norm of the gradient of rhoa and of rhob separatedly
      IF (needs%norm_drho_spin) THEN
        NULLIFY(rho_set%norm_drhoa,rho_set%norm_drhob)
        CALL reallocate(rho_set%norm_drhoa,1,na,1,1,1,1)
        CALL reallocate(rho_set%norm_drhob,1,na,1,1,1,1)
        rho_set%owns%norm_drho_spin=.TRUE.
        rho_set%has%norm_drho_spin=.FALSE.
      END IF
!
      IF (needs%drhoa_drhob) THEN
        NULLIFY(rho_set%drhoa_drhob)
        CALL reallocate(rho_set%drhoa_drhob,1,na,1,1,1,1)
        rho_set%owns%drhoa_drhob=.TRUE.
        rho_set%has%drhoa_drhob=.FALSE.
      END IF
!     Allocate the storage space for the components of the gradient for the total rho
      IF (needs%drho) THEN
        DO idir = 1,3
          NULLIFY(rho_set%drho(idir)%array)
          CALL reallocate(rho_set%drho(idir)%array,1,na,1,1,1,1)
        END DO
        rho_set%owns%drho=.TRUE.        
        rho_set%has%drho=.FALSE.
      END IF
!     Allocate the storage space for the components of the gradient for rhoa and rhob
      IF (needs%drho_spin) THEN
        DO idir = 1,3
          NULLIFY(rho_set%drhoa(idir)%array,rho_set%drhob(idir)%array)
          CALL reallocate(rho_set%drhoa(idir)%array,1,na,1,1,1,1)
          CALL reallocate(rho_set%drhob(idir)%array,1,na,1,1,1,1)
        END DO
        rho_set%owns%drho_spin=.TRUE.        
        rho_set%has%drho_spin=.FALSE.
      END IF
!
    END SELECT

    ! tau part
    IF (needs%tau.OR.needs%tau_spin) THEN
       CALL stop_program(routine_name,module_name,__LINE__,&
                        "GAPW + metafunctionals not implemented yet")
    END IF
    IF (needs%tau) THEN
      NULLIFY(rho_set%tau)
      CALL reallocate(rho_set%tau,1,na,1,1,1,1)
      rho_set%owns%tau=.TRUE.
    END IF
    IF (needs%tau_spin) THEN
      NULLIFY(rho_set%tau_a,rho_set%tau_b)
      CALL reallocate(rho_set%tau_a,1,na,1,1,1,1)
      CALL reallocate(rho_set%tau_b,1,na,1,1,1,1)
      rho_set%owns%tau_spin=.TRUE.
      rho_set%has%tau_spin=.FALSE.
    END IF


  END SUBROUTINE xc_rho_set_atom_update

!   ***************************************************************************

  SUBROUTINE fill_rho_set(rho_set,lsd,tddft,needs,rho,drho,oodrho,w,na,ir,error)

    ! ARGUMENTS
    TYPE(xc_rho_set_type), POINTER                :: rho_set
    LOGICAL, INTENT(IN)                           :: lsd, tddft
    TYPE(xc_rho_cflags_type), INTENT(in)          :: needs
    REAL(dp), DIMENSION(:,:),  POINTER            :: rho 
    REAL(dp), DIMENSION(:,:,:,:),  POINTER        :: drho 
    REAL(dp), DIMENSION(:,:),  POINTER            :: oodrho 
    REAL(dp), DIMENSION(:),  POINTER              :: w 
    INTEGER, INTENT(IN)                           :: na, ir
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                    :: error

    ! LOCALS
    INTEGER                                       :: ia, idir, nspins
    LOGICAL                                       :: failure
    LOGICAL                                       :: gradient_f, tddft_split
    REAL(dp), DIMENSION(:,:), POINTER             :: tau

    ! PARAMETERS
    REAL(KIND=dp), PARAMETER                      :: f13 = (1.0_dp/3.0_dp)
    CHARACTER(LEN=*), PARAMETER                   :: routine_name = "fill_rho_set",&
                                                     routineP = module_name//':'//routine_name
        

    NULLIFY(tau)
    failure = .FALSE.
    tddft_split = tddft.AND.(.NOT.lsd)
    nspins = 1
    IF (lsd .OR. tddft_split) nspins = 2

    ! some checks
    IF (lsd) THEN
       CPPrecondition(SIZE(rho,2)==2,cp_failure_level,routineP,error,failure)
    ELSE
       CPPrecondition(SIZE(rho,2)==1,cp_failure_level,routineP,error,failure)
    END IF
    SELECT CASE(nspins)
    CASE(1)
       CPPrecondition(.NOT.needs%rho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%norm_drho_spin,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%drhoa_drhob,cp_failure_level,routineP,error,failure)
       CPPrecondition(.NOT.needs%rho_spin_1_3,cp_failure_level,routineP,error,failure)
    CASE(2)
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    gradient_f=(needs%drho_spin.OR.needs%norm_drho_spin.OR.&
                needs%drhoa_drhob.OR.needs%drho.OR.needs%norm_drho)
    
    SELECT CASE(nspins)
    CASE(1)
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
             rho_set%norm_drho(ia,1,1) = drho(4,ia,ir,1)
          END DO
          rho_set%owns%norm_drho=.TRUE.
          rho_set%has%norm_drho=.TRUE.
       END IF
       ! Give the three components of the gradient of the density
       IF (needs%drho) THEN
          DO idir = 1,3
             DO ia = 1,na
                rho_set%drho(idir)%array(ia,1,1) = drho(idir,ia,ir,1)
             END DO
          END DO
          rho_set%owns%drho=.TRUE.
          rho_set%has%drho=.TRUE.
       END IF
    CASE(2)
       ! Give the total density
       IF (needs%rho) THEN
          ! this should never be the case unless you use LDA functionals with LSD
          IF (.NOT.tddft_split) THEN
             DO ia = 1,na
                rho_set%rho(ia,1,1) = rho(ia,1)+rho(ia,2)
             END DO
          ELSE
             DO ia = 1,na
                rho_set%rho(ia,1,1) = rho(ia,1)
             END DO
          END IF
          rho_set%owns%rho=.TRUE.
          rho_set%has%rho=.TRUE.
       END IF
       ! Give the total density to 1/3  
       IF (needs%rho_1_3) THEN
          IF (.NOT.tddft_split) THEN
             DO ia = 1,na
                rho_set%rho_1_3(ia,1,1) = MAX(rho(ia,1)+rho(ia,2),0.0_dp)**f13
             END DO
          ELSE
             DO ia = 1,na
                rho_set%rho_1_3(ia,1,1) = MAX(rho(ia,1),0.0_dp)**f13
             END DO
          END IF
          rho_set%owns%rho_1_3=.TRUE.
          rho_set%has%rho_1_3=.TRUE.
       END IF
       ! Give the spin densities to 1/3  
       IF (needs%rho_spin_1_3) THEN
          IF (.NOT.tddft_split) THEN
             DO ia = 1,na
                rho_set%rhoa_1_3(ia,1,1) = MAX(rho(ia,1),0.0_dp)**f13
                rho_set%rhob_1_3(ia,1,1) = MAX(rho(ia,2),0.0_dp)**f13
             END DO
          ELSE
             DO ia = 1,na
                rho_set%rhoa_1_3(ia,1,1) = MAX(0.5_dp*rho(ia,1),0.0_dp)**f13
                rho_set%rhob_1_3(ia,1,1) = rho_set%rhoa_1_3(ia,1,1)
             END DO
          END IF
          rho_set%owns%rho_spin_1_3=.TRUE.
          rho_set%has%rho_spin_1_3=.TRUE.
       END IF
       ! Give the spin densities rhoa and rhob
       IF (needs%rho_spin) THEN
          IF (.NOT.tddft_split) THEN
             DO ia = 1,na
                rho_set%rhoa(ia,1,1) = rho(ia,1)
                rho_set%rhob(ia,1,1) = rho(ia,2)
             END DO
          ELSE
             DO ia = 1,na
                rho_set%rhoa(ia,1,1) = 0.5_dp*rho(ia,1)
                rho_set%rhob(ia,1,1) = rho_set%rhoa(ia,1,1)
             END DO
          END IF
          rho_set%owns%rho_spin=.TRUE.
          rho_set%has%rho_spin=.TRUE.
       END IF
       ! Give the norm of the gradient of the total density
       IF (needs%norm_drho) THEN
          IF (.NOT.tddft_split) THEN
             DO ia = 1,na
                rho_set%norm_drho(ia,1,1) = SQRT(&
                     (drho(1,ia,ir,1)+drho(1,ia,ir,2))**2+&
                     (drho(2,ia,ir,1)+drho(2,ia,ir,2))**2+&
                     (drho(3,ia,ir,1)+drho(3,ia,ir,2))**2)
             END DO
          ELSE
             DO ia = 1,na
                rho_set%norm_drho(ia,1,1) = SQRT(&
                     (drho(1,ia,ir,1))**2+&
                     (drho(2,ia,ir,1))**2+&
                     (drho(3,ia,ir,1))**2)
             END DO
          END IF
          rho_set%owns%norm_drho=.TRUE.
          rho_set%has%norm_drho=.TRUE.
       END IF
       ! Give the norm of the gradient of rhoa and of rhob separatedly
       IF (needs%norm_drho_spin) THEN
          IF (.NOT.tddft_split) THEN
             DO ia = 1,na
                rho_set%norm_drhoa(ia,1,1) = drho(4,ia,ir,1)
                rho_set%norm_drhob(ia,1,1) = drho(4,ia,ir,2)
             END DO
          ELSE
             DO ia = 1,na
                rho_set%norm_drhoa(ia,1,1) = 0.5_dp*drho(4,ia,ir,1)
                rho_set%norm_drhob(ia,1,1) = rho_set%norm_drhoa(ia,1,1)
             END DO
          END IF
          rho_set%owns%norm_drho_spin=.TRUE.
          rho_set%has%norm_drho_spin=.TRUE.
       END IF
       !
       IF (needs%drhoa_drhob) THEN
          IF (.NOT.tddft_split) THEN
             DO ia =1,na
                rho_set%drhoa_drhob(ia,1,1) =&
                     (drho(1,ia,ir,1)*drho(1,ia,ir,2))+&
                     (drho(2,ia,ir,1)*drho(2,ia,ir,2))+&
                     (drho(3,ia,ir,1)*drho(3,ia,ir,2))   
             END DO
          ELSE
             DO ia =1,na
                rho_set%drhoa_drhob(ia,1,1) =&
                     (0.25_dp*drho(1,ia,ir,1)**2)+&
                     (0.25_dp*drho(2,ia,ir,1)**2)+&
                     (0.25_dp*drho(3,ia,ir,1)**2)   
             END DO
          END IF
          rho_set%owns%drhoa_drhob=.TRUE.
          rho_set%has%drhoa_drhob=.TRUE.
       END IF
       ! Give the components of the gradient for the total rho
       IF (needs%drho) THEN
          IF (.NOT.tddft_split) THEN
             DO idir = 1,3
                DO ia = 1,na
                   rho_set%drho(idir)%array(ia,1,1) = drho(idir,ia,ir,1)+drho(idir,ia,ir,2)
                END DO
             END DO
          ELSE
             DO idir = 1,3
                DO ia = 1,na
                   rho_set%drho(idir)%array(ia,1,1) = drho(idir,ia,ir,1)
                END DO
             END DO
          END IF
          rho_set%owns%drho=.TRUE.        
          rho_set%has%drho=.TRUE.
       END IF
       ! Give the components of the gradient for rhoa and rhob
       IF (needs%drho_spin) THEN
          IF (.NOT.tddft_split) THEN
             DO idir = 1,3
                DO ia = 1,na
                   rho_set%drhoa(idir)%array(ia,1,1) = drho(idir,ia,ir,1)
                   rho_set%drhob(idir)%array(ia,1,1) = drho(idir,ia,ir,2)
                END DO
             END DO
          ELSE
             DO idir = 1,3
                DO ia = 1,na
                   rho_set%drhoa(idir)%array(ia,1,1) = 0.5_dp*drho(idir,ia,ir,1)
                   rho_set%drhob(idir)%array(ia,1,1) = rho_set%drhoa(idir)%array(ia,1,1)
                END DO
             END DO
          END IF
          rho_set%owns%drho_spin=.TRUE.        
          rho_set%has%drho_spin=.TRUE.
       END IF
       !
    END SELECT

    ! tau part
    IF (needs%tau.OR.needs%tau_spin) THEN
       CALL stop_program(routine_name,module_name,__LINE__,&
            "GAPW + metafunctionals not implemented yet")
       CPPrecondition(ASSOCIATED(tau),cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(tau,2)==nspins,cp_failure_level,routineP,error,failure)
    END IF
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
       DO ia = 1,na
          rho_set%tau_a(ia,1,1) = tau(ia,1)
          rho_set%tau_b(ia,1,1) = tau(ia,2)
       END DO
       rho_set%owns%tau_spin=.TRUE.
       rho_set%has%tau_spin=.TRUE.
    END IF


  END SUBROUTINE fill_rho_set

END MODULE xc_atom
