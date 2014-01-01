!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations for DFTB
!> \par History
!>      Started as a copy from the GPW file
!> \author JHU (03.2013)
! *****************************************************************************
MODULE qs_scf_post_dftb
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: twopi
  USE moments_utils,                   ONLY: get_reference_point
  USE mulliken,                        ONLY: mulliken_charges
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE qs_dftb_types,                   ONLY: qs_dftb_atom_type
  USE qs_dftb_utils,                   ONLY: get_dftb_atom_param
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_dftb'
  PUBLIC :: scf_post_calculation_dftb

CONTAINS

! *****************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!> \note
!> \par History
!>      03.2013 copy of scf_post_gpw
!> \author JHU
! *****************************************************************************
  SUBROUTINE scf_post_calculation_dftb(dft_section, scf_env, qs_env, error)

    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_dftb', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: ana
    CHARACTER(LEN=default_string_length)     :: aname
    INTEGER                                  :: handle, iat, iatom, ikind, &
                                                istat, nat, natom, nkind, &
                                                nspins, unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: zeff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: mcharge
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: charges
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_dftb_atom_type), POINTER         :: dftb_kind
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input, moments_section, &
                                                print_key

    CALL timeset(routineN,handle)

    failure=.FALSE.

    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    ! Here we start with data that needs a postprocessing...
    IF (.NOT. failure) THEN
       NULLIFY(dft_control,input,rho,para_env,matrix_s,matrix_p)
       CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
            dft_control=dft_control,input=input,rho=rho,&
            natom=natom,para_env=para_env,matrix_s=matrix_s,error=error)
       nspins = dft_control%nspins 
       matrix_p => rho%rho_ao
       ! Mulliken charges
       ALLOCATE(charges(natom,nspins),mcharge(natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       !
       CALL mulliken_charges(matrix_p,matrix_s(1)%matrix,para_env,charges,error=error)
       !
       nkind = SIZE(atomic_kind_set)
       DO ikind=1,nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,dftb_parameter=dftb_kind)
          CALL get_dftb_atom_param(dftb_kind,zeff=zeff)
          DO iatom=1,nat
             iat = atomic_kind%atom_list(iatom)
             mcharge(iat) = zeff - SUM(charges(iat,1:nspins))
          END DO
       END DO

       print_key => section_vals_get_subs_vals(input,"DFT%PRINT%MULLIKEN", error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MULLIKEN",extension=".mulliken",&
                                       middle_name="",log_filename=.FALSE.,error=error)
          IF (unit_nr > 0) THEN
             WRITE (UNIT=unit_nr,FMT="(/,/,T2,A)") "MULLIKEN POPULATION ANALYSIS"
             IF (nspins == 1) THEN
                WRITE (UNIT=unit_nr,FMT="(/,T2,A,T70,A)")&
                       " # Atom   Element   Kind        Atomic population"," Net charge"
                DO ikind=1,nkind
                   atomic_kind => atomic_kind_set(ikind)
                   CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,dftb_parameter=dftb_kind)
                   CALL get_dftb_atom_param(dftb_kind,name=aname)
                   ana = ADJUSTR(TRIM(ADJUSTL(aname)))
                   DO iatom=1,nat
                      iat = atomic_kind%atom_list(iatom)
                      WRITE (UNIT=unit_nr,&
                             FMT="(T2,I7,6X,A2,3X,I6,T39,F12.6,T69,F12.6)")&
                             iat,ana,ikind,charges(iat,1),mcharge(iat)
                   END DO
                END DO
                WRITE (UNIT=unit_nr,&
                       FMT="(T2,A,T39,F12.6,T69,F12.6,/)")&
                           "# Total charge",SUM(charges(:,1)),SUM(mcharge(:))
             ELSE
                WRITE (UNIT=unit_nr,FMT="(/,T2,A)")&
                       "# Atom  Element  Kind  Atomic population (alpha,beta)   Net charge  Spin moment"
                DO ikind=1,nkind
                   atomic_kind => atomic_kind_set(ikind)
                   CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,dftb_parameter=dftb_kind)
                   CALL get_dftb_atom_param(dftb_kind,name=aname)
                   ana = ADJUSTR(TRIM(ADJUSTL(aname)))
                   DO iatom=1,nat
                      iat = atomic_kind%atom_list(iatom)
                      WRITE (UNIT=unit_nr,&
                             FMT="(T2,I6,5X,A2,2X,I6,T29,4(1X,F12.6))")&
                             iat,ana,ikind,charges(iat,1:2),mcharge(iat),charges(iat,1)-charges(iat,2)
                   END DO
                END DO
                WRITE (UNIT=unit_nr,&
                       FMT="(T2,A,T29,4(1X,F12.6),/)")&
                           "# Total charge and spin",SUM(charges(:,1)),SUM(charges(:,2)),SUM(mcharge(:))
             END IF
             CALL m_flush(unit_nr)
          END IF
          CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)
       END IF

       ! Dipole Moments
       print_key => section_vals_get_subs_vals(input,"DFT%PRINT%MOMENTS", error=error)
       IF(BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file))THEN
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MOMENTS",&
               extension=".data",middle_name="dftb_dipole",log_filename=.FALSE.,error=error)
          moments_section => section_vals_get_subs_vals(input,"DFT%PRINT%MOMENTS",error=error)
          CALL dftb_dipole(qs_env, moments_section, unit_nr, mcharge, error)
          CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)
       END IF

       DEALLOCATE(charges,mcharge,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

    END IF
    CALL timestop(handle)

  END SUBROUTINE scf_post_calculation_dftb

! *****************************************************************************

  SUBROUTINE dftb_dipole(qs_env, input, unit_nr, charges, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: input
    INTEGER, INTENT(in)                      :: unit_nr
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: charges
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dftb_dipole', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description, dipole_type
    COMPLEX(KIND=dp)                         :: dzeta, dzphase(3), zeta, &
                                                zphase(3)
    COMPLEX(KIND=dp), DIMENSION(3)           :: dggamma, ggamma
    INTEGER                                  :: i, iat, ikind, j, nat, &
                                                reference
    LOGICAL                                  :: do_berry
    REAL(KIND=dp) :: charge_tot, ci(3), dci(3), dipole(3), dipole_deriv(3), &
      drcc(3), dria(3), dtheta, gvec(3), q, rcc(3), ria(3), theta, tmp(3), &
      via(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ref_point
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    NULLIFY(atomic_kind_set,cell)
    CALL get_qs_env(qs_env,atomic_kind_set=atomic_kind_set,&
                    particle_set=particle_set,cell=cell,error=error)

    ! Reference point
    reference = section_get_ival(input,keyword_name="REFERENCE",error=error)
    NULLIFY(ref_point)
    description='[DIPOLE]'
    CALL section_vals_val_get(input,"REF_POINT",r_vals=ref_point,error=error)
    CALL section_vals_val_get(input,"PERIODIC",l_val=do_berry,error=error)

    CALL get_reference_point(rcc,drcc,qs_env=qs_env,reference=reference,ref_point=ref_point,error=error)

    ! Dipole deriv will be the derivative of the Dipole(dM/dt=\sum e_j v_j)
    dipole_deriv = 0.0_dp
    dipole       = 0.0_dp
    IF(do_berry)THEN
       dipole_type="[BERRY PHASE]"
       rcc = pbc(rcc,cell)
       charge_tot = 0._dp
       charge_tot=SUM(charges)
       ria  = twopi * MATMUL(cell%h_inv, rcc)
       zphase  = CMPLX(COS(ria), SIN(ria), dp)**charge_tot

       dria = twopi * MATMUL(cell%h_inv, drcc)
       dzphase = charge_tot * CMPLX(-SIN(ria), COS(ria), dp)**(charge_tot-1.0_dp)* dria

       ggamma  = CMPLX(1.0_dp,0.0_dp,KIND=dp)
       dggamma = CMPLX(0.0_dp,0.0_dp,KIND=dp)
       DO ikind = 1, SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat)
          DO i = 1, nat
             iat = atomic_kind%atom_list(i)
             ria = particle_set(iat)%r(:)
             ria = pbc(ria,cell)
             via = particle_set(iat)%v(:)
             q=charges(iat)
             DO j = 1, 3
                gvec  = twopi*cell%h_inv(j,:)
                theta = SUM(ria(:)*gvec(:))
                dtheta= SUM(via(:)*gvec(:))
                zeta  = CMPLX(COS(theta),SIN(theta),KIND=dp)**(-q)
                dzeta = -q * CMPLX(-SIN(theta),COS(theta),KIND=dp)**(-q-1.0_dp)* dtheta
                dggamma(j)= dggamma(j)* zeta + ggamma(j)* dzeta
                ggamma(j) = ggamma(j)* zeta
             END DO
          ENDDO
       END DO
       dggamma = dggamma * zphase + ggamma * dzphase
       ggamma = ggamma * zphase
       IF(ALL(REAL(ggamma,KIND=dp)/= 0.0_dp))THEN
          tmp = AIMAG(ggamma)/REAL(ggamma,KIND=dp)
          ci  = ATAN(tmp)
          dci =(1.0_dp /(1.0_dp + tmp**2))* &
              (AIMAG(dggamma)*REAL(ggamma,KIND=dp)-AIMAG(ggamma)*REAL(dggamma,KIND=dp))/(REAL(ggamma,KIND=dp))**2
          dipole       = MATMUL(cell%hmat, ci)/twopi
          dipole_deriv = MATMUL(cell%hmat, dci)/twopi
       END IF
    ELSE
       dipole_type="[Non Periodic]"
       DO i = 1, SIZE(particle_set)
          ! no pbc(particle_set(i)%r(:),cell) so that the total dipole is the sum of the molecular dipoles
          ria = particle_set(i)%r(:) 
          q=charges(i)
          dipole = dipole - q *(ria-rcc)
          dipole_deriv(:)= dipole_deriv(:)- q *(particle_set(i)%v(:)- drcc)
       END DO
    END IF
    IF(unit_nr>0)THEN
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"DIPOLE "//TRIM(dipole_type)//"(A.U.)|",dipole
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"DIPOLE "//TRIM(dipole_type)//"(Debye)|",dipole*debye
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"DIPOLE "//TRIM(dipole_type)//" DERIVATIVE(A.U.)|",dipole_deriv
    END IF

  END SUBROUTINE dftb_dipole

! *****************************************************************************

END MODULE qs_scf_post_dftb
