!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief contains the master routines for the Hartree-Fock RI Method   
!> \author JGH July 2009    
! *****************************************************************************
MODULE hfx_ri_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE basis_set_types,                 ONLY: geminal_basis_set_type,&
                                             get_geminal_basis_set
  USE kinds,                           ONLY: dp
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_geminals,                     ONLY: geminal_coulomb,&
                                             geminal_gto_coulomb
  USE ri_environment_types,            ONLY: &
       get_ri_env, ri_env_create, ri_environment_type, ri_vector_add, &
       ri_vector_copy, ri_vector_create, ri_vector_dot, ri_vector_norm, &
       ri_vector_release, ri_vector_scale, ri_vector_type, set_ri_env
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *****************************************************************************

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri_methods'

  PUBLIC :: hfx_ri_energy_potential, hfx_ri_env_create

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief creates and initializes an ri_env
!> \param ri_env the ri_environment you want to create
!>
!> \author JGH
! *****************************************************************************
  SUBROUTINE hfx_ri_env_create (qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, istat, nat, &
                                                ngem, nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natoms, nbasis
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(geminal_basis_set_type), POINTER    :: geminal_basis_set
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_coeff, ri_rhs

    CALL timeset(routineN,handle)

    CALL ri_env_create(ri_env,error)

    CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,error=error)
    nkind = SIZE(atomic_kind_set)
    ALLOCATE(nbasis(nkind),natoms(nkind),stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    nbasis=0
    natoms=0
    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       NULLIFY(geminal_basis_set)
       CALL get_atomic_kind(atomic_kind=atomic_kind,natom=nat,geminal_basis_set=geminal_basis_set)
       IF (.NOT.ASSOCIATED(geminal_basis_set)) CYCLE
       CALL get_geminal_basis_set(geminal_basis_set=geminal_basis_set,ngeminals=ngem)
       natoms(ikind) = nat
       nbasis(ikind) = ngem
    END DO

    NULLIFY(ri_coeff,ri_rhs)
    CALL ri_vector_create(ri_coeff, nkind, natoms, nbasis, error)
    CALL ri_vector_create(ri_rhs, nkind, natoms, nbasis, error)

    CALL set_ri_env(ri_env,coeff=ri_coeff,rhs=ri_rhs)
    CALL set_qs_env(qs_env=qs_env,hfx_ri_env=ri_env,error=error)

    DEALLOCATE(nbasis,natoms,stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    CALL timestop(handle)

  END SUBROUTINE hfx_ri_env_create

! *****************************************************************************
!> \brief creates and initializes an ri_env
!> \param ri_env the ri_environment you want to create
!> 
!> \author JGH
! *****************************************************************************
  SUBROUTINE hfx_ri_energy_potential (qs_env,calculate_force,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_energy_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: calculate_virial
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)

    CALL optimize_ri_coeff (qs_env,error)

    CALL get_qs_env(qs_env=qs_env,virial=virial,error=error)
    calculate_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)

    CALL ri_energy_potential(qs_env,calculate_force,calculate_virial,error)

    CALL timestop(handle)

  END SUBROUTINE hfx_ri_energy_potential

! *****************************************************************************

  SUBROUTINE optimize_ri_coeff (qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'optimize_ri_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, istat, iter, &
                                                niter, nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natoms, nbasis
    LOGICAL                                  :: converged
    REAL(KIND=dp)                            :: alpha, beta, energy, eps_ri, &
                                                papk, rnorm, rrk, rro
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_apk, ri_coeff, ri_pk, &
                                                ri_res, ri_rhs

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,hfx_ri_env=ri_env,error=error)
    CALL get_ri_env(ri_env=ri_env,coeff=ri_coeff,rhs=ri_rhs)

    ! calculate rhs of linear equation
    CALL geminal_gto_coulomb(qs_env,.FALSE.,.FALSE.,.FALSE.,.FALSE.,error)

    !allocate scratch vectors
    nkind = SIZE(ri_coeff%vector)
    ALLOCATE(nbasis(nkind),natoms(nkind),stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    DO ikind=1,nkind
      natoms(ikind) = ri_coeff%vector(ikind)%natom
      nbasis(ikind) = ri_coeff%vector(ikind)%nbasis
    END DO

    NULLIFY(ri_res,ri_pk,ri_apk)
    CALL ri_vector_create(ri_res, nkind, natoms, nbasis, error)
    CALL ri_vector_create(ri_pk, nkind, natoms, nbasis, error)
    CALL ri_vector_create(ri_apk, nkind, natoms, nbasis, error)

    DEALLOCATE(nbasis,natoms,stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    niter = 20
    eps_ri = 1.e-5_dp

    converged = .FALSE.
    ! calculate residual
    CALL geminal_coulomb(ri_coeff,ri_res,qs_env,energy,.FALSE.,.FALSE.,.FALSE.,error)
    CALL ri_vector_add(-1._dp, ri_rhs, ri_res, error)
    CALL ri_vector_scale(-1._dp, ri_res, error)
    CALL ri_vector_copy(ri_res, ri_pk, error)
    CALL ri_vector_dot(rrk, ri_res, ri_res, error)
    DO iter=1,niter
       CALL geminal_coulomb(ri_pk,ri_apk,qs_env,energy,.FALSE.,.FALSE.,.FALSE.,error)
       CALL ri_vector_dot(papk, ri_pk, ri_apk, error)
       alpha = rrk/papk
       CALL ri_vector_add(alpha, ri_pk, ri_coeff, error)
       CALL ri_vector_add(-alpha, ri_apk, ri_res, error)
       CALL ri_vector_norm(rnorm, ri_res, error)
       IF ( rnorm < eps_ri ) THEN
          converged = .TRUE.
          EXIT
       END IF
       rro = rrk
       CALL ri_vector_dot(rrk, ri_res, ri_res, error)
       beta = rrk/rro
       CALL ri_vector_scale(beta, ri_pk, error)
       CALL ri_vector_add(1._dp, ri_res, ri_pk, error)
    END DO

    CALL ri_vector_release(ri_res, error)
    CALL ri_vector_release(ri_pk, error)
    CALL ri_vector_release(ri_apk, error)

    CALL timestop(handle)

  END SUBROUTINE optimize_ri_coeff

! *****************************************************************************

  SUBROUTINE ri_energy_potential (qs_env,calculate_force,calculate_virial,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: calculate_force, &
                                                calculate_virial
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ri_energy_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, istat, nkind
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natoms, nbasis
    LOGICAL                                  :: calculate_energy, &
                                                calculate_fock
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(ri_environment_type), POINTER       :: ri_env
    TYPE(ri_vector_type), POINTER            :: ri_coeff, ri_res, ri_rhs

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,hfx_ri_env=ri_env,energy=energy,error=error)
    CALL get_ri_env(ri_env=ri_env,coeff=ri_coeff,rhs=ri_rhs)
    nkind = SIZE(ri_coeff%vector)
    ALLOCATE(nbasis(nkind),natoms(nkind),stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
    DO ikind=1,nkind
      natoms(ikind) = ri_coeff%vector(ikind)%natom
      nbasis(ikind) = ri_coeff%vector(ikind)%nbasis
    END DO

    NULLIFY(ri_res)
    CALL ri_vector_create(ri_res, nkind, natoms, nbasis, error)

    DEALLOCATE(nbasis,natoms,stat=istat)
    CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)

    calculate_energy = .TRUE.
    calculate_fock = .TRUE.

    CALL geminal_coulomb(ri_coeff,ri_res,qs_env,ri_env%ehfx1,calculate_energy,calculate_force,&
                         calculate_virial,error)
    CALL geminal_gto_coulomb(qs_env,calculate_fock,calculate_energy,&
                             calculate_force,calculate_virial,error)

    energy%ex = 0.5_dp*(ri_env%ehfx1 + ri_env%ehfx2)

    CALL ri_vector_release(ri_res, error)

    CALL timestop(handle)

  END SUBROUTINE ri_energy_potential

! *****************************************************************************

END MODULE hfx_ri_methods

