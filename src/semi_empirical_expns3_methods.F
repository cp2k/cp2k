!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Methods for handling the 1/R^3 residual integral part
!> \author Teodoro Laino (12.2008) [tlaino]
! *****************************************************************************
MODULE semi_empirical_expns3_methods
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cp_control_types,                ONLY: semi_empirical_control_type
  USE f77_blas
  USE input_constants,                 ONLY: do_method_undef
  USE kinds,                           ONLY: dp
  USE semi_empirical_expns3_types,     ONLY: semi_empirical_expns3_create
  USE semi_empirical_int3_utils,       ONLY: coeff_int_3,&
                                             ijkl_low_3
  USE semi_empirical_int_arrays,       ONLY: indexa,&
                                             l_index
  USE semi_empirical_types,            ONLY: semi_empirical_type
  USE semi_empirical_utils,            ONLY: get_se_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PARAMETER, PRIVATE          :: debug_this_module   = .FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'semi_empirical_expns3_methods'

  PUBLIC :: semi_empirical_expns3_setup

CONTAINS
! *****************************************************************************
!> \brief Setup the quantity necessary to handle the slowly convergent
!>        residual integral term 1/R^3
!>
!> \author Teodoro Laino [tlaino]
!> \date 12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE semi_empirical_expns3_setup(atomic_kind_set, se_control, method_id, error)
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(semi_empirical_control_type), &
      POINTER                                :: se_control
    INTEGER, INTENT(IN)                      :: method_id
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'semi_empirical_expns3_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, itype, j, nkinds, stat
    LOGICAL                                  :: check, failure
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_i, atomic_kind_j
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj

    failure     = .FALSE.
    IF (se_control%do_ewald_r3) THEN
       NULLIFY(sepi, sepj, atomic_kind_i, atomic_kind_j)
       nkinds  = SIZE(atomic_kind_set)
       DO i = 1, nkinds
          atomic_kind_i => atomic_kind_set(i)
          CALL get_atomic_kind(atomic_kind_i, se_parameter=sepi)
          check = .NOT.ASSOCIATED(sepi%expns3_int)
          CPPostcondition(check,cp_failure_level,routineP,error,failure)
          ALLOCATE(sepi%expns3_int(nkinds),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO j = 1, nkinds
             NULLIFY(sepi%expns3_int(j)%expns3)
             CALL semi_empirical_expns3_create(sepi%expns3_int(j)%expns3,error=error)
          END DO
       END DO

       itype = get_se_type(method_id)
       DO i = 1, nkinds
          atomic_kind_i => atomic_kind_set(i)
          CALL get_atomic_kind(atomic_kind_i, se_parameter=sepi)
          DO j = 1, nkinds
             atomic_kind_j => atomic_kind_set(j)
             CALL get_atomic_kind(atomic_kind_j, se_parameter=sepj)
             CALL setup_c3_coeff(sepi, sepj, i, j, itype, error)
          END DO
       END DO
    END IF
  END SUBROUTINE semi_empirical_expns3_setup

! *****************************************************************************
!> \brief For any given semi-empirical pair i,j evaluates the coefficient of
!>        the integral residual part ( 1/r^3 term )
!>        The integral expression, unfortunately, does not allow any kind of
!>        separability. It is, therefore, mandatory to compute this coefficient
!>        as a pair term, instead of as an atomic quantity.
!>
!> \author Teodoro Laino [tlaino]
!> \date 12.2008 [tlaino]
! *****************************************************************************
  SUBROUTINE setup_c3_coeff(sepi, sepj, ikind, jkind, itype, error)
    TYPE(semi_empirical_type), POINTER       :: sepi, sepj
    INTEGER, INTENT(IN)                      :: ikind, jkind, itype
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'setup_c3_coeff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ij, j, kl, kr, li, lk
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: core_core, e1b(9), e2a(9), r, &
                                                zi, zj

    failure = .FALSE.
    ! Set the distance to 0 (the coefficient is anyway independent of the atomic
    ! position)
    r   = 0.0_dp
    ! Nuclei-Nuclei     contribution
    ij  = indexa(1, 1)
    zi  = -sepi%zeff
    zj  = -sepj%zeff
    core_core = ijkl_low_3(sepi, sepj, ij, ij, 0, 0, 0, 0, -1, r, itype, coeff_int_3, error) * zi * zj

    ! Electron(i)-Nuclei(j)   contribution
    kl = indexa(1, 1)
    e1b(1) = ijkl_low_3 (sepi, sepj, kl, ij, 0, 0, 0, 0, 2, r, itype, coeff_int_3, error) * zj
    IF (sepi%natorb > 1) THEN
       kl = indexa(2, 2)
       e1b(2) = ijkl_low_3 (sepi, sepj, kl, ij, 1, 1, 0, 0, 2, r, itype, coeff_int_3, error) * zj
       kl = indexa(3, 3)
       e1b(3) = ijkl_low_3 (sepi, sepj, kl, ij, 1, 1, 0, 0, 2, r, itype, coeff_int_3, error) * zj
       kl = indexa(4, 4)
       e1b(4) = ijkl_low_3 (sepi, sepj, kl, ij, 1, 1, 0, 0, 2, r, itype, coeff_int_3, error) * zj
       ! Consistency check
       CPPostcondition(e1b(2)==e1b(3),cp_failure_level,routineP,error,failure)
       CPPostcondition(e1b(3)==e1b(4),cp_failure_level,routineP,error,failure)
       IF (sepi%dorb) THEN
          kl = indexa(5, 5)
          e1b(5) = ijkl_low_3 (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, r, itype, coeff_int_3, error) * zj
          kl = indexa(6, 6)
          e1b(6) = ijkl_low_3 (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, r, itype, coeff_int_3, error) * zj
          kl = indexa(7, 7)
          e1b(7) = ijkl_low_3 (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, r, itype, coeff_int_3, error) * zj
          kl = indexa(8, 8)
          e1b(8) = ijkl_low_3 (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, r, itype, coeff_int_3, error) * zj
          kl = indexa(9, 9)
          e1b(9) = ijkl_low_3 (sepi, sepj, kl, ij, 2, 2, 0, 0, 2, r, itype, coeff_int_3, error) * zj
          ! Consistency check
          CPPostcondition(e1b(5)==e1b(6),cp_failure_level,routineP,error,failure)
          CPPostcondition(e1b(6)==e1b(7),cp_failure_level,routineP,error,failure)
          CPPostcondition(e1b(7)==e1b(8),cp_failure_level,routineP,error,failure)
          CPPostcondition(e1b(8)==e1b(9),cp_failure_level,routineP,error,failure)
       END IF
    END IF

    ! Electron(j)-Nuclei(i)   contribution
    kl = indexa(1, 1)
    e2a(1) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 0, 0, 1, r, itype, coeff_int_3, error) * zi
    IF (sepj%natorb > 1) THEN
       kl = indexa(2, 2)
       e2a(2) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 1, 1, 1, r, itype, coeff_int_3, error) * zi
       kl = indexa(3, 3)
       e2a(3) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 1, 1, 1, r, itype, coeff_int_3, error) * zi
       kl = indexa(4, 4)
       e2a(4) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 1, 1, 1, r, itype, coeff_int_3, error) * zi
       ! Consistency check
       CPPostcondition(e2a(2)==e2a(3),cp_failure_level,routineP,error,failure)
       CPPostcondition(e2a(3)==e2a(4),cp_failure_level,routineP,error,failure)
       IF (sepj%dorb) THEN
          kl = indexa(5, 5)
          e2a(5) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, r, itype, coeff_int_3, error) * zi
          kl = indexa(6, 6)
          e2a(6) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, r, itype, coeff_int_3, error) * zi
          kl = indexa(7, 7)
          e2a(7) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, r, itype, coeff_int_3, error) * zi
          kl = indexa(8, 8)
          e2a(8) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, r, itype, coeff_int_3, error) * zi
          kl = indexa(9, 9)
          e2a(9) = ijkl_low_3 (sepi, sepj, ij, kl, 0, 0, 2, 2, 1, r, itype, coeff_int_3, error) * zi
          ! Consistency check
          CPPostcondition(e2a(5)==e2a(6),cp_failure_level,routineP,error,failure)
          CPPostcondition(e2a(6)==e2a(7),cp_failure_level,routineP,error,failure)
          CPPostcondition(e2a(7)==e2a(8),cp_failure_level,routineP,error,failure)
          CPPostcondition(e2a(8)==e2a(9),cp_failure_level,routineP,error,failure)
       END IF
    END IF

    ! Copy info into the semi-empirical type (i)
    sepi%expns3_int(jkind)%expns3%core_core           = core_core
    sepi%expns3_int(jkind)%expns3%e1b(1:sepi%natorb)  = e1b(1:sepi%natorb)
    sepi%expns3_int(jkind)%expns3%e2a(1:sepj%natorb)  = e2a(1:sepj%natorb)
    ! Copy info into the semi-empirical type (j)
    sepj%expns3_int(ikind)%expns3%core_core           = core_core
    sepj%expns3_int(ikind)%expns3%e1b(1:sepj%natorb)  = e2a(1:sepj%natorb)
    sepj%expns3_int(ikind)%expns3%e2a(1:sepi%natorb)  = e1b(1:sepi%natorb)

    ! Electron-Electron contribution - sepi/sepj
    kr = 0
    DO i = 1, sepi%natorb
       li = l_index(i)
       ij = indexa(i,i)
       DO j = 1, sepj%natorb
          lk = l_index(j)
          kl = indexa(j,j)
          kr = kr + 1
          sepi%expns3_int(jkind)%expns3%w(kr) = &
               ijkl_low_3 (sepi, sepj, ij, kl, li, li, lk, lk, 0, r, do_method_undef, coeff_int_3, error)
       END DO
    END DO

    ! Electron-Electron contribution - sepj/sepi
    kr = 0
    DO i = 1, sepj%natorb
       li = l_index(i)
       ij = indexa(i,i)
       DO j = 1, sepi%natorb
          lk = l_index(j)
          kl = indexa(j,j)
          kr = kr + 1
          sepj%expns3_int(ikind)%expns3%w(kr) = &
               ijkl_low_3 (sepj, sepi, ij, kl, li, li, lk, lk, 0, r, do_method_undef, coeff_int_3, error)
       END DO
    END DO

  END SUBROUTINE setup_c3_coeff

END MODULE semi_empirical_expns3_methods
