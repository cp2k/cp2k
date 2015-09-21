!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Debug the derivatives of the the rotational matrices
!>
!> \param sepi ...
!> \param sepj ...
!> \param rjiv ...
!> \param ij_matrix ...
!> \param do_invert ...
!> \date 04.2008 [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
SUBROUTINE check_rotmat_der( sepi, sepj, rjiv, ij_matrix, do_invert)
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_utils,        ONLY: rotmat
  USE semi_empirical_types,            ONLY: rotmat_create,&
                                             rotmat_release,&
                                             rotmat_type,&
                                             semi_empirical_type,&
                                             se_int_control_type
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(semi_empirical_type), POINTER       :: sepi, sepj
  REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rjiv
  TYPE(rotmat_type), POINTER               :: ij_matrix
  LOGICAL, INTENT(IN)                      :: do_invert

  CHARACTER(len=*), PARAMETER :: moduleN='semi_empirical_int_debug',&
       routineN = 'check_rotmat_der', routineP = moduleN//':'//routineN

  LOGICAL                                  :: failure, check_value
  REAL(KIND=dp)                            :: dx, r, r0(3), x(3)
  TYPE(rotmat_type), POINTER               :: matrix, matrix_m, matrix_n, &
       matrix_p
  INTEGER                                  :: imap(3),i,j,k,l

  NULLIFY(matrix_m, matrix_p, matrix_n, matrix)
  failure = .FALSE.
  CALL rotmat_create(matrix_p)
  CALL rotmat_create(matrix_m)
  CALL rotmat_create(matrix_n)
  dx = 1.0E-6_dp
  imap(1) = 1
  imap(2) = 2
  imap(3) = 3
  IF (do_invert) THEN
     imap(1) = 3
     imap(2) = 2
     imap(3) = 1
  END IF
  ! Check derivatives: analytical VS numerical
  WRITE(*,*)"DEBUG::"//routineP
  DO j = 1, 3
     x          = 0.0_dp
     x(imap(j)) = dx
     DO i = 1, 2
        IF (i==1) matrix => matrix_p
        IF (i==2) matrix => matrix_m
        r0 = rjiv + (-1.0_dp)**(i-1) * x
        r  = SQRT(DOT_PRODUCT(r0,r0))
        CALL rotmat (sepi, sepj, r0, r, matrix, do_derivatives=.FALSE., debug_invert=do_invert)
     END DO
     ! SP
     matrix_n%sp_d(j, :, :) = (matrix_p%sp - matrix_m%sp)/(2.0_dp*dx)
     DO i = 1, 3
        DO k = 1, 3
           IF (.NOT.check_value(matrix_n%sp_d(j,k,i), ij_matrix%sp_d(j,k,i), dx, 0.1_dp)) THEN
              WRITE(*,*)"ERROR for SP rotation matrix derivative SP(j,k,i), j,k,i::",j,k,i
              CPABORT("")
           END IF
        END DO
     END DO
     ! PP
     matrix_n%pp_d(j, :, :, :) = (matrix_p%pp - matrix_m%pp)/(2.0_dp*dx)
     DO i = 1, 3
        DO k = 1, 3
           DO l = 1, 6
              IF (.NOT.check_value(matrix_n%pp_d(j,l,k,i), ij_matrix%pp_d(j,l,k,i), dx, 0.1_dp)) THEN
                 WRITE(*,*)"ERROR for PP rotation matrix derivative PP(j,l,k,i), j,l,k,i::",j,l,k,i
                 CPABORT("")
              END IF
           END DO
        END DO
     END DO
     ! d-orbitals debug
     IF (sepi%dorb.OR.sepj%dorb) THEN
        ! SD
        matrix_n%sd_d(j, :, :) = (matrix_p%sd - matrix_m%sd)/(2.0_dp*dx)
        DO i = 1, 5
           DO k = 1, 5
              IF (.NOT.check_value(matrix_n%sd_d(j,k,i), ij_matrix%sd_d(j,k,i), dx, 0.1_dp)) THEN
                 WRITE(*,*)"ERROR for SD rotation matrix derivative SD(j,k,i), j,k,i::",j,k,i
                 CPABORT("")
              END IF
           END DO
        END DO
        ! DP
        matrix_n%pd_d(j, :, :, :) = (matrix_p%pd - matrix_m%pd)/(2.0_dp*dx)
        DO i = 1, 3
           DO k = 1, 5
              DO l = 1, 15
                 IF (.NOT.check_value(matrix_n%pd_d(j,l,k,i), ij_matrix%pd_d(j,l,k,i), dx, 0.1_dp)) THEN
                    WRITE(*,*)"ERROR for DP rotation matrix derivative DP(j,l,k,i), j,l,k,i::",j,l,k,i
                    CPABORT("")
                 END IF
              END DO
           END DO
        END DO
        ! DD
        matrix_n%dd_d(j, :, :, :) = (matrix_p%dd - matrix_m%dd)/(2.0_dp*dx)
        DO i = 1, 5
           DO k = 1, 5
              DO l = 1, 15
                 IF (.NOT.check_value(matrix_n%dd_d(j,l,k,i), ij_matrix%dd_d(j,l,k,i), dx, 0.1_dp)) THEN
                    WRITE(*,*)"ERROR for DD rotation matrix derivative DD(j,l,k,i), j,l,k,i::",j,l,k,i
                    CPABORT("")
                 END IF
              END DO
           END DO
        END DO
     END IF
  END DO
  CALL rotmat_release(matrix_p)
  CALL rotmat_release(matrix_m)
  CALL rotmat_release(matrix_n)
END SUBROUTINE check_rotmat_der

! *****************************************************************************
!> \brief Check Numerical Vs Analytical
!> \param sepi ...
!> \param sepj ...
!> \param rijv ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param invert ...
!> \param ii ...
!> \param kk ...
!> \param v_d ...
!> \par History
!>      04.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
!> \note
!>      Debug routine
! *****************************************************************************
SUBROUTINE rot_2el_2c_first_debug(sepi, sepj, rijv, se_int_control, se_taper, invert, ii, kk, v_d)
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_utils,        ONLY: rotmat
  USE semi_empirical_int_arrays,       ONLY: indexb
  USE semi_empirical_int_num,          ONLY: terep_num
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             rotmat_type,&
                                             rotmat_create,&
                                             rotmat_release,&
                                             se_int_control_type,&
                                             se_taper_type
  USE semi_empirical_int_utils,        ONLY: rot_2el_2c_first
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(semi_empirical_type), POINTER       :: sepi, sepj
  REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rijv
  TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
  LOGICAL, INTENT(IN)                      :: invert
  TYPE(se_taper_type), POINTER             :: se_taper
  INTEGER, INTENT(IN)                      :: ii, kk
  REAL(KIND=dp), DIMENSION(3, 45, 45), &
    INTENT(IN)                             :: v_d

  CHARACTER(len=*), PARAMETER :: moduleN='semi_empirical_int_debug',&
       routineN = 'rot_2el_2c_first', routineP = moduleN//':'//routineN

  REAL(KIND=dp), DIMENSION(491)            :: rep
  LOGICAL, DIMENSION(45, 45)               :: logv
  REAL(KIND=dp), DIMENSION(45, 45)         :: v_n, v_p, v_m

  LOGICAL                                  :: failure, check_value
  REAL(KIND=dp)                            :: dx, r, r0(3), x(3)
  TYPE(rotmat_type), POINTER               :: matrix
  INTEGER                                  :: imap(3),i,j,k, limkl

  NULLIFY(matrix)
  failure = .FALSE.
  dx = 1.0E-6_dp
  imap(1) = 1
  imap(2) = 2
  imap(3) = 3
  IF (invert) THEN
     imap(1) = 3
     imap(2) = 2
     imap(3) = 1
  END IF
  limkl = indexb(kk, kk)
  ! Check derivatives: analytical VS numerical
  WRITE(*,*)"DEBUG::"//routineP
  DO j = 1, 3
     x          = 0.0_dp
     x(imap(j)) = dx
     DO i = 1, 2
        r0 = rijv + (-1.0_dp)**(i-1) * x
        r  = SQRT(DOT_PRODUCT(r0,r0))

        CALL rotmat_create(matrix)
        CALL rotmat (sepi, sepj, r0, r, matrix, do_derivatives=.FALSE., debug_invert=invert)

        ! Compute integrals in diatomic frame
        CALL terep_num(sepi,sepj,r,rep,se_taper=se_taper,se_int_control=se_int_control)

        IF (i==1) THEN
           CALL rot_2el_2c_first(sepi, sepj, r0, se_int_control, se_taper, invert, ii, kk, rep, logv, matrix,&
                v_p, lgrad=.FALSE.)
        END IF
        IF (i==2) THEN
           CALL rot_2el_2c_first(sepi, sepj, r0, se_int_control, se_taper, invert, ii, kk, rep, logv, matrix,&
                v_m, lgrad=.FALSE.)
        END IF
        CALL rotmat_release(matrix)
     END DO
     ! Check numerical VS analytical
     DO i = 1, 45
        DO k = 1, limkl
           ! Compute the numerical derivative
           v_n(i,k) = (v_p(i,k)-v_m(i,k))/(2.0_dp*dx)
           IF (.NOT.check_value(v_d(j,i,k), v_n(i,k), dx, 0.1_dp)) THEN
              WRITE(*,*)"ERROR for  rot_2el_2c_first derivative V_D(j,i,k), j,i,k::",j,i,k
              CPABORT("")
           END IF
        END DO
     END DO
  END DO
END SUBROUTINE rot_2el_2c_first_debug

! *****************************************************************************
!> \brief Check Numerical Vs Analytical ssss
!> \param sepi ...
!> \param sepj ...
!> \param r ...
!> \param dssss ...
!> \param itype ...
!> \param se_int_control ...
!> \param se_taper ...
!> \par History
!>      04.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
!> \note
!>      Debug routine
! *****************************************************************************
SUBROUTINE check_dssss_nucint_ana (sepi,sepj,r,dssss,itype,se_int_control,se_taper)
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_num,          ONLY: ssss_nucint_num
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             se_int_control_type,&
                                             se_taper_type
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(semi_empirical_type), POINTER       :: sepi, sepj
  REAL(dp), INTENT(IN)                     :: r
  REAL(dp), INTENT(IN)                     :: dssss
  INTEGER, INTENT(IN)                      :: itype
  TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
  TYPE(se_taper_type), POINTER                :: se_taper

  CHARACTER(len=*), PARAMETER :: moduleN='semi_empirical_int_debug',&
       routineN = 'check_dssss_nucint_ana', routineP = moduleN//':'//routineN

  REAL(dp)                                 :: delta, nssss, od, rn, ssssm, ssssp
  LOGICAL                                  :: check_value

  delta = 1.0E-8_dp
  od = 0.5_dp/delta
  rn = r + delta
  CALL ssss_nucint_num(sepi,sepj,rn,ssssp,itype,se_taper,se_int_control)
  rn = r - delta
  CALL ssss_nucint_num(sepi,sepj,rn,ssssm,itype,se_taper,se_int_control)
  nssss = od * (ssssp - ssssm)
  ! check
  WRITE(*,*)"DEBUG::"//routineP
  IF (.NOT.check_value(nssss, dssss, delta, 0.1_dp)) THEN
     WRITE(*,*)"ERROR for SSSS derivative SSSS"
     CPABORT("")
  END IF
END SUBROUTINE check_dssss_nucint_ana

! *****************************************************************************
!> \brief Check Numerical Vs Analytical core
!> \param sepi ...
!> \param sepj ...
!> \param r ...
!> \param dcore ...
!> \param itype ...
!> \param se_int_control ...
!> \param se_taper ...
!> \par History
!>      04.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
!> \note
!>      Debug routine
! *****************************************************************************
SUBROUTINE check_dcore_nucint_ana (sepi,sepj,r,dcore,itype,se_int_control,se_taper)
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_num,          ONLY: core_nucint_num
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             se_int_control_type,&
                                             se_taper_type
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(semi_empirical_type), POINTER       :: sepi, sepj
  REAL(dp), INTENT(IN)                     :: r
  REAL(dp), DIMENSION(10, 2), INTENT(IN)   :: dcore
  INTEGER, INTENT(IN)                      :: itype
  TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
  TYPE(se_taper_type), POINTER             :: se_taper

  CHARACTER(len=*), PARAMETER :: moduleN='semi_empirical_int_debug',&
       routineN = 'check_dcore_nucint_ana', routineP = moduleN//':'//routineN

  INTEGER                                  :: i, j
  REAL(dp)                                 :: delta, od, rn
  REAL(dp), DIMENSION(10, 2)               :: corem, corep, ncore
  LOGICAL                                  :: check_value

  delta = 1.0E-8_dp
  od = 0.5_dp/delta
  rn = r + delta
  CALL core_nucint_num(sepi,sepj,rn,corep,itype,se_taper,se_int_control)
  rn = r - delta
  CALL core_nucint_num(sepi,sepj,rn,corem,itype,se_taper,se_int_control)
  ncore = od * (corep - corem)
  ! check
  WRITE(*,*)"DEBUG::"//routineP
  DO i = 1, 2
     DO j = 1, 10
        IF (.NOT.check_value(ncore(j,i), dcore(j,i), delta, 0.1_dp)) THEN
           WRITE(*,*)"ERROR for CORE derivative CORE(j,i), j,i::",j,i
           CPABORT("")
        END IF
     END DO
  END DO
END SUBROUTINE check_dcore_nucint_ana

! *****************************************************************************
!> \brief Low level comparison function between numberical and analytical values
!> \param num ...
!> \param ana ...
!> \param minval ...
!> \param thrs ...
!> \retval passed ...
!> \par History
!>      04.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
!> \note
!>      Debug routine
! *****************************************************************************
FUNCTION check_value(num, ana, minval, thrs) RESULT(passed)
  USE kinds,                           ONLY: dp
  IMPLICIT NONE
  REAL(KIND=dp)                            :: num, ana, thrs, minval
  LOGICAL                                  :: passed

  passed = .TRUE.

  IF      ((ABS(num)<minval).AND.(ABS(ana)>minval)) THEN
     WRITE(*,*)"WARNING ---> ",num, ana, thrs
  ELSE IF ((ABS(num)>minval).AND.(ABS(ana)<minval)) THEN
     WRITE(*,*)"WARNING ---> ",num, ana, thrs
  ELSE IF ((ABS(num)<minval).AND.(ABS(ana)<minval)) THEN
     ! skip..
     RETURN
  END IF
  IF (ABS((num-ana)/num*100._dp)>thrs) THEN
     WRITE(*,*) ABS(num-ana)/num*100._dp, thrs
     passed = .FALSE.
  END IF
  IF (.NOT.passed) THEN
     WRITE(*,*)"ANALYTICAL ::",ana
     WRITE(*,*)"NUMERICAL  ::",num
  END IF
END FUNCTION check_value

! *****************************************************************************
!> \brief Check Numerical Vs Analytical
!> \param sepi ...
!> \param sepj ...
!> \param rijv ...
!> \param itype ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param e1b ...
!> \param e2a ...
!> \param de1b ...
!> \param de2a ...
!> \par History
!>      04.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
!> \note
!>      Debug routine
! *****************************************************************************
SUBROUTINE check_drotnuc_ana(sepi, sepj, rijv, itype, se_int_control, se_taper, e1b, e2a, de1b, de2a)
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_num,          ONLY: rotnuc_num,&
                                             drotnuc_num
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             se_int_control_type,&
                                             se_taper_type
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(semi_empirical_type), POINTER       :: sepi, sepj
  REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
  INTEGER, INTENT(IN)                      :: itype
  TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
  TYPE(se_taper_type), POINTER             :: se_taper
  REAL(dp), DIMENSION(45), INTENT(IN), &
    OPTIONAL                               :: e1b, e2a
  REAL(dp), DIMENSION(3, 45), &
    INTENT(IN), OPTIONAL                   :: de1b, de2a

  CHARACTER(len=*), PARAMETER :: moduleN = 'semi_empirical_int_debug', &
    routineN = 'check_drotnuc_ana', routineP = moduleN//':'//routineN

  INTEGER                                  :: i, j
  LOGICAL                                  :: l_de1b, l_de2a, l_e1b, l_e2a,&
                                              lgrad, check_value
  REAL(dp)                                 :: delta
  REAL(KIND=dp), DIMENSION(45)             :: e1b2, e2a2
  REAL(KIND=dp), DIMENSION(3, 45)          :: de1b2, de2a2

  l_e1b  = PRESENT(e1b)
  l_e2a  = PRESENT(e2a)
  l_de1b = PRESENT(de1b)
  l_de2a = PRESENT(de2a)
  lgrad  = l_de1b.OR.l_de2a
  delta = 1.0E-5_dp
  ! Check value of integrals
  WRITE(*,*)"DEBUG::"//routineP
  CALL rotnuc_num(sepi,sepj,rijv,e1b2,e2a2,itype,se_int_control,se_taper=se_taper)
  IF (l_e1b) THEN
     DO j = 1, 45
        IF (.NOT.check_value(e1b2(j), e1b(j), delta, 0.1_dp)) THEN
           WRITE(*,*)"ERROR for E1B value E1B(j), j::",j
           CPABORT("")
        END IF
     END DO
  END IF
  IF (l_e2a) THEN
     DO J = 1, 45
        IF (.NOT.check_value(e2a2(j), e2a(j), delta, 0.1_dp)) THEN
           WRITE(*,*)"ERROR for E2A value E2A(j), j::",j
           CPABORT("")
        END IF
     END DO
  END IF

  ! Check derivatives
  IF (lgrad) THEN
     CALL drotnuc_num(sepi,sepj,rijv,de1b2,de2a2,itype,delta=delta,&
          se_int_control=se_int_control,se_taper=se_taper)
     IF (l_de1b) THEN
        DO i = 1, 3
           DO j = 1, 45
              ! Additional check on the value of the integral before checking derivatives
              IF (ABS(e1b2(j))>delta) THEN
                 IF (.NOT.check_value(de1b2(i,j), de1b(i,j), delta, 0.1_dp)) THEN
                    WRITE(*,*)"ERROR for derivative of E1B.  DE1B(i,j), i,j::",i,j
                    CPABORT("")
                 END IF
              END IF
           END DO
        END DO
     END IF
     IF (l_de2a) THEN
        DO i = 1, 3
           DO j = 1, 45
              ! Additional check on the value of the integral before checking derivatives
              IF (ABS(e2a2(j))>delta) THEN
                 IF (.NOT.check_value(de2a2(i,j), de2a(i,j), delta, 0.1_dp)) THEN
                    WRITE(*,*)"ERROR for derivative of E2A.  DE2A(i,j), i,j::",i,j
                    CPABORT("")
                 END IF
              END IF
           END DO
        END DO
     END IF
  END IF
END SUBROUTINE check_drotnuc_ana

! *****************************************************************************
!> \brief Check Numerical Vs Analytical CORE-CORE
!> \param sepi ...
!> \param sepj ...
!> \param rijv ...
!> \param itype ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param enuc ...
!> \param denuc ...
!> \par History
!>      04.2007 created [tlaino]
!> \author Teodoro Laino - Zurich University
!> \note
!>      Debug routine
! *****************************************************************************
SUBROUTINE check_dcorecore_ana(sepi, sepj, rijv, itype,se_int_control, se_taper, enuc, denuc)
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_num,          ONLY: corecore_num,&
                                             dcorecore_num
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             se_int_control_type,&
                                             se_taper_type
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(semi_empirical_type), POINTER       :: sepi, sepj
  REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
  INTEGER, INTENT(IN)                      :: itype
  REAL(dp), INTENT(IN), OPTIONAL           :: enuc
  REAL(dp), DIMENSION(3), INTENT(IN), &
       OPTIONAL                            :: denuc
  TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
  TYPE(se_taper_type), POINTER             :: se_taper

  CHARACTER(len=*), PARAMETER :: moduleN = 'semi_empirical_int_debug', &
       routineN = 'check_dcorecore_ana', routineP = moduleN//':'//routineN

  LOGICAL                                  :: check_value
  INTEGER                                  :: j
  REAL(dp)                                 :: enuc_num, delta
  REAL(dp), DIMENSION(3)                   :: denuc_num

  WRITE(*,*)"DEBUG::"//routineP
  delta = 1.0E-7_dp
  ! check
  IF (PRESENT(enuc)) THEN
     CALL  corecore_num(sepi,sepj,rijv, enuc_num,itype,se_int_control,se_taper)
     IF (.NOT.check_value(enuc, enuc_num, delta, 0.001_dp)) THEN
        WRITE(*,*)"ERROR for CORE-CORE energy value (numerical different from analytical)!!"
        CPABORT("")
     END IF
  END IF
  IF (PRESENT(denuc)) THEN
     CALL dcorecore_num(sepi,sepj,rijv,denuc_num,itype,delta,se_int_control,se_taper)
     DO j = 1, 3
        IF (.NOT.check_value(denuc(j), denuc_num(j), delta, 0.001_dp)) THEN
           WRITE(*,*)"ERROR for CORE-CORE energy derivative value (numerical different from analytical). DENUC(j), j::", j
           CPABORT("")
        END IF
     END DO
  END IF
END SUBROUTINE check_dcorecore_ana

! *****************************************************************************
!> \brief Check Numerical Vs Analytical DTEREP_ANA
!> \param sepi ...
!> \param sepj ...
!> \param r ...
!> \param ri ...
!> \param dri ...
!> \param se_int_control ...
!> \param se_taper ...
!> \param lgrad ...
!> \par History
!>      04.2007 created [tlaino]
!> \author Teodoro Laino - Zurich University
!> \note
!>      Debug routine
! *****************************************************************************
SUBROUTINE check_dterep_ana (sepi,sepj,r,ri,dri,se_int_control,se_taper,lgrad)
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_num,          ONLY: terep_num
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             se_int_control_type,&
                                             se_taper_type
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(semi_empirical_type), POINTER       :: sepi, sepj
  REAL(dp), INTENT(IN)                     :: r
  REAL(dp), DIMENSION(491), INTENT(IN)     :: ri, dri
  TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
  LOGICAL, INTENT(IN)                      :: lgrad
  TYPE(se_taper_type), POINTER             :: se_taper

  CHARACTER(len=*), PARAMETER :: moduleN = 'semi_empirical_int_debug', &
       routineN = 'check_dterep_ana', routineP = moduleN//':'//routineN

  LOGICAL                                  :: check_value
  INTEGER                                  :: j, i
  REAL(dp)                                 :: delta, od, rn
  REAL(dp), DIMENSION(491)                 :: nri, ri0, rim, rip

  delta = 1.0E-8_dp
  od = 0.5_dp/delta
  rn = r
  CALL terep_num(sepi,sepj,rn,ri0,se_taper,se_int_control)
  IF (lgrad) THEN
     rn = r + delta
     CALL terep_num(sepi,sepj,rn,rip,se_taper,se_int_control)
     rn = r - delta
     CALL terep_num(sepi,sepj,rn,rim,se_taper,se_int_control)
     nri = od * (rip - rim)
  END IF
  ! check
  WRITE(*,*)"DEBUG::"//routineP
  DO j = 1, 491
     IF  (ABS(ri(j)-ri0(j))>EPSILON(0.0_dp)) THEN
        WRITE(*,*)"Error in value of the integral RI",j,ri(j),ri0(j)
        CPABORT("")
     END IF
     IF (lgrad) THEN
        IF (.NOT.check_value(dri(j), nri(j), delta*100.0_dp, 0.1_dp)) THEN
           WRITE(*,*)"ERROR for derivative of RI integral, RI(j), j::",j
           WRITE(*,*)"FULL SET OF INTEGRALS: INDX  ANAL  NUM   DIFF"
           DO i = 1, 491
              WRITE(*,'(I5,3F15.9)')i, dri(i), nri(i), nri(i)-dri(i)
           END DO
           CPABORT("")
        END IF
     END IF
  END DO
END SUBROUTINE check_dterep_ana

! *****************************************************************************
!> \brief Check Numerical Vs Analytical ROTINT_ANA
!> \param sepi ...
!> \param sepj ...
!> \param rijv ...
!> \param w ...
!> \param dw ...
!> \param se_int_control ...
!> \param se_taper ...
!> \par History
!>      04.2008 created [tlaino]
!> \author Teodoro Laino - Zurich University
!> \note
!>      Debug routine
! *****************************************************************************
SUBROUTINE check_rotint_ana(sepi,sepj,rijv,w,dw,se_int_control,se_taper)
  
  USE kinds,                           ONLY: dp
  USE semi_empirical_int_num,          ONLY: rotint_num,&
                                             drotint_num
  USE semi_empirical_types,            ONLY: semi_empirical_type,&
                                             se_int_control_type,&
                                             se_taper_type
#include "./common/cp_common_uses.f90"
  IMPLICIT NONE
  TYPE(semi_empirical_type), POINTER       :: sepi, sepj
  REAL(dp), DIMENSION(3), INTENT(IN)       :: rijv
  REAL(dp), DIMENSION(2025), INTENT(IN), &
    OPTIONAL                               :: w
  REAL(dp), DIMENSION(3, 2025), &
    INTENT(IN), OPTIONAL                   :: dw
  TYPE(se_int_control_type), INTENT(IN)    :: se_int_control
  TYPE(se_taper_type), POINTER             :: se_taper

  CHARACTER(len=*), PARAMETER :: moduleN = 'semi_empirical_int_debug', &
       routineN = 'rotint_ana', routineP = moduleN//':'//routineN

  REAL(dp), DIMENSION(2025)                :: w2
  REAL(dp), DIMENSION(3, 2025)             :: dw2
  LOGICAL                                  :: check_value
  INTEGER                                  :: j, i
  REAL(KIND=dp)                            :: delta

  delta = 1.0E-6_dp
  WRITE(*,*)"DEBUG::"//routineP
  IF (PRESENT(w)) THEN
     w2=0.0_dp
     CALL rotint_num(sepi,sepj,rijv,w2,se_int_control,se_taper=se_taper)
     DO j = 1, 2025
        IF (.NOT.check_value(w(j), w2(j), delta, 0.1_dp)) THEN
           WRITE(*,*)"ERROR for integral value W(j), j::",j
           CPABORT("")
        END IF
     END DO
  END IF
  IF (PRESENT(dw)) THEN
     ! Numerical derivatives are obviosly a big problem..
     ! First of all let's decide if the value we get for delta is compatible
     ! with a reasonable value of the integral.. (compatible if the value of the
     ! integral is greater than 1.0E-6)
     CALL drotint_num(sepi,sepj,rijv,dw2,delta=delta,se_int_control=se_int_control,se_taper=se_taper)
     CALL rotint_num(sepi,sepj,rijv,w2,se_int_control=se_int_control,se_taper=se_taper)
     DO i = 1, 3
        DO j = 1, 2025
           IF ((ABS(w2(j))>delta).AND.(ABS(dw2(i,j))>delta*10)) THEN
              IF (.NOT.check_value(dw(i,j), dw2(i,j), delta, 0.1_dp)) THEN
                 WRITE(*,*)"ERROR for derivative of the integral value W(j). DW(i,j) i,j::",i,j
                 CPABORT("")
              END IF
           END IF
        END DO
     END DO
  END IF
END SUBROUTINE check_rotint_ana
