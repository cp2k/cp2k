!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate RI-RPA energy and Sigma correction to the RPA energies
!>         using the cubic spline based on eigen values of Q(w).
!> \par History
! **************************************************************************************************
MODULE rpa_sigma_functional
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE input_constants,                 ONLY: sigma_PBE0_S1,&
                                              sigma_PBE0_S2,&
                                              sigma_PBE_S1,&
                                              sigma_PBE_S2,&
                                              sigma_none
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_sigma_functional'

   PUBLIC :: rpa_sigma_matrix_spectral, rpa_sigma_create, rpa_sigma_type, finalize_rpa_sigma

   TYPE rpa_sigma_type
      PRIVATE
      REAL(KIND=dp)                                      :: e_sigma_corr = 0.0_dp
      REAL(KIND=dp)                                      :: e_rpa_by_eig_val = 0.0_dp
      INTEGER                                            :: sigma_param = 0
      TYPE(cp_fm_type)                                   :: mat_Q_diagonal = cp_fm_type()
      TYPE(cp_fm_type)                                   :: fm_evec = cp_fm_type()
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: sigma_eigenvalue
      INTEGER                                            :: dimen_RI_red = 0
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ... Collect the Q(w) (fm_mat_Q) matrix to create rpa_sigma a derived type variable.
!>             and write out the choosen parametrization for the cubic spline.
!> \param rpa_sigma ...
!> \param sigma_param ...
!> \param fm_mat_Q ...
!> \param unit_nr ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE rpa_sigma_create(rpa_sigma, sigma_param, fm_mat_Q, unit_nr, para_env)

      TYPE(rpa_sigma_type), INTENT(OUT)                  :: rpa_sigma
      INTEGER                                            :: sigma_param
      TYPE(cp_fm_type)                                   :: fm_mat_Q
      INTEGER                                            :: unit_nr

      CLASS(mp_comm_type), INTENT(IN)                    :: para_env

      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct

      ! Getting information about the Q matrix and creating initializing two matrices to pass it to the diagonalising driver.
      CALL cp_fm_get_info(fm_mat_Q, matrix_struct=matrix_struct, nrow_global=rpa_sigma%dimen_RI_red)

      ALLOCATE (rpa_sigma%sigma_eigenvalue(rpa_sigma%dimen_RI_red))

      CALL cp_fm_create(rpa_sigma%fm_evec, matrix_struct)
      CALL cp_fm_create(rpa_sigma%mat_Q_diagonal, matrix_struct)

      rpa_sigma%sigma_param = sigma_param

      SELECT CASE (rpa_sigma%sigma_param)

      CASE (sigma_none)
         ! There is nothing to do
      CASE DEFAULT
         CPABORT("Unknown parameterization")

      CASE (sigma_PBE0_S1)
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3, A)") &
            "SIGMA_INFO| Sigma eigenvalues parameterized with PBE0_S1 reference"

      CASE (sigma_PBE0_S2)
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3, A)") &
            "SIGMA_INFO| Sigma eigenvalues parameterized with PBE0_S2 reference"

      CASE (sigma_PBE_S1)
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3, A)") &
            "SIGMA_INFO| Sigma eigenvalues parameterized with PBE_S1 reference"

      CASE (sigma_PBE_S2)
         IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3, A)") &
            "SIGMA_INFO| Sigma eigenvalues parameterized with PBE_S2 reference"
      END SELECT
      IF (unit_nr > 0) CALL m_flush(unit_nr)
      CALL para_env%sync()

   END SUBROUTINE rpa_sigma_create

! **************************************************************************************************
!> \brief ... Memory cleanup routine
!> \param rpa_sigma ...
! **************************************************************************************************
   SUBROUTINE rpa_sigma_cleanup(rpa_sigma)

      TYPE(rpa_sigma_type), INTENT(INOUT)                :: rpa_sigma

      CALL cp_fm_release(rpa_sigma%mat_Q_diagonal)
      CALL cp_fm_release(rpa_sigma%fm_evec)
      DEALLOCATE (rpa_sigma%sigma_eigenvalue)

   END SUBROUTINE rpa_sigma_cleanup

! **************************************************************************************************
!> \brief ... Diagonalize and store the eigenvalues of fm_mat_Q in rpa_sigma%sigma_eigenvalue.
!> \param rpa_sigma ...
!> \param fm_mat_Q ...
!> \param wj ...
!> \param para_env_RPA ...
! **************************************************************************************************
   SUBROUTINE rpa_sigma_matrix_spectral(rpa_sigma, fm_mat_Q, wj, para_env_RPA)

      TYPE(rpa_sigma_type)                               :: rpa_sigma
      TYPE(cp_fm_type)                                   :: fm_mat_Q
      REAL(KIND=dp)                                      :: wj
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_RPA

      ! copy the Q matrix into the dummy matrix to avoid changing it.
      CALL cp_fm_to_fm(fm_mat_Q, rpa_sigma%mat_Q_diagonal)

      !diagonalising driver
      CALL choose_eigv_solver(rpa_sigma%mat_Q_diagonal, rpa_sigma%fm_evec, rpa_sigma%sigma_eigenvalue)

      ! Computing the integration to calculate the sigma correction.
      CALL compute_e_sigma_corr_by_freq_int(rpa_sigma, wj, para_env_RPA)

   END SUBROUTINE rpa_sigma_matrix_spectral
! **************************************************************************************************

! **************************************************************************************************
!> \brief ... To compute the e_sigma_corr and e_rpa_by_eig_val by freq integration over the eigenvalues of Q(w)
!>            e_sigma_corr = - H(sigma) &  e_rpa_by_eig_val = log(1+sigma)-sigma
!> \param rpa_sigma ...
!> \param wj ...
!> \param para_env_RPA ...
! **************************************************************************************************
   SUBROUTINE compute_e_sigma_corr_by_freq_int(rpa_sigma, wj, para_env_RPA)
      TYPE(rpa_sigma_type), INTENT(INOUT)                :: rpa_sigma
      REAL(KIND=dp), INTENT(IN)                          :: wj
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_RPA

      INTEGER                                            :: iaux
      REAL(KIND=dp)                                      :: dedw_rpa, dedw_sigma

      dedw_sigma = 0.0_dp
      dedw_rpa = 0.0_dp

      ! Loop which  take each eigenvalue to: i) get E_RPA & ii) integrates the spline to get E_c correction
      DO iaux = 1, rpa_sigma%dimen_RI_red
         IF (rpa_sigma%sigma_eigenvalue(iaux) > 0.0_dp) THEN
            IF (MODULO(iaux, para_env_RPA%num_pe) /= para_env_RPA%mepos) CYCLE
            dedw_rpa = dedw_rpa + LOG(1.0_dp + rpa_sigma%sigma_eigenvalue(iaux)) - rpa_sigma%sigma_eigenvalue(iaux)
            IF (MODULO(iaux, para_env_RPA%num_pe) /= para_env_RPA%mepos) CYCLE
            dedw_sigma = dedw_sigma - cubic_spline_integr(rpa_sigma%sigma_eigenvalue(iaux), rpa_sigma%sigma_param)
         END IF
      END DO

      ! (use 2.0_dp its better for compilers)
      rpa_sigma%e_sigma_corr = rpa_sigma%e_sigma_corr + (wj*dedw_sigma/(2.0_dp*pi*2.0_dp))
      rpa_sigma%e_rpa_by_eig_val = rpa_sigma%e_rpa_by_eig_val + (wj*dedw_rpa/(2.0_dp*pi*2.0_dp))

   END SUBROUTINE compute_e_sigma_corr_by_freq_int

! **************************************************************************************************
!> \brief ... Save the calculated value of E_c correction to the global variable and  memory clean.
!> \param rpa_sigma ...
!> \param unit_nr ...
!> \param e_sigma_corr ...
!> \param para_env ...
!> \param do_minimax_quad ...
! **************************************************************************************************
   SUBROUTINE finalize_rpa_sigma(rpa_sigma, unit_nr, e_sigma_corr, para_env, do_minimax_quad)
      TYPE(rpa_sigma_type), INTENT(INOUT)                :: rpa_sigma
      INTEGER                                            :: unit_nr
      REAL(KIND=dp), INTENT(OUT)                         :: e_sigma_corr
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      LOGICAL, INTENT(IN)                                :: do_minimax_quad

      IF (do_minimax_quad) rpa_sigma%e_rpa_by_eig_val = rpa_sigma%e_rpa_by_eig_val/2.0_dp
      CALL para_env%sum(rpa_sigma%e_rpa_by_eig_val)
      e_sigma_corr = rpa_sigma%e_sigma_corr
      IF (unit_nr > 0) WRITE (unit_nr, '(T3,A,T56,F25.14)') &
         'RI-RPA energy from eigenvalues of Q(w)  = ', &
         rpa_sigma%e_rpa_by_eig_val

      CALL rpa_sigma_cleanup(rpa_sigma)

   END SUBROUTINE finalize_rpa_sigma

! **************************************************************************************************
!> \brief ... integrates cubic spline to get eigenvalue sigma based E_c correction.
!> \param sigma ...
!> \param sigma_param ...
!> \return ... Integration value H(sigma) wrt to coupling constant eq 35 in Trushin et al  JCP 2021
! **************************************************************************************************
   FUNCTION cubic_spline_integr(sigma, sigma_param) RESULT(integral)
      REAL(KIND=dp), INTENT(in)                          :: sigma
      INTEGER                                            :: sigma_param
      REAL(KIND=dp)                                      :: integral

      INTEGER                                            :: i, m, n
      REAL(KIND=dp)                                      :: h
      REAL(KIND=dp), ALLOCATABLE                         :: coeff(:, :), x_coord(:)

      SELECT CASE (sigma_param)
      CASE (sigma_PBE0_S1)
         n = 21
         ALLOCATE (x_coord(n))
         ALLOCATE (coeff(4, n))

         coeff(1, 1) = 0.000000000000D+00
         coeff(1, 2) = 0.000000000000D+00
         coeff(1, 3) = -0.149500660756D-03
         coeff(1, 4) = -0.353017276233D-02
         coeff(1, 5) = -0.109810247734D-01
         coeff(1, 6) = -0.231246943777D-01
         coeff(1, 7) = -0.268999962858D-01
         coeff(1, 8) = -0.634751994007D-03
         coeff(1, 9) = 0.118792892470D-01
         coeff(1, 10) = -0.473431931326D-01
         coeff(1, 11) = -0.817589390539D-01
         coeff(1, 12) = 0.125726011069D-01
         coeff(1, 13) = 0.108028492092D+00
         coeff(1, 14) = 0.193548206759D+00
         coeff(1, 15) = 0.358395561305D-01
         coeff(1, 16) = -0.497714974829D-01
         coeff(1, 17) = 0.341059348835D-01
         coeff(1, 18) = 0.341050720155D-01
         coeff(1, 19) = 0.785549033229D-01
         coeff(1, 20) = 0.000000000000D+00
         coeff(1, 21) = 0.000000000000D+00
         coeff(2, 1) = 0.000000000000D+00
         coeff(2, 2) = 0.000000000000D+00
         coeff(2, 3) = -0.208376539581D+01
         coeff(2, 4) = -0.469755869285D+01
         coeff(2, 5) = -0.565503803415D+01
         coeff(2, 6) = -0.135502867642D+01
         coeff(2, 7) = 0.000000000000D+00
         coeff(2, 8) = 0.284340701746D+01
         coeff(2, 9) = 0.000000000000D+00
         coeff(2, 10) = -0.342695931351D+01
         coeff(2, 11) = 0.000000000000D+00
         coeff(2, 12) = 0.358739081268D+01
         coeff(2, 13) = 0.203368806130D+01
         coeff(2, 14) = 0.000000000000D+00
         coeff(2, 15) = -0.901387663218D+00
         coeff(2, 16) = 0.000000000000D+00
         coeff(2, 17) = 0.000000000000D+00
         coeff(2, 18) = 0.000000000000D+00
         coeff(2, 19) = 0.000000000000D+00
         coeff(2, 20) = 0.000000000000D+00
         coeff(2, 21) = 0.000000000000D+00
         coeff(3, 1) = -0.000000000000D+00
         coeff(3, 2) = -0.322176662524D+05
         coeff(3, 3) = -0.267090835643D+04
         coeff(3, 4) = -0.373532067350D+04
         coeff(3, 5) = -0.797121299000D+03
         coeff(3, 6) = 0.111299540119D+03
         coeff(3, 7) = 0.299284621116D+04
         coeff(3, 8) = -0.319333485618D+02
         coeff(3, 9) = -0.140910103454D+04
         coeff(3, 10) = -0.848330431187D+01
         coeff(3, 11) = 0.435025012278D+03
         coeff(3, 12) = -0.700327539634D+01
         coeff(3, 13) = 0.545486142353D+01
         coeff(3, 14) = -0.453346282407D+02
         coeff(3, 15) = 0.371921027910D+00
         coeff(3, 16) = 0.464101795796D+01
         coeff(3, 17) = -0.190069531714D-04
         coeff(3, 18) = 0.514345336660D-01
         coeff(3, 19) = -0.431543078188D-02
         coeff(3, 20) = -0.000000000000D+00
         coeff(3, 21) = 0.000000000000D+00
         coeff(4, 1) = 0.000000000000D+00
         coeff(4, 2) = 0.152897717268D+09
         coeff(4, 3) = 0.902815532735D+06
         coeff(4, 4) = 0.191760493084D+07
         coeff(4, 5) = 0.445372471512D+06
         coeff(4, 6) = 0.188362654331D+04
         coeff(4, 7) = -0.383203258784D+06
         coeff(4, 8) = -0.170027418959D+05
         coeff(4, 9) = 0.819629330224D+05
         coeff(4, 10) = 0.560228610945D+04
         coeff(4, 11) = -0.108203002413D+05
         coeff(4, 12) = -0.363378668069D+03
         coeff(4, 13) = -0.260332257619D+03
         coeff(4, 14) = 0.291068208088D+03
         coeff(4, 15) = 0.122322834276D+02
         coeff(4, 16) = -0.132875656470D+02
         coeff(4, 17) = 0.343356030115D-04
         coeff(4, 18) = -0.212958640167D-01
         coeff(4, 19) = 0.389311916174D-03
         coeff(4, 20) = 0.000000000000D+00
         coeff(4, 21) = 0.000000000000D+00
         x_coord(1) = 0.000000000000D+00
         x_coord(2) = 0.100000000000D-04
         x_coord(3) = 0.100000000000D-03
         x_coord(4) = 0.100000000000D-02
         x_coord(5) = 0.215443469000D-02
         x_coord(6) = 0.464158883000D-02
         x_coord(7) = 0.100000000000D-01
         x_coord(8) = 0.146779926762D-01
         x_coord(9) = 0.215443469003D-01
         x_coord(10) = 0.316227766017D-01
         x_coord(11) = 0.464158883361D-01
         x_coord(12) = 0.681292069058D-01
         x_coord(13) = 0.100000000000D+00
         x_coord(14) = 0.158489319246D+00
         x_coord(15) = 0.251188643151D+00
         x_coord(16) = 0.398107170553D+00
         x_coord(17) = 0.630957344480D+00
         x_coord(18) = 0.100000000000D+01
         x_coord(19) = 0.261015721568D+01
         x_coord(20) = 0.100000000000D+02
         x_coord(21) = 0.215443469000D+02

      CASE (sigma_PBE0_S2)
         n = 21
         ALLOCATE (x_coord(n))
         ALLOCATE (coeff(4, n))

         coeff(1, 1) = 0.000000000000D+00
         coeff(1, 2) = 0.000000000000D+00
         coeff(1, 3) = -0.431405252048D-04
         coeff(1, 4) = -0.182874853131D-02
         coeff(1, 5) = -0.852003132762D-02
         coeff(1, 6) = -0.218177403992D-01
         coeff(1, 7) = -0.305777654735D-01
         coeff(1, 8) = -0.870882903969D-02
         coeff(1, 9) = 0.137878988102D-01
         coeff(1, 10) = -0.284352007440D-01
         coeff(1, 11) = -0.798812002431D-01
         coeff(1, 12) = -0.334010771574D-02
         coeff(1, 13) = 0.934182748715D-01
         coeff(1, 14) = 0.204960802253D+00
         coeff(1, 15) = 0.213204380281D-01
         coeff(1, 16) = -0.401220283037D-01
         coeff(1, 17) = 0.321629738336D-01
         coeff(1, 18) = 0.321618301891D-01
         coeff(1, 19) = 0.808763912948D-01
         coeff(1, 20) = 0.000000000000D+00
         coeff(1, 21) = 0.000000000000D+00
         coeff(2, 1) = 0.000000000000D+00
         coeff(2, 2) = 0.000000000000D+00
         coeff(2, 3) = -0.661870777583D+00
         coeff(2, 4) = -0.289752912590D+01
         coeff(2, 5) = -0.558979946652D+01
         coeff(2, 6) = -0.267765704540D+01
         coeff(2, 7) = 0.000000000000D+00
         coeff(2, 8) = 0.389592612611D+01
         coeff(2, 9) = 0.000000000000D+00
         coeff(2, 10) = -0.382296397421D+01
         coeff(2, 11) = 0.000000000000D+00
         coeff(2, 12) = 0.327772498106D+01
         coeff(2, 13) = 0.239633724310D+01
         coeff(2, 14) = 0.000000000000D+00
         coeff(2, 15) = -0.726304793204D+00
         coeff(2, 16) = 0.000000000000D+00
         coeff(2, 17) = 0.000000000000D+00
         coeff(2, 18) = 0.000000000000D+00
         coeff(2, 19) = 0.000000000000D+00
         coeff(2, 20) = 0.000000000000D+00
         coeff(2, 21) = 0.000000000000D+00
         coeff(3, 1) = -0.000000000000D+00
         coeff(3, 2) = -0.862385254713D+04
         coeff(3, 3) = -0.192306222883D+04
         coeff(3, 4) = -0.520047462362D+04
         coeff(3, 5) = -0.877473657666D+03
         coeff(3, 6) = 0.841408344046D+02
         coeff(3, 7) = 0.216516760964D+04
         coeff(3, 8) = 0.296702212913D+03
         coeff(3, 9) = -0.867733655494D+03
         coeff(3, 10) = -0.188410055380D+03
         coeff(3, 11) = 0.336084151111D+03
         coeff(3, 12) = 0.489746728744D+01
         coeff(3, 13) = 0.158746877181D+02
         coeff(3, 14) = -0.562764882273D+02
         coeff(3, 15) = 0.134759277149D+01
         coeff(3, 16) = 0.399959778866D+01
         coeff(3, 17) = -0.251917983154D-04
         coeff(3, 18) = 0.563694092760D-01
         coeff(3, 19) = -0.444296223097D-02
         coeff(3, 20) = -0.000000000000D+00
         coeff(3, 21) = 0.000000000000D+00
         coeff(4, 1) = 0.000000000000D+00
         coeff(4, 2) = 0.366429086790D+08
         coeff(4, 3) = 0.504466528222D+06
         coeff(4, 4) = 0.232980923705D+07
         coeff(4, 5) = 0.392124287301D+06
         coeff(4, 6) = 0.206173887726D+05
         coeff(4, 7) = -0.249217659838D+06
         coeff(4, 8) = -0.563519876566D+05
         coeff(4, 9) = 0.448530826095D+05
         coeff(4, 10) = 0.143140667434D+05
         coeff(4, 11) = -0.800144415404D+04
         coeff(4, 12) = -0.391685311241D+03
         coeff(4, 13) = -0.414433988077D+03
         coeff(4, 14) = 0.376550449117D+03
         coeff(4, 15) = 0.510124747789D+01
         coeff(4, 16) = -0.114511339236D+02
         coeff(4, 17) = 0.455083767664D-04
         coeff(4, 18) = -0.233390912502D-01
         coeff(4, 19) = 0.400817027790D-03
         coeff(4, 20) = 0.000000000000D+00
         coeff(4, 21) = 0.000000000000D+00
         x_coord(1) = 0.000000000000D+00
         x_coord(2) = 0.100000000000D-04
         x_coord(3) = 0.100000000000D-03
         x_coord(4) = 0.100000000000D-02
         x_coord(5) = 0.215443469000D-02
         x_coord(6) = 0.464158883000D-02
         x_coord(7) = 0.100000000000D-01
         x_coord(8) = 0.146779926762D-01
         x_coord(9) = 0.215443469003D-01
         x_coord(10) = 0.316227766017D-01
         x_coord(11) = 0.464158883361D-01
         x_coord(12) = 0.681292069058D-01
         x_coord(13) = 0.100000000000D+00
         x_coord(14) = 0.158489319246D+00
         x_coord(15) = 0.251188643151D+00
         x_coord(16) = 0.398107170553D+00
         x_coord(17) = 0.630957344480D+00
         x_coord(18) = 0.100000000000D+01
         x_coord(19) = 0.261015721568D+01
         x_coord(20) = 0.100000000000D+02
         x_coord(21) = 0.215443469000D+02

      CASE (sigma_PBE_S1)
         n = 22
         ALLOCATE (x_coord(n))
         ALLOCATE (coeff(4, n))

         coeff(1, 1) = 0.000000000000D+00
         coeff(1, 2) = 0.000000000000D+00
         coeff(1, 3) = -0.493740326815D-04
         coeff(1, 4) = -0.136110637329D-02
         coeff(1, 5) = -0.506905111755D-02
         coeff(1, 6) = -0.127411222930D-01
         coeff(1, 7) = -0.220144968504D-01
         coeff(1, 8) = -0.239939034695D-01
         coeff(1, 9) = -0.436386416290D-01
         coeff(1, 10) = -0.117890214262D+00
         coeff(1, 11) = -0.141123921668D+00
         coeff(1, 12) = 0.865524876740D-01
         coeff(1, 13) = 0.179390274565D+00
         coeff(1, 14) = 0.269368658116D+00
         coeff(1, 15) = 0.785040456996D-01
         coeff(1, 16) = 0.490248637276D-01
         coeff(1, 17) = -0.111571911794D+00
         coeff(1, 18) = -0.197712184164D-01
         coeff(1, 19) = -0.197716870218D-01
         coeff(1, 20) = -0.372253617253D-01
         coeff(1, 21) = 0.000000000000D+00
         coeff(1, 22) = 0.000000000000D+00
         coeff(2, 1) = 0.000000000000D+00
         coeff(2, 2) = 0.000000000000D+00
         coeff(2, 3) = -0.709484897949D+00
         coeff(2, 4) = -0.197447407686D+01
         coeff(2, 5) = -0.315478745349D+01
         coeff(2, 6) = -0.229603163128D+01
         coeff(2, 7) = -0.670801534786D+00
         coeff(2, 8) = -0.704199644986D+00
         coeff(2, 9) = -0.400987325224D+01
         coeff(2, 10) = -0.269982990241D+01
         coeff(2, 11) = 0.000000000000D+00
         coeff(2, 12) = 0.472814414167D+01
         coeff(2, 13) = 0.207638470052D+01
         coeff(2, 14) = 0.000000000000D+00
         coeff(2, 15) = -0.389846972557D+00
         coeff(2, 16) = -0.298496119087D+00
         coeff(2, 17) = 0.000000000000D+00
         coeff(2, 18) = 0.000000000000D+00
         coeff(2, 19) = -0.601781536636D-06
         coeff(2, 20) = 0.000000000000D+00
         coeff(2, 21) = 0.000000000000D+00
         coeff(2, 22) = 0.000000000000D+00
         coeff(3, 1) = -0.000000000000D+00
         coeff(3, 2) = -0.104035132381D+05
         coeff(3, 3) = -0.108777473624D+04
         coeff(3, 4) = -0.219328637518D+04
         coeff(3, 5) = -0.260711341283D+03
         coeff(3, 6) = 0.132509852177D+02
         coeff(3, 7) = 0.165970301474D+03
         coeff(3, 8) = -0.460909893146D+03
         coeff(3, 9) = -0.112939707971D+04
         coeff(3, 10) = 0.465035067500D+02
         coeff(3, 11) = 0.123097490767D+04
         coeff(3, 12) = -0.876616265219D+02
         coeff(3, 13) = 0.790484996078D+01
         coeff(3, 14) = -0.624281400584D+02
         coeff(3, 15) = 0.324152775194D+01
         coeff(3, 16) = -0.632212496608D+01
         coeff(3, 17) = 0.202215332970D+01
         coeff(3, 18) = -0.308693235932D-06
         coeff(3, 19) = -0.495067060383D-02
         coeff(3, 20) = 0.116855980641D-02
         coeff(3, 21) = -0.000000000000D+00
         coeff(3, 22) = 0.000000000000D+00
         coeff(4, 1) = 0.000000000000D+00
         coeff(4, 2) = 0.478661516427D+08
         coeff(4, 3) = 0.285187385316D+06
         coeff(4, 4) = 0.971371823345D+06
         coeff(4, 5) = 0.116156741398D+06
         coeff(4, 6) = 0.172191903906D+05
         coeff(4, 7) = -0.241613612898D+05
         coeff(4, 8) = 0.213790845631D+05
         coeff(4, 9) = 0.790063233314D+05
         coeff(4, 10) = 0.201667888760D+04
         coeff(4, 11) = -0.344519214370D+05
         coeff(4, 12) = 0.963471669433D+03
         coeff(4, 13) = -0.292417702205D+03
         coeff(4, 14) = 0.433842720035D+03
         coeff(4, 15) = -0.132982468090D+02
         coeff(4, 16) = 0.199358142858D+02
         coeff(4, 17) = -0.365297127483D+01
         coeff(4, 18) = 0.434041376596D-07
         coeff(4, 19) = 0.101490424907D-02
         coeff(4, 20) = -0.796902275213D-04
         coeff(4, 21) = 0.000000000000D+00
         coeff(4, 22) = 0.000000000000D+00
         x_coord(1) = 0.000000000000D+00
         x_coord(2) = 0.100000000000D-04
         x_coord(3) = 0.100000000000D-03
         x_coord(4) = 0.100000000000D-02
         x_coord(5) = 0.215443469000D-02
         x_coord(6) = 0.464158883000D-02
         x_coord(7) = 0.100000000000D-01
         x_coord(8) = 0.146779926762D-01
         x_coord(9) = 0.215443469003D-01
         x_coord(10) = 0.316227766017D-01
         x_coord(11) = 0.464158883361D-01
         x_coord(12) = 0.681292069058D-01
         x_coord(13) = 0.100000000000D+00
         x_coord(14) = 0.158489319246D+00
         x_coord(15) = 0.251188643151D+00
         x_coord(16) = 0.398107170553D+00
         x_coord(17) = 0.630957344480D+00
         x_coord(18) = 0.100000000000D+01
         x_coord(19) = 0.237137370566D+01
         x_coord(20) = 0.562341325000D+01
         x_coord(21) = 0.153992652606D+02
         x_coord(22) = 0.316227766000D+02

      CASE (sigma_PBE_S2)
         n = 22
         ALLOCATE (x_coord(n))
         ALLOCATE (coeff(4, n))

         coeff(1, 1) = 0.000000000000D+00
         coeff(1, 2) = 0.000000000000D+00
         coeff(1, 3) = -0.156157535801D-03
         coeff(1, 4) = -0.365199003270D-02
         coeff(1, 5) = -0.108302033233D-01
         coeff(1, 6) = -0.203436953346D-01
         coeff(1, 7) = -0.214330355346D-01
         coeff(1, 8) = 0.109617244934D-03
         coeff(1, 9) = 0.813969827075D-02
         coeff(1, 10) = -0.701367130014D-01
         coeff(1, 11) = -0.162002361715D+00
         coeff(1, 12) = 0.337288711362D-01
         coeff(1, 13) = 0.140348429629D+00
         coeff(1, 14) = 0.271234417677D+00
         coeff(1, 15) = 0.780732751240D-01
         coeff(1, 16) = 0.436066976238D-01
         coeff(1, 17) = -0.106097689688D+00
         coeff(1, 18) = -0.133141637069D-01
         coeff(1, 19) = -0.133143525246D-01
         coeff(1, 20) = -0.430994711278D-01
         coeff(1, 21) = 0.000000000000D+00
         coeff(1, 22) = 0.000000000000D+00
         coeff(2, 1) = 0.000000000000D+00
         coeff(2, 2) = 0.000000000000D+00
         coeff(2, 3) = -0.217211651544D+01
         coeff(2, 4) = -0.473638379726D+01
         coeff(2, 5) = -0.487821808504D+01
         coeff(2, 6) = -0.433631413905D+00
         coeff(2, 7) = 0.000000000000D+00
         coeff(2, 8) = 0.193813387881D+01
         coeff(2, 9) = 0.000000000000D+00
         coeff(2, 10) = -0.695060290528D+01
         coeff(2, 11) = 0.000000000000D+00
         coeff(2, 12) = 0.502541925806D+01
         coeff(2, 13) = 0.273498669354D+01
         coeff(2, 14) = 0.000000000000D+00
         coeff(2, 15) = -0.448708826169D+00
         coeff(2, 16) = -0.332102918195D+00
         coeff(2, 17) = 0.000000000000D+00
         coeff(2, 18) = 0.000000000000D+00
         coeff(2, 19) = -0.242488141082D-06
         coeff(2, 20) = 0.000000000000D+00
         coeff(2, 21) = 0.000000000000D+00
         coeff(2, 22) = 0.000000000000D+00
         coeff(3, 1) = -0.000000000000D+00
         coeff(3, 2) = -0.337014964214D+05
         coeff(3, 3) = -0.285795351280D+04
         coeff(3, 4) = -0.372723918347D+04
         coeff(3, 5) = -0.516689374427D+03
         coeff(3, 6) = 0.480322803175D+02
         coeff(3, 7) = 0.253894893657D+04
         coeff(3, 8) = -0.535684993409D+02
         coeff(3, 9) = -0.162223464755D+04
         coeff(3, 10) = -0.319667723139D+03
         coeff(3, 11) = 0.101401359817D+04
         coeff(3, 12) = -0.862770702569D+02
         coeff(3, 13) = 0.212578002151D+02
         coeff(3, 14) = -0.625949163782D+02
         coeff(3, 15) = 0.357838438707D+01
         coeff(3, 16) = -0.543078279308D+01
         coeff(3, 17) = 0.204380282001D+01
         coeff(3, 18) = -0.124376927880D-06
         coeff(3, 19) = -0.844892173480D-02
         coeff(3, 20) = 0.135295689023D-02
         coeff(3, 21) = -0.000000000000D+00
         coeff(3, 22) = 0.000000000000D+00
         coeff(4, 1) = 0.000000000000D+00
         coeff(4, 2) = 0.160253203309D+09
         coeff(4, 3) = 0.106174857663D+07
         coeff(4, 4) = 0.211694319531D+07
         coeff(4, 5) = 0.377995031873D+06
         coeff(4, 6) = -0.941771032564D+03
         coeff(4, 7) = -0.332306990184D+06
         coeff(4, 8) = -0.850176312292D+04
         coeff(4, 9) = 0.844978829514D+05
         coeff(4, 10) = 0.249933770676D+05
         coeff(4, 11) = -0.275803550484D+05
         coeff(4, 12) = 0.105308484492D+04
         coeff(4, 13) = -0.508788318128D+03
         coeff(4, 14) = 0.432758845019D+03
         coeff(4, 15) = -0.144367801061D+02
         coeff(4, 16) = 0.175904487062D+02
         coeff(4, 17) = -0.369208055752D+01
         coeff(4, 18) = 0.174842962107D-07
         coeff(4, 19) = 0.173203285755D-02
         coeff(4, 20) = -0.922652326542D-04
         coeff(4, 21) = 0.000000000000D+00
         coeff(4, 22) = 0.000000000000D+00
         x_coord(1) = 0.000000000000D+00
         x_coord(2) = 0.100000000000D-04
         x_coord(3) = 0.100000000000D-03
         x_coord(4) = 0.100000000000D-02
         x_coord(5) = 0.215443469000D-02
         x_coord(6) = 0.464158883000D-02
         x_coord(7) = 0.100000000000D-01
         x_coord(8) = 0.146779926762D-01
         x_coord(9) = 0.215443469003D-01
         x_coord(10) = 0.316227766017D-01
         x_coord(11) = 0.464158883361D-01
         x_coord(12) = 0.681292069058D-01
         x_coord(13) = 0.100000000000D+00
         x_coord(14) = 0.158489319246D+00
         x_coord(15) = 0.251188643151D+00
         x_coord(16) = 0.398107170553D+00
         x_coord(17) = 0.630957344480D+00
         x_coord(18) = 0.100000000000D+01
         x_coord(19) = 0.237137370566D+01
         x_coord(20) = 0.562341325000D+01
         x_coord(21) = 0.153992652606D+02
         x_coord(22) = 0.316227766000D+02

      END SELECT

      ! determine to which interval sigma eigenvalue belongs
      m = intervalnum(x_coord, n, sigma)

      ! Numerically evaluate integral
      integral = 0.0_dp
      IF (m == 1) THEN
         integral = 0.5_dp*coeff(2, 1)*sigma
      END IF
      IF ((m > 1) .AND. (m < n)) THEN
         h = sigma - x_coord(m)
         integral = 0.5_dp*coeff(2, 1)*x_coord(2)**2/sigma &
                    + (coeff(1, m)*h + coeff(2, m)/2.0_dp*h**2 + &
                       coeff(3, m)/3.0_dp*h**3 + coeff(4, m)/4.0_dp*h**4)/sigma
         DO i = 2, m - 1
            h = x_coord(i + 1) - x_coord(i)
            integral = integral &
                       + (coeff(1, i)*h + coeff(2, i)/2.0_dp*h**2 + &
                          coeff(3, i)/3.0_dp*h**3 + coeff(4, i)/4.0_dp*h**4)/sigma
         END DO
      END IF
      IF (m == n) THEN
         integral = 0.5_dp*coeff(2, 1)*x_coord(2)**2/sigma
         DO i = 2, m - 1
            h = x_coord(i + 1) - x_coord(i)
            integral = integral &
                       + (coeff(1, i)*h + coeff(2, i)/2.0_dp*h**2 &
                          + coeff(3, i)/3.0_dp*h**3 + coeff(4, i)/4.0_dp*h**4)/sigma
         END DO
         integral = integral + coeff(1, n)*(1.0_dp - x_coord(n)/sigma)
      END IF
      integral = integral*sigma

      DEALLOCATE (x_coord, coeff)

   END FUNCTION cubic_spline_integr

! **************************************************************************************************
!> \brief ... Determine the interval which contains the sigma eigenvalue.
!> \param x_coord ...
!> \param n ...
!> \param sigma ...
!> \return ... an integer m
! **************************************************************************************************
   INTEGER FUNCTION intervalnum(x_coord, n, sigma) RESULT(inum)

      INTEGER                                            :: n
      REAL(KIND=dp), INTENT(in)                          :: x_coord(n), sigma

      INTEGER                                            :: i

      IF (sigma <= 0.0_dp) CPABORT('intervalnum: sigma should be positive')

      inum = -1
      IF (sigma > x_coord(n)) inum = n
      DO i = 1, n - 1
         IF ((sigma > x_coord(i)) .AND. (sigma <= x_coord(i + 1))) inum = i
      END DO

      IF (inum == -1) CPABORT('interval: something was wrong')

   END FUNCTION intervalnum

END MODULE rpa_sigma_functional
