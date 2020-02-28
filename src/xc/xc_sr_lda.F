!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates the short range correlation LDA energy (Improved version of Paola Gori-Giorgi's Code
!> \par History
!>      18-MAR-2002, TCH, working version
!>      fawzi (04.2004)  : adapted to the new xc interface
!> \see functionals_utilities
! **************************************************************************************************
MODULE xc_sr_lda
#:include "xc_perdew_wang.fypp"

   USE bibliography, ONLY: Paziani2006, &
                           Toulouse2004, &
                           cite_reference
   USE input_section_types, ONLY: section_vals_type, &
                                  section_vals_val_get
   USE kinds, ONLY: dp
   USE mathconstants, ONLY: pi
   USE xc_input_constants, ONLY: pw_dmc, &
                                 pw_orig, &
                                 pw_vmc
   USE xc_derivative_set_types, ONLY: xc_derivative_set_type, &
                                      xc_dset_get_derivative
   USE xc_derivative_types, ONLY: xc_derivative_get, &
                                  xc_derivative_type
   USE xc_rho_cflags_types, ONLY: xc_rho_cflags_type
   USE xc_rho_set_types, ONLY: xc_rho_set_get, &
                               xc_rho_set_type
   USE xc_functionals_utilities, ONLY: set_util
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_sr_lda'

@:global_var_pw92()
   REAL(KIND=dp), PARAMETER, PRIVATE :: &
      Acoul = 2._dp*(LOG(2._dp) - 1._dp)/pi**2, &
      aQ2 = 5.84605_dp, &
      cQ2 = 3.91744_dp, &
      dQ2 = 3.44851_dp, &
      bQ2 = dQ2 - 3._dp/(2._dp*pi*Acoul)*(4._dp/(9._dp*pi))**(1._dp/3._dp), &
      f02 = 4._dp/(9._dp*(2._dp**(1._dp/3._dp) - 1._dp)), &
      alpha = (4._dp/9._dp/pi)**(1._dp/3._dp), &
      cf = (9._dp*pi/4._dp)**(1._dp/3._dp), &
      p2p = 0.04_dp, &
      p3p = 0.4319_dp, &
      Cg0 = 0.0819306_dp, &
      Fg0 = 0.752411_dp, &
      Dg0 = -0.0127713_dp, &
      Eg0 = 0.00185898_dp, &
      Bg0 = 0.7317_dp - Fg0, &
      adib = 0.784949_dp, &
      q1a = -0.388_dp, &
      q2a = 0.676_dp, &
      q3a = 0.547_dp, &
      t1a = -4.95_dp, &
      t2a = 1._dp, &
      t3a = 0.31_dp

   PUBLIC :: sr_lda_info, sr_lda_eval, sr_lsd_eval

CONTAINS

! **************************************************************************************************
!> \brief Return some info on the functionals.
!> \param reference full reference
!> \param shortform short reference
!> \param lsd ...
!> \param needs ...
!> \param max_deriv ...
! **************************************************************************************************
   SUBROUTINE sr_lda_info(reference, shortform, lsd, needs, max_deriv)
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      LOGICAL, INTENT(IN), OPTIONAL                      :: lsd
      TYPE(xc_rho_cflags_type), INTENT(INOUT), OPTIONAL  :: needs
      INTEGER, INTENT(OUT), OPTIONAL                     :: max_deriv

      CHARACTER(len=*), PARAMETER :: routineN = 'sr_lda_info', routineP = moduleN//':'//routineN

      CALL cite_reference(Toulouse2004)
      CALL cite_reference(Paziani2006)
      IF (PRESENT(reference)) THEN
         reference = "J. Toulouse, A. Savin, and H.-J. Flad," &
                     //" Int. J. Quantum Chem. 100, 1074-1056 (2004)"
      END IF
      IF (PRESENT(shortform)) THEN
         shortform = "J. Toulouse et al., IJQC 100, 1074-1056 (2004)"
      END IF
      IF (PRESENT(needs)) THEN
         IF (lsd) THEN
            needs%rho_spin = .TRUE.
         ELSE
            needs%rho = .TRUE.
         END IF
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 1

   END SUBROUTINE sr_lda_info

@:init_pw92()

! **************************************************************************************************
!> \brief Calculate the correlation energy and its derivatives
!>      wrt to rho (the electron density) up to 3rd order. This
!>      is the short-range LDA version of the Perdew-Wang correlation energy
!>      If no order argument is given, then the routine calculates
!>      just the energy.
!> \param rho_set ...
!> \param deriv_set ...
!> \param order order of derivatives to calculate
!>      order must lie between -2 and 2. If it is negative then only
!>      that order will be calculated, otherwise all derivatives up to
!>      that order will be calculated.
!> \param sr_section ...
! **************************************************************************************************
   SUBROUTINE sr_lda_eval(rho_set, deriv_set, order, sr_section)

      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      INTEGER, INTENT(in)                                :: order
      TYPE(section_vals_type), POINTER                   :: sr_section

      CHARACTER(len=*), PARAMETER :: routineN = 'sr_lda_eval', routineP = moduleN//':'//routineN

      INTEGER                                            :: npoints, handle, method
      INTEGER, DIMENSION(:, :), POINTER                  :: bo
      REAL(KIND=dp)                                      :: omega, rho_cutoff, sc, sx
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dummy, e_0, e_rho, rho
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(sr_section, 'SCALE_X', r_val=sx)
      CALL section_vals_val_get(sr_section, 'SCALE_C', r_val=sc)
      CALL section_vals_val_get(sr_section, 'OMEGA', r_val=omega)
      CALL section_vals_val_get(sr_section, 'PARAMETRIZATION', i_val=method)

      NULLIFY (bo, rho, e_0, e_rho, dummy)
      CPASSERT(ASSOCIATED(rho_set))
      CPASSERT(rho_set%ref_count > 0)
      CPASSERT(ASSOCIATED(deriv_set))
      CPASSERT(deriv_set%ref_count > 0)
      CALL xc_rho_set_get(rho_set, rho=rho, &
                          local_bounds=bo, rho_cutoff=rho_cutoff)

      CALL perdew_wang_init(method, rho_cutoff)

      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      dummy => rho

      e_0 => dummy
      e_rho => dummy

      IF (order >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, "", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)
      END IF
      IF (order >= 1 .OR. order == -1) THEN
         deriv => xc_dset_get_derivative(deriv_set, "(rho)", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_rho)
      END IF
      IF (order > 1 .OR. order < -1) THEN
         CPABORT("derivatives bigger than 1 not implemented")
      END IF

      CALL sr_lda_calc(rho, omega, rho_cutoff, e_0, e_rho, &
                       npoints, order, sx, sc)

      CALL timestop(handle)

   END SUBROUTINE sr_lda_eval

! **************************************************************************************************
!> \brief ...
!> \param rho ...
!> \param omega ...
!> \param rho_cutoff ...
!> \param e_0 ...
!> \param e_rho ...
!> \param npoints ...
!> \param order ...
!> \param sx ...
!> \param sc ...
! **************************************************************************************************
   SUBROUTINE sr_lda_calc(rho, omega, rho_cutoff, e_0, e_rho, npoints, order, sx, sc)
      !FM low level calc routine
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rho
      REAL(KIND=dp), INTENT(IN)                          :: omega, rho_cutoff
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_0, e_rho
      INTEGER, INTENT(IN)                                :: npoints, order
      REAL(KIND=dp), INTENT(IN)                          :: sx, sc

      CHARACTER(len=*), PARAMETER :: routineN = 'sr_lda_calc', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, k
      REAL(KIND=dp)                                      :: my_rho, rs
      REAL(KIND=dp), DIMENSION(0:1)                      :: ed

      CALL timeset(routineN, handle)

      IF (sc /= 0.0_dp .OR. sx /= 0.0_dp) THEN
!$OMP PARALLEL DO PRIVATE (k, ed, my_rho, rs) DEFAULT(NONE)&
!$OMP SHARED(npoints,rho,rho_cutoff,omega,e_0,e_rho,order,sc,sx)
         DO k = 1, npoints

            my_rho = rho(k)
            IF (rho(k) > rho_cutoff) THEN
               rs = (3.0_dp/4.0_dp/pi/my_rho)**(1.0_dp/3.0_dp)

               CALL ldasr(rs, omega, ed(0), ed(1), sx, sc)

               IF (order >= 0) THEN
                  e_0(k) = e_0(k) + ed(0)*my_rho
               END IF
               IF (order >= 1 .OR. order == -1) THEN
                  e_rho(k) = e_rho(k) + ed(1)
               END IF
            END IF

         END DO
!$OMP END PARALLEL DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE sr_lda_calc

! **************************************************************************************************
!> \brief Calculate the correlation energy and its derivatives
!>      wrt to rho (the electron density) up to 3rd order. This
!>      is the short-range LSD version of the Perdew-Wang correlation energy
!>      If no order argument is given, then the routine calculates
!>      just the energy.
!> \param rho_set ...
!> \param deriv_set ...
!> \param order order of derivatives to calculate
!>      order must lie between -3 and 3. If it is negative then only
!>      that order will be calculated, otherwise all derivatives up to
!>      that order will be calculated.
!> \param sr_section ...
! **************************************************************************************************
   SUBROUTINE sr_lsd_eval(rho_set, deriv_set, order, sr_section)

      TYPE(xc_rho_set_type), POINTER                     :: rho_set
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      INTEGER, INTENT(IN), OPTIONAL                      :: order
      TYPE(section_vals_type), POINTER                   :: sr_section

      CHARACTER(len=*), PARAMETER :: routineN = 'sr_lsd_eval', routineP = moduleN//':'//routineN

      INTEGER                                            :: npoints, handle, method
      INTEGER, DIMENSION(:, :), POINTER                  :: bo
      REAL(KIND=dp)                                      :: omega, rho_cutoff, sc, sx
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: a, b, dummy, e_0, ea, eb
      TYPE(xc_derivative_type), POINTER                  :: deriv

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(sr_section, 'SCALE_X', r_val=sx)
      CALL section_vals_val_get(sr_section, 'SCALE_C', r_val=sc)
      CALL section_vals_val_get(sr_section, 'OMEGA', r_val=omega)
      CALL section_vals_val_get(sr_section, 'PARAMETRIZATION', i_val=method)

      NULLIFY (bo, a, b, e_0, ea, eb)
      CPASSERT(ASSOCIATED(rho_set))
      CPASSERT(rho_set%ref_count > 0)
      CPASSERT(ASSOCIATED(deriv_set))
      CPASSERT(deriv_set%ref_count > 0)
      CALL xc_rho_set_get(rho_set, rhoa=a, rhob=b, &
                          local_bounds=bo, rho_cutoff=rho_cutoff)

      CALL perdew_wang_init(method, rho_cutoff)

      npoints = (bo(2, 1) - bo(1, 1) + 1)*(bo(2, 2) - bo(1, 2) + 1)*(bo(2, 3) - bo(1, 3) + 1)

      ! meaningful default for the arrays we don't need: let us make compiler
      ! and debugger happy...
      dummy => a

      e_0 => dummy
      ea => dummy; eb => dummy

      IF (order >= 0) THEN
         deriv => xc_dset_get_derivative(deriv_set, "", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=e_0)
      END IF
      IF (order >= 1 .OR. order == -1) THEN
         deriv => xc_dset_get_derivative(deriv_set, "(rhoa)", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=ea)
         deriv => xc_dset_get_derivative(deriv_set, "(rhob)", &
                                         allocate_deriv=.TRUE.)
         CALL xc_derivative_get(deriv, deriv_data=eb)
      END IF
      IF (order > 1 .OR. order < -1) THEN
         CPABORT("derivatives bigger than 1 not implemented")
      END IF

      CALL sr_lsd_calc(a, b, omega, rho_cutoff, e_0, ea, eb, npoints, order, sx, sc)

      CALL timestop(handle)

   END SUBROUTINE sr_lsd_eval

! **************************************************************************************************
!> \brief ...
!> \param rhoa ...
!> \param rhob ...
!> \param omega ...
!> \param rho_cutoff ...
!> \param e_0 ...
!> \param ea ...
!> \param eb ...
!> \param npoints ...
!> \param order ...
!> \param sx ...
!> \param sc ...
! **************************************************************************************************
   SUBROUTINE sr_lsd_calc(rhoa, rhob, omega, rho_cutoff, e_0, ea, eb, npoints, order, sx, sc)
      !FM low-level computation routine
      REAL(KIND=dp), DIMENSION(*), INTENT(IN)            :: rhoa, rhob
      REAL(KIND=dp), INTENT(IN)                          :: omega, rho_cutoff
      REAL(KIND=dp), DIMENSION(*), INTENT(INOUT)         :: e_0, ea, eb
      INTEGER, INTENT(IN)                                :: npoints, order
      REAL(KIND=dp), INTENT(IN)                          :: sx, sc

      CHARACTER(len=*), PARAMETER :: routineN = 'sr_lsd_calc', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, k
      REAL(KIND=dp)                                      :: my_rhoa, my_rhob, rho, rs, zeta
      REAL(KIND=dp), DIMENSION(0:5)                      :: ed

      CALL timeset(routineN, handle)

      IF (sc /= 0.0_dp .OR. sx /= 0.0_dp) THEN
!$OMP PARALLEL DO PRIVATE (k, rho, ed, my_rhoa, my_rhob, rs, zeta) DEFAULT(NONE)&
!$OMP SHARED(npoints,rhoa,rhob,rho_cutoff,omega,order,e_0,ea,eb,sx,sc)
         DO k = 1, npoints

            my_rhoa = rhoa(k)
            my_rhob = rhob(k)
            rho = my_rhoa + my_rhob
            IF (rho > rho_cutoff) THEN
               rs = (3.0_dp/4.0_dp/pi/rho)**(1.0_dp/3.0_dp)
               zeta = (my_rhoa - my_rhob)/rho

               CALL lsdsr(rs, zeta, omega, ed(0), ed(1), ed(2), sx, sc)
               IF (order >= 0) THEN
                  e_0(k) = e_0(k) + ed(0)*rho
               END IF
               IF (order >= 1 .OR. order == -1) THEN
                  ea(k) = ea(k) + ed(1)
                  eb(k) = eb(k) + ed(2)
               END IF
            END IF

         END DO
!$OMP END PARALLEL DO

         CALL timestop(handle)

      END IF

   END SUBROUTINE sr_lsd_calc

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param mu ...
!> \param excsr ...
!> \param vxcsr ...
!> \param sx ...
!> \param sc ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE ldasr(rs, mu, excsr, vxcsr, sx, sc)
      REAL(KIND=dp), INTENT(IN)                          :: rs, mu
      REAL(KIND=dp), INTENT(OUT)                         :: excsr, vxcsr
      REAL(KIND=dp), INTENT(IN)                          :: sx, sc

      CHARACTER(len=*), PARAMETER :: routineN = 'ldasr', routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: ec, ecd, eclr, ex, exlr, vc, vclr, vx, &
                                                            vxlr

      IF (sx /= 0.0_dp) THEN
         ex = -3._dp*cf/rs/4._dp/pi
         vx = -(3._dp/2._dp/pi)**(2._dp/3._dp)/rs

         CALL exchangelr_lda(rs, mu, exlr, vxlr)
      ELSE
         ex = 0.0_dp
         vx = 0.0_dp

         exlr = 0.0_dp
         vxlr = 0.0_dp
      END IF

      IF (sc /= 0.0_dp) THEN
         CALL ecPW_lda(rs, ec, ecd)
         vc = ec - rs/3._dp*ecd

         CALL ecorrlr_lda(rs, mu, eclr, vclr, ec, ecd)
      ELSE
         ec = 0.0_dp
         vc = 0.0_dp

         eclr = 0.0_dp
         vclr = 0.0_dp
      END IF

      excsr = sx*ex + sc*ec - (sx*exlr + sc*eclr)
      vxcsr = sx*vx + sc*vc - (sx*vxlr + sc*vclr)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param mu ...
!> \param eclr ...
!> \param vclr ...
!> \param ec ...
!> \param ecd ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE ecorrlr_lda(rs, mu, eclr, vclr, ec, ecd)
      REAL(KIND=dp), INTENT(IN)                          :: rs, mu
      REAL(KIND=dp), INTENT(OUT)                         :: eclr, vclr
      REAL(KIND=dp), INTENT(IN)                          :: ec, ecd

      CHARACTER(len=*), PARAMETER :: routineN = 'ecorrlr_lda', routineP = moduleN//':'//routineN

      REAL(KIND=dp) :: a1, a1rs, a2, a2rs, a3, a3rs, a4, a4rs, a5, a5rs, b0, coe2, coe2rs, coe3, &
                       coe3rs, coe4, coe4rs, coe5, coe5rs, d2anti, d2antid, d3anti, d3antid, eclrrs, x, z

      b0 = adib*rs
      z = 0._dp

      d2anti = (q1a + q2a*rs)*EXP(-q3a*rs)/rs
      d3anti = (t1a + t2a*rs)*EXP(-t3a*rs)/rs**2

      d2antid = -((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/rs**2)*EXP(-q3a*rs)
      d3antid = -((rs*t2a*(1._dp + rs*t3a) + t1a*(2._dp + rs*t3a))/rs**3)*EXP(-rs*t3a)

      coe2 = -3._dp/8._dp/rs**3*(g0(rs) - 0.5_dp)
      coe2rs = -3._dp/8._dp/rs**3*g0d(rs) + 9._dp/8._dp/rs**4*(g0(rs) - 0.5_dp)

      coe3 = -g0(rs)/SQRT(2._dp*pi)/rs**3
      coe3rs = -g0d(rs)/SQRT(2._dp*pi)/rs**3 + 3._dp*g0(rs)/SQRT(2._dp*pi)/rs**4

      coe4 = -9._dp/64._dp/rs**3*(.5_dp*dpol(rs*2._dp**(1._dp/3._dp)) + d2anti - cf**2/5._dp/rs**2)
      coe4rs = -3._dp/rs*coe4 - 9._dp/64._dp/rs**3*(((1._dp + z)/2._dp)**(5._dp/3._dp)*dpold(rs*(2._dp/(1._dp + z))** &
                                                                              (1._dp/3._dp)) + ((1._dp - z)/2._dp)**(5._dp/3._dp)* &
                                                    dpold(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) + (1._dp - z**2)*d2antid &
                                                    + cf**2/5._dp*((1._dp + z)**(8._dp/3._dp) + (1._dp - z)**(8._dp/3._dp))/rs**3)

      coe5 = -9._dp/40._dp/SQRT(2._dp*pi)/rs**3*(((1._dp + z)/2._dp)**2*dpol(rs*(2._dp/(1._dp + z))**(1._dp/3._dp)) &
                                                 + ((1._dp - z)/2._dp)**2*dpol(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) &
                                                 + (1._dp - z**2)*d3anti)
      coe5rs = -3._dp/rs*coe5 - 9._dp/(40._dp*SQRT(2._dp*pi)*rs**3)*( &
               ((1._dp + z)/2._dp)**(5._dp/3._dp)*dpold(rs*(2._dp/(1._dp + z))**(1._dp/3._dp)) &
               + ((1._dp - z)/2._dp)**(5._dp/3._dp)*dpold(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) + (1._dp - z**2)*d3antid)

      a1 = 4._dp*b0**6*coe3 + b0**8*coe5
      a1rs = 24._dp*adib*b0**5*coe3 + 4._dp*b0**6*coe3rs + 8._dp*adib*b0**7*coe5 + b0**8*coe5rs

      a2 = 4._dp*b0**6*coe2 + b0**8*coe4 + 6._dp*b0**4*ec
      a2rs = 24._dp*adib*b0**5*coe2 + 4._dp*b0**6*coe2rs + 8._dp*adib*b0**7*coe4 + b0**8*coe4rs &
             + 24._dp*adib*b0**3*ec + 6._dp*b0**4*ecd

      a3 = b0**8*coe3
      a3rs = 8._dp*adib*b0**7*coe3 + b0**8*coe3rs

      a4 = b0**6*(b0**2*coe2 + 4._dp*ec)
      a4rs = 8._dp*adib*b0**7*coe2 + b0**8*coe2rs + 24._dp*adib*b0**5*ec + 4._dp*b0**6*ecd

      a5 = b0**8*ec
      a5rs = 8._dp*adib*b0**7*ec + b0**8*ecd

      x = mu*SQRT(rs)

      eclr = (Qrpa(x) + mu**3*(a1 + mu*(a2 + mu*(a3 + mu*(a4 + a5*mu**2)))))/((1._dp + b0**2*mu**2)**4)

      eclrrs = -8._dp*adib/(1._dp + b0**2*mu**2)*b0*mu**2*eclr + &
               1._dp/((1._dp + b0**2*mu**2)**4)*(mu/(2._dp*SQRT(rs))*Qrpad(x) + &
                                                 mu**3*(a1rs + mu*(a2rs + mu*(a3rs + mu*(a4rs + a5rs*mu**2)))))

      vclr = eclr - rs/3._dp*eclrrs

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param z ...
!> \param mu ...
!> \param vxlrup ...
!> \param vxlrdown ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE exchangelr_lda(rs, mu, exlr, vxlr)
      REAL(KIND=dp), INTENT(IN)                          :: rs, mu
      REAL(KIND=dp), INTENT(OUT)                         :: exlr, vxlr

      CHARACTER(len=*), PARAMETER :: routineN = 'exchangelr_lda', &
                                     routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: derrs, fx, fx1, y

      y = alpha/2._dp*mu*rs
      fx = -((y*(-3._dp + 4._dp*y**2 + (2._dp - 4._dp*y**2)*EXP(-.25_dp/y**2)) + SQRT(pi)*ERF(.5_dp/y))/pi)
      exlr = mu*fx
      fx1 = (3._dp*(1._dp + (-4._dp + 4._dp*EXP(-.25_dp/y**2))*y**2))/pi
      derrs = alpha/4._dp*mu**2*fx1
      vxlr = 2._dp/3._dp*rs*derrs

      vxlr = exlr - vxlr

   END SUBROUTINE

! **************************************************************************************************
!> \brief PW92 energy functional
!> \param rs ...
!> \param ec ...
!> \param ecd ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE ecPW_lda(rs, ec, ecd)
      REAL(KIND=dp), INTENT(IN)                          :: rs
      REAL(KIND=dp), INTENT(OUT)                         :: ec, ecd

      REAL(KIND=dp) :: G(0:1)

      CALL calc_g(rs, 0, G, 1)
      ec = G(0)
      ecd = G(1)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param z ...
!> \param mu ...
!> \param excsr ...
!> \param vxcsrup ...
!> \param vxcsrdown ...
!> \param sx ...
!> \param sc ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE lsdsr(rs, z, mu, excsr, vxcsrup, vxcsrdown, sx, sc)
      REAL(KIND=dp), INTENT(IN)                          :: rs, z, mu
      REAL(KIND=dp), INTENT(OUT)                         :: excsr, vxcsrup, vxcsrdown
      REAL(KIND=dp), INTENT(IN)                          :: sx, sc

      CHARACTER(len=*), PARAMETER :: routineN = 'lsdsr', routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: ec, ecd, eclr, ecz, ex, exlr, vcdown, &
                                                            vclrdown, vclrup, vcup, vxdown, &
                                                            vxlrdown, vxlrup, vxup

      IF (sx /= 0.0_dp) THEN
         ex = -3._dp*cf/rs/8._dp/pi*((1._dp + z)**(4._dp/3._dp) + &
                                     (1._dp - z)**(4._dp/3._dp))

         vxup = -(1._dp + z)**(1._dp/3._dp)*(3._dp/2._dp/pi)**(2._dp/3._dp)/rs
         vxdown = -(1._dp - z)**(1._dp/3._dp)*(3._dp/2._dp/pi)**(2._dp/3._dp)/rs

         CALL exchangelr_lsd(rs, z, mu, exlr, vxlrup, vxlrdown)
      ELSE
         ex = 0.0_dp
         vxup = 0.0_dp
         vxdown = 0.0_dp

         exlr = 0.0_dp
         vxlrup = 0.0_dp
         vxlrdown = 0.0_dp
      END IF

      IF (sc /= 0.0_dp) THEN
         CALL ecPW_lsd(rs, z, ec, ecd, ecz)
         vcup = ec - rs/3._dp*ecd - (z - 1._dp)*ecz
         vcdown = ec - rs/3._dp*ecd - (z + 1._dp)*ecz

         CALL ecorrlr_lsd(rs, z, mu, eclr, vclrup, vclrdown, ec, ecd, ecz)
      ELSE
         ec = 0.0_dp
         vcup = 0.0_dp
         vcdown = 0.0_dp

         eclr = 0.0_dp
         vclrup = 0.0_dp
         vclrdown = 0.0_dp
      END IF

      excsr = sx*ex + sc*ec - (sx*exlr + sc*eclr)
      vxcsrup = sx*vxup + sc*vcup - (sx*vxlrup + sc*vclrup)
      vxcsrdown = sx*vxdown + sc*vcdown - (sx*vxlrdown + sc*vclrdown)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param z ...
!> \param mu ...
!> \param eclr ...
!> \param ec ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE ecorrlr_lsd(rs, z, mu, eclr, vclrup, vclrdown, ec, ecd, ecz)
      REAL(KIND=dp), INTENT(IN)                          :: rs, z, mu
      REAL(KIND=dp), INTENT(OUT)                         :: eclr, vclrup, vclrdown
      REAL(KIND=dp), INTENT(IN)                          :: ec, ecd, ecz

      CHARACTER(len=*), PARAMETER :: routineN = 'ecorrlr_lsd', routineP = moduleN//':'//routineN

      REAL(KIND=dp) :: a1, a1rs, a1z, a2, a2rs, a2z, a3, a3rs, a3z, a4, a4rs, a4z, a5, a5rs, a5z, &
                       b0, coe2, coe2rs, coe2z, coe3, coe3rs, coe3z, coe4, coe4rs, coe4z, coe5, coe5rs, coe5z, &
                       d2anti, d2antid, d3anti, d3antid, eclrrs, eclrz, phi, x

      phi = ((1._dp + z)**(2._dp/3._dp) + (1._dp - z)**(2._dp/3._dp))/2._dp

      b0 = adib*rs

      d2anti = (q1a + q2a*rs)*EXP(-q3a*rs)/rs
      d3anti = (t1a + t2a*rs)*EXP(-t3a*rs)/rs**2

      d2antid = -((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/rs**2)*EXP(-q3a*rs)
      d3antid = -((rs*t2a*(1._dp + rs*t3a) + t1a*(2._dp + rs*t3a))/rs**3)*EXP(-rs*t3a)

      coe2 = -3._dp/8._dp/rs**3*(1._dp - z**2)*(g0(rs) - 0.5_dp)
      coe2rs = -3._dp/8._dp/rs**3*(1._dp - z**2)*g0d(rs) + 9._dp/8._dp/rs**4*(1._dp - z**2)*(g0(rs) - 0.5_dp)
      coe2z = -3._dp/8._dp/rs**3*(-2._dp*z)*(g0(rs) - 0.5_dp)

      coe3 = -(1._dp - z**2)*g0(rs)/SQRT(2._dp*pi)/rs**3
      coe3rs = -(1._dp - z**2)*g0d(rs)/SQRT(2._dp*pi)/rs**3 + 3._dp*(1._dp - z**2)*g0(rs)/SQRT(2._dp*pi)/rs**4
      coe3z = 2._dp*z*g0(rs)/(SQRT(2._dp*pi)*rs**3)

      IF (ABS(z) >= 1._dp) THEN

         coe4 = -9._dp/64._dp/rs**3*(dpol(rs) - cf**2*2**(5._dp/3._dp)/5._dp/rs**2)
         coe4rs = -3._dp/rs*coe4 - 9._dp/64._dp/rs**3*(dpold(rs) + 2._dp*cf**2*2**(5._dp/3._dp)/5._dp/rs**3)
         coe4z = -9._dp/64._dp/rs**3*(dpol(rs) - rs/6._dp*dpold(rs) - 2._dp*d2anti &
                                      - 4._dp/15._dp*cf**2*2._dp**(5._dp/3._dp)/rs**2)*z

         coe5 = -9._dp/40._dp/SQRT(2._dp*pi)/rs**3*dpol(rs)
         coe5rs = -3._dp/rs*coe5 - 9._dp/40._dp/SQRT(2._dp*pi)/rs**3*dpold(rs)
         coe5z = -9._dp/40._dp/SQRT(2._dp*pi)/rs**3*(dpol(rs) - rs/6._dp*dpold(rs) - 2._dp*d3anti)*z

      ELSE

         coe4 = -9._dp/64._dp/rs**3*(((1._dp + z)/2._dp)**2* &
                                     dpol(rs*(2._dp/(1._dp + z))**(1._dp/3._dp)) + ((1._dp - z)/2._dp)**2 &
                                     *dpol(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) + &
                                     (1._dp - z**2)*d2anti - cf**2/10._dp*((1._dp + z)**(8._dp/3._dp) &
                                                                           + (1._dp - z)**(8._dp/3._dp))/rs**2)
         coe4rs = -3._dp/rs*coe4 - 9._dp/64._dp/rs**3*(((1._dp + z)/2._dp)**(5._dp/3._dp)*dpold(rs*(2._dp/(1._dp + z))** &
                                                                              (1._dp/3._dp)) + ((1._dp - z)/2._dp)**(5._dp/3._dp)* &
                                                       dpold(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) + (1._dp - z**2)*d2antid &
                                                      + cf**2/5._dp*((1._dp + z)**(8._dp/3._dp) + (1._dp - z)**(8._dp/3._dp))/rs**3)
         coe4z = -9._dp/64._dp/rs**3*(1._dp/2._dp*(1._dp + z)*dpol(rs*(2._dp/(1._dp + z))**(1._dp/3._dp)) &
                                      - 1._dp/2._dp*(1._dp - z)*dpol(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) &
                                      - rs/6._dp*((1._dp + z)/2._dp)**(2._dp/3._dp)*dpold(rs*(2/(1._dp + z))**(1._dp/3._dp)) &
                                      + rs/6._dp*((1._dp - z)/2._dp)**(2._dp/3._dp)*dpold(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) &
                                      - 2._dp*z*d2anti - 4._dp/15._dp*cf**2/rs**2*((1._dp + z)**(5._dp/3._dp) &
                                                                                   - (1._dp - z)**(5._dp/3._dp)))

         coe5 = -9._dp/40._dp/SQRT(2._dp*pi)/rs**3*(((1._dp + z)/2._dp)**2*dpol(rs*(2._dp/(1._dp + z))**(1._dp/3._dp)) &
                                                    + ((1._dp - z)/2._dp)**2*dpol(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) &
                                                    + (1._dp - z**2)*d3anti)
         coe5rs = -3._dp/rs*coe5 - 9._dp/(40._dp*SQRT(2._dp*pi)*rs**3)*( &
                  ((1._dp + z)/2._dp)**(5._dp/3._dp)*dpold(rs*(2._dp/(1._dp + z))**(1._dp/3._dp)) &
                  + ((1._dp - z)/2._dp)**(5._dp/3._dp)*dpold(rs*(2._dp/(1._dp - z))**(1._dp/3._dp)) + (1._dp - z**2)*d3antid)
         coe5z = -9._dp/40._dp/SQRT(2._dp*pi)/rs**3*(1._dp/2._dp*(1._dp + z)*dpol(rs*(2/(1._dp + z))**(1._dp/3._dp)) &
                                                     - 1._dp/2._dp*(1._dp - z)*dpol(rs*(2/(1._dp - z))**(1._dp/3._dp)) &
                                                     - rs/6._dp*((1._dp + z)/2._dp)**(2._dp/3._dp)*dpold(rs*(2/(1._dp + z)) &
                                                                    **(1._dp/3._dp)) + rs/6._dp*((1._dp - z)/2._dp)**(2._dp/3._dp) &
                                                     *dpold(rs*(2/(1._dp - z))**(1._dp/3._dp)) - 2._dp*z*d3anti)

      END IF

      a1 = 4._dp*b0**6*coe3 + b0**8*coe5
      a1rs = 24._dp*adib*b0**5*coe3 + 4._dp*b0**6*coe3rs + 8._dp*adib*b0**7*coe5 + b0**8*coe5rs
      a1z = 4._dp*b0**6*coe3z + b0**8*coe5z

      a2 = 4._dp*b0**6*coe2 + b0**8*coe4 + 6._dp*b0**4*ec
      a2rs = 24._dp*adib*b0**5*coe2 + 4._dp*b0**6*coe2rs + 8._dp*adib*b0**7*coe4 + b0**8*coe4rs &
             + 24._dp*adib*b0**3*ec + 6._dp*b0**4*ecd
      a2z = 4._dp*b0**6*coe2z + b0**8*coe4z + 6._dp*b0**4*ecz

      a3 = b0**8*coe3
      a3rs = 8._dp*adib*b0**7*coe3 + b0**8*coe3rs
      a3z = b0**8*coe3z

      a4 = b0**6*(b0**2*coe2 + 4._dp*ec)
      a4rs = 8._dp*adib*b0**7*coe2 + b0**8*coe2rs + 24._dp*adib*b0**5*ec + 4._dp*b0**6*ecd
      a4z = b0**6*(b0**2*coe2z + 4._dp*ecz)

      a5 = b0**8*ec
      a5rs = 8._dp*adib*b0**7*ec + b0**8*ecd
      a5z = b0**8*ecz

      x = mu*SQRT(rs)/phi

      eclr = (phi**3*Qrpa(x) + mu**3*(a1 + mu*(a2 + mu*(a3 + mu*(a4 + a5*mu**2)))))/((1._dp + b0**2*mu**2)**4)

      eclrrs = -8._dp*adib/(1._dp + b0**2*mu**2)*b0*mu**2*eclr + &
               1._dp/((1._dp + b0**2*mu**2)**4)*(phi**2*mu/(2._dp*SQRT(rs))*Qrpad(x) + &
                                                 mu**3*(a1rs + mu*(a2rs + mu*(a3rs + mu*(a4rs + a5rs*mu**2)))))

      IF (z >= 1._dp) THEN
         vclrup = eclr - rs/3._dp*eclrrs
         vclrdown = 0._dp
      ELSE IF (z <= -1._dp) THEN
         vclrup = 0._dp
         vclrdown = eclr - rs/3._dp*eclrrs
      ELSE

         eclrz = (phi**2*((1._dp + z)**(-1._dp/3._dp) - (1._dp - z)**(-1._dp/3._dp)) &
                  *Qrpa(x) - phi*Qrpad(x)*mu*SQRT(rs)*((1._dp + z)**(-1._dp/3._dp) &
                                                       - (1._dp - z)**(-1._dp/3._dp))/3._dp + &
                  mu**3*(a1z + mu*(a2z + mu*(a3z + mu*(a4z + a5z*mu**2)))))/((1._dp + b0**2*mu**2)**4)

         vclrup = eclr - rs/3._dp*eclrrs - (z - 1._dp)*eclrz
         vclrdown = eclr - rs/3._dp*eclrrs - (z + 1._dp)*eclrz
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION g0(rs)
      REAL(KIND=dp), INTENT(IN)                          :: rs
      REAL(KIND=dp)                                      :: g0

      g0 = (1._dp - (0.7317_dp - Fg0)*rs + Cg0*rs**2 + Dg0*rs**3 + Eg0*rs**4)*EXP(-ABS(Fg0)*rs)/2._dp

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION g0d(rs)
      REAL(KIND=dp), INTENT(IN)                          :: rs
      REAL(KIND=dp)                                      :: g0d

      g0d = (-Bg0 + 2.0_dp*Cg0*rs + 3.0_dp*Dg0*rs**2 + 4.0_dp*Eg0*rs**3)/2._dp*EXP(-Fg0*rs) &
            - (Fg0*(1.0_dp - Bg0*rs + Cg0*rs**2 + Dg0*rs**3 + Eg0*rs**4))/ &
            2._dp*EXP(-Fg0*rs)

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION dpol(rs)
      REAL(KIND=dp), INTENT(IN)                          :: rs
      REAL(KIND=dp)                                      :: dpol

      dpol = 2._dp**(5._dp/3._dp)/5._dp*cf**2/rs**2*(1._dp + (p3p - 0.454555_dp)*rs) &
             /(1._dp + p3p*rs + p2p*rs**2)

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION dpold(rs)
      REAL(KIND=dp), INTENT(IN)                          :: rs
      REAL(KIND=dp)                                      :: dpold

      dpold = 2._dp**(5._dp/3._dp)/5._dp*cf**2* &
              (-2._dp + (0.454555 - 4._dp*p3p)*rs + &
               (-4._dp*p2p + (0.90911 - 2.*p3p)*p3p)*rs**2 &
               + p2p*(1.363665 - 3._dp*p3p)*rs**3)/ &
              (rs**3*(1._dp + p3p*rs + p2p*rs**2)**2)

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION Qrpa(x)
      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: Qrpa

      Qrpa = Acoul*LOG((1._dp + aQ2*x + bQ2*x**2 + cQ2*x**3)/(1._dp + aQ2*x + dQ2*x**2))

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param x ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION Qrpad(x)
      REAL(KIND=dp), INTENT(IN)                          :: x
      REAL(KIND=dp)                                      :: Qrpad

      Qrpad = Acoul*((x*(bQ2*(2._dp + aQ2*x) + cQ2*x*(3._dp + 2._dp*aQ2*x) + dQ2*(-2._dp - aQ2*x + cQ2*x**3)))/ &
                     ((1._dp + aQ2*x + dQ2*x**2)*(1._dp + aQ2*x + bQ2*x**2 + cQ2*x**3)))

   END FUNCTION

! **************************************************************************************************
!> \brief ...
!> \param rs ...
!> \param z ...
!> \param mu ...
!> \param vxlrup ...
!> \param vxlrdown ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE exchangelr_lsd(rs, z, mu, exlr, vxlrup, vxlrdown)
      REAL(KIND=dp), INTENT(IN)                          :: rs, z, mu
      REAL(KIND=dp), INTENT(OUT)                         :: exlr, vxlrup, vxlrdown

      CHARACTER(len=*), PARAMETER :: routineN = 'exchangelr_lsd', &
                                     routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: derrs, derz, fx, fx1, x, y

      IF (z >= 1._dp) THEN
         x = rs*alpha*mu
         y = .5_dp**(4._dp/3._dp)*x
         fx = -((y*(-3._dp + 4._dp*y**2 + (2._dp - 4._dp*y**2)*EXP(-.25/y**2)) + SQRT(pi)*ERF(.5_dp/y))/pi)
         exlr = mu*fx
         vxlrup = mu*(x/(2._dp**(1._dp/3._dp)*pi) - x/(2._dp**(1._dp/3._dp)*pi)* &
                      EXP(-2._dp**(2._dp/3._dp)/x**2) - &
                      ERF(2._dp**(1._dp/3._dp)/x)/SQRT(pi))
         vxlrdown = 0._dp
      ELSE IF (z <= -1._dp) THEN
         x = rs*alpha*mu
         y = .5_dp**(4._dp/3._dp)*x
         fx = -((y*(-3._dp + 4._dp*y**2 + (2._dp - 4._dp*y**2)*EXP(-.25/y**2)) + SQRT(pi)*ERF(.5_dp/y))/pi)
         exlr = mu*fx
         vxlrdown = mu*(x/(2._dp**(1._dp/3._dp)*pi) - x/(2._dp**(1._dp/3._dp)*pi)* &
                        EXP(-2._dp**(2._dp/3._dp)/x**2) - &
                        ERF(2._dp**(1._dp/3._dp)/x)/SQRT(pi))
         vxlrup = 0._dp
      ELSE
         y = alpha/2._dp/(1.+z)**(1._dp/3._dp)*mu*rs
         fx = -((y*(-3._dp + 4._dp*y**2 + (2._dp - 4._dp*y**2)*EXP(-.25_dp/y**2)) + &
                 SQRT(pi)*ERF(.5_dp/y))/pi)
         exlr = (1._dp + z)*mu*fx/2._dp
         fx1 = (3._dp*(1._dp + (-4._dp + 4._dp*EXP(-.25_dp/y**2))*y**2))/pi
         derrs = alpha/4._dp*(1._dp + z)**(2._dp/3._dp)*mu**2*fx1
         derz = 1._dp/2._dp*mu*fx - 1._dp/6._dp*fx1*mu*y
         vxlrup = rs/3._dp*derrs + (z - 1._dp)*derz
         vxlrdown = rs/3._dp*derrs + (z + 1._dp)*derz

         y = alpha/2._dp/(1.-z)**(1._dp/3._dp)*mu*rs
         fx = -((y*(-3._dp + 4._dp*y**2 + (2._dp - 4._dp*y**2)*EXP(-.25_dp/y**2)) + &
                 SQRT(pi)*ERF(.5_dp/y))/pi)
         exlr = exlr + (1._dp - z)*mu*fx/2._dp
         fx1 = (3._dp*(1._dp + (-4._dp + 4._dp*EXP(-.25_dp/y**2))*y**2))/pi
         derrs = alpha/4._dp*(1._dp - z)**(2._dp/3._dp)*mu**2*fx1
         derz = -1._dp/2._dp*mu*fx + 1._dp/6._dp*fx1*mu*y
         vxlrup = vxlrup + rs/3._dp*derrs + (z - 1._dp)*derz
         vxlrdown = vxlrdown + rs/3._dp*derrs + (z + 1._dp)*derz

         vxlrup = exlr - vxlrup
         vxlrdown = exlr - vxlrdown
      ENDIF

   END SUBROUTINE

! **************************************************************************************************
!> \brief PW92 energy functional
!> \param rs ...
!> \param z ...
!> \param ec ...
!> \param ecd ...
!> \param ecz ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE ecPW_lsd(rs, z, ec, ecd, ecz)
      REAL(KIND=dp), INTENT(IN)                          :: rs, z
      REAL(KIND=dp), INTENT(OUT)                         :: ec, ecd, ecz

      CHARACTER(len=*), PARAMETER :: routineN = 'ecPW_lsd', routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: alfac(0:1), ec0(0:1), ec1(0:1), ff

      IF (ABS(z) >= 1._dp) THEN
         CALL calc_g(rs, 0, ec0, 0)
         CALL calc_g(rs, 1, ec1, 1)
         CALL calc_g(rs, -1, alfac, 0)
         alfac = -alfac

         ec = ec1(0)
         ecd = ec1(1)
         ecz = SIGN(-4._dp/f02*alfac(0) + (ec1(0) - ec0(0))*(4._dp + 2._dp**(4._dp/3._dp)/3._dp/ &
                                                             (2._dp**(1._dp/3._dp) - 1._dp)), z)
      ELSE
         ff = ((1._dp + z)**(4._dp/3._dp) + (1._dp - z)**(4._dp/3._dp) - &
               2._dp)/(2._dp**(4._dp/3._dp) - 2._dp)

         CALL calc_g(rs, 0, ec0, 1)
         CALL calc_g(rs, 1, ec1, 1)
         CALL calc_g(rs, -1, alfac, 1)
         alfac = -alfac

         ec = ec0(0) + alfac(0)*ff/f02*(1._dp - z**4) + (ec1(0) - ec0(0))*ff*z**4
         ecd = ec0(1) + alfac(1)*ff/f02*(1._dp - z**4) + (ec1(1) - ec0(1))*ff*z**4
         ecz = alfac(0)*(-4._dp*z**3)*ff/f02 + alfac(0)*(1._dp - z**4)/f02* &
               4._dp/3._dp*((1._dp + z)**(1._dp/3._dp) - (1._dp - z)**(1._dp/3._dp))/ &
               (2._dp**(4._dp/3._dp) - 2._dp) + (ec1(0) - ec0(0))*(4._dp*z**3*ff + &
                                                            4._dp/3._dp*((1._dp + z)**(1._dp/3._dp) - (1._dp - z)**(1._dp/3._dp))/ &
                                                                   (2._dp**(4._dp/3._dp) - 2._dp)*z**4)
      END IF

   END SUBROUTINE

@:calc_g()

END MODULE xc_sr_lda
