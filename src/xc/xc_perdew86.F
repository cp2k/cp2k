!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculate the Perdew Correlation from 1986
!> \par History
!>      JGH (26.02.2003) : OpenMP enabled
!>      fawzi (04.2004)  : adapted to the new xc interface
!> \author JGH (03.03.2002)
! *****************************************************************************
MODULE xc_perdew86

  
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_get_derivative
  USE xc_derivative_types,             ONLY: xc_derivative_get,&
                                             xc_derivative_type
  USE xc_functionals_utilities,        ONLY: calc_rs_pw,&
                                             set_util
  USE xc_input_constants,              ONLY: pz_orig
  USE xc_perdew_zunger,                ONLY: pz_lda_eval
  USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
  USE xc_rho_set_types,                ONLY: xc_rho_set_get,&
                                             xc_rho_set_type
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  REAL(KIND=dp), PARAMETER :: pi = 3.14159265358979323846264338_dp
  REAL(KIND=dp), PARAMETER :: f13 = 1.0_dp/3.0_dp, &
                          f23 = 2.0_dp*f13, &
                          f43 = 4.0_dp*f13, &
                          f53 = 5.0_dp*f13, &
                          f76 = 7.0_dp/6.0_dp, &
                          frs = 1.6119919540164696407_dp, &
                          fpe = 0.19199566167376364_dp

  PUBLIC :: p86_lda_info, p86_lda_eval

  REAL(KIND=dp) :: eps_rho
  LOGICAL :: debug_flag

  REAL(KIND=dp), PARAMETER :: a = 0.023266_dp, &
                          b = 7.389e-6_dp, &
                          c = 8.723_dp, &
                          d = 0.472_dp, &
                          pc1 = 0.001667_dp, &
                          pc2 = 0.002568_dp, &
                          pci = pc1+pc2
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_perdew86'

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param cutoff ...
!> \param debug ...
! *****************************************************************************
  SUBROUTINE p86_init ( cutoff, debug )

    REAL(KIND=dp), INTENT(IN)                :: cutoff
    LOGICAL, INTENT(IN), OPTIONAL            :: debug

    eps_rho = cutoff
    CALL set_util ( cutoff )

    IF ( PRESENT ( debug ) ) THEN
       debug_flag = debug
    ELSE
       debug_flag = .FALSE.
    ENDIF

  END SUBROUTINE p86_init

! *****************************************************************************
!> \brief ...
!> \param reference ...
!> \param shortform ...
!> \param needs ...
!> \param max_deriv ...
! *****************************************************************************
  SUBROUTINE p86_lda_info ( reference, shortform, needs, max_deriv)
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: reference, shortform
    TYPE(xc_rho_cflags_type), &
      INTENT(inout), OPTIONAL                :: needs
    INTEGER, INTENT(out), OPTIONAL           :: max_deriv

    IF ( PRESENT ( reference ) ) THEN
      reference = "J. P. Perdew, Phys. Rev. B, 33, 8822 (1986) {LDA version}"
    END IF
    IF ( PRESENT ( shortform ) ) THEN
      shortform = "Perdew 1986 correlation energy functional {LDA}"
    END IF
    IF ( PRESENT(needs) ) THEN
       needs%rho=.TRUE.
       needs%norm_drho=.TRUE.
    END IF
    IF (PRESENT(max_deriv)) max_deriv=3

  END SUBROUTINE p86_lda_info

! *****************************************************************************
!> \brief ...
!> \param rho_set ...
!> \param deriv_set ...
!> \param order ...
!> \param p86_params ...
! *****************************************************************************
  SUBROUTINE p86_lda_eval ( rho_set, deriv_set, order, p86_params)


    TYPE(xc_rho_set_type), POINTER           :: rho_set
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    INTEGER, INTENT(IN)                      :: order
    TYPE(section_vals_type), POINTER         :: p86_params

    CHARACTER(len=*), PARAMETER :: routineN = 'p86_lda_eval', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, m, npoints, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bo
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: drho_cutoff, rho_cutoff
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: rs
    REAL(KIND=dp), DIMENSION(:, :, :), POINTER :: e_0, e_ndrho, &
      e_ndrho_ndrho, e_ndrho_ndrho_ndrho, e_rho, e_rho_ndrho, &
      e_rho_ndrho_ndrho, e_rho_rho, e_rho_rho_ndrho, e_rho_rho_rho, grho, rho
    TYPE(xc_derivative_type), POINTER        :: deriv

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(bo,rho,e_0,e_rho,e_ndrho,&
         e_rho_rho, e_rho_ndrho, e_ndrho_ndrho,&
         e_rho_rho_rho, e_rho_rho_ndrho, e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho)
    CPASSERT(ASSOCIATED(rho_set))
    CPASSERT(rho_set%ref_count>0)
    CPASSERT(ASSOCIATED(deriv_set))
    CPASSERT(deriv_set%ref_count>0)

    ! calculate the perdew_zunger correlation
    CALL pz_lda_eval( pz_orig, rho_set, deriv_set, order, p86_params)

    CALL xc_rho_set_get(rho_set,rho=rho,&
         norm_drho=grho,local_bounds=bo,rho_cutoff=rho_cutoff,&
         drho_cutoff=drho_cutoff)
    npoints=(bo(2,1)-bo(1,1)+1)*(bo(2,2)-bo(1,2)+1)*(bo(2,3)-bo(1,3)+1)
    CALL p86_init(rho_cutoff)
    m = ABS(order)

    ALLOCATE ( rs(npoints), STAT=stat )
    CPASSERT(stat==0)

    CALL calc_rs_pw ( rho, rs, npoints )
    IF ( order>=0 ) THEN
       deriv => xc_dset_get_derivative(deriv_set,"",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_0)

       CALL p86_u_0 ( rho, rs, grho, e_0, npoints)
    END IF
    IF ( order>=1 .OR. order==-1 ) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho)
       deriv => xc_dset_get_derivative(deriv_set,"(norm_drho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho)

       CALL p86_u_1 ( rho, grho, rs, e_rho,&
            e_ndrho, npoints)
    END IF
    IF ( order>=2 .OR. order==-2 ) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho)
       deriv => xc_dset_get_derivative(deriv_set,"(rho)(norm_drho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drho)(norm_drho)", allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho)

       CALL p86_u_2 ( rho, grho, rs, e_rho_rho,&
            e_rho_ndrho, e_ndrho_ndrho, npoints)
    END IF
    IF ( order>=3 .OR. order ==-3 ) THEN
       deriv => xc_dset_get_derivative(deriv_set,"(rho)(rho)(rho)",&
            allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_rho)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(rho)(rho)(norm_drho)",allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_rho_ndrho)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(rho)(norm_drho)(norm_drho)",allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_rho_ndrho_ndrho)
       deriv => xc_dset_get_derivative(deriv_set,&
            "(norm_drho)(norm_drho)(norm_drho)", allocate_deriv=.TRUE.)
       CALL xc_derivative_get(deriv,deriv_data=e_ndrho_ndrho_ndrho)

       CALL p86_u_3 ( rho, grho, rs, e_rho_rho_rho,&
            e_rho_rho_ndrho, e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho,&
            npoints)
    END IF
    IF ( order>3.OR.order<-3) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="derivatives bigger than 3 not implemented", &
            error_level=cp_failure_level)
    END IF
    DEALLOCATE ( rs, STAT=stat )
    CPASSERT(stat==0)
    CALL timestop(handle)

  END SUBROUTINE p86_lda_eval

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param rs ...
!> \param grho ...
!> \param e_0 ...
!> \param npoints ...
! *****************************************************************************
  SUBROUTINE p86_u_0 ( rho, rs, grho, e_0, npoints)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, rs, grho
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_0
    INTEGER, INTENT(in)                      :: npoints

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: cr, ep, g, or, phi, r, x

!$omp parallel do private(ip,g,r,x,or,cr,phi,ep) default(none)&
!$omp shared(npoints,rho,eps_rho,grho,rs,e_0)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         g = grho(ip)
         r = rs(ip)
         x = r*frs
         or = 1.0_dp/rho(ip)
         cr= pc1 + (pc2+a*r+b*r*r)/(1.0_dp+c*r+d*r*r+1.e4_dp*b*r*r*r)
         phi = fpe*pci/cr * g * SQRT(x)*or
         ep = EXP(-phi)
         e_0(ip) = e_0(ip) + x*or*g*g*cr*ep
      END IF

    END DO

  END SUBROUTINE p86_u_0

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param rs ...
!> \param e_rho ...
!> \param e_ndrho ...
!> \param npoints ...
! *****************************************************************************
  SUBROUTINE p86_u_1 ( rho, grho, rs, e_rho, e_ndrho, npoints)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, grho, rs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho, e_ndrho
    INTEGER, INTENT(in)                      :: npoints

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: cr, dcr, dphig, dphir, dpv, &
                                                dq, ep, ff, g, or, p, phi, q, &
                                                r, x

!$omp parallel do private(ip,g,r,x,or,p,dpv,q,dq,cr,dcr,dphig,phi,dphir,ep,ff) default(none)&
!$omp shared(npoints,rho,eps_rho,grho,rs,e_rho,e_ndrho)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         g = grho(ip)
         r = rs(ip)
         x = r*frs
         or = 1.0_dp/rho(ip)
         p = pc2+a*r+b*r*r
         dpv = a+2.0_dp*b*r
         q = 1.0_dp+c*r+d*r*r+1.e4_dp*b*r*r*r
         dq = c+2.0_dp*d*r+3.e4_dp*b*r*r
         cr= pc1 + p/q
         dcr = (dpv*q-p*dq)/(q*q) * (-f13*r*or)
         dphig = fpe*pci/cr * SQRT(x)*or
         phi = dphig * g
         dphir = -phi*(dcr/cr + f76*or)
         ep = EXP(-phi)
         ff = x*or*g*ep
         e_rho(ip) = e_rho(ip) + ff*g*dcr-ff*g*cr*dphir-ff*g*cr*f43*or
         e_ndrho(ip) = e_ndrho(ip) + ff*cr*(2.0_dp-g*dphig)
      END IF

    END DO

  END SUBROUTINE p86_u_1

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param rs ...
!> \param e_rho_rho ...
!> \param e_rho_ndrho ...
!> \param e_ndrho_ndrho ...
!> \param npoints ...
! *****************************************************************************
  SUBROUTINE p86_u_2 ( rho, grho, rs, e_rho_rho, e_rho_ndrho, e_ndrho_ndrho,&
       npoints)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, grho, rs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho_rho, e_rho_ndrho, &
                                                e_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints

    INTEGER                                  :: ip
    REAL(KIND=dp)                            :: cr, d2cr, d2p, d2phir, d2q, &
                                                dcr, dphig, dphigr, dphir, &
                                                dpv, dq, ep, g, or, p, phi, &
                                                q, r, x

!$omp parallel do private(ip,x,r,cr,phi,ep,g,or,p,q,dpv,dq,dphir,dcr) &
!$omp             private(dphig,dphigr,d2phir,d2cr,d2p,d2q) default(none) &
!$omp shared(npoints,rho,eps_rho,grho,rs,e_rho_rho,e_rho_ndrho,e_ndrho_ndrho)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         g = grho(ip)
         r = rs(ip)
         x = r*frs
         or = 1.0_dp/rho(ip)
         p = pc2+a*r+b*r*r
         dpv = a+2.0_dp*b*r
         d2p = 2.0_dp*b
         q = 1.0_dp+c*r+d*r*r+1.e4_dp*b*r*r*r
         dq = c+2.0_dp*d*r+3.e4_dp*b*r*r
         d2q = 2.0_dp*d+6.e4_dp*b*r
         cr= pc1 + p/q
         dcr = (dpv*q-p*dq)/(q*q) * (-f13*r*or)
         d2cr = (d2p*q*q-p*q*d2q-2*dpv*dq*q+2*p*dq*dq)/(q*q*q)*(f13*r*or)**2 + &
                (dpv*q-p*dq)/(q*q) * f13*f43*r*or*or
         dphig = fpe*pci/cr * SQRT(x)*or
         phi = dphig * g
         dphir = -phi*(dcr/cr + f76*or)
         d2phir = -dphir*(dcr/cr + f76*or)-&
                   phi*((d2cr*cr-dcr*dcr)/(cr*cr) - f76*or*or)
         dphigr = -dphig*(dcr/cr + f76*or)
         ep = EXP(-phi)
         e_rho_rho(ip) = e_rho_rho(ip) + x*or*ep*g*g*&
              (-f43*or*dcr+d2cr-dcr*dphir + &
              f43*or*cr*dphir-dcr*dphir-cr*d2phir+cr*dphir*dphir + &
              f43*or*(7.*f13*or*cr-dcr+cr*dphir))
         e_rho_ndrho(ip) = e_rho_ndrho(ip) + x*or*ep*g*&
              (-2*f43*cr*or+2*dcr-2*cr*dphir+f43*or*g*cr*dphig-&
              g*dcr*dphig+g*cr*dphir*dphig-g*cr*dphigr)
         e_ndrho_ndrho(ip) = e_ndrho_ndrho(ip) + x*or*ep*cr*&
              (2.0_dp-4.0_dp*g*dphig+g*g*dphig*dphig)
      END IF

    END DO

  END SUBROUTINE p86_u_2

! *****************************************************************************
!> \brief ...
!> \param rho ...
!> \param grho ...
!> \param rs ...
!> \param e_rho_rho_rho ...
!> \param e_rho_rho_ndrho ...
!> \param e_rho_ndrho_ndrho ...
!> \param e_ndrho_ndrho_ndrho ...
!> \param npoints ...
! *****************************************************************************
  SUBROUTINE p86_u_3 ( rho, grho, rs, e_rho_rho_rho,&
       e_rho_rho_ndrho, e_rho_ndrho_ndrho, e_ndrho_ndrho_ndrho,&
       npoints)

    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: rho, grho, rs
    REAL(KIND=dp), DIMENSION(*), &
      INTENT(INOUT)                          :: e_rho_rho_rho, &
                                                e_rho_rho_ndrho, &
                                                e_rho_ndrho_ndrho, &
                                                e_ndrho_ndrho_ndrho
    INTEGER, INTENT(in)                      :: npoints

    INTEGER                                  :: ip
    REAL(KIND=dp) :: cr, d2cr, d2p, d2phir, d2phirg, d2pq, d2q, d2z, d3cr, &
      d3phir, d3pq, d3q, d3z, dcr, dphig, dphigr, dphir, dpq, dpv, dq, dz, &
      ep, g, or, oz, p, phi, pq, q, r, x

!$omp parallel do private(ip,x, r, cr, phi, ep, g, or, p, q, dpv, dq, dphir, dcr, dphig) &
!$omp             private(dphigr, d2phir, d3phir, d2cr, d3cr, d2p, d2q, d2phirg, d3q) &
!$omp             private(pq, dpq, d2pq, d3pq, oz, dz, d2z, d3z) default(none) &
!$omp             shared(npoints,rho,eps_rho,grho,e_rho_rho_rho,e_rho_rho_ndrho,e_rho_ndrho_ndrho,e_ndrho_ndrho_ndrho)
    DO ip = 1, npoints

      IF ( rho(ip) > eps_rho ) THEN
         g = grho(ip)
         r = rs(ip)
         x = r*frs
         or = 1.0_dp/rho(ip)
         p = pc2+a*r+b*r*r
         dpv = a+2.0_dp*b*r
         d2p = 2.0_dp*b
         q = 1.0_dp+c*r+d*r*r+1.e4_dp*b*r*r*r
         dq = c+2.0_dp*d*r+3.e4_dp*b*r*r
         d2q = 2.0_dp*d+6.e4*b*r
         d3q = 6.e4*b
         pq = p/q
         dpq = (dpv*q-p*dq)/(q*q)
         d2pq = (d2p*q*q-2*dpv*dq*q+2*p*dq*dq-p*d2q*q)/(q*q*q)
         d3pq = -(3*d2p*dq*q*q-6*dpv*dq*dq*q+3*dpv*d2q*q*q+6*p*dq*dq*dq-6*p*dq*d2q*q &
                  +p*d3q*q*q)/(q*q*q*q)
         cr= pc1 + pq
         dcr = dpq*(-f13*r*or)
         d2cr = d2pq*f13*f13*r*r*or*or + dpq*f13*f43*r*or*or
         d3cr = d3pq*(-f13*r*or)**3 + 3*d2pq*(-f13*f13*f43*r*r*or*or*or) + &
                dpq*(-f13*f43*f13*7*r*or*or*or)
         oz = SQRT(x)*or/cr
         dz = dcr/cr+f76*or
         d2z = d2cr/cr+2*f76*dcr/cr*or+f76/6.*or*or
         d3z = d3cr/cr+3*f76*d2cr/cr*or+3*f76/6.*dcr/cr*or*or-5*f76/36.*or*or*or
         dphig = fpe*pci*oz
         phi = dphig * g
         dphir = -phi*dz
         dphigr = -dphig*dz
         d2phir = -phi*(d2z-2*dz*dz)
         d3phir = -phi*(d3z-6*d2z*dz+6*dz*dz*dz)
         d2phirg = -dphigr*dz-&
                   dphig*((d2cr*cr-dcr*dcr)/(cr*cr) - f76*or*or)
         ep = EXP(-phi)
         e_rho_rho_rho(ip) = e_rho_rho_rho(ip)&
              + g*g*x*or*ep*(-280./27.*or*or*or*cr+3*28./9.*or*or*dcr+&
              3*28./9.*or*or*cr*(-dphir)-4*or*d2cr-8*or*dcr*(-dphir)-&
              4*or*cr*(-d2phir+dphir*dphir)+d3cr+3*d2cr*(-dphir)+&
              3*dcr*(-d2phir+dphir*dphir)+cr*(-d3phir+3*dphir*d2phir-&
              dphir**3))
         e_rho_rho_ndrho(ip)=e_rho_rho_ndrho(ip)&
              + 2.*x*or*ep*g*(-f43*or*dcr+d2cr-dcr*dphir + &
              f43*or*cr*dphir-dcr*dphir-cr*d2phir+cr*dphir*dphir + &
              f43*or*(7.*f13*or*cr-dcr+cr*dphir)) - &
              dphig*x*or*ep*g*g*(-f43*or*dcr+d2cr-dcr*dphir + &
              f43*or*cr*dphir-dcr*dphir-cr*d2phir+cr*dphir*dphir + &
              f43*or*(7.*f13*or*cr-dcr+cr*dphir)) + &
              x*or*ep*g*g*(-dcr*dphigr+f43*or*cr*dphigr-dcr*dphigr-cr*d2phirg+&
              2.*cr*dphigr*dphir+f43*or*cr*dphigr)
         e_rho_ndrho_ndrho(ip) = e_rho_ndrho_ndrho(ip)&
              + x*or*ep*(-2*f43*cr*or+2*dcr-2*cr*dphir+f43*or*g*cr*dphig-&
              g*dcr*dphig+g*cr*dphir*dphig-g*cr*dphigr) + &
              x*or*ep*g*(-2*cr*dphigr+f43*or*cr*dphig-&
              dcr*dphig+cr*dphir*dphig+g*cr*dphigr*dphig-cr*dphigr)-&
              x*or*ep*g*dphig*(-2*f43*cr*or+2*dcr-2*cr*dphir+f43*or*g*cr*dphig-&
              g*dcr*dphig+g*cr*dphir*dphig-g*cr*dphigr)
         e_ndrho_ndrho_ndrho(ip) = e_ndrho_ndrho_ndrho(ip)&
              + x*or*ep*cr*dphig*(-6.0_dp+6.0_dp*g*dphig-g*g*dphig*dphig)
      END IF

    END DO

  END SUBROUTINE p86_u_3

END MODULE xc_perdew86

