!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Includes all necessary routines, functions and parameters from
!>        libxc. Provides CP2K routines/functions where the LibXC calling list
!>        is version dependent. The naming convention for such 
!>        routines/functions is xc_f90_XXX --> 'xc_libxc_wrap_XXX'. All version
!>        independent routines/functions are just bypassed to higher level
!>        module file 'xc_libxc'.
!>        
!> \note For LibXC versions 2.1.x and above.
!>       Once the LibXC-API is stable, remove all 'xc_libxc_wrap_XXX'
!>       routines/functions and use 'xc_f90_lib_m' directly in 'xc_libxc'.
!>       Marques, Oliveira, Burnus, CPC 183, 2272 (2012)).
!>
!> \par History
!>      08.2015 created [A. Gloess]
!> \author A. Gloessa (agloess)
! *****************************************************************************
MODULE xc_libxc_wrap
#if defined(__LIBXC2) || defined(__LIBXC3)
#define __LIBXC
#endif
#if defined (__LIBXC)
#include <xc_version.h>
! check for LibXC version
#if (XC_MAJOR_VERSION < 2) || ((XC_MAJOR_VERSION == 2) && (XC_MINOR_VERSION < 1))
  This version of CP2K ONLY works with libxc versions 2.1.x and above.
  Furthermore, -I[LIBXC_DIR]/INCLUDE needs to be added to FCFLAGS.
#else
  USE libxc_funcs_m
  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_f90_lib_m,                    ONLY: xc_f90_func_end,&
                                             xc_f90_func_init,&
                                             xc_f90_functional_get_number,&
                                             xc_f90_family_from_id,&
!
                                             xc_f90_gga_exc,&
                                             xc_f90_gga_exc_vxc,&
                                             xc_f90_gga_fxc,&
                                             xc_f90_gga_lb_set_par,&
                                             xc_f90_gga_vxc,&
                                             xc_f90_gga_x_hjs_set_par,&
                                             xc_f90_gga_x_wpbeh_set_par,&
!
                                             xc_f90_hyb_gga_xc_hse_set_par,&
!
                                             xc_f90_info_family,&
                                             xc_f90_info_kind,&
                                             xc_f90_info_name,&
                                             xc_f90_info_refs,&
                                             xc_f90_info_flags,&
!
                                             xc_f90_lda,&
                                             xc_f90_lda_c_1d_csc_set_par,&
                                             xc_f90_lda_c_2d_prm_set_par,&
                                             xc_f90_lda_c_xalpha_set_par,&
                                             xc_f90_lda_exc,&
                                             xc_f90_lda_exc_vxc,&
                                             xc_f90_lda_fxc,&
                                             xc_f90_lda_kxc,&
                                             xc_f90_lda_vxc,&
                                             xc_f90_lda_x_1d_set_par,&
                                             xc_f90_lda_x_set_par,&
!
                                             xc_f90_mgga,&
                                             xc_f90_mgga_exc,&
                                             xc_f90_mgga_exc_vxc,&
                                             xc_f90_mgga_fxc,&
                                             xc_f90_mgga_vxc,&
                                             xc_f90_mgga_x_tb09_set_par,&
!
                                             xc_f90_pointer_t,&
!
                                             XC_FAMILY_UNKNOWN,&
                                             XC_FAMILY_NONE,&
                                             XC_FAMILY_LDA,&
                                             XC_FAMILY_GGA,&
                                             XC_FAMILY_MGGA,&
                                             XC_FAMILY_LCA,&
                                             XC_FAMILY_OEP,&
                                             XC_FAMILY_HYB_GGA,&
                                             XC_FAMILY_HYB_MGGA,&
!
                                             XC_UNPOLARIZED,&
                                             XC_POLARIZED,&
                                             XC_NON_RELATIVISTIC,&
                                             XC_RELATIVISTIC,&
!
                                             XC_EXCHANGE,&
                                             XC_CORRELATION,&
                                             XC_EXCHANGE_CORRELATION,&
                                             XC_KINETIC,&
!
                                             XC_FLAGS_HAVE_EXC,&
                                             XC_FLAGS_HAVE_VXC,&
                                             XC_FLAGS_HAVE_FXC,&
                                             XC_FLAGS_HAVE_KXC,&
                                             XC_FLAGS_HAVE_LXC,&
                                             XC_FLAGS_1D,&
                                             XC_FLAGS_2D,&
                                             XC_FLAGS_3D,&
                                             XC_FLAGS_STABLE,&
                                             XC_FLAGS_DEVELOPMENT,&
                                             XC_GGA_XC_LB,&
                                             XC_GGA_K_ABSR1

#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_libxc_wrap'

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: libxc_version = XC_VERSION

  PUBLIC :: xc_f90_pointer_t
  PUBLIC :: xc_f90_func_init, xc_f90_func_end
  PUBLIC :: xc_f90_info_family, xc_f90_info_kind, xc_f90_info_name
  PUBLIC :: xc_f90_gga_exc, xc_f90_gga_exc_vxc, xc_f90_gga_fxc,&
            xc_f90_gga_lb_set_par, xc_f90_gga_vxc, xc_f90_gga_x_hjs_set_par,&
            xc_f90_gga_x_wpbeh_set_par, xc_f90_hyb_gga_xc_hse_set_par
  PUBLIC :: xc_f90_lda, xc_f90_lda_c_1d_csc_set_par, xc_f90_lda_c_2d_prm_set_par,&
            xc_f90_lda_c_xalpha_set_par, xc_f90_lda_exc, xc_f90_lda_exc_vxc,&
            xc_f90_lda_fxc, xc_f90_lda_kxc, xc_f90_lda_vxc,&
            xc_f90_lda_x_1d_set_par, xc_f90_lda_x_set_par
  PUBLIC :: xc_f90_mgga, xc_f90_mgga_exc, xc_f90_mgga_exc_vxc, xc_f90_mgga_fxc,&
            xc_f90_mgga_vxc, xc_f90_mgga_x_tb09_set_par

  PUBLIC :: XC_FAMILY_UNKNOWN, XC_FAMILY_NONE,&
            XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_LCA,&
            XC_FAMILY_OEP, XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA

  PUBLIC :: XC_UNPOLARIZED, XC_POLARIZED, XC_NON_RELATIVISTIC, XC_RELATIVISTIC

  PUBLIC :: XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC

  PUBLIC :: XC_FLAGS_HAVE_EXC, XC_FLAGS_HAVE_VXC, XC_FLAGS_HAVE_FXC,&
            XC_FLAGS_HAVE_KXC, XC_FLAGS_HAVE_LXC,&
            XC_FLAGS_1D, XC_FLAGS_2D, XC_FLAGS_3D,&
            XC_FLAGS_STABLE, XC_FLAGS_DEVELOPMENT,&
            XC_GGA_XC_LB, XC_GGA_K_ABSR1

! wrappers for routines where interface has changed between versions
  PUBLIC :: xc_libxc_wrap_info_refs,&
            xc_libxc_wrap_version,&
            xc_libxc_wrap_functional_get_number,&
            xc_libxc_wrap_needs_laplace,&
            xc_libxc_wrap_functional_set_params,&
            xc_libxc_wrap_functional_buggy

CONTAINS

! *****************************************************************************
!> \brief Provides the reference(s) for this functional.
!> \param xc_info ...
!> \param polarized ...
!> \param sc ...
!> \param reference ...
!>
!> \author A. Gloessa (agloess)
! *****************************************************************************
  SUBROUTINE xc_libxc_wrap_info_refs(xc_info, polarized, sc, reference)
    TYPE(xc_f90_pointer_t), INTENT(IN) :: xc_info
    INTEGER, INTENT(IN)                :: polarized
    REAL(KIND=dp), INTENT(IN)          :: sc
    CHARACTER(LEN=*), INTENT(OUT)      :: reference

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_info_refs', &
      routineP = moduleN//':'//routineN

    INTEGER, PARAMETER                 :: maxlen = 67 ! max line length in output
    INTEGER                            :: handle, i_ref, first, last, empty
    CHARACTER(LEN=200)                 :: ref_string ! string for one reference
    CHARACTER(LEN=200)                 :: tmp_string ! modified string for one reference
#if (XC_MAJOR_VERSION == 2)
    TYPE(xc_f90_pointer_t)             :: str    ! this will hold a (char **) pointer
#endif

    CALL timeset(routineN, handle)

    i_ref=0
#if (XC_MAJOR_VERSION == 2)
    CALL xc_f90_info_refs(xc_info, i_ref, str, ref_string)
#else
    CALL xc_f90_info_refs(xc_info, i_ref, ref_string)
#endif
    first = 1
    DO WHILE (i_ref >= 0)
       WRITE(tmp_string, '(a1,i1,a2,a)') '[', i_ref, '] ', TRIM(ref_string)
       last = first + LEN_TRIM(tmp_string) - 1
       reference(first:last) = TRIM(tmp_string)
       first = last + 1
       empty = last - MOD(last, maxlen) + maxlen
       IF (empty /= last) THEN
          ! fill with 'spaces'
          reference(first:empty) = ' '
          first = empty + 1
       END IF
#if (XC_MAJOR_VERSION == 2)
       CALL xc_f90_info_refs(xc_info, i_ref, str, ref_string)
#else
       CALL xc_f90_info_refs(xc_info, i_ref, ref_string)
#endif
    END DO
    SELECT CASE (polarized)
     CASE (XC_UNPOLARIZED)
        WRITE(tmp_string, "('{scale=',f5.3,', spin-unpolarized}')") sc
     CASE (XC_POLARIZED)
        WRITE(tmp_string, "('{scale=',f5.3,', spin-polarized}')") sc
     CASE default
        CALL stop_program(routineN,moduleN,__LINE__,&
           "Unsupported value for variable 'polarized'.")
    END SELECT
    last = first + LEN_TRIM(tmp_string) - 1
    reference(first:last) = TRIM(tmp_string)
    first = last + 1
    ! fill with 'spaces'
    reference(first:LEN(reference)) = ' '

    CALL timestop(handle)

  END SUBROUTINE xc_libxc_wrap_info_refs

! *****************************************************************************
!> \brief Provides a version string.
!> \param version ...
!> \author A. Gloessa (agloess)
!> \note Minor and micro version could be defined as character for SVN trunk
!>       version (e.g. 3.x.x)!
!>
! *****************************************************************************
  SUBROUTINE xc_libxc_wrap_version(version)
    CHARACTER(LEN=*), INTENT(OUT)            :: version

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_version', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

! the string that is output
!#if (XC_MAJOR_VERSION >= 3)
!    INTEGER                            :: vmajor, vminor, vmicro
!#endif

    CALL timeset(routineN, handle)

!#if (XC_MAJOR_VERSION >= 3)
!    CALL xc_f90_version(vmajor, vminor, vmicro)
!    write(version(1:5),'(i1,a,i1,a,i1)') vmajor, '.', vminor, '.', vmicro
!#else
!    CALL xc_f90_version(vmajor, vminor)
!    write(version(1:5),'(i1,a,i1,a)') vmajor, '.', vminor, '.?'
!#endif
    version = TRIM(libxc_version)

    CALL timestop(handle)

  END SUBROUTINE xc_libxc_wrap_version

! *****************************************************************************
!> \brief Provides the functional ID.
!> \param func_string ...
!> \retval func_id ...
!> \author A. Gloessa (agloess)
!> \note Remove prefix to keep compatibility, functionals can be specified (in
!>       LIBXC section) as:
!>       GGA_X_...  or  XC_GGA_X_...
!>       Starting from version 2.2.0 both name conventions are allowed, before
!>       the 'XC_' prefix was necessary.
!>
! *****************************************************************************
  INTEGER FUNCTION xc_libxc_wrap_functional_get_number(func_string) RESULT(func_id)
    CHARACTER(LEN=*), INTENT(IN)       :: func_string

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_version', &
      routineP = moduleN//':'//routineN

    INTEGER                            :: handle

    CALL timeset(routineN, handle)

    IF (func_string(1:3) == "XC_") THEN
      func_id = xc_f90_functional_get_number(func_string(4:LEN_TRIM(func_string)))
    ELSE
      func_id = xc_f90_functional_get_number(func_string(1:LEN_TRIM(func_string)))
    END IF
    IF (func_id == -1) THEN
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_string)//": wrong functional name")
    END IF

    CALL timestop(handle)

  END FUNCTION xc_libxc_wrap_functional_get_number

! *****************************************************************************
!> \brief Wrapper for functionals that need the Laplacian, all others can use
!>        a dummy array.
!> \param func_id ...
!>
!> \retval xc_libxc_wrap_needs_laplace ...
!> \author A. Gloessa (agloess)
! *****************************************************************************
  LOGICAL FUNCTION xc_libxc_wrap_needs_laplace(func_id)
    ! Only some MGGA functionals need 'r->u' and 'r->us', which are calculated
    ! by default in 'src/work_mgga_x.c' and 'src/work_mgga_c.c' as 'r.u' and
    ! 'r.us'.
    ! All other MGGA functionals one could use a dummy (of the right size) for
    ! 'lapl'.
    ! 
    ! It would be best, if LibXC would provide the information:
    !   need%lapl
    !   need%tau
    !   nedd%...
    ! as FLAGS.
    ! As a work-around we call the corresponding routine twice with different
    ! values for 'lapl' and check if the results are changing.
    ! 
    ! [last checked version: 3.x.x (18.08.2015)]
    ! 
    ! 'r->u':
    ! -------
    ! # src/mgga_x_br89.c
    ! XC_MGGA_X_BR89              206 /* Becke-Roussel 89  */
    ! XC_MGGA_X_BJ06              207 /* Becke & Johnson correction to Becke-Roussel 89  */
    ! XC_MGGA_X_TB09              208 /* Tran & Blaha correction to Becke & Johnson  */
    ! XC_MGGA_X_RPP09             209 /* Rasanen, Pittalis, and Proetto correction to Becke & Johnson  */
    ! 
    ! # src/mgga_x_2d_prhg07.c
    ! XC_MGGA_X_2D_PRHG07         210   /* Pittalis, Rasanen, Helbig, Gross Exchange Functional */
    ! XC_MGGA_X_2D_PRHG07_PRP10   211   /* PRGH07 with PRP10 correction */
    ! 
    ! # src/mgga_x_mk00.c
    ! XC_MGGA_X_MK00              230 /* Exchange for accurate virtual orbital energies */
    ! XC_MGGA_X_MK00B             243 /* Exchange for accurate virtual orbital energies (v. B) */
    ! 
    ! 
    ! 'r->us':
    ! -------
    ! # src/mgga_xc_zlp.c
    ! XC_MGGA_XC_ZLP               42 /* Zhao, Levy & Parr, Eq. (21) */
    ! 
    ! # src/mgga_c_cs.c
    ! XC_MGGA_C_CS                 72 /* Colle and Salvetti */
    ! 
    ! # src/mgga_c_cc06.c
    ! XC_MGGA_C_CC06              229 /* Cancio and Chou 2006 */
    !
    INTEGER, INTENT(IN)                        :: func_id

    CHARACTER(LEN=*), PARAMETER :: routineN = 'xc_libxc_wrap_needs_laplace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    REAL(KIND=dp)                            :: rho, norm_drho, laplace_rho, &
                                                my_tau
    REAL(KIND=dp), DIMENSION(:), &
      ALLOCATABLE                            :: dummy
    REAL(KIND=dp), DIMENSION(:,:), &
      ALLOCATABLE                            :: val
    TYPE(xc_f90_pointer_t)                   :: xc_func, xc_info

    CALL timeset(routineN, handle)

    xc_libxc_wrap_needs_laplace = .FALSE.

    ! Some MGGa need the laplace explicit and some just need an arbitrary array
    ! of the correct size. Here we call the calculation routine twice for two
    ! different values of 'laplace_rho' and check if the result has changed.
    !
    ! Assumption (.true. in v2.1.0 - v3.x.x): 
    !             if 
    !                functional is Laplace-dependent for XC_UNPOLARIZED
    !             then
    !                functional will be Laplace-dependent for XC_POLARIZED too.
    !
!$OMP CRITICAL(libxc_init)
    CALL xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
!$OMP END CRITICAL(libxc_init)
!$OMP BARRIER
    IF (xc_f90_info_family(xc_info) == XC_FAMILY_MGGA) THEN
       IF (IAND(xc_f90_info_flags(xc_info), XC_FLAGS_HAVE_EXC) == XC_FLAGS_HAVE_EXC) THEN
          ALLOCATE(val(1,2))
          rho = 2.0e-3_dp; norm_drho = 2.0e-3_dp; my_tau = 5.0e-1_dp
          laplace_rho = -1.0e-1_dp
          CALL xc_f90_mgga_exc(xc_func, 1, rho, norm_drho, &
                               laplace_rho, my_tau, val(1,1))
          laplace_rho = -3.0e-1_dp
          CALL xc_f90_mgga_exc(xc_func, 1, rho, norm_drho, &
                               laplace_rho, my_tau, val(1,2))
          IF (val(1,1) /= val(1,2)) xc_libxc_wrap_needs_laplace = .TRUE.
          DEALLOCATE(val)
       END IF 
       IF (IAND(xc_f90_info_flags(xc_info), XC_FLAGS_HAVE_VXC) == XC_FLAGS_HAVE_VXC) THEN
          ALLOCATE(val(1,2), dummy(3))
          rho = 2.0e-3_dp; norm_drho = 2.0e-3_dp; my_tau = 5.0e-1_dp
          laplace_rho = -1.0e-1_dp
          CALL xc_f90_mgga_vxc(xc_func, 1, rho, norm_drho,&
                     laplace_rho, my_tau, dummy(1), dummy(2), val(1,1), dummy(3))
          laplace_rho = -3.0e-1_dp
          CALL xc_f90_mgga_vxc(xc_func, 1, rho, norm_drho,&
                     laplace_rho, my_tau, dummy(1), dummy(2), val(1,2), dummy(3))
          IF (val(1,1) /= val(1,2)) xc_libxc_wrap_needs_laplace = .TRUE.
          DEALLOCATE(val, dummy)
       END IF
       IF (IAND(xc_f90_info_flags(xc_info), XC_FLAGS_HAVE_FXC) == XC_FLAGS_HAVE_FXC) THEN
          ALLOCATE(val(4,2),dummy(6))
          rho = 2.0e-3_dp; norm_drho = 2.0e-3_dp; my_tau = 5.0e-1_dp
          laplace_rho = -1.0e-1_dp
          CALL xc_f90_mgga_fxc(xc_func, 1, rho, norm_drho,&
                     laplace_rho, my_tau,&
                     dummy(1), dummy(2), val(1,1), dummy(3), dummy(4), val(2,1),&
                     dummy(5), val(3,1), dummy(6), val(4,1))
          laplace_rho = -3.0e-1_dp
          CALL xc_f90_mgga_fxc(xc_func, 1, rho, norm_drho,&
                     laplace_rho, my_tau,&
                     dummy(1), dummy(2), val(1,2), dummy(3), dummy(4), val(2,2),&
                     dummy(5), val(3,2), dummy(6), val(4,2))
          IF (val(1,1) /= val(1,2) .OR. val(2,1) /= val(2,2) .OR. &
              val(3,1) /= val(3,2) .OR. val(4,1) /= val(4,2)) &
             xc_libxc_wrap_needs_laplace = .TRUE.
          DEALLOCATE(val, dummy)
       END IF
    END IF
    CALL xc_f90_func_end(xc_func)

!    IF (func_id == XC_MGGA_X_BR89)            xc_libxc_wrap_needs_laplace = .TRUE.
!    IF (func_id == XC_MGGA_X_BJ06)            xc_libxc_wrap_needs_laplace = .TRUE.
!    IF (func_id == XC_MGGA_X_TB09)            xc_libxc_wrap_needs_laplace = .TRUE.
!    IF (func_id == XC_MGGA_X_RPP09)           xc_libxc_wrap_needs_laplace = .TRUE.
!    IF (func_id == XC_MGGA_X_2D_PRHG07)       xc_libxc_wrap_needs_laplace = .TRUE.
!    IF (func_id == XC_MGGA_X_2D_PRHG07_PRP10) xc_libxc_wrap_needs_laplace = .TRUE.
!    IF (func_id == XC_MGGA_X_MK00)            xc_libxc_wrap_needs_laplace = .TRUE.
!    IF (func_id == XC_MGGA_X_MK00B)           xc_libxc_wrap_needs_laplace = .TRUE.
!#if (XC_MAJOR_VERSION == 3)
!    IF (func_id == XC_MGGA_XC_ZLP)            xc_libxc_wrap_needs_laplace = .TRUE.
!#endif
!    IF (func_id == XC_MGGA_C_CS)              xc_libxc_wrap_needs_laplace = .TRUE.
!    IF (func_id == XC_MGGA_C_CC06)            xc_libxc_wrap_needs_laplace = .TRUE.

    CALL timestop(handle)

  END FUNCTION xc_libxc_wrap_needs_laplace

! *****************************************************************************
!> \brief Wrapper for functionals that need special parameters.
!> \param xc_func ...
!> \param xc_info ...
!> \param func_id ...
!> \param params ...
!> \param no_exc ...
!>
!> \author A. Gloessa (agloess)
! *****************************************************************************
  SUBROUTINE xc_libxc_wrap_functional_set_params(xc_func, xc_info, func_id, params, no_exc)
    TYPE(xc_f90_pointer_t), INTENT(INOUT)    :: xc_func
    TYPE(xc_f90_pointer_t), INTENT(IN)       :: xc_info
    INTEGER, INTENT(IN)                      :: func_id
    REAL(KIND=dp), DIMENSION(*), INTENT(IN)  :: params
    LOGICAL, INTENT(INOUT)                   :: no_exc

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'xc_libxc_wrap_functional_set_params', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN, handle)

!ToDo: This is dangerous, since wrong parameter ordering or missing ones
!      result in wrong numbers or a crash. 
!      --> LibXC needs to provide more informations via a function.
    IF (params(1) < 0.99e20_dp) THEN
       SELECT CASE (func_id)
        CASE (XC_LDA_X)
           CALL xc_f90_lda_x_set_par(xc_func,params(1),NINT(params(2)),params(3))
        CASE (XC_LDA_C_XALPHA)
           CALL xc_f90_lda_c_xalpha_set_par(xc_func,params(1))
        CASE (XC_LDA_C_2D_PRM)
           CALL xc_f90_lda_c_2d_prm_set_par(xc_func,params(1))
        CASE (XC_LDA_C_1D_CSC)
           CALL xc_f90_lda_c_1d_csc_set_par(xc_func,NINT(params(1)),params(2))
        CASE (XC_LDA_X_1D)
           CALL xc_f90_lda_x_1d_set_par(xc_func,NINT(params(1)),params(2))
        CASE (XC_GGA_X_LB)
           CALL xc_f90_gga_lb_set_par(xc_func,NINT(params(1)),params(2),params(3),params(4))
        CASE (XC_MGGA_X_TB09)
           CALL xc_f90_mgga_x_tb09_set_par(xc_func,params(1))
        CASE (XC_HYB_GGA_XC_HSE03, XC_HYB_GGA_XC_HSE06)
           CALL xc_f90_hyb_gga_xc_hse_set_par(xc_func,params(1),params(2))
        CASE (XC_HYB_GGA_XC_HJS_PBE, XC_HYB_GGA_XC_HJS_PBE_SOL, &
              XC_HYB_GGA_XC_HJS_B88, XC_HYB_GGA_XC_HJS_B97X, &
              XC_GGA_X_HJS_PBE, XC_GGA_X_HJS_PBE_SOL, &
              XC_GGA_X_HJS_B88, XC_GGA_X_HJS_B97X)
           CALL xc_f90_gga_x_hjs_set_par(xc_func,params(1))
        CASE (XC_GGA_X_WPBEH)
           CALL xc_f90_gga_x_wpbeh_set_par(xc_func,params(1))
       END SELECT
    END IF

    IF (IAND(xc_f90_info_flags(xc_info), XC_FLAGS_HAVE_EXC) == XC_FLAGS_HAVE_EXC) THEN
       no_exc = .FALSE.
    ELSE
       no_exc = .TRUE.
    END IF

    CALL timestop(handle)

  END SUBROUTINE xc_libxc_wrap_functional_set_params

! *****************************************************************************
!> \brief Wrapper for known buggy functionals.
!> \param func_id ...
!> \param grad_deriv ...
!> \author A. Gloessa (agloess)
! *****************************************************************************
  SUBROUTINE xc_libxc_wrap_functional_buggy(func_id, grad_deriv)
!
! BugFix information was taken from:
!
! http://www.tddft.org/programs/octopus/wiki/index.php/Libxc_changes
!
!
! corrected with: 2.1.1 and 2.2.1
!    XC_GGA_X_SOGGA11, XC_HYB_GGA_X_SOGGA11_X, XC_HYB_GGA_XC_CAM_B3LYP,
!    XC_HYB_GGA_XC_TUNED_CAM_B3LYP
!
! corrected with: 2.1.2 and 2.2.2
!    XC_GGA_X_LB, XC_GGA_X_LBM, XC_MGGA_C_REVTPSS, XC_MGGA_C_M08_HX,
!    XC_MGGA_C_M08_SO, XC_MGGA_C_M11, XC_MGGA_C_M11_L 
!  + all XC_HYB_GGA and XC_HYB_MGGA_*
!
! corrected with: 2.1.3 and 2.2.3
!    XC_GGA_X_N12, XC_GGA_C_N12, XC_GGA_C_N12_SX, XC_HYB_GGA_X_N12_SX,
!    XC_HYB_GGA_XC_O3LYP, XC_HYB_GGA_XC_X3LYP, XC_GGA_X_B88,
!    XC_GGA_X_OPTB88_VDW, XC_GGA_X_MB88, XC_GGA_K_LLP, XC_GGA_K_FR_B88,
!    XC_GGA_K_THAKKAR
! + all XC_GGA_X_HJS_*
! + all XC_HYB_GGA_XC_HJS_*
!    second order derivatives (XC_GGA_X_B88, XC_GGA_X_OPTB88_VDW,
!     XC_GGA_X_MB88, XC_GGA_K_LLP, XC_GGA_K_FR_B88, XC_GGA_K_THAKKAR)
!
! corrected with: 3.x.x
!    XC_MGGA_X_M11, XC_MGGA_X_M11_L, XC_HYB_MGGA_X_MS2H
!
!
! Note, some variables were redefined between version 2.1.x and 3.x.x!
!
! XC_HYB_MGGA_X_M11 --> XC_MGGA_X_M11
! XC_MGGA_X_MS2H    --> XC_HYB_MGGA_X_MS2H
!
    INTEGER, INTENT(IN)                      :: func_id
    INTEGER, INTENT(IN), OPTIONAL            :: grad_deriv

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'xc_libxc_wrap_functional_buggy', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=256)                       :: func_name
    INTEGER                                  :: handle, i, my_grad_deriv, &
                                                nbuggy
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: buggy
    LOGICAL                                  :: bug

    CALL timeset(routineN, handle)

    bug = .FALSE.
    my_grad_deriv = 0
    IF (PRESENT(grad_deriv)) my_grad_deriv = ABS(grad_deriv)
    SELECT CASE (libxc_version)
#if ( XC_MAJOR_VERSION == 2)
     CASE ("2.1.0", "2.2.0")
        IF (my_grad_deriv == 2) THEN
           nbuggy = 34; ALLOCATE(buggy(nbuggy))
           buggy(:) = (/ XC_GGA_X_SOGGA11, XC_HYB_GGA_X_SOGGA11_X, &
            XC_HYB_GGA_XC_CAM_B3LYP, XC_HYB_GGA_XC_TUNED_CAM_B3LYP, XC_GGA_X_LB, &
            XC_GGA_X_LBM, XC_MGGA_C_REVTPSS, XC_MGGA_C_M08_HX, XC_MGGA_C_M08_SO, &
            XC_MGGA_C_M11, XC_MGGA_C_M11_L, XC_GGA_X_N12, XC_GGA_C_N12, &
            XC_GGA_C_N12_SX,XC_HYB_GGA_X_N12_SX, XC_HYB_GGA_XC_O3LYP, &
            XC_HYB_GGA_XC_X3LYP, XC_GGA_X_B88, XC_GGA_X_OPTB88_VDW, &
            XC_GGA_X_MB88, XC_GGA_K_LLP, XC_GGA_K_FR_B88, XC_GGA_K_THAKKAR, &
            XC_HYB_MGGA_X_M11, XC_MGGA_X_M11_L, XC_MGGA_X_MS2H, &
            XC_GGA_X_HJS_PBE, XC_GGA_X_HJS_PBE_SOL, &
            XC_GGA_X_HJS_B88, XC_GGA_X_HJS_B97X, XC_HYB_GGA_XC_HJS_PBE, &
            XC_HYB_GGA_XC_HJS_PBE_SOL, XC_HYB_GGA_XC_HJS_B88, XC_HYB_GGA_XC_HJS_B97X /)
        ELSE
           nbuggy = 28; ALLOCATE(buggy(nbuggy))
           buggy(:) = (/ XC_GGA_X_SOGGA11, XC_HYB_GGA_X_SOGGA11_X, &
            XC_HYB_GGA_XC_CAM_B3LYP, XC_HYB_GGA_XC_TUNED_CAM_B3LYP, XC_GGA_X_LB, &
            XC_GGA_X_LBM, XC_MGGA_C_REVTPSS, XC_MGGA_C_M08_HX, XC_MGGA_C_M08_SO, &
            XC_MGGA_C_M11, XC_MGGA_C_M11_L, XC_GGA_X_N12, XC_GGA_C_N12, &
            XC_GGA_C_N12_SX,XC_HYB_GGA_X_N12_SX, XC_HYB_GGA_XC_O3LYP, &
            XC_HYB_GGA_XC_X3LYP, &
            XC_HYB_MGGA_X_M11, XC_MGGA_X_M11_L, XC_MGGA_X_MS2H, &
            XC_GGA_X_HJS_PBE, XC_GGA_X_HJS_PBE_SOL, &
            XC_GGA_X_HJS_B88, XC_GGA_X_HJS_B97X, XC_HYB_GGA_XC_HJS_PBE, &
            XC_HYB_GGA_XC_HJS_PBE_SOL, XC_HYB_GGA_XC_HJS_B88, XC_HYB_GGA_XC_HJS_B97X /)
        END IF
        IF ((xc_f90_family_from_id(func_id) == XC_FAMILY_HYB_GGA) .OR. &
            (xc_f90_family_from_id(func_id) == XC_FAMILY_HYB_MGGA)) THEN
           bug = .TRUE.
        END IF
     CASE ("2.1.1", "2.2.1")
        IF (my_grad_deriv == 2) THEN
           nbuggy = 30; ALLOCATE(buggy(nbuggy))
           buggy(:) = (/ XC_GGA_X_LB, &
            XC_GGA_X_LBM, XC_MGGA_C_REVTPSS, XC_MGGA_C_M08_HX, XC_MGGA_C_M08_SO, &
            XC_MGGA_C_M11, XC_MGGA_C_M11_L, XC_GGA_X_N12, XC_GGA_C_N12, &
            XC_GGA_C_N12_SX,XC_HYB_GGA_X_N12_SX, XC_HYB_GGA_XC_O3LYP, &
            XC_HYB_GGA_XC_X3LYP, XC_GGA_X_B88, XC_GGA_X_OPTB88_VDW, &
            XC_GGA_X_MB88, XC_GGA_K_LLP, XC_GGA_K_FR_B88, XC_GGA_K_THAKKAR, &
            XC_HYB_MGGA_X_M11, XC_MGGA_X_M11_L, XC_MGGA_X_MS2H, &
            XC_GGA_X_HJS_PBE, XC_GGA_X_HJS_PBE_SOL, &
            XC_GGA_X_HJS_B88, XC_GGA_X_HJS_B97X, XC_HYB_GGA_XC_HJS_PBE, &
            XC_HYB_GGA_XC_HJS_PBE_SOL, XC_HYB_GGA_XC_HJS_B88, XC_HYB_GGA_XC_HJS_B97X /)
        ELSE
           nbuggy = 24; ALLOCATE(buggy(nbuggy))
           buggy(:) = (/ XC_GGA_X_LB, &
            XC_GGA_X_LBM, XC_MGGA_C_REVTPSS, XC_MGGA_C_M08_HX, XC_MGGA_C_M08_SO, &
            XC_MGGA_C_M11, XC_MGGA_C_M11_L, XC_GGA_X_N12, XC_GGA_C_N12, &
            XC_GGA_C_N12_SX,XC_HYB_GGA_X_N12_SX, XC_HYB_GGA_XC_O3LYP, &
            XC_HYB_GGA_XC_X3LYP, &
            XC_HYB_MGGA_X_M11, XC_MGGA_X_M11_L, XC_MGGA_X_MS2H, &
            XC_GGA_X_HJS_PBE, XC_GGA_X_HJS_PBE_SOL, &
            XC_GGA_X_HJS_B88, XC_GGA_X_HJS_B97X, XC_HYB_GGA_XC_HJS_PBE, &
            XC_HYB_GGA_XC_HJS_PBE_SOL, XC_HYB_GGA_XC_HJS_B88, XC_HYB_GGA_XC_HJS_B97X /)
        END IF
        IF ((xc_f90_family_from_id(func_id) == XC_FAMILY_HYB_GGA) .OR. &
            (xc_f90_family_from_id(func_id) == XC_FAMILY_HYB_MGGA)) THEN
           bug = .TRUE.
        END IF
     CASE ("2.1.2", "2.2.2")
        IF (my_grad_deriv == 2) THEN
           nbuggy = 23; ALLOCATE(buggy(nbuggy))
           buggy(:) = (/ XC_GGA_X_N12, XC_GGA_C_N12, &
            XC_GGA_C_N12_SX,XC_HYB_GGA_X_N12_SX, XC_HYB_GGA_XC_O3LYP, &
            XC_HYB_GGA_XC_X3LYP, XC_GGA_X_B88, XC_GGA_X_OPTB88_VDW, &
            XC_GGA_X_MB88, XC_GGA_K_LLP, XC_GGA_K_FR_B88, XC_GGA_K_THAKKAR, &
            XC_HYB_MGGA_X_M11, XC_MGGA_X_M11_L, XC_MGGA_X_MS2H, &
            XC_GGA_X_HJS_PBE, XC_GGA_X_HJS_PBE_SOL, &
            XC_GGA_X_HJS_B88, XC_GGA_X_HJS_B97X, XC_HYB_GGA_XC_HJS_PBE, &
            XC_HYB_GGA_XC_HJS_PBE_SOL, XC_HYB_GGA_XC_HJS_B88, XC_HYB_GGA_XC_HJS_B97X /)
        ELSE
           nbuggy = 17; ALLOCATE(buggy(nbuggy))
           buggy(:) = (/ XC_GGA_X_N12, XC_GGA_C_N12, &
            XC_GGA_C_N12_SX,XC_HYB_GGA_X_N12_SX, XC_HYB_GGA_XC_O3LYP, &
            XC_HYB_GGA_XC_X3LYP, &
            XC_HYB_MGGA_X_M11, XC_MGGA_X_M11_L, XC_MGGA_X_MS2H, &
            XC_GGA_X_HJS_PBE, XC_GGA_X_HJS_PBE_SOL, &
            XC_GGA_X_HJS_B88, XC_GGA_X_HJS_B97X, XC_HYB_GGA_XC_HJS_PBE, &
            XC_HYB_GGA_XC_HJS_PBE_SOL, XC_HYB_GGA_XC_HJS_B88, XC_HYB_GGA_XC_HJS_B97X /)
        END IF
     CASE ("2.1.3", "2.2.3")
        nbuggy = 3; ALLOCATE(buggy(nbuggy))
        buggy(:) = (/ XC_HYB_MGGA_X_M11, XC_MGGA_X_M11_L, XC_MGGA_X_MS2H /)
#endif
     CASE ("3.x.x")
        nbuggy = 0
     CASE default
        nbuggy = 0
           CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                " This version ("//TRIM(libxc_version)//") of LibXC is new or unknown."//&
                " ======== Please check the results carefully. ========    "//&
                " More informations on bugfixes can be found at: "//&
                " http://www.tddft.org/programs/octopus/wiki/index.php/Libxc_changes"//&
                " Calculation continues using unsupported LibXC version. "//&
                CPSourceFileRef,&
                only_ionode=.TRUE.)
    END SELECT

    DO i=1,nbuggy
       IF (func_id == buggy(i)) THEN
          bug = .TRUE.
       END IF
    END DO

    IF (bug) THEN
       CALL xc_f90_functional_get_name(func_id, func_name)
       CALL stop_program(routineN,moduleN,__LINE__,&
          TRIM(func_name)//": deactivated since buggy in version "//TRIM(libxc_version)//" of LibXC.")
    END IF

    IF (ALLOCATED(buggy)) DEALLOCATE(buggy)
    
    CALL timestop(handle)

  END SUBROUTINE xc_libxc_wrap_functional_buggy
#endif
#endif
END MODULE xc_libxc_wrap
