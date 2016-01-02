!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2016  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE eri_mme_test

  USE cell_types,                      ONLY: cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             init_cell
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE eri_mme_methods,                 ONLY: eri_mme_set_params
  USE eri_mme_methods_low,             ONLY: eri_mme_2c_integrate_low
  USE eri_mme_types,                   ONLY: create_eri_mme_section,&
                                             dir_lattice_space,&
                                             eri_mme_init_read_input,&
                                             eri_mme_param,&
                                             eri_mme_release,&
                                             lmax,&
                                             rec_lattice_space
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE input_val_types,                 ONLY: real_t
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_sum
  USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                             ncoset
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eri_mme_test'

  PUBLIC :: create_eri_mme_test_section, eri_mme_perf_acc_test

  CONTAINS

! *****************************************************************************
!> \brief Create input section for unit testing
!> \param section ...
! *****************************************************************************
  SUBROUTINE create_eri_mme_test_section(section)
    TYPE(section_type), INTENT(INOUT), &
      POINTER                                :: section

    CHARACTER(len=*), PARAMETER :: routineN = 'create_eri_mme_test_section', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    NULLIFY(keyword, subsection)

    CPASSERT(.NOT.ASSOCIATED(section))
    CALL section_create(section, "ERI_MME_TEST", &
                        description = "Parameters for testing the ERI_MME method for electron repulsion integrals. "//&
                                      "Testing w.r.t. performance and accuracy. ", &
                        n_keywords = 5, n_subsections = 1)


    CALL create_eri_mme_section(subsection)
    CALL section_add_subsection(section, subsection)
    CALL section_release(subsection)

    CALL keyword_create(keyword, &
                        name = "_SECTION_PARAMETERS_", &
                        description = "Controls the activation the ERI_MME test. ", &
                        default_l_val = .FALSE., &
                        lone_keyword_l_val = .TRUE.)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name = "ABC", &
                        description = "Specify the lengths of the cell vectors A, B, and C. ", &
                        usage = "ABC 10.000 10.000 10.000", unit_str = "angstrom", &
                        n_var = 3, type_of_var = real_t)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name = "NAB_MIN", &
                        description = "Minimum number of atomic distances to consider. ", &
                        default_i_val = 8)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name = "NREP", &
                        description = "Number of repeated calculation of each integral. ", &
                        default_i_val = 1)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name = "CHECK_ACCURACY", &
                        description = "Whether integrals should be compared to reference values, "//&
                                      "created on the fly by exact method (G-space sum on grid without "//&
                                      "minimax approximation). Note: only feasible for not too many "//&
                                      "integrals and maximum exponent around 10.0. ", &
                        default_l_val = .TRUE.)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name = "LMAX", &
                        description = "Maximum total angular momentum quantum number. ", &
                        default_i_val = 6)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name = "ZET_MIN", &
                        description = "Minimum exponent. ", &
                        default_r_val = 0.001_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name = "ZET_MAX", &
                        description = "Maximum exponent. ", &
                        default_r_val = 1.0_dp)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

    CALL keyword_create(keyword, name = "NZET", &
                        description = "Number of exponents (logarithmic partition between ZET_MIN and ZET_MAX). ", &
                        default_i_val = 4)
    CALL section_add_keyword(section, keyword)
    CALL keyword_release(keyword)

 END SUBROUTINE create_eri_mme_test_section

! *****************************************************************************
!> \brief Unit test for performance and accuracy (orthorhombic cell only)
!> \param para_env ...
!> \param iw ...
!> \param eri_mme_test_section ...
! *****************************************************************************
  SUBROUTINE eri_mme_perf_acc_test(para_env, iw, eri_mme_test_section)
    TYPE(cp_para_env_type), INTENT(IN), &
      POINTER                                :: para_env
    INTEGER, INTENT(IN)                      :: iw
    TYPE(section_vals_type), INTENT(INOUT), &
      POINTER                                :: eri_mme_test_section

    CHARACTER(len=*), PARAMETER :: routineN = 'eri_mme_perf_acc_test', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: count_ab, i, iab, irep, izet, &
                                                jab, kab, l, lmax, nab, &
                                                nab_min, nab_xyz, nrep, nzet
    LOGICAL                                  :: acc_check, test_accuracy
    REAL(KIND=dp)                            :: acc_G, acc_R, G_res, t0, t1, &
                                                t2, zet_fac, zetmax, zetmin
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: zet
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: rab, time_G, time_R
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: I_diff_G, I_diff_R, I_ref, &
                                                I_test_G, I_test_R
    REAL(KIND=dp), DIMENSION(3, 3)           :: ht
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cell_par
    TYPE(cell_type), POINTER                 :: box
    TYPE(eri_mme_param)                      :: param
    TYPE(section_vals_type), POINTER         :: eri_mme_section

    NULLIFY(box, eri_mme_section, cell_par)

    eri_mme_section => section_vals_get_subs_vals(eri_mme_test_section, "ERI_MME")
    CALL eri_mme_init_read_input(eri_mme_section, param)
    CALL section_vals_val_get(eri_mme_test_section, "ABC", r_vals = cell_par)
    CALL section_vals_val_get(eri_mme_test_section, "NAB_MIN", i_val = nab_min)
    CALL section_vals_val_get(eri_mme_test_section, "NREP", i_val = nrep)
    CALL section_vals_val_get(eri_mme_test_section, "CHECK_ACCURACY", l_val = test_accuracy)
    CALL section_vals_val_get(eri_mme_test_section, "LMAX", i_val = lmax)
    CALL section_vals_val_get(eri_mme_test_section, "ZET_MIN", r_val = zetmin)
    CALL section_vals_val_get(eri_mme_test_section, "ZET_MAX", r_val = zetmax)
    CALL section_vals_val_get(eri_mme_test_section, "NZET", i_val = nzet)
    IF(nzet .LE. 0)&
       CPABORT("Number of exponents NZET must be greater than 0.")


    CALL init_orbital_pointers(lmax)

    ! Create ranges of zet to be tested
    ALLOCATE(zet(nzet))

    zet(1) = zetmin
    IF(nzet .GT. 1) THEN
      zet_fac = (zetmax/zetmin)**(1.0_dp/(nzet - 1))
      DO i = 1, nzet - 1
        zet(i + 1) = zet(i)*zet_fac
      ENDDO
    ENDIF

    ! Create range of rab (atomic distances) to be tested
    nab_xyz = CEILING(REAL(nab_min, KIND = dp)**(1.0_dp/3.0_dp) - 1.0E-06)
    nab = nab_xyz**3

    ALLOCATE(rab(3, nab))
    count_ab = 0
    DO iab = 1, nab_xyz
      DO jab = 1, nab_xyz
        DO kab = 1, nab_xyz
          count_ab = count_ab + 1
          rab(:, count_ab) = [iab*ABS(cell_par(1)), jab*ABS(cell_par(2)), kab*ABS(cell_par(3))]/nab_xyz
        ENDDO
      ENDDO
    ENDDO

    ! initialize cell
    CALL cell_create(box)
    box%hmat = 0.0_dp
    box%hmat(1, 1) = cell_par(1)
    box%hmat(2, 2) = cell_par(2)
    box%hmat(3, 3) = cell_par(3)
    CALL init_cell ( box )

    ! initialize MME method
    CALL eri_mme_set_params(param, box, MINVAL(zet), [MAXVAL(zet)], [lmax], lmax, para_env)
    IF(iw>0) CALL m_flush(iw)

    ! Calculate reference values (Exact expression in G space converged to high precision)
    IF (test_accuracy) THEN
      ht = twopi*TRANSPOSE(param%h_inv)
      G_res = 0.5_dp*param%G_min

      ALLOCATE(I_ref(ncoset(lmax), ncoset(lmax), nab, nzet))
      I_ref(:, :, :, :) = 0.0_dp

      DO izet = 1, nzet
        DO iab = 1, nab
          CALL eri_mme_2c_integrate_low(param, 0, lmax, 0, lmax, zet(izet), zet(izet), rab(:, iab), &
                                     I_ref(:, :, iab, izet), 0, 0, normalize = .TRUE., exact_method = .TRUE.)

        ENDDO
      ENDDO
    ENDIF

    ! test performance and accuracy of MME method
    ALLOCATE(I_test_G(ncoset(lmax), ncoset(lmax), nab, nzet))
    ALLOCATE(I_test_R(ncoset(lmax), ncoset(lmax), nab, nzet))
    ALLOCATE(I_diff_G(ncoset(lmax), ncoset(lmax), nab, nzet))
    ALLOCATE(I_diff_R(ncoset(lmax), ncoset(lmax), nab, nzet))

    ALLOCATE(time_G(0:lmax, nzet))
    ALLOCATE(time_R(0:lmax, nzet))
    DO l = 0, lmax
      DO izet = 1, nzet
        CALL CPU_TIME(t0)
        DO irep = 1, nrep
          DO iab = 1, nab
            CALL eri_mme_2c_integrate_low(param, 0, l, 0, l, zet(izet), zet(izet), rab(:, iab), &
                                       I_test_G(:, :, iab, izet), 0, 0, normalize = .TRUE., &
                                       lattice_space = rec_lattice_space)
          ENDDO
        ENDDO
        CALL CPU_TIME(t1)
        DO irep = 1, nrep
          DO iab = 1, nab
            CALL eri_mme_2c_integrate_low(param, 0, l, 0, l, zet(izet), zet(izet), rab(:, iab), &
                                       I_test_R(:, :, iab, izet), 0, 0, normalize = .TRUE., &
                                       lattice_space = dir_lattice_space)
          ENDDO
        ENDDO
        CALL CPU_TIME(t2)
        time_G(l, izet) = t1 - t0
        time_R(l, izet) = t2 - t1
      ENDDO
    ENDDO

    CALL mp_sum(time_G, para_env%group)
    CALL mp_sum(time_R, para_env%group)

    IF(test_accuracy) THEN
      I_diff_G(:, :, :, :) = ABS(I_test_G - I_ref)
      I_diff_R(:, :, :, :) = ABS(I_test_R - I_ref)
    ENDIF

    IF(iw>0) THEN
      WRITE(iw, '(T2, A, T61, I20)') "ERI_MME| Number of atomic distances:", nab
      WRITE(iw, '(/T2, A)') "ERI_MME| Test results for time and accuracy in G and R space"
      WRITE(iw, '(T11, A)') "l, zet, cpu time G, cpu time R, accuracy G, accuracy R"

      DO l = 0, lmax
        DO izet = 1, nzet
          IF(test_accuracy) THEN
            acc_G = MAXVAL(I_diff_G(ncoset(l - 1) + 1:ncoset(l), ncoset(l - 1) + 1:ncoset(l), :, izet))
            acc_R = MAXVAL(I_diff_R(ncoset(l - 1) + 1:ncoset(l), ncoset(l - 1) + 1:ncoset(l), :, izet))
          ELSE
            acc_G = 0.0_dp
            acc_R = 0.0_dp
          ENDIF

          WRITE(iw, '(T11, I1, 1X, ES9.2, 1X, ES9.2, 1X, ES9.2, 1X, ES9.2, 1X, ES9.2)') &
                   l, zet(izet), time_G(l, izet), time_R(l, izet), acc_G, acc_R
        ENDDO
      ENDDO

      IF(test_accuracy) THEN
        WRITE(iw, '(/T2, A, 32X, ES9.2, 1X, ES9.2, 1X)') "ERI_MME| Maximum error G, R:", &
          MAXVAL(I_diff_G), MAXVAL(I_diff_R)

        acc_check =  param%error_tot .GE. MAX(MAXVAL(I_diff_G), MAXVAL(I_diff_R))
        IF(.NOT.acc_check)&
           CPABORT("Actual error greater than upper bound estimate.")

      ENDIF
    ENDIF

    CALL cell_release(box)
    CALL eri_mme_release(param)

  END SUBROUTINE eri_mme_perf_acc_test

END MODULE eri_mme_test
