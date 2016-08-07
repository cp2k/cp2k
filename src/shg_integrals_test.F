!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates 2-center integrals for different r12 operators comparing the Solid harmonic
!>        Gaussian integral scheme to the Obara-Saika (OS) scheme
!> \author  Dorothea Golze [05.2016]
! **************************************************************************************************
MODULE shg_integrals_test

   USE basis_set_types,                 ONLY: allocate_gto_basis_set,&
                                              deallocate_gto_basis_set,&
                                              gto_basis_set_type,&
                                              init_orb_basis_set,&
                                              read_gto_basis_set
   USE generic_os_integrals,            ONLY: int_overlap_ab_os,&
                                              int_overlap_aba_os,&
                                              int_overlap_abb_os,&
                                              two_center_os_integrals
   USE generic_shg_integrals,           ONLY: int_overlap_ab_shg,&
                                              int_overlap_aba_shg,&
                                              int_overlap_abb_shg,&
                                              precalc_angular_shg_part,&
                                              two_center_shg_integrals
   USE generic_shg_integrals_init,      ONLY: contraction_matrix_shg,&
                                              contraction_matrix_shg_mix,&
                                              get_clebsch_gordon_coefficients
   USE input_constants,                 ONLY: operator_coulomb,&
                                              operator_gauss,&
                                              operator_verf,&
                                              operator_verfc,&
                                              operator_vgauss
   USE input_cp2k_subsys,               ONLY: create_basis_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: &
        section_add_keyword, section_add_subsection, section_create, section_release, &
        section_type, section_vals_get, section_vals_get_subs_vals, section_vals_type, &
        section_vals_val_get
   USE input_val_types,                 ONLY: real_t
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE orbital_pointers,                ONLY: init_orbital_pointers
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'shg_integrals_test'

   PUBLIC :: create_shg_integrals_test_section, shg_integrals_perf_acc_test

CONTAINS

! **************************************************************************************************
!> \brief Create input section for unit testing
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_shg_integrals_test_section(section)
      TYPE(section_type), INTENT(INOUT), POINTER         :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_shg_integrals_test_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      NULLIFY (keyword, subsection)

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "SHG_INTEGRALS_TEST", &
                          description="Parameters for testing the SHG 2-center integrals for "// &
                          "different r12 operators. Test w.r.t. performance and accurarcy.", &
                          n_keywords=4, n_subsections=1)

      CALL create_basis_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL keyword_create(keyword, &
                          name="_SECTION_PARAMETERS_", &
                          description="Controls the activation the SHG integral test. ", &
                          default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="ABC", &
                          description="Specify the lengths of the cell vectors A, B, and C. ", &
                          usage="ABC 10.000 10.000 10.000", unit_str="angstrom", &
                          n_var=3, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NAB_MIN", &
                          description="Minimum number of atomic distances to consider. ", &
                          default_i_val=8)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NREP", &
                          description="Number of repeated calculation of each integral. ", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="CHECK_ACCURACY", &
                          description="Causes abortion when SHG and OS integrals differ "// &
                          "more what's given by ACCURACY_LEVEL.", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="ACCURACY_LEVEL", &
                          description="Level of accuracy for comparison of SHG and OS "// &
                          "integrals.", &
                          default_r_val=1.0E-8_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="CALCULATE_DERIVATIVES", &
                          description="Calculates also the derivatives of the integrals.", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_OVERLAP", &
                          description="Calculates the integrals (a|b).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_release(keyword)
      CALL keyword_create(keyword, name="TEST_COULOMB", &
                          description="Calculates the integrals (a|1/r12|b).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_VERF", &
                          description="Calculates the integrals (a|erf(omega*r12)/r12|b).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_VERFC", &
                          description="Calculates the integrals (a|erfc(omega*r12)/r12|b).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_VGAUSS", &
                          description="Calculates the integrals (a|exp(omega*r12^2)/r12|b).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_GAUSS", &
                          description="Calculates the integrals (a|exp(omega*r12^2)|b).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_OVERLAP_ABA", &
                          description="Calculates the integrals (a|b|b).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="TEST_OVERLAP_ABB", &
                          description="Calculates the integrals (a|b|b).", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_shg_integrals_test_section

! **************************************************************************************************
!> \brief Unit test for performance and accuracy of the SHG integrals
!> \param iw output unit
!> \param shg_integrals_test_section ...
! **************************************************************************************************
   SUBROUTINE shg_integrals_perf_acc_test(iw, shg_integrals_test_section)
      INTEGER, INTENT(IN)                                :: iw
      TYPE(section_vals_type), INTENT(INOUT), POINTER    :: shg_integrals_test_section

      CHARACTER(len=*), PARAMETER :: routineN = 'shg_integrals_perf_acc_test', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: basis_type
      INTEGER                                            :: count_ab, handle, iab, jab, kab, lamax, &
                                                            lbmax, lcamax, lcbmax, lmax, nab, &
                                                            nab_min, nab_xyz, nrep, nrep_bas
      LOGICAL                                            :: acc_check, calc_derivatives, &
                                                            test_overlap_aba, test_overlap_abb
      REAL(KIND=dp)                                      :: acc_param
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cell_par
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scona_shg, sconb_shg
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb, oba, obb
      TYPE(section_vals_type), POINTER                   :: basis_section

      CALL timeset(routineN, handle)
      NULLIFY (oba, obb, fba, fbb, basis_section, cell_par)
      CALL section_vals_val_get(shg_integrals_test_section, "ABC", r_vals=cell_par)
      CALL section_vals_val_get(shg_integrals_test_section, "NAB_MIN", i_val=nab_min)
      CALL section_vals_val_get(shg_integrals_test_section, "NREP", i_val=nrep)
      CALL section_vals_val_get(shg_integrals_test_section, "CHECK_ACCURACY", l_val=acc_check)
      CALL section_vals_val_get(shg_integrals_test_section, "ACCURACY_LEVEL", r_val=acc_param)
      CALL section_vals_val_get(shg_integrals_test_section, "CALCULATE_DERIVATIVES", l_val=calc_derivatives)
      CALL section_vals_val_get(shg_integrals_test_section, "TEST_OVERLAP_ABA", l_val=test_overlap_aba)
      CALL section_vals_val_get(shg_integrals_test_section, "TEST_OVERLAP_ABB", l_val=test_overlap_abb)

      !*** Read the basis set information
      basis_section => section_vals_get_subs_vals(shg_integrals_test_section, "BASIS")
      CALL section_vals_get(basis_section, n_repetition=nrep_bas)
      IF (.NOT. (nrep_bas == 2 .OR. nrep_bas == 3)) THEN
         CALL cp_abort(__LOCATION__, &
                       "Provide basis sets")
      END IF
      CALL allocate_gto_basis_set(oba)
      CALL read_gto_basis_set(TRIM("A"), basis_type, oba, basis_section, irep=1)
      lamax = MAXVAL(oba%lmax)
      CALL allocate_gto_basis_set(obb)
      CALL read_gto_basis_set(TRIM("B"), basis_type, obb, basis_section, irep=2)
      lbmax = MAXVAL(obb%lmax)
      lmax = MAX(lamax, lbmax)
      IF (test_overlap_aba) THEN
         CALL allocate_gto_basis_set(fba)
         CALL read_gto_basis_set(TRIM("CA"), basis_type, fba, basis_section, irep=3)
         lcamax = MAXVAL(fba%lmax)
         lmax = MAX(lamax+lcamax, lbmax)
      ENDIF
      IF (test_overlap_abb) THEN
         CALL allocate_gto_basis_set(fbb)
         CALL read_gto_basis_set(TRIM("CB"), basis_type, fbb, basis_section, irep=3)
         lcbmax = MAXVAL(fbb%lmax)
         lmax = MAX(lamax, lbmax+lcbmax)
      ENDIF
      IF (test_overlap_aba .AND. test_overlap_abb) THEN
         lmax = MAX(MAX(lamax+lcamax, lbmax), MAX(lamax, lbmax+lcbmax))
      ENDIF
      !*** Initialize basis set information
      CALL init_orbital_pointers(lmax+1)
      CALL init_spherical_harmonics(lmax, output_unit=-100)
      oba%norm_type = 2
      CALL init_orb_basis_set(oba)
      obb%norm_type = 2
      CALL init_orb_basis_set(obb)
      IF (test_overlap_aba) THEN
         fba%norm_type = 2
         CALL init_orb_basis_set(fba)
      ENDIF
      IF (test_overlap_abb) THEN
         fbb%norm_type = 2
         CALL init_orb_basis_set(fbb)
      ENDIF
      ! if shg integrals are later actually used in the code, contraction_matrix_shg should be
      ! moved to init_orb_basis_set and scon_shg should become an element of gto_basis_set_type
      CALL contraction_matrix_shg(oba, scona_shg)
      CALL contraction_matrix_shg(obb, sconb_shg)

      !*** Create range of rab (atomic distances) to be tested
      nab_xyz = CEILING(REAL(nab_min, KIND=dp)**(1.0_dp/3.0_dp)-1.0E-06)
      nab = nab_xyz**3

      ALLOCATE (rab(3, nab))
      count_ab = 0
      DO iab = 1, nab_xyz
         DO jab = 1, nab_xyz
            DO kab = 1, nab_xyz
               count_ab = count_ab+1
               rab(:, count_ab) = [iab*ABS(cell_par(1)), jab*ABS(cell_par(2)), kab*ABS(cell_par(3))]/nab_xyz
            ENDDO
         ENDDO
      ENDDO

      !*** Calculate the SHG integrals

      CALL test_shg_operator12_integrals(oba, obb, rab, nrep, scona_shg, sconb_shg, &
                                         shg_integrals_test_section, acc_check, &
                                         acc_param, calc_derivatives, iw)

      CALL test_shg_overlap_integrals(oba, obb, rab, nrep, scona_shg, sconb_shg, &
                                      shg_integrals_test_section, acc_check, &
                                      acc_param, calc_derivatives, iw)

      CALL test_shg_overlap_aba_integrals(oba, obb, fba, fbb, rab, nrep, scona_shg, sconb_shg, &
                                          shg_integrals_test_section, acc_check, &
                                          acc_param, calc_derivatives, iw)

      DEALLOCATE (scona_shg, sconb_shg, rab)
      CALL deallocate_gto_basis_set(oba)
      CALL deallocate_gto_basis_set(obb)
      IF (test_overlap_aba) CALL deallocate_gto_basis_set(fba)
      IF (test_overlap_abb) CALL deallocate_gto_basis_set(fbb)

      CALL timestop(handle)

   END SUBROUTINE shg_integrals_perf_acc_test

! **************************************************************************************************
!> \brief ...
!> \param oba ...
!> \param obb ...
!> \param rab ...
!> \param nrep ...
!> \param scona_shg ...
!> \param sconb_shg ...
!> \param shg_integrals_test_section ...
!> \param acc_check ...
!> \param acc_param ...
!> \param calc_derivatives ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE test_shg_operator12_integrals(oba, obb, rab, nrep, scona_shg, sconb_shg, &
                                            shg_integrals_test_section, acc_check, &
                                            acc_param, calc_derivatives, iw)
      TYPE(gto_basis_set_type), POINTER                  :: oba, obb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rab
      INTEGER, INTENT(IN)                                :: nrep
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scona_shg, sconb_shg
      TYPE(section_vals_type), INTENT(IN), POINTER       :: shg_integrals_test_section
      LOGICAL, INTENT(IN)                                :: acc_check
      REAL(KIND=dp), INTENT(IN)                          :: acc_param
      LOGICAL, INTENT(IN)                                :: calc_derivatives
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'test_shg_operator12_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iab, irep, nab, nfa, nfb
      LOGICAL                                            :: test_any, test_coulomb, test_gauss, &
                                                            test_verf, test_verfc, test_vgauss
      REAL(KIND=dp) :: ddmax_coulomb, ddmax_gauss, ddmax_verf, ddmax_verfc, ddmax_vgauss, ddtemp, &
         dmax_coulomb, dmax_gauss, dmax_verf, dmax_verfc, dmax_vgauss, dtemp, omega
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vab_os, vab_shg
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dvab_os, dvab_shg

      CALL section_vals_val_get(shg_integrals_test_section, "TEST_COULOMB", l_val=test_coulomb)
      CALL section_vals_val_get(shg_integrals_test_section, "TEST_VERF", l_val=test_verf)
      CALL section_vals_val_get(shg_integrals_test_section, "TEST_VERFC", l_val=test_verfc)
      CALL section_vals_val_get(shg_integrals_test_section, "TEST_VGAUSS", l_val=test_vgauss)
      CALL section_vals_val_get(shg_integrals_test_section, "TEST_GAUSS", l_val=test_gauss)

      test_any = (test_coulomb .OR. test_verf .OR. test_verfc .OR. test_vgauss .OR. test_gauss)

      IF (test_any) THEN
         nfa = oba%nsgf
         nfb = obb%nsgf
         ALLOCATE (vab_shg(nfa, nfb), dvab_shg(nfa, nfb, 3))
         ALLOCATE (vab_os(nfa, nfb), dvab_os(nfa, nfb, 3))
         omega = 2.3_dp
         dmax_coulomb = 0.0_dp
         ddmax_coulomb = 0.0_dp
         dmax_verf = 0.0_dp
         ddmax_verf = 0.0_dp
         dmax_verfc = 0.0_dp
         ddmax_verfc = 0.0_dp
         dmax_vgauss = 0.0_dp
         ddmax_vgauss = 0.0_dp
         dmax_gauss = 0.0_dp
         ddmax_gauss = 0.0_dp

         nab = SIZE(rab, 2)
         DO irep = 1, nrep
            DO iab = 1, nab
               !*** Coulomb: (a|1/r12|b)
               IF (test_coulomb) THEN
                  CALL two_center_shg_integrals(operator_coulomb, vab_shg, dvab_shg, rab(:, iab), &
                                                oba, obb, scona_shg, sconb_shg, &
                                                calculate_forces=calc_derivatives)
                  CALL two_center_os_integrals(operator_coulomb, vab_os, dvab_os, rab(:, iab), oba, obb, &
                                               calculate_forces=calc_derivatives)
                  CALL calculate_deviation_ab(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
                  dmax_coulomb = MAX(dmax_coulomb, dtemp)
                  ddmax_coulomb = MAX(ddmax_coulomb, ddtemp)
               ENDIF
               !*** verf: (a|erf(omega*r12)/r12|b)
               IF (test_verf) THEN
                  CALL two_center_shg_integrals(operator_verf, vab_shg, dvab_shg, rab(:, iab), &
                                                oba, obb, scona_shg, sconb_shg, omega, &
                                                calc_derivatives)
                  CALL two_center_os_integrals(operator_verf, vab_os, dvab_os, rab(:, iab), oba, obb, &
                                               omega, calc_derivatives)
                  CALL calculate_deviation_ab(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
                  dmax_verf = MAX(dmax_verf, dtemp)
                  ddmax_verf = MAX(ddmax_verf, ddtemp)
               ENDIF
               !*** verfc: (a|erfc(omega*r12)/r12|b)
               IF (test_verfc) THEN
                  CALL two_center_shg_integrals(operator_verfc, vab_shg, dvab_shg, rab(:, iab), &
                                                oba, obb, scona_shg, sconb_shg, omega, &
                                                calc_derivatives)
                  CALL two_center_os_integrals(operator_verfc, vab_os, dvab_os, rab(:, iab), oba, obb, &
                                               omega, calc_derivatives)
                  CALL calculate_deviation_ab(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
                  dmax_verfc = MAX(dmax_verfc, dtemp)
                  ddmax_verfc = MAX(ddmax_verfc, ddtemp)
               ENDIF
               !*** vgauss: (a|exp(omega*r12^2)/r12|b)
               IF (test_vgauss) THEN
                  CALL two_center_shg_integrals(operator_vgauss, vab_shg, dvab_shg, rab(:, iab), &
                                                oba, obb, scona_shg, sconb_shg, omega, &
                                                calc_derivatives)
                  CALL two_center_os_integrals(operator_vgauss, vab_os, dvab_os, rab(:, iab), oba, obb, &
                                               omega, calc_derivatives)
                  CALL calculate_deviation_ab(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
                  dmax_vgauss = MAX(dmax_vgauss, dtemp)
                  ddmax_vgauss = MAX(ddmax_vgauss, ddtemp)
               ENDIF
               !*** gauss: (a|exp(omega*r12^2)|b)
               IF (test_gauss) THEN
                  CALL two_center_shg_integrals(operator_gauss, vab_shg, dvab_shg, rab(:, iab), &
                                                oba, obb, scona_shg, sconb_shg, omega, &
                                                calc_derivatives)
                  CALL two_center_os_integrals(operator_gauss, vab_os, dvab_os, rab(:, iab), oba, obb, &
                                               omega, calc_derivatives)
                  CALL calculate_deviation_ab(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
                  dmax_gauss = MAX(dmax_gauss, dtemp)
                  ddmax_gauss = MAX(ddmax_gauss, ddtemp)
               ENDIF
            ENDDO
         ENDDO

         IF (iw > 0) THEN
            WRITE (iw, FMT="(/,T2,A)") "TEST INFO FOR 2-CENTER SHG and OS INTEGRALS:"
            WRITE (iw, FMT="(T2,A)") "Maximal deviation between SHG and OS integrals and their derivatives"
            IF (test_coulomb) THEN
               WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|1/r12|b]", &
                  dmax_coulomb, ddmax_coulomb
            ENDIF
            IF (test_verf) THEN
               WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|erf(omega*r12)/r12|b]", &
                  dmax_verf, ddmax_verf
            ENDIF
            IF (test_verfc) THEN
               WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|erfc(omega*r12)/r12|b]", &
                  dmax_verfc, ddmax_verfc
            ENDIF
            IF (test_vgauss) THEN
               WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|exp(-omega*r12^2)/r12|b]", &
                  dmax_vgauss, ddmax_vgauss
            ENDIF
            IF (test_gauss) THEN
               WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|exp(-omega*r12^2)|b]", &
                  dmax_gauss, ddmax_gauss
            ENDIF

            IF (acc_check) THEN
               IF ((dmax_coulomb >= acc_param) .OR. (ddmax_coulomb >= acc_param)) THEN
                  CPABORT("[a|1/r12|b]: Deviation between SHG and OS larger than 1.0E-08")
               ENDIF
               IF ((dmax_verf >= acc_param) .OR. (ddmax_verf >= acc_param)) THEN
                  CPABORT("[a|erf(omega*r12)/r12|b]: Dev. between SHG and OS larger than 1.0E-08")
               ENDIF
               IF ((dmax_verfc >= acc_param) .OR. (ddmax_verfc >= acc_param)) THEN
                  CPABORT("[a|erfc(omega*r12)/r12|b]: Dev. between SHG and OS larger than 1.0E-08")
               ENDIF
               IF ((dmax_vgauss >= acc_param) .OR. (ddmax_vgauss >= acc_param)) THEN
                  CPABORT("[a|exp(-omega*r12^2)/r12|b]: Dev. between SHG and OS larger than 1.0E-08")
               ENDIF
               IF ((dmax_gauss >= acc_param) .OR. (ddmax_gauss >= acc_param)) THEN
                  CPABORT("[a|exp(-omega*r12^2)|b]: Dev. between SHG and OS larger than 1.0E-08")
               ENDIF
            ENDIF
         ENDIF
         DEALLOCATE (vab_shg, vab_os, dvab_shg, dvab_os)
      ENDIF

   END SUBROUTINE test_shg_operator12_integrals

! **************************************************************************************************
!> \brief ...
!> \param oba ...
!> \param obb ...
!> \param rab ...
!> \param nrep ...
!> \param scona_shg ...
!> \param sconb_shg ...
!> \param shg_integrals_test_section ...
!> \param acc_check ...
!> \param acc_param ...
!> \param calc_derivatives ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE test_shg_overlap_integrals(oba, obb, rab, nrep, scona_shg, sconb_shg, &
                                         shg_integrals_test_section, acc_check, &
                                         acc_param, calc_derivatives, iw)
      TYPE(gto_basis_set_type), POINTER                  :: oba, obb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rab
      INTEGER, INTENT(IN)                                :: nrep
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scona_shg, sconb_shg
      TYPE(section_vals_type), INTENT(IN), POINTER       :: shg_integrals_test_section
      LOGICAL, INTENT(IN)                                :: acc_check
      REAL(KIND=dp), INTENT(IN)                          :: acc_param
      LOGICAL, INTENT(IN)                                :: calc_derivatives
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'test_shg_overlap_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iab, irep, la_max, lb_max, nab, nfa, nfb
      LOGICAL                                            :: test_overlap
      REAL(KIND=dp)                                      :: ddmax_overlap, ddtemp, dmax_overlap, &
                                                            dtemp, dummy, ra(3), rb(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Waux_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sab_os, sab_shg
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dsab_os, dsab_shg

      CALL section_vals_val_get(shg_integrals_test_section, "TEST_OVERLAP", &
                                l_val=test_overlap)
      IF (test_overlap) THEN
         !effectively switch off screening; makes no sense for the tests
         oba%set_radius(:) = 1.0E+09_dp
         obb%set_radius(:) = 1.0E+09_dp
         oba%pgf_radius(:, :) = 1.0E+09_dp
         obb%pgf_radius(:, :) = 1.0E+09_dp
         nfa = oba%nsgf
         nfb = obb%nsgf
         la_max = MAXVAL(oba%lmax)
         lb_max = MAXVAL(obb%lmax)
         ra(:) = 0.0_dp
         rb(:) = 0.0_dp
         dummy = 0.0_dp
         dmax_overlap = 0.0_dp
         ddmax_overlap = 0.0_dp
         ALLOCATE (sab_shg(nfa, nfb), dsab_shg(nfa, nfb, 3))
         ALLOCATE (sab_os(nfa, nfb), dsab_os(nfa, nfb, 3))
         nab = SIZE(rab, 2)
         DO irep = 1, nrep
            DO iab = 1, nab
               CALL precalc_angular_shg_part(la_max, lb_max, rab(:, iab), Waux_mat, dWaux_mat, &
                                             calc_derivatives)
               CALL int_overlap_ab_shg(sab_shg, dsab_shg, rab(:, iab), oba, obb, scona_shg, &
                                       sconb_shg, Waux_mat, dWaux_mat, calc_derivatives, &
                                       contraction_high=.TRUE.)
               CALL int_overlap_ab_os(sab_os, dsab_os, ra, rb, rab(:, iab), oba, obb, &
                                      calc_derivatives, debug=.FALSE., dmax=dummy)
               CALL calculate_deviation_ab(sab_shg, sab_os, dsab_shg, dsab_os, dtemp, ddtemp)
               dmax_overlap = MAX(dmax_overlap, dtemp)
               ddmax_overlap = MAX(ddmax_overlap, ddtemp)
            ENDDO
         ENDDO

         IF (iw > 0) THEN
            WRITE (iw, FMT="(/,T2,A)") "TEST INFO FOR 2-CENTER OVERLAP SHG and OS INTEGRALS:"
            WRITE (iw, FMT="(T2,A)") "Maximal deviation between SHG and OS integrals and their derivatives"
            WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|b]", &
               dmax_overlap, ddmax_overlap
         ENDIF
         IF (acc_check) THEN
            IF ((dmax_overlap >= acc_param) .OR. (ddmax_overlap >= acc_param)) THEN
               CPABORT("[a|b]: Deviation between SHG and OS larger than 1.0E-08")
            ENDIF
         ENDIF
         DEALLOCATE (Waux_mat, dWaux_mat)
         DEALLOCATE (sab_shg, sab_os, dsab_shg, dsab_os)
      ENDIF

   END SUBROUTINE test_shg_overlap_integrals

! **************************************************************************************************
!> \brief ...
!> \param oba ...
!> \param obb ...
!> \param fba ...
!> \param fbb ...
!> \param rab ...
!> \param nrep ...
!> \param scon_oba ...
!> \param scon_obb ...
!> \param shg_integrals_test_section ...
!> \param acc_check ...
!> \param acc_param ...
!> \param calc_derivatives ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE test_shg_overlap_aba_integrals(oba, obb, fba, fbb, rab, nrep, scon_oba, scon_obb, &
                                             shg_integrals_test_section, acc_check, &
                                             acc_param, calc_derivatives, iw)
      TYPE(gto_basis_set_type), POINTER                  :: oba, obb, fba, fbb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: rab
      INTEGER, INTENT(IN)                                :: nrep
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scon_oba, scon_obb
      TYPE(section_vals_type), INTENT(IN), POINTER       :: shg_integrals_test_section
      LOGICAL, INTENT(IN)                                :: acc_check
      REAL(KIND=dp), INTENT(IN)                          :: acc_param
      LOGICAL, INTENT(IN)                                :: calc_derivatives
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'test_shg_overlap_aba_integrals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iab, irep, la_max, laa_max, lb_max, &
                                                            lbb_max, maxl_orb, maxl_ri, nab, nba, &
                                                            nbb, nfa, nfb
      INTEGER, DIMENSION(:, :), POINTER                  :: ncg_none0
      INTEGER, DIMENSION(:, :, :), POINTER               :: cg_none0_list, fba_index, fbb_index, &
                                                            oba_index, obb_index
      LOGICAL                                            :: test_overlap_aba, test_overlap_abb
      REAL(KIND=dp)                                      :: ddmax_overlap_aba, ddmax_overlap_abb, &
                                                            ddtemp, dmax_overlap_aba, &
                                                            dmax_overlap_abb, dtemp, dummy, ra(3), &
                                                            rb(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: Waux_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: dWaux_mat
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: cg_coeff, saba_os, saba_shg, sabb_os, &
                                                            sabb_shg
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: dsaba_os, dsaba_shg, dsabb_os, &
                                                            dsabb_shg, scona_mix, sconb_mix

      CALL section_vals_val_get(shg_integrals_test_section, "TEST_OVERLAP_ABA", &
                                l_val=test_overlap_aba)
      CALL section_vals_val_get(shg_integrals_test_section, "TEST_OVERLAP_ABB", &
                                l_val=test_overlap_abb)
      IF (test_overlap_aba .OR. test_overlap_abb) THEN
         !effectively switch off screening; makes no sense for the tests
         oba%set_radius(:) = 1.0E+09_dp
         obb%set_radius(:) = 1.0E+09_dp
         oba%pgf_radius(:, :) = 1.0E+09_dp
         obb%pgf_radius(:, :) = 1.0E+09_dp
         nba = oba%nsgf
         nbb = obb%nsgf
         maxl_orb = MAX(MAXVAL(oba%lmax), MAXVAL(obb%lmax))
         la_max = MAXVAL(oba%lmax)
         lb_max = MAXVAL(obb%lmax)
         IF (test_overlap_aba) THEN
            fba%set_radius(:) = 1.0E+09_dp
            fba%pgf_radius(:, :) = 1.0E+09_dp
            nfa = fba%nsgf
            maxl_ri = MAXVAL(fba%lmax)
            laa_max = MAXVAL(oba%lmax)+MAXVAL(fba%lmax)
            ALLOCATE (saba_shg(nba, nbb, nfa), dsaba_shg(nba, nbb, nfa, 3))
            ALLOCATE (saba_os(nba, nbb, nfa), dsaba_os(nba, nbb, nfa, 3))
            CALL contraction_matrix_shg_mix(oba, fba, oba_index, fba_index, scona_mix)
         ENDIF
         IF (test_overlap_abb) THEN
            fbb%set_radius(:) = 1.0E+09_dp
            fbb%pgf_radius(:, :) = 1.0E+09_dp
            nfb = fbb%nsgf
            maxl_ri = MAXVAL(fbb%lmax)
            lbb_max = MAXVAL(obb%lmax)+MAXVAL(fbb%lmax)
            ALLOCATE (sabb_shg(nba, nbb, nfb), dsabb_shg(nba, nbb, nfb, 3))
            ALLOCATE (sabb_os(nba, nbb, nfb), dsabb_os(nba, nbb, nfb, 3))
            CALL contraction_matrix_shg_mix(obb, fbb, obb_index, fbb_index, sconb_mix)
         ENDIF
         ra(:) = 0.0_dp
         rb(:) = 0.0_dp
         dummy = 0.0_dp
         dmax_overlap_aba = 0.0_dp
         ddmax_overlap_aba = 0.0_dp
         dmax_overlap_abb = 0.0_dp
         ddmax_overlap_abb = 0.0_dp
         CALL get_clebsch_gordon_coefficients(cg_coeff, cg_none0_list, ncg_none0, maxl_orb, maxl_ri)
         nab = SIZE(rab, 2)
         IF (test_overlap_aba) THEN
            DO irep = 1, nrep
               DO iab = 1, nab
                  saba_shg = 0.0_dp
                  dsaba_shg = 0.0_dp
                  CALL precalc_angular_shg_part(laa_max, lb_max, rab(:, iab), Waux_mat, dWaux_mat, &
                                                calc_derivatives)
                  CALL int_overlap_aba_shg(saba_shg, dsaba_shg, rab(:, iab), oba, obb, fba, &
                                           scon_obb, scona_mix, oba_index, fba_index, &
                                           cg_coeff, cg_none0_list, ncg_none0, &
                                           Waux_mat, dWaux_mat, calc_derivatives)
                  CALL int_overlap_aba_os(saba_os, dsaba_os, ra, rb, rab(:, iab), oba, obb, fba, &
                                          calc_derivatives, debug=.FALSE., dmax=dummy)
                  CALL calculate_deviation_abx(saba_shg, saba_os, dsaba_shg, dsaba_os, dtemp, ddtemp)
                  dmax_overlap_aba = MAX(dmax_overlap_aba, dtemp)
                  ddmax_overlap_aba = MAX(ddmax_overlap_aba, ddtemp)
                  DEALLOCATE (Waux_mat, dWaux_mat)
               ENDDO
            ENDDO
            DEALLOCATE (oba_index, fba_index, scona_mix)
            DEALLOCATE (saba_shg, saba_os, dsaba_shg, dsaba_os)
         ENDIF
         IF (test_overlap_abb) THEN
            DO irep = 1, nrep
               DO iab = 1, nab
                  sabb_shg = 0.0_dp
                  dsabb_shg = 0.0_dp
                  CALL precalc_angular_shg_part(lbb_max, la_max, rab(:, iab), Waux_mat, dWaux_mat, &
                                                calc_derivatives)
                  CALL int_overlap_abb_shg(sabb_shg, dsabb_shg, rab(:, iab), oba, obb, fbb, &
                                           scon_oba, sconb_mix, obb_index, fbb_index, &
                                           cg_coeff, cg_none0_list, ncg_none0, &
                                           Waux_mat, dWaux_mat, calc_derivatives)
                  CALL int_overlap_abb_os(sabb_os, dsabb_os, ra, rb, rab(:, iab), oba, obb, fbb, &
                                          calc_derivatives, debug=.FALSE., dmax=dummy)
                  CALL calculate_deviation_abx(sabb_shg, sabb_os, dsabb_shg, dsabb_os, dtemp, ddtemp)
                  dmax_overlap_abb = MAX(dmax_overlap_abb, dtemp)
                  ddmax_overlap_abb = MAX(ddmax_overlap_abb, ddtemp)
                  DEALLOCATE (Waux_mat, dWaux_mat)
               ENDDO
            ENDDO
            DEALLOCATE (obb_index, fbb_index, sconb_mix)
            DEALLOCATE (sabb_shg, sabb_os, dsabb_shg, dsabb_os)
         ENDIF
         IF (iw > 0) THEN
            WRITE (iw, FMT="(/,T2,A)") "TEST INFO [a|b|x] OVERLAP SHG and OS INTEGRALS:"
            WRITE (iw, FMT="(T2,A)") "Maximal deviation between SHG and OS integrals and their derivatives"
            IF (test_overlap_aba) THEN
               WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|b|a]", &
                  dmax_overlap_aba, ddmax_overlap_aba
            ENDIF
            IF (test_overlap_abb) THEN
               WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|b|b]", &
                  dmax_overlap_abb, ddmax_overlap_abb
            ENDIF
         ENDIF
         IF (acc_check) THEN
            IF ((dmax_overlap_aba >= acc_param) .OR. (ddmax_overlap_aba >= acc_param)) THEN
               CPABORT("[a|b|a]: Deviation between SHG and OS larger than 1.0E-08")
            ENDIF
            IF ((dmax_overlap_abb >= acc_param) .OR. (ddmax_overlap_abb >= acc_param)) THEN
               CPABORT("[a|b|b]: Deviation between SHG and OS larger than 1.0E-08")
            ENDIF
         ENDIF
         DEALLOCATE (cg_coeff, cg_none0_list, ncg_none0)
      ENDIF

   END SUBROUTINE test_shg_overlap_aba_integrals

! **************************************************************************************************
!> \brief Calculation of the deviation between SHG and OS integrals
!> \param vab_shg integral matrix obtained from the SHG scheme
!> \param vab_os integral matrix obtained from the OS scheme
!> \param dvab_shg derivative of the integrals, SHG
!> \param dvab_os derivative of the integrals, OS
!> \param dmax maximal deviation of vab matrices
!> \param ddmax maximal deviation of dvab matrices
! **************************************************************************************************
   SUBROUTINE calculate_deviation_ab(vab_shg, vab_os, dvab_shg, dvab_os, dmax, ddmax)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: vab_shg, vab_os
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: dvab_shg, dvab_os
      REAL(KIND=dp), INTENT(OUT)                         :: dmax, ddmax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_deviation_ab', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: diff

      dmax = 0.0_dp
      ddmax = 0.0_dp

      ! integrals vab
      DO j = 1, SIZE(vab_shg, 2)
         DO i = 1, SIZE(vab_shg, 1)
            diff = ABS(vab_shg(i, j)-vab_os(i, j))
            dmax = MAX(dmax, diff)
         ENDDO
      ENDDO

      ! derivatives dvab
      DO k = 1, 3
         DO j = 1, SIZE(dvab_shg, 2)
            DO i = 1, SIZE(dvab_shg, 1)
               diff = ABS(dvab_shg(i, j, k)-dvab_os(i, j, k))
               ddmax = MAX(ddmax, diff)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE calculate_deviation_ab

! **************************************************************************************************
!> \brief Calculation of the deviation between SHG and OS integrals
!> \param vab_shg integral matrix obtained from the SHG scheme
!> \param vab_os integral matrix obtained from the OS scheme
!> \param dvab_shg derivative of the integrals, SHG
!> \param dvab_os derivative of the integrals, OS
!> \param dmax maximal deviation of vab matrices
!> \param ddmax maximal deviation of dvab matrices
! **************************************************************************************************
   SUBROUTINE calculate_deviation_abx(vab_shg, vab_os, dvab_shg, dvab_os, dmax, ddmax)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: vab_shg, vab_os
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN)   :: dvab_shg, dvab_os
      REAL(KIND=dp), INTENT(OUT)                         :: dmax, ddmax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_deviation_abx', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, k, l
      REAL(KIND=dp)                                      :: diff

      dmax = 0.0_dp
      ddmax = 0.0_dp

      ! integrals vab
      DO k = 1, SIZE(vab_shg, 3)
         DO j = 1, SIZE(vab_shg, 2)
            DO i = 1, SIZE(vab_shg, 1)
               diff = ABS(vab_shg(i, j, k)-vab_os(i, j, k))
               dmax = MAX(dmax, diff)
            ENDDO
         ENDDO
      ENDDO

      ! derivatives dvab
      DO l = 1, 3
         DO k = 1, SIZE(dvab_shg, 3)
            DO j = 1, SIZE(dvab_shg, 2)
               DO i = 1, SIZE(dvab_shg, 1)
                  diff = ABS(dvab_shg(i, j, k, l)-dvab_os(i, j, k, l))
                  ddmax = MAX(ddmax, diff)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE calculate_deviation_abx
END MODULE shg_integrals_test
