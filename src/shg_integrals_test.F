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
   USE generic_os_integrals,            ONLY: two_center_os_integrals
   USE generic_shg_integrals,           ONLY: two_center_shg_integrals
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
   USE mathconstants,                   ONLY: fac,&
                                              pi
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
                          "more than 1.0E-08.", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
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
      INTEGER                                            :: count_ab, iab, irep, jab, kab, lamax, &
                                                            lbmax, lmax, nab, nab_min, nab_xyz, &
                                                            nfa, nfb, nrep, nrep_bas
      LOGICAL                                            :: acc_check
      REAL(KIND=dp) :: acc_param, ddmax_coulomb, ddmax_gauss, ddmax_verf, ddmax_verfc, &
         ddmax_vgauss, ddtemp, dmax_coulomb, dmax_gauss, dmax_verf, dmax_verfc, dmax_vgauss, &
         dtemp, omega
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cell_par
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vab_os, vab_shg
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dvab_os, dvab_shg, scona_shg, sconb_shg
      TYPE(gto_basis_set_type), POINTER                  :: fba, fbb
      TYPE(section_vals_type), POINTER                   :: basis_section

      NULLIFY (fba, fbb, basis_section, cell_par)
      CALL section_vals_val_get(shg_integrals_test_section, "ABC", r_vals=cell_par)
      CALL section_vals_val_get(shg_integrals_test_section, "NAB_MIN", i_val=nab_min)
      CALL section_vals_val_get(shg_integrals_test_section, "NREP", i_val=nrep)
      CALL section_vals_val_get(shg_integrals_test_section, "CHECK_ACCURACY", l_val=acc_check)

      !*** Read the basis set information
      basis_section => section_vals_get_subs_vals(shg_integrals_test_section, "BASIS")
      CALL section_vals_get(basis_section, n_repetition=nrep_bas)
      IF (nrep_bas /= 2) RETURN
      CALL allocate_gto_basis_set(fba)
      CALL read_gto_basis_set(TRIM("A"), basis_type, fba, basis_section, irep=1)
      lamax = MAXVAL(fba%lmax)
      CALL allocate_gto_basis_set(fbb)
      CALL read_gto_basis_set(TRIM("B"), basis_type, fbb, basis_section, irep=2)
      lbmax = MAXVAL(fbb%lmax)
      lmax = MAX(lamax, lbmax)

      !*** Initialize basis set information
      CALL init_orbital_pointers(lmax+1)
      CALL init_spherical_harmonics(lmax, output_unit=-100)
      fba%norm_type = 2
      CALL init_orb_basis_set(fba)
      fbb%norm_type = 2
      CALL init_orb_basis_set(fbb)
      ! if shg integrals are later actually used in the code, contraction_matrix_shg should be
      ! moved to init_orb_basis_set and scon_shg should become an element of gto_basis_set_type
      CALL contraction_matrix_shg(fba, scona_shg)
      CALL contraction_matrix_shg(fbb, sconb_shg)

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
      nfa = fba%nsgf
      nfb = fbb%nsgf
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
      DO irep = 1, nrep
         DO iab = 1, nab
            !*** Coulomb: (a|1/r12|b)
            CALL two_center_shg_integrals(operator_coulomb, vab_shg, dvab_shg, rab(:, iab), &
                                          fba, fbb, scona_shg, sconb_shg, calculate_forces=.TRUE.)
            CALL two_center_os_integrals(operator_coulomb, vab_os, dvab_os, rab(:, iab), fba, fbb, &
                                         calculate_forces=.TRUE.)
            CALL calculate_deviation(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
            dmax_coulomb = MAX(dmax_coulomb, dtemp)
            ddmax_coulomb = MAX(ddmax_coulomb, ddtemp)
            !*** verf: (a|erf(omega*r12)/r12|b)
            CALL two_center_shg_integrals(operator_verf, vab_shg, dvab_shg, rab(:, iab), &
                                          fba, fbb, scona_shg, sconb_shg, omega, &
                                          calculate_forces=.TRUE.)
            CALL two_center_os_integrals(operator_verf, vab_os, dvab_os, rab(:, iab), fba, fbb, &
                                         omega, calculate_forces=.TRUE.)
            CALL calculate_deviation(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
            dmax_verf = MAX(dmax_verf, dtemp)
            ddmax_verf = MAX(ddmax_verf, ddtemp)
            !*** verfc: (a|erfc(omega*r12)/r12|b)
            CALL two_center_shg_integrals(operator_verfc, vab_shg, dvab_shg, rab(:, iab), &
                                          fba, fbb, scona_shg, sconb_shg, omega, &
                                          calculate_forces=.TRUE.)
            CALL two_center_os_integrals(operator_verfc, vab_os, dvab_os, rab(:, iab), fba, fbb, &
                                         omega, calculate_forces=.TRUE.)
            CALL calculate_deviation(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
            dmax_verfc = MAX(dmax_verfc, dtemp)
            ddmax_verfc = MAX(ddmax_verfc, ddtemp)
            !*** vgauss: (a|exp(omega*r12^2)/r12|b)
            CALL two_center_shg_integrals(operator_vgauss, vab_shg, dvab_shg, rab(:, iab), &
                                          fba, fbb, scona_shg, sconb_shg, omega, &
                                          calculate_forces=.TRUE.)
            CALL two_center_os_integrals(operator_vgauss, vab_os, dvab_os, rab(:, iab), fba, fbb, &
                                         omega, calculate_forces=.TRUE.)
            CALL calculate_deviation(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
            dmax_vgauss = MAX(dmax_vgauss, dtemp)
            ddmax_vgauss = MAX(ddmax_vgauss, ddtemp)
            !*** gauss: (a|exp(omega*r12^2)|b)
            CALL two_center_shg_integrals(operator_gauss, vab_shg, dvab_shg, rab(:, iab), &
                                          fba, fbb, scona_shg, sconb_shg, omega, &
                                          calculate_forces=.TRUE.)
            CALL two_center_os_integrals(operator_gauss, vab_os, dvab_os, rab(:, iab), fba, fbb, &
                                         omega, calculate_forces=.TRUE.)
            CALL calculate_deviation(vab_shg, vab_os, dvab_shg, dvab_os, dtemp, ddtemp)
            dmax_gauss = MAX(dmax_gauss, dtemp)
            ddmax_gauss = MAX(ddmax_gauss, ddtemp)
         ENDDO
      ENDDO

      IF (iw > 0) THEN
         WRITE (iw, FMT="(/,T2,A)") "TEST INFO FOR 2-CENTER SHG and OS INTEGRALS:"
         WRITE (iw, FMT="(T2,A)") "Maximal deviation between SHG and OS integrals and their derivatives"
         WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|1/r12|b]", &
            dmax_coulomb, ddmax_coulomb
         WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|erf(omega*r12)/r12|b]", &
            dmax_verf, ddmax_verf
         WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|erfc(omega*r12)/r12|b]", &
            dmax_verfc, ddmax_verfc
         WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|exp(-omega*r12^2)/r12|b]", &
            dmax_vgauss, ddmax_vgauss
         WRITE (iw, FMT="(T2,A,T53,ES12.5,4X,ES12.5)") "SHG_INTEGRALS | [a|exp(-omega*r12^2)|b]", &
            dmax_gauss, ddmax_gauss

         IF (acc_check) THEN
            acc_param = 1.0E-08_dp
            IF ((dmax_coulomb >= acc_param) .OR. (ddmax_coulomb >= acc_param)) THEN
               CPABORT("[a|1/r12|b]: Deviation between SHG and OS larger than 1.0E-08")
            ENDIF
            IF ((dmax_verf >= acc_param) .OR. (ddmax_verf >= acc_param)) THEN
               CPABORT("[a|erf(omega*r12)/r12|b]: Deviation between SHG and OS larger than 1.0E-08")
            ENDIF
            IF ((dmax_verfc >= acc_param) .OR. (ddmax_verfc >= acc_param)) THEN
               CPABORT("[a|erfc(omega*r12)/r12|b]: Deviation between SHG and OS larger than 1.0E-08")
            ENDIF
            IF ((dmax_vgauss >= acc_param) .OR. (ddmax_vgauss >= acc_param)) THEN
               CPABORT("[a|exp(-omega*r12^2)/r12|b]: Dev. between SHG and OS larger than 1.0E-08")
            ENDIF
            IF ((dmax_gauss >= acc_param) .OR. (ddmax_gauss >= acc_param)) THEN
               CPABORT("[a|exp(-omega*r12^2)|b]: Deviation between SHG and OS larger than 1.0E-08")
            ENDIF
         ENDIF
      ENDIF
      DEALLOCATE (scona_shg, sconb_shg, rab, vab_shg, dvab_shg, vab_os, dvab_os)
      CALL deallocate_gto_basis_set(fba)
      CALL deallocate_gto_basis_set(fbb)

   END SUBROUTINE shg_integrals_perf_acc_test

! **************************************************************************************************
!> \brief contraction matrix for the SHG integrals. If these integrals are actually used,
!>        this routine should be moved, e.g. in basis_set_types.F, where the basis set
!>        information are initialized
!> \param basis ...
!> \param scon_shg contraction matrix
! **************************************************************************************************
   SUBROUTINE contraction_matrix_shg(basis, scon_shg)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: scon_shg

      CHARACTER(len=*), PARAMETER :: routineN = 'contraction_matrix_shg', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, ishell, l, maxpgf, maxshell, &
                                                            nset
      INTEGER, DIMENSION(:), POINTER                     :: npgf, nshell
      REAL(KIND=dp)                                      :: aif, gcc, norm_w
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: norm, zet

      nset = basis%nset
      npgf => basis%npgf
      nshell => basis%nshell
      zet => basis%zet

      maxpgf = SIZE(basis%gcc, 1)
      maxshell = SIZE(basis%gcc, 2)
      ALLOCATE (norm(basis%nset, maxshell))
      ALLOCATE (scon_shg(maxpgf, maxshell, nset))
      scon_shg = 0.0_dp

      CALL basis_norm_shg(basis, norm)

      DO iset = 1, nset
         DO ishell = 1, nshell(iset)
            l = basis%l(ishell, iset)
            norm_w = norm(iset, ishell)
            DO ipgf = 1, npgf(iset)
               aif = 1.0_dp/((2._dp*zet(ipgf, iset))**l)
               gcc = basis%gcc(ipgf, ishell, iset)
               scon_shg(ipgf, ishell, iset) = norm(iset, ishell)*gcc*aif
            END DO
         END DO
      END DO

      DEALLOCATE (norm)

   END SUBROUTINE contraction_matrix_shg

!*****************************************************************************
!> \brief normalization solid harmonic Gaussians (SHG)
!> \param basis ...
!> \param norm ...
! **************************************************************************************************
   SUBROUTINE basis_norm_shg(basis, norm)

      TYPE(gto_basis_set_type), POINTER                  :: basis
      REAL(dp), DIMENSION(:, :), POINTER                 :: norm

      CHARACTER(len=*), PARAMETER :: routineN = 'basis_norm_shg', routineP = moduleN//':'//routineN

      INTEGER                                            :: ipgf, iset, ishell, jpgf, l
      REAL(KIND=dp)                                      :: aai, aaj, cci, ccj, expa, ppl

      norm = 0._dp

      DO iset = 1, basis%nset
         DO ishell = 1, basis%nshell(iset)
            l = basis%l(ishell, iset)
            expa = 0.5_dp*REAL(2*l+3, dp)
            ppl = fac(2*l+2)*pi**(1.5_dp)/fac(l+1)
            ppl = ppl/(2._dp**REAL(2*l+1, dp))
            ppl = ppl/REAL(2*l+1, dp)
            DO ipgf = 1, basis%npgf(iset)
               cci = basis%gcc(ipgf, ishell, iset)
               aai = basis%zet(ipgf, iset)
               DO jpgf = 1, basis%npgf(iset)
                  ccj = basis%gcc(jpgf, ishell, iset)
                  aaj = basis%zet(jpgf, iset)
                  norm(iset, ishell) = norm(iset, ishell)+cci*ccj*ppl/(aai+aaj)**expa
               END DO
            END DO
            norm(iset, ishell) = 1.0_dp/SQRT(norm(iset, ishell))
         END DO
      END DO

   END SUBROUTINE basis_norm_shg

! **************************************************************************************************
!> \brief Calculation of the deviation between SHG and OS integrals
!> \param vab_shg integral matrix obtained from the SHG scheme
!> \param vab_os integral matrix obtained from the OS scheme
!> \param dvab_shg derivative of the integrals, SHG
!> \param dvab_os derivative of the integrals, OS
!> \param dmax maximal deviation of vab matrices
!> \param ddmax maximal deviation of dvab matrices
! **************************************************************************************************
   SUBROUTINE calculate_deviation(vab_shg, vab_os, dvab_shg, dvab_os, dmax, ddmax)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: vab_shg, vab_os
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: dvab_shg, dvab_os
      REAL(KIND=dp), INTENT(OUT)                         :: dmax, ddmax

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_deviation', &
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

   END SUBROUTINE
END MODULE shg_integrals_test
